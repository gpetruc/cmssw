#include "L1TriggerScouting/TauTagging/plugins/alpaka/TransformKernel.h"

#include "HeterogeneousCore/AlpakaInterface/interface/HistoContainer.h"
#include "HeterogeneousCore/AlpakaInterface/interface/radixSort.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"


namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels {

  using namespace cms::alpakatools;

  template <typename TAcc, typename T>
  ALPAKA_FN_ACC T px(const TAcc& acc, T pt, T phi) {
    return pt * alpaka::math::cos(acc, phi);
  }

  template <typename TAcc, typename T>
  ALPAKA_FN_ACC T py(const TAcc& acc, T pt, T phi) {
    return pt * alpaka::math::sin(acc, phi);
  }

  template <typename TAcc, typename T>
  ALPAKA_FN_ACC T pz(const TAcc& acc, T pt, T eta) {
    return pt * alpaka::math::sinh(acc, eta);
  }

  template <typename TAcc, typename T>
  ALPAKA_FN_ACC T energy(const TAcc& acc, T px, T py, T pz, T mass = 0.13957f) {
    return alpaka::math::sqrt(acc, 
        alpaka::math::pow(acc, px, 2) + 
        alpaka::math::pow(acc, py, 2) + 
        alpaka::math::pow(acc, pz, 2) + 
        alpaka::math::pow(acc, mass, 2));
  }

  template <typename TAcc, typename T>
  ALPAKA_FN_ACC T jet_pt(const TAcc& acc, T Px, T Py) {
    return alpaka::math::sqrt(acc, alpaka::math::pow(acc, Px, 2) + alpaka::math::pow(acc, Py, 2));
  }

  template <typename TAcc, typename T>
  ALPAKA_FN_ACC T jet_eta(const TAcc& acc, T Pt, T Pz) {
    return (Pt > 0.0) ? alpaka::math::asinh(acc, Pz / Pt) : 0.0; // ask simone later 
  }

  template <typename TAcc, typename T>
  ALPAKA_FN_ACC T jet_phi(const TAcc& acc, T Px, T Py) {
    return alpaka::math::atan2(acc, Py, Px);
  }

  ALPAKA_FN_ACC float phi(Acc1D const& acc, float phi, float Phi) {
    auto kPi = alpaka::math::constants::pi;
    return alpaka::math::remainder(acc, phi - Phi + kPi, 2.0 * kPi) - kPi;
  }

  template <typename TAcc>
  ALPAKA_FN_ACC float charge(const TAcc& acc, int pdgid) {
    auto charge = (alpaka::math::abs(acc, pdgid) == 211) ? ((pdgid > 0) ? +1 : -1) : ((pdgid > 0) ? -1 : +1);
    return charge;
  }

  class ComputeClueTauFeaturesKernel {
  public:
    ALPAKA_FN_ACC void operator()(
        Acc1D const& acc,
        PFCandidateDeviceCollection::ConstView pf,
        LongIndexSoA::ConstView indexes, 
        ClusterOffsetsSoA::ConstView offsets,
        SoftTauInputDeviceTensor::View clue_taus) const {
      for (auto block_idx: independent_groups(acc, offsets.metadata().size() - 1)) {
        auto begin = offsets.offsets()[block_idx];
        auto end = offsets.offsets()[block_idx + 1];
        auto block_dim = end - begin;
        if (block_dim == 0)
          continue;

        auto E = 0.0f;
        auto Px = 0.0f;
        auto Py = 0.0f;
        auto Pz = 0.0f;
        // jet kinematics
        auto Pt = 0.0f;
        auto Eta = 0.0f;
        auto Phi = 0.0f;
        // auto mass = 0.0f;
        if (once_per_block(acc)) {
          for (auto i = 0; i < block_dim; i++) {
            auto idx = indexes.indexes()[i+begin];
            auto px_v = px(acc, pf.pt()[idx], pf.phi()[idx]);
            auto py_v = py(acc, pf.pt()[idx], pf.phi()[idx]);
            auto pz_v = pz(acc, pf.pt()[idx], pf.eta()[idx]);
            auto e_v = energy(acc, px_v, py_v, pz_v);
            Px += px_v;
            Py += py_v;
            Pz += pz_v;
            E += e_v;
          }

          Pt = jet_pt(acc, Px, Py);
          Eta = jet_eta(acc, Pt, Pz);
          Phi = jet_phi(acc, Py, Px);
        }

        auto clue_tau = clue_taus[block_idx];
        auto cluster_size = (block_dim > JetFeatures::RowsAtCompileTime) ? JetFeatures::RowsAtCompileTime : block_dim;
        for (auto tid : independent_group_elements(acc, cluster_size)) {
          auto thread_idx = tid + begin; 
          auto index = indexes.indexes()[thread_idx];

          auto pdgid_abs = alpaka::math::abs(acc, static_cast<int>(pf.pdgid()[index]));
          // features
          clue_tau.features()(tid, 0) = pf.pt()[index];
          clue_tau.features()(tid, 1) = pf.eta()[index] - Eta;
          clue_tau.features()(tid, 2) = phi(acc, pf.phi()[index], Phi);
          clue_tau.features()(tid, 3) = charge(acc, pf.pdgid()[index]);
          clue_tau.features()(tid, 4) = pf.z0()[index];
          // one hot-encoding from pdgid
          clue_tau.features()(tid, 5) = (pdgid_abs == 211) ? 1.0f : 0.0f;
          clue_tau.features()(tid, 6) = (pdgid_abs == 130) ? 1.0f : 0.0f;
          clue_tau.features()(tid, 7) = (pdgid_abs == 11) ? 1.0f : 0.0f;
          clue_tau.features()(tid, 8) = (pdgid_abs == 13) ? 1.0f : 0.0f;
          clue_tau.features()(tid, 9) = (pdgid_abs == 22) ? 1.0f : 0.0f;

          // pad mask
          clue_tau.pad_mask()(tid) = 1.0f;
        }
      }
    }
  };

  SoftTauInputDeviceTensor transform(Queue& queue, 
                 const PFCandidateDeviceCollection& pf,
                 const CandsClusterBxDeviceCollection& candsClusterBx) {
    const auto kNumClusters = candsClusterBx.const_view<ClusterOffsetsSoA>().metadata().size() - 1;
    auto input_tensor = SoftTauInputDeviceTensor(kNumClusters, queue);
    input_tensor.zeroInitialise(queue);

    std::cout << "INITIALIZED SoftTauInputDeviceTensor" << std::endl;

    alpaka::exec<Acc1D>(queue, 
      make_workdiv<Acc1D>(kNumClusters, 128), 
      ComputeClueTauFeaturesKernel{}, 
      pf.const_view(),
      candsClusterBx.view<LongIndexSoA>(),
      candsClusterBx.view<ClusterOffsetsSoA>(),
      input_tensor.view());

    return input_tensor;
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels
