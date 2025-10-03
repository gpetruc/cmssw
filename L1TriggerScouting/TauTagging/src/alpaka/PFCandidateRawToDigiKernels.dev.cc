#include "L1TriggerScouting/TauTagging/interface/alpaka/PFCandidateRawToDigiKernels.h"

#include "DataFormats/L1TParticleFlow/interface/PFCandidate.h"
#include "HeterogeneousCore/AlpakaInterface/interface/prefixScan.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels {

  PFCandidateRawToDigiKernels::PFCandidateRawToDigiKernels(Queue& queue) { initialize(queue); }

  // Initialize device constant memory for the kernels.
  // Called only once (thread-safe)
  void PFCandidateRawToDigiKernels::initialize(Queue& queue) {
    std::call_once(init_flag_, [&]() {
      // pdgid mapping
      constexpr int16_t host_pdgid[8] = {130, 22, -211, 211, 11, -11, 13, -13};
      auto view = createView(cms::alpakatools::host(), host_pdgid, Vec1D{8});
      alpaka::memcpy(queue, kPdgid<Acc1D>, view);

      // hw to float conversion: 3.14 / 720.0
      constexpr float host_pi_720 = alpaka::math::constants::pi / 720.0f;
      auto view_var = createView(cms::alpakatools::host(), &host_pi_720, Vec1D{1});
      alpaka::memcpy(queue, kPi720<Acc1D>, view_var);

      // mass mapping
      constexpr float host_masses[8] = {0.5, 0.0, 0.13, 0.13, 0.0005, 0.0005, 0.105, 0.105};
      auto view_mass = createView(cms::alpakatools::host(), host_masses, Vec1D{8});
      alpaka::memcpy(queue, kMass<Acc1D>, view_mass);

      // charge mapping
      // 0: neutral hadron, 0: photon, -1 pi-, +1 pi+, -1 e-, +1 e+, -1 mu-, +1 mu+
      constexpr int8_t host_charges[8] = {0, 0, -1, +1, -1, +1, -1, +1};
      auto view_charge = createView(cms::alpakatools::host(), host_charges, Vec1D{8});
      alpaka::memcpy(queue, kCharge<Acc1D>, view_charge);

      // type mapping
      constexpr uint8_t host_types[8] = {static_cast<uint8_t>(::l1t::PFCandidate::ParticleType::NeutralHadron),
                                         static_cast<uint8_t>(::l1t::PFCandidate::ParticleType::Photon),
                                         static_cast<uint8_t>(::l1t::PFCandidate::ParticleType::ChargedHadron),
                                         static_cast<uint8_t>(::l1t::PFCandidate::ParticleType::ChargedHadron),
                                         static_cast<uint8_t>(::l1t::PFCandidate::ParticleType::Electron),
                                         static_cast<uint8_t>(::l1t::PFCandidate::ParticleType::Electron),
                                         static_cast<uint8_t>(::l1t::PFCandidate::ParticleType::Muon),
                                         static_cast<uint8_t>(::l1t::PFCandidate::ParticleType::Muon)};
      auto view_type = createView(cms::alpakatools::host(), host_types, Vec1D{8});
      alpaka::memcpy(queue, kType<Acc1D>, view_type);
    });
  }

  class RawToDigiKernel {
  public:
    ALPAKA_FN_ACC void operator()(Acc1D const& acc, data_t* p_data, PFCandidateDeviceCollection::View pf) const {
      for (int32_t idx : cms::alpakatools::uniform_elements(acc, pf.metadata().size())) {
        uint64_t data = p_data[idx];

        // hardware values
        auto hwPt = decodeBits<uint16_t>(data, 0, 14);
        auto hwEta = decodeBitsSigned<int16_t>(data, 14, 12);
        auto hwPhi = decodeBitsSigned<int16_t>(data, 26, 11);
        auto pid = decodeBits<uint8_t>(data, 37, 3);

        // hw values
        pf.hwPt()[idx] = hwPt;
        pf.hwEta()[idx] = hwEta;
        pf.hwPhi()[idx] = hwPhi;

        // convert to real values
        pf.pt()[idx] = hwPt * 0.25f;
        pf.eta()[idx] = hwEta * kPi720<Acc1D>.get();
        pf.phi()[idx] = hwPhi * kPi720<Acc1D>.get();

        // pid mapping
        pf.mass()[idx] = kMass<Acc1D>.get()[pid];
        pf.puppiw()[idx] = 1.0f;
        pf.charge()[idx] = kCharge<Acc1D>.get()[pid];
        pf.type()[idx] = kType<Acc1D>.get()[pid];
        pf.pdgid()[idx] = kPdgid<Acc1D>.get()[pid];

        if (pid > 1) {
          auto hwZ0 = decodeBitsSigned<int16_t>(data, 40, 10);
          auto hwDxy = decodeBitsSigned<int16_t>(data, 50, 8);
          auto hwQual = decodeBits<uint8_t>(data, 58, 3);

          // hw values
          pf.hwZ0()[idx] = hwZ0;
          pf.hwDxy()[idx] = hwDxy;
          pf.hwPuppiw()[idx] = 1 << 8;  // default (-1 = 255)
          pf.hwQual()[idx] = hwQual;

          // convert to real values
          pf.z0()[idx] = hwZ0 * 0.05f;
          pf.dxy()[idx] = hwDxy * 0.05f;

        } else {
          auto hwPuppiw = decodeBits<uint16_t>(data, 40, 10);
          auto hwQual = decodeBits<uint8_t>(data, 50, 6);

          // hw values
          pf.hwZ0()[idx] = 0;
          pf.hwDxy()[idx] = 0;
          pf.hwPuppiw()[idx] = hwPuppiw;
          pf.hwQual()[idx] = hwQual;

          // convert to real values
          pf.z0()[idx] = 0.0f;
          pf.dxy()[idx] = 0.0f;
        }
      }
    }
  };

  void decode(Queue& queue, data_t* p_data, PFCandidateDeviceCollection& pf) {
    // move host residing data to device memory space
    auto extent = Vec1D{pf.const_view().metadata().size()};
    auto p_data_device = alpaka::allocAsyncBuf<data_t, Idx>(queue, extent);
    alpaka::memcpy(queue, p_data_device, createView(cms::alpakatools::host(), p_data, extent));

    // grid dims can be tuned for performance
    uint32_t threads_per_block = 1024;
    uint32_t blocks_per_grid = cms::alpakatools::divide_up_by(pf.const_view().metadata().size(), threads_per_block);
    auto grid = cms::alpakatools::make_workdiv<Acc1D>(blocks_per_grid, threads_per_block);

    // decode particles features
    alpaka::exec<Acc1D>(queue, grid, RawToDigiKernel{}, p_data_device.data(), pf.view());
  }

  void decode(Queue& queue, data_t* h_data, BxLookupDeviceCollection& bx_lookup) {
    // move host residing data to device memory space
    auto extent = Vec1D(bx_lookup.const_view().metadata().size());
    auto h_data_device = alpaka::allocAsyncBuf<data_t, Idx>(queue, extent);
    alpaka::memcpy(queue, h_data_device, createView(cms::alpakatools::host(), h_data, extent));

    // grid dims can be tuned for performance
    uint32_t threads_per_block = 1024;
    uint32_t blocks_per_grid =
        cms::alpakatools::divide_up_by(bx_lookup.const_view().metadata().size(), threads_per_block);
    auto grid = cms::alpakatools::make_workdiv<Acc1D>(blocks_per_grid, threads_per_block);

    // accumulate buffer with events sizes
    alpaka::exec<Acc1D>(
        queue,
        grid,
        [] ALPAKA_FN_ACC(Acc1D const& acc, data_t* data, BxIndexSoA::View bx_index, OffsetsSoA::View offsets) {
          if (cms::alpakatools::once_per_grid(acc))
            offsets.offsets()[0] = 0;

          for (int32_t idx : cms::alpakatools::uniform_elements(acc, offsets.metadata().size() - 1)) {
            auto range = decodeBits<uint32_t>(data[idx], 0, 12);
            offsets.offsets()[idx + 1] = range;
            bx_index.bx()[idx] = decodeBits<uint32_t>(data[idx], 12, 12);
          }
        },
        h_data_device.data(),
        bx_lookup.view<BxIndexSoA>(),
        bx_lookup.view<OffsetsSoA>());

    // prefix sum to build association map used for span extraction and batching
    auto pc = alpaka::allocAsyncBuf<int32_t, Idx>(queue, Vec1D{1});
    alpaka::memset(queue, pc, 0x00);
    alpaka::exec<Acc1D>(queue,
                        grid,
                        cms::alpakatools::multiBlockPrefixScan<uint32_t>{},
                        bx_lookup.view<OffsetsSoA>().offsets(),
                        bx_lookup.view<OffsetsSoA>().offsets(),
                        bx_lookup.view<OffsetsSoA>().metadata().size(),
                        blocks_per_grid,
                        pc.data(),
                        alpaka::getPreferredWarpSize(alpaka::getDev(queue)));
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels