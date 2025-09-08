#include "L1TriggerScouting/Phase2/interface/alpaka/L1TScPhase2W3PiKernels.h"

#include "alpaka/alpaka.hpp"
#include "HeterogeneousCore/AlpakaInterface/interface/host.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels {

  using namespace cms::alpakatools;

  ALPAKA_FN_HOST uint32_t threadsPerBlockUpperBound(uint32_t val) {
    if (val <= 0)
      return 1;
    return std::pow(2, std::ceil(std::log2(val)));
  }

  template <typename TAcc>
  ALPAKA_FN_ACC int8_t charge(TAcc const& acc, int16_t cls) {
    return alpaka::math::abs(acc, static_cast<int>(cls)) == 11 ? (cls > 0 ? -1 : +1) : (cls > 0 ? +1 : -1);
  }

  template <typename TAcc>
  ALPAKA_FN_ACC float energy(TAcc const& acc, float pt, float eta, float mass) {
    float pz = pt * alpaka::math::sinh(acc, eta);
    float p = alpaka::math::sqrt(acc, pt * pt + pz * pz);
    return alpaka::math::sqrt(acc, p * p + mass * mass);
  }

  template <typename TAcc>
  ALPAKA_FN_ACC float massInvariant(
      TAcc const& acc, PuppiDeviceCollection::ConstView data, uint32_t i, uint32_t j, uint32_t k) {
    auto px1 = data.pt()[i] * alpaka::math::cos(acc, data.phi()[i]);
    auto py1 = data.pt()[i] * alpaka::math::sin(acc, data.phi()[i]);
    auto pz1 = data.pt()[i] * alpaka::math::sinh(acc, data.eta()[i]);
    auto e1 = energy(acc, data.pt()[i], data.eta()[i], 0.1396);

    auto px2 = data.pt()[j] * alpaka::math::cos(acc, data.phi()[j]);
    auto py2 = data.pt()[j] * alpaka::math::sin(acc, data.phi()[j]);
    auto pz2 = data.pt()[j] * alpaka::math::sinh(acc, data.eta()[j]);
    auto e2 = energy(acc, data.pt()[j], data.eta()[j], 0.1396);

    auto px3 = data.pt()[k] * alpaka::math::cos(acc, data.phi()[k]);
    auto py3 = data.pt()[k] * alpaka::math::sin(acc, data.phi()[k]);
    auto pz3 = data.pt()[k] * alpaka::math::sinh(acc, data.eta()[k]);
    auto e3 = energy(acc, data.pt()[k], data.eta()[k], 0.1396);

    auto t_energy = e1 + e2 + e3;
    auto t_px = px1 + px2 + px3;
    auto t_py = py1 + py2 + py3;
    auto t_pz = pz1 + pz2 + pz3;

    auto t_momentum = t_px * t_px + t_py * t_py + t_pz * t_pz;
    auto invariant_mass = t_energy * t_energy - t_momentum;

    return invariant_mass > 0 ? alpaka::math::sqrt(acc, invariant_mass) : 0.0f;
  }

  template <typename TAcc>
  ALPAKA_FN_ACC float deltaPhi(TAcc const& acc, float phi1, float phi2) {
    const float M_PI_CONST = 3.14159265358979323846;
    const float M_2_PI_CONST = 2.0 * M_PI_CONST;
    auto r = alpaka::math::fmod(acc, phi2 - phi1, M_2_PI_CONST);
    if (r < -M_PI_CONST)
      return r + M_2_PI_CONST;
    if (r > M_PI_CONST)
      return r - M_2_PI_CONST;
    return r;
  }

  template <typename TAcc>
  ALPAKA_FN_ACC bool angularSeparation(TAcc const& acc,
                                       PuppiDeviceCollection::ConstView data,
                                       uint32_t pidx,
                                       uint32_t idx) {
    static constexpr float ang_sep_lower_bound = 0.5 * 0.5;
    float delta_eta = data.eta()[pidx] - data.eta()[idx];
    float delta_phi = deltaPhi(acc, data.phi()[pidx], data.phi()[idx]);
    float ang_sep = delta_eta * delta_eta + delta_phi * delta_phi;
    if (ang_sep < ang_sep_lower_bound)
      return false;
    return true;
  }

  template <typename TAcc>
  ALPAKA_FN_ACC bool coneIsolation(TAcc const& acc,
                                   PuppiDeviceCollection::ConstView data,
                                   uint32_t thread_idx,
                                   uint32_t span_begin,
                                   uint32_t span_end) {
    const float min_threshold = 0.01 * 0.01;
    const float max_threshold = 0.25 * 0.25;
    const float max_isolation_threshold = 2.0;

    float accumulated = 0.0f;
    for (auto idx = span_begin; idx < span_end; idx++) {
      if (thread_idx == idx)
        continue;
      auto delta_eta = data.eta()[thread_idx] - data.eta()[idx];
      auto delta_phi = deltaPhi(acc, data.phi()[thread_idx], data.phi()[idx]);

      float th_value = delta_eta * delta_eta + delta_phi * delta_phi;
      if (th_value >= min_threshold && th_value <= max_threshold) {
        accumulated += data.pt()[idx];
      }
    }
    return accumulated <= max_isolation_threshold * data.pt()[thread_idx];
  }

  // TODO: - move constants to accelerator device space constant memory
  //       - define traits for shared memory to be queried dynamically from hardware device
  //       - move cuts to portable struct and MoveToDeviceCache<> to persist in multi-gpu environments
  //       - is it possible to vectorize innermost loops to reduce register pressure
  class W3PiKernel {
  public:
    template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
    ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                  PuppiDeviceCollection::View data,
                                  OrbitEventIndexMapDeviceCollection::View map) const {
      const uint8_t SHARED_MEM_BLOCK = 128;
      const uint8_t min_threshold = 7;
      const uint8_t int_threshold = 12;
      const uint8_t high_threshold = 15;
      const float invariant_mass_upper_bound = 100.0;
      const float invariant_mass_lower_bound = 60.0;

      auto& size = alpaka::declareSharedVar<int, __COUNTER__>(acc);
      auto& mask = alpaka::declareSharedVar<int[SHARED_MEM_BLOCK], __COUNTER__>(acc);
      auto& best_score = alpaka::declareSharedVar<float, __COUNTER__>(acc);

      if (once_per_block(acc)) {
        size = 0;
        best_score = 0.0f;
        for (auto idx = 0; idx < SHARED_MEM_BLOCK; idx++)
          mask[idx] = 0;
      }

      uint32_t grid_dim = alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0];
      for (uint32_t block_idx : independent_groups(acc, grid_dim)) {
        uint32_t begin = map.offsets()[block_idx];
        uint32_t end = map.offsets()[block_idx + 1];
        if (end == 0xFFFFFFFF)
          continue;
        if (end - begin == 0)
          continue;
        uint32_t block_dim = end - begin;
        for (uint32_t tid : independent_group_elements(acc, block_dim)) {
          auto thread_idx = tid + begin;  // global index
          auto cls = alpaka::math::abs(acc, static_cast<int>(data.pdgid()[thread_idx]));
          if (cls != 211 && cls != 11)
            continue;
          auto pt = data.pt()[thread_idx];
          if (pt < min_threshold)
            continue;
          alpaka::atomicAdd(acc, &mask[tid], 1);
        }

        alpaka::syncBlockThreads(acc);

        if (once_per_block(acc)) {
          for (auto idx = 0; idx < SHARED_MEM_BLOCK; idx++) {
            if (mask[idx] == 1)
              alpaka::atomicAdd(acc, &size, 1);
          }
          if (size < 3)
            continue;
        }

        for (uint32_t tid : independent_group_elements(acc, block_dim)) {
          auto thread_idx = tid + begin;  // global index
          if (mask[tid] == 0)
            continue;
          if (data.pt()[thread_idx] < high_threshold)
            continue;
          if (!coneIsolation(acc, data, thread_idx, begin, end))
            continue;
          for (uint32_t i = 0; i < block_dim; i++) {
            auto global_i_idx = i + begin;
            if (mask[i] == 0)
              continue;
            if (global_i_idx == thread_idx || data.pt()[global_i_idx] < int_threshold)
              continue;
            if (data.pt()[global_i_idx] > data.pt()[thread_idx] ||
                (data.pt()[global_i_idx] == data.pt()[thread_idx] && global_i_idx < thread_idx))
              continue;
            if (!angularSeparation(acc, data, thread_idx, global_i_idx))
              continue;
            for (uint32_t j = 0; j < block_dim; j++) {
              auto global_j_idx = j + begin;
              if (mask[j] == 0)
                continue;
              if (global_j_idx == thread_idx || global_j_idx == global_i_idx || data.pt()[global_i_idx] < min_threshold)
                continue;
              if (data.pt()[global_j_idx] > data.pt()[thread_idx] ||
                  (data.pt()[j] == data.pt()[thread_idx] && global_j_idx < thread_idx))
                continue;
              if (data.pt()[global_j_idx] > data.pt()[global_i_idx] ||
                  (data.pt()[global_j_idx] == data.pt()[global_i_idx] && global_j_idx < global_i_idx))
                continue;
              if (alpaka::math::abs(
                      acc,
                      static_cast<int>(charge(acc, data.pdgid()[thread_idx]) + charge(acc, data.pdgid()[global_i_idx]) +
                                       charge(acc, data.pdgid()[global_j_idx]))) != 1)
                continue;
              auto mass = massInvariant(acc, data, thread_idx, global_i_idx, global_j_idx);
              if (mass < invariant_mass_lower_bound || mass > invariant_mass_upper_bound)
                continue;
              if (angularSeparation(acc, data, thread_idx, global_j_idx) &&
                  angularSeparation(acc, data, global_i_idx, global_j_idx)) {
                if (coneIsolation(acc, data, global_i_idx, begin, end) &&
                    coneIsolation(acc, data, global_j_idx, begin, end)) {
                  alpaka::atomicAdd(acc, &data.selection()[thread_idx], static_cast<uint32_t>(1));
                  alpaka::atomicAdd(acc, &data.selection()[global_i_idx], static_cast<uint32_t>(1));
                  alpaka::atomicAdd(acc, &data.selection()[global_j_idx], static_cast<uint32_t>(1));
                }
              }
            }
          }
        }
      }
    }
  };

  void runW3Pi(Queue& queue, PuppiDeviceCollection& puppi, OrbitEventIndexMapDeviceCollection& orbit_association_map) {
    uint32_t threads_per_block = threadsPerBlockUpperBound(128);  // particles per single processing block on hardware.
    uint32_t blocks_per_grid = orbit_association_map.view().metadata().size() - 1;
    auto grid = make_workdiv<Acc1D>(blocks_per_grid, threads_per_block);

    alpaka::exec<Acc1D>(queue, grid, W3PiKernel{}, puppi.view(), orbit_association_map.view());
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels