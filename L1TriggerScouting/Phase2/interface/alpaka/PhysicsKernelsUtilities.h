#ifndef L1TriggerScouting_Phase2_interface_alpaka_PhysicsKernelsUtilities_h
#define L1TriggerScouting_Phase2_interface_alpaka_PhysicsKernelsUtilities_h

#include <alpaka/alpaka.hpp>

#include "DataFormats/L1ScoutingSoA/interface/alpaka/PuppiDeviceCollection.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels {

  ALPAKA_FN_HOST inline uint32_t make_threads_per_block(uint32_t val) {
    if (val <= 0)
      return 1;
    return std::pow(2, std::ceil(std::log2(val)));
  }

  template <typename TAcc>
  ALPAKA_FN_ACC inline void printBits(TAcc const& acc, uint64_t bits) {
    for (int i = std::numeric_limits<uint64_t>::digits - 1; i >= 0; --i) {
      printf("%llu", (bits >> i) & 1ULL);
      if (i % 8 == 0)
        printf(" ");
    }
    printf("\n");
  }

  template <typename TAcc>
  ALPAKA_FN_ACC inline int8_t charge(TAcc const& acc, int16_t cls) {
    return alpaka::math::abs(acc, static_cast<int>(cls)) == 11 ? (cls > 0 ? -1 : +1) : (cls > 0 ? +1 : -1);
  }

  template <typename TAcc>
  ALPAKA_FN_ACC inline float energy(TAcc const& acc, float pt, float eta, float mass) {
    float pz = pt * alpaka::math::sinh(acc, eta);
    float p = alpaka::math::sqrt(acc, pt * pt + pz * pz);
    return alpaka::math::sqrt(acc, p * p + mass * mass);
  }

  template <typename TAcc>
  ALPAKA_FN_ACC inline float deltaPhi(TAcc const& acc, float phi1, float phi2) {
    const float M_2_PI_CONST = 2.0 * alpaka::math::constants::pi;
    auto r = alpaka::math::fmod(acc, phi2 - phi1, M_2_PI_CONST);
    if (r < -alpaka::math::constants::pi)
      return r + M_2_PI_CONST;
    if (r > alpaka::math::constants::pi)
      return r - M_2_PI_CONST;
    return r;
  }

  template <typename TAcc>
  ALPAKA_FN_ACC inline float massInvariant(
      TAcc const& acc, PuppiDeviceCollection::ConstView puppi, uint32_t i, uint32_t j, uint32_t k) {
    auto px1 = puppi.pt()[i] * alpaka::math::cos(acc, puppi.phi()[i]);
    auto py1 = puppi.pt()[i] * alpaka::math::sin(acc, puppi.phi()[i]);
    auto pz1 = puppi.pt()[i] * alpaka::math::sinh(acc, puppi.eta()[i]);
    auto e1 = energy(acc, puppi.pt()[i], puppi.eta()[i], 0.1396);

    auto px2 = puppi.pt()[j] * alpaka::math::cos(acc, puppi.phi()[j]);
    auto py2 = puppi.pt()[j] * alpaka::math::sin(acc, puppi.phi()[j]);
    auto pz2 = puppi.pt()[j] * alpaka::math::sinh(acc, puppi.eta()[j]);
    auto e2 = energy(acc, puppi.pt()[j], puppi.eta()[j], 0.1396);

    auto px3 = puppi.pt()[k] * alpaka::math::cos(acc, puppi.phi()[k]);
    auto py3 = puppi.pt()[k] * alpaka::math::sin(acc, puppi.phi()[k]);
    auto pz3 = puppi.pt()[k] * alpaka::math::sinh(acc, puppi.eta()[k]);
    auto e3 = energy(acc, puppi.pt()[k], puppi.eta()[k], 0.1396);

    auto t_energy = e1 + e2 + e3;
    auto t_px = px1 + px2 + px3;
    auto t_py = py1 + py2 + py3;
    auto t_pz = pz1 + pz2 + pz3;

    auto t_momentum = t_px * t_px + t_py * t_py + t_pz * t_pz;
    auto invariant_mass = t_energy * t_energy - t_momentum;

    return invariant_mass > 0 ? alpaka::math::sqrt(acc, invariant_mass) : 0.0f;
  }

  template <typename TAcc>
  ALPAKA_FN_ACC inline bool angularSeparation(TAcc const& acc,
                                              PuppiDeviceCollection::ConstView puppi,
                                              uint32_t pidx,
                                              uint32_t idx,
                                              const float ang_sep_lower_bound) {
    float delta_eta = puppi.eta()[pidx] - puppi.eta()[idx];
    float delta_phi = deltaPhi(acc, puppi.phi()[pidx], puppi.phi()[idx]);
    float ang_sep = delta_eta * delta_eta + delta_phi * delta_phi;
    if (ang_sep < ang_sep_lower_bound)
      return false;
    return true;
  }

  template <typename TAcc>
  ALPAKA_FN_ACC inline bool coneIsolation(TAcc const& acc,
                                          PuppiDeviceCollection::ConstView puppi,
                                          uint32_t thread_idx,
                                          uint32_t span_begin,
                                          uint32_t span_end,
                                          const float min_threshold,
                                          const float max_threshold,
                                          const float max_isolation_threshold) {
    float accumulated = 0.0f;
    for (auto idx = span_begin; idx < span_end; idx++) {
      if (thread_idx == idx)
        continue;
      auto delta_eta = puppi.eta()[thread_idx] - puppi.eta()[idx];
      auto delta_phi = deltaPhi(acc, puppi.phi()[thread_idx], puppi.phi()[idx]);

      float th_value = delta_eta * delta_eta + delta_phi * delta_phi;
      if (th_value >= min_threshold && th_value <= max_threshold) {
        accumulated += puppi.pt()[idx];
      }
    }
    return accumulated <= max_isolation_threshold * puppi.pt()[thread_idx];
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels

#endif  // L1TriggerScouting_Phase2_interface_alpaka_PhysicsKernelsUtilities_h