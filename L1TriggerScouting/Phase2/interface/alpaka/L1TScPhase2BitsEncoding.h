#ifndef L1TriggerScouting_Phase2_interface_alpaka_L1TScPhase2BitsEncoding_h
#define L1TriggerScouting_Phase2_interface_alpaka_L1TScPhase2BitsEncoding_h

#include <cstdint>
#include <alpaka/alpaka.hpp>

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels {

  template <typename T, unsigned int start, unsigned int width>
  ALPAKA_FN_ACC T decodeBits(uint64_t word) {
    static_assert(std::is_integral<T>::value, "extract_unsigned_bits expects integral types");
    constexpr uint64_t mask = (width < 64) ? ((1ULL << width) - 1) : ~0ULL;
    return static_cast<T>((word >> start) & mask);
  }

  template <typename T, unsigned int start, unsigned int width>
  ALPAKA_FN_ACC T decodeBitsSigned(uint64_t word) {
    static_assert(std::is_integral<T>::value && std::is_signed<T>::value,
                  "extract_signed_bits expects signed integral types");
    auto sdata = static_cast<int64_t>(word << (64 - width - start));
    return static_cast<T>(sdata >> (64 - width));  // arithmetic right shift
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels

#endif  // L1TriggerScouting_Phase2_interface_alpaka_L1TScPhase2BitsEncoding_h