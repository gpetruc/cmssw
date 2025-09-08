#ifndef L1TriggerScouting_Phase2_interface_alpaka_L1TScPhase2BitsEncoding_h
#define L1TriggerScouting_Phase2_interface_alpaka_L1TScPhase2BitsEncoding_h

#include <cstdint>
#include <alpaka/alpaka.hpp>

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels {

  template <typename T>
  ALPAKA_FN_ACC T decodeBits(uint64_t word, unsigned int start, unsigned int width) {
    static_assert(std::is_integral<T>::value, "extract_unsigned_bits expects integral types");
    uint64_t mask = (width < 64) ? ((1ULL << width) - 1) : ~0ULL;
    return static_cast<T>((word >> start) & mask);
  }

  template <typename T>
  ALPAKA_FN_ACC T decodeBitsSigned(uint64_t word, unsigned int start, unsigned int width) {
    static_assert(std::is_integral<T>::value && std::is_signed<T>::value,
                  "extract_signed_bits expects signed integral types");
    uint64_t raw = (word >> start) & ((1ULL << width) - 1);
    if (raw & (1ULL << (width - 1)))
      raw |= (~0ULL << width);  // manual sign extension
    return static_cast<T>(raw);
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels

#endif  // L1TriggerScouting_Phase2_interface_alpaka_L1TScPhase2BitsEncoding_h