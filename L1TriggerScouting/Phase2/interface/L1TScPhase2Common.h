#ifndef L1TriggerScouting_Phase2_interface_L1TScPhase2Common_h
#define L1TriggerScouting_Phase2_interface_L1TScPhase2Common_h

#include <cstdint>

namespace l1sc {

  using data_t = uint64_t;
  constexpr uint32_t kBitFieldSize =std::numeric_limits<data_t>::digits;
  constexpr uint32_t kOrbitSize = 3564;

}  // namespace l1sc

#endif  // L1TriggerScouting_Phase2_interface_L1TScPhase2Common_h