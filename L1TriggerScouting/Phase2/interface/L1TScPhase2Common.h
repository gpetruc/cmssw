#ifndef L1TriggerScouting_Phase2_interface_L1TScPhase2Common_h
#define L1TriggerScouting_Phase2_interface_L1TScPhase2Common_h

#include <compare>
#include <cstdint>

namespace l1sc {

  using data_t = uint64_t;

  enum class Environment : int { kProduction = 0, kDevelopment = 1, kTest = 2, kDebug = 3 };

  constexpr std::strong_ordering operator<=>(Environment t, Environment u) {
    return static_cast<int>(t) <=> static_cast<int>(u);
  }

}  // namespace l1sc

#endif  // L1TriggerScouting_Phase2_interface_L1TScPhase2Common_h