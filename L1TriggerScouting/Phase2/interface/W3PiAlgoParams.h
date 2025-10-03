#ifndef L1TriggerScouting_Phase2_interface_W3PiAlgoParams_h
#define L1TriggerScouting_Phase2_interface_W3PiAlgoParams_h

#include <cstdint>

struct W3PiAlgoParams {
  // pT thresholds
  uint8_t pT_min;
  uint8_t pT_int;
  uint8_t pT_max;

  // mass bounds
  float invariant_mass_lower_bound;
  float invariant_mass_upper_bound;

  // isolation limits
  float min_deltar_threshold;
  float max_deltar_threshold;
  float max_isolation_threshold;

  // angular separation bound
  float ang_sep_lower_bound;
};

#endif  // L1TriggerScouting_Phase2_interface_W3PiAlgoParams_h