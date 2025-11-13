#ifndef DataFormats_L1ScoutingSoA_interface_LongOffsetsSoA_h
#define DataFormats_L1ScoutingSoA_interface_LongOffsetsSoA_h

#include "DataFormats/SoATemplate/interface/SoALayout.h"

namespace l1sc {

  GENERATE_SOA_LAYOUT(LongOffsetsLayout, SOA_COLUMN(int32_t, offsets))

  using LongOffsetsSoA = LongOffsetsLayout<>;

}  // namespace l1sc

#endif  // DataFormats_L1ScoutingSoA_interface_LongOffsetsSoA_h