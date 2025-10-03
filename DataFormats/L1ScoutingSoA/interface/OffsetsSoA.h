#ifndef DataFormats_L1ScoutingSoA_interface_OffsetsSoA_h
#define DataFormats_L1ScoutingSoA_interface_OffsetsSoA_h

#include "DataFormats/SoATemplate/interface/SoALayout.h"

namespace l1sc {

  GENERATE_SOA_LAYOUT(OffsetsLayout, SOA_COLUMN(uint32_t, offsets))

  using OffsetsSoA = OffsetsLayout<>;

}  // namespace l1sc

#endif  // DataFormats_L1ScoutingSoA_interface_OffsetsSoA_h