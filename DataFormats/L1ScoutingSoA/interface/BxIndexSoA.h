#ifndef DataFormats_L1ScoutingSoA_interface_BxIndexSoA_h
#define DataFormats_L1ScoutingSoA_interface_BxIndexSoA_h

#include "DataFormats/SoATemplate/interface/SoALayout.h"

namespace l1sc {

  GENERATE_SOA_LAYOUT(BxIndexLayout, SOA_COLUMN(uint32_t, bx))

  using BxIndexSoA = BxIndexLayout<>;

}  // namespace l1sc

#endif  // DataFormats_L1ScoutingSoA_interface_BxIndexSoA_h