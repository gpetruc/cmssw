#ifndef DataFormats_L1ScoutingSoA_interface_W3PiTable_h
#define DataFormats_L1ScoutingSoA_interface_W3PiTable_h

#include "DataFormats/SoATemplate/interface/SoALayout.h"

namespace l1sc {

  GENERATE_SOA_LAYOUT(W3PiTableLayout, SOA_COLUMN(uint32_t, i), SOA_COLUMN(uint32_t, j), SOA_COLUMN(uint32_t, k))

  using W3PiTable = W3PiTableLayout<>;

}  // namespace l1sc

#endif  // DataFormats_L1ScoutingSoA_interface_W3PiTable_h