#ifndef DataFormats_L1ScoutingSoA_interface_IndexSoA_h
#define DataFormats_L1ScoutingSoA_interface_IndexSoA_h

#include "DataFormats/SoATemplate/interface/SoALayout.h"

namespace l1sc {

  GENERATE_SOA_LAYOUT(IndexLayout, SOA_COLUMN(uint32_t, indexes))

  using IndexSoA = IndexLayout<>;

}  // namespace l1sc

#endif  // DataFormats_L1ScoutingSoA_interface_IndexSoA_h