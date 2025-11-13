#ifndef DataFormats_L1ScoutingSoA_interface_LongIndexSoA_h
#define DataFormats_L1ScoutingSoA_interface_LongIndexSoA_h

#include "DataFormats/SoATemplate/interface/SoALayout.h"

namespace l1sc {

  GENERATE_SOA_LAYOUT(LongIndexLayout, SOA_COLUMN(int32_t, indexes))

  using LongIndexSoA = LongIndexLayout<>;

}  // namespace l1sc

#endif  // DataFormats_L1ScoutingSoA_interface_IndexSoA_h