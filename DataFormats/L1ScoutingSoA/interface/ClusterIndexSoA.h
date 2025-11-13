#ifndef DataFormats_L1ScoutingSoA_interface_ClusterIndexSoA_h
#define DataFormats_L1ScoutingSoA_interface_ClusterIndexSoA_h

#include "DataFormats/SoATemplate/interface/SoALayout.h"

namespace l1sc {

  GENERATE_SOA_LAYOUT(ClusterIndexLayout, SOA_COLUMN(int32_t, indexes))

  using ClusterIndexSoA = ClusterIndexLayout<>;

}  // namespace l1sc

#endif  // DataFormats_L1ScoutingSoA_interface_IndexSoA_h