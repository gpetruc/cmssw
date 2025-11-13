#ifndef DataFormats_L1ScoutingSoA_interface_ClusterOffsetsSoA_h
#define DataFormats_L1ScoutingSoA_interface_ClusterOffsetsSoA_h

#include "DataFormats/SoATemplate/interface/SoALayout.h"

namespace l1sc {

  GENERATE_SOA_LAYOUT(ClusterOffsetsLayout, SOA_COLUMN(int32_t, offsets))

  using ClusterOffsetsSoA = ClusterOffsetsLayout<>;

}  // namespace l1sc

#endif  // DataFormats_L1ScoutingSoA_interface_OffsetsSoA_h