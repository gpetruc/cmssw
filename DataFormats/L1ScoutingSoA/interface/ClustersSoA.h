#ifndef DataFormats_L1ScoutingSoA_interface_ClustersSoA_h
#define DataFormats_L1ScoutingSoA_interface_ClustersSoA_h

#include "DataFormats/SoATemplate/interface/SoALayout.h"

namespace l1sc {

  GENERATE_SOA_LAYOUT(ClustersLayout, SOA_COLUMN(int, cluster), SOA_COLUMN(int, is_seed))

  using ClustersSoA = ClustersLayout<>;

}  // namespace l1sc

#endif  // DataFormats_L1ScoutingSoA_interface_ClustersSoA_h