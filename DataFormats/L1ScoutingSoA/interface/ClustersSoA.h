#ifndef DataFormats_L1ScoutingSoA_interface_ClustersSoA_h
#define DataFormats_L1ScoutingSoA_interface_ClustersSoA_h

#include "DataFormats/SoATemplate/interface/SoALayout.h"

namespace l1sc {

  GENERATE_SOA_LAYOUT(ClustersLayout, SOA_COLUMN(int, cluster), SOA_COLUMN(int, is_seed))

  GENERATE_SOA_LAYOUT(
      ClusterObjLayout, SOA_COLUMN(float, pt), SOA_COLUMN(float, eta), SOA_COLUMN(float, phi), SOA_COLUMN(int, cluster))

  using ClustersSoA = ClustersLayout<>;
  using ClusterObjSoA = ClusterObjLayout<>;

}  // namespace l1sc

#endif  // DataFormats_L1ScoutingSoA_interface_ClustersSoA_h