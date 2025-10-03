#ifndef DataFormats_L1ScoutingSoA_interface_ClustersHostCollection_h
#define DataFormats_L1ScoutingSoA_interface_ClustersHostCollection_h

#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/ClustersSoA.h"

namespace l1sc {

  using ClustersHostCollection = PortableHostCollection<ClustersSoA>;

}  // namespace l1sc

#endif  // DataFormats_L1ScoutingSoA_interface_ClustersHostCollection_h