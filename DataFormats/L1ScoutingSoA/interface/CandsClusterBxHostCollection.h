#ifndef DataFormats_L1ScoutingSoA_interface_CandsClusterBxHostCollection_h
#define DataFormats_L1ScoutingSoA_interface_CandsClusterBxHostCollection_h

#include <alpaka/alpaka.hpp>

#include "DataFormats/Portable/interface/PortableCollection.h"
#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/LongIndexSoA.h"
#include "DataFormats/L1ScoutingSoA/interface/LongOffsetsSoA.h"
#include "DataFormats/L1ScoutingSoA/interface/ClusterOffsetsSoA.h"
#include "DataFormats/L1ScoutingSoA/interface/ClusterIndexSoA.h"

namespace l1sc {

  using CandsClusterBxHostCollection = PortableMultiCollection<alpaka::DevCpu, LongIndexSoA, ClusterOffsetsSoA, ClusterIndexSoA, LongOffsetsSoA>;

}  // namespace l1sc

#endif  // DataFormats_L1ScoutingSoA_interface_BxLookupHostCollection_h