#ifndef DataFormats_L1ScoutingSoA_interface_BxLookupHostCollection_h
#define DataFormats_L1ScoutingSoA_interface_BxLookupHostCollection_h

#include <alpaka/alpaka.hpp>

#include "DataFormats/Portable/interface/PortableCollection.h"
#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/BxIndexSoA.h"
#include "DataFormats/L1ScoutingSoA/interface/OffsetsSoA.h"

namespace l1sc {

  using BxLookupHostCollection = PortableMultiCollection<alpaka::DevCpu, BxIndexSoA, OffsetsSoA>;

}  // namespace l1sc

#endif  // DataFormats_L1ScoutingSoA_interface_BxLookupHostCollection_h