#ifndef DataFormats_L1ScoutingSoA_interface_LongAssociationMapHost_h
#define DataFormats_L1ScoutingSoA_interface_LongAssociationMapHost_h

#include <alpaka/alpaka.hpp>

#include "DataFormats/Portable/interface/PortableCollection.h"
#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/LongIndexSoA.h"
#include "DataFormats/L1ScoutingSoA/interface/LongOffsetsSoA.h"

namespace l1sc {

  using LongAssociationMapHost = PortableMultiCollection<alpaka::DevCpu, LongIndexSoA, LongOffsetsSoA>;

}  // namespace l1sc

#endif  // DataFormats_L1ScoutingSoA_interface_LongAssociationMapHost_h