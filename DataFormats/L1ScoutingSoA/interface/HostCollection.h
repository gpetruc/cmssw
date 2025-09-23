#ifndef DataFormats_L1ScoutingSoA_interface_HostCollection_h
#define DataFormats_L1ScoutingSoA_interface_HostCollection_h

#include "DataFormats/Portable/interface/PortableCollection.h"
#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/SoA.h"

namespace l1sc {

  using PuppiHostCollection = PortableHostCollection<PuppiSoA>;
  using NbxMapHostCollection = PortableMultiCollection<alpaka::DevCpu, NbxSoA, OffsetsSoA>;
  using W3PiPuppiTableHostCollection = PortableHostCollection<W3PiPuppiTableSoA>;

}  // namespace l1sc

#endif  // DataFormats_L1ScoutingSoA_interface_HostCollection_h