#ifndef DataFormats_L1ScoutingSoA_interface_PuppiHostCollection_h
#define DataFormats_L1ScoutingSoA_interface_PuppiHostCollection_h

#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/PuppiSoA.h"

namespace l1sc {

  using PuppiHostCollection = PortableHostCollection<PuppiSoA>;

}  // namespace l1sc

#endif  // DataFormats_L1ScoutingSoA_interface_PuppiHostCollection_h