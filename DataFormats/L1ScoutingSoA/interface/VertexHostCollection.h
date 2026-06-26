#ifndef DataFormats_L1ScoutingSoA_interface_VertexHostCollection_h
#define DataFormats_L1ScoutingSoA_interface_VertexHostCollection_h

#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/VertexSoA.h"

namespace l1sc {

  using VertexHostCollection = PortableHostCollection<VertexSoA>;

}  // namespace l1sc

#endif  // DataFormats_L1ScoutingSoA_interface_VertexHostCollection_h