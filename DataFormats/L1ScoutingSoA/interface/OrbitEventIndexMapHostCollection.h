#ifndef DataFormats_L1ScoutingSoA_interface_OrbitEventIndexMapHostCollection_h
#define DataFormats_L1ScoutingSoA_interface_OrbitEventIndexMapHostCollection_h

#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/OrbitEventIndexMapSoA.h"

namespace l1sc {

  using OrbitEventIndexMapHostCollection = PortableHostCollection<OrbitEventIndexMapSoA>;

}  // namespace l1sc

#endif  // DataFormats_L1ScoutingSoA_interface_OrbitEventIndexMapHostCollection_h