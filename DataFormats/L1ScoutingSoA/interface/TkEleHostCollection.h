#ifndef DataFormats_L1ScoutingSoA_interface_TkEleHostCollection_h
#define DataFormats_L1ScoutingSoA_interface_TkEleHostCollection_h

#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/TkEleSoA.h"

namespace l1sc {

  using TkEleHostCollection = PortableHostCollection<TkEleSoA>;

}  // namespace l1sc

#endif  // DataFormats_L1ScoutingSoA_interface_TkEleHostCollection_h