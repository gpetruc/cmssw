#ifndef DataFormats_L1ScoutingSoA_interface_TkEmHostCollection_h
#define DataFormats_L1ScoutingSoA_interface_TkEmHostCollection_h

#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/TkEmSoA.h"

namespace l1sc {

  using TkEmHostCollection = PortableHostCollection<TkEmSoA>;

}  // namespace l1sc

#endif  // DataFormats_L1ScoutingSoA_interface_TkEmHostCollection_h