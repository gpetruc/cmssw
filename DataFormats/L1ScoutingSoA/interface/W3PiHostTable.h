#ifndef DataFormats_L1ScoutingSoA_interface_W3PiHostTable_h
#define DataFormats_L1ScoutingSoA_interface_W3PiHostTable_h

#include <alpaka/alpaka.hpp>

#include "DataFormats/Portable/interface/PortableCollection.h"
#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/W3PiTable.h"

namespace l1sc {

  using W3PiHostTable = PortableHostCollection<W3PiTable>;

}  // namespace l1sc

#endif  // DataFormats_L1ScoutingSoA_interface_W3PiHostTable_h