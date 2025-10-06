#ifndef DataFormats_L1ScoutingSoA_interface_TkEleSoA_h
#define DataFormats_L1ScoutingSoA_interface_TkEleSoA_h

#include "DataFormats/SoATemplate/interface/SoALayout.h"

namespace l1sc {

  GENERATE_SOA_LAYOUT(TkEleLayout,
                      SOA_COLUMN(float, pt),
                      SOA_COLUMN(float, eta),
                      SOA_COLUMN(float, phi),
                      SOA_COLUMN(uint8_t, quality),
                      SOA_COLUMN(float, isolation));

  using TkEleSoA = TkEleLayout<>;

}  // namespace l1sc

#endif  // DataFormats_L1ScoutingSoA_interface_TkEleSoA_h