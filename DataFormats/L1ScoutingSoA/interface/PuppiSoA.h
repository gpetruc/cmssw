#ifndef DataFormats_L1ScoutingSoA_interface_PuppiSoA_h
#define DataFormats_L1ScoutingSoA_interface_PuppiSoA_h

#include "DataFormats/SoATemplate/interface/SoALayout.h"

namespace l1sc {

  GENERATE_SOA_LAYOUT(PuppiLayout,
                      SOA_COLUMN(float, pt),
                      SOA_COLUMN(float, eta),
                      SOA_COLUMN(float, phi),
                      SOA_COLUMN(float, z0),
                      SOA_COLUMN(float, dxy),
                      SOA_COLUMN(float, puppiw),
                      SOA_COLUMN(uint8_t, quality),
                      SOA_COLUMN(int16_t, pdgid));

  using PuppiSoA = PuppiLayout<>;

}  // namespace l1sc

#endif  // DataFormats_L1ScoutingSoA_interface_PuppiSoA_h