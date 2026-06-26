#ifndef DataFormats_L1ScoutingSoA_interface_VertexSoA_h
#define DataFormats_L1ScoutingSoA_interface_VertexSoA_h

#include "DataFormats/SoATemplate/interface/SoALayout.h"

namespace l1sc {

  GENERATE_SOA_LAYOUT(VertexLayout,
                      SOA_COLUMN(bool, valid),
                      SOA_COLUMN(float, z0),
                      SOA_COLUMN(uint16_t, multIn),
                      SOA_COLUMN(float, sumPt),
                      SOA_COLUMN(uint8_t, quality),
                      SOA_COLUMN(uint16_t, multOut),
                      SOA_COLUMN(uint16_t, unassigned));

  using VertexSoA = VertexLayout<>;

}  // namespace l1sc

#endif  // DataFormats_L1ScoutingSoA_interface_VertexSoA_h