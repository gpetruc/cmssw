#ifndef DataFormats_L1ScoutingSoA_interface_PFCandidateSoA_h
#define DataFormats_L1ScoutingSoA_interface_PFCandidateSoA_h

#include "DataFormats/SoATemplate/interface/SoALayout.h"

namespace l1sc {

  GENERATE_SOA_LAYOUT(PFCandidateLayout,
                      // real features
                      SOA_COLUMN(float, pt),
                      SOA_COLUMN(float, eta),
                      SOA_COLUMN(float, phi),
                      SOA_COLUMN(float, mass),
                      SOA_COLUMN(float, z0),
                      SOA_COLUMN(float, dxy),
                      SOA_COLUMN(float, puppiw),
                      SOA_COLUMN(int8_t, charge),
                      SOA_COLUMN(uint8_t, type),
                      SOA_COLUMN(int16_t, pdgid),
                      // hw features
                      SOA_COLUMN(uint16_t, hwPt),
                      SOA_COLUMN(int16_t, hwEta),
                      SOA_COLUMN(int16_t, hwPhi),
                      SOA_COLUMN(int16_t, hwZ0),
                      SOA_COLUMN(int16_t, hwDxy),
                      SOA_COLUMN(int8_t, hwQual),
                      SOA_COLUMN(int16_t, hwPuppiw));

  using PFCandidateSoA = PFCandidateLayout<>;
  using PFCandidateSoAView = PFCandidateSoA::View;
  using PFCandidateSoAConstView = PFCandidateSoA::ConstView;

}  // namespace l1sc

#endif  // DataFormats_L1ScoutingSoA_interface_PFCandidateSoA_h