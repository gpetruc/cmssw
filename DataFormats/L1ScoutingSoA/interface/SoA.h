#ifndef DataFormats_L1ScoutingSoA_interface_SoA_h
#define DataFormats_L1ScoutingSoA_interface_SoA_h

#include "DataFormats/SoATemplate/interface/SoACommon.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"
#include "DataFormats/SoATemplate/interface/SoAView.h"

namespace l1sc {

  // TODO: selection has to be moved to separate memory blob
  GENERATE_SOA_LAYOUT(PuppiLayout,
                      SOA_COLUMN(float, pt),
                      SOA_COLUMN(float, eta),
                      SOA_COLUMN(float, phi),
                      SOA_COLUMN(float, z0),
                      SOA_COLUMN(float, dxy),
                      SOA_COLUMN(float, puppiw),
                      SOA_COLUMN(uint8_t, quality),
                      SOA_COLUMN(int16_t, pdgid),
                      SOA_COLUMN(uint64_t, selection));

  using PuppiSoA = PuppiLayout<>;
  using PuppiSoAView = PuppiSoA::View;
  using PuppiSoAConstView = PuppiSoA::ConstView;

  GENERATE_SOA_LAYOUT(NbxLayout, SOA_COLUMN(uint32_t, bx), SOA_COLUMN(uint32_t, selected))

  using NbxSoA = NbxLayout<>;
  using NbxSoAView = NbxSoA::View;
  using NbxSoAConstView = NbxSoA::ConstView;

  GENERATE_SOA_LAYOUT(OffsetsLayout, SOA_COLUMN(uint32_t, offsets))

  using OffsetsSoA = OffsetsLayout<>;
  using OffsetsSoAView = OffsetsSoA::View;
  using OffsetsSoAConstView = OffsetsSoA::ConstView;

  GENERATE_SOA_LAYOUT(W3PiPuppiTableLayout, SOA_COLUMN(uint32_t, i), SOA_COLUMN(uint32_t, j), SOA_COLUMN(uint32_t, k))

  using W3PiPuppiTableSoA = W3PiPuppiTableLayout<>;
  using W3PiPuppiTableSoAView = W3PiPuppiTableSoA::View;
  using W3PiPuppiTableSoAConstView = W3PiPuppiTableSoA::ConstView;

}  // namespace l1sc

#endif  // DataFormats_L1ScoutingSoA_interface_SoA_h