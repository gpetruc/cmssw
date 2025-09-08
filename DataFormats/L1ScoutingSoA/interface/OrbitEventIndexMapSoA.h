#ifndef DataFormats_L1ScoutingSoA_interface_OrbitEventIndexMapSoA_h
#define DataFormats_L1ScoutingSoA_interface_OrbitEventIndexMapSoA_h

#include "DataFormats/Common/interface/StdArray.h"
#include "DataFormats/SoATemplate/interface/SoACommon.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"
#include "DataFormats/SoATemplate/interface/SoAView.h"

namespace l1sc {

  GENERATE_SOA_LAYOUT(OrbitEventIndexMapLayout, SOA_COLUMN(uint32_t, offsets))

  using OrbitEventIndexMapSoA = OrbitEventIndexMapLayout<>;
  using OrbitEventIndexMapSoAView = OrbitEventIndexMapSoA::View;
  using OrbitEventIndexMapSoAConstView = OrbitEventIndexMapSoA::ConstView;

}  // namespace l1sc

#endif  // DataFormats_L1ScoutingSoA_interface_OrbitEventIndexMapSoA_h