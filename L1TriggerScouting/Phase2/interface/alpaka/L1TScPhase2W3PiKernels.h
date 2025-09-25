#ifndef L1TriggerScouting_Phase2_plugins_alpaka_L1TScPhase2W3PiKernels_h
#define L1TriggerScouting_Phase2_plugins_alpaka_L1TScPhase2W3PiKernels_h

#include <alpaka/alpaka.hpp>
#include "DataFormats/L1ScoutingSoA/interface/alpaka/DeviceCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels {

  void runW3Pi(Queue& queue, PuppiDeviceCollection& puppi, NbxMapDeviceCollection& nbx_map);
  W3PiPuppiTableHostCollection makeW3PiPuppiTable(Queue& queue, PuppiHostCollection& puppi_host, NbxMapHostCollection& nbx_map_host);

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels

#endif  // L1TriggerScouting_Phase2_plugins_alpaka_L1TScPhase2W3PiKernels_h