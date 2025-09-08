#ifndef L1TriggerScouting_Phase2_plugins_alpaka_L1TScPhase2W3PiKernels_h
#define L1TriggerScouting_Phase2_plugins_alpaka_L1TScPhase2W3PiKernels_h

#include <alpaka/alpaka.hpp>
#include "DataFormats/L1ScoutingSoA/interface/alpaka/OrbitEventIndexMapDeviceCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/PuppiDeviceCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels {

  void runW3Pi(Queue& queue, PuppiDeviceCollection& puppi, OrbitEventIndexMapDeviceCollection& orbit_association_map);

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels

#endif  // L1TriggerScouting_Phase2_plugins_alpaka_L1TScPhase2W3PiKernels_h