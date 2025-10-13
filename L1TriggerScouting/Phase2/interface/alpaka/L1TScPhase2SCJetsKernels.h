#ifndef L1TriggerScouting_Phase2_plugins_alpaka_L1TScPhase2SCJetsKernels_h
#define L1TriggerScouting_Phase2_plugins_alpaka_L1TScPhase2SCJetsKernels_h

#include <alpaka/alpaka.hpp>

#include "DataFormats/L1ScoutingSoA/interface/alpaka/BxLookupDeviceCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/ClustersDeviceCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/PuppiDeviceCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels {

  using namespace ::l1sc;

  class L1TScPhase2SCJetsKernels {
  public:
    explicit L1TScPhase2SCJetsKernels();

    void run(Queue& queue,
             const PuppiDeviceCollection& src,
             const BxLookupDeviceCollection& bx_lookup,
             float R2,
             ClustersDeviceCollection& clusters,
             ClusterObjDeviceCollection& jets) const;

    void run(Queue& queue,
             const PuppiDeviceCollection& src,
             const BxLookupDeviceCollection& bx_lookup,
             float R2,
             unsigned int nJets,
             ClustersDeviceCollection& clusters,
             ClusterObjDeviceCollection& jets) const;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels

#endif  // L1TriggerScouting_Phase2_plugins_alpaka_L1TScPhase2SCJetsKernels_h