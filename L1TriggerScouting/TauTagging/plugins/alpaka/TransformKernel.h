#ifndef L1TriggerScouting_TauTagging_plugins_alpaka_TransformKernel_h
#define L1TriggerScouting_TauTagging_plugins_alpaka_TransformKernel_h

#include <alpaka/alpaka.hpp>

#include <cstdio>

#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "DataFormats/Portable/interface/PortableHostObject.h"
#include "DataFormats/Portable/interface/alpaka/PortableObject.h"
#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/SoftTauDeviceTensor.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/AssociationMapDevice.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/BxLookupDeviceCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/ClustersDeviceCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/PFCandidateDeviceCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

// #define __DEBUG__

// These definitions are not stored in DataFormats/L1ScoutingSoA/
// since are designed to be used as helper types in the kernels
// and not meant to be stored in the FW event struct
namespace l1sc {

  struct PortableCounter {
    int value;
  };
  using CounterHost = PortableHostObject<PortableCounter>;

}  // namespace l1sc

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc {

  using namespace ::l1sc;
  using CounterDevice = PortableObject<PortableCounter>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels {

  using namespace ::l1sc;

  SoftTauInputDeviceTensor transform(Queue& queue, 
                 const PFCandidateDeviceCollection& pf, 
                 const AssociationMapDevice& clusterCandsMap);

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels

#endif  // L1TriggerScouting_TauTagging_plugins_alpaka_TransformKernel_h