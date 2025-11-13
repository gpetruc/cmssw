#ifndef DataFormats_L1ScoutingSoA_interface_alpaka_CandsClusterBxDeviceCollection_h
#define DataFormats_L1ScoutingSoA_interface_alpaka_CandsClusterBxDeviceCollection_h

#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/CandsClusterBxHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/LongIndexSoA.h"
#include "DataFormats/L1ScoutingSoA/interface/LongOffsetsSoA.h"
#include "DataFormats/L1ScoutingSoA/interface/ClusterOffsetsSoA.h"
#include "DataFormats/L1ScoutingSoA/interface/ClusterIndexSoA.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/CopyToHost.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc {

  // make the names from the top-level `l1sc` namespace visible for unqualified lookup
  // inside the `ALPAKA_ACCELERATOR_NAMESPACE::l1sc` namespace
  using namespace ::l1sc;

  using CandsClusterBxDeviceCollection = PortableCollection4<LongIndexSoA, ClusterOffsetsSoA, ClusterIndexSoA, LongOffsetsSoA>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc

ASSERT_DEVICE_MATCHES_HOST_COLLECTION(l1sc::CandsClusterBxDeviceCollection, l1sc::CandsClusterBxHostCollection);

#endif  // DataFormats_L1ScoutingSoA_interface_alpaka_BxLookupDeviceCollection_h