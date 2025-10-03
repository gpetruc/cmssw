#ifndef DataFormats_L1ScoutingSoA_interface_alpaka_ClustersDeviceCollection_h
#define DataFormats_L1ScoutingSoA_interface_alpaka_ClustersDeviceCollection_h

#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/ClustersHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/ClustersSoA.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/CopyToHost.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc {

  using namespace ::l1sc;

  using ClustersDeviceCollection = PortableCollection<ClustersSoA>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc

ASSERT_DEVICE_MATCHES_HOST_COLLECTION(l1sc::ClustersDeviceCollection, l1sc::ClustersHostCollection);

#endif  // DataFormats_L1ScoutingSoA_interface_alpaka_ClustersDeviceCollection_h