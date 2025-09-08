#ifndef DataFormats_L1ScoutingSoA_interface_alpaka_OrbitEventIndexMapDeviceCollection
#define DataFormats_L1ScoutingSoA_interface_alpaka_OrbitEventIndexMapDeviceCollection

#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/OrbitEventIndexMapHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/OrbitEventIndexMapSoA.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/CopyToHost.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc {

  /**
   * make the names from the top-level `l1sc` namespace visible for unqualified lookup
   * inside the `ALPAKA_ACCELERATOR_NAMESPACE::l1sc` namespace
   */
  using namespace ::l1sc;

  using OrbitEventIndexMapDeviceCollection = PortableCollection<OrbitEventIndexMapSoA>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc

ASSERT_DEVICE_MATCHES_HOST_COLLECTION(l1sc::OrbitEventIndexMapDeviceCollection, l1sc::OrbitEventIndexMapHostCollection);

#endif  // DataFormats_L1ScoutingSoA_interface_alpaka_OrbitEventIndexMapDeviceCollection