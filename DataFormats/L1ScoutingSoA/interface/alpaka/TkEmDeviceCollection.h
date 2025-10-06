#ifndef DataFormats_L1ScoutingSoA_interface_alpaka_TkEmDeviceCollection_h
#define DataFormats_L1ScoutingSoA_interface_alpaka_TkEmDeviceCollection_h

#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/TkEmHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/TkEmSoA.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/CopyToHost.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc {

  // make the names from the top-level `l1sc` namespace visible for unqualified lookup
  // inside the `ALPAKA_ACCELERATOR_NAMESPACE::l1sc` namespace
  using namespace ::l1sc;

  using TkEmDeviceCollection = PortableCollection<TkEmSoA>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc

ASSERT_DEVICE_MATCHES_HOST_COLLECTION(l1sc::TkEmDeviceCollection, l1sc::TkEmHostCollection);

#endif  // DataFormats_L1ScoutingSoA_interface_alpaka_TkEmDeviceCollection_h