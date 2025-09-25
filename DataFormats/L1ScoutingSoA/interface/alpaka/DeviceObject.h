#ifndef DataFormats_PortableTestObjects_interface_alpaka_DeviceObject_h
#define DataFormats_PortableTestObjects_interface_alpaka_DeviceObject_h

#include "DataFormats/Portable/interface/alpaka/PortableObject.h"
#include "DataFormats/L1ScoutingSoA/interface/HostObject.h"
#include "DataFormats/L1ScoutingSoA/interface/Struct.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/CopyToHost.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc {

  /**
   * make the names from the top-level `l1sc` namespace visible for unqualified lookup
   * inside the `ALPAKA_ACCELERATOR_NAMESPACE::l1sc` namespace
   */
  using namespace ::l1sc;

  using W3PiTripletDeviceObject = PortableObject<W3PiTriplet>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc

ASSERT_DEVICE_MATCHES_HOST_COLLECTION(l1sc::W3PiTripletDeviceObject, l1sc::W3PiTripletHostObject);

#endif  // DataFormats_PortableTestObjects_interface_alpaka_DeviceObject_h