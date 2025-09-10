#ifndef DataFormats_L1ScoutingSoA_interface_alpaka_CounterDevice_h
#define DataFormats_L1ScoutingSoA_interface_alpaka_CounterDevice_h

#include "DataFormats/Portable/interface/alpaka/PortableObject.h"
#include "DataFormats/Portable/interface/PortableHostObject.h"
#include "DataFormats/L1ScoutingSoA/interface/CounterHost.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/CopyToHost.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc {

  /**
   * make the names from the top-level `l1sc` namespace visible for unqualified lookup
   * inside the `ALPAKA_ACCELERATOR_NAMESPACE::l1sc` namespace
   */
  using namespace ::l1sc;

  using CounterDevice = PortableObject<unsigned int>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc

ASSERT_DEVICE_MATCHES_HOST_COLLECTION(l1sc::CounterDevice, l1sc::CounterHost);
#endif  // DataFormats_L1ScoutingSoA_interface_alpaka_CounterDevice_h