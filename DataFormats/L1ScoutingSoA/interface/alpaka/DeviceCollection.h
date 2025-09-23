#ifndef DataFormats_L1ScoutingSoA_interface_alpaka_DeviceCollection_h
#define DataFormats_L1ScoutingSoA_interface_alpaka_DeviceCollection_h

#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/HostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/SoA.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/CopyToHost.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc {

  /**
   * make the names from the top-level `l1sc` namespace visible for unqualified lookup
   * inside the `ALPAKA_ACCELERATOR_NAMESPACE::l1sc` namespace
   */
  using namespace ::l1sc;

  using PuppiDeviceCollection = PortableCollection<PuppiSoA>;
  using NbxMapDeviceCollection = PortableCollection2<NbxSoA, OffsetsSoA>;
  using W3PiPuppiTableDeviceCollection = PortableCollection<W3PiPuppiTableSoA>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc

ASSERT_DEVICE_MATCHES_HOST_COLLECTION(l1sc::PuppiDeviceCollection, l1sc::PuppiHostCollection);
ASSERT_DEVICE_MATCHES_HOST_COLLECTION(l1sc::NbxMapDeviceCollection, l1sc::NbxMapHostCollection);
ASSERT_DEVICE_MATCHES_HOST_COLLECTION(l1sc::W3PiPuppiTableDeviceCollection, l1sc::W3PiPuppiTableHostCollection);

#endif  // DataFormats_L1ScoutingSoA_interface_alpaka_DeviceCollection_h