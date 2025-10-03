#ifndef DataFormats_L1ScoutingSoA_interface_alpaka_PuppiDeviceCollection_h
#define DataFormats_L1ScoutingSoA_interface_alpaka_PuppiDeviceCollection_h

#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/PuppiHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/PuppiSoA.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/CopyToHost.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc {

  // make the names from the top-level `l1sc` namespace visible for unqualified lookup
  // inside the `ALPAKA_ACCELERATOR_NAMESPACE::l1sc` namespace
  using namespace ::l1sc;

  using PuppiDeviceCollection = PortableCollection<PuppiSoA>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc

ASSERT_DEVICE_MATCHES_HOST_COLLECTION(l1sc::PuppiDeviceCollection, l1sc::PuppiHostCollection);

#endif  // DataFormats_L1ScoutingSoA_interface_alpaka_PuppiDeviceCollection_h