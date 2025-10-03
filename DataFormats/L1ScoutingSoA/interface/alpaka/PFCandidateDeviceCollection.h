#ifndef DataFormats_L1ScoutingSoA_interface_alpaka_PFCandidateDeviceCollection_h
#define DataFormats_L1ScoutingSoA_interface_alpaka_PFCandidateDeviceCollection_h

#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/PFCandidateHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/PFCandidateSoA.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/CopyToHost.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc {

  using namespace ::l1sc;

  using PFCandidateDeviceCollection = PortableCollection<PFCandidateSoA>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc

ASSERT_DEVICE_MATCHES_HOST_COLLECTION(l1sc::PFCandidateDeviceCollection, l1sc::PFCandidateHostCollection);

#endif  // DataFormats_L1ScoutingSoA_interface_alpaka_PFCandidateDeviceCollection_h