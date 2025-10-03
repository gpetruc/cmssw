#ifndef L1TriggerScouting_Phase2_plugins_alpaka_L1TScPhase2W3PiKernels_h
#define L1TriggerScouting_Phase2_plugins_alpaka_L1TScPhase2W3PiKernels_h

#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "DataFormats/Portable/interface/PortableHostObject.h"
#include "DataFormats/Portable/interface/alpaka/PortableObject.h"
#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/BxLookupDeviceCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/PuppiDeviceCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/SelectedBxDeviceCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/W3PiDeviceTable.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/CopyToHost.h"
#include "L1TriggerScouting/Phase2/interface/W3PiAlgoParams.h"

// These definitions are not stored in DataFormats/L1ScoutingSoA/
// since are designed to be used as helper types in the kernels
// and not meant to be stored in the FW event struct
namespace l1sc {

  struct PortableCounter {
    int value;
  };
  using CounterHost = PortableHostObject<PortableCounter>;

  GENERATE_SOA_LAYOUT(BufferLayout, SOA_COLUMN(uint32_t, value));
  using Buffer = BufferLayout<>;
  using BufferHost = PortableHostCollection<Buffer>;

  GENERATE_SOA_LAYOUT(SelectionBitmaskLayout, SOA_COLUMN(uint64_t, bits));
  using SelectionBitmaskSoA = SelectionBitmaskLayout<>;
  using SelectionBitmaskHostCollection = PortableHostCollection<SelectionBitmaskSoA>;

}  // namespace l1sc

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc {

  using namespace ::l1sc;
  using SelectionBitmaskDeviceCollection = PortableCollection<SelectionBitmaskSoA>;
  using BufferDevice = PortableCollection<Buffer>;
  using CounterDevice = PortableObject<PortableCounter>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels {

  std::tuple<SelectedBxDeviceCollection, W3PiDeviceTable> runW3Pi(Queue& queue,
                                                                  const PuppiDeviceCollection& puppi,
                                                                  const BxLookupDeviceCollection& bx_lookup,
                                                                  const W3PiAlgoParams* params,
                                                                  const bool fast_path = false);

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels

#endif  // L1TriggerScouting_Phase2_plugins_alpaka_L1TScPhase2W3PiKernels_h