#ifndef L1TriggerScouting_Phase2_plugins_alpaka_L1TScPhase2VertexRawToDigiKernels_h
#define L1TriggerScouting_Phase2_plugins_alpaka_L1TScPhase2VertexRawToDigiKernels_h

#include <alpaka/alpaka.hpp>

#include "DataFormats/L1ScoutingSoA/interface/alpaka/BxLookupDeviceCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/VertexDeviceCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "L1TriggerScouting/Phase2/interface/L1TScPhase2Common.h"
#include "L1TriggerScouting/Phase2/plugins/alpaka/L1TScPhase2BitsEncoding.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels {

  class L1TScPhase2VertexRawToDigiKernels {
  public:
    L1TScPhase2VertexRawToDigiKernels() = default;
    explicit L1TScPhase2VertexRawToDigiKernels(Queue& queue);

    void initialize(Queue& queue);

  private:
    inline static std::once_flag init_flag_;
  };

  void decode(Queue& queue, data_t* p_data, VertexDeviceCollection& Vertex);

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels

#endif  // L1TriggerScouting_Phase2_plugins_alpaka_L1TScPhase2VertexRawToDigiKernels_h