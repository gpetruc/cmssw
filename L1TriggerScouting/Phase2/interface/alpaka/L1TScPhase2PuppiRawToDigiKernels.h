#ifndef L1TriggerScouting_Phase2_interface_alpaka_L1TScPhase2PuppiRawToDigiKernels_h
#define L1TriggerScouting_Phase2_interface_alpaka_L1TScPhase2PuppiRawToDigiKernels_h

#include <alpaka/alpaka.hpp>

#include "DataFormats/L1ScoutingSoA/interface/alpaka/OrbitEventIndexMapDeviceCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/PuppiDeviceCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "L1TriggerScouting/Phase2/interface/L1TScPhase2Common.h"
#include "L1TriggerScouting/Phase2/interface/alpaka/L1TScPhase2BitsEncoding.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels {

  /**
   * Device constant memory constructs.
   */
  ALPAKA_STATIC_ACC_MEM_CONSTANT alpaka::DevGlobal<Acc1D, const int16_t[8]> kPdgid;
  ALPAKA_STATIC_ACC_MEM_CONSTANT alpaka::DevGlobal<Acc1D, const float> kPi720;

  class L1TScPhase2PuppiRawToDigiKernels {
  public:
    L1TScPhase2PuppiRawToDigiKernels() = default;
    explicit L1TScPhase2PuppiRawToDigiKernels(Queue& queue);

    void initialize(Queue& queue);

  private:
    static std::once_flag init_flag_;
  };

  void rawToDigi(Queue& queue, data_t* p_data, PuppiDeviceCollection& puppi);
  void associateOrbitEventIndex(Queue& queue,
                                data_t* h_data,
                                OrbitEventIndexMapDeviceCollection& orbit_association_map);

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels

#endif  // L1TriggerScouting_Phase2_interface_alpaka_L1TScPhase2PuppiRawToDigiKernels_h