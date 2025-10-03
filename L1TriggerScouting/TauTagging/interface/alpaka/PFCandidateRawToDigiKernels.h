#ifndef L1TriggerScouting_TauTagging_interface_alpaka_PFCandidateRawToDigiKernels_h
#define L1TriggerScouting_TauTagging_interface_alpaka_PFCandidateRawToDigiKernels_h

#include <alpaka/alpaka.hpp>

#include "DataFormats/L1ScoutingSoA/interface/alpaka/BxLookupDeviceCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/PFCandidateDeviceCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "L1TriggerScouting/Phase2/interface/L1TScPhase2Common.h"
#include "L1TriggerScouting/Phase2/interface/alpaka/L1TScPhase2BitsEncoding.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels {

  // Device constant memory constructs.
  ALPAKA_STATIC_ACC_MEM_CONSTANT alpaka::DevGlobal<Acc1D, const int16_t[8]> kPdgid;
  ALPAKA_STATIC_ACC_MEM_CONSTANT alpaka::DevGlobal<Acc1D, const float[8]> kMass;
  ALPAKA_STATIC_ACC_MEM_CONSTANT alpaka::DevGlobal<Acc1D, const int8_t[8]> kCharge;
  ALPAKA_STATIC_ACC_MEM_CONSTANT alpaka::DevGlobal<Acc1D, const uint8_t[8]> kType;
  ALPAKA_STATIC_ACC_MEM_CONSTANT alpaka::DevGlobal<Acc1D, const float> kPi720;

  class PFCandidateRawToDigiKernels {
  public:
    PFCandidateRawToDigiKernels() = default;
    explicit PFCandidateRawToDigiKernels(Queue& queue);

    void initialize(Queue& queue);

  private:
    inline static std::once_flag init_flag_;
  };

  void decode(Queue& queue, data_t* p_data, PFCandidateDeviceCollection& pf);
  void decode(Queue& queue, data_t* h_data, BxLookupDeviceCollection& bx_lookup);

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels

#endif  // L1TriggerScouting_TauTagging_interface_alpaka_PFCandidateRawToDigiKernels_h