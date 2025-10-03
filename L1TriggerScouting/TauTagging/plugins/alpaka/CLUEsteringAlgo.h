#ifndef L1TriggerScouting_TauTagging_plugins_alpaka_CLUEsteringAlgo_h
#define L1TriggerScouting_TauTagging_plugins_alpaka_CLUEsteringAlgo_h

#include <alpaka/alpaka.hpp>

#include "CLUEstering/CLUEstering.hpp"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/BxLookupDeviceCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/ClustersDeviceCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/PFCandidateDeviceCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels {

  using namespace ::l1sc;

  constexpr size_t kDims = 2;

  class CLUEsteringAlgo {
  public:
    explicit CLUEsteringAlgo(float dc, float rhoc, float dm, bool wrap_coords);

    void run(Queue& queue, const PFCandidateDeviceCollection& pf, ClustersDeviceCollection& clusters) const;
    void run(Queue& queue,
             const PFCandidateDeviceCollection& pf,
             const BxLookupDeviceCollection& bx_lookup,
             ClustersDeviceCollection& clusters) const;

  private:
    float dc_, rhoc_, dm_;
    bool wrap_coords_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels

#endif  // L1TriggerScouting_TauTagging_plugins_alpaka_CLUEsteringAlgo_h