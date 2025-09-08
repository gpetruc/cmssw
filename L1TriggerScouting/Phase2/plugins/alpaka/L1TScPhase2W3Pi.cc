#include "DataFormats/L1ScoutingSoA/interface/alpaka/OrbitEventIndexMapDeviceCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/PuppiDeviceCollection.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDGetToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/Event.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EventSetup.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/EDProducer.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "L1TriggerScouting/Phase2/interface/alpaka/L1TScPhase2W3PiKernels.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc {

  using namespace ::l1sc;

  /**
   * @class L1TScPhase2W3Pi
   * @brief Heterogeneous implementation of cut based W -> 3pi algorithm. 
   */
  class L1TScPhase2W3Pi : public stream::EDProducer<> {
  public:
    L1TScPhase2W3Pi(const edm::ParameterSet &params);

    void produce(device::Event &event, const device::EventSetup &event_setup) override;
    static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);
    void beginStream(edm::StreamID) override;
    void endStream() override;
    void logDebugMessage(int event_id, const PuppiHostCollection &puppi, const OrbitEventIndexMapHostCollection &map);

  private:
    const device::EDGetToken<PuppiDeviceCollection> puppi_in_token_;                            // PUPPI candidates
    const device::EDGetToken<OrbitEventIndexMapDeviceCollection> orbit_association_map_token_;  // orbit association map
    const device::EDPutToken<PuppiDeviceCollection> puppi_out_token_;                           // PUPPI candidates
    const bool verbose_;                                                                        // verbose output
    const int verbose_level_;                                                                   // verbose level
    uint32_t w3pi_ct_ = 0;
    uint32_t w3pi_candidates_ = 0;
    std::chrono::steady_clock::time_point start_, end_;
    std::chrono::steady_clock::time_point ev_start_, ev_end_;
  };

  // __________________________________________________________________________________________________________________
  // IMPLEMENTATION

  L1TScPhase2W3Pi::L1TScPhase2W3Pi(const edm::ParameterSet &params)
      : EDProducer<>(params),
        puppi_in_token_{consumes(params.getParameter<edm::InputTag>("src"))},
        orbit_association_map_token_{consumes(params.getParameter<edm::InputTag>("src"))},
        puppi_out_token_{produces()},
        verbose_(params.getUntrackedParameter<bool>("verbose")),
        verbose_level_(params.getUntrackedParameter<int>("verboseLevel")) {}

  void L1TScPhase2W3Pi::produce(device::Event &event, const device::EventSetup &event_setup) {
    // TODO: remove const_cast operation and replace with separate buffer for association map
    auto &puppi = const_cast<PuppiDeviceCollection &>(event.get(puppi_in_token_));
    auto &orbit_association_map =
        const_cast<OrbitEventIndexMapDeviceCollection &>(event.get(orbit_association_map_token_));

    kernels::runW3Pi(event.queue(), puppi, orbit_association_map);

    // debug log to stdout
    if (verbose_) {
      w3pi_candidates_ += orbit_association_map.const_view().metadata().size() - 1;
      auto puppi_host = PuppiHostCollection(puppi.view().metadata().size(), event.queue());
      auto orbit_association_map_host =
          OrbitEventIndexMapHostCollection(orbit_association_map.view().metadata().size(), event.queue());
      alpaka::memcpy(event.queue(), puppi_host.buffer(), puppi.buffer());
      alpaka::memcpy(event.queue(), orbit_association_map_host.buffer(), orbit_association_map.buffer());
      alpaka::wait(event.queue());  // explicitly synchronize the device
      logDebugMessage(event.id().event(), puppi_host, orbit_association_map_host);
    }

    event.emplace(puppi_out_token_, std::move(puppi));
  }

  void L1TScPhase2W3Pi::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src");
    desc.addUntracked<bool>("verbose", false);
    desc.addUntracked<int>("verboseLevel", 0);
    descriptions.addWithDefaultLabel(desc);
  }

  void L1TScPhase2W3Pi::logDebugMessage(int event_id,
                                        const PuppiHostCollection &puppi,
                                        const OrbitEventIndexMapHostCollection &map) {
    auto cache_w3pi_ct = w3pi_ct_;
    for (int32_t idx = 0; idx < map.view().metadata().size() - 1; ++idx) {
      auto span_begin = map.view().offsets()[idx];
      auto span_end = map.view().offsets()[idx + 1];
      if (span_end - span_begin == 0)
        continue;
      int ct = 0;
      for (uint32_t i = span_begin; i < span_end; ++i) {
        ct += puppi.view().selection()[i];
      }
      if (ct > 0) {
        w3pi_ct_ = w3pi_ct_ + static_cast<uint32_t>(ct / 3);
      }
    }
    fmt::print("[DEBUG] l1sc::L1TScPhase2W3Pi: OK -> {} ({})\n", w3pi_ct_ - cache_w3pi_ct, event_id);
  }

  void L1TScPhase2W3Pi::beginStream(edm::StreamID) {
    if (verbose_) {
      fmt::print("============================================================\n");
      start_ = std::chrono::steady_clock::now();
    }
  }

  void L1TScPhase2W3Pi::endStream() {
    if (verbose_) {
      end_ = std::chrono::steady_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_ - start_);
      fmt::print("============================================================\n");
      fmt::print("[DEBUG] l1sc::L1TScPhase2W3Pi: {} -> {} ({} ms)\n", w3pi_candidates_, w3pi_ct_, duration.count());
      fmt::print("============================================================\n");
    }
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc

DEFINE_FWK_ALPAKA_MODULE(l1sc::L1TScPhase2W3Pi);