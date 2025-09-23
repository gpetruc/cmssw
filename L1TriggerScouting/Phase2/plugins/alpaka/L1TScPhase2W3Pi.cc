#include "DataFormats/L1ScoutingSoA/interface/alpaka/DeviceCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/DeviceObject.h"
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

  private:
    const device::EDGetToken<PuppiDeviceCollection> puppi_in_token_;                            // PUPPI candidates
    const device::EDGetToken<NbxMapDeviceCollection> nbx_map_in_token_;                            // orbit association map
    const device::EDPutToken<PuppiDeviceCollection> puppi_out_token_;                           // PUPPI candidates
    const device::EDPutToken<NbxMapDeviceCollection> nbx_map_out_token_;
    const device::EDPutToken<W3PiPuppiTableDeviceCollection> w3pi_table_out_token_;
    const bool verbose_;                                                                        // verbose output
  };

  // __________________________________________________________________________________________________________________
  // IMPLEMENTATION

  L1TScPhase2W3Pi::L1TScPhase2W3Pi(const edm::ParameterSet &params)
      : EDProducer<>(params),
        puppi_in_token_{consumes(params.getParameter<edm::InputTag>("src"))},
        nbx_map_in_token_{consumes(params.getParameter<edm::InputTag>("src"))},
        puppi_out_token_{produces()},
        nbx_map_out_token_{produces()},
        w3pi_table_out_token_{produces()},
        verbose_(params.getUntrackedParameter<bool>("verbose")) {}

  void L1TScPhase2W3Pi::produce(device::Event &event, const device::EventSetup &event_setup) {
    auto timestamp = std::chrono::steady_clock::now();
    // TODO: remove const_cast operation and replace with separate buffer for association map
    auto &puppi = const_cast<PuppiDeviceCollection &>(event.get(puppi_in_token_));
    auto &nbx_map = const_cast<NbxMapDeviceCollection &>(event.get(nbx_map_in_token_));

    auto timestamp2 = std::chrono::steady_clock::now();
    kernels::runW3Pi(event.queue(), puppi, nbx_map);
    alpaka::wait(event.queue());
    auto elapsed2 = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - timestamp2);

    // TODO: analysis should not have explicit device sync to produce table, implement kernels for this
    std::array<int32_t, 2> const sizes{{nbx_map.view<NbxSoA>().metadata().size(), nbx_map.view<OffsetsSoA>().metadata().size()}};
    auto puppi_host = PuppiHostCollection(puppi.view().metadata().size(), event.queue());
    auto nbx_map_host = NbxMapHostCollection(sizes, event.queue());
    alpaka::memcpy(event.queue(), puppi_host.buffer(), puppi.buffer());
    alpaka::memcpy(event.queue(), nbx_map_host.buffer(), nbx_map.buffer());
    alpaka::wait(event.queue());  // explicitly synchronize the device

    auto timestamp3 = std::chrono::steady_clock::now();
    auto table_host = kernels::makeW3PiPuppiTable(event.queue(), puppi_host, nbx_map_host);
    auto table = W3PiPuppiTableDeviceCollection(table_host.view().metadata().size(), event.queue());
    alpaka::memcpy(event.queue(), table.buffer(), table_host.buffer());
    alpaka::wait(event.queue());
    auto elapsed3 = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - timestamp3);
    // END TODO

    // debug log
    // if (verbose_) {
    //   fmt::print("[DEBUG] l1sc::L1TScPhase2W3Pi: OK");
    // }

    event.emplace(puppi_out_token_, std::move(puppi));
    event.emplace(nbx_map_out_token_, std::move(nbx_map));
    event.emplace(w3pi_table_out_token_, std::move(table));

    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - timestamp);
    // if (verbose_) {
    fmt::print("[DEBUG] l1sc::L1TScPhase2W3Pi: OK {} us\n", elapsed.count());
    fmt::print("\tl1sc::runW3Pi: OK {} us\n", elapsed2.count());
    fmt::print("\tl1sc::makeW3PiPuppiTable: OK {} us\n", elapsed3.count());
    // }
  }

  void L1TScPhase2W3Pi::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src");
    desc.addUntracked<bool>("verbose", false);
    descriptions.addWithDefaultLabel(desc);
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc

DEFINE_FWK_ALPAKA_MODULE(l1sc::L1TScPhase2W3Pi);