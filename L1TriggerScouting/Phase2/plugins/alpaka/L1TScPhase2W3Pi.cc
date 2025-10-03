#include "DataFormats/L1ScoutingSoA/interface/alpaka/BxLookupDeviceCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/PuppiDeviceCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/SelectedBxDeviceCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/W3PiDeviceTable.h"
#include "DataFormats/Portable/interface/PortableObject.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "HeterogeneousCore/AlpakaCore/interface/MoveToDeviceCache.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDGetToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/Event.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EventSetup.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/EDProducer.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "L1TriggerScouting/Phase2/interface/L1TScPhase2Common.h"
#include "L1TriggerScouting/Phase2/interface/W3PiAlgoParams.h"
#include "L1TriggerScouting/Phase2/interface/alpaka/SynchronizingTimer.h"
#include "L1TriggerScouting/Phase2/interface/alpaka/L1TScPhase2W3PiKernels.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc {

  using namespace ::l1sc;

  // Heterogeneous implementation of cut based W -> 3pi algorithm.
  class L1TScPhase2W3Pi : public stream::EDProducer<> {
  public:
    L1TScPhase2W3Pi(const edm::ParameterSet &params)
        : EDProducer<>(params),
          puppi_token_{consumes(params.getParameter<edm::InputTag>("src"))},
          bx_lookup_token_{consumes(params.getParameter<edm::InputTag>("src"))},
          selected_bx_token_{produces()},
          w3pi_table_token_{produces()},
          environment_{static_cast<Environment>(params.getUntrackedParameter<int>("environment"))},
          fast_path_{params.getParameter<bool>("fast_path")},
          sync_timer_(std::in_place, "L1TScPhase2W3Pi", environment_),
          w3pi_params_{PortableHostObject<W3PiAlgoParams>{
              cms::alpakatools::host(),
              W3PiAlgoParams{
                  .pT_min = static_cast<uint8_t>(params.getParameter<double>("pT_min")),
                  .pT_int = static_cast<uint8_t>(params.getParameter<double>("pT_int")),
                  .pT_max = static_cast<uint8_t>(params.getParameter<double>("pT_max")),
                  .invariant_mass_lower_bound =
                      static_cast<float>(params.getParameter<double>("invariant_mass_lower_bound")),
                  .invariant_mass_upper_bound =
                      static_cast<float>(params.getParameter<double>("invariant_mass_upper_bound")),
                  .min_deltar_threshold = static_cast<float>(params.getParameter<double>("min_deltar_threshold")),
                  .max_deltar_threshold = static_cast<float>(params.getParameter<double>("max_deltar_threshold")),
                  .max_isolation_threshold = static_cast<float>(params.getParameter<double>("max_isolation_threshold")),
                  .ang_sep_lower_bound = static_cast<float>(params.getParameter<double>("ang_sep_lower_bound"))}}} {}

    void produce(device::Event &event, const device::EventSetup &event_setup) override {
      // debug / test
      sync_timer_.value().start(event.queue());

      // query input data (consumed on device-side)
      const auto &puppi = event.get(puppi_token_);
      const auto &bx_lookup = event.get(bx_lookup_token_);

      // get control params from device cache
      const W3PiAlgoParams *w3pi_params = w3pi_params_.get(event.queue()).data();

      // run w -> 3pi algorithm impl
      auto [selected_bxs, w3pi_table] = kernels::runW3Pi(event.queue(), puppi, bx_lookup, w3pi_params, fast_path_);

      // store device-side products
      event.emplace(selected_bx_token_, std::move(selected_bxs));
      event.emplace(w3pi_table_token_, std::move(w3pi_table));

      // debug / test end
      sync_timer_.value().sync(event.queue());
    }

    static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<edm::InputTag>("src");
      // params
      desc.add<double>("pT_min", 7);
      desc.add<double>("pT_int", 12);
      desc.add<double>("pT_max", 15);
      desc.add<double>("invariant_mass_lower_bound", 40.0);
      desc.add<double>("invariant_mass_upper_bound", 150.0);
      desc.add<double>("min_deltar_threshold", 0.01 * 0.01);
      desc.add<double>("max_deltar_threshold", 0.25 * 0.25);
      desc.add<double>("max_isolation_threshold", 2.0);
      desc.add<double>("ang_sep_lower_bound", 0.5 * 0.5);
      // fast mode
      desc.add<bool>("fast_path", false);
      // debug
      desc.addUntracked<int>("environment", static_cast<int>(Environment::kProduction));
      descriptions.addWithDefaultLabel(desc);
    }

  private:
    // consume device-side data
    const device::EDGetToken<PuppiDeviceCollection> puppi_token_;
    const device::EDGetToken<BxLookupDeviceCollection> bx_lookup_token_;

    // produce device-side data
    const device::EDPutToken<SelectedBxDeviceCollection> selected_bx_token_;
    const device::EDPutToken<W3PiDeviceTable> w3pi_table_token_;

    // utility members
    const Environment environment_;
    const bool fast_path_;

    // debug / test stats
    std::optional<SynchronizingTimer> sync_timer_;

    // control params
    cms::alpakatools::MoveToDeviceCache<Device, PortableHostObject<W3PiAlgoParams>> w3pi_params_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc

DEFINE_FWK_ALPAKA_MODULE(l1sc::L1TScPhase2W3Pi);