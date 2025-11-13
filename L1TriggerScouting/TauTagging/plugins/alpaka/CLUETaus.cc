#include "DataFormats/L1ScoutingSoA/interface/alpaka/AssociationMapDevice.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/CandsClusterBxDeviceCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/ClustersDeviceCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/PFCandidateDeviceCollection.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/Event.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EventSetup.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/EDProducer.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "L1TriggerScouting/Phase2/interface/L1TScPhase2Common.h"
#include "L1TriggerScouting/TauTagging/plugins/alpaka/CLUEsteringAlgo.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc {

  class CLUETaus : public stream::EDProducer<> {
  public:
    explicit CLUETaus(const edm::ParameterSet &params)
        : EDProducer<>(params),
          pf_candidates_token_{consumes(params.getParameter<edm::InputTag>("src"))},
          bx_lookup_token_{consumes(params.getParameter<edm::InputTag>("src"))},
          cluestering_token_{produces()},
          candsclusterbx_token_{produces()},
          clustering_(static_cast<float>(params.getParameter<double>("dc")),
                      static_cast<float>(params.getParameter<double>("rhoc")),
                      static_cast<float>(params.getParameter<double>("dm")),
                      params.getParameter<bool>("wrapCoords")) {}

    void produce(device::Event &event, const device::EventSetup &event_setup) override {
      // get collection from device memory space (implicit copy done by framework)
      const auto &pf = event.get(pf_candidates_token_);
      const auto n_points = pf.const_view().metadata().size();

      // allocate buffer for the index of the cluster for each pf candidate
      auto clusters = ClustersDeviceCollection(n_points, event.queue());

      // run CLUEstering algo
      const auto &bx_lookup = event.get(bx_lookup_token_);
      auto candsclusterbx_map = clustering_.run(event.queue(), pf, bx_lookup, clusters);
      event.emplace(candsclusterbx_token_, std::move(candsclusterbx_map));

      // move clustering results to event storage
      event.emplace(cluestering_token_, std::move(clusters));
    }

    static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<edm::InputTag>("src");
      desc.add<double>("dc");
      desc.add<double>("rhoc");
      desc.add<double>("dm");
      desc.add<bool>("wrapCoords");
      descriptions.addWithDefaultLabel(desc);
    }

  private:
    // get device pf data
    const device::EDGetToken<PFCandidateDeviceCollection> pf_candidates_token_;
    // get association map if runScouting=True
    const device::EDGetToken<BxLookupDeviceCollection> bx_lookup_token_;
    // put device clustering data
    const device::EDPutToken<ClustersDeviceCollection> cluestering_token_;
    const device::EDPutToken<CandsClusterBxDeviceCollection> candsclusterbx_token_;
    // algorithm
    const kernels::CLUEsteringAlgo clustering_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc

DEFINE_FWK_ALPAKA_MODULE(l1sc::CLUETaus);