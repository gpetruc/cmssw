#include "DataFormats/L1ScoutingSoA/interface/alpaka/BxLookupDeviceCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/PuppiDeviceCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/ClustersDeviceCollection.h"
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
#include "L1TriggerScouting/Phase2/interface/alpaka/SynchronizingTimer.h"
#include "L1TriggerScouting/Phase2/interface/alpaka/L1TScPhase2SCJetsKernels.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc {

  using namespace ::l1sc;

  class L1TScPhase2SCJets : public stream::EDProducer<> {
  public:
    L1TScPhase2SCJets(const edm::ParameterSet &params)
        : EDProducer<>(params),
          src_candidates_token_{consumes(params.getParameter<edm::InputTag>("src"))},
          bx_lookup_token_{consumes(params.getParameter<edm::InputTag>("src"))},
          clusters_token_{produces()},
          jets_token_{produces()},
          R2_{std::pow(params.getParameter<double>("rParam"), 2)},
          nJets_{params.getParameter<unsigned int>("nJets")} {}

    void produce(device::Event &event, const device::EventSetup &event_setup) override {
      // get raw data input
      const auto &src = event.get(src_candidates_token_);
      const auto &bx_lookup = event.get(bx_lookup_token_);

      // allocate buffer
      const auto nsrc = src.const_view().metadata().size();
      auto clusters = ClustersDeviceCollection(nsrc, event.queue());
      auto jets = ClusterObjDeviceCollection(nsrc, event.queue());

      // run
      if (nJets_ == 0) {
        kernels_.run(event.queue(), src, bx_lookup, R2_, clusters, jets);
      } else {
        kernels_.run(event.queue(), src, bx_lookup, R2_, nJets_, clusters, jets);
      }

      // move clustering results to event storage
      event.emplace(clusters_token_, std::move(clusters));
      event.emplace(jets_token_, std::move(jets));
    };

    static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<edm::InputTag>("src");
      desc.add<double>("rParam", 0.4);
      desc.add<unsigned int>("nJets", 0);
      descriptions.addWithDefaultLabel(desc);
    };

  private:
    // get device pf data
    const device::EDGetToken<PuppiDeviceCollection> src_candidates_token_;
    // get association map if runScouting=False
    const device::EDGetToken<BxLookupDeviceCollection> bx_lookup_token_;
    // put device clustering data
    const device::EDPutToken<ClustersDeviceCollection> clusters_token_;
    const device::EDPutToken<ClusterObjDeviceCollection> jets_token_;

    // kernel
    kernels::L1TScPhase2SCJetsKernels kernels_;

    // params
    double R2_;
    unsigned int nJets_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc

DEFINE_FWK_ALPAKA_MODULE(l1sc::L1TScPhase2SCJets);