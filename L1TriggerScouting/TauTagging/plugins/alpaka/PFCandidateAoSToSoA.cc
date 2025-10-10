#include "DataFormats/L1ScoutingSoA/interface/PuppiHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/PuppiDeviceCollection.h"
#include "DataFormats/L1TParticleFlow/interface/PFCandidate.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/EDPutToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/Event.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EventSetup.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/EDProducer.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "L1TriggerScouting/Phase2/interface/L1TScPhase2Common.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc {

  using namespace ::l1sc;

  class PFCandidateAoSToSoA : public stream::EDProducer<> {
  public:
    PFCandidateAoSToSoA(const edm::ParameterSet &params)
        : EDProducer<>(params),
          pf_candidates_aos_token_{consumes(params.getParameter<edm::InputTag>("src"))},
          pf_candidates_soa_token_{produces()},
          environment_{static_cast<Environment>(params.getUntrackedParameter<int>("environment"))} {}

    void produce(device::Event &event, const device::EventSetup &event_setup) override {
      // grab PF candidates
      const auto &pf_candidates_aos = event.get(pf_candidates_aos_token_);
      // filter out eta domain and estimate mem block size
      size_t size = std::count_if(pf_candidates_aos.begin(), pf_candidates_aos.end(), [](l1t::PFCandidate c) {
        return c.eta() > -2.4 && c.eta() < 2.4;
      });
      // allocate buffer to store converted soa
      auto pf_candidates_soa = PuppiHostCollection(size, event.queue());
      alpaka::wait(event.queue());

      // convert aos to soa
      size_t idx_target = 0;
      for (size_t idx = 0; idx < pf_candidates_aos.size(); ++idx) {
        const auto &pf_candidate = pf_candidates_aos[idx];
        if (std::abs(pf_candidate.eta()) >= 2.4)
          continue;
        auto pf_view = pf_candidates_soa.view()[idx_target];

        pf_view.pt() = static_cast<float>(pf_candidate.pt());
        pf_view.eta() = static_cast<float>(pf_candidate.eta());
        pf_view.phi() = static_cast<float>(pf_candidate.phi());
        pf_view.z0() = static_cast<float>(pf_candidate.z0());
        pf_view.dxy() = static_cast<float>(pf_candidate.dxy());
        pf_view.puppiw() = static_cast<float>(pf_candidate.puppiWeight());
        pf_view.quality() = static_cast<int16_t>(pf_candidate.hwQual());
        pf_view.pdgid() = static_cast<uint8_t>(pf_candidate.pdgId());

        ++idx_target;
      }

      // move converted soa to event storage
      event.emplace(pf_candidates_soa_token_, std::move(pf_candidates_soa));
    }

    static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<edm::InputTag>("src");
      desc.addUntracked<int>("environment", static_cast<int>(Environment::kProduction));
      descriptions.addWithDefaultLabel(desc);
    }

  private:
    // get host data
    const edm::EDGetTokenT<std::vector<l1t::PFCandidate>> pf_candidates_aos_token_;
    // produce host data that will be automatically copied to device
    const edm::EDPutTokenT<PuppiHostCollection> pf_candidates_soa_token_;
    // debug / test
    const Environment environment_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc

DEFINE_FWK_ALPAKA_MODULE(l1sc::PFCandidateAoSToSoA);