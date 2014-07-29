// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

class TrackCounter : public edm::EDAnalyzer {
   public:
      explicit TrackCounter(const edm::ParameterSet&);
      ~TrackCounter();

   private:
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

      // ----------member data ---------------------------
      std::string name_;
      edm::EDGetTokenT<reco::TrackCollection> tkToken_;
      std::vector<std::pair<std::string,StringCutObjectSelector<reco::Track>>> cuts_;
      std::map<std::string,unsigned int> overall_;
};

TrackCounter::TrackCounter(const edm::ParameterSet& iConfig):
    name_(iConfig.getParameter<std::string>("@module_label")),
    tkToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks")))
{
    edm::ParameterSet cuts = iConfig.getParameter<edm::ParameterSet>("cuts");
    std::vector<std::string> names = cuts.getParameterNamesForType<std::string>();
    std::sort(names.begin(),names.end());
    for (auto n : names) {
        cuts_.push_back(std::make_pair(n, StringCutObjectSelector<reco::Track>(cuts.getParameter<std::string>(n))));
    }
}

TrackCounter::~TrackCounter()
{
    printf("--- overall track report from %s ---\n", name_.c_str());
    for (const auto &p : overall_) {
        printf("%-30s %6u\n", p.first.c_str(), p.second);
    }
    printf("\n");
}


void
TrackCounter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edm::Handle<reco::TrackCollection> tracks;
    iEvent.getByToken(tkToken_, tracks);
   
    printf("--- track report from %s ---\n", name_.c_str());
    printf("%-30s %6d\n", "total tracks", int(tracks->size()));
    overall_["total"] += tracks->size();
    for (const auto &p : cuts_) {
        unsigned int pass = std::count_if(tracks->begin(), tracks->end(), p.second);
        printf("%-30s %6d\n", p.first.c_str(), int(pass));
        overall_[p.first] += pass;
    }
    printf("\n");
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackCounter);
