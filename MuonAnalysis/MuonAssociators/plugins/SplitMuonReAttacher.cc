
//
// $Id: SplitMuonReAttacher.cc,v 1.5 2010/07/12 20:56:11 gpetrucc Exp $
//

/**
  \class    pat::SplitMuonReAttacher SplitMuonReAttacher.h "MuonAnalysis/MuonAssociators/interface/SplitMuonReAttacher.h"
  \brief    Matcher of reconstructed objects to other reconstructed objects using the tracks inside them 
            
  \author   Giovanni Petrucciani
  \version  $Id: SplitMuonReAttacher.cc,v 1.5 2010/07/12 20:56:11 gpetrucc Exp $
*/


#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "RecoTracker/TrackProducer/interface/TrackMerger.h"


// template-related workaround for bug in OwnVector+Ptr
namespace edm { using std::advance; }

namespace pat {

  class SplitMuonReAttacher : public edm::EDProducer {
    public:
      explicit SplitMuonReAttacher(const edm::ParameterSet & iConfig);
      virtual ~SplitMuonReAttacher() { }

      virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);

    private:
      /// Labels for input collections
      edm::InputTag src_;

      /// Use reco info to determine who's inner and who's outer 
      bool useReco_;
    
      /// Use segment references for matching
      bool useSegmentReferences_;

      TrackMerger merger_;
  };
  
} // namespace

pat::SplitMuonReAttacher::SplitMuonReAttacher(const edm::ParameterSet & iConfig) :
    src_(iConfig.getParameter<edm::InputTag>("src")),
    useReco_(iConfig.existsAs<bool>("useReco") ? iConfig.getParameter<bool>("useReco") : false),
    merger_(iConfig)
{
    produces<std::vector<TrackCandidate> >(); 
    if (!useReco_) throw cms::Exception("Configuration") << "Sorry, RECO is needed this morning\n";
}

void 
pat::SplitMuonReAttacher::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {
    using namespace edm;
    using namespace std;

    Handle<View<reco::Muon> > src;
    iEvent.getByLabel(src_, src);

    auto_ptr<vector<TrackCandidate> > out(new vector<TrackCandidate>());
    unsigned int nsrc = src->size();
    out->reserve(nsrc);

    merger_.init(iSetup);

    for (unsigned int i = 0; i < nsrc; ++i) {
        const reco::Muon &mu1 = (*src)[i];
        if (mu1.track().isNull()) continue;
        for (unsigned int j = i+1; j < nsrc; ++j) {
            const reco::Muon &mu2 = (*src)[j];
            if (mu2.track().isNull()) continue;
            // see if the two match
            if (deltaR(mu1,mu2) > 0.03) continue;
            if (mu1.charge() != mu2.charge()) continue;
            // figure out who's inner and who's outer
            reco::TrackRef inner = mu1.track(), outer = mu2.track();
            if (useReco_) {
                if (mu1.track()->outerPosition().R() > mu2.track()->outerPosition().R()) {
                    std::swap(inner, outer);
                }
            }
            std::cout << "Attempting to merge muons.\n";
            std::cout << "    Inner pt " << inner->pt() << ", eta " << inner->eta() << ", phi " << inner->phi() << ",  hits " << inner->hitPattern().numberOfHits() << std::endl;
            std::cout << "    Outer pt " << outer->pt() << ", eta " << outer->eta() << ", phi " << outer->phi() << ",  hits " << outer->hitPattern().numberOfHits() << std::endl;
            out->push_back(merger_.merge(*inner, *outer));
            std::cout << "    output candidate has  " << (out->back().recHits().second - out->back().recHits().first) << " hits" << std::endl;
        }
    }
    iEvent.put(out);
}

#include "FWCore/Framework/interface/MakerMacros.h"
using namespace pat;
DEFINE_FWK_MODULE(SplitMuonReAttacher);
