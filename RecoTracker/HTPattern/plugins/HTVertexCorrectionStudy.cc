
//

/**
  \class    HTVertexCorrectionStudy HTVertexCorrectionStudy.h "RecoTracker/HTPattern/interface/HTVertexCorrectionStudy.h"
            
  \author   Giovanni Petrucciani
*/


#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

class HTVertexCorrectionStudy : public edm::EDProducer {
    public:
      explicit HTVertexCorrectionStudy(const edm::ParameterSet & iConfig);
      virtual ~HTVertexCorrectionStudy() { }

      virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup) override;

     
    private:
      /// Labels for input collections
      edm::EDGetTokenT<std::vector<reco::Track> > tracks_;

      // Configurables
      StringCutObjectSelector<reco::Track> trackSelection_;

      Float_t dz_, dxy_, phi_;
};

HTVertexCorrectionStudy::HTVertexCorrectionStudy(const edm::ParameterSet & iConfig) :
    tracks_(consumes<std::vector<reco::Track> >(iConfig.getParameter<edm::InputTag>("tracks"))),
    trackSelection_(iConfig.getParameter<std::string>("trackSelection")),
    dz_(iConfig.getParameter<double>("dz")),
    dxy_(iConfig.getParameter<double>("dxy")),
    phi_(iConfig.getParameter<double>("phi"))
{
    produces<std::vector<reco::Vertex> >();
}

void 
HTVertexCorrectionStudy::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {
    double x = 0, y = 0, z = 0; int ntracks = 0;
    edm::Handle<std::vector<reco::Track> > tracks;
    iEvent.getByToken(tracks_, tracks);
    for (const reco::Track & tk : *tracks) {
        if (!trackSelection_(tk)) continue;
        x += tk.vx(); 
        y += tk.vy(); 
        z += tk.vz(); 
        ntracks++;
    }
    std::auto_ptr<std::vector<reco::Vertex> > out(new std::vector<reco::Vertex>());
    out->push_back(reco::Vertex(
        reco::Vertex::Point(x/ntracks - dxy_ * std::cos(phi_), 
                            y/ntracks - dxy_ * std::sin(phi_), 
                            z/ntracks - dz_),
        reco::Vertex::Error(),
        0.0, 10., ntracks));

    for (const reco::Track & tk : *tracks) {
        if (!trackSelection_(tk)) continue;
        printf("track pt = %7.2f, eta %+5.3f, phi %+5.3f:    dz = %+7.4f, dxy = %+7.4f\n",
                tk.pt(), tk.eta(), tk.phi(),
                tk.dz(out->back().position()), tk.dxy(out->back().position()));
    }
    iEvent.put(out);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HTVertexCorrectionStudy);
