
//

/**
  \class    HTBinningStudies HTBinningStudies.h "RecoTracker/HTPattern/interface/HTBinningStudies.h"
            
  \author   Giovanni Petrucciani
*/


#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "RecoTracker/HTPattern/interface/HTHits.h"
#include "RecoTracker/HTPattern/interface/HTHitMap.h"
#include "RecoTracker/HTPattern/interface/HTDebugger.h"
#include "RecoTracker/HTPattern/interface/TrackCandidateBuilderFromCluster.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include <TTree.h>

class HTBinningStudies : public edm::EDProducer {
    public:
      explicit HTBinningStudies(const edm::ParameterSet & iConfig);
      virtual ~HTBinningStudies() { }

      virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup) override;

     
    private:
      /// Labels for input collections
      edm::EDGetTokenT<std::vector<reco::Track> > tracks_;

      // Configurables
      StringCutObjectSelector<reco::Track> trackSelection_;

      // EventSetup stuff
      edm::ESHandle<TrackerGeometry> geometry_;
      edm::ESHandle<TrackerTopology> tTopo_;
      edm::ESHandle<MagneticField> bfield_;

      TTree *tree_;
      Float_t  trackPt_, trackEta_;
      Float_t  deta_, dphi_; Int_t layer_;
};

HTBinningStudies::HTBinningStudies(const edm::ParameterSet & iConfig) :
    tracks_(consumes<std::vector<reco::Track> >(iConfig.getParameter<edm::InputTag>("tracks"))),
    trackSelection_(iConfig.getParameter<std::string>("trackSelection"))
{
    edm::Service<TFileService> fs;
    tree_ = fs->make<TTree>("t","t");
    tree_->Branch("pt",&trackPt_,"pt/F");
    tree_->Branch("eta",&trackEta_,"eta/F");
    tree_->Branch("deta",&deta_,"deta/F");
    tree_->Branch("dphi",&dphi_,"dphi/F");
    tree_->Branch("layer",&layer_,"layer/I");
}

void 
HTBinningStudies::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {
    //Get the geometry
    iSetup.get<TrackerDigiGeometryRecord>().get(geometry_);

    //Retrieve tracker topology from geometry
    iSetup.get<IdealGeometryRecord>().get(tTopo_);

    //Retrieve magnetic field
    iSetup.get<IdealMagneticFieldRecord>().get(bfield_);
    double bfield = bfield_->inTesla(GlobalPoint(0.,0.,0.)).z();

    edm::Handle<std::vector<reco::Track> > tracks;
    iEvent.getByToken(tracks_, tracks);
    for (const reco::Track & tk : *tracks) {
        if (!trackSelection_(tk)) continue;
        trackPt_ = tk.pt(); trackEta_ = tk.eta();
        GlobalPoint p0(tk.vx(), tk.vy(), tk.vz());
        float alpha = 0.5 * 0.003 * bfield / trackPt_;

        for (trackingRecHit_iterator ihit = tk.recHitsBegin(), ehit = tk.recHitsEnd(); ihit != ehit; ++ihit) {
            const TrackingRecHit *hit = &**ihit;
            if (!hit->isValid()) continue;
            DetId id = hit->geographicalId(); 
            const GeomDet* geomDet = geometry_->idToDet(id);
            if (geomDet == 0) continue;
            layer_ = tTopo_->layer(id);
            if (id.subdetId() > 2) layer_ += 3; // strip
            if (id.subdetId() == StripSubdetector::TOB || id.subdetId() == StripSubdetector::TEC) layer_ += 4;
            GlobalVector vec = geomDet->surface().toGlobal( hit->localPosition() ) - p0;
            deta_ = fabs(vec.eta() - tk.eta());
            dphi_ = fabs(deltaPhi(vec.phi()-alpha*tk.phi(), tk.phi()));
            printf("hit on layer %d: dphi = %7.5f\n", layer_, dphi_);
            tree_->Fill();
        }
    }
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HTBinningStudies);
