
//

/**
  \class    HTFinalOutcomeStudy HTFinalOutcomeStudy.h "RecoTracker/HTPattern/interface/HTFinalOutcomeStudy.h"
            
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
#include "RecoTracker/DebugTools/interface/TrackMixingAssociator.h"

class HTFinalOutcomeStudy : public edm::EDProducer {
    public:
      explicit HTFinalOutcomeStudy(const edm::ParameterSet & iConfig);
      virtual ~HTFinalOutcomeStudy() { }

      virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup) override;

     
    private:
      /// Labels for input collections
      edm::EDGetTokenT<std::vector<reco::Track> > tracksCkf_;
      edm::EDGetTokenT<std::vector<reco::Track> > tracksHT_;

      // Configurables
      StringCutObjectSelector<reco::Track> trackSelectionCkf_;

      // EventSetup stuff
      edm::ESHandle<TrackerGeometry> geometry_;
      edm::ESHandle<TrackerTopology> tTopo_;
      edm::ESHandle<MagneticField> bfield_;
};

HTFinalOutcomeStudy::HTFinalOutcomeStudy(const edm::ParameterSet & iConfig) :
    tracksCkf_(consumes<std::vector<reco::Track> >(iConfig.getParameter<edm::InputTag>("tracksCkf"))),
    tracksHT_(consumes<std::vector<reco::Track> >(iConfig.getParameter<edm::InputTag>("tracksHT"))),
    trackSelectionCkf_(iConfig.getParameter<std::string>("trackSelectionCkf"))
{
}

void 
HTFinalOutcomeStudy::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {
    edm::Handle<std::vector<reco::Track> > tracksHT;
    iEvent.getByToken(tracksHT_, tracksHT);
    edm::Handle<std::vector<reco::Track> > tracksCkf;
    iEvent.getByToken(tracksCkf_, tracksCkf);

    TrackMixingAssociator associator;
    associator.registerTrackEvent(1, *tracksCkf);
    associator.registerTrackEvent(2, *tracksHT);
    std::vector<TrackMixingAssociator::TrackAssociation> assoc;

    printf("\n\n========================================== CKF -> HT SUMMARY ==================================");
    for (const reco::Track & tk : *tracksCkf) {
        if (!trackSelectionCkf_(tk)) continue;
        associator.associateToTracks(tk, assoc);

        printf("\nTrack pt %8.2f +/- %7.2f, eta %+5.3f, phi %+5.3f, q %+2d, hits %2d, layers %2d, 3d layers %2d, vz %+6.3f, algo %2d, %4s (rho = %5.1f, z = %+5.1f --> rho = %5.1f, z = %+5.1f)\n", 
                tk.pt(), tk.ptError(), tk.eta(), tk.phi(), tk.charge(), tk.numberOfValidHits(), tk.hitPattern().trackerLayersWithMeasurement(), tk.hitPattern().pixelLayersWithMeasurement() + tk.hitPattern().numberOfValidStripLayersWithMonoAndStereo(), 
                tk.vz(), tk.algo(), tk.quality(reco::Track::qualitySize) ? "good" : "fake", tk.innerPosition().Rho(), tk.innerPosition().Z(), tk.outerPosition().Rho(), tk.outerPosition().Z()); 

        for (const auto &a : assoc) {
            if (a.sharedHits <= 2 || a.eventId == 1) continue;
            printf("   -> pt %8.2f +/- %7.2f, eta %+5.3f, phi %+5.3f, q %+2d, hits %2d, layers %2d, 3d layers %2d, vz %+6.3f: associated layers %2d (%3.0f%%), hits %2d (%3.0f%%)\n", 
                a.track->pt(), a.track->ptError(), a.track->eta(), a.track->phi(), a.track->charge(), a.track->hitPattern().numberOfValidHits(), a.track->hitPattern().trackerLayersWithMeasurement(), a.track->hitPattern().pixelLayersWithMeasurement() + a.track->hitPattern().numberOfValidStripLayersWithMonoAndStereo(), 
                a.track->vz(), a.sharedHitPattern.trackerLayersWithMeasurement(), a.sharedHitPattern.trackerLayersWithMeasurement()*100./tk.hitPattern().trackerLayersWithMeasurement(), a.sharedHits, (100.0*a.sharedHits)/tk.found());
        }
    }

    printf("\n\n========================================== HT -> CKF SUMMARY ==================================");
    for (const reco::Track & tk : *tracksHT) {
        associator.associateToTracks(tk, assoc);

        printf("\nTrack pt %8.2f +/- %7.2f, eta %+5.3f, phi %+5.3f, q %+2d, hits %2d, layers %2d, 3d layers %2d, vz %+6.3f, algo %2d, %4s (rho = %5.1f, z = %+5.1f --> rho = %5.1f, z = %+5.1f)\n", 
                tk.pt(), tk.ptError(), tk.eta(), tk.phi(), tk.charge(), tk.numberOfValidHits(), tk.hitPattern().trackerLayersWithMeasurement(), tk.hitPattern().pixelLayersWithMeasurement() + tk.hitPattern().numberOfValidStripLayersWithMonoAndStereo(), 
                tk.vz(), tk.algo(), tk.quality(reco::Track::qualitySize) ? "good" : "fake", tk.innerPosition().Rho(), tk.innerPosition().Z(), tk.outerPosition().Rho(), tk.outerPosition().Z()); 

        for (const auto &a : assoc) {
            if (a.sharedHits <= 2 || a.eventId == 2) continue;
            printf("   -> pt %8.2f +/- %7.2f, eta %+5.3f, phi %+5.3f, q %+2d, hits %2d, layers %2d, 3d layers %2d, vz %+6.3f: associated layers %2d (%3.0f%%), hits %2d (%3.0f%%)\n", 
                a.track->pt(), a.track->ptError(), a.track->eta(), a.track->phi(), a.track->charge(), a.track->hitPattern().numberOfValidHits(), a.track->hitPattern().trackerLayersWithMeasurement(), a.track->hitPattern().pixelLayersWithMeasurement() + a.track->hitPattern().numberOfValidStripLayersWithMonoAndStereo(), 
                a.track->vz(), a.sharedHitPattern.trackerLayersWithMeasurement(), a.sharedHitPattern.trackerLayersWithMeasurement()*100./a.track->hitPattern().trackerLayersWithMeasurement(), a.sharedHits, (100.0*a.sharedHits)/a.track->found());
        }
    }
    printf("\n\n");

}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HTFinalOutcomeStudy);
