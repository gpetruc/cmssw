
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
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorBase.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"


class HTFinalOutcomeStudy : public edm::EDProducer {
    public:
      explicit HTFinalOutcomeStudy(const edm::ParameterSet & iConfig);
      virtual ~HTFinalOutcomeStudy() { }

      virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup) override;

     
    private:
      /// Labels for input collections
      edm::EDGetTokenT<std::vector<reco::Track> > tracksCkf_;
      edm::EDGetTokenT<edm::View<reco::Track> > tracksHT_;

      // Configurables
      StringCutObjectSelector<reco::Track> trackSelectionCkf_;
      StringCutObjectSelector<reco::Track> trackSelectionHT_;

      // MC Configurables   
      bool doMC_;
      edm::EDGetTokenT<std::vector<TrackingParticle> > tracksSim_;
      std::string simAssociatorLabel_;

      // EventSetup stuff
      edm::ESHandle<TrackerGeometry> geometry_;
      edm::ESHandle<TrackerTopology> tTopo_;
      edm::ESHandle<MagneticField> bfield_;
      edm::ESHandle<TrackAssociatorBase> simAssociator_;
};

HTFinalOutcomeStudy::HTFinalOutcomeStudy(const edm::ParameterSet & iConfig) :
    tracksCkf_(consumes<std::vector<reco::Track> >(iConfig.getParameter<edm::InputTag>("tracksCkf"))),
    tracksHT_(consumes<edm::View<reco::Track> >(iConfig.getParameter<edm::InputTag>("tracksHT"))),
    trackSelectionCkf_(iConfig.getParameter<std::string>("trackSelectionCkf")),
    trackSelectionHT_(iConfig.getParameter<std::string>("trackSelectionHT")),
    doMC_(iConfig.getParameter<bool>("doMC"))
{
    if (doMC_)  {
        tracksSim_ = consumes<std::vector<TrackingParticle> >(iConfig.getParameter<edm::InputTag>("tracksSim"));
        simAssociatorLabel_ = iConfig.getParameter<std::string>("simAssociator");
    }
}

void 
HTFinalOutcomeStudy::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {
    edm::Handle<edm::View<reco::Track> > tracksHT;
    iEvent.getByToken(tracksHT_, tracksHT);
    edm::Handle<std::vector<reco::Track> > tracksCkf;
    iEvent.getByToken(tracksCkf_, tracksCkf);

    TrackMixingAssociator associator;
    associator.registerTrackEvent(1, *tracksCkf);
    associator.registerTrackEvent(2, *tracksHT);
    std::vector<TrackMixingAssociator::TrackAssociation> assoc;

    if (doMC_) {
        printf("\n\n========================================== CKF -> HT SUMMARY ==================================");
        unsigned int ngood = 0, foundgood = 0, founddup = 0;
        iSetup.get<TrackAssociatorRecord>().get(simAssociatorLabel_, simAssociator_);
        edm::Handle<TrackingParticleCollection> tracksSim;
        iEvent.getByToken(tracksSim_, tracksSim);
        reco::SimToRecoCollection simRecColl = simAssociator_->associateSimToReco(tracksHT, tracksSim, &iEvent, &iSetup);
        //edm::ESHandle<ParametersDefinerForTP> parametersDefinerTP; 
        //iSetup.get<TrackAssociatorRecord>().get("LhcParametersDefinerForTP", parametersDefinerTP);    
        for (unsigned int i = 0, n = tracksSim->size(); i < n; ++i) {
            TrackingParticleRef tpr(tracksSim, i);
            const TrackingParticle &tp = *tpr;
            TrackingParticle::LorentzVector decay;
            if (!tp.decayVertices().empty()) decay = tp.decayVertices()[0]->position();
            printf("\nTrPar pt %9.4f +/- %8.4f, eta %+5.3f, phi %+5.3f, q %+2d, hits %2d, layers %2d, 3d layers %2d, vz %+6.3f, algo %2d, %4s (rho = %5.1f, z = %+5.1f --> rho = %5.1f, z = %+5.1f)\n", 
                    tp.pt(), 0.0, tp.eta(), tp.phi(), tp.charge(), tp.numberOfTrackerHits(), tp.numberOfTrackerLayers(), -1, tp.vz(), 0, "sim.", 
                    tp.vertex().Rho(), tp.vertex().Z(), decay.Rho(), decay.Z()); 
            auto match = simRecColl.find(tpr);
            bool found = false, dup = false;
            if (match != simRecColl.end()) {
                for (const auto & a : match->val) {
                    const reco::Track &ta = *a.first;
                    printf("   -> pt %9.4f +/- %8.4f, eta %+5.3f, phi %+5.3f, q %+2d, hits %2d, layers %2d, 3d layers %2d, vz %+6.3f: association quality %7.3f\n", 
                            ta.pt(), ta.ptError(), ta.eta(), ta.phi(), ta.charge(), ta.hitPattern().numberOfValidHits(), ta.hitPattern().trackerLayersWithMeasurement(), ta.hitPattern().pixelLayersWithMeasurement() + ta.hitPattern().numberOfValidStripLayersWithMonoAndStereo(), 
                            ta.vz(), a.second);
                    if (found) dup = true;
                    found = true;
                }
            }
            ngood++; if (found) foundgood++; if (dup) founddup++;
        }
        printf("\nREPORT SIM: good tracking particles %d  found %d  duplicated %d\n", ngood, foundgood, founddup);
    }

    unsigned int ngood = 0, foundgood = 0, foundfake = 0, founddup = 0;
    printf("\n\n========================================== CKF -> HT SUMMARY ==================================");
    for (const reco::Track & tk : *tracksCkf) {
        if (!trackSelectionCkf_(tk)) continue;
        associator.associateToTracks(tk, assoc);

        printf("\nTrack pt %9.4f +/- %8.4f, eta %+5.3f, phi %+5.3f, q %+2d, hits %2d, layers %2d, 3d layers %2d, vz %+6.3f, algo %2d, %4s (rho = %5.1f, z = %+5.1f --> rho = %5.1f, z = %+5.1f)\n", 
                tk.pt(), tk.ptError(), tk.eta(), tk.phi(), tk.charge(), tk.numberOfValidHits(), tk.hitPattern().trackerLayersWithMeasurement(), tk.hitPattern().pixelLayersWithMeasurement() + tk.hitPattern().numberOfValidStripLayersWithMonoAndStereo(), 
                tk.vz(), tk.algo(), tk.quality(reco::Track::qualitySize) ? "good" : "fake", tk.innerPosition().Rho(), tk.innerPosition().Z(), tk.outerPosition().Rho(), tk.outerPosition().Z()); 

        bool found = false, dup = false;
        for (const auto &a : assoc) {
            if (a.sharedHits <= 2 || a.eventId == 1) continue;
            printf("   -> pt %9.4f +/- %8.4f, eta %+5.3f, phi %+5.3f, q %+2d, hits %2d, layers %2d, 3d layers %2d, vz %+6.3f: associated layers %2d (%3.0f%%), hits %2d (%3.0f%%)\n", 
                a.track->pt(), a.track->ptError(), a.track->eta(), a.track->phi(), a.track->charge(), a.track->hitPattern().numberOfValidHits(), a.track->hitPattern().trackerLayersWithMeasurement(), a.track->hitPattern().pixelLayersWithMeasurement() + a.track->hitPattern().numberOfValidStripLayersWithMonoAndStereo(), 
                a.track->vz(), a.sharedHitPattern.trackerLayersWithMeasurement(), a.sharedHitPattern.trackerLayersWithMeasurement()*100./tk.hitPattern().trackerLayersWithMeasurement(), a.sharedHits, (100.0*a.sharedHits)/tk.found());
            int lay  = tk.hitPattern().trackerLayersWithMeasurement();
            int slay = a.sharedHitPattern.trackerLayersWithMeasurement();
            if (slay > 3 && (slay >= 8 || slay >= 0.6*lay)) {
                if (found) dup = true;
                found = true;
            }
        }
        if (tk.quality(reco::Track::qualitySize)) {
            ngood++;
            if (found) foundgood++;
            if (dup)   founddup++;
        }
    }

    printf("\n\n========================================== HT -> CKF SUMMARY ==================================");
    for (const reco::Track & tk : *tracksHT) {
        associator.associateToTracks(tk, assoc);

        printf("\nTrack pt %9.4f +/- %8.4f, eta %+5.3f, phi %+5.3f, q %+2d, hits %2d, layers %2d, 3d layers %2d, vz %+6.3f, algo %2d, %4s (rho = %5.1f, z = %+5.1f --> rho = %5.1f, z = %+5.1f)\n", 
                tk.pt(), tk.ptError(), tk.eta(), tk.phi(), tk.charge(), tk.numberOfValidHits(), tk.hitPattern().trackerLayersWithMeasurement(), tk.hitPattern().pixelLayersWithMeasurement() + tk.hitPattern().numberOfValidStripLayersWithMonoAndStereo(), 
                tk.vz(), tk.algo(), tk.quality(reco::Track::qualitySize) ? "good" : "fake", tk.innerPosition().Rho(), tk.innerPosition().Z(), tk.outerPosition().Rho(), tk.outerPosition().Z()); 

        for (const auto &a : assoc) {
            if (a.sharedHits <= 2 || a.eventId == 2) continue;
            printf("   -> pt %9.4f +/- %8.4f, eta %+5.3f, phi %+5.3f, q %+2d, hits %2d, layers %2d, 3d layers %2d, vz %+6.3f: associated layers %2d (%3.0f%%), hits %2d (%3.0f%%)\n", 
                a.track->pt(), a.track->ptError(), a.track->eta(), a.track->phi(), a.track->charge(), a.track->hitPattern().numberOfValidHits(), a.track->hitPattern().trackerLayersWithMeasurement(), a.track->hitPattern().pixelLayersWithMeasurement() + a.track->hitPattern().numberOfValidStripLayersWithMonoAndStereo(), 
                a.track->vz(), a.sharedHitPattern.trackerLayersWithMeasurement(), a.sharedHitPattern.trackerLayersWithMeasurement()*100./a.track->hitPattern().trackerLayersWithMeasurement(), a.sharedHits, (100.0*a.sharedHits)/a.track->found());
        }
        if (trackSelectionHT_(tk) && !tk.quality(reco::Track::qualitySize)) foundfake++;
    }
    printf("\n\n");

    printf("REPORT: good tracks %d  found %d  duplicated %d  fakes %d\n", ngood, foundgood, founddup, foundfake);

}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HTFinalOutcomeStudy);
