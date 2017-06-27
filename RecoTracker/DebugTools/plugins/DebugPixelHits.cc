
// system include files
#include <memory>
#include <cmath>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "TrackingTools/TrackRefitter/interface/TrackTransformer.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/TrackerRecHit2D/interface/TrackerSingleRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "RecoTracker/MeasurementDet/interface/MeasurementTrackerEvent.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimatorBase.h"
#include "TrackingTools/MeasurementDet/interface/TempMeasurements.h"


//#include "Geometry/TrackerGeometryBuilder/interface/GluedGeomDet.h"
#include "Geometry/CommonTopologies/interface/Topology.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
//#include "CalibFormats/SiPixelObjects/interface/SiPixelQuality.h"
//#include "CalibTracker/Records/interface/SiPixelQualityRcd.h"
#include "RecoTracker/MeasurementDet/src/TkMeasurementDetSet.h"


class DebugPixelHits : public edm::EDProducer {
    public:
        explicit DebugPixelHits(const edm::ParameterSet&);
        ~DebugPixelHits();

    private:
        virtual void produce(edm::Event&, const edm::EventSetup&);

        // ----------member data ---------------------------
        const edm::EDGetTokenT<edm::View<reco::Candidate>> pairs_;    
        const edm::EDGetTokenT<edmNew::DetSetVector<SiPixelCluster>> pixelClusterLabel_;
        edm::EDGetTokenT<MeasurementTrackerEvent> tracker_;
        /// Track Transformer
        TrackTransformer refitter_;
        const std::string propagatorOpposite_;
        //std::string pixelQualityLabel_;
        const std::string estimatorName_;
        const double rescaleError_;

        edm::ESHandle<TrackerTopology> theTrkTopo;
        edm::ESHandle<Propagator> thePropagatorOpposite;
        edm::ESHandle<Chi2MeasurementEstimatorBase> theEstimator;

        //edm::ESHandle<SiPixelQuality> thePixelQuality;
};

//
// constructors and destructor
//
DebugPixelHits::DebugPixelHits(const edm::ParameterSet& iConfig):
    pairs_(consumes<edm::View<reco::Candidate>>(iConfig.getParameter<edm::InputTag>("pairs"))),
    //pixelClusterLabel_(consumes<edmNew::DetSetVector<SiPixelCluster>>(iConfig.getParameter<edm::InputTag>("pixelClusters"))),
    tracker_(consumes<MeasurementTrackerEvent>(iConfig.getParameter<edm::InputTag>("tracker"))),
    refitter_(iConfig),
    propagatorOpposite_(iConfig.getParameter<std::string>("PropagatorOpposite")),
    estimatorName_(iConfig.getParameter<std::string>("Chi2MeasurementEstimator")),
    rescaleError_(iConfig.getParameter<double>("rescaleError"))
    //pixelQualityLabel_(iConfig.getParameter<std::string>("SiPixelQuality"))
{
}


DebugPixelHits::~DebugPixelHits()
{
}

    void
DebugPixelHits::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;

    // read input
    Handle<View<reco::Candidate> > pairs;
    iEvent.getByToken(pairs_, pairs);

    if (pairs->empty()) return;

    iSetup.get<TrackerTopologyRcd>().get(theTrkTopo);
    refitter_.setServices(iSetup);

    iSetup.get<TrackingComponentsRecord>().get(estimatorName_, theEstimator);
    iSetup.get<TrackingComponentsRecord>().get(propagatorOpposite_, thePropagatorOpposite);
    //iSetup.get<SiPixelQualityRcd>().get(pixelQualityLabel_, thePixelQuality);

    //edm::Handle<edmNew::DetSetVector<SiPixelCluster> > pixelC; 
    //iEvent.getByToken(pixelClusterLabel_, pixelC); 

    Handle<MeasurementTrackerEvent> tracker;
    iEvent.getByToken(tracker_, tracker);
    //const PxMeasurementConditionSet & pixelConds = tracker->measurementTracker().pixelDetConditions();

    for (const reco::Candidate & pair : *pairs) {
        const reco::Muon &mu = dynamic_cast<const reco::Muon &>(*pair.daughter(1)->masterClone());
        if (mu.innerTrack().isNull()) continue;
        const reco::Track & mutk = *mu.innerTrack();
        bool PXB1Valid =  mutk.hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 1);
        bool PXB2Valid =  mutk.hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 2);
        std::cout << "\n\nMuon with pt " << mu.pt() << " eta " << mu.eta() << 
            ", found PXB hits: " << mutk.hitPattern().numberOfValidPixelBarrelHits() <<
            ", lost PXB inner hits along: "  << mutk.hitPattern().numberOfLostPixelBarrelHits(reco::HitPattern::TRACK_HITS) << 
            ", lost PXB inner hits before: " << mutk.hitPattern().numberOfLostPixelBarrelHits(reco::HitPattern::MISSING_INNER_HITS) << 
            ", PXB1 hit: " << PXB1Valid << 
            ", PXB2 hit: " << PXB2Valid << 
            ", PXB3 hit: " << mutk.hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 3) << std::endl;
        int nhits = mutk.recHitsSize();
        std::vector<Trajectory> traj  = refitter_.transform(mutk);
        if (traj.size() != 1) continue; 
        std::cout << "Refitted ok" << std::endl;
        const SiPixelCluster *PXB1Cluster = nullptr;
        if (PXB1Valid) {
            for (int i = 0; i < nhits; ++i) {
                const TrackingRecHit *hit = &* mutk.recHit(i);
                if (hit->isValid() && hit->geographicalId().subdetId() == PixelSubdetector::SubDetector::PixelBarrel && 
                        theTrkTopo->layer(hit->geographicalId()) == 1) { 
                    const auto * pixhit = dynamic_cast<const SiPixelRecHit*>(hit->hit());
                    if (!pixhit) throw cms::Exception("CorruptData", "Valid PXB1 hit that is not a SiPixelRecHit");
                    auto clustref = pixhit->cluster();
                    if (clustref.isNull()) throw cms::Exception("CorruptData", "Valid PXB1 SiPixelRecHit with null cluster ref");
                    PXB1Cluster = &*clustref;
                    break;
                }
            }
            if (!PXB1Cluster) std::cout << "WARNING: did not find Cluster for PXB1 Hit" << std::endl;
        }
        TrajectoryStateOnSurface tsosPXB2;
        if (PXB2Valid) {
            for (const auto &tm : traj.front().measurements()) {
                if (tm.recHit().get() && tm.recHitR().isValid()) {
                    DetId where = tm.recHitR().geographicalId();
                    if (where.subdetId() == PixelSubdetector::SubDetector::PixelBarrel && theTrkTopo->layer(where) == 2) {
                        tsosPXB2 = tm.updatedState().isValid() ? tm.updatedState() : tm.backwardPredictedState();
                        break;
                    }
                }
            }
            if (!tsosPXB2.isValid()) std::cout << "WARNING: did not find state for PXB2 Hit" << std::endl;
        }
        if (!tsosPXB2.isValid()) continue; // for now
        tsosPXB2.rescaleError(rescaleError_);
        const GeometricSearchTracker * gst = tracker->geometricSearchTracker();
        //const auto & PXBLs = gst->pixelBarrelLayers();
        //for (const auto * PXBLayer : PXBLs) { std::cout << "PXB Layer with radius = " << PXBLayer->specificSurface().radius() << std::endl; }
        const auto *pxbLayer1 = gst->pixelBarrelLayers().front();
        auto compDets = pxbLayer1->compatibleDets(tsosPXB2, *thePropagatorOpposite, *theEstimator);
        for (const auto & detAndState : compDets) {
            std::cout << "Compatible module " << detAndState.first->geographicalId().rawId() << std::endl; 
            const auto &mdet = tracker->idToDet(detAndState.first->geographicalId());
            printf("   module center: rho = % 6.2f  z = %+6.3f  phi = %+6.3f\n", mdet.position().perp(),  mdet.position().z(), mdet.position().phi().value());
            printf("   track state:   rho = % 6.2f  z = %+6.3f  phi = %+6.3f\n", detAndState.second.globalPosition().perp(), detAndState.second.globalPosition().z(), detAndState.second.globalPosition().phi().value());
            const auto & tkpos = detAndState.second.localPosition();
            const auto & tkerr = detAndState.second.localError().positionError();
            printf("              local x = %+6.3f +- %5.3f   y = %+6.3f +- %5.3f \n",  tkpos.x(), std::sqrt(tkerr.xx()),  tkpos.y(), std::sqrt(tkerr.yy()));
            //auto rechits = mdet.recHits(detAndState.second);
            //for (const auto & hit : rechits) {
            LocalVector localDir = detAndState.second.localMomentum()/detAndState.second.localMomentum().mag();
            float chargeCorr = std::abs(localDir.z());
            std::vector<std::pair<float, TrackingRecHit::ConstRecHitPointer>> rechitsAndChi2 ;
            for (const auto & hitAndChi2 : mdet.fastMeasurements(detAndState.second, tsosPXB2, *thePropagatorOpposite, *theEstimator)) {
                if (hitAndChi2.recHit()->isValid()) rechitsAndChi2.emplace_back(hitAndChi2.estimate(), hitAndChi2.recHit());
            } 
            if (rechitsAndChi2.empty()) {
                for (const auto & hit : mdet.recHits(detAndState.second)) { 
                    if (hit->isValid() && std::max(std::abs(hit->localPosition().x()-tkpos.x()), std::abs(hit->localPosition().y()-tkpos.y())) < 0.2) {
                         rechitsAndChi2.emplace_back(-99, hit);
                    }
                }
            }
            std::sort(rechitsAndChi2.begin(), rechitsAndChi2.end());
            for (const auto & hitAndChi2 : rechitsAndChi2) {
                const auto & hit = hitAndChi2.second;
                const auto &hitpos = hit->localPosition();
                const auto &hiterr = hit->localPositionError();
                const auto * pixhit = dynamic_cast<const SiPixelRecHit*>(&*hit);    if (!pixhit) throw cms::Exception("CorruptData", "Valid PXB1 hit that is not a SiPixelRecHit");
                auto clustref = pixhit->cluster();                                  if (clustref.isNull()) throw cms::Exception("CorruptData", "Valid PXB1 SiPixelRecHit with null cluster ref");
                printf("             rechit x = %+6.3f +- %5.3f   y = %+6.3f +- %5.3f     dx = %+6.3f +- %5.3f   dy = %+6.3f +- %5.3f   chi2 = %7.2f    ontrack = %c  ",  
                            hitpos.x(), std::sqrt(hiterr.xx()), hitpos.y(), std::sqrt(hiterr.yy()),
                            hitpos.x()-tkpos.x(), std::sqrt(hiterr.xx()+tkerr.xx()), hitpos.y()-tkpos.y(), std::sqrt(hiterr.yy()+tkerr.yy()), hitAndChi2.first,
                            (clustref.get() == PXB1Cluster ? 'Y' : 'n'));
                printf(" size: %3d (%2d x %2d), raw charge %8d   corr charge %8.1f\n", clustref->size(), clustref->sizeX(), clustref->sizeY(), clustref->charge(), clustref->charge()*chargeCorr);
            }
        }
#if 0
        for (const GeomDet *det : gdets) {
            where = det->geographicalId(); mdet = tracker->idToDet(where);
            int denseIndex = pixelConds.find(where());
            if (denseIndex == pixelConds.nDet()) { std::cout << "Module missing in pixel conditions set" << std::endl; continue; }
            int nPixels = pixelConds.totalPixels(denseIndex);
            std::cout << "Analyzing module at " << where() << ", isActive? " << mdet.isActive() << ", pixels " << nPixels << std::endl;
            float utraj = 0, uerr = 0; unsigned int uapv = 0; bool pred = false, hascluster = false, hasdigi = false, anycluster = false, anydigi = false; int cmode = -1;
            if (tsosBefore.isValid()) {
                TrajectoryStateOnSurface tsos = thePropagator->propagate(tsosBefore, det->surface());
                if (tsos.isValid()) {  
                    pred = true;
                    utraj = det->topology().measurementPosition( tsos.localPosition() ).x();
                    uerr  = std::sqrt( det->topology().measurementError( tsos.localPosition(), tsos.localError().positionError() ).uu() ); 
                    uapv = std::min<unsigned int>(nPixels-1,std::max<float>(0,utraj))/128;
                    std::cout << "  Searching around pixel " << utraj << " +/- " << uerr << "    APV: " << uapv << std::endl;
                } else {
                    std::cout << "  Failed to propagate??" << std::endl;
                }
            }
            if (!mdet.isActive()) {
                std::cout << "  Detector is inactive" << std::endl;
                continue;
            }
            std::cout << "  Bad components on the detector" << std::endl;
            std::cout << "  APVs (or fibers): ";
            for (unsigned int iapv = 0; iapv < 5; ++iapv) {
                if (thePixelQuality->IsApvBad(where(), iapv) || thePixelQuality->IsFiberBad(where(), iapv/2)) std::cout << iapv << " " ;
            } 
            std::cout << std::endl;
            std::cout << "  Pixels: ";
            SiPixelQuality::Range range = thePixelQuality->getRange(where());
            for (auto pixel = range.first; pixel < range.second; ++pixel) std::cout << (*pixel) << " " ;
            std::cout << std::endl;
            auto cl_iter = pixelC->find(where);
            if (cl_iter == pixelC->end()) {
                std::cout << "  ... no pixel clusters on this detid" << std::endl;
            } else {
                edmNew::DetSet<SiPixelCluster> clusters = *cl_iter;
                for (const SiPixelCluster &cluster : clusters) {
                    std::cout << "  Cluster of " << cluster.amplitudes().size() << " pixels: " << std::endl;
                    const std::vector<uint8_t> & amps = cluster.amplitudes();
                    for (unsigned int s = cluster.firstPixel(), i = 0, e  = amps.size(); i < e; ++s, ++i) {  
                        std::cout << "   " << std::setw(4) << s << " | " << (s/128) << " | "; bar(amps[i], 2);
                        if (pred && std::abs(s-utraj) < 5) { hascluster = true; anycluster = true; }
                        if (pred && (s/128) == uapv) anycluster = true;
                    }
                }
            }
            auto di_iter = pixelD->find(where);
            if (di_iter == pixelD->end()) {
                std::cout << "  ... no pixel digis on this detid" << std::endl;
            } else {
                std::cout << "  Digis on this detid" << std::endl;
                const edm::DetSet<SiPixelDigi> & digis = *di_iter;
                for (unsigned int idigi = 0, ndigi = digis.size(); idigi < ndigi; ++idigi) {
                    if (idigi > 0 && (digis[idigi].pixel() > digis[idigi-1].pixel()+1)) std::cout << "      ---------------------" << std::endl;
                    std::cout << "   " << std::setw(4) << digis[idigi].pixel() << " | " << (digis[idigi].pixel()/128) << " | "; bar(digis[idigi].adc(), 2, 1024);
                    if (pred && std::abs(digis[idigi].pixel()-utraj) < 5) { hasdigi = true; anydigi = true; }
                    if (pred && (digis[idigi].pixel()/128) == uapv) anydigi = true;
                }
            }
            auto cm_iter = pixelCM->find(where);
            if (cm_iter == pixelCM->end()) {
                std::cout << "  ... no pixel common mode on this detid" << std::endl;
            } else {
                std::cout << "  Common mode on this detid" << std::endl;
                const edm::DetSet<SiPixelRawDigi> & apvs = *cm_iter;
                for (unsigned int iapv = 0, napv = apvs.size(); iapv < napv; ++iapv) {
                    std::cout << "   " << std::setw(4) << "..." << " | " << (iapv) << " | "; bar(apvs[iapv].adc(), 2, 1024);
                    if (pred && (iapv == uapv)) cmode = apvs[iapv].adc();
                }
            }
            if (pred) {
                std::cout << "  Summary: " << ( hascluster  ? " cluster" : (anycluster ? " " : " no-clusters")) << 
                    ( hasdigi  ? " digi" : (anydigi ? " " : " no-digis"));
                if (cmode >= 0) std::cout <<  " CM=" << cmode;
                if (utraj < 0 || utraj > nPixels) std::cout << " maybe-outside" ;
                std::cout << std::endl;
            }

        }
        for (const auto &tm : traj.front().measurements()) {
            if (tm.recHit().get() && !tm.recHit()->isValid()) {
                DetId where =  tm.recHitR().geographicalId();
                int subdet = where.subdetId(), layer = theTrkTopo->layer(where);
                if (!layersToDebug_.empty()) {
                    if (std::find(layersToDebug_.begin(), layersToDebug_.end(), sDETS[subdet]+sLAYS[layer]) == layersToDebug_.end()) {
                        continue;
                    }
                }
                std::cout << " missing hit on " << sDETS[subdet] << " layer " << layer << " type = " << tm.recHitR().getType() << ", detid " << where() << std::endl;
            }
        }
#endif
    } 
}


//define this as a plug-in
DEFINE_FWK_MODULE(DebugPixelHits);
