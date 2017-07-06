
// system include files
#include <memory>
#include <cmath>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

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
#include "../interface/BadComponents.h"

#include "DataFormats/Scalers/interface/LumiScalers.h"
#include "DataFormats/Luminosity/interface/LumiDetails.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include <TTree.h>

class DebugPixelHits : public edm::one::EDAnalyzer<edm::one::SharedResources> {
    public:
        explicit DebugPixelHits(const edm::ParameterSet&);
        ~DebugPixelHits();

    private:
        virtual void analyze(const edm::Event&, const edm::EventSetup&) ;

        // ----------member data ---------------------------
        const edm::EDGetTokenT<edm::View<reco::Candidate>> pairs_;    
        const edm::EDGetTokenT<edmNew::DetSetVector<SiPixelCluster>> pixelClusterLabel_;
        const edm::EDGetTokenT<MeasurementTrackerEvent> tracker_;
        const edm::EDGetTokenT<LumiScalersCollection> lumiScaler_;
        const edm::EDGetTokenT<std::vector<reco::Vertex>> vertices_; 
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
        BadComponents theBadComponents;

        TTree *tree_;
        int run_, lumi_, bx_, instLumi_, npv_; uint64_t event_; 
        float pair_mass_, track_pt_, track_eta_, track_phi_;
        float pv_ndof_, pv_z0_, tag_dz_, tag_dxy_, pv_xyErr_, pv_zErr_, track_dz_, track_dxy_, track_dzErr_, track_dxyErr_;
        int source_det_, source_layer_;
        int detid_, roc_, hitFound_, hitOnTrack_, detIsActive_, maybeBadROC_, trackHasHit_, trackHasLostHit_, hitInRandomWindow_;
        float track_global_phi_, track_global_z_, track_local_x_, track_local_y_, track_exp_sizeX_, track_exp_sizeY_, track_exp_charge_;
        float hit_global_phi_, hit_global_z_, hit_local_x_, hit_local_y_, hit_sizeX_, hit_sizeY_, hit_firstpixel_x_, hit_firstpixel_y_, hit_chi2_, hit_charge_, hitInRandomWindowDistance_;

        int debug_;
};

//
// constructors and destructor
//
DebugPixelHits::DebugPixelHits(const edm::ParameterSet& iConfig):
    pairs_(consumes<edm::View<reco::Candidate>>(iConfig.getParameter<edm::InputTag>("pairs"))),
    //pixelClusterLabel_(consumes<edmNew::DetSetVector<SiPixelCluster>>(iConfig.getParameter<edm::InputTag>("pixelClusters"))),
    tracker_(consumes<MeasurementTrackerEvent>(iConfig.getParameter<edm::InputTag>("tracker"))),
    lumiScaler_(consumes<LumiScalersCollection>(iConfig.getParameter<edm::InputTag>("lumiScalers"))),
    vertices_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices"))),
    refitter_(iConfig),
    propagatorOpposite_(iConfig.getParameter<std::string>("PropagatorOpposite")),
    estimatorName_(iConfig.getParameter<std::string>("Chi2MeasurementEstimator")),
    rescaleError_(iConfig.getParameter<double>("rescaleError")),
    debug_(iConfig.getUntrackedParameter<int>("debug",0))
    //pixelQualityLabel_(iConfig.getParameter<std::string>("SiPixelQuality"))
{
    if (iConfig.existsAs<std::string>("badComponentsFile")) {
        theBadComponents.init(iConfig.getParameter<std::string>("badComponentsFile"));
    }

    usesResource("TFileService");
    edm::Service<TFileService> fs;
    tree_ = fs->make<TTree>("tree","tree");
    tree_->Branch("run",  &run_,  "run/i");
    tree_->Branch("lumi", &lumi_, "lumi/i");
    tree_->Branch("bx", &bx_, "bx/i");
    tree_->Branch("instLumi", &instLumi_, "instLumi/i");
    tree_->Branch("npv", &npv_, "npv/i");
    tree_->Branch("pair_mass", &pair_mass_, "pair_mass/F");
    tree_->Branch("track_pt", &track_pt_, "track_pt/F");
    tree_->Branch("track_eta", &track_eta_, "track_eta/F");
    tree_->Branch("track_phi", &track_phi_, "track_phi/F");
    tree_->Branch("track_global_phi", &track_global_phi_, "track_global_phi/F");
    tree_->Branch("track_global_z", &track_global_z_, "track_global_z/F");
    tree_->Branch("track_local_x", &track_local_x_, "track_local_x/F");
    tree_->Branch("track_local_y", &track_local_y_, "track_local_y/F");
    tree_->Branch("track_exp_sizeX", &track_exp_sizeX_, "track_exp_sizeX/F");
    tree_->Branch("track_exp_sizeY", &track_exp_sizeY_, "track_exp_sizeY/F");
    tree_->Branch("track_exp_charge", &track_exp_charge_, "track_exp_charge/F");
    tree_->Branch("hit_global_phi", &hit_global_phi_, "hit_global_phi/F");
    tree_->Branch("hit_global_z", &hit_global_z_, "hit_global_z/F");
    tree_->Branch("hit_local_x", &hit_local_x_, "hit_local_x/F");
    tree_->Branch("hit_local_y", &hit_local_y_, "hit_local_y/F");
    tree_->Branch("hit_firstpixel_x", &hit_firstpixel_x_, "hit_firstpixel_x/F");
    tree_->Branch("hit_firstpixel_y", &hit_firstpixel_y_, "hit_firstpixel_y/F");
    tree_->Branch("hit_sizeX", &hit_sizeX_, "hit_sizeX/F");
    tree_->Branch("hit_sizeY", &hit_sizeY_, "hit_sizeY/F");
    tree_->Branch("hit_chi2", &hit_chi2_, "hit_chi2/F");
    tree_->Branch("hit_charge", &hit_charge_, "hit_charge/F");

    tree_->Branch("pv_ndof", &pv_ndof_, "pv_ndof/F");
    tree_->Branch("pv_z0", &pv_z0_, "pv_z0/F");
    tree_->Branch("tag_dz", &tag_dz_, "tag_dz/F");
    tree_->Branch("tag_dxy", &tag_dxy_, "tag_dxy/F");
    tree_->Branch("pv_xyErr", &pv_xyErr_, "pv_xyErr/F");
    tree_->Branch("pv_zErr", &pv_zErr_, "pv_zErr/F");
    tree_->Branch("track_dz", &track_dz_, "track_dz/F");
    tree_->Branch("track_dxy", &track_dxy_, "track_dxy/F");
    tree_->Branch("track_dzErr", &track_dzErr_, "track_dzErr/F");
    tree_->Branch("track_dxyErr", &track_dxyErr_, "track_dxyErr/F");

    tree_->Branch("detid", &detid_, "detid/I");
    tree_->Branch("roc", &roc_, "roc/I");
    tree_->Branch("source_det", &source_det_, "source_det/I");
    tree_->Branch("source_layer", &source_layer_, "source_layer/I");
    tree_->Branch("hitFound", &hitFound_, "hitFound/I");
    tree_->Branch("hitOnTrack", &hitOnTrack_, "hitOnTrack/I");
    tree_->Branch("hitInRandomWindow", &hitInRandomWindow_, "hitInRandomWindow/I");
    tree_->Branch("hitInRandomWindowDistance", &hitInRandomWindowDistance_, "hitInRandomWindowDistance/F");
    tree_->Branch("detIsActive", &detIsActive_, "detIsActive/I");
    tree_->Branch("maybeBadROC", &maybeBadROC_, "maybeBadROC/I");
    tree_->Branch("trackHasHit", &trackHasHit_, "trackHasHit/I");
    tree_->Branch("trackHasLostHit", &trackHasLostHit_, "trackHasLostHit/I");

}


DebugPixelHits::~DebugPixelHits()
{
}

    void
DebugPixelHits::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    run_  = iEvent.id().run();
    lumi_ = iEvent.id().luminosityBlock();
    bx_   = iEvent.bunchCrossing();

    edm::Handle<LumiScalersCollection> lumiScaler;
    iEvent.getByToken(lumiScaler_, lumiScaler);
    instLumi_ = (!lumiScaler->empty()) ? lumiScaler->front().instantLumi() : -99;

    edm::Handle<std::vector<reco::Vertex>> vertices;
    iEvent.getByToken(vertices_, vertices);
    npv_ = vertices->size();

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

    std::vector<int> badROCs;
    for (const reco::Candidate & pair : *pairs) {
        pair_mass_ = pair.mass();
        const reco::Muon &tag = dynamic_cast<const reco::Muon &>(*pair.daughter(0)->masterClone());
        if (tag.innerTrack().isNull()) continue;
        const reco::Muon &mu = dynamic_cast<const reco::Muon &>(*pair.daughter(1)->masterClone());
        if (mu.innerTrack().isNull()) continue;
        const reco::Vertex *bestPV = &vertices->front();
        float bestDZ = std::abs(tag.innerTrack()->dz(bestPV->position()));
        for (const reco::Vertex &pv : *vertices) {
            float thisDZ = std::abs(tag.innerTrack()->dz(pv.position()));
            if (thisDZ < bestDZ) { bestDZ = thisDZ; bestPV = &pv; }
        }
        pv_ndof_ = bestPV->ndof();
        pv_z0_ = bestPV->z();
        pv_xyErr_ = std::hypot(bestPV->xError(), bestPV->yError());
        pv_zErr_ = bestPV->zError();
        tag_dz_ = tag.innerTrack()->dz(bestPV->position());
        tag_dxy_ = tag.innerTrack()->dxy(bestPV->position());
        const reco::Track & mutk = *mu.innerTrack();
        track_pt_ = mutk.pt();
        track_eta_ = mutk.eta();
        track_phi_ = mutk.phi();
        track_dz_ = mutk.dz(bestPV->position());
        track_dxy_ = mutk.dxy(bestPV->position());
        track_dzErr_ = mutk.dzError();
        track_dxyErr_ = mutk.dxyError();
        trackHasHit_ = mutk.hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 1);
        trackHasLostHit_ = mutk.hitPattern().numberOfLostPixelBarrelHits(reco::HitPattern::MISSING_INNER_HITS);
        bool PXB1Valid =  mutk.hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 1);
        bool PXB2Valid =  mutk.hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 2);
        if (debug_) printf("\n\nMuon with pt %f eta %f, found PXB hits: %d, lost PXB inner hits along: %d, lost PXB inner hits before: %d, PXB1 hit: %d, PXB2 hit: %d, PXB3 hit: %d, PXF1 hit: %d, PXF2 hit: %d\n",
                            mu.pt(),  mu.eta(), mutk.hitPattern().numberOfValidPixelBarrelHits(), 
                            mutk.hitPattern().numberOfLostPixelBarrelHits(reco::HitPattern::TRACK_HITS), mutk.hitPattern().numberOfLostPixelBarrelHits(reco::HitPattern::MISSING_INNER_HITS), 
                            PXB1Valid, PXB2Valid, mutk.hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 3),
                            mutk.hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelEndcap, 1), 
                            mutk.hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelEndcap, 2));
        int nhits = mutk.recHitsSize();
        std::vector<Trajectory> traj  = refitter_.transform(mutk);
        if (traj.size() != 1) continue; 
        std::vector<const SiPixelCluster *> PXB1Clusters;
        if (PXB1Valid) {
            for (int i = 0; i < nhits; ++i) {
                const TrackingRecHit *hit = &* mutk.recHit(i);
                if (hit->isValid() && hit->geographicalId().subdetId() == PixelSubdetector::SubDetector::PixelBarrel && 
                        theTrkTopo->layer(hit->geographicalId()) == 1) { 
                    const auto * pixhit = dynamic_cast<const SiPixelRecHit*>(hit->hit());
                    if (!pixhit) throw cms::Exception("CorruptData", "Valid PXB1 hit that is not a SiPixelRecHit");
                    auto clustref = pixhit->cluster();
                    if (clustref.isNull()) throw cms::Exception("CorruptData", "Valid PXB1 SiPixelRecHit with null cluster ref");
                    PXB1Clusters.push_back(&*clustref);
                }
            }
            if (PXB1Clusters.empty()) std::cout << "WARNING: did not find Cluster for PXB1 Hit" << std::endl;
        }
        TrajectoryStateOnSurface tsosPXB2;
        for (const auto &tm : traj.front().measurements()) {
            if (tm.recHit().get() && tm.recHitR().isValid()) {
                DetId where = tm.recHitR().geographicalId();
                source_det_ = where.subdetId();
                source_layer_ = theTrkTopo->layer(where);
                if (source_det_ != PixelSubdetector::SubDetector::PixelBarrel ||  source_layer_ != 1) {
                    tsosPXB2 = tm.updatedState().isValid() ? tm.updatedState() : tm.backwardPredictedState();
                    if (debug_) printf("starting state on det %d, layer %d, r = %5.1f, z = %+6.1f\n", 
                                        source_det_, source_layer_, tsosPXB2.globalPosition().perp(), tsosPXB2.globalPosition().z());
                    break;
                }
            }
        }
        if (!tsosPXB2.isValid()) std::cout << "WARNING: did not find state for PXB2 Hit" << std::endl;
        if (!tsosPXB2.isValid()) continue; // for now
        tsosPXB2.rescaleError(rescaleError_);
        const GeometricSearchTracker * gst = tracker->geometricSearchTracker();
        //const auto & PXBLs = gst->pixelBarrelLayers();
        //for (const auto * PXBLayer : PXBLs) { std::cout << "PXB Layer with radius = " << PXBLayer->specificSurface().radius() << std::endl; }
        const auto *pxbLayer1 = gst->pixelBarrelLayers().front();
        auto compDets = pxbLayer1->compatibleDets(tsosPXB2, *thePropagatorOpposite, *theEstimator);
        for (const auto & detAndState : compDets) {
            bool found = false;
            if (debug_) printf("Compatible module %u\n", detAndState.first->geographicalId().rawId()); 
            detid_ = detAndState.first->geographicalId().rawId();
            const auto &mdet = tracker->idToDet(detAndState.first->geographicalId());
            detIsActive_ = mdet.isActive();
            if (debug_) printf("   module center: rho = % 6.2f  z = %+6.3f  phi = %+6.3f\n", mdet.position().perp(),  mdet.position().z(), mdet.position().phi().value());
            if (debug_) printf("   track state:   rho = % 6.2f  z = %+6.3f  phi = %+6.3f\n", detAndState.second.globalPosition().perp(), detAndState.second.globalPosition().z(), detAndState.second.globalPosition().phi().value());
            const auto & tkpos = detAndState.second.localPosition();
            const auto & tkerr = detAndState.second.localError().positionError();
            if (debug_) printf("              local x = %+6.3f +- %5.3f   y = %+6.3f +- %5.3f \n",  tkpos.x(), std::sqrt(tkerr.xx()),  tkpos.y(), std::sqrt(tkerr.yy()));
            LocalVector localDir = detAndState.second.localMomentum()/detAndState.second.localMomentum().mag();
            float chargeCorr = std::abs(localDir.z());
            track_global_phi_ = detAndState.second.globalPosition().phi();
            track_global_z_   = detAndState.second.globalPosition().z();
            track_local_x_ = tkpos.x();
            track_local_y_ = tkpos.y();
            roc_ = (tkpos.x() > 0)*8 + std::min(std::max<int>(0,std::floor(tkpos.y()+3.24)),7);
            track_exp_sizeX_ = 1.5f;
            track_exp_sizeY_ = std::max<float>(1.f, std::hypot(1.3f, 1.9f*mutk.pz()/mutk.pt()));
            track_exp_charge_ = std::sqrt(1.08f + std::pow(mutk.pz()/mutk.pt(),2)) * 26000;
            theBadComponents.fetch(iEvent.id().run(), iEvent.id().luminosityBlock(),  detAndState.first->geographicalId().rawId(), badROCs);
            maybeBadROC_ = false;
            for (int badROC: badROCs) {
                float rocx = ((badROC/8)-0.5)*82*0.0100; // one ROC is 80 pixels, each 100um wide, but there are big pixels
                float rocy = ((badROC%8)-3.5)*54*0.0150; // one ROC is 52 pixels, each 150um wide, but there are big pixels
                if (debug_) printf("      bad ROC local x = %+6.3f            y = %+6.3f    (gio) \n",  rocx, rocy);
                //MeasurementPoint rocixiy((badROC/8)*80+40, (badROC%8)*52 + 26);
                //LocalPoint rocxy = (dynamic_cast<const PixelGeomDetUnit*>(detAndState.first))->specificTopology().localPosition(MeasurementPoint(rocixiy));
                //printf("      bad ROC local x = %+6.3f            y = %+6.3f    (top) \n",  rocxy.x(), rocxy.y());
                if (std::abs(rocx - tkpos.x()) < (40*0.0100 + 2*0.0100 + 5*std::sqrt(tkerr.xx()))) maybeBadROC_ = true; // 1 half-ROC + 2 pixels + 5 sigma
                if (std::abs(rocy - tkpos.y()) < (26*0.0150 + 2*0.0150 + 5*std::sqrt(tkerr.yy()))) maybeBadROC_ = true; // 1 half-ROC + 2 pixels + 5 sigma
            }
            //auto rechits = mdet.recHits(detAndState.second);
            //for (const auto & hit : rechits) {
            std::vector<std::pair<float, TrackingRecHit::ConstRecHitPointer>> rechitsAndChi2 ;
            for (const auto & hitAndChi2 : mdet.fastMeasurements(detAndState.second, tsosPXB2, *thePropagatorOpposite, *theEstimator)) {
                if (hitAndChi2.recHit()->isValid()) rechitsAndChi2.emplace_back(hitAndChi2.estimate(), hitAndChi2.recHit());
            } 
            const auto & allhits = mdet.recHits(detAndState.second);
            if (rechitsAndChi2.empty()) {
                for (const auto & hit : allhits) { 
                    float distance = std::max(std::abs(hit->localPosition().x()-tkpos.x()), std::abs(hit->localPosition().y()-tkpos.y()));
                    if (hit->isValid() && distance < 0.2) {
                         rechitsAndChi2.emplace_back(99 + distance, hit);
                    }
                }
            }
            std::sort(rechitsAndChi2.begin(), rechitsAndChi2.end());
            // now we make the matching around a fake position, in the same TBM; go inwards by one ROC if one is not in the innermost ROC
            bool innermostROC = std::abs(tkpos.y()) < 0.8;
            float fakey = tkpos.y() + std::copysign(tkpos.y(), 0.8)*(innermostROC ? +1 : -1); 
            hitInRandomWindow_ = false; hitInRandomWindowDistance_ = 0.2;
            for (const auto & hit : allhits) { 
                float distance = std::max(std::abs(hit->localPosition().x()-tkpos.x()), std::abs(hit->localPosition().y()-fakey));
                if (hit->isValid() && distance < hitInRandomWindowDistance_) {
                    hitInRandomWindow_ = true; hitInRandomWindowDistance_ = distance;
                }
            }
            for (const auto & hitAndChi2 : rechitsAndChi2) {
                const auto & hit = hitAndChi2.second;
                const auto &hitpos = hit->localPosition();
                const auto &hiterr = hit->localPositionError();
                const auto * pixhit = dynamic_cast<const SiPixelRecHit*>(&*hit);    if (!pixhit) throw cms::Exception("CorruptData", "Valid PXB1 hit that is not a SiPixelRecHit");
                auto clustref = pixhit->cluster();                                  if (clustref.isNull()) throw cms::Exception("CorruptData", "Valid PXB1 SiPixelRecHit with null cluster ref");
                auto hitgpos = detAndState.first->toGlobal(hitpos);
                hit_global_phi_ = hitgpos.phi();
                hit_global_z_   = hitgpos.z();
                hit_local_x_ = hitpos.x();
                hit_local_y_ = hitpos.y();
                hit_firstpixel_x_ = clustref->minPixelRow();
                hit_firstpixel_y_ = clustref->minPixelCol();
                hit_chi2_ = hitAndChi2.first;
                hit_charge_ = clustref->charge();
                hit_sizeX_ = clustref->sizeX();
                hit_sizeY_ = clustref->sizeY();
                hitFound_ = true;
                hitOnTrack_ = (std::find(PXB1Clusters.begin(), PXB1Clusters.end(), clustref.get()) != PXB1Clusters.end());
                if (!found) { tree_->Fill(); found = true; }
                if (debug_) printf("             rechit x = %+6.3f +- %5.3f   y = %+6.3f +- %5.3f     dx = %+6.3f +- %5.3f   dy = %+6.3f +- %5.3f   chi2 = %7.2f    ontrack = %c  ",  
                                    hitpos.x(), std::sqrt(hiterr.xx()), hitpos.y(), std::sqrt(hiterr.yy()),
                                    hitpos.x()-tkpos.x(), std::sqrt(hiterr.xx()+tkerr.xx()), hitpos.y()-tkpos.y(), std::sqrt(hiterr.yy()+tkerr.yy()), hitAndChi2.first,
                                    (hitOnTrack_ ? 'Y' : 'n'));
                if (debug_) printf(" size: %3d (%2d x %2d), raw charge %8d   corr charge %8.1f\n", clustref->size(), clustref->sizeX(), clustref->sizeY(), clustref->charge(), clustref->charge()*chargeCorr);
                if (debug_) {
                    for (const auto &P : clustref->pixels()) {
                        printf("      pixel at x = %3d y = %3d \n", P.x, P.y);
                    }
                }
            }
            if (!found) {
                hitFound_ = false;
                tree_->Fill();
            }
        }
    } 

    if (debug_ > 0) debug_--;
}


//define this as a plug-in
DEFINE_FWK_MODULE(DebugPixelHits);
