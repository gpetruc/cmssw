#ifndef TrackCandidateBuilderFromCluster_h
#define TrackCandidateBuilderFromCluster_h

#include "RecoTracker/HTPattern/interface/HTHitMap.h"
#include "RecoTracker/TkSeedingLayers/interface/SeedComparitor.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimatorBase.h"
#include "TrackingTools/DetLayers/interface/NavigationSetter.h"
#include "TrackingTools/DetLayers/interface/NavigationSchool.h"
#include "RecoTracker/CkfPattern/interface/BaseCkfTrajectoryBuilder.h"
#include "TrackingTools/TrajectoryCleaning/interface/TrajectoryCleaner.h"
#include "RecoTracker/Record/interface/NavigationSchoolRecord.h"
#include "RecoTracker/MeasurementDet/interface/MeasurementTrackerEvent.h"
#include "RecoTracker/CkfPattern/interface/TransientInitialStateEstimator.h"


class TrackCandidateBuilderFromCluster {
    public:
        TrackCandidateBuilderFromCluster(const edm::ParameterSet &iConfig) ;
        ~TrackCandidateBuilderFromCluster() { delete initialStateEstimator_; }

        void init(const edm::EventSetup& es, const MeasurementTrackerEvent &evt, const SeedComparitor *ifilter) ;
        void done();
        void setHits(const HTHitsSpher  &hitsHiRes, const HTHitsSpher &hitsLowRes, 
                     const HTHitMap     &mapHiRes,  const HTHitMap    &mapLowRes,
                     std::vector<bool>  &maskHiRes,  std::vector<bool>     &maskLowRes,
                     int etashift, int phishift, bool useLowRes) {
            hitsHiRes_ = &hitsHiRes; hitsLowRes_ = &hitsLowRes;
            mapHiRes_ = &mapHiRes;  mapLowRes_ = &mapLowRes;
            maskHiRes_ = &maskHiRes; maskLowRes_ = &maskLowRes;
            etashift_ = etashift; phishift_ = phishift;
            useLowRes_ = useLowRes;
            bfieldAtOrigin_ = bfield_->inTesla(GlobalPoint(hitsHiRes.x0(), hitsHiRes.y0(), hitsHiRes.z0())).z();
        }
        void setAlphaRange(float alphaMin, float alphaMax) {
            alphaMin_ = alphaMin; alphaMax_ = alphaMax;
            if (alphaMax_ < alphaMin_) std::swap(alphaMax_,alphaMin_);
        }
        void run(const HTCluster &cluster, unsigned int nseedlayersCut, unsigned int nlayersCut, TrackCandidateCollection & tcCollection, TrajectorySeedCollection & seedCollection, TrajectorySeedCollection *seedsFromAllClusters=0) ;

        struct ClusteredHit {
            ClusteredHit(const HTHitsSpher &hits, unsigned int i, unsigned int ly, float etaErr=1) :
                id(i), hit(hits.hit(i)), layer(ly), eta(hits.eta(i)), etaerr(etaErr), phi(hits.phi(i)), rho(hits.rho(i)), z(hits.z(i)) {};
            unsigned int id;
            const TrackingRecHit *hit;
            int layer;
            float eta, etaerr, phi, rho, z;
            float dist2d(const ClusteredHit &other) {
                return dphi(other) + deta(other);
            }
            float dist2d(float eta0, float phi0, float alpha, float beta) {
                return dphi(phi0, alpha) + deta(eta0, beta);
            }
            float dphi(const ClusteredHit &other) {
                return dphi(other.phi, 0);
            }
            float dphi(float phi0, float alpha) {
                float ret = std::abs(phi-alpha*rho-phi0);
                while (ret > float(2*M_PI)) ret -= float(2*M_PI);
                if (ret > float(M_PI)) ret = float(2*M_PI) - ret;
                return ret;
            }
            float deta(const ClusteredHit &other) {
                return std::abs(eta-other.eta);
            }
            float deta(float eta0, float beta) {
                return std::abs(eta - beta/rho - eta0);
            }
            float detaE(float eta0, float beta) {
                return std::abs(eta - beta/rho - eta0)*etaerr;
            }
        };

        static int hitid(const TrackingRecHit *hit) ;
    protected:
        // --- Parameters ---
        std::string propagatorLabel_, propagatorOppositeLabel_, hitBuilderLabel_, estimatorLabel_, navigationSchoolLabel_, trajectoryBuilderLabel_, trajectoryCleanerLabel_;
        bool cleanTrajectoryAfterInOut_;
        edm::ParameterSet initialStateEstimatorPSet_;
        bool useHitsSplitting_, splitSeedHits_;
        float dphiCutPair_,    detaCutPair_;
        bool  pairsOnSeedCellOnly_;
        float dphiCutHits_[2], detaCutHits_[2];
        unsigned int minHits_, minHitsIfNoCkf_;
        unsigned int maxFailMicroClusters_;
        std::vector<double> startingCovariance_;

        // --- ES Stuff ---
        const SeedComparitor *filter = nullptr;
        edm::ESHandle<TrackerGeometry> tracker_;
        edm::ESHandle<Propagator> propagator_, propagatorOpposite_, anyPropagator_;
        edm::ESHandle<MagneticField> bfield_;
        edm::ESHandle<TransientTrackingRecHitBuilder> hitBuilder_;
        edm::ESHandle<Chi2MeasurementEstimatorBase>   estimator_;
        edm::ESHandle<NavigationSchool>  navigationSchool_;
        edm::ESHandle<TrajectoryBuilder> trajectoryBuilderTemplate_;
        edm::ESHandle<TrajectoryCleaner> trajectoryCleaner_;
        TransientInitialStateEstimator*  initialStateEstimator_;

        const HTHitsSpher  *hitsHiRes_, *hitsLowRes_;
        const HTHitMap *mapHiRes_,  *mapLowRes_;
        std::vector<bool> *maskHiRes_, *maskLowRes_;
        unsigned int etashift_, phishift_;
        bool useLowRes_;
        float bfieldAtOrigin_;
        float alphaMin_, alphaMax_;

        //std::vector<bool> *maskStripClusters_, *maskPixelClusters_;
        
        // per-event locals
        std::auto_ptr<NavigationSetter> navigationSetter_;
        std::auto_ptr<BaseCkfTrajectoryBuilder> trajectoryBuilder_;

        // --- Helper Functions ---
        typedef std::vector<std::pair<uint32_t,unsigned int>> HitIndex;
        FreeTrajectoryState startingState(float eta0, float phi0, float alpha) const ;
        void refitAlphaBetaCorr(const std::vector<ClusteredHit> &cluster, const int *ihits, unsigned int nhits, float &alphacorr, float &betacorr, float &phi0, float &eta0) const ; 
        void dumpAsSeed(const std::vector<ClusteredHit> &cluster, TrajectorySeedCollection &seedCollection) const ;
        const TrajectorySeed * makeSeed(const std::vector<ClusteredHit> &hits, std::pair<int,int> seedhits, const std::array<int,20> &layerhits,  float eta0, float phi0, float alpha, TrajectorySeedCollection &out) const ;
        const TrajectorySeed * makeSeed2Way(const std::vector<ClusteredHit> &hits, std::pair<int,int> seedhits, const std::array<int,20> &layerhits,  float eta0, float phi0, float alpha, TrajectorySeedCollection &out) const ;
        void makeTrajectories(const TrajectorySeed &seed, std::vector<Trajectory> &out) const ;
        void indexHits(const std::vector<ClusteredHit> &hits, HitIndex &index) const ;
        int  findHit(const TrackingRecHit *hit, const std::vector<ClusteredHit> &hits, HitIndex &index) const ;
        void saveCandidates(const std::vector<Trajectory> &cands,  TrackCandidateCollection & tcCollection) const ;  
        void dumpTraj(const Trajectory &traj) const ;
        bool prefilterCluster(const HTCluster &cluster) const ;
};

#endif
