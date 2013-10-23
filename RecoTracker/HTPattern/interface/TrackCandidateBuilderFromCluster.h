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
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"


class TrackCandidateBuilderFromCluster {
    public:
        TrackCandidateBuilderFromCluster(const edm::ParameterSet &iConfig) ;

        void init(const edm::EventSetup& es, const SeedComparitor *ifilter) ;
        void setHits(const HTHitsSpher  &hitsHiRes, const HTHitsSpher &hitsLowRes, 
                     const HTHitMap     &mapHiRes,  const HTHitMap    &mapLowRes,
                     std::vector<bool>  &maskHiRes,  std::vector<bool>     &maskLowRes,
                     int etashift, int phishift) {
            hitsHiRes_ = &hitsHiRes; hitsLowRes_ = &hitsLowRes;
            mapHiRes_ = &mapHiRes;  mapLowRes_ = &mapLowRes;
            maskHiRes_ = &maskHiRes; maskLowRes_ = &maskLowRes;
            etashift_ = etashift; phishift_ = phishift;
        }
        void run(const HTCluster &cluster, TrajectorySeedCollection & seedCollection, TrajectorySeedCollection *seedsFromAllClusters=0) ;

        struct ClusteredHit {
            ClusteredHit(const HTHitsSpher &hits, unsigned int i, unsigned int ly, float etaErr=0.02) :
                id(i), hit(hits.hit(i)), layer(ly), eta(hits.eta(i)), etaerr(etaErr), phi(hits.phi(i)), rho(hits.rho(i)), z(hits.z(i)) {};
            unsigned int id;
            const TrackingRecHit *hit;
            int layer;
            float eta, etaerr, phi, rho, z;
            float dist2d(const ClusteredHit &other) {
                float dphi = std::abs(phi-other.phi); 
                if (dphi > float(M_PI)) dphi = float(2*M_PI) - dphi;
                return dphi + std::abs(eta-other.eta);
            }
            float dist2d(float eta0, float phi0, float alpha, float beta) {
                float dphi = std::abs(phi-alpha*rho-phi0); 
                if (dphi > float(M_PI)) dphi = float(2*M_PI) - dphi;
                return dphi + std::abs(eta-beta*rho - eta0);
            }
        };
    protected:
        // --- Parameters ---
        std::string propagatorLabel_, hitBuilderLabel_;
        float dCut_, dcCut_;
        unsigned int minHits_;
        std::vector<double> startingCovariance_;

        // --- ES Stuff ---
        const SeedComparitor *filter = nullptr;
        edm::ESHandle<TrackerGeometry> tracker_;
        edm::ESHandle<Propagator> propagator_;
        edm::ESHandle<MagneticField> bfield_;
        edm::ESHandle<TransientTrackingRecHitBuilder> hitBuilder_;

        const HTHitsSpher  *hitsHiRes_, *hitsLowRes_;
        const HTHitMap *mapHiRes_,  *mapLowRes_;
        std::vector<bool> *maskHiRes_, *maskLowRes_;
        unsigned int etashift_, phishift_;

        // --- Helper Functions ---
        FreeTrajectoryState startingState(float eta0, float phi0, float alpha) const ;
        void refitAlphaBetaCorr(const std::vector<ClusteredHit> &cluster, const int *ihits, unsigned int nhits, float &alphacorr, float &betacorr, float &phi0, float &eta0) const ; 
        void dumpAsSeed(const std::vector<ClusteredHit> &cluster, TrajectorySeedCollection &seedCollection) const ;
         
};

#endif
