#ifndef TrackCandidateBuilderFromCluster_h
#define TrackCandidateBuilderFromCluster_h

#include "RecoTracker/HTPattern/interface/HTHitMap.h"

class TrackCandidateBuilderFromCluster {
    public:
        TrackCandidateBuilderFromCluster() {}
        void setHits(const HTHitsSpher  &hitsHiRes, const HTHitsSpher &hitsLowRes, 
                     const HTHitMap     &mapHiRes,  const HTHitMap    &mapLowRes,
                     std::vector<bool>  &maskHiRes,  std::vector<bool>     &maskLowRes,
                     int etashift, int phishift) {
            hitsHiRes_ = &hitsHiRes; hitsLowRes_ = &hitsLowRes;
            mapHiRes_ = &mapHiRes;  mapLowRes_ = &mapLowRes;
            maskHiRes_ = &maskHiRes; maskLowRes_ = &maskLowRes;
            etashift_ = etashift; phishift_ = phishift;
        }
        void run(const HTCluster &cluster) ;

    protected:
        // --- ES Stuff ---
        // geometry
        // propagator
        // updator
        // filter
        // estimator
        // ttrhbuilder

        const HTHitsSpher  *hitsHiRes_, *hitsLowRes_;
        const HTHitMap *mapHiRes_,  *mapLowRes_;
        std::vector<bool> *maskHiRes_, *maskLowRes_;
        unsigned int etashift_, phishift_;
    private:
         
};

#endif
