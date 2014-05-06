#include "RecoPixelVertexing/PixelLowPtUtilities/interface/ClusterShapeFilterMeasurementEstimator.h"
#include "RecoPixelVertexing/PixelLowPtUtilities/interface/ClusterShapeHitFilter.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/ProjectedSiStripRecHit2D.h"
#include <cassert>

ClusterShapeFilterMeasurementEstimator::HitReturnType
ClusterShapeFilterMeasurementEstimator::estimate( const TrajectoryStateOnSurface& ts, const TransientTrackingRecHit& hit) const 
{
    if (!compatibleHit(hit,ts)) return HitReturnType(false,9e9);
    return base_ ? base_->estimate(ts,hit) : HitReturnType(true,0);
}

bool 
ClusterShapeFilterMeasurementEstimator::compatibleHit(const TransientTrackingRecHit &hit, const TrajectoryStateOnSurface& tsos) const 
{
    assert(hit.isValid() && tsos.isValid());
    GlobalVector direction = tsos.globalDirection();
    if (hit.geographicalId().subdetId() <= 2) {
        if (!filterPixelHits_) return true;    
        const SiPixelRecHit *pixhit = dynamic_cast<const SiPixelRecHit *>(hit.hit());
        if (pixhit == 0) throw cms::Exception("LogicError", "Found a valid hit on the pixel detector which is not a SiPixelRecHit\n");
        //printf("Cheching hi hit on detid %10d, local direction is x = %9.6f, y = %9.6f, z = %9.6f\n", hit.geographicalId().rawId(), direction.x(), direction.y(), direction.z());
        return filter_->isCompatible(*pixhit, direction);
    } else {
        if (!filterStripHits_) return true;
        GlobalPoint position = tsos.globalPosition();
        const std::type_info &tid = typeid(*hit.hit());
        if (tid == typeid(SiStripMatchedRecHit2D)) {
            const SiStripMatchedRecHit2D* matchedHit = dynamic_cast<const SiStripMatchedRecHit2D *>(hit.hit());
            assert(matchedHit != 0);
            return (filter_->isCompatible(DetId(matchedHit->monoId()), matchedHit->monoCluster(), position, direction) &&
                    filter_->isCompatible(DetId(matchedHit->stereoId()), matchedHit->stereoCluster(), position, direction));
        } else if (tid == typeid(SiStripRecHit2D)) {
            const SiStripRecHit2D* recHit = dynamic_cast<const SiStripRecHit2D *>(hit.hit());
            assert(recHit != 0);
            return filter_->isCompatible(*recHit, position, direction);
        } else if (tid == typeid(ProjectedSiStripRecHit2D)) {
            const ProjectedSiStripRecHit2D* precHit = dynamic_cast<const ProjectedSiStripRecHit2D *>(hit.hit());
            assert(precHit != 0);
            return filter_->isCompatible(precHit->originalHit(), position, direction);
        } else {
            //printf("Questo e' un %s, che ci fo?\n", tid.name());
            return true;
        }
    }
}
