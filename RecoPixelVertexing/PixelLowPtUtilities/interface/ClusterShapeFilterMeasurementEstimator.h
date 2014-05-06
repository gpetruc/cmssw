#ifndef ClusterShapeFilterMeasurementEstimator_h
#define ClusterShapeFilterMeasurementEstimator_h

#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimatorBase.h"
#include "RecoPixelVertexing/PixelLowPtUtilities/interface/ClusterShapeHitFilter.h"

class ClusterShapeFilterMeasurementEstimator : public Chi2MeasurementEstimatorBase {
    public:
        ClusterShapeFilterMeasurementEstimator(const ClusterShapeHitFilter *filter, 
                const Chi2MeasurementEstimatorBase  *base, bool filterPixelHits, bool filterStripHits) :
            Chi2MeasurementEstimatorBase(base ? base->chiSquaredCut() : 0, base ? base->nSigmaCut() : 0),
            filter_(filter), base_(base), filterPixelHits_(filterPixelHits), filterStripHits_(filterStripHits) {}

        /** Returns pair( true, value) if the TrajectoryStateOnSurface is compatible
         *  with the RecHit, and pair( false, value) if it is not compatible.
         *  The TrajectoryStateOnSurface must be on the same Surface as the RecHit. 
         *  For an estimator where there is no value computed, e.g. fixed
         *  window estimator, only the first(bool) part is of interest.
         */
        virtual HitReturnType estimate( const TrajectoryStateOnSurface& ts, 
                const TransientTrackingRecHit& hit) const override;

        /** Returns true if the TrajectoryStateOnSurface is compatible with the
         *  Plane, false otherwise.
         *  The TrajectoryStateOnSurface must be on the plane.  */
        virtual SurfaceReturnType estimate( const TrajectoryStateOnSurface& ts, 
                const Plane& plane) const override { 
            return base_ ? base_->estimate(ts,plane) : true; 
        }

        virtual ClusterShapeFilterMeasurementEstimator* clone() const override {
            return new ClusterShapeFilterMeasurementEstimator(*this);
        }
    protected:
        const ClusterShapeHitFilter *filter_;
        const Chi2MeasurementEstimatorBase *base_;
        bool  filterPixelHits_, filterStripHits_;

        /// check if the hit is compatible
        bool  compatibleHit(const TransientTrackingRecHit &hit, const TrajectoryStateOnSurface& ts) const ;
};

#endif
