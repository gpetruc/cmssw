#ifndef SimplisticTrajectoryBuilder_h
#define SimplisticTrajectoryBuilder_h

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
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "TrackingTools/PatternTools/interface/TempTrajectory.h"
#include <TStopwatch.h>

#ifdef  HTDEBUG
#define DEBUG HTDEBUG
#else
#define DEBUG 0.5
#endif
#define DEBUG1_printf  if (DEBUG>=1) printf
#define DEBUG2_printf  if (DEBUG>=2) printf
#define DEBUG3_printf  if (DEBUG>=3) printf
#define DEBUG4_printf  if (DEBUG>=4) printf

namespace reco { struct Track; }

class SimplisticTrajectoryBuilder {
    public:
        SimplisticTrajectoryBuilder(const edm::ParameterSet &iConfig) ;
        ~SimplisticTrajectoryBuilder() ;

        void setEvent(const MeasurementTrackerEvent &iEvent, const edm::EventSetup &iSetup) ;

        bool run(const std::vector<const TrackingRecHit *> & hits, const TrajectoryStateOnSurface & stateOnFirstHit, PropagationDirection startingDirection, std::vector<Trajectory> &out) const ;

        void debugSetTrack(const reco::Track *target) const { targetTrack_ = target; }
    protected:
        enum LayerPolicy { anyLayer, pixelLayer, noLayer };
        static LayerPolicy str2policy(const std::string &p) { if (p == "anyLayer") return anyLayer; if (p == "pixelLayer") return pixelLayer; return noLayer; }

        // --- Parameters ---
        std::string propagatorLabel_, propagatorOppositeLabel_, hitBuilderLabel_, estimatorLabel_, estimatorRebuildLabel_, trajectoryCleanerLabel_, trajectoryFilterLabel_, trajectoryFilterStartLabel_, trajectoryFilterRebuildLabel_;
        double foundHitBonus_, lostHitPenalty_;
        LayerPolicy searchHits_, useGrouped_;

        // --- Event ---
        const MeasurementTrackerEvent * mtEvent_;

        // --- ES Stuff ---
        edm::ESHandle<TrackerGeometry> tracker_;
        edm::ESHandle<Propagator> propagator_, propagatorOpposite_;
        edm::ESHandle<MagneticField> bfield_;
        edm::ESHandle<TransientTrackingRecHitBuilder> hitBuilder_;
        edm::ESHandle<Chi2MeasurementEstimatorBase>   estimator_, estimatorRebuild_;
        edm::ESHandle<TrajectoryCleaner> trajectoryCleaner_;
        edm::ESHandle<TrajectoryFilter>  trajectoryFilter_, trajectoryFilterStart_, trajectoryFilterRebuild_;
        edm::ESHandle<TrackerTopology>   tTopo_;

        // --- Debug stuff ---
        mutable const reco::Track * targetTrack_; // the pointee is const, the pointer is mutable

        void fitStartingHits(const std::vector<const TrackingRecHit *> & hits, const TrajectoryStateOnSurface & stateOnFirstHit, TempTrajectory & traj) const ;
        void continueTrajectory(TempTrajectory & traj) const ;
        void rebuildTrajectory(const TempTrajectory &src, TempTrajectory & out) const ;
        bool nextLayer(TempTrajectory & traj, const Chi2MeasurementEstimatorBase &est, const TransientTrackingRecHit *hint) const ;
};

#endif
