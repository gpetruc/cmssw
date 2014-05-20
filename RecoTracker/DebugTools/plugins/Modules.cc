#include "RecoTracker/DebugTools/interface/CkfDebugTrackCandidateMaker.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ModuleFactory.h"

using cms::CkfDebugTrackCandidateMaker;

//DEFINE_FWK_MODULE(CkfDebugTrackCandidateMaker);

// TrajectoryFilter
#include "FWCore/Utilities/interface/typelookup.h"

#include "TrackingTools/TrajectoryFiltering/interface/TrajectoryFilterFactory.h"
#include "RecoTracker/DebugTools/interface/ClusterShapeDebugTrajectoryFilter.h"
DEFINE_EDM_PLUGIN(TrajectoryFilterFactory, ClusterShapeDebugTrajectoryFilter, "ClusterShapeDebugTrajectoryFilter");

