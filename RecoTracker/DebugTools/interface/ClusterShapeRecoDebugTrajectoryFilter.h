#ifndef _ClusterShapeRecoDebugTrajectoryFilter_h_
#define _ClusterShapeRecoDebugTrajectoryFilter_h_

#include "RecoTracker/DebugTools/interface/ClusterShapeDebugTrajectoryFilter.h"
#include "RecoTracker/DebugTools/interface/TrackMixingAssociator.h"

class ClusterShapeRecoDebugTrajectoryFilter : public ClusterShapeDebugTrajectoryFilter {
    public:
        ClusterShapeRecoDebugTrajectoryFilter(const edm::ParameterSet &iConfig, edm::ConsumesCollector& iC);

    private:
        virtual void fillAssociations(const TrackingRecHit *hit, std::vector<Id2> &out) const  override;
        virtual void initAssociator(const edm::Event &, const edm::EventSetup &) override;

        edm::EDGetTokenT<std::vector<reco::Track>> theTracks_;
        edm::Handle<std::vector<reco::Track>>  theTracksHandle_;
        std::unique_ptr<TrackMixingAssociator> theRecoAssociator_;
};

#endif

