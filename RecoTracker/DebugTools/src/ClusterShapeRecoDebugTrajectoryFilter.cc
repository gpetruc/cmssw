#include "RecoTracker/DebugTools/interface/ClusterShapeRecoDebugTrajectoryFilter.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

ClusterShapeRecoDebugTrajectoryFilter::ClusterShapeRecoDebugTrajectoryFilter(const edm::ParameterSet &iConfig, edm::ConsumesCollector& iC) :
    ClusterShapeDebugTrajectoryFilter(iConfig,iC),
    theTracks_(iC.consumes<std::vector<reco::Track>>(theConfig.getParameter<edm::InputTag>("referenceTracks")))
{
}

void ClusterShapeRecoDebugTrajectoryFilter::fillAssociations(const TrackingRecHit *hit, std::vector<Id2> &out) const
{
    out.clear();
    std::vector<TrackMixingAssociator::TrackAssociation> assoc;
    theRecoAssociator->associateToTracks(*hit, assoc);
    for (const auto & a : assoc) {
        out.push_back(Id2(a.eventId, a.track - &theTracksHandle_->front()));
    }
}

void ClusterShapeRecoDebugTrajectoryFilter::initAssociator(const edm::Event &e, const edm::EventSetup &) 
{
    e.getByToken(theTracks_, theTracksHandle_);
    theRecoAssociator->registerTrackEvent(1,*theTracksHandle_);
}
