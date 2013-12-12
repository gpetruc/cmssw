#include "RecoTracker/HTPattern/interface/SimplisticTrajectoryBuilder.h"
#include "RecoTracker/HTPattern/interface/TrackCandidateBuilderFromCluster.h"
#include "TrackingTools/KalmanUpdators/interface/KFUpdator.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateAccessor.h"
#include "RecoTracker/Record/interface/CkfComponentsRecord.h"
#include "TrackingTools/TrajectoryFiltering/interface/TrajectoryFilter.h"
#include "RecoTracker/CkfPattern/interface/GroupedTrajCandLess.h"
#include "TrackingTools/MeasurementDet/interface/LayerMeasurements.h"
#include "RecoTracker/HTPattern/src/TrajectorySegmentBuilder.h"
#include "RecoTracker/HTPattern/interface/HTDebugger.h"

SimplisticTrajectoryBuilder::SimplisticTrajectoryBuilder(const edm::ParameterSet &iConfig) :
    propagatorLabel_(iConfig.getParameter<std::string>("propagator")),
    propagatorOppositeLabel_(iConfig.getParameter<std::string>("propagatorOpposite")),
    hitBuilderLabel_(iConfig.getParameter<std::string>("TTRHBuilder")),
    estimatorLabel_(iConfig.getParameter<std::string>("chi2MeasurementEstimator")),
    estimatorRebuildLabel_(iConfig.getParameter<std::string>("chi2MeasurementEstimatorForRebuild")),
    trajectoryCleanerLabel_(iConfig.getParameter<std::string>("trajectoryCleaner")),
    trajectoryFilterLabel_(iConfig.getParameter<std::string>("trajectoryFilter")),
    trajectoryFilterStartLabel_(iConfig.getParameter<std::string>("trajectoryFilterStart")),
    trajectoryFilterRebuildLabel_(iConfig.getParameter<std::string>("trajectoryFilterRebuild")),
    foundHitBonus_(iConfig.getParameter<double>("foundHitBonus")),
    lostHitPenalty_(iConfig.getParameter<double>("lostHitPenalty")),
    searchHits_(str2policy(iConfig.getParameter<std::string>("searchHits"))),
    useGrouped_(str2policy(iConfig.getParameter<std::string>("useGrouped")))
{
}

SimplisticTrajectoryBuilder::~SimplisticTrajectoryBuilder()
{
}

void
SimplisticTrajectoryBuilder::setEvent(const MeasurementTrackerEvent &iEvent, const edm::EventSetup &iSetup) 
{
    mtEvent_ = & iEvent;

    // get tracker
    iSetup.get<TrackerDigiGeometryRecord>().get(tracker_);
    // get propagator
    iSetup.get<TrackingComponentsRecord>().get(propagatorLabel_, propagator_);
    iSetup.get<TrackingComponentsRecord>().get(propagatorOppositeLabel_, propagatorOpposite_);
    // get 
    iSetup.get<TransientRecHitRecord>().get(hitBuilderLabel_,hitBuilder_);
    
    iSetup.get<IdealMagneticFieldRecord>().get(bfield_);
    iSetup.get<TrackingComponentsRecord>().get(estimatorLabel_,estimator_);  
    iSetup.get<TrackingComponentsRecord>().get(estimatorRebuildLabel_,estimatorRebuild_);  

    iSetup.get<TrajectoryCleaner::Record>().get(trajectoryCleanerLabel_, trajectoryCleaner_);

    iSetup.get<TrajectoryFilter::Record>().get(trajectoryFilterLabel_, trajectoryFilter_);
    iSetup.get<TrajectoryFilter::Record>().get(trajectoryFilterRebuildLabel_, trajectoryFilterRebuild_);
    iSetup.get<TrajectoryFilter::Record>().get(trajectoryFilterStartLabel_, trajectoryFilterStart_);

    iSetup.get<IdealGeometryRecord>().get(tTopo_);
}

bool
SimplisticTrajectoryBuilder::run(const std::vector<const TrackingRecHit *> & hits, const TrajectoryStateOnSurface & stateOnFirstHit, PropagationDirection startingDirection, std::vector<Trajectory> &out) const
{
    if (DEBUG >= 2) targetTrack_ = HTDebugger::matchCandidate(hits);

    TempTrajectory traj(startingDirection);
    fitStartingHits(hits, stateOnFirstHit, traj);
    if (!trajectoryFilterStart_->qualityFilter(traj)) {
        DEBUG3_printf("\tTrajectory fails the post-seeding filter.\n");
        return false;
    }
    
    DEBUG3_printf("\t\nNow doing pattern reco.\n");
    continueTrajectory(traj);
    if (!trajectoryFilter_->qualityFilter(traj)) {
        DEBUG3_printf("\tTrajectory fails the filter after forward pattern reco.\n");
        return false;
    }

    DEBUG3_printf("\t\nNow doing backwards pattern reco.\n");
    TempTrajectory rebuilt;
    rebuildTrajectory(traj, rebuilt);
    if (!trajectoryFilterRebuild_->qualityFilter(rebuilt)) {
        DEBUG3_printf("\tTrajectory fails the filter after backwards pattern reco.\n");
        return false;
    }

    out.push_back(rebuilt.toTrajectory());
    // reverse, as it should match the seed direction
    out.back().reverse();
    // put in a dummy but valid seed
    edm::OwnVector<TrackingRecHit> seedHits;  seedHits.push_back(*hits.front());
    PTrajectoryStateOnDet const & PTraj = trajectoryStateTransform::persistentState(stateOnFirstHit, hits.front()->geographicalId().rawId());
    out.back().setSharedSeed(boost::shared_ptr<TrajectorySeed>(new TrajectorySeed(PTraj,std::move(seedHits),startingDirection)));
    return true;
}

void
SimplisticTrajectoryBuilder::fitStartingHits(const std::vector<const TrackingRecHit *> & hits, const TrajectoryStateOnSurface & stateOnFirstHit, TempTrajectory & traj) const
{
    const Propagator * propagator = (traj.direction() == alongMomentum ? propagator_ : propagatorOpposite_).product();
    TrajectoryStateOnSurface tsos = stateOnFirstHit;
    unsigned nhits = hits.size(), myskip = 0, myhits = 0, mylost = 0, mypixelhits = 0;
    bool finalfilter = false;
    const TrackingRecHit* lasthit = 0;
    KFUpdator updator;
    for (unsigned int ihit = 0; ihit < nhits; ++ihit) {
        const TrackingRecHit* hit = hits[ihit];
        bool firstOfMatch  = (ihit < nhits-1 && abs(hits[ihit+1]->geographicalId().rawId()-hit->geographicalId().rawId()) == 1);
        TrajectoryStateOnSurface state = (ihit > 0 ? propagator->propagate(tsos, tracker_->idToDet(hit->geographicalId())->surface()) : tsos);
        if (!state.isValid()) { DEBUG3_printf ("\tfailed propagation\n");
            if (tsos.isValid()) { 
                DEBUG4_printf("\t\t\tstarting state: rho = %9.4f, z = %+9.4f +/- %7.4f, phi = %+5.3f, eta = %+5.3f\n",
                        tsos.globalPosition().perp(), tsos.globalPosition().z(), std::sqrt(tsos.cartesianError().matrix()(2,2)), float(tsos.globalPosition().phi()), tsos.globalPosition().eta());
                TransientTrackingRecHit::RecHitPointer tth = hitBuilder_->build(hit);
                DEBUG4_printf("\t\t\trechit prefit : rho = %9.4f, z = %+9.4f +/- %7.4f, phi = %+5.3f, eta = %+5.3f (detid %10d)\n",
                        tth->globalPosition().perp(), tth->globalPosition().z(), std::sqrt(tth->globalPositionError().czz()), float(tth->globalPosition().phi()), tth->globalPosition().eta(), tth->geographicalId().rawId());
            }
            // for matched pairs, we mark the failure only if both fail 
            if (!firstOfMatch && (lasthit == 0 || abs(lasthit->geographicalId().rawId()-hit->geographicalId().rawId()) > 1)) { 
                myskip++; 
                mylost++; 
            }
            if (myskip >= 2) break;
            continue; 
        }
        DEBUG4_printf("\t\t\tpropag.  state: rho = %6.2f, z = %+7.2f +/- %5.2f, phi = %+5.3f, eta = %+5.3f, u = %+7.3f +/- %6.3f, v = %+7.3f +/- %6.3f\n",
                state.globalPosition().perp(), state.globalPosition().z(), std::sqrt(state.cartesianError().matrix()(2,2)), float(state.globalPosition().phi()), state.globalPosition().eta(),
                state.localParameters().position().x(), std::sqrt(state.localError().matrix()(3,3)), state.localParameters().position().y(), std::sqrt(state.localError().matrix()(4,4)));
        TransientTrackingRecHit::RecHitPointer tth = hitBuilder_->build(hit);
        DEBUG4_printf("\t\t\trechit prefit : rho = %6.2f, z = %+7.2f +/- %5.2f, phi = %+5.3f, eta = %+5.3f, u = %+7.3f +/- %6.3f, v = %+7.3f +/- %6.3f\n",
                tth->globalPosition().perp(), tth->globalPosition().z(), std::sqrt(tth->globalPositionError().czz()), float(tth->globalPosition().phi()), tth->globalPosition().eta(),
                tth->localPosition().x(), std::sqrt(tth->parametersError()(1,1)), tth->localPosition().y(), std::sqrt(tth->parametersError()(2,2)));
        tth = tth->clone(state);
        DEBUG4_printf("\t\t\trechit postfit: rho = %6.2f, z = %+7.2f +/- %5.2f, phi = %+5.3f, eta = %+5.3f, u = %+7.3f +/- %6.3f, v = %+7.3f +/- %6.3f\n",
                tth->globalPosition().perp(), tth->globalPosition().z(), std::sqrt(tth->globalPositionError().czz()), float(tth->globalPosition().phi()), tth->globalPosition().eta(),
                tth->localPosition().x(), std::sqrt(tth->parametersError()(1,1)), tth->localPosition().y(), std::sqrt(tth->parametersError()(2,2)));
        std::pair<bool,double> chi2 = estimator_->estimate(state, *tth);
        double chi2cut = estimator_->chiSquaredCut();
        if (myskip == 1) chi2cut *= 0.6;
        if (mylost > 0 ) chi2cut *= 0.6;
        if (finalfilter) chi2cut *= 0.6;
        if (chi2.second > chi2cut) {
            // for matched pairs, we mark the failure only if both fail 
            if (!firstOfMatch && (lasthit == 0 || abs(lasthit->geographicalId().rawId()-hit->geographicalId().rawId()) > 2)) { 
                myskip++; 
                mylost++; 
            }
            DEBUG3_printf("\tSkipping hit on detid %10d/%7d due to bad chi2 = %7.1f (threshold %6.1f, %d consec. skipped hits)\n", hit->geographicalId().rawId(), TrackCandidateBuilderFromCluster::hitid(hit), chi2.second, chi2cut, myskip);
            if (myskip >= 2) break;
            continue;
       } else {
            myskip = 0;
        }
        TrajectoryStateOnSurface updated = updator.update(state, *tth);
        if (!updated.isValid()) { DEBUG3_printf("\tfailed update\n"); break; } 
        traj.push(TrajectoryMeasurement(state, updated, tth, chi2.second));
        tsos = updated; lasthit = hit; 
        myhits++; if (hit->geographicalId().subdetId() <= 2) mypixelhits++;
        DEBUG3_printf("\tKF-fitted %2d hits (last chi2 = %7.1f; pt = %7.2f +/- %7.2f, q = %+1d; last detid %10d/%7d, matched %d)\n", myhits, chi2.second, tsos.globalMomentum().perp(), TrajectoryStateAccessor(*tsos.freeState()).inversePtError()*tsos.globalMomentum().perp2(), tsos.charge(),lasthit->geographicalId().rawId(),TrackCandidateBuilderFromCluster::hitid(lasthit), HTDebugger::isMatchedToTrack(hit, targetTrack_));
        DEBUG4_printf("\t\t local q/p = %+8.4f +/- %8.4f, dxdz = %+8.4f +/- %8.4f, dydz = %+8.4f +/- %8.4f, x = %+8.5f +/- %8.5f, y = %+8.5f +/- %8.5f\n",
                            tsos.localParameters().qbp(),  std::sqrt(tsos.localError().matrix()(0,0)), tsos.localParameters().dxdz(), std::sqrt(tsos.localError().matrix()(1,1)), tsos.localParameters().dydz(), std::sqrt(tsos.localError().matrix()(2,2)),
                            tsos.localParameters().position().x(), std::sqrt(tsos.localError().matrix()(3,3)), tsos.localParameters().position().y(), std::sqrt(tsos.localError().matrix()(4,4)));
        if (!trajectoryFilterStart_->toBeContinued(traj)) {
            DEBUG3_printf("\tTrajectory filter requests to stop the trajectory here.\n");
            break;
        }
        if (!finalfilter) {
            // we want to know if we already qualify for the end, so that we can tighten the chi2 cut
            finalfilter = trajectoryFilterStart_->qualityFilter(traj);
        }
    }
}

void 
SimplisticTrajectoryBuilder::continueTrajectory(TempTrajectory & traj) const
{
    do {
        bool hasnext = nextLayer(traj, *estimator_, 0);
        if (!hasnext) break;
    } while (trajectoryFilter_->toBeContinued(traj));

    while (!traj.empty() && !traj.lastMeasurement().recHit()->isValid()) {
        traj.pop();
    }
}

void
SimplisticTrajectoryBuilder::rebuildTrajectory(const TempTrajectory &src, TempTrajectory & out) const
{
    out = TempTrajectory(src.direction() == alongMomentum ? oppositeToMomentum : alongMomentum);
    const TrajectoryMeasurement & last = src.lastMeasurement();
    TrajectoryStateOnSurface tsos = last.updatedState();
    tsos.rescaleError(10.);
    TrajectoryStateOnSurface upd = KFUpdator().update(tsos, *last.recHit());
    if (!tsos.isValid()) { DEBUG1_printf("ERROR updating back the last hit?\n"); return; }
    out.push(TrajectoryMeasurement(tsos, upd, last.recHit(), 0.)); 
    TempTrajectory::DataContainer::iterator ilast = src.measurements().begin(), ifirst = src.measurements().end(); 
    do {
        const TransientTrackingRecHit *hint = 0;
        if (out.lastMeasurement().recHit()->isValid() && ilast != ifirst) {
            TempTrajectory::DataContainer::iterator it = ilast; 
            for (++it; it != ifirst; ++it) {
                if (!it->recHit()->isValid()) continue;
                if (it->recHit()->geographicalId() == out.lastMeasurement().recHit()->geographicalId()) {
                    ilast = it;
                    ++it;
                    if (it != ifirst && it->recHit()->isValid()) {
                        hint = &*it->recHit();
                    }
                    break;
                }
            }
        }
        bool hasnext = nextLayer(out, *estimatorRebuild_, hint);
        if (!hasnext) break;
    } while (trajectoryFilterRebuild_->toBeContinued(out));

    while (!out.empty() && !out.lastMeasurement().recHit()->isValid()) {
        out.pop();
    }
}

bool
SimplisticTrajectoryBuilder::nextLayer(TempTrajectory & traj, const Chi2MeasurementEstimatorBase &est, const TransientTrackingRecHit *hint) const
{
    TrajectoryStateOnSurface startingState = traj.lastMeasurement().updatedState();
    const Propagator * propagator = (traj.direction() == alongMomentum ? propagator_ : propagatorOpposite_).product();
    const DetLayer *startLayer = traj.lastLayer();
    if (startLayer == 0) { 
        startLayer = mtEvent_->geometricSearchTracker()->idToLayer(traj.lastMeasurement().recHit()->geographicalId());
        if (startLayer == 0) { 
            DEBUG1_printf("Missing layer!\n"); 
            return false; 
        }
    }
    std::vector<const DetLayer*> nextLayers = startLayer->nextLayers( *startingState.freeState(), traj.direction() );
    DEBUG2_printf("\tStarting from layer %3d (subdetector %2d), found %ld compatible layers\n", startLayer->seqNum(), startLayer->subDetector(), nextLayers.size()); 
    LayerMeasurements layerMeasurements(mtEvent_->measurementTracker(), *mtEvent_);
    KFUpdator updator;
    bool foundSegments = false;
    std::vector<TempTrajectory> nextValid; std::vector<TempTrajectory> nextInvalid;
    const DetLayer *hintLayer = (hint != 0 && hint->isValid() ? mtEvent_->geometricSearchTracker()->idToLayer(hint->geographicalId()) : 0);
    if (hint != 0 && hint->isValid() && hintLayer != 0) {
        DEBUG2_printf("\t\thint for hit %10d/%7d on layer %3d\n", hint->geographicalId().rawId(), TrackCandidateBuilderFromCluster::hitid(hint->hit()), hintLayer->seqNum());
    }
    if (nextLayers.size() == 1 && nextLayers.front() == hintLayer && 
            (searchHits_ == noLayer || (searchHits_ == pixelLayer && hintLayer->subDetector() <= 2))) {
        DEBUG2_printf("\t\t layer matches with hint (detid %10d)\n", hint->geographicalId().rawId());
        TrajectoryStateOnSurface prop = propagator->propagate(startingState, tracker_->idToDet(hint->geographicalId())->surface());
        if (prop.isValid()) {
            TransientTrackingRecHit::RecHitPointer hintrefit = hint->clone(prop);
            std::pair<bool,double> chi2 = est.estimate( prop, *hintrefit );
            TrajectoryStateOnSurface updated = updator.update(prop, *hintrefit);
            if (chi2.first) {
                traj.push(TrajectoryMeasurement(prop, updated, hintrefit, chi2.second));
                DEBUG2_printf("\t\t using rechit from hint: %10d/%7d, chi2 %7.1f\n", hint->geographicalId().rawId(), TrackCandidateBuilderFromCluster::hitid(hint->hit()), chi2.second);
                return true;
            }
        }
    }

    for (const DetLayer * newLayer : nextLayers) {
        bool grouped = (useGrouped_ == anyLayer);
        if (useGrouped_ == pixelLayer && newLayer->subDetector() <= 2) grouped = true;

        if (grouped) {
            TrajectorySegmentBuilder layerBuilder(&layerMeasurements, *newLayer, *propagator, updator, est, true, true);
            std::vector<TempTrajectory> segments = layerBuilder.segments(startingState);

            DEBUG2_printf("\t => going to layer %3d (subdetector %2d), number of segments = %ld\n", newLayer->seqNum(), newLayer->subDetector(), segments.size()); 

            if ( !segments.empty() )  foundSegments = true;
            for ( std::vector<TempTrajectory>::iterator is=segments.begin();
                    is!=segments.end(); is++ ) {

                // assume "invalid hit only" segment is last in list
                const TempTrajectory::DataContainer & measurements = is->measurements();
                if ( measurements.size()==1 && (measurements.front().recHit()->getType() == TrackingRecHit::missing) ) {
                    if (!nextValid.empty()) continue;
                    if (nextInvalid.empty()) {
                        DetId where = measurements.front().recHit()->geographicalId();
                        nextInvalid.push_back(traj);
                        nextInvalid.back().push(measurements.front());
                        DEBUG2_printf("\t\tinvalid hit will be added on subdet %d, layer %d (if needed)\n", where.subdetId(), tTopo_->layer(where));
                        //DEBUG2_printf("\t\tinvalid hit will be added on subdet %d, layer %d", measurements.front().recHit()->geographicalId().subdetId(), tTopo_->layer(measurements.front().recHit()->geographicalId()));
                    }
                    continue;
                }
                else if (DEBUG >= 3) {
                    DEBUG3_printf("\t\tsegment with %d hits\n", measurements.size());
                    for (const auto &tm : measurements) {
                        const TrackingRecHit *hit = tm.recHit()->hit();
                        DEBUG3_printf("\t\t\thit with chi2 = %7.1f on detid %10d/%7d, matched %d (%s)\n", tm.estimate(), hit->geographicalId().rawId(),TrackCandidateBuilderFromCluster::hitid(hit), HTDebugger::isMatchedToTrack(hit, targetTrack_), typeid(*hit).name());
                    }
                }
                if (is->foundHits() == 0) {
                    DEBUG3_printf("\t\tthis segment has no valid hits.\n");
                    if (!nextInvalid.empty() && nextInvalid.back().lastMeasurement().recHit()->getType() == TrackingRecHit::missing) {
                        DEBUG3_printf("\t\tclearning existing vector of invalid hits.\n");
                        nextInvalid.clear();
                    }
                    if (nextInvalid.empty()) {
                        DetId where = measurements.front().recHit()->geographicalId();
                        nextInvalid.push_back(traj);
                        nextInvalid.back().push(measurements.front());
                        DEBUG2_printf("\t\tinvalid hit will be added on subdet %d, layer %d (if needed)\n", where.subdetId(), where.rawId() ? tTopo_->layer(where) : 0);
                    }
                    continue; 
                }

                //----  avoid to add the same hits more than once in the trajectory ----
                bool toBeRejected(false);
                for (const TempTrajectory::DataContainer::const_iterator revIt = measurements.rbegin(); revIt!=measurements.rend(); --revIt) {
                    for (const TempTrajectory::DataContainer::const_iterator newTrajMeasIt = traj.measurements().rbegin(); 
                            newTrajMeasIt != traj.measurements().rend(); --newTrajMeasIt) {
                        if (revIt->recHitR().geographicalId()==newTrajMeasIt->recHitR().geographicalId() 
                                && (revIt->recHitR().geographicalId() != DetId(0)) ){
                            toBeRejected = true;
                            goto rejected; //break;  // see http://stackoverflow.com/questions/1257744/can-i-use-break-to-exit-multiple-nested-for-loops
                        }
                    }
                }
                rejected:;    // http://xkcd.com/292/
                if (toBeRejected) { 
                    continue; //Are we sure about this????
                }
                //------------------------

                // create new candidate
                TempTrajectory newTraj(traj);
                newTraj.join(*is);
                nextValid.push_back(newTraj);
            }
        } else {
            std::vector<TrajectoryMeasurement> meas = layerMeasurements.measurements( *newLayer, startingState, *propagator, est );
            foundSegments = true;
            for (const TrajectoryMeasurement & m : meas) {
                if (m.recHit()->isValid()) {
                    const TrackingRecHit *hit = m.recHit()->hit();
                    DEBUG3_printf("\t\t\thit with chi2 = %7.1f on detid %10d/%7d, matched %d (%s)\n", m.estimate(), hit->geographicalId().rawId(),TrackCandidateBuilderFromCluster::hitid(hit), HTDebugger::isMatchedToTrack(hit, targetTrack_), typeid(*hit).name());
                    TempTrajectory newTraj(traj);
                    newTraj.push(m);
                    nextValid.push_back(newTraj);
                } else if (nextInvalid.empty() || (nextInvalid.back().lastMeasurement().recHit()->getType() == TrackingRecHit::missing && 
                                                   m.recHit()->getType() != TrackingRecHit::missing))  {
                    nextInvalid.clear();
                    TempTrajectory newTraj(traj);
                    newTraj.push(m);
                    nextValid.push_back(newTraj);
                }
                break;  // they are already sorted
            }
        }
    }
    if (foundSegments) {
        if (nextValid.empty()) {
            traj = nextInvalid.front();
        } else {
            std::sort(nextValid.begin(), nextValid.end(), GroupedTrajCandLess(foundHitBonus_, lostHitPenalty_));
            traj = nextValid.front();
        }
    }
    return foundSegments;
}
