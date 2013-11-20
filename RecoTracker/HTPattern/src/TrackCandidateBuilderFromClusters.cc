#include "RecoTracker/HTPattern/interface/TrackCandidateBuilderFromCluster.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "TrackingTools/KalmanUpdators/interface/KFUpdator.h"
#include "RecoTracker/TkSeedGenerator/interface/FastHelix.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h" 
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateAccessor.h"
#include "RecoTracker/Record/interface/CkfComponentsRecord.h"
#include "FWCore/Utilities/interface/isFinite.h"
#include "DataFormats/TrackerRecHit2D/interface/TrackerSingleRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Math/GenVector/CylindricalEta3D.h"
#include "RecoTracker/HTPattern/interface/HTDebugger.h"
#include "TrackingTools/TransientTrackingRecHit/interface/InvalidTransientRecHit.h"

#ifdef  HTDEBUG
#define DEBUG HTDEBUG
#else
#define DEBUG 0.5
#endif
#define DEBUG1_printf  if (DEBUG>=1) printf
#define DEBUG2_printf  if (DEBUG>=2) printf
#define DEBUG3_printf  if (DEBUG>=3) printf
#define DEBUG4_printf  if (DEBUG>=4) printf

namespace {
    unsigned int mask2layer(unsigned int layermask) {
        unsigned int layer = 0;
        while ((layermask & 1) == 0) { layermask >>= 1; layer++; }
        return layer;
    }
    unsigned int countbits(unsigned int layermask) {
        unsigned int layers = 0;
        while (layermask != 0) { layers += (layermask & 1); layermask >>= 1; }
        return layers;
    }
};

TrackCandidateBuilderFromCluster::TrackCandidateBuilderFromCluster(const edm::ParameterSet &iConfig) :
    propagatorLabel_(iConfig.getParameter<std::string>("propagator")),
    propagatorOppositeLabel_(iConfig.getParameter<std::string>("propagatorOpposite")),
    hitBuilderLabel_(iConfig.getParameter<std::string>("TTRHBuilder")),
    estimatorLabel_(iConfig.getParameter<std::string>("chi2MeasurementEstimator")),
    navigationSchoolLabel_(iConfig.getParameter<std::string>("NavigationSchool")),
    trajectoryBuilderLabel_(iConfig.getParameter<std::string>("TrajectoryBuilder")),
    trajectoryCleanerLabel_(iConfig.getParameter<std::string>("TrajectoryCleaner")),
    cleanTrajectoryAfterInOut_(iConfig.getParameter<bool>("cleanTrajectoryAfterInOut")),
    initialStateEstimatorPSet_(iConfig.getParameter<edm::ParameterSet>("transientInitialStateEstimatorParameters")),
    useHitsSplitting_(iConfig.getParameter<bool>("useHitsSplitting")),
    splitSeedHits_(iConfig.getParameter<bool>("splitSeedHits")),
    dphiCutPair_(iConfig.getParameter<double>("dphiCutPair")),
    detaCutPair_(iConfig.getParameter<double>("detaCutPair")),
    pairsOnSeedCellOnly_(iConfig.getParameter<bool>("pairsOnSeedCellOnly")),
    minHits_(iConfig.getParameter<uint32_t>("minHits")),
    minHitsIfNoCkf_(iConfig.getParameter<uint32_t>("minHitsIfNoCkf")),
    maxFailMicroClusters_(iConfig.getParameter<uint32_t>("maxFailMicroClusters")),
    startingCovariance_(iConfig.getParameter<std::vector<double>>("startingCovariance")),
    initialStateEstimator_(0)
{
    std::vector<double> dphiCutHits = iConfig.getParameter<std::vector<double>>("dphiCutHits");
    std::vector<double> detaCutHits = iConfig.getParameter<std::vector<double>>("detaCutHits");
    dphiCutHits_[0] = dphiCutHits[0]; dphiCutHits_[1] = dphiCutHits[1]; 
    detaCutHits_[0] = detaCutHits[0]; detaCutHits_[1] = detaCutHits[1]; 
}

void 
TrackCandidateBuilderFromCluster::init(const edm::EventSetup& es, const MeasurementTrackerEvent &evt, const SeedComparitor *ifilter)
{
    filter = ifilter;
    // get tracker
    es.get<TrackerDigiGeometryRecord>().get(tracker_);
    // get propagator
    es.get<TrackingComponentsRecord>().get(propagatorLabel_, propagator_);
    es.get<TrackingComponentsRecord>().get(propagatorOppositeLabel_, propagatorOpposite_);
    es.get<TrackingComponentsRecord>().get("AnyDirectionAnalyticalPropagator", anyPropagator_);
    //
    es.get<TransientRecHitRecord>().get(hitBuilderLabel_,hitBuilder_);
    
    es.get<IdealMagneticFieldRecord>().get(bfield_);
    es.get<TrackingComponentsRecord>().get(estimatorLabel_,estimator_);  

    es.get<TrajectoryCleaner::Record>().get(trajectoryCleanerLabel_, trajectoryCleaner_);

    es.get<CkfComponentsRecord>().get(trajectoryBuilderLabel_, trajectoryBuilderTemplate_);    
    const BaseCkfTrajectoryBuilder *builder = dynamic_cast<const BaseCkfTrajectoryBuilder *>(trajectoryBuilderTemplate_.product());
    assert(builder != 0);
    trajectoryBuilder_.reset(builder->clone(&evt)); 
    //maskStripClusters_ = & evt.stripClustersToSkip();
    //maskPixelClusters_ = & evt.pixelClustersToSkip();

    es.get<NavigationSchoolRecord>().get(navigationSchoolLabel_, navigationSchool_);
    navigationSetter_.reset(new NavigationSetter(*navigationSchool_));

    if (initialStateEstimator_ == 0) {
        initialStateEstimator_ = new TransientInitialStateEstimator(es, initialStateEstimatorPSet_);
    }
    initialStateEstimator_->setEventSetup(es);
}

void 
TrackCandidateBuilderFromCluster::done() 
{
    trajectoryBuilder_.reset();
    navigationSetter_.reset();
}

void
TrackCandidateBuilderFromCluster::run(const HTCluster &cluster, unsigned int nseedlayersCut, unsigned int nlayersCut, TrackCandidateCollection & tcCollection, TrajectorySeedCollection & seedCollection, TrajectorySeedCollection *seedsFromAllClusters) 
{
    std::vector<ClusteredHit> hits; hits.reserve(50);
    HitIndex hitindex;
    std::vector<Trajectory> allTrajectories; allTrajectories.reserve(5);
    const HTCell &seed = mapHiRes_->get(cluster.ieta(), cluster.iphi());
    DEBUG1_printf("\n --- Cluster #%03d at eta %+5.3f, phi %+5.3f with %d seed layers, %d hi-res layers, %d layers; front detid %d --- \n", 
                        seed.icluster(), hitsHiRes_->eta(seed.hits().front()), hitsHiRes_->phi(seed.hits().front()),
                        cluster.nseedlayers(), cluster.nlayers(), cluster.nmorelayers(),
                        hitsHiRes_->hit(seed.hits().front())->geographicalId().rawId());
    unsigned int seedlayermask = 0, layermask = 0;
    for (int hit : seed.hits()) {
        if ((*maskHiRes_)[hit]) continue;
        seedlayermask |= hitsHiRes_->layermask(hit);
        hits.push_back(ClusteredHit(*hitsHiRes_, hit, mask2layer(hitsHiRes_->layermask(hit))));
    }
    DEBUG2_printf("\tpost-mask seed layers: %d\n", countbits(seedlayermask));
    if (seedlayermask == 0) return; // nothing to do
    if (countbits(seedlayermask) < nseedlayersCut) return;
    if (DEBUG >= 2) {
        DEBUG2_printf("tracks associated to seed cell hits only:\n");
        TrajectorySeedCollection tmp;
        dumpAsSeed(hits, tmp);
        HTDebugger::printAssociatedTracks(tmp.back());
    }
    unsigned int nseedhits = hits.size();
    layermask = seedlayermask;
    std::array<const HTCell *, 8> cells; unsigned int ncells;
    mapHiRes_->getNeighbours(cluster.ieta(), cluster.iphi(), cells, ncells);
    for (unsigned int ic = 0; ic < ncells; ++ic) {
        for (int hit : cells[ic]->hits()) {
            if ((*maskHiRes_)[hit]) continue;
            layermask |= hitsHiRes_->layermask(hit);
            hits.push_back(ClusteredHit(*hitsHiRes_, hit, mask2layer(hitsHiRes_->layermask(hit))));
        }
    }
    DEBUG2_printf("\tpost-mask hi-res layers: %d\n", countbits(layermask));
    unsigned int nhireshits = hits.size();
    if (useLowRes_ && hitsLowRes_->size()) {
        const HTCell &lowres = mapLowRes_->get(cluster.ieta()>>etashift_, cluster.iphi()>>phishift_);
        for (int hit : lowres.hits()) {
            if ((*maskLowRes_)[hit]) continue;
            int il = mask2layer(hitsLowRes_->layermask(hit)); if (il < 3) continue;
            const TrackingRecHit * thishit = hitsLowRes_->hit(hit);
            bool found = false;
            for (const ClusteredHit &other : hits) {
                if (il == other.layer && other.hit->geographicalId() == thishit->geographicalId() && other.hit->sharesInput(thishit, TrackingRecHit::some)) {
                    found = true; break;
                }
            }
            if (!found) {
                hits.push_back(ClusteredHit(*hitsLowRes_, hit, mask2layer(hitsLowRes_->layermask(hit)), 0.2));
                layermask |= hitsLowRes_->layermask(hit);
            }
        }
    }
    DEBUG2_printf("\tpost-mask total layers: %d\n", countbits(layermask));
    if (countbits(layermask) < nlayersCut) return;
    unsigned int nhits = hits.size();
    DEBUG2_printf("found %d seed hits, %d hi-res hits, %d hits\n", nseedhits, nhireshits, nhits);
    if (seedsFromAllClusters) if (DEBUG>0)  dumpAsSeed(hits, *seedsFromAllClusters);
    if (seedsFromAllClusters) if (DEBUG>=1) HTDebugger::printAssociatedTracks(seedsFromAllClusters->back());
    if (seedsFromAllClusters) if (DEBUG>=2) HTDebugger::beginLoggingCluster(seedsFromAllClusters->back(), nseedhits, nhireshits, nhits, hitsHiRes_->alpha());

    typedef std::pair<uint16_t,uint16_t> hitpair;
    typedef std::pair<float, hitpair> distpair;
    std::vector<distpair> pairs;
    for (unsigned int i1 = 0; i1 < nseedhits; ++i1) {
        for (unsigned int i2 = i1+1; i2 < nhits; ++i2) {
            if (pairsOnSeedCellOnly_ && i2 >= nseedhits) break;
            if (hits[i2].layer == hits[i1].layer) continue; // no same-layer pairs in seeding
            if (std::abs(hits[i1].layer - hits[i2].layer) <= 2) { // to too-far-away-pairs
                float dphi = hits[i1].dphi(hits[i2]);
                float deta = hits[i1].deta(hits[i2]);
                //printf("\t pair %3d, %3d: (%+5.3f,%+5.3f, ly %2d)  (%+5.3f,%+5.3f ly %2d): d2d %.4f\n", i1,i2, hits[i1].eta,hits[i1].phi,hits[i1].layer, hits[i2].eta,hits[i2].phi,hits[i2].layer, d2d);
                if (dphi < dphiCutPair_ && deta < detaCutPair_) pairs.push_back(distpair(deta, hitpair(i1,i2)));
            }
        }
    }
    DEBUG2_printf("\t pairs %d\n", int(pairs.size()));
    std::sort(pairs.begin(), pairs.end());
    std::vector<uint8_t> microcluster(nhits, 0);
    uint8_t microclusters = 0;
    std::array<int, 20>   layerhits;   
    std::array<float, 20> layerdist;
    unsigned int failMicroClusters = 0;
    for (unsigned int ip = 0, np = pairs.size(); ip < np && microclusters < 100; ++ip) {
        unsigned int i1 = pairs[ip].second.first; 
        unsigned int i2 = pairs[ip].second.second; 
        if (microcluster[i1] != 0 || microcluster[i2] != 0) {
            continue;
        }
        ++microclusters;
        if (++failMicroClusters > maxFailMicroClusters_) { DEBUG2_printf("Reached limit of consecutive failing microclusters.\n"); break; }
        DEBUG2_printf("starting microcluster %d with pair %3d, %3d: (%+5.3f,%+5.3f, ly %2d, detid %10d/%7d)  (%+5.3f,%+5.3f ly %2d, detid %10d/%7d), deta = %5.3f, dr = %5.3f\n", 
                            microclusters, i1,i2, hits[i1].eta,hits[i1].phi,hits[i1].layer,hits[i1].hit->geographicalId().rawId(),hitid(hits[i1].hit), 
                            hits[i2].eta,hits[i2].phi,hits[i2].layer,hits[i2].hit->geographicalId().rawId(),hitid(hits[i2].hit), hits[i1].deta(hits[i2]), hits[i1].dist2d(hits[i2]));
        if (DEBUG>=2) HTDebugger::printAssociatedTracks(hits[i1].hit, hits[i2].hit);
        // apply pt correction from two hits
        float alphacorr = (hits[i1].phi - hits[i2].phi)/(hits[i1].rho - hits[i2].rho);
        // apply eta correction (approximate: if dz is the displacement wrt z0, then eta = eta0 + dz / (rho cosh(eta)) =: eta0 + beta / rho
        float betacorr  = (hits[i1].eta - hits[i2].eta)/(1/hits[i1].rho - 1/hits[i2].rho);
        float phi0 = hits[i1].phi - alphacorr*hits[i1].rho;
        float eta0 = hits[i1].eta - betacorr/hits[i1].rho; 
        std::array<int, 8> mchits4refit;
        mchits4refit[0] = i1; mchits4refit[1] = i2;
        if (DEBUG>=2) HTDebugger::logSeedingPair(hits[i1].hit, hits[i1].layer, i1 < nseedhits ? 0 : (i1 < nhireshits ? 1 : 2),
                       hits[i2].hit, hits[i2].layer, i2 < nseedhits ? 0 : (i2 < nhireshits ? 1 : 2),
                       pairs[ip].first, hits[i1].dphi(hits[i2]), hits[i1].deta(hits[i2]));
        unsigned int microclusternhits = 2;
        for (int iteration = 0; iteration <= 1; ++iteration) {
            DEBUG2_printf("   iteration %d: eta0,phi0 = %+5.3f,%+5.3f     alphacorr = %+8.5f (pt %8.4f), betacorr = %+8.5f\n", iteration+1, eta0, phi0, alphacorr, 0.5f*0.003f*bfieldAtOrigin_/std::max(std::abs(alphacorr-hitsHiRes_->alpha()),0.0000001f), betacorr);
            microclusternhits = 2;
            std::fill(layerhits.begin(), layerhits.end(), -1);
            std::fill(layerdist.begin(), layerdist.end(), 999.f);
            layerhits[hits[i1].layer] = i1;
            layerhits[hits[i2].layer] = i2;
            layerdist[hits[i1].layer] = 0.f; // always keep
            layerdist[hits[i2].layer] = 0.f; // the seed hits
            for (unsigned int i3 = 0; i3 < nhits; ++i3) {
                if (i3 == i1 || i3 == i2 || microcluster[i3]) continue;
                float dphi = hits[i3].dphi(phi0,alphacorr);
                float deta = hits[i3].detaE(eta0,betacorr); // include the scale factor due to the resolution of low-res hits
                float dc   = deta + dphi;
                if (DEBUG >= 4 || (dphi <= dphiCutHits_[iteration] && deta <= detaCutHits_[iteration])) {
                    DEBUG3_printf("\thit %2d %c (%+5.3f,%+5.3f, ly %2d): d2d %.4f, d2dc %.4f    deta %.4f detac %.4f detacE %.4f   dphi %.4f dphic %.4f    detid %10d/%7d %s\n",i3, 
                                  i3 < nseedhits ? 'S' : (i3 < nhireshits ? 'H' : 'L'),
                                  hits[i3].eta,hits[i3].phi,hits[i3].layer, hits[i3].dist2d(eta0,phi0,0.f,0.f), dc, hits[i3].deta(eta0,0), hits[i3].deta(eta0,betacorr), deta, hits[i3].dphi(phi0,0), hits[i3].dphi(phi0,alphacorr),
                                  hits[i3].hit->geographicalId().rawId(),hitid(hits[i3].hit),
                                  HTDebugger::thirdHitOnSameTrack(hits[i1].hit, hits[i2].hit, hits[i3].hit) ? " (good)" : " (not good)");
                    if (DEBUG>=3) HTDebugger::logThirdHit(hits[i3].hit, hits[i3].layer, i3 < nseedhits ? 0 : (i3 < nhireshits ? 1 : 2),
                                        dc, hits[i3].dphi(phi0,alphacorr), hits[i3].deta(eta0,betacorr), 
                                        iteration+1, hitsHiRes_->alpha()+alphacorr, 0);
                }
                if (dphi > dphiCutHits_[iteration] || deta > detaCutHits_[iteration]) continue;
                if (dc > layerdist[hits[i3].layer]) continue;
                if (layerhits[hits[i3].layer] == -1) {
                    microclusternhits++;
                }
                layerdist[hits[i3].layer] = dc;
                layerhits[hits[i3].layer] = i3;
            }
            if (iteration == 0) {
                if (microclusternhits == 2) break; // nothing to refit
                refitAlphaBetaCorr(hits, &layerhits[0], layerhits.size(), alphacorr, betacorr, phi0, eta0);
            }
        }

        DEBUG2_printf("microcluster %d (%d hits)\n", microclusters, microclusternhits);
        const reco::Track *match = 0;
        const char *matchstr[3] = { "(not good)", "(mixed)", "(good)" };
        if (DEBUG>=2) {
            std::vector<const TrackingRecHit *> phits;
            for (unsigned int il = 0; il < layerhits.size(); ++il) {
                if (layerhits[il] != -1) phits.push_back(hits[layerhits[il]].hit);
            } 
            match = HTDebugger::matchCandidate(phits);
        }
        for (unsigned int il = 0; il < layerhits.size(); ++il) {
            if (layerhits[il] == -1) continue;
            DEBUG2_printf("\ton layer %2d (seed: %c, prec %c): %3d, dist %7.4f deta %.4f dphi %.4f (rho %5.1f) detid %10d/%7d %s\n", il, 
                (layerhits[il] == int(i1) || layerhits[il] == int(i2) ? 'Y' : 'n'),  
                layerhits[il] < int(nseedhits) ? 'S' : (layerhits[il] < int(nhireshits) ? 'H' : 'L'),
                layerhits[il], layerdist[il], hits[layerhits[il]].detaE(eta0,betacorr), hits[layerhits[il]].dphi(phi0,alphacorr),
                hits[layerhits[il]].rho, hits[layerhits[il]].hit->geographicalId().rawId(),hitid(hits[layerhits[il]].hit),
                matchstr[HTDebugger::isMatchedToTrack(hits[layerhits[il]].hit, match)]);
                //HTDebugger::thirdHitOnSameTrack(hits[i1].hit, hits[i2].hit, hits[layerhits[il]].hit) ? " (good)" : " (not good)");
            if (layerhits[il] != int(i1) && layerhits[il] != int(i2)) {
                if (DEBUG >= 2) HTDebugger::logThirdHit(hits[layerhits[il]].hit, hits[layerhits[il]].layer, layerhits[il] < int(nseedhits) ? 0 : (layerhits[il] < int(nhireshits) ? 1 : 2),
                                        hits[layerhits[il]].dist2d(eta0,phi0,alphacorr,betacorr), hits[layerhits[il]].dphi(phi0,alphacorr), hits[layerhits[il]].deta(eta0,betacorr), 
                                        3, hitsHiRes_->alpha()+alphacorr, 1);
            }
        }
        if (microclusternhits < minHits_) continue;
        const TrajectorySeed *seed = makeSeed(hits, std::make_pair(i1,i2), layerhits, eta0,phi0,hitsHiRes_->alpha()+alphacorr, seedCollection);
        if (seed == 0) {
            seed = 0; // makeSeed2Way(hits, std::make_pair(i1,i2), layerhits, eta0,phi0,hitsHiRes_->alpha()+alphacorr, seedCollection); // FIXME
        }
        if (seed == 0) continue;
        if (DEBUG>=2) HTDebugger::printAssociatedTracks(*seed);

        std::vector<Trajectory> trajectories;
        makeTrajectories(*seed, trajectories);
        
        if (trajectories.empty()) continue;
        for (Trajectory &traj : trajectories) {
            if (!traj.isValid() || traj.measurements().empty()) continue;
            DEBUG2_printf("built trajectory with %d found hits, %d lost hits\n", traj.foundHits(), traj.lostHits());
            if (DEBUG >= 2) HTDebugger::printAssociatedTracks(traj);
            failMicroClusters = 0;
            for (const TrajectoryMeasurement &tm : traj.measurements()) {
                if (!tm.recHit()->isValid()) continue;
                int ihit = findHit(tm.recHit()->hit(), hits, hitindex);
                if (ihit == -1) continue;
                microcluster[ihit] = microclusters;
                std::vector<bool> *mask = ihit < int(nhireshits) ?  maskHiRes_ : maskLowRes_;
                (*mask)[hits[ihit].id] = true;
            }
            allTrajectories.push_back(std::move(traj));
        }
    }
    if (!allTrajectories.empty()) {
        DEBUG1_printf("Final list of trajectories from this cluster:\n");
        trajectoryCleaner_->clean(allTrajectories);
        saveCandidates(allTrajectories, tcCollection);
    }
}

FreeTrajectoryState
TrackCandidateBuilderFromCluster::startingState(float eta0, float phi0, float alpha) const 
{
    // alpha = 0.5 * 0.003 * B[T] / pT[GeV] -->
    float pt = std::abs(alpha ? 0.5 * 0.003f * bfieldAtOrigin_ / alpha : 100.0f); 
    ROOT::Math::CylindricalEta3D<double> p3d(pt, eta0, phi0);
    GlobalPoint  x(hitsHiRes_->x0(), hitsHiRes_->y0(), hitsHiRes_->z0());
    GlobalVector p(p3d.x(), p3d.y(), p3d.z());
    GlobalTrajectoryParameters gtp(x,p, alpha > 0 ? -1 : +1, &*bfield_);
    DEBUG2_printf("Cluster for pT %8.4f, eta %+5.3f, phi %+5.3f, q %+2d\n", pt, eta0, phi0, alpha > 0 ? -1 : +1);

    AlgebraicSymMatrix55 cov;
    for (unsigned int i = 0; i < 5; ++i) cov(i,i) = startingCovariance_[i];

    return FreeTrajectoryState(gtp, CurvilinearTrajectoryError(cov));
}

void
TrackCandidateBuilderFromCluster::dumpAsSeed(const std::vector<ClusteredHit> &cluster, TrajectorySeedCollection &seedCollection) const 
{
    edm::OwnVector<TrackingRecHit> seedHits; 
    seedHits.reserve(cluster.size());
    for (const ClusteredHit &hit : cluster) seedHits.push_back(*hit.hit);
    FreeTrajectoryState state(startingState(cluster.front().eta, cluster.front().phi, hitsHiRes_->alpha()));
    TrajectoryStateOnSurface tsos(state, tracker_->idToDet(cluster.front().hit->geographicalId())->surface());
    PTrajectoryStateOnDet const & PTraj = trajectoryStateTransform::persistentState(tsos, cluster.front().hit->geographicalId().rawId());
    TrajectorySeed seed(PTraj,std::move(seedHits),alongMomentum);
    seedCollection.push_back(seed);
}

void 
TrackCandidateBuilderFromCluster::refitAlphaBetaCorr(const std::vector<ClusteredHit> &cluster, const int *ihits, unsigned int nhits, float &alphacorr, float &betacorr, float &phi0, float &eta0) const {
    // Linear fits:
    //    eta - eta0 = p[0] + p[1]/rho
    //    phi - phi0 = p[0] + p[1]*rho
    AlgebraicSymMatrix22 rhoij, rhoinvij; AlgebraicVector2 etai, phii;
    for (unsigned int i = 0; i < nhits; ++i) {
        if (ihits[i] == -1) continue;
        const ClusteredHit &hit = cluster[ihits[i]];
        rhoij(0,0) += 1;
        rhoij(0,1) += hit.rho;
        rhoij(1,1) += hit.rho*hit.rho;
        phii(0) += (hit.phi-phi0);
        phii(1) += (hit.phi-phi0)*hit.rho;
        rhoinvij(0,0) += 1;
        rhoinvij(0,1) += 1/hit.rho;
        rhoinvij(1,1) += 1/(hit.rho*hit.rho);
        etai(0) += (hit.eta-eta0);
        etai(1) += (hit.eta-eta0)/hit.rho;
    }
    if (rhoij.Invert()) {
        AlgebraicVector2 dphi = rhoij*phii; 
        phi0 += dphi[0];
        alphacorr = dphi[1];
    }
    if (rhoinvij.Invert()) {
        AlgebraicVector2 deta = rhoinvij*etai; 
        eta0 += deta[0];
        betacorr = deta[1];
    }
    if (DEBUG>=2) {
        std::vector<const TrackingRecHit *> hitv;
        for (unsigned int i = 0; i < nhits; ++i) {
            if (ihits[i] == -1) continue;
            hitv.push_back(cluster[ihits[i]].hit);
        }
        HTDebugger::debugHelixParameters(hitv, hitsHiRes_->alpha() + alphacorr, betacorr, phi0, eta0, bfieldAtOrigin_, GlobalPoint(hitsHiRes_->x0(), hitsHiRes_->y0(), hitsHiRes_->z0()));
    }
}

const TrajectorySeed *
TrackCandidateBuilderFromCluster::makeSeed(const std::vector<ClusteredHit> &hits, std::pair<int,int> seedhits, const std::array<int,20> &layerhits,  float eta0, float phi0, float alpha, TrajectorySeedCollection &out) const 
{
    // sort hits along R, since the layer number sorting is not correct in the endcaps
    typedef std::pair<int,const TrackingRecHit *> HitP;
    std::vector<std::pair<float,HitP>> hitsSorted; std::list<SiStripRecHit2D> splittedHits; 
    for (unsigned int il = 0; il < layerhits.size(); ++il) {
        if (layerhits[il] == -1) continue;
        if (splitSeedHits_ && (typeid(*hits[layerhits[il]].hit) == typeid(SiStripMatchedRecHit2D))) {
            const SiStripMatchedRecHit2D &matched = dynamic_cast<const SiStripMatchedRecHit2D &>(*hits[layerhits[il]].hit);
            const GeomDet *  monoDet = tracker_->idToDet(DetId(matched.monoId()));
            const GeomDet *stereoDet = tracker_->idToDet(DetId(matched.stereoId()));
            float epsilon = 0; 
            if (dynamic_cast<const GeomDetUnit *>(monoDet)->type().isBarrel()) {
                epsilon = monoDet->position().perp() - stereoDet->position().perp();
            } else {
                epsilon = monoDet->position().z() - stereoDet->position().z();
                if (eta0 < 0) epsilon = -epsilon;
            }
            splittedHits.push_back(matched.monoHit());
            hitsSorted.push_back(std::make_pair(hypot(hits[layerhits[il]].rho, hits[layerhits[il]].z - hitsHiRes_->z0())+epsilon, HitP(layerhits[il],&splittedHits.back())));
            splittedHits.push_back(matched.stereoHit());
            hitsSorted.push_back(std::make_pair(hypot(hits[layerhits[il]].rho, hits[layerhits[il]].z - hitsHiRes_->z0())-epsilon, HitP(layerhits[il],&splittedHits.back())));
        } else {
            hitsSorted.push_back(std::make_pair(hypot(hits[layerhits[il]].rho, hits[layerhits[il]].z - hitsHiRes_->z0()), HitP(layerhits[il],hits[layerhits[il]].hit)));
        }
    }
    std::sort(hitsSorted.begin(), hitsSorted.end());

    // ok now we build the seed
    FreeTrajectoryState fts(startingState(eta0,phi0,alpha));
    TrajectoryStateOnSurface tsos;
    edm::OwnVector<TrackingRecHit> seedHits; seedHits.reserve(hitsSorted.size());
    KFUpdator updator;
    const TrackingRecHit* lasthit = 0;  
    unsigned int myhits = 0, myskip = 0;
    for (unsigned int ipair = 0, npair = hitsSorted.size(); ipair < npair; ++ipair) {
        const TrackingRecHit* hit = hitsSorted[ipair].second.second;
        bool firstOfMatch  = (ipair < npair-1 && abs(hitsSorted[ipair+1].second.second->geographicalId().rawId()-hit->geographicalId().rawId()) == 1);
        TrajectoryStateOnSurface state = (lasthit == 0) ?
            propagator_->propagate(fts,  tracker_->idToDet(hit->geographicalId())->surface())
            : propagator_->propagate(tsos, tracker_->idToDet(hit->geographicalId())->surface());
        if (!state.isValid()) { DEBUG3_printf ("\tfailed propagation\n");
            if (tsos.isValid()) { 
                DEBUG4_printf("\t\t\tstarting state: rho = %9.4f, z = %+9.4f +/- %7.4f, phi = %+5.3f, eta = %+5.3f\n",
                        tsos.globalPosition().perp(), tsos.globalPosition().z(), std::sqrt(tsos.cartesianError().matrix()(2,2)), float(tsos.globalPosition().phi()), tsos.globalPosition().eta());
                TransientTrackingRecHit::RecHitPointer tth = hitBuilder_->build(hit);
                DEBUG4_printf("\t\t\trechit prefit : rho = %9.4f, z = %+9.4f +/- %7.4f, phi = %+5.3f, eta = %+5.3f (detid %10d)\n",
                        tth->globalPosition().perp(), tth->globalPosition().z(), std::sqrt(tth->globalPositionError().czz()), float(tth->globalPosition().phi()), tth->globalPosition().eta(), tth->geographicalId().rawId());
            }
            // for matched pairs, we mark the failure only if both fail 
            if (!firstOfMatch && (lasthit == 0 || abs(lasthit->geographicalId().rawId()-hit->geographicalId().rawId()) > 1)) myskip++;  
            if (myskip >= 2) break;
            continue; 
        }
        //printf("\t\t propagated state: 1/p = %+9.5f +/- %9.5f,  pt = %7.2f +/- %7.2f\n", 
        //                state.globalParameters().signedInverseMomentum(), std::sqrt(state.curvilinearError().matrix()(0,0)),
        //                state.globalMomentum().perp(), TrajectoryStateAccessor(*state.freeState()).inversePtError()*state.globalMomentum().perp2());
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
        if (!chi2.first) {
            // for matched pairs, we mark the failure only if both fail 
            if (!firstOfMatch && (lasthit == 0 || abs(lasthit->geographicalId().rawId()-hit->geographicalId().rawId()) > 2)) myskip++; 
            DEBUG3_printf("\tSkipping hit on detid %10d/%7d due to bad chi2 = %7.1f (%d consec. skipped hits)\n", hit->geographicalId().rawId(), hitid(hit), chi2.second, myskip);
            if (myskip >= 2) break;
            continue;
        } else if (myskip == 1 && chi2.second >= 0.6*estimator_->chiSquaredCut()) {
            // for matched pairs, we mark the failure only if both fail 
            if (!firstOfMatch && (lasthit == 0 || abs(lasthit->geographicalId().rawId()-hit->geographicalId().rawId()) > 1)) myskip++; 
            DEBUG3_printf("\tSkipping hit on detid %10d/%7d due to bad chi2 = %7.1f (tighter threshold of %6.1f used after one lost hit; %d consec. skipped hits)\n", hit->geographicalId().rawId(), hitid(hit), chi2.second, 0.6*estimator_->chiSquaredCut(), myskip);
            if (myskip >= 2) break;
            continue;
        } else {
            myskip = 0;
        }
        TrajectoryStateOnSurface updated = updator.update(state, *tth);
        if (!updated.isValid()) { DEBUG3_printf("\tfailed update\n"); break; } 
        tsos = updated;
        //printf("\t\t updated    state: 1/p = %+9.5f +/- %9.5f,  pt = %7.2f +/- %7.2f\n", 
        //                tsos.globalParameters().signedInverseMomentum(), std::sqrt(tsos.curvilinearError().matrix()(0,0)),
        //                tsos.globalMomentum().perp(), TrajectoryStateAccessor(*tsos.freeState()).inversePtError()*tsos.globalMomentum().perp2());
        seedHits.push_back(tth->hit()->clone());
        lasthit = hit; 
        myhits++;
        if (myhits > 1) DEBUG3_printf("\tKF-fitted %2d hits (last chi2 = %7.1f; pt = %7.2f +/- %7.2f, q = %+1d; last detid %10d/%7d)\n", myhits, chi2.second, tsos.globalMomentum().perp(), TrajectoryStateAccessor(*tsos.freeState()).inversePtError()*tsos.globalMomentum().perp2(), tsos.charge(),lasthit->geographicalId().rawId(),hitid(lasthit));
    }
    if (myhits >= minHits_) {
        PTrajectoryStateOnDet const & PTraj = trajectoryStateTransform::persistentState(tsos, lasthit->geographicalId().rawId());
        TrajectorySeed seed(PTraj,std::move(seedHits),alongMomentum);
        out.push_back(seed);
        DEBUG3_printf("--> saved as trajectory seed!\n");
        return &out.back();
    }
    return 0;
}

const TrajectorySeed *
TrackCandidateBuilderFromCluster::makeSeed2Way(const std::vector<ClusteredHit> &hits, std::pair<int,int> seedhits, const std::array<int,20> &layerhits,  float eta0, float phi0, float alpha, TrajectorySeedCollection &out) const 
{
    // sort hits along R, since the layer number sorting is not correct in the endcaps
    typedef std::pair<int,const TrackingRecHit *> HitP;
    std::vector<std::pair<float,HitP>> hitsSorted;
    for (unsigned int il = 0; il < layerhits.size(); ++il) {
        if (layerhits[il] == -1) continue;
        hitsSorted.push_back(std::make_pair(hypot(hits[layerhits[il]].rho, hits[layerhits[il]].z - hitsHiRes_->z0()), HitP(layerhits[il],hits[layerhits[il]].hit)));
    }
    std::sort(hitsSorted.begin(), hitsSorted.end());

    FreeTrajectoryState fts(startingState(eta0,phi0,alpha));
    TrajectoryStateOnSurface tsos;
    KFUpdator updator;

    // get the index of the innermost and outermost seed hit
    int isin = -1, isout = -1;
    for (unsigned int i = 0, n = hitsSorted.size(); i < n; ++i) {
        if (hitsSorted[i].second.first == seedhits.first || hitsSorted[i].second.first == seedhits.second) {
            if (isin == -1) {
                isin = i;
            } else { 
                isout = i; 
                break; 
            }
        }
    }
    assert(isin != -1 && isout != -1);
    const TrackingRecHit* lasthit = 0; 
    unsigned int istart = isin;
    unsigned int myhits = 0, myskip = 0;
    if (isin != 0 && isout != -1) {
        DEBUG3_printf("\tSeed hits: %d, %d. Will attempt back-propagation.\n", isin, isout);
        for (int iback = isout; iback >= 0; --iback) {
            const TrackingRecHit *hit = hitsSorted[iback].second.second;
            TrajectoryStateOnSurface state = (iback == isout) 
                        ?         propagator_->propagate(fts,  tracker_->idToDet(hit->geographicalId())->surface())
                        : propagatorOpposite_->propagate(tsos, tracker_->idToDet(hit->geographicalId())->surface());
            if (!state.isValid()) { 
                DEBUG3_printf ("\tfailed back-propagation\n"); 
                break; 
            }
            TransientTrackingRecHit::RecHitPointer tth = hitBuilder_->build(hit);
            tth = tth->clone(state);
            std::pair<bool,double> chi2 = estimator_->estimate(state, *tth);
            if (!chi2.first) {
                DEBUG3_printf("\tSkipping backwards hit on detid %10d/%7d due to bad chi2 = %7.1f\n", hit->geographicalId().rawId(), hitid(hit), chi2.second);
                myskip++; 
                if (myskip == 2) { break; }
                continue;
            } else {
                myskip = 0;
            }
            TrajectoryStateOnSurface updated = updator.update(state, *tth);
            if (!updated.isValid()) { DEBUG3_printf("\tfailed back-update\n"); break; } 
            tsos = updated;
            lasthit = hit;
            istart = iback;
            DEBUG3_printf("\tKF-back-fitted hits (last chi2 = %7.1f; pt = %7.2f +/- %7.2f, q = %+1d; last detid %10d/%7d)\n", chi2.second, tsos.globalMomentum().perp(), TrajectoryStateAccessor(*tsos.freeState()).inversePtError()*tsos.globalMomentum().perp2(), tsos.charge(),lasthit->geographicalId().rawId(),hitid(lasthit));
        }
        if (istart > unsigned(isin)) {
            DEBUG3_printf("\tDone back-fitting but failed to reach even the innermost seed hit. will discard the attempt.\n");
            tsos = TrajectoryStateOnSurface();
        } else {
            DEBUG3_printf("\tDone back-fitting successfully from hit %d to hit %d (innermost seed hit %d); now will start forward-fitting from there.\n",istart,isout,isin);
            //tsos.rescaleError(10.0); 
            tsos = TrajectoryStateOnSurface();
        }
    }
    // ok now we build the seed
    edm::OwnVector<TrackingRecHit> seedHits; seedHits.reserve(hitsSorted.size());
    myhits = 0; myskip = 0; lasthit = 0;
    for (unsigned int i = istart, n = hitsSorted.size(); i < n; ++i) { 
        const TrackingRecHit* hit = hitsSorted[i].second.second;
        TrajectoryStateOnSurface state = (lasthit == 0) ?
              propagator_->propagate(fts,  tracker_->idToDet(hit->geographicalId())->surface())
            : propagator_->propagate(tsos, tracker_->idToDet(hit->geographicalId())->surface());
        if (!state.isValid()) { DEBUG3_printf ("\tfailed propagation on %10d/%7d\n",hit->geographicalId().rawId(),hitid(hit)); break; }
        TransientTrackingRecHit::RecHitPointer tth = hitBuilder_->build(hit);
        tth = tth->clone(state);
        std::pair<bool,double> chi2 = estimator_->estimate(state, *tth);
        if (!chi2.first) {
            DEBUG3_printf("\tSkipping hit on detid %10d/%7d due to bad chi2 = %7.1f\n", hit->geographicalId().rawId(), hitid(hit), chi2.second);
            myskip++; if (myskip == 2) { break; }
            continue;
        } else {
            myskip = 0;
        }
        TrajectoryStateOnSurface updated = updator.update(state, *tth);
        if (!updated.isValid()) { DEBUG3_printf("\tfailed update\n"); break; } 
        tsos = updated;
        seedHits.push_back(tth->hit()->clone());
        lasthit = hit;
        myhits++;
        DEBUG3_printf("\tKF-fitted %2d hits (last chi2 = %7.1f; pt = %7.2f +/- %7.2f, q = %+1d; last detid %10d/%7d)\n", myhits, chi2.second, tsos.globalMomentum().perp(), TrajectoryStateAccessor(*tsos.freeState()).inversePtError()*tsos.globalMomentum().perp2(), tsos.charge(),lasthit->geographicalId().rawId(),hitid(lasthit));
    }
    if (myhits >= minHits_) {
        PTrajectoryStateOnDet const & PTraj = trajectoryStateTransform::persistentState(tsos, lasthit->geographicalId().rawId());
        TrajectorySeed seed(PTraj,std::move(seedHits),alongMomentum);
        out.push_back(seed);
        DEBUG3_printf("--> saved as trajectory seed!\n");
        return &out.back();
    }
    return 0;
}


void
TrackCandidateBuilderFromCluster::makeTrajectories(const TrajectorySeed &seed, std::vector<Trajectory> &out) const 
{
    out.clear();
    auto const & startTraj = trajectoryBuilder_->buildTrajectories(seed, out, nullptr);
    DEBUG3_printf("\t created %lu trajectories in-out (before cleaning) \n", out.size());
    if (cleanTrajectoryAfterInOut_) {
        trajectoryCleaner_->clean(out);
        for (Trajectory &t : out) { if (t.measurements().empty()) t.invalidate(); }
        out.erase(std::remove_if(out.begin(),out.end(), std::not1(std::mem_fun_ref(&Trajectory::isValid))), out.end());
    }
    DEBUG3_printf("\t created %lu trajectories in-out (after cleaning) \n", out.size());
    if (DEBUG>=3) for (const Trajectory &t : out) { if (t.isValid()) dumpTraj(t); }

    if (out.empty()) return;

    // make a copy of the in-out trajectory before rebuild
    std::vector<Trajectory> safe(out);

    trajectoryBuilder_->rebuildTrajectories(startTraj, seed, out);
    DEBUG3_printf("\t created %lu trajectories out-in (before cleaning) \n", out.size());
    trajectoryCleaner_->clean(out);
    for (Trajectory &t : out) { if (t.measurements().empty()) t.invalidate(); }
    out.erase(std::remove_if(out.begin(),out.end(), std::not1(std::mem_fun_ref(&Trajectory::isValid))), out.end());
    DEBUG3_printf("\t created %lu trajectories out-in (after cleaning)\n", out.size());
    if (DEBUG>=3) for (const Trajectory &t : out) { if (t.isValid()) dumpTraj(t); }
    if (!safe.empty() && out.empty()) {
        DEBUG3_printf("\t Somehow we lost the trajectory in the backwards fit. trying to understand why\n");
        KFUpdator updator;
        for (Trajectory &traj : safe) {
            Trajectory refit(traj.sharedSeed(), traj.direction()); 
            refit.reverse();
            Trajectory::DataContainer const & data = traj.measurements();
            refit.reserve(data.size());

            // determine where to start from
            TrajectoryStateOnSurface state; int lastState = -1;
            for (lastState = data.size()-2; lastState >= 0; --lastState) {
                if (data[lastState].updatedState().isValid() && data[lastState].recHitR().isValid()) {
                    state = data[lastState].updatedState();
                    state.rescaleError(100.);
                    DEBUG3_printf("\t\tstarting state: %d on detid %10d/%7d\n", lastState, data[lastState].recHit()->geographicalId().rawId(), hitid(data[lastState].recHit()->hit()));
                    break;
                }
            }
            if (lastState == -1) { DEBUG3_printf("\t\tcouldn't find a state to restart the refit\n"); continue; }

            unsigned int nfail = 0;
            for (int istate = lastState; istate >= 0; --istate) {
                TransientTrackingRecHit::ConstRecHitPointer hit = data[istate].recHit();
                if (hit->det() == 0) {
                    DEBUG3_printf("\t\tstate with no det: no idea what to do\n"); 
                    continue;
                }
                DEBUG3_printf("\t\t\tstarting state: rho = %6.2f, z = %+7.2f phi = %+5.3f, eta = %+5.3f\n", state.globalPosition().perp(), state.globalPosition().z(), float(state.globalPosition().phi()), state.globalPosition().eta());
                TrajectoryStateOnSurface tsos = propagatorOpposite_->propagate(state, hit->det()->surface());
                if (!tsos.isValid()) {
                    DEBUG3_printf("\t\tfailed propagation to hit %10d/%7d, will skip it\n", hit->geographicalId().rawId(), hitid(hit->hit()));
                    nfail++;
                    if (nfail > 2) break; else continue;
                }
                DEBUG3_printf("\t\t\tpropag.  state: rho = %6.2f, z = %+7.2f +/- %5.2f, phi = %+5.3f, eta = %+5.3f, u = %+7.3f +/- %6.3f, v = %+7.3f +/- %6.3f\n", 
                                    tsos.globalPosition().perp(), tsos.globalPosition().z(), std::sqrt(tsos.cartesianError().matrix()(2,2)), float(tsos.globalPosition().phi()), tsos.globalPosition().eta(),
                                    tsos.localParameters().position().x(), std::sqrt(tsos.localError().matrix()(3,3)), tsos.localParameters().position().y(), std::sqrt(tsos.localError().matrix()(4,4)));
                if (!hit->isValid()) {
                    DEBUG3_printf("\t\tstate has an invalid hit: adding it and carrying on\n");
                    refit.push(TrajectoryMeasurement(tsos,hit));
                    state = tsos; nfail = 0;
                } else {
                    hit = hit->clone(tsos);
                    DEBUG3_printf("\t\t\trechit pos.   : rho = %6.2f, z = %+7.2f +/- %5.2f, phi = %+5.3f, eta = %+5.3f, u = %+7.3f +/- %6.3f, v = %+7.3f +/- %6.3f\n", 
                                        hit->globalPosition().perp(), hit->globalPosition().z(), std::sqrt(hit->globalPositionError().czz()), float(hit->globalPosition().phi()), hit->globalPosition().eta(),
                                        hit->localPosition().x(), std::sqrt(hit->parametersError()(1,1)), hit->localPosition().y(), std::sqrt(hit->parametersError()(2,2)));
                    std::pair<bool,double> chi2 = estimator_->estimate(tsos, *hit);
                    if (chi2.first) {
                        TrajectoryStateOnSurface updated = updator.update(tsos, *hit); 
                        if (!updated.isValid()) {
                            DEBUG3_printf("\t\tupdate failed on hit on detid %10d/%7d with chi2 = %7.1f, replace with invalid hit\n", hit->geographicalId().rawId(), hitid(hit->hit()), chi2.second);
                            hit = InvalidTransientRecHit::build(hit->det());
                            refit.push(TrajectoryMeasurement(tsos,hit));
                            state = tsos;
                            nfail++;
                            if (nfail > 2) break;
                        } else {
                            DEBUG3_printf("\t\tgood hit on detid %10d/%7d with chi2 = %7.1f\n", hit->geographicalId().rawId(), hitid(hit->hit()), chi2.second);
                            refit.push(TrajectoryMeasurement(tsos,updated,hit,chi2.second));
                            nfail = 0;
                            state = updated;
                        }
                    } else {
                        DEBUG3_printf("\t\tstate has hit detid %10d/%7d with bad chi2 = %7.1f; replace with invalid hit\n", hit->geographicalId().rawId(), hitid(hit->hit()), chi2.second);
                        hit = InvalidTransientRecHit::build(hit->det());
                        refit.push(TrajectoryMeasurement(tsos,hit));
                        state = tsos;
                        nfail++;
                        if (nfail > 2) break;
                    }
                }
            }
            while (!refit.empty() && !refit.lastMeasurement().recHit()->isValid()) {
                refit.pop();
            }
            refit.reverse();
            DEBUG3_printf("\tManually refitted trajectory:\n");
            if (DEBUG>=3) dumpTraj(refit);
            out.push_back(std::move(refit));
        }    
    } else {
        bool dropped = false;
        for (Trajectory &t : out) { 
            if (!t.isValid()) continue;
            if (t.foundHits() < int( t.chiSquared() == 0 ? minHitsIfNoCkf_ : minHits_)) {
                t.invalidate();
                dropped = true;
           }
        }
        if (dropped) {
            out.erase(std::remove_if(out.begin(),out.end(), std::not1(std::mem_fun_ref(&Trajectory::isValid))), out.end());
            DEBUG3_printf("\t created %lu trajectories out-in (after my filtering)\n", out.size());
        }
    }
}


void
TrackCandidateBuilderFromCluster::saveCandidates(const std::vector<Trajectory> &cands,  TrackCandidateCollection & tcCollection) const 
{
    for (const Trajectory & traj : cands) {
        if (!traj.isValid() || traj.measurements().empty()) continue;

        Trajectory::RecHitContainer thits;
        traj.recHitsV(thits,useHitsSplitting_);
        edm::OwnVector<TrackingRecHit> recHits;
        recHits.reserve(thits.size());
        for (const auto & hit : thits) {
            recHits.push_back( hit->hit()->clone() );
        }

        TransientTrackingRecHit::ConstRecHitPointer firsthit = traj.firstMeasurement().recHit();
        if (!firsthit->isValid()) {
            printf("\tFirst hit of trajectory is not valid - makes no sense!!!\n");
            continue;
        }
        TrajectoryStateOnSurface tsos = traj.firstMeasurement().updatedState();
        if (!tsos.isValid()) tsos = traj.firstMeasurement().backwardPredictedState();
        if (!tsos.isValid()) tsos = traj.firstMeasurement().forwardPredictedState();
        TrajectoryStateOnSurface tsos2 = traj.lastMeasurement().updatedState();
        if (!tsos2.isValid()) tsos2 = traj.lastMeasurement().backwardPredictedState();
        if (!tsos2.isValid()) tsos2 = traj.lastMeasurement().forwardPredictedState();
        if (tsos.isValid()) { 
            float phipos = tsos.globalMomentum().phi(); if (phipos < 0) phipos += 2*M_PI;
            DEBUG1_printf("\tTrajectory pt %8.4f +/- %7.2f, eta %+5.3f, phi %+5.3f, q %+2d, hits %2d (rho = %5.1f, z = %+5.1f --> rho = %5.1f, z = %+5.1f)\n",  
                tsos.globalMomentum().perp(), TrajectoryStateAccessor(*tsos.freeState()).inversePtError()*tsos.globalMomentum().perp2(), 
                tsos.globalMomentum().eta(), phipos, 
                tsos.charge(), traj.foundHits(),
                tsos.globalPosition().perp(), tsos.globalPosition().z(),
                tsos2.isValid() ? tsos2.globalPosition().perp() : 0, tsos2.isValid() ? tsos2.globalPosition().z() : 0);
            if (DEBUG>=1) HTDebugger::printAssociatedTracks(traj,(DEBUG==1?3:1));
        } else {
            printf("\tinitial state is not valid?\n");
            continue;
        }

        PTrajectoryStateOnDet state;
        if(useHitsSplitting_ && firsthit->geographicalId() != thits.front()->geographicalId() ){
            TrajectoryStateOnSurface propagated = anyPropagator_->propagate(tsos,thits.front()->det()->surface());
            if (!propagated.isValid()) { printf("\tcannot fix initial state by propagation??\n"); continue; }
            tsos = propagated;
            firsthit = thits.front();
        }
        state = trajectoryStateTransform::persistentState(tsos, firsthit->geographicalId().rawId());
        tcCollection.push_back(TrackCandidate(recHits,traj.seed(),state,traj.seedRef(),traj.nLoops() ) );

    }
}

void
TrackCandidateBuilderFromCluster::indexHits(const std::vector<ClusteredHit> &hits, std::vector<std::pair<uint32_t,unsigned int>> &index) const 
{
    unsigned int n = hits.size();
    index.resize(n);
    for (unsigned int i = 0; i < n; ++i) {
        index.push_back(std::make_pair(hits[i].hit->geographicalId().rawId(), i));
    }
    std::sort(index.begin(), index.end());
}

int 
TrackCandidateBuilderFromCluster::findHit(const TrackingRecHit *hit, const std::vector<ClusteredHit> &hits, HitIndex &index) const 
{
    //printf("\t\t find hit %u in %lu hits\n", detid, hits.size());
    if (index.empty()) indexHits(hits, index);
    uint32_t detid = hit->geographicalId().rawId();
    HitIndex::const_iterator match = std::lower_bound(index.begin(), index.end(), std::make_pair(detid,0u));
    for (;match != index.end() && match->first == detid; ++match) {
        if (hits[match->second].hit->sharesInput(hit, TrackingRecHit::all)) return match->second;
    } 
    return -1;
}

void
TrackCandidateBuilderFromCluster::dumpTraj(const Trajectory &traj) const
{
    printf("\t\tpropagation direction: %s\n", (traj.direction() == alongMomentum ? "along" : "opposite"));
    for (const TrajectoryMeasurement &tm : traj.measurements()) {
        TrajectoryStateOnSurface tsos;
        if (tm.updatedState().isValid()) {
            tsos = tm.updatedState();
            printf("\t\tstate [updated]: ");
        } else if (tm.forwardPredictedState().isValid()) {
            tsos = tm.forwardPredictedState();
            printf("\t\tstate [fw pred]: ");
        } else if (tm.backwardPredictedState().isValid()) {
            tsos = tm.backwardPredictedState();
            printf("\t\tstate [bw pred]: ");
        } else {
            printf("\t\t no state?");
        }
        if (tsos.isValid()) { 
            float phipos = tsos.globalMomentum().phi(); if (phipos < 0) phipos += 2*M_PI;
            printf("rho = %5.1f, z = %+6.1f, eta = %+5.3f, phi = %+5.3f, pt = %7.2f +/- %7.2f, q = %+1d   ", 
                tsos.globalPosition().perp(), tsos.globalPosition().z(),
                tsos.globalMomentum().eta(), phipos, 
                tsos.globalMomentum().perp(), TrajectoryStateAccessor(*tsos.freeState()).inversePtError()*tsos.globalMomentum().perp2(), tsos.charge());
 
        }
        if (tm.recHit()->isValid()) {
            float phipos = tm.recHit()->globalPosition().phi(); if (phipos < 0) phipos += 2*M_PI;
            printf("hit (eta %+5.3f, phi %+5.3f, detid %10d/%6d, chi2 %5.1f)", tm.recHit()->globalPosition().eta(), phipos, tm.recHit()->geographicalId().rawId(),hitid(tm.recHit()->hit()),tm.estimate());
        } else {
            printf("no hit");
        }
        printf("\n");
    }
}

int TrackCandidateBuilderFromCluster::hitid(const TrackingRecHit *hit)
{
    if (!hit->isValid()) return 0;
    if (typeid(*hit) == typeid(SiPixelRecHit)) {
        return static_cast<const SiPixelRecHit *>(hit)->cluster().key();
    } else if (typeid(*hit) == typeid(SiStripMatchedRecHit2D)) {
        return static_cast<const SiStripMatchedRecHit2D *>(hit)->monoClusterRef().key();
    } else {
        const TrackerSingleRecHit *sihit = dynamic_cast<const TrackerSingleRecHit *>(hit);
        if (sihit) return sihit->omniClusterRef().key(); 
    }
    return 666666;
}
