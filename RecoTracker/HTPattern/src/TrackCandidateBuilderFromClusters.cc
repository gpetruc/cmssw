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
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit1D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Math/GenVector/CylindricalEta3D.h"
#include "RecoTracker/HTPattern/interface/HTDebugger.h"
#include "RecoTracker/HTPattern/interface/SimplisticTrajectoryBuilder.h"
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
    useSimpleTB_(iConfig.getParameter<bool>("useSimpleTB")),
    simpleTB_(useSimpleTB_ ? new SimplisticTrajectoryBuilder(iConfig.getParameter<edm::ParameterSet>("simpleTBParameters")) : 0),
    useHitsSplitting_(iConfig.getParameter<bool>("useHitsSplitting")),
    splitSeedHits_(iConfig.getParameter<bool>("splitSeedHits")),
    dphiCutPair_(iConfig.getParameter<double>("dphiCutPair")),
    detaCutPair_(iConfig.getParameter<double>("detaCutPair")),
    pairsOnSeedCellOnly_(iConfig.getParameter<bool>("pairsOnSeedCellOnly")),
    minHits_(iConfig.getParameter<uint32_t>("minHits")),
    maxSeedHits_(iConfig.getParameter<uint32_t>("maxSeedHits")),
    minHitsIfNoCkf_(iConfig.getParameter<uint32_t>("minHitsIfNoCkf")),
    popOneLayer_(iConfig.getParameter<bool>("popOneLayer")),
    keepPixelTriplets_(iConfig.getParameter<bool>("keepPixelTriplets")),
    maxFailMicroClusters_(iConfig.getParameter<uint32_t>("maxFailMicroClusters")),
    startingCovariance_(iConfig.getParameter<std::vector<double>>("startingCovariance")),
    initialStateEstimator_(0),
    alphaMin_(-9e9),alphaMax_(9e9)
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

    es.get<IdealGeometryRecord>().get(tTopo_);

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

    if (useSimpleTB_) simpleTB_->setEvent(evt,es);

    allTrajectories_.reserve(200);
    timerAll_.Reset();
    timerSeed_.Reset();
    timerTraj_.Reset();
}

void 
TrackCandidateBuilderFromCluster::done(TrackCandidateCollection &tcCollection) 
{
    if (!allTrajectories_.empty()) {
        DEBUG1_printf("Final list of trajectories from HT:\n");
        trajectoryCleaner_->clean(allTrajectories_);
        saveCandidates(allTrajectories_, tcCollection);
    }
    allTrajectories_.clear();
    trajectoryBuilder_.reset();
    navigationSetter_.reset();
    timerAll_.Stop();
    timerSeed_.Stop();
    timerTraj_.Stop();
    printf("Time: %6.2f s total, of which %6.2f s for seeding, %6.2f for trajectory building\n", timerAll_.RealTime(), timerSeed_.RealTime(), timerTraj_.RealTime());
}

bool 
TrackCandidateBuilderFromCluster::prefilterCluster(const HTCluster &cluster) const 
{
    const HTCell &seed = mapHiRes_->get(cluster.ieta(), cluster.iphi());
    for (unsigned int ii1 = 0, n = seed.hits().size(); ii1 < n-1; ++ii1) {
        unsigned int i1 = seed.hits()[ii1];
        if ((*maskHiRes_)[i1]) continue;
        for (unsigned int ii2 = ii1+1; ii2 < n; ++ii2) {
            unsigned int i2 = seed.hits()[ii2];
            if ((*maskHiRes_)[i2]) continue;
            if (hitsHiRes_->layermask(i1) == hitsHiRes_->layermask(i2)) continue;
            if (std::abs(hitsHiRes_->eta(i1) - hitsHiRes_->eta(i2)) > detaCutPair_) continue;
            float dphi = std::abs(hitsHiRes_->phi(i1) - hitsHiRes_->phi(i2)); if (dphi > M_PI) dphi = 2*M_PI-dphi;
            if (dphi > dphiCutPair_) continue;
            if (abs(mask2layer(hitsHiRes_->layermask(i1)) - mask2layer(hitsHiRes_->layermask(i2))) > 2) continue;
            return true;
        }
    }
    return false;
}

void
TrackCandidateBuilderFromCluster::run(const HTCluster &cluster, unsigned int nseedlayersCut, unsigned int nlayersCut, TrajectorySeedCollection & seedCollection, TrajectorySeedCollection *seedsFromAllClusters) 
{
    timerAll_.Start(false);
    //if (pairsOnSeedCellOnly_ && !prefilterCluster(cluster)) return;
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
    if (countbits(layermask) < nlayersCut) {
        timerAll_.Stop();
        return;
    }
    unsigned int nhits = hits.size();
    DEBUG2_printf("found %d seed hits, %d hi-res hits, %d hits\n", nseedhits, nhireshits, nhits);
    if (seedsFromAllClusters) if (DEBUG>0)  dumpAsSeed(hits, *seedsFromAllClusters);
    if (seedsFromAllClusters) if (DEBUG>=1) HTDebugger::printAssociatedTracks(seedsFromAllClusters->back());
    if (seedsFromAllClusters) if (DEBUG>=2) HTDebugger::beginLoggingCluster(seedsFromAllClusters->back(), nseedhits, nhireshits, nhits, hitsHiRes_->alpha());

    typedef std::pair<uint16_t,uint16_t> hitpair;
    typedef std::pair<float, hitpair> distpair;
    std::vector<distpair> pairs;
    //float dphidrMax = 0.5*0.003*std::abs(bfieldAtOrigin_)/(alphaMax_-alphaMin_); 
    for (unsigned int i1 = 0; i1 < nseedhits; ++i1) {
        for (unsigned int i2 = i1+1; i2 < nhits; ++i2) {
            if (pairsOnSeedCellOnly_ && i2 >= nseedhits) break;
            if (hits[i2].layer == hits[i1].layer) continue; // no same-layer pairs in seeding
            if (std::abs(hits[i1].layer - hits[i2].layer) <= 2) { // to too-far-away-pairs
                float dphi = hits[i1].dphi(hits[i2]);
                float deta = hits[i1].deta(hits[i2]);
                DEBUG3_printf("\t pair %3d, %3d: (%+5.3f,%+5.3f, ly %2d)  (%+5.3f,%+5.3f ly %2d): dphi %.4f, deta %.4f, %s\n", i1,i2, 
                            hits[i1].eta,hits[i1].phi,hits[i1].layer, hits[i2].eta,hits[i2].phi,hits[i2].layer, dphi, deta, HTDebugger::twoHitsOnSameTrack(hits[i1].hit, hits[i2].hit) ? "(good)" : "(not good)");
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
        //DEBUG2_printf("\t alpha uncropped = %+8.5f (pt %8.4f)\n", alphacorr+hitsHiRes_->alpha(), 0.5f*0.003f*bfieldAtOrigin_/std::max(std::abs(alphacorr+hitsHiRes_->alpha()),0.0000001f));
        if (alphacorr + hitsHiRes_->alpha() > alphaMax_) alphacorr = alphaMax_ - hitsHiRes_->alpha();
        if (alphacorr + hitsHiRes_->alpha() < alphaMin_) alphacorr = alphaMin_ - hitsHiRes_->alpha();
        //DEBUG2_printf("\t alpha   cropped = %+8.5f (pt %8.4f)\n", alphacorr+hitsHiRes_->alpha(), 0.5f*0.003f*bfieldAtOrigin_/std::max(std::abs(alphacorr+hitsHiRes_->alpha()),0.0000001f));
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
            DEBUG2_printf("   iteration %d: eta0,phi0 = %+5.3f,%+5.3f     alphacorr = %+8.5f (pt %8.4f), betacorr = %+8.5f\n", iteration+1, eta0, phi0, alphacorr, 0.5f*0.003f*bfieldAtOrigin_/std::max(std::abs(alphacorr+hitsHiRes_->alpha()),0.0000001f), betacorr);
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
                if (DEBUG >= 5 || (dphi <= dphiCutHits_[iteration] && deta <= detaCutHits_[iteration])) {
                    DEBUG4_printf("\thit %2d %c (%+5.3f,%+5.3f, ly %2d): d2d %.4f, d2dc %.4f    deta %.4f detac %.4f detacE %.4f   dphi %.4f dphic %.4f    detid %10d/%7d %s\n",i3, 
                                  i3 < nseedhits ? 'S' : (i3 < nhireshits ? 'H' : 'L'),
                                  hits[i3].eta,hits[i3].phi,hits[i3].layer, hits[i3].dist2d(eta0,phi0,0.f,0.f), dc, hits[i3].deta(eta0,0), hits[i3].deta(eta0,betacorr), deta, hits[i3].dphi(phi0,0), hits[i3].dphi(phi0,alphacorr),
                                  hits[i3].hit->geographicalId().rawId(),hitid(hits[i3].hit),
                                  HTDebugger::thirdHitOnSameTrack(hits[i1].hit, hits[i2].hit, hits[i3].hit) ? " (good)" : " (not good)");
                    if (DEBUG>=4) HTDebugger::logThirdHit(hits[i3].hit, hits[i3].layer, i3 < nseedhits ? 0 : (i3 < nhireshits ? 1 : 2),
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
                //DEBUG2_printf("\t alpha uncropped = %+8.5f (pt %8.4f)\n", alphacorr+hitsHiRes_->alpha(), 0.5f*0.003f*bfieldAtOrigin_/std::max(std::abs(alphacorr+hitsHiRes_->alpha()),0.0000001f));
                if (alphacorr + hitsHiRes_->alpha() > alphaMax_) alphacorr = alphaMax_ - hitsHiRes_->alpha();
                if (alphacorr + hitsHiRes_->alpha() < alphaMin_) alphacorr = alphaMin_ - hitsHiRes_->alpha();
                //DEBUG2_printf("\t alpha   cropped = %+8.5f (pt %8.4f)\n", alphacorr+hitsHiRes_->alpha(), 0.5f*0.003f*bfieldAtOrigin_/std::max(std::abs(alphacorr+hitsHiRes_->alpha()),0.0000001f));
            }
        }
        unsigned int microclusterpixelhits = (layerhits[0] != -1) + (layerhits[1] != -1) + (layerhits[2] != -1);
        DEBUG2_printf("microcluster %d (%d hits, %d pixel hits)\n", microclusters, microclusternhits, microclusterpixelhits);
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
        if (microclusternhits < minHits_ && (!keepPixelTriplets_ || microclusterpixelhits < 3)) continue;

        std::vector<Trajectory> trajectories;
        if (!useSimpleTB_) {
            const TrajectorySeed *seed = makeSeed(hits, std::make_pair(i1,i2), layerhits, eta0,phi0,hitsHiRes_->alpha()+alphacorr, seedCollection);
            if (seed == 0) {
                seed = 0; // makeSeed2Way(hits, std::make_pair(i1,i2), layerhits, eta0,phi0,hitsHiRes_->alpha()+alphacorr, seedCollection); // FIXME
            }
            if (seed == 0) continue;
            if (DEBUG>=2) HTDebugger::printAssociatedTracks(*seed);

            makeTrajectories(*seed, trajectories);
        } else {
            makeSimpleTrajectories(hits, std::make_pair(i1,i2), layerhits, eta0,phi0,hitsHiRes_->alpha()+alphacorr, trajectories);
        }
        
        if (trajectories.empty()) continue;
        for (Trajectory &traj : trajectories) {
            if (!traj.isValid() || traj.measurements().empty()) continue;
            int quality = qualityCut(traj); 
            if (quality == -1) continue;
            DEBUG2_printf("built trajectory with %d found hits, %d lost hits, quality %d\n", traj.foundHits(), traj.lostHits(), quality);
            if (DEBUG >= 2) HTDebugger::printAssociatedTracks(traj);
            if (quality > 0) failMicroClusters = 0;
            for (const TrajectoryMeasurement &tm : traj.measurements()) {
                if (!tm.recHit()->isValid()) continue;
                int ihit = findHit(tm.recHit()->hit(), hits, hitindex);
                if (ihit == -1) continue;
                // for quality >= 1, mask hits in the micro-cluster
                if (quality > 0) microcluster[ihit] = microclusters;
                // for quality >= 2, mask hits in the full tracking 
                std::vector<bool> *mask = ihit < int(nhireshits) ?  maskHiRes_ : maskLowRes_;
                if (quality > 1) (*mask)[hits[ihit].id] = true;
            }
            allTrajectories.push_back(std::move(traj));
        }
        unsigned int usedhits = 0;
        for (unsigned int i = 0; i < nhits; ++i) if (microcluster[i]) usedhits++;
        DEBUG3_printf("number of masked hits on the cluster: %d\n",usedhits);
    }
    if (!allTrajectories.empty()) {
        DEBUG1_printf("Final list of trajectories from this cluster:\n");
        trajectoryCleaner_->clean(allTrajectories);
        for (Trajectory &t : allTrajectories) {
            if (t.isValid()) {
                DEBUG1_printf("\ttrajectory with %d found hits, %d lost hits\n", t.foundHits(), t.lostHits());
                if (DEBUG >= 1) HTDebugger::printAssociatedTracks(t);
                allTrajectories_.push_back(std::move(t));
            }
        }
    }
    timerAll_.Stop();
}

FreeTrajectoryState
TrackCandidateBuilderFromCluster::startingState(float pseudoeta0, float phi0, float alpha) const 
{
    float eta0 = HTHitsSpher::etaFromPseudoEta(pseudoeta0);
    
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

TrajectoryStateOnSurface
TrackCandidateBuilderFromCluster::startingState(float pseudoeta0, float phi0, float alpha, const TrackingRecHit *hit) const 
{
    float eta0 = HTHitsSpher::etaFromPseudoEta(pseudoeta0);
    const GeomDet *det = tracker_->idToDet(DetId(hit->geographicalId()));
    //TransientTrackingRecHit::RecHitPointer tth = hitBuilder_->build(hit);
    LocalPoint  lx(0,0,0); // tth->localPosition());
    GlobalPoint  x = det->toGlobal(lx); //tth->globalPosition();
    float pt = std::abs(alpha ? 0.5 * 0.003f * bfieldAtOrigin_ / alpha : 100.0f);
    math::RhoEtaPhiVectorF p(pt, eta0, phi0 + 2*alpha*x.perp());
    LocalVector lp = det->toLocal(GlobalVector(p.X(), p.Y(), p.Z()));
    LocalTrajectoryParameters ltp(lx, lp, alpha > 0 ? -1 : +1);
    AlgebraicSymMatrix55 cov;

    float ptinv2 = std::max(0.25f, 1/(pt*pt));
    float pscale = (1/lp.mag2())/ptinv2; // convert from 1/pt^2 to 1/p^2
    cov(0,0) = startingCovariance_[0] * pscale;
    for (unsigned int i = 1; i < 5; ++i) cov(i,i) = startingCovariance_[i];

    return TrajectoryStateOnSurface(ltp, LocalTrajectoryError(cov), det->surface(), &*bfield_);
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
    if (DEBUG>=4) {
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
    timerSeed_.Start(false);
    // sort hits along R, since the layer number sorting is not correct in the endcaps
    typedef std::pair<int,const TrackingRecHit *> HitP;
    std::vector<std::pair<float,HitP>> hitsSorted; std::list<SiStripRecHit2D> splittedHits2D; std::list<SiStripRecHit1D> splittedHits1D; 
    int minlayer = 20, maxlayer = 0;
    for (unsigned int il = 0; il < layerhits.size(); ++il) {
        if (layerhits[il] == -1) continue;
        if (hits[layerhits[il]].layer > maxlayer) maxlayer = hits[layerhits[il]].layer;
        if (hits[layerhits[il]].layer < minlayer) minlayer = hits[layerhits[il]].layer;
        if (splitSeedHits_ && (typeid(*hits[layerhits[il]].hit) == typeid(SiStripMatchedRecHit2D))) {
            const SiStripMatchedRecHit2D &matched = dynamic_cast<const SiStripMatchedRecHit2D &>(*hits[layerhits[il]].hit);
            const GeomDet *  monoDet = tracker_->idToDet(DetId(matched.monoId()));
            const GeomDet *stereoDet = tracker_->idToDet(DetId(matched.stereoId()));
            float epsilon = 0; 
            if (dynamic_cast<const GeomDetUnit *>(monoDet)->type().isBarrel()) {
                epsilon = monoDet->position().perp() - stereoDet->position().perp();
                splittedHits1D.push_back(SiStripRecHit1D(LocalPoint(),LocalError(),DetId(matched.monoId()),matched.monoClusterRef()));
                hitsSorted.push_back(std::make_pair(hypot(hits[layerhits[il]].rho, hits[layerhits[il]].z - hitsHiRes_->z0())+epsilon, HitP(layerhits[il],&splittedHits1D.back())));
                splittedHits1D.push_back(SiStripRecHit1D(LocalPoint(),LocalError(),DetId(matched.stereoId()),matched.stereoClusterRef()));
                hitsSorted.push_back(std::make_pair(hypot(hits[layerhits[il]].rho, hits[layerhits[il]].z - hitsHiRes_->z0())-epsilon, HitP(layerhits[il],&splittedHits1D.back())));
            } else {
                epsilon = monoDet->position().z() - stereoDet->position().z();
                if (eta0 < 0) epsilon = -epsilon;
                splittedHits2D.push_back(matched.monoHit());
                hitsSorted.push_back(std::make_pair(hypot(hits[layerhits[il]].rho, hits[layerhits[il]].z - hitsHiRes_->z0())+epsilon, HitP(layerhits[il],&splittedHits2D.back())));
                splittedHits2D.push_back(matched.stereoHit());
                hitsSorted.push_back(std::make_pair(hypot(hits[layerhits[il]].rho, hits[layerhits[il]].z - hitsHiRes_->z0())-epsilon, HitP(layerhits[il],&splittedHits2D.back())));
            }
        } else {
            hitsSorted.push_back(std::make_pair(hypot(hits[layerhits[il]].rho, hits[layerhits[il]].z - hitsHiRes_->z0()), HitP(layerhits[il],hits[layerhits[il]].hit)));
        }
    }
    std::sort(hitsSorted.begin(), hitsSorted.end());

    float seedlayer = 0.5*(hits[seedhits.first].layer + hits[seedhits.second].layer);
    DEBUG2_printf("\tLayers: %d - %d  (avg seed layer: %4.1f)\n", minlayer,maxlayer,seedlayer);

    PropagationDirection pdir = alongMomentum;
    const Propagator *propagator = &*propagator_;
    if (seedlayer >= 7 && ((maxlayer-seedlayer) < (seedlayer-minlayer))) {
        pdir = oppositeToMomentum;
        propagator = &*propagatorOpposite_;
        std::reverse(hitsSorted.begin(), hitsSorted.end());
    }

    TrajectoryStateOnSurface tsos = startingState(eta0, phi0, alpha, hitsSorted.front().second.second);
    edm::OwnVector<TrackingRecHit> seedHits; seedHits.reserve(hitsSorted.size());
    KFUpdator updator;
    const TrackingRecHit* lasthit = 0;  
    unsigned int myhits = 0, myskip = 0, mypixelhits = 0, mylost = 0;
    int lastlayerIndex = 0, lastlayerDet = 0, lastlayerLayer = 0; TrajectoryStateOnSurface lastlayerState;
    for (unsigned int ipair = 0, npair = hitsSorted.size(); ipair < npair; ++ipair) {
        const TrackingRecHit* hit = hitsSorted[ipair].second.second;
        bool firstOfMatch  = (ipair < npair-1 && abs(hitsSorted[ipair+1].second.second->geographicalId().rawId()-hit->geographicalId().rawId()) == 1);
        TrajectoryStateOnSurface state = propagator->propagate(tsos, tracker_->idToDet(hit->geographicalId())->surface());
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
        double chi2cut = estimator_->chiSquaredCut();
        if (myskip == 1)       chi2cut *= 0.6;
        if (mylost > 0 )       chi2cut *= 0.6;
        if (myhits > minHits_) chi2cut *= 0.6;
        if (chi2.second > chi2cut) {
            // for matched pairs, we mark the failure only if both fail 
            if (!firstOfMatch && (lasthit == 0 || abs(lasthit->geographicalId().rawId()-hit->geographicalId().rawId()) > 2)) { 
                myskip++; 
                mylost++; 
            }
            DEBUG3_printf("\tSkipping hit on detid %10d/%7d due to bad chi2 = %7.1f (threshold %6.1f, %d consec. skipped hits)\n", hit->geographicalId().rawId(), hitid(hit), chi2.second, chi2cut, myskip);
            if (myskip >= 2) break;
            continue;
       } else {
            myskip = 0;
        }
        TrajectoryStateOnSurface updated = updator.update(state, *tth);
        if (!updated.isValid()) { DEBUG3_printf("\tfailed update\n"); break; } 
        int layer = tTopo_->layer(hit->geographicalId());
        if (lastlayerDet != hit->geographicalId().subdetId() || lastlayerLayer != layer) {
            lastlayerIndex = seedHits.size();
            lastlayerDet   = hit->geographicalId().subdetId();
            lastlayerLayer = layer;
            lastlayerState = tsos;
        }
        tsos = updated;
        seedHits.push_back(tth->hit()->clone());
        lasthit = hit; 
        myhits++; if (hit->geographicalId().subdetId() <= 2) mypixelhits++;
        DEBUG3_printf("\tKF-fitted %2d hits (last chi2 = %7.1f; pt = %7.2f +/- %7.2f, q = %+1d; last detid %10d/%7d)\n", myhits, chi2.second, tsos.globalMomentum().perp(), TrajectoryStateAccessor(*tsos.freeState()).inversePtError()*tsos.globalMomentum().perp2(), tsos.charge(),lasthit->geographicalId().rawId(),hitid(lasthit));
        DEBUG4_printf("\t\t local q/p = %+8.4f +/- %8.4f, dxdz = %+8.4f +/- %8.4f, dydz = %+8.4f +/- %8.4f, x = %+8.5f +/- %8.5f, y = %+8.5f +/- %8.5f\n",
                            tsos.localParameters().qbp(),  std::sqrt(tsos.localError().matrix()(0,0)), tsos.localParameters().dxdz(), std::sqrt(tsos.localError().matrix()(1,1)), tsos.localParameters().dydz(), std::sqrt(tsos.localError().matrix()(2,2)),
                            tsos.localParameters().position().x(), std::sqrt(tsos.localError().matrix()(3,3)), tsos.localParameters().position().y(), std::sqrt(tsos.localError().matrix()(4,4)));
        if (alphaMin_ > -99) {
            float alpha    = -0.5 * 0.003f * bfieldAtOrigin_ * tsos.globalParameters().signedInverseTransverseMomentum(); 
            float errAlpha =  std::abs(0.5 * 0.003f * bfieldAtOrigin_ * TrajectoryStateAccessor(*tsos.freeState()).inversePtError());
            if (alpha + 4*errAlpha < alphaMin_ || alpha -4*errAlpha > alphaMax_) {
                DEBUG3_printf("\talpha = %+7.4f +/- %6.4f outside the range [ %+7.4f, %+7.4f ] --> rejecting seed\n", alpha, errAlpha, alphaMin_, alphaMax_);
                timerSeed_.Stop();
                return 0;
            }
        }
        
        if (!firstOfMatch && myhits >= maxSeedHits_) {
            DEBUG3_printf("\tThat's enough for now, let's call it a seed...\n");
            break;
        }
    }
    if (myhits >= minHits_ || (keepPixelTriplets_ && mypixelhits == 3)) { 
        if (popOneLayer_ && lastlayerDet != 0) {
            while (int(seedHits.size()) > lastlayerIndex) seedHits.pop_back();
            tsos = lastlayerState;
            lasthit = & seedHits.back();
            DEBUG3_printf("\tRemoving one layer. Last hit will be %10d/%7d\n", lasthit->geographicalId().rawId(),hitid(lasthit));
        }
        PTrajectoryStateOnDet const & PTraj = trajectoryStateTransform::persistentState(tsos, lasthit->geographicalId().rawId());
        TrajectorySeed seed(PTraj,std::move(seedHits),pdir);
        out.push_back(seed);
        DEBUG3_printf("--> saved as trajectory seed!\n");
        timerSeed_.Stop();
        return &out.back();
    }
    timerSeed_.Stop();
    return 0;
}

void
TrackCandidateBuilderFromCluster::makeSimpleTrajectories(const std::vector<ClusteredHit> &hits, std::pair<int,int> seedhits, const std::array<int,20> &layerhits,  float eta0, float phi0, float alpha, std::vector<Trajectory> &out) const 
{
    timerSeed_.Start(false);
    // sort hits along R, since the layer number sorting is not correct in the endcaps
    typedef std::pair<int,const TrackingRecHit *> HitP;
    std::vector<std::pair<float,HitP>> hitsSorted; std::list<SiStripRecHit2D> splittedHits2D; std::list<SiStripRecHit1D> splittedHits1D; 
    int minlayer = 20, maxlayer = 0;
    for (unsigned int il = 0; il < layerhits.size(); ++il) {
        if (layerhits[il] == -1) continue;
        if (hits[layerhits[il]].layer > maxlayer) maxlayer = hits[layerhits[il]].layer;
        if (hits[layerhits[il]].layer < minlayer) minlayer = hits[layerhits[il]].layer;
        if (splitSeedHits_ && (typeid(*hits[layerhits[il]].hit) == typeid(SiStripMatchedRecHit2D))) {
            const SiStripMatchedRecHit2D &matched = dynamic_cast<const SiStripMatchedRecHit2D &>(*hits[layerhits[il]].hit);
            const GeomDet *  monoDet = tracker_->idToDet(DetId(matched.monoId()));
            const GeomDet *stereoDet = tracker_->idToDet(DetId(matched.stereoId()));
            float epsilon = 0; 
            if (dynamic_cast<const GeomDetUnit *>(monoDet)->type().isBarrel()) {
                epsilon = monoDet->position().perp() - stereoDet->position().perp();
                splittedHits1D.push_back(SiStripRecHit1D(LocalPoint(),LocalError(),DetId(matched.monoId()),matched.monoClusterRef()));
                hitsSorted.push_back(std::make_pair(hypot(hits[layerhits[il]].rho, hits[layerhits[il]].z - hitsHiRes_->z0())+epsilon, HitP(layerhits[il],&splittedHits1D.back())));
                splittedHits1D.push_back(SiStripRecHit1D(LocalPoint(),LocalError(),DetId(matched.stereoId()),matched.stereoClusterRef()));
                hitsSorted.push_back(std::make_pair(hypot(hits[layerhits[il]].rho, hits[layerhits[il]].z - hitsHiRes_->z0())-epsilon, HitP(layerhits[il],&splittedHits1D.back())));
            } else {
                epsilon = monoDet->position().z() - stereoDet->position().z();
                if (eta0 < 0) epsilon = -epsilon;
                splittedHits2D.push_back(matched.monoHit());
                hitsSorted.push_back(std::make_pair(hypot(hits[layerhits[il]].rho, hits[layerhits[il]].z - hitsHiRes_->z0())+epsilon, HitP(layerhits[il],&splittedHits2D.back())));
                splittedHits2D.push_back(matched.stereoHit());
                hitsSorted.push_back(std::make_pair(hypot(hits[layerhits[il]].rho, hits[layerhits[il]].z - hitsHiRes_->z0())-epsilon, HitP(layerhits[il],&splittedHits2D.back())));
            }
        } else {
            hitsSorted.push_back(std::make_pair(hypot(hits[layerhits[il]].rho, hits[layerhits[il]].z - hitsHiRes_->z0()), HitP(layerhits[il],hits[layerhits[il]].hit)));
        }
    }
    std::sort(hitsSorted.begin(), hitsSorted.end());

    float seedlayer = 0.5*(hits[seedhits.first].layer + hits[seedhits.second].layer);
    DEBUG2_printf("\tLayers: %d - %d  (avg seed layer: %4.1f)\n", minlayer,maxlayer,seedlayer);

    PropagationDirection pdir = alongMomentum;
    if (seedlayer >= 7 && ((maxlayer-seedlayer) < (seedlayer-minlayer))) {
        pdir = oppositeToMomentum;
        std::reverse(hitsSorted.begin(), hitsSorted.end());
    }
    std::vector<const TrackingRecHit *> vhits(hitsSorted.size());
    for (unsigned int i = 0, n = hitsSorted.size(); i < n; ++i) {
        vhits[i] = hitsSorted[i].second.second;
    }

    TrajectoryStateOnSurface tsos = startingState(eta0, phi0, alpha, hitsSorted.front().second.second);
    timerSeed_.Stop();
    timerTraj_.Start(false);
    out.clear();
    simpleTB_->run(vhits, tsos, pdir, out);
    timerTraj_.Stop();
}



#if 0
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
#endif

void
TrackCandidateBuilderFromCluster::makeTrajectories(const TrajectorySeed &seed, std::vector<Trajectory> &out) const 
{
    timerTraj_.Start(false);
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
   
    if (popOneLayer_) { 
        for (Trajectory &t : out) { if (t.chiSquared() == 0) t.invalidate(); }
        out.erase(std::remove_if(out.begin(),out.end(), std::not1(std::mem_fun_ref(&Trajectory::isValid))), out.end());
        DEBUG3_printf("\t created %lu trajectories in-out (after removing non-fitted ones) \n", out.size());
    }
    if (out.empty()) {
        timerTraj_.Stop();
        return;
    }

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
        if (DEBUG>3) {
            DEBUG4_printf("\t Somehow we lost the trajectory in the backwards fit. trying to understand why\n");
            std::vector<Trajectory> dummy;
            for (Trajectory &traj : safe) {
                recoverBackwardsFit(traj, dummy, 1);
            }
        } 
    } 

    timerTraj_.Stop();
    /*
    std::vector<Trajectory> more;
    for (Trajectory &traj : out) {
        if (traj.chiSquared() == 0 && traj.foundHits() >= int(minHitsIfNoCkf_)) {
            if (recoverBackwardsFit(traj, more, 0)) {
                traj.invalidate();
            }
        }
    }
    out.insert(out.end(), more.begin(), more.end());
    */
}

bool 
TrackCandidateBuilderFromCluster::recoverBackwardsFit(const Trajectory &traj, std::vector<Trajectory> &out, int skipHits) const 
{
    Trajectory refit(traj.sharedSeed(), traj.direction()); 
    refit.reverse();
    Trajectory::DataContainer const & data = traj.measurements();
    refit.reserve(data.size());

    // determine where to start from, skipping the last skipHits hits
    TrajectoryStateOnSurface state; int lastState = -1;
    for (lastState = data.size()-1-skipHits; lastState >= 0; --lastState) {
        if (data[lastState].updatedState().isValid() && data[lastState].recHitR().isValid()) {
            state = data[lastState].updatedState();
            state.rescaleError(100.);
            DEBUG3_printf("\t\tstarting state: %d on detid %10d/%7d\n", lastState, data[lastState].recHit()->geographicalId().rawId(), hitid(data[lastState].recHit()->hit()));
            break;
        }
    }
    if (lastState == -1) { DEBUG3_printf("\t\tcouldn't find a state to restart the refit\n"); return false; }

    KFUpdator updator;
    const Propagator * propagator = &* (traj.direction() == alongMomentum ? propagatorOpposite_ : propagator_);
    unsigned int nfail = 0;
    for (int istate = lastState; istate >= 0; --istate) {
        TransientTrackingRecHit::ConstRecHitPointer hit = data[istate].recHit();
        if (hit->det() == 0) {
            DEBUG3_printf("\t\tstate with no det: no idea what to do\n"); 
            continue;
        }
        DEBUG3_printf("\t\t\tstarting state: rho = %6.2f, z = %+7.2f phi = %+5.3f, eta = %+5.3f\n", state.globalPosition().perp(), state.globalPosition().z(), float(state.globalPosition().phi()), state.globalPosition().eta());
        TrajectoryStateOnSurface tsos = propagator->propagate(state, hit->det()->surface());
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
    if (refit.measurements().size() < 2) return false;
    if (refit.direction() != refit.seed().direction()) {
        refit.reverse();
    }
    DEBUG3_printf("\tManually refitted trajectory:\n");
    if (DEBUG>=3) dumpTraj(refit);
    out.push_back(std::move(refit));
    return true;
}

int 
TrackCandidateBuilderFromCluster::qualityCut(const Trajectory &traj) const 
{
    if (!traj.isValid() || traj.measurements().size() < 3) return -1;
    unsigned int pixelhits = 0, validhits = 0;
    for (const TrajectoryMeasurement &m : traj.measurements()) { 
        if (m.recHit()->isValid()) {
            validhits++;
            if (m.recHit()->geographicalId().subdetId() <= 2) pixelhits++;
        }
    }
    if (validhits < minHits_ && (!keepPixelTriplets_ || pixelhits < 3)) return -1;
    if (traj.chiSquared() == 0) return (validhits < minHitsIfNoCkf_ ? -1 : 0);
    if (validhits + pixelhits > 8) return 2;
    return 1; 
}

void
TrackCandidateBuilderFromCluster::saveCandidates(std::vector<Trajectory> &cands,  TrackCandidateCollection & tcCollection) const 
{
    for (Trajectory & traj : cands) {
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
            if (DEBUG>=2) dumpTraj(traj);
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
        if (typeid(*hits[i].hit) == typeid(SiStripMatchedRecHit2D)) {
            const SiStripMatchedRecHit2D &matched = static_cast<const SiStripMatchedRecHit2D &>(*hits[i].hit);
            index.push_back(std::make_pair(matched.monoId(), i));
            index.push_back(std::make_pair(matched.stereoId(), i));
        }
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
    const reco::Track *match = (DEBUG > 0 ? HTDebugger::matchCandidate(traj) : 0); 
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
            const char *matchstr[3] = { "(not good)", "(mixed)", "(good)" };
            float phipos = tm.recHit()->globalPosition().phi(); if (phipos < 0) phipos += 2*M_PI;
            printf("hit (eta %+5.3f, phi %+5.3f, detid %10d/%6d, chi2 %5.1f) %s", tm.recHit()->globalPosition().eta(), phipos, tm.recHit()->geographicalId().rawId(),hitid(tm.recHit()->hit()),tm.estimate(),
                   match == 0 ? "" : matchstr[HTDebugger::isMatchedToTrack(tm.recHit()->hit(), match)]);
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
