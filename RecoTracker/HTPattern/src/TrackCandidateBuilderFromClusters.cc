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
#include "Math/GenVector/CylindricalEta3D.h"


namespace {
    unsigned int mask2layer(unsigned int layermask) {
        unsigned int layer = 0;
        while ((layermask & 1) == 0) { layermask >>= 1; layer++; }
        return layer;
    }
};

TrackCandidateBuilderFromCluster::TrackCandidateBuilderFromCluster(const edm::ParameterSet &iConfig) :
    propagatorLabel_(iConfig.getParameter<std::string>("propagator")),
    hitBuilderLabel_(iConfig.getParameter<std::string>("TTRHBuilder")),
    dCut_(iConfig.getParameter<double>("dist2dCut")),
    dcCut_(iConfig.getParameter<double>("dist2dCorrCut")),
    minHits_(iConfig.getParameter<uint32_t>("minHits")),
    startingCovariance_(iConfig.getParameter<std::vector<double>>("startingCovariance"))
{
}

void 
TrackCandidateBuilderFromCluster::init(const edm::EventSetup& es, const SeedComparitor *ifilter)
{
    filter = ifilter;
    // get tracker
    es.get<TrackerDigiGeometryRecord>().get(tracker_);
    // get propagator
    es.get<TrackingComponentsRecord>().get(propagatorLabel_, propagator_);
    //
    es.get<TransientRecHitRecord>().get(hitBuilderLabel_,hitBuilder_);
    
    es.get<IdealMagneticFieldRecord>().get(bfield_);
}

void
TrackCandidateBuilderFromCluster::run(const HTCluster &cluster, TrajectorySeedCollection & seedCollection, TrajectorySeedCollection *seedsFromAllClusters) 
{
    printf("\n --- Cluster with %d seed layers, %d hi-res layers, %d layers --- \n", cluster.nseedlayers(), cluster.nlayers(), cluster.nmorelayers());
    std::vector<ClusteredHit> hits; hits.reserve(50);
    const HTCell &seed = mapHiRes_->get(cluster.ieta(), cluster.iphi());
    for (int hit : seed.hits()) {
        if ((*maskHiRes_)[hit]) continue;
        hits.push_back(ClusteredHit(*hitsHiRes_, hit, mask2layer(hitsHiRes_->layermask(hit))));
    }
    unsigned int nseedhits = hits.size();
    std::array<const HTCell *, 8> cells; unsigned int ncells;
    mapHiRes_->getNeighbours(cluster.ieta(), cluster.iphi(), cells, ncells);
    for (unsigned int ic = 0; ic < ncells; ++ic) {
        for (int hit : cells[ic]->hits()) {
            if ((*maskHiRes_)[hit]) continue;
            hits.push_back(ClusteredHit(*hitsHiRes_, hit, mask2layer(hitsHiRes_->layermask(hit))));
        }
    }
    unsigned int nhireshits = hits.size();
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
        if (!found) hits.push_back(ClusteredHit(*hitsLowRes_, hit, mask2layer(hitsLowRes_->layermask(hit)), 0.2));
    }
    unsigned int nhits = hits.size();
    printf("found %d seed hits, %d hi-res hits, %d hits\n", nseedhits, nhireshits, nhits);
    if (seedsFromAllClusters) dumpAsSeed(hits, *seedsFromAllClusters);

    typedef std::pair<uint16_t,uint16_t> hitpair;
    typedef std::pair<float, hitpair> distpair;
    std::vector<distpair> pairs;
    for (unsigned int i1 = 0; i1 < nseedhits; ++i1) {
        for (unsigned int i2 = i1+1; i2 < nhits; ++i2) {
            if (hits[i2].layer == hits[i1].layer) continue; // no same-layer pairs in seeding
            if (std::abs(hits[i1].layer - hits[i2].layer) <= 2) { // to too-far-away-pairs
                float d2d = hits[i1].dist2d(hits[i2]);
                //printf("\t pair %3d, %3d: (%+5.3f,%+5.3f, ly %2d)  (%+5.3f,%+5.3f ly %2d): d2d %.4f\n", i1,i2, hits[i1].eta,hits[i1].phi,hits[i1].layer, hits[i2].eta,hits[i2].phi,hits[i2].layer, d2d);
                if (d2d < dCut_) pairs.push_back(distpair(d2d, hitpair(i1,i2)));
            }
        }
    }
    printf("\t pairs %d\n", int(pairs.size()));
    std::sort(pairs.begin(), pairs.end());
    std::vector<uint8_t> microcluster(nhits, 0);
    uint8_t microclusters = 0;
    std::array<int, 20>   layerhits;   
    std::array<float, 20> layerdist;
    for (unsigned int ip = 0, np = pairs.size(); ip < np && microclusters < 100 && pairs[ip].first < dCut_; ++ip) {
        unsigned int i1 = pairs[ip].second.first; 
        unsigned int i2 = pairs[ip].second.second; 
        if (microcluster[i1] != 0 || microcluster[i2] != 0) {
            continue;
        }
        ++microclusters;
        printf("starting microcluster %d with pair %3d, %3d: (%+5.3f,%+5.3f, ly %2d)  (%+5.3f,%+5.3f ly %2d)\n", microclusters, i1,i2, hits[i1].eta,hits[i1].phi,hits[i1].layer, hits[i2].eta,hits[i2].phi,hits[i2].layer);
        microcluster[i1] = microcluster[i1] = microclusters;
        // apply pt correction from two hits
        double alphacorr = (hits[i1].phi - hits[i2].phi)/(hits[i1].rho - hits[i2].rho);
        // apply eta correction (approximate)
        double betacorr  = (hits[i1].eta - hits[i2].eta)/(hits[i1].rho - hits[i2].rho);
        double phi0 = hits[i1].phi - alphacorr*hits[i1].phi;
        double eta0 = hits[i1].eta - betacorr*hits[i1].eta; 
        unsigned int microclusternhits = 2;
        std::fill(layerhits.begin(), layerhits.end(), -1);
        std::fill(layerdist.begin(), layerdist.end(), dcCut_);
        layerhits[hits[i1].layer] = i1;
        layerhits[hits[i2].layer] = i2;
        layerdist[hits[i1].layer] = 0.f;
        layerdist[hits[i2].layer] = 0.f;
        for (unsigned int i3 = 0; i3 < nhits; ++i3) {
            if (microcluster[i3]) continue;
            float dc = hits[i3].dist2d(eta0,phi0,alphacorr,betacorr);
            //printf("\thit %2d (%+5.3f,%+5.3f, ly %2d): d2d(0) %.4f, d2dc(0) %.4f\n",i3, hits[i3].eta,hits[i3].phi,hits[i3].layer, hits[i3].dist2d(eta0,phi0,0.f,0.f), dc);
            if (dc > layerdist[hits[i3].layer]) continue;
            if (layerhits[i3] == -1) microclusternhits++;
            layerdist[hits[i3].layer] = dc;
            layerhits[hits[i3].layer] = i3;
        }
        printf("microcluster %d (%d hits)\n", microclusters, microclusternhits);
        for (unsigned int il = 0; il < layerhits.size(); ++il) {
            if (layerhits[il] == -1) continue;
            printf("\ton layer %2d (seed: %c): %3d, dist %7.4f (rho %5.1f)\n", il, (layerhits[il] == int(i1) || layerhits[il] == int(i2) ? 'Y' : 'n'),  layerhits[il], layerdist[il], hits[layerhits[il]].rho);
            microcluster[layerhits[il]] = microclusters;
        }
        if (microclusternhits < minHits_) continue;
        // sort hits along R, since the layer number sorting is not correct in the endcaps
        std::vector<std::pair<float,const TrackingRecHit *> > hitsSorted;
        for (unsigned int il = 0; il < layerhits.size(); ++il) {
            if (layerhits[il] == -1) continue;
            hitsSorted.push_back(std::make_pair(hypot(hits[layerhits[il]].rho, hits[layerhits[il]].z - hitsHiRes_->z0()), hits[layerhits[il]].hit));
        }
        std::sort(hitsSorted.begin(), hitsSorted.end());
        // ok now we build the seed
        FreeTrajectoryState fts(startingState(eta0,phi0,hitsHiRes_->alpha()+alphacorr));
        TrajectoryStateOnSurface tsos;
        edm::OwnVector<TrackingRecHit> seedHits; seedHits.reserve(nhits);
        KFUpdator updator;
        const TrackingRecHit* lasthit = 0;
        unsigned int myhits = 0;
        for (const std::pair<float,const TrackingRecHit *> & pair: hitsSorted) {
            const TrackingRecHit* hit = pair.second;
            TrajectoryStateOnSurface state = (lasthit == 0) ?
                  propagator_->propagate(fts,  tracker_->idToDet(hit->geographicalId())->surface())
                : propagator_->propagate(tsos, tracker_->idToDet(hit->geographicalId())->surface());
            if (!state.isValid()) { printf ("\tfailed propagation\n"); break; }
            TransientTrackingRecHit::RecHitPointer tth = hitBuilder_->build(hit);
            tth = tth->clone(state);
            TrajectoryStateOnSurface updated = updator.update(state, *tth);
            if (!updated.isValid()) { printf("\tfailed update\n"); break; }
            tsos = updated;
            seedHits.push_back(tth->hit()->clone());
            lasthit = hit;
            myhits++;
            if (myhits > 1) printf("\tKF-fitted %d hits (pt = %7.2f +/- %7.2f, q = %+1d)\n", myhits, tsos.globalMomentum().perp(), TrajectoryStateAccessor(*tsos.freeState()).inversePtError()*tsos.globalMomentum().perp2(), tsos.charge());
        }
        if (myhits >= minHits_) {
            PTrajectoryStateOnDet const & PTraj = trajectoryStateTransform::persistentState(tsos, lasthit->geographicalId().rawId());
            TrajectorySeed seed(PTraj,std::move(seedHits),alongMomentum);
            seedCollection.push_back(seed);
            for (unsigned int il = 0; il < layerhits.size(); ++il) {
                if (layerhits[il] == -1) continue;
                std::vector<bool> *mask = layerhits[il] < int(nhireshits) ?  maskHiRes_ : maskLowRes_;
                (*mask)[hits[layerhits[il]].id] = true;
            }
            printf("--> saved as trajectory seed!\n");
        }
    }
}

FreeTrajectoryState
TrackCandidateBuilderFromCluster::startingState(float eta0, float phi0, float alpha) const 
{
    // alpha = 0.003 * B[T] / pT[GeV] -->
    float pt = std::abs(alpha ? 0.003f * bfield_->nominalValue() / alpha : 100.0f); 
    ROOT::Math::CylindricalEta3D<double> p3d(pt, eta0, phi0);
    GlobalPoint  x(hitsHiRes_->x0(), hitsHiRes_->y0(), hitsHiRes_->z0());
    GlobalVector p(p3d.x(), p3d.y(), p3d.z());
    GlobalTrajectoryParameters gtp(x,p, alpha > 0 ? +1 : -1, &*bfield_);
    printf("Cluster for pT %8.4f, eta %+5.3f, phi %+5.3f, q %+2d\n", pt, eta0, phi0, alpha > 0 ? +1 : -1);

    AlgebraicSymMatrix55 cov;
    for (unsigned int i = 0; i < 5; ++i) cov(i,i) = startingCovariance_[i];

    return FreeTrajectoryState(gtp, CurvilinearTrajectoryError(cov));
}

void
TrackCandidateBuilderFromCluster::dumpAsSeed(const std::vector<ClusteredHit> &cluster, TrajectorySeedCollection &seedCollection)
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
