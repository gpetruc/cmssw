#include "RecoTracker/HTPattern/interface/TrackCandidateBuilderFromCluster.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"

namespace {
    unsigned int mask2layer(unsigned int layermask) {
        unsigned int layer = 0;
        while ((layermask & 1) == 0) { layermask >>= 1; layer++; }
        return layer;
    }
};

void
TrackCandidateBuilderFromCluster::run(const HTCluster &cluster) 
{
    std::array<std::vector<const TrackingRecHit *>, 20> hits;
    std::array<std::vector<int>,                    20> hitids;
    for (auto v: hits) v.clear();

    const HTCell &seed = mapHiRes_->get(cluster.ieta(), cluster.iphi());
    for (int hit : seed.hits()) {
        //printf("[cell @%p, hit %d, layermaks %x, detid %d]\n", (void*)(&seed), hit, hitsHiRes_->layermask(hit), hitsHiRes_->hit(hit)->geographicalId().rawId());
        hitids[mask2layer(hitsHiRes_->layermask(hit))].push_back(hit);
        hits[mask2layer(hitsHiRes_->layermask(hit))].push_back(hitsHiRes_->hit(hit));
    }

    printf("\nseed hits on cluster (%2d seed layers = %2d)\n", cluster.nseedlayers(), seed.nlayers());
    for (unsigned int il = 0; il < hits.size(); ++il) {
       if (!hits[il].empty()) printf("on layer %2d: %2ld seed hits\n", il, hits[il].size());
       for (unsigned int ih = 0, nh = hits[il].size(); ih < nh; ++ih) {
            printf("\t %d (eta %+6.4f, phi %+6.4f)\n", hits[il][ih]->geographicalId().rawId(), hitsHiRes_->eta(hitids[il][ih]), hitsHiRes_->phi(hitids[il][ih]));
       }
    }

    std::array<const HTCell *, 8> cells; unsigned int ncells;
    mapHiRes_->getNeighbours(cluster.ieta(), cluster.iphi(), cells, ncells);
    for (unsigned int ic = 0; ic < ncells; ++ic) {
        for (int hit : cells[ic]->hits()) {
            hits[mask2layer(hitsHiRes_->layermask(hit))].push_back(hitsHiRes_->hit(hit));
            hitids[mask2layer(hitsHiRes_->layermask(hit))].push_back(hit);
        }
    }

    printf("\nhits after including neighbouring cells (%2d layers)\n", cluster.nlayers());
    for (unsigned int il = 0; il < hits.size(); ++il) {
       if (!hits[il].empty()) printf("on layer %2d: %2ld hits\n", il, hits[il].size());
       for (unsigned int ih = 0, nh = hits[il].size(); ih < nh; ++ih) {
            printf("\t %d (eta %+6.4f, phi %+6.4f)\n", hits[il][ih]->geographicalId().rawId(), hitsHiRes_->eta(hitids[il][ih]), hitsHiRes_->phi(hitids[il][ih]));
       }
    }
    
    const HTCell &lowres = mapLowRes_->get(cluster.ieta()>>etashift_, cluster.iphi()>>phishift_);
    for (int hit : lowres.hits()) {
        int il = mask2layer(hitsLowRes_->layermask(hit)); if (il < 3) continue;
        std::vector<const TrackingRecHit *> & myhits = hits[il];
        const TrackingRecHit * thishit = hitsLowRes_->hit(hit);
        bool found = false;
        for (const TrackingRecHit *other : myhits) {
            if ((other->geographicalId() == thishit->geographicalId() && other->sharesInput(thishit, TrackingRecHit::some))) {
                found = true; break;
            }
        }
        if (!found) myhits.push_back(thishit);
    }

    printf("\nhits on cluster, including low-resolution (%2d layers)\n", cluster.nmorelayers());
    for (unsigned int il = 0; il < hits.size(); ++il) {
       if (!hits[il].empty()) printf("on layer %2d: %2ld hits\n", il, hits[il].size());
       for (auto *h : hits[il]) printf("\t %012d\n", h->geographicalId().rawId());
    }


}
