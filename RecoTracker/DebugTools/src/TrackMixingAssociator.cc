#include "RecoTracker/DebugTools/interface/TrackMixingAssociator.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit1D.h"
#include "DataFormats/TrackerRecHit2D/interface/ProjectedSiStripRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2D.h"
#include <algorithm>

//#define DEBUG_TrackMixingAssociator
#ifdef DEBUG_TrackMixingAssociator
#define DEBUG_printf  printf
namespace {
    int hitid(const TrackingRecHit *hit)
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
}
#else
#define DEBUG_printf  if (0) printf
#define hitid(x) 0
#endif

void 
TrackMixingAssociator::registerTrackEvent(int eid, const reco::TrackCollection &tracks) 
{
    for (reco::TrackCollection::const_iterator it = tracks.begin(), ed = tracks.end(); it != ed; ++it) {
        const reco::Track &tk = *it;
        DEBUG_printf("registering track of pt %.3f, eta %+.3f\n", tk.pt(), tk.eta());
        for (trackingRecHit_iterator ithit = tk.recHitsBegin(), edhit = tk.recHitsEnd(); ithit != edhit; ++ithit) {
            const TrackingRecHit *hit = &**ithit;
            uint32_t detid = hit->geographicalId().rawId();
            DEBUG_printf("  registered hit %10d/%7d, %s\n", hit->geographicalId().rawId(), hitid(hit), typeid(*hit).name());
            allHits_[detid].push_back(HitRecord(eid, hit, &tk));
        }
    }
}

void 
TrackMixingAssociator::registerSeedEvent(int eid, const TrajectorySeedCollection &tracks) 
{
    for (TrajectorySeedCollection::const_iterator it = tracks.begin(), ed = tracks.end(); it != ed; ++it) {
        const TrajectorySeed &tk = *it;
        for (TrajectorySeed::range r = tk.recHits(); r.first != r.second; ++r.first) {
            const TrackingRecHit *hit = &*r.first;
            uint32_t detid = hit->geographicalId().rawId();
            if (typeid(*hit) == typeid(SiStripMatchedRecHit2D)) {
                const SiStripMatchedRecHit2D &mh = static_cast<const SiStripMatchedRecHit2D &>(*hit);
                allSeedHits_[mh.monoId()  ].push_back(SeedRecord(eid, &mh, &tk));
                allSeedHits_[mh.stereoId()].push_back(SeedRecord(eid, &mh, &tk));
            } else if (typeid(*hit) == typeid(ProjectedSiStripRecHit2D)) {
                const ProjectedSiStripRecHit2D &ph = static_cast<const ProjectedSiStripRecHit2D &>(*hit);
                allSeedHits_[ph.originalHit().geographicalId().rawId()].push_back(SeedRecord(eid, &ph, &tk));
            } else {
                allSeedHits_[detid].push_back(SeedRecord(eid, hit, &tk));
            }
        }
    }
}

void 
TrackMixingAssociator::registerTrajEvent(int eid, const std::vector<Trajectory> &tracks)
{
    for (const Trajectory &tk : tracks) {
        for (const TrajectoryMeasurement &tm : tk.measurements()) {
            if (!tm.recHit()->isValid()) continue;
            const TrackingRecHit *hit = tm.recHit()->hit();
            uint32_t detid = hit->geographicalId().rawId();
            if (typeid(*hit) == typeid(SiStripMatchedRecHit2D)) {
                const SiStripMatchedRecHit2D &mh = static_cast<const SiStripMatchedRecHit2D &>(*hit);
                allTrajHits_[mh.monoId()  ].push_back(TrajRecord(eid, &mh, &tk));
                allTrajHits_[mh.stereoId()].push_back(TrajRecord(eid, &mh, &tk));
            } else if (typeid(*hit) == typeid(ProjectedSiStripRecHit2D)) {
                const ProjectedSiStripRecHit2D &ph = static_cast<const ProjectedSiStripRecHit2D &>(*hit);
                allTrajHits_[ph.originalHit().geographicalId().rawId()].push_back(TrajRecord(eid, &ph, &tk));
            } else {
                allTrajHits_[detid].push_back(TrajRecord(eid, hit, &tk));
            }
        }
    }
}

void 
TrackMixingAssociator::registerTrackCandEvent(int eid, const std::vector<TrackCandidate> &tracks)
{
    for (const TrackCandidate &tk : tracks) {
        for (TrajectorySeed::range r = tk.recHits(); r.first != r.second; ++r.first) {
            const TrackingRecHit *hit = &*r.first;
            uint32_t detid = hit->geographicalId().rawId();
            if (typeid(*hit) == typeid(SiStripMatchedRecHit2D)) {
                const SiStripMatchedRecHit2D &mh = static_cast<const SiStripMatchedRecHit2D &>(*hit);
                allTrackCandHits_[mh.monoId()  ].push_back(TrackCandRecord(eid, &mh, &tk));
                allTrackCandHits_[mh.stereoId()].push_back(TrackCandRecord(eid, &mh, &tk));
            } else if (typeid(*hit) == typeid(ProjectedSiStripRecHit2D)) {
                const ProjectedSiStripRecHit2D &ph = static_cast<const ProjectedSiStripRecHit2D &>(*hit);
                allTrackCandHits_[ph.originalHit().geographicalId().rawId()].push_back(TrackCandRecord(eid, &ph, &tk));
            } else {
                allTrackCandHits_[detid].push_back(TrackCandRecord(eid, hit, &tk));
            }
        }
    }
}




void 
TrackMixingAssociator::registerClusterEvent(int id, const SiPixelClusterCollection &allPixels, const SiStripClusterCollection &allStrips) 
{
   registerClusters(id, allStrips); 
   registerClusters(id, allPixels); 
}


template<typename T>
void
TrackMixingAssociator::associateHitToRec(const TrackingRecHit &hit, std::vector<RecAssociation<T> > &out, boost::unordered_map<uint32_t, std::vector<RecRecord<T> > > &record) {
    const std::vector<RecRecord<T> > & hitRecords = record[ hit.geographicalId().rawId() ];
    DEBUG_printf("   associateHitToRec(%10d/%7d, %s): %d records\n", hit.geographicalId().rawId(), hitid(&hit), typeid(hit).name(), int(hitRecords.size()));
    for (typename std::vector<RecRecord<T> >::const_iterator itr = hitRecords.begin(), edr = hitRecords.end(); itr != edr; ++itr) {
        if (matchesRecord(hit, *itr)) {
            DEBUG_printf("      matched, track @%p!\n", (const void*)itr->track);
            typename std::vector<RecAssociation<T> >::iterator as = std::find(out.begin(), out.end(), itr->track);
            if (as == out.end()) { out.push_back(RecAssociation<T>(itr->eventId, itr->track)); as = out.end() - 1; }
            as->addSharedHit(hit);
        }
    }
}

template<typename T>
void
TrackMixingAssociator::associateToRec(const reco::Track &tk, std::vector<RecAssociation<T> > &out, boost::unordered_map<uint32_t, std::vector<RecRecord<T> > > &record) 
{
    out.clear();
    for (trackingRecHit_iterator ithit = tk.recHitsBegin(), edhit = tk.recHitsEnd(); ithit != edhit; ++ithit) {
        const TrackingRecHit &hit = **ithit;
        associateHitToRec<T>(hit,out,record);
        
    }
    std::sort(out.begin(), out.end());
}

template<typename T>
void
TrackMixingAssociator::associateToRec(const TrajectorySeed &tk, std::vector<RecAssociation<T> > &out, boost::unordered_map<uint32_t, std::vector<RecRecord<T> > > &record) 
{
    out.clear();
    for (TrajectorySeed::range r = tk.recHits(); r.first != r.second; ++r.first) {
        const TrackingRecHit &hit = *r.first;
        DEBUG_printf("associateToRec, hit %10d/%7d, %s:\n", hit.geographicalId().rawId(), hitid(&hit), typeid(hit).name());
        if (typeid(hit) == typeid(SiStripMatchedRecHit2D)) {
            const SiStripMatchedRecHit2D &mh = static_cast<const SiStripMatchedRecHit2D &>(hit);
            associateHitToRec<T>(mh.monoHit(),out,record);
            associateHitToRec<T>(mh.stereoHit(),out,record);
        } else if (typeid(hit) == typeid(ProjectedSiStripRecHit2D)) {
            const ProjectedSiStripRecHit2D &ph = static_cast<const ProjectedSiStripRecHit2D &>(hit);
            associateHitToRec<T>(ph.originalHit(),out,record);
        } else {
            associateHitToRec<T>(hit,out,record);
        }
    }
    std::sort(out.begin(), out.end());
}

template<typename T>
void
TrackMixingAssociator::associateToRec(const TrackCandidate &tk, std::vector<RecAssociation<T> > &out, boost::unordered_map<uint32_t, std::vector<RecRecord<T> > > &record) 
{
    out.clear();
    for (TrackCandidate::range r = tk.recHits(); r.first != r.second; ++r.first) {
        const TrackingRecHit &hit = *r.first;
        DEBUG_printf("associateToRec, hit %10d/%7d, %s:\n", hit.geographicalId().rawId(), hitid(&hit), typeid(hit).name());
        if (typeid(hit) == typeid(SiStripMatchedRecHit2D)) {
            const SiStripMatchedRecHit2D &mh = static_cast<const SiStripMatchedRecHit2D &>(hit);
            associateHitToRec<T>(mh.monoHit(),out,record);
            associateHitToRec<T>(mh.stereoHit(),out,record);
        } else if (typeid(hit) == typeid(ProjectedSiStripRecHit2D)) {
            const ProjectedSiStripRecHit2D &ph = static_cast<const ProjectedSiStripRecHit2D &>(hit);
            associateHitToRec<T>(ph.originalHit(),out,record);
        } else {
            associateHitToRec<T>(hit,out,record);
        }
    }
    std::sort(out.begin(), out.end());
}



template<typename T>
void
TrackMixingAssociator::associateToRec(const Trajectory &tk, std::vector<RecAssociation<T> > &out, boost::unordered_map<uint32_t, std::vector<RecRecord<T> > > &record) 
{
    out.clear();
    for (const TrajectoryMeasurement &tm : tk.measurements()) {
        if (!tm.recHit()->isValid()) continue;
        const TrackingRecHit &hit = *tm.recHit()->hit();
        DEBUG_printf("associateToRec, hit %10d/%7d, %s:\n", hit.geographicalId().rawId(), hitid(&hit), typeid(hit).name());
        if (typeid(hit) == typeid(SiStripMatchedRecHit2D)) {
            const SiStripMatchedRecHit2D &mh = static_cast<const SiStripMatchedRecHit2D &>(hit);
            associateHitToRec<T>(mh.monoHit(),out,record);
            associateHitToRec<T>(mh.stereoHit(),out,record);
        } else if (typeid(hit) == typeid(ProjectedSiStripRecHit2D)) {
            const ProjectedSiStripRecHit2D &ph = static_cast<const ProjectedSiStripRecHit2D &>(hit);
            associateHitToRec<T>(ph.originalHit(),out,record);
        } else {
            associateHitToRec<T>(hit,out,record);
        }
    }
    std::sort(out.begin(), out.end());
}

template<typename T>
void
TrackMixingAssociator::associateToRec(const std::vector<const TrackingRecHit *> &tk, std::vector<RecAssociation<T> > &out, boost::unordered_map<uint32_t, std::vector<RecRecord<T> > > &record) 
{
    out.clear();
    for (const TrackingRecHit *phit : tk) {
        const TrackingRecHit &hit = *phit;
        DEBUG_printf("associateToRec, hit %10d/%7d, %s:\n", hit.geographicalId().rawId(), hitid(&hit), typeid(hit).name());
        if (typeid(hit) == typeid(SiStripMatchedRecHit2D)) {
            const SiStripMatchedRecHit2D &mh = static_cast<const SiStripMatchedRecHit2D &>(hit);
            associateHitToRec<T>(mh.monoHit(),out,record);
            associateHitToRec<T>(mh.stereoHit(),out,record);
        } else if (typeid(hit) == typeid(ProjectedSiStripRecHit2D)) {
            const ProjectedSiStripRecHit2D &ph = static_cast<const ProjectedSiStripRecHit2D &>(hit);
            associateHitToRec<T>(ph.originalHit(),out,record);
        } else {
            associateHitToRec<T>(hit,out,record);
        }
    }
    std::sort(out.begin(), out.end());
}

template<typename T>
void
TrackMixingAssociator::associateToRec(const TrackingRecHit &hit, std::vector<RecAssociation<T> > &out, boost::unordered_map<uint32_t, std::vector<RecRecord<T> > > &record) 
{
    out.clear();
    DEBUG_printf("associateToRec, hit %10d/%7d, %s:\n", hit.geographicalId().rawId(), hitid(&hit), typeid(hit).name());
    if (typeid(hit) == typeid(SiStripMatchedRecHit2D)) {
        const SiStripMatchedRecHit2D &mh = static_cast<const SiStripMatchedRecHit2D &>(hit);
        associateHitToRec<T>(mh.monoHit(),out,record);
        associateHitToRec<T>(mh.stereoHit(),out,record);
    } else if (typeid(hit) == typeid(ProjectedSiStripRecHit2D)) {
        const ProjectedSiStripRecHit2D &ph = static_cast<const ProjectedSiStripRecHit2D &>(hit);
        associateHitToRec<T>(ph.originalHit(),out,record);
    } else {
        associateHitToRec<T>(hit,out,record);
    }
    std::sort(out.begin(), out.end());
}



void
TrackMixingAssociator::associateToTracks(const reco::Track &tk, std::vector<TrackAssociation> &out) {
    return associateToRec<reco::Track>(tk,out,allHits_); 
}

void
TrackMixingAssociator::associateToTracks(const TrajectorySeed &tk, std::vector<TrackAssociation> &out) { 
    return associateToRec<reco::Track>(tk,out,allHits_);
}

void
TrackMixingAssociator::associateToTracks(const Trajectory &tk, std::vector<TrackAssociation> &out) { 
    return associateToRec<reco::Track>(tk,out,allHits_);
}

void
TrackMixingAssociator::associateToTracks(const TrackCandidate &tk, std::vector<TrackAssociation> &out) { 
    return associateToRec<reco::Track>(tk,out,allHits_);
}

void
TrackMixingAssociator::associateToTracks(const std::vector<const TrackingRecHit *> &tk, std::vector<TrackAssociation> &out) {
    return associateToRec<reco::Track>(tk,out,allHits_);
}

void
TrackMixingAssociator::associateToTracks(const TrackingRecHit &tk, std::vector<TrackAssociation> &out) {
    return associateToRec<reco::Track>(tk,out,allHits_);
}

void
TrackMixingAssociator::associateToSeed(const reco::Track &tk, std::vector<SeedAssociation> &out) { 
    return associateToRec<TrajectorySeed>(tk,out,allSeedHits_); 
}

void
TrackMixingAssociator::associateToSeed(const TrajectorySeed &tk, std::vector<SeedAssociation> &out) { 
    return associateToRec<TrajectorySeed>(tk,out,allSeedHits_); 
} 

void
TrackMixingAssociator::associateToTrajs(const reco::Track &tk, std::vector<TrajAssociation> &out) {
    return associateToRec<Trajectory>(tk,out,allTrajHits_); 
}

void
TrackMixingAssociator::associateToTrajs(const TrajectorySeed &tk, std::vector<TrajAssociation> &out) { 
    return associateToRec<Trajectory>(tk,out,allTrajHits_);
}

void
TrackMixingAssociator::associateToTrajs(const Trajectory &tk, std::vector<TrajAssociation> &out) { 
    return associateToRec<Trajectory>(tk,out,allTrajHits_);
}

void
TrackMixingAssociator::associateToTrajs(const std::vector<const TrackingRecHit *> &tk, std::vector<TrajAssociation> &out) {
    return associateToRec<Trajectory>(tk,out,allTrajHits_);
}


void
TrackMixingAssociator::associateToTrackCands(const reco::Track &tk, std::vector<TrackCandAssociation> &out) {
    return associateToRec<TrackCandidate>(tk,out,allTrackCandHits_); 
}

void
TrackMixingAssociator::associateToTrackCands(const TrajectorySeed &tk, std::vector<TrackCandAssociation> &out) { 
    return associateToRec<TrackCandidate>(tk,out,allTrackCandHits_);
}

void
TrackMixingAssociator::associateToTrackCands(const Trajectory &tk, std::vector<TrackCandAssociation> &out) { 
    return associateToRec<TrackCandidate>(tk,out,allTrackCandHits_);
}

void
TrackMixingAssociator::associateToTrackCands(const std::vector<const TrackingRecHit *> &tk, std::vector<TrackCandAssociation> &out) {
    return associateToRec<TrackCandidate>(tk,out,allTrackCandHits_);
}

void 
TrackMixingAssociator::associateHitToClusters(const TrackingRecHit &hit, std::vector<TrackClusterAssociation> &out) 
{
    const std::vector<ClusterRecord> & clusterRecords = allClusters_[ hit.geographicalId().rawId() ];
    std::vector<int> eventIds;
    int nmatches = 0;
    for (std::vector<ClusterRecord>::const_iterator itr = clusterRecords.begin(), edr = clusterRecords.end(); itr != edr; ++itr) {
        if (matches(hit, *itr)) { nmatches++; eventIds.push_back(itr->eventId); }
    }
    for (std::vector<int>::const_iterator iti = eventIds.begin(), edi = eventIds.end(); iti != edi; ++iti) {
        std::vector<TrackClusterAssociation>::iterator as = std::find(out.begin(), out.end(), *iti);
        if (as == out.end()) { out.push_back(TrackClusterAssociation(*iti)); as = out.end() - 1; }
        as->addHit(hit, nmatches == 1);
    }
}

void 
TrackMixingAssociator::associateToClusters(const reco::Track &tk, std::vector<TrackClusterAssociation> &out) 
{
    out.clear();
    for (trackingRecHit_iterator ithit = tk.recHitsBegin(), edhit = tk.recHitsEnd(); ithit != edhit; ++ithit) {
        const TrackingRecHit &hit = **ithit;
        associateHitToClusters(hit,out);
    }
    std::sort(out.begin(), out.end());
}


void 
TrackMixingAssociator::associateToClusters(const TrajectorySeed &tk, std::vector<TrackClusterAssociation> &out) 
{
    out.clear();
    for (TrajectorySeed::range r = tk.recHits(); r.first != r.second; ++r.first) {
        const TrackingRecHit &hit = *r.first;
        if (typeid(hit) == typeid(SiStripMatchedRecHit2D)) {
            const SiStripMatchedRecHit2D &mh = static_cast<const SiStripMatchedRecHit2D &>(hit);
            associateHitToClusters(mh.monoHit(),out);
            associateHitToClusters(mh.stereoHit(),out);
        } else if (typeid(hit) == typeid(ProjectedSiStripRecHit2D)) {
            const ProjectedSiStripRecHit2D &ph = static_cast<const ProjectedSiStripRecHit2D &>(hit);
            associateHitToClusters(ph.originalHit(),out);
        } else {
            associateHitToClusters(hit,out);
        }
    }
    std::sort(out.begin(), out.end());
}


bool
TrackMixingAssociator::matches(const SiStripCluster &c1, const SiStripCluster &c2) 
{
    // either: 
    //   - c1 begins after c2 begins but before c2 ends
    //   - c2 begins after c1 begins but before c1 ends
    // the rightmost beginning must be on the left of the leftmost ending
    return std::max(c1.firstStrip(),c2.firstStrip()) <= std::min(c1.firstStrip()+c1.amplitudes().size(),c2.firstStrip()+c2.amplitudes().size());
}

bool
TrackMixingAssociator::matches(const SiPixelCluster &c1, const SiPixelCluster &c2)
{
    return ( std::max(c1.minPixelRow(), c2.minPixelRow()) <= std::min(c1.maxPixelRow(), c2.maxPixelRow()) ) && 
           ( std::max(c1.minPixelCol(), c2.minPixelCol()) <= std::min(c1.maxPixelCol(), c2.maxPixelCol()) ) ;
}

template<typename T, typename Record>
bool TrackMixingAssociator::matchesRecordSpecific(const TrackingRecHit &hit, const Record &record) 
{
    if (typeid(hit) == typeid(T) && typeid(T) == typeid(*record.hit)) {
        return matches(*(static_cast<const T &>(hit)).cluster(),
                       *(static_cast<const T &>(*record.hit)).cluster());
    }
    return false;
}
template<typename Record>
bool TrackMixingAssociator::matchesRecordTracker1D(const TrackingRecHit &hit, const Record &record) 
{
    const TrackerSingleRecHit *sihit = dynamic_cast<const TrackerSingleRecHit *>(&hit);
    if (sihit == 0) return false;
    const TrackerSingleRecHit *rshit = dynamic_cast<const TrackerSingleRecHit *>(record.hit);
    if (rshit != 0) {
        return matches(sihit->stripCluster(), rshit->stripCluster());
    }
    return false;
}



template<typename Record>
bool 
TrackMixingAssociator::matchesRecord(const TrackingRecHit &hit, const Record &record) 
{
    if (typeid(hit) != typeid(SiStripMatchedRecHit2D)) {
        return typeid(hit) == typeid(SiPixelRecHit) ? 
                matchesRecordSpecific<SiPixelRecHit,Record>(hit,record) :
                matchesRecordTracker1D<Record>(hit,record);
    } else {
        const SiStripMatchedRecHit2D &mhit = static_cast<const SiStripMatchedRecHit2D &>(hit);
        return matchesRecordTracker1D<Record>(mhit.monoHit(),  record) ||
               matchesRecordTracker1D<Record>(mhit.stereoHit(),record);
    }
}

bool
TrackMixingAssociator::matches(const TrackingRecHit &hit, const ClusterRecord &record) 
{
    if (record.stripCluster != 0 && typeid(hit) == typeid(SiStripRecHit2D)) {
         return matches( *(static_cast<const SiStripRecHit2D &>(hit)).cluster(), *record.stripCluster);
    } else if (record.stripCluster != 0 && typeid(hit) == typeid(SiStripRecHit1D)) {
         return matches( *(static_cast<const SiStripRecHit1D &>(hit)).cluster(), *record.stripCluster);
    } else if (record.stripCluster != 0 && typeid(hit) == typeid(SiStripMatchedRecHit2D)) {
         return matches( (static_cast<const SiStripMatchedRecHit2D &>(hit)).monoCluster(),   *record.stripCluster) ||
                matches( (static_cast<const SiStripMatchedRecHit2D &>(hit)).stereoCluster(), *record.stripCluster);
    } else if (record.pixelCluster != 0 && typeid(hit) == typeid(SiPixelRecHit)) {
         return matches( *(static_cast<const SiPixelRecHit &>(hit)).cluster(), *record.pixelCluster);
    }
    return false;
}

template<typename T>
void TrackMixingAssociator::registerClusters(int eventId, const edmNew::DetSetVector<T> &clusters)
{
    for (typename edmNew::DetSetVector<T>::const_iterator itds = clusters.begin(), edds = clusters.end(); itds != edds; ++itds) {
        typename edmNew::DetSet<T> ds = *itds;
        std::vector<ClusterRecord> & theseClusters = allClusters_[ds.detId()];
        for (typename edmNew::DetSet<T>::const_iterator it = ds.begin(), ed = ds.end(); it != ed; ++it) {
            theseClusters.push_back(ClusterRecord(eventId, *it));
        }
    }
}

void TrackMixingAssociator::getAssociatedClusters(const TrackingRecHit &hit, std::vector<std::pair<int, const SiStripCluster *>> &out) 
{
    if (typeid(hit) == typeid(SiStripMatchedRecHit2D)) {
        std::vector<std::pair<int, const SiStripCluster *>> o1, o2;
        const SiStripMatchedRecHit2D &mh = static_cast<const SiStripMatchedRecHit2D &>(hit);
        getAssociatedClusters(mh.monoHit(),o1);
        getAssociatedClusters(mh.stereoHit(),o2);
        out.swap(o1);
        out.insert(out.end(), o2.begin(), o2.end());
    } else if (typeid(hit) == typeid(ProjectedSiStripRecHit2D)) {
        const ProjectedSiStripRecHit2D &ph = static_cast<const ProjectedSiStripRecHit2D &>(hit);
        getAssociatedClusters(ph.originalHit(),out);
    } else {
        const TrackerSingleRecHit *sihit = dynamic_cast<const TrackerSingleRecHit *>(&hit);
        assert(sihit != 0);
        getAssociatedClusters(hit.geographicalId().rawId(), sihit->stripCluster(), out);
    }
}

void TrackMixingAssociator::getAssociatedClusters(const TrackingRecHit &hit, std::vector<std::pair<int, const SiPixelCluster *>> &out) 
{
    assert(typeid(hit) == typeid(SiPixelRecHit));
    getAssociatedClusters(hit.geographicalId().rawId(), *(static_cast<const SiPixelRecHit &>(hit)).cluster(), out);
}



void TrackMixingAssociator::getAssociatedClusters(uint32_t detid, const SiStripCluster &hit, std::vector<std::pair<int, const SiStripCluster *>> &out)
{
    out.clear();
    const std::vector<ClusterRecord> & clusterRecords = allClusters_[ detid ];
    for (std::vector<ClusterRecord>::const_iterator itr = clusterRecords.begin(), edr = clusterRecords.end(); itr != edr; ++itr) {
        if (itr->stripCluster != 0 && matches(hit, *itr->stripCluster)) {
            out.push_back(std::make_pair(itr->eventId, itr->stripCluster));
        }
    }
}

void TrackMixingAssociator::getAssociatedClusters(uint32_t detid, const SiPixelCluster &hit, std::vector<std::pair<int, const SiPixelCluster *>> &out)
{
    out.clear();
    const std::vector<ClusterRecord> & clusterRecords = allClusters_[ detid ];
    for (std::vector<ClusterRecord>::const_iterator itr = clusterRecords.begin(), edr = clusterRecords.end(); itr != edr; ++itr) {
        if (itr->pixelCluster != 0 && matches(hit, *itr->pixelCluster)) {
            out.push_back(std::make_pair(itr->eventId, itr->pixelCluster));
        }
    }
}


