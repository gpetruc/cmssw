#ifndef UserCode_TrackerTrackMixing_interface_TrackMixingAssociator_h
#define UserCode_TrackerTrackMixing_interface_TrackMixingAssociator_h

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidate.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"

#include <boost/unordered_map.hpp>

class TrackMixingAssociator {
public:
    typedef edmNew::DetSetVector<SiStripCluster> SiStripClusterCollection;
    typedef edmNew::DetSetVector<SiPixelCluster> SiPixelClusterCollection;

    template<typename T>
    struct RecAssociation {
        /// Simple structure to hold the association between one track and another track

        /// Id of the event (as provided in 'registreTrackEvent')
        int eventId;

        /// Pointer to the other track
        const T *track;

        /// Hits shared
        unsigned int sharedHits;

        /// HitPattern filled with shared hits, so that their location in the tracker is also available
        reco::HitPattern sharedHitPattern;

        /// Empty constructor, then one should add the shared hits one by one.
        RecAssociation(int eventId1, const T *track1, int sharedHits1=0) : eventId(eventId1), track(track1), sharedHits(sharedHits1) {}

        /// Include one shared hit
        void addSharedHit(const TrackingRecHit &hit) {
            sharedHitPattern.set(hit, sharedHits);
            sharedHits++;
        }

        /// To sort the matches by relevance. Sorts by shared hits (desc), then by event id (asc), then 'randomly'
        bool operator<(const RecAssociation<T> &other) const { 
            if (sharedHits != other.sharedHits) return sharedHits > other.sharedHits;
            if (eventId    != other.eventId)    return eventId    < other.eventId;
            return track < other.track;
        }
        bool operator==(const RecAssociation<T> &other)      const { return (other.track == track); }
        bool operator==(const T           *otherTrack) const { return (otherTrack  == track); }
    };
    typedef RecAssociation<reco::Track> TrackAssociation;
    typedef RecAssociation<TrajectorySeed> SeedAssociation;
    typedef RecAssociation<Trajectory> TrajAssociation;
    typedef RecAssociation<TrackCandidate> TrackCandAssociation;

    struct TrackClusterAssociation {
        /** Association between a track and a collection of clusters
            It contains:
             - the event id of the cluster collection (as per 'registerClusterEvent')
             - the total number and hitpattern for the hits from a given collection
             - the total number and hitpattern for the exclusive hits from a given collection 
               (hits that don't have overlap with any other cluster collection)
             - the total number and hitpattern for the shared hits from a given collection
               (hits that overlap also with other cluster collections beyond this one) 
         */

        /// Id of the event (as provided in 'registreTrackEvent')
        int eventId;

        unsigned int exclusiveHits, sharedHits, allHits;
        reco::HitPattern   exclusiveHitPattern, sharedHitPattern, allHitPattern;

        /// Create the empty association, then add the hits
        TrackClusterAssociation(int evId) : eventId(evId), 
                                            exclusiveHits(0), sharedHits(0), allHits(0), 
                                            exclusiveHitPattern(), sharedHitPattern(), allHitPattern() { }

        /// Add one hit; if exclusive is false, it assumes the hit is shared
        void addHit(const TrackingRecHit &hit, bool exclusive) {
            allHitPattern.set(hit, allHits); allHits++;
            if (exclusive) {
                exclusiveHitPattern.set(hit, exclusiveHits); exclusiveHits++;
            } else {
                sharedHitPattern.set(hit, sharedHits); sharedHits++;
            }
        }

        /// Sort (by event id, asc)
        bool operator<(const TrackClusterAssociation &other) const { return eventId < other.eventId; }
        bool operator==(int otherEventId) const { return eventId == otherEventId; }
        bool operator==(const TrackClusterAssociation &other) const { return eventId == other.eventId; }
    };
    
    /// Register this collection of tracks coming from event id 'id'
    void registerTrackEvent(int id, const reco::TrackCollection &tracks) ;

    /// Register this collection of tracks coming from event id 'id'
    void registerSeedEvent(int id, const TrajectorySeedCollection &tracks) ;
    
    /// Register this collection of tracks coming from event id 'id'
    void registerTrajEvent(int id, const std::vector<Trajectory> &tracks) ;
    
    /// Register this collection of tracks coming from event id 'id'
    void registerTrackCandEvent(int id, const std::vector<TrackCandidate> &tracks) ;

    /// Register this collection of clusters coming from event id 'id'
    void registerClusterEvent(int id, const SiPixelClusterCollection &allPixels, const SiStripClusterCollection &allStrips) ;

    /// Try to associate track 'tk' with the registered track collections, fill in the vector of TrackAssociation.
    void associateToTracks(const reco::Track &tk, std::vector<TrackAssociation> &out) ;

    /// Try to associate track 'tk' with the registered track collections, fill in the vector of TrackAssociation.
    void associateToTracks(const TrajectorySeed &tk, std::vector<TrackAssociation> &out) ;

    /// Try to associate track 'tk' with the registered track collections, fill in the vector of TrackAssociation.
    void associateToTracks(const Trajectory &tk, std::vector<TrackAssociation> &out) ;

    /// Try to associate track 'tk' with the registered track collections, fill in the vector of TrackAssociation.
    void associateToTracks(const TrackCandidate &tk, std::vector<TrackAssociation> &out) ;

    /// Try to associate track 'tk' with the registered track collections, fill in the vector of TrackAssociation.
    void associateToTracks(const std::vector<const TrackingRecHit *> &tk, std::vector<TrackAssociation> &out) ;

    /// Try to associate track 'tk' with the registered track collections, fill in the vector of TrackAssociation.
    void associateToTracks(const TrackingRecHit &hit, std::vector<TrackAssociation> &out) ;


    /// Try to associate track 'tk' with the registered track collections, fill in the vector of TrackAssociation.
    void associateToSeed(const reco::Track &tk, std::vector<SeedAssociation> &out) ;

    /// Try to associate track 'tk' with the registered track collections, fill in the vector of TrackAssociation.
    void associateToSeed(const TrajectorySeed &tk, std::vector<SeedAssociation> &out) ;

    /// Try to associate track 'tk' with the registered track collections, fill in the vector of TrajAssociation.
    void associateToTrajs(const reco::Track &tk, std::vector<TrajAssociation> &out) ;

    /// Try to associate track 'tk' with the registered track collections, fill in the vector of TrajAssociation.
    void associateToTrajs(const TrajectorySeed &tk, std::vector<TrajAssociation> &out) ;

    /// Try to associate track 'tk' with the registered track collections, fill in the vector of TrajAssociation.
    void associateToTrajs(const Trajectory &tk, std::vector<TrajAssociation> &out) ;

    /// Try to associate track 'tk' with the registered track collections, fill in the vector of TrajAssociation.
    void associateToTrajs(const std::vector<const TrackingRecHit *> &tk, std::vector<TrajAssociation> &out) ;

    /// Try to associate track 'tk' with the registered track collections, fill in the vector of TrackCandAssociation.
    void associateToTrackCands(const reco::Track &tk, std::vector<TrackCandAssociation> &out) ;

    /// Try to associate track 'tk' with the registered track collections, fill in the vector of TrackCandAssociation.
    void associateToTrackCands(const TrajectorySeed &tk, std::vector<TrackCandAssociation> &out) ;

    /// Try to associate track 'tk' with the registered track collections, fill in the vector of TrackCandAssociation.
    void associateToTrackCands(const Trajectory &tk, std::vector<TrackCandAssociation> &out) ;

    /// Try to associate track 'tk' with the registered track collections, fill in the vector of TrackCandAssociation.
    void associateToTrackCands(const std::vector<const TrackingRecHit *> &tk, std::vector<TrackCandAssociation> &out) ;

    /// Try to associate track 'tk' with the registered cluster collections, fill in the vector of TrackClusterAssociation.
    void associateToClusters(const reco::Track &tk, std::vector<TrackClusterAssociation> &out) ;

    /// Try to associate track 'tk' with the registered cluster collections, fill in the vector of TrackClusterAssociation.
    void associateToClusters(const TrajectorySeed &tk, std::vector<TrackClusterAssociation> &out) ;

    /// Get clusters associated to a cluster
    void getAssociatedClusters(const TrackingRecHit &hit, std::vector<std::pair<int, const SiStripCluster *>> &out) ;
    void getAssociatedClusters(uint32_t detid, const SiStripCluster &hit, std::vector<std::pair<int, const SiStripCluster *>> &out);
    void getAssociatedClusters(const TrackingRecHit &hit, std::vector<std::pair<int, const SiPixelCluster *>> &out) ;
    void getAssociatedClusters(uint32_t detid, const SiPixelCluster &hit, std::vector<std::pair<int, const SiPixelCluster *>> &out);

private:
    template<typename T>
    struct RecRecord {
        int   eventId;
        const TrackingRecHit *hit;
        const T              *track;
        RecRecord(int eventID, const TrackingRecHit *h, const T *t) : eventId(eventID), hit(h), track(t) {}
    };
    typedef RecRecord<reco::Track> HitRecord;
    typedef RecRecord<TrajectorySeed> SeedRecord;
    typedef RecRecord<Trajectory> TrajRecord;
    typedef RecRecord<TrackCandidate> TrackCandRecord;
    boost::unordered_map<uint32_t, std::vector<HitRecord> > allHits_;
    boost::unordered_map<uint32_t, std::vector<SeedRecord> > allSeedHits_;
    boost::unordered_map<uint32_t, std::vector<TrajRecord> > allTrajHits_;
    boost::unordered_map<uint32_t, std::vector<TrackCandRecord> > allTrackCandHits_;

    struct ClusterRecord {
        int   eventId;
        const SiStripCluster * stripCluster;
        const SiPixelCluster * pixelCluster;
        ClusterRecord(int id, const SiStripCluster &clust) : eventId(id), stripCluster(&clust), pixelCluster(0)      {}
        ClusterRecord(int id, const SiPixelCluster &clust) : eventId(id), stripCluster(0),      pixelCluster(&clust) {}
    };
    boost::unordered_map<uint32_t, std::vector<ClusterRecord> > allClusters_;

    // --- register helpers ---
    template<typename T>
    void registerClusters(int eventId, const edmNew::DetSetVector<T> &clusters) ;

    // --- match helpers ---
    static bool matches(const SiStripCluster &c1, const SiStripCluster &c2) ;
    static bool matches(const SiPixelCluster &c1, const SiPixelCluster &c2) ;

    template<typename T, typename Record>
    static bool matchesRecordSpecific(const TrackingRecHit &hit, const Record &record) ;

    template<typename Record>
    static bool matchesRecordTracker1D(const TrackingRecHit &hit, const Record &record) ;

    template<typename Record>
    static bool matchesRecord(const TrackingRecHit &hit, const Record &record) ;

    static bool matches(const TrackingRecHit &hit, const ClusterRecord &record) ;

    template<typename T>
    void associateToRec(const reco::Track &tk, std::vector<RecAssociation<T> > &out, boost::unordered_map<uint32_t, std::vector<RecRecord<T> > > &record) ;

    template<typename T>
    void associateToRec(const TrajectorySeed &tk, std::vector<RecAssociation<T> > &out, boost::unordered_map<uint32_t, std::vector<RecRecord<T> > > &record) ;

    template<typename T>
    void associateToRec(const TrackCandidate &tk, std::vector<RecAssociation<T> > &out, boost::unordered_map<uint32_t, std::vector<RecRecord<T> > > &record) ;

    template<typename T>
    void associateToRec(const Trajectory &tk, std::vector<RecAssociation<T> > &out, boost::unordered_map<uint32_t, std::vector<RecRecord<T> > > &record) ;

    template<typename T>
    void associateToRec(const std::vector<const TrackingRecHit *> &tk, std::vector<RecAssociation<T> > &out, boost::unordered_map<uint32_t, std::vector<RecRecord<T> > > &record) ;

    template<typename T>
    void associateToRec(const TrackingRecHit &tk, std::vector<RecAssociation<T> > &out, boost::unordered_map<uint32_t, std::vector<RecRecord<T> > > &record) ;

    template<typename T>
    void associateHitToRec(const TrackingRecHit &hit, std::vector<RecAssociation<T> > &out, boost::unordered_map<uint32_t, std::vector<RecRecord<T> > > &record) ;

    void associateHitToClusters(const TrackingRecHit &hit, std::vector<TrackClusterAssociation> &out) ;
};

#endif
