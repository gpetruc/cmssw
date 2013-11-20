class Trajectory;
class TrackCandidate;
class MagneticField;
class TrackingGeometry;
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "RecoTracker/HTPattern/interface/HTHits.h"
#include "RecoTracker/HTPattern/interface/HTHitMap.h"
#include "DataFormats/TrackReco/interface/Track.h"

namespace HTDebugger {
    void dumpTracks(const std::string &name, const std::vector<reco::Track> & tracks, double pT=0) ;
    void dumpHTHits3D(const std::string &name, const HTHits3D &hits, double zref) ;
    void dumpHTHitsSpher(const std::string &name, const HTHitsSpher &hits) ;
    void dumpHTHitMap(const std::string &name, const HTHitMap &map);
    void dumpHTClusters(const std::string &name, const HTHitMap &map, const HTHitsSpher &hits, unsigned int minlayers, unsigned int minmorelayers=0) ;

    void registerTracks(const std::vector<reco::Track> & tracks) ;
    void registerClustersAsSeeds(int step, const TrajectorySeedCollection &seed);
    void registerTrackCandidates(const std::vector<TrackCandidate> &tcs);
    void printAssociatedTracks(const TrajectorySeed &cluster);
    void printAssociatedTracks(const Trajectory &traj, int minHits=1);
    void printAssociatedTracks(const TrackingRecHit *hit1, const TrackingRecHit *hit2);
    bool thirdHitOnSameTrack(const TrackingRecHit *hit1, const TrackingRecHit *hit2, const TrackingRecHit *hit3);
    void printAssociatedSeeds(const std::vector<reco::Track> & tracks) ;
    void printBackAssociation(const std::vector<reco::Track> & tracks, const std::vector<TrackCandidate> &tcs, const TrackingGeometry &g, const MagneticField &mf, float ptMin) ;
    void beginLoggingCluster(const TrajectorySeed &cluster, unsigned int seedhits, unsigned int hits, unsigned int morehits, double alpha) ;
    void logSeedingPair(const TrackingRecHit *hit1, int layer1, int class1, const TrackingRecHit *hit2, int layer2, int class2, float dr, float dphi, float deta) ;
    void logThirdHit(const TrackingRecHit *hit3, int layer3, int class3, float d2c, float dphic, float detac, int pass, float thisalpha, bool finalsel);
    void debugHelixParameters(const std::vector<const TrackingRecHit *> & hits, float alpha, float beta, float phi0, float eta0, float bfieldAtOrigin, GlobalPoint vertex) ;
    const reco::Track * matchCandidate(const std::vector<const TrackingRecHit *> & hits) ;
    int  isMatchedToTrack(const TrackingRecHit *hit, const reco::Track *tk) ;

};
