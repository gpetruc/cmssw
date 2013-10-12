#include "RecoTracker/HTPattern/interface/HTHits.h"
#include "RecoTracker/HTPattern/interface/HTHitMap.h"
#include "DataFormats/TrackReco/interface/Track.h"

namespace HTDebugger {
    void dumpTracks(const std::string &name, const std::vector<reco::Track> & tracks) ;
    void dumpHTHits3D(const std::string &name, const HTHits3D &hits, double zref) ;
    void dumpHTHitsSpher(const std::string &name, const HTHitsSpher &hits) ;
    void dumpHTHitMap(const std::string &name, const HTHitMap &map);
    void dumpHTClusters(const std::string &name, const HTHitMap &map, const HTHitsSpher &hits) ;
};
