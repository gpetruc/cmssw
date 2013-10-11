#include "RecoTracker/HTPattern/interface/HTHits.h"
#include "RecoTracker/HTPattern/interface/HTHitMap.h"

namespace HTDebugger {
    void dumpHTHits3D(const std::string &name, const HTHits3D &hits, double zref) ;
    void dumpHTHitsSpher(const std::string &name, const HTHitsSpher &hits) ;
    void dumpHTHitMap(const std::string &name, const HTHitMap &map);
};
