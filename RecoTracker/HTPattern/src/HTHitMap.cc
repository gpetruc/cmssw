#include "RecoTracker/HTPattern/interface/HTHitMap.h"

void HTHitMap::clusterize(unsigned int ieta, unsigned int iphi, HTCell &cell) {
    cell.setCluster(nclusters_);
    clusters_.push_back(HTCluster(ieta,iphi,cell));
    if (ieta > 0)       clusters_.back().addCell(get(ieta-1,iphi));
    if (ieta < neta_-1) clusters_.back().addCell(get(ieta+1,iphi));
    clusters_.back().addCell(get(ieta, (iphi > 0 ? iphi : nphi_) - 1 ));
    clusters_.back().addCell(get(ieta, (iphi + 1 < nphi_ ? iphi+1 : 0)));
    nclusters_++;
}
