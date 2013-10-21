#include "RecoTracker/HTPattern/interface/HTHitMap.h"

#if 0
void HTHitMap::getNeighbours(unsigned int ieta, unsigned int iphi, std::array<HTCell *, 8> & cells, unsigned int &ncells) {
    ncells = 0;
    unsigned int iphim = (iphi > 0 ? iphi : nphi_) - 1;
    unsigned int iphip = (iphi + 1 < nphi_ ? iphi+1 : 0);
    cells[ncells++] = & get(iphim, ieta);
    cells[ncells++] = & get(iphip, ieta);
    if (ieta > 0) {
        cells[ncells++] = & get(iphim, ieta-1);
        cells[ncells++] = & get(iphi , ieta-1);
        cells[ncells++] = & get(iphip, ieta-1);
    }
    if (ieta < neta_-1) {
        cells[ncells++] = & get(iphim, ieta+1);
        cells[ncells++] = & get(iphi , ieta+1);
        cells[ncells++] = & get(iphip, ieta+1);
    }
}
#endif
void HTHitMap::getNeighbours(unsigned int ieta, unsigned int iphi, std::array<const HTCell *, 8> & cells, unsigned int &ncells) const {
    ncells = 0;
    unsigned int iphim = (iphi > 0 ? iphi : nphi_) - 1;
    unsigned int iphip = (iphi + 1 < nphi_ ? iphi+1 : 0);
    cells[ncells++] = & get(ieta, iphim);
    cells[ncells++] = & get(ieta, iphip);
    if (ieta > 0) {
        cells[ncells++] = & get(ieta-1, iphim);
        cells[ncells++] = & get(ieta-1, iphi );
        cells[ncells++] = & get(ieta-1, iphip);
    }
    if (ieta < neta_-1) {
        cells[ncells++] = & get(ieta+1, iphim);
        cells[ncells++] = & get(ieta+1, iphi );
        cells[ncells++] = & get(ieta+1, iphip);
    }
}



void HTHitMap::clusterize() {
    std::array<const HTCell *, 8> cells; unsigned int ncells;
    for (HTCluster & cluster : clusters_) {
        cluster.updateSeed(get(cluster.ieta(),cluster.iphi()));
        getNeighbours(cluster.ieta(),cluster.iphi(), cells, ncells);
        for (unsigned int i = 0; i < ncells; ++i) {
            cluster.addCell(*cells[i]);
        }
    }
}
