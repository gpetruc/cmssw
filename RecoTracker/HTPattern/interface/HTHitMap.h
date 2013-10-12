#ifndef RecoTracker_HTPattern_HTHitMap 
#define RecoTracker_HTPattern_HTHitMap

#include "HTHits.h"

class HTCell {   
    public:
        HTCell() : nlayers_(0), layermask_(0), nhits_(0), cluster_(-1) {}
        unsigned int addHit(unsigned int hit, unsigned int layermask) {
            hits_.push_back(hit); nhits_++;
            if ((layermask & layermask_) == 0) {
                layermask_ += layermask;
                nlayers_++;
                return nlayers_;
            }
            return 0;
        }
        void clear() { if (nhits_) { nlayers_ = 0; layermask_ = 0; nhits_ = 0; cluster_ = -1; hits_.clear(); } }
        unsigned int nhits() const { return nhits_; }
        unsigned int nlayers() const { return nlayers_; }
        unsigned int layermask() const { return layermask_; }
        int icluster() const { return cluster_; }
        void setCluster(int cluster) { cluster_ = cluster; }
        const std::vector<unsigned int> & hits() const  { return hits_; }
    protected:
        unsigned int nlayers_, layermask_, nhits_;
        int cluster_;
        std::vector<unsigned int> hits_; 
};

class HTCluster {
    public:
        HTCluster(unsigned int ieta, unsigned int iphi, const HTCell &cell) :
            ieta_(ieta), iphi_(iphi), nlayers_(cell.nlayers()), layermask_(cell.layermask()) {}
        unsigned int ieta() const { return ieta_; }
        unsigned int iphi() const { return iphi_; }
        unsigned int nlayers() const { return nlayers_; }
        void addCell(const HTCell &cell) {
            unsigned int addmask = (cell.layermask() & (~layermask_));
            for (; addmask > 0; addmask >>= 1) {
                nlayers_ += (addmask & 1);
            }
        }
    protected:
        unsigned int ieta_, iphi_;
        unsigned int nlayers_, layermask_;
};

class HTHitMap {
    public:
        HTHitMap(unsigned int neta, unsigned int nphi) :
            neta_(neta), nphi_(nphi), cells_(neta*nphi), nclusters_(0) {}

        //HTCell & get(unsigned int ieta, unsigned int iphi) { 
        //   return cells_[ieta*nphi_+iphi];
        //}
        const HTCell & get(unsigned int ieta, unsigned int iphi) const { 
           return cells_[ieta*nphi_+iphi];
        }
        void clear(unsigned int ieta, unsigned int iphi) {
           return cells_[ieta*nphi_+iphi].clear();
        }
        void clearClusters() { clusters_.clear(); nclusters_ = 0; }
        void addHit(unsigned int ieta, unsigned int iphi, unsigned int hit, unsigned int layermask, unsigned int layerSeedCut) {
           HTCell &cell = cells_[ieta*nphi_+iphi];
           if (cell.addHit(hit,layermask) > layerSeedCut && cell.icluster() == -1) {
            clusterize(ieta,iphi,cell);
           }   
        }
        unsigned int etabins() const { return neta_; }
        unsigned int phibins() const { return nphi_; }
        const std::vector<HTCluster> & clusters() const { return clusters_; }
    protected:
        void clusterize(unsigned int ieta, unsigned int iphi, HTCell &cell);
        unsigned int neta_, nphi_;
        std::vector<HTCell> cells_; 
        unsigned int nclusters_;
        std::vector<HTCluster> clusters_;
};
#endif
