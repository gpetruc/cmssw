#ifndef RecoTracker_HTPattern_HTHitMap 
#define RecoTracker_HTPattern_HTHitMap

#include "HTHits.h"
#include <array>
#include <cstdio>
#include <cassert>

class HTCell {   
    public:
        HTCell() : nlayers_(0), layermask_(0), nhits_(0), cluster_(-1) {}
        unsigned int addHit(unsigned int hit, unsigned int layermask) {
            hits_.push_back(hit); nhits_++;
            //printf("adding hit %d with layermask %x to cell @%p with layermask %x, layers %d: ", hit, layermask, (void*)(this), layermask_, nlayers_);
            if ((layermask & layermask_) == 0) {
                layermask_ += layermask;
                nlayers_++;
                //printf(" ---> layermask %x, layers %d\n", layermask_, nlayers_);
                return nlayers_;
            }
            //printf(" ---> layermask %x, layers %d\n", layermask_, nlayers_);
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
            ieta_(ieta), iphi_(iphi), 
            nseedlayers_(cell.nlayers()), seedlayermask_(cell.layermask()),
            nlayers_(cell.nlayers()), layermask_(cell.layermask()), 
            nmorelayers_(cell.nlayers()), morelayermask_(cell.layermask())
            { }
            //{ printf("created cluster from cell $%p: layers %d, layermask %d\n", (void*)(&cell), nseedlayers_, seedlayermask_); }
        unsigned int ieta() const { return ieta_; }
        unsigned int iphi() const { return iphi_; }
        unsigned int nlayers() const { return nlayers_; }
        unsigned int nseedlayers() const { return nseedlayers_; }
        unsigned int nmorelayers() const { return nmorelayers_; }
        void updateSeed(const HTCell &cell) {
            nseedlayers_ = cell.nlayers(); seedlayermask_ = cell.layermask();
            nlayers_ = cell.nlayers();     layermask_ = cell.layermask(); 
            nmorelayers_ = cell.nlayers(); morelayermask_ = cell.layermask();
        }
        void addCell(const HTCell &cell) {
            addCell(cell,layermask_,nlayers_);
            morelayermask_ = layermask_;
            nmorelayers_   = nlayers_;
        }
        void addMoreCell(const HTCell &cell) {
            addCell(cell,morelayermask_,nmorelayers_);
        }
        void addCell(const HTCell &cell, unsigned int &layermask, unsigned int &nlayers) {
            unsigned int addmask = (cell.layermask() & (~layermask));
            //printf("Adding a cell: starting mask %x (%d), cell %x, addmask %x", layermask, nlayers, cell.layermask(), addmask);
            for (; addmask > 0; addmask >>= 1) {
                nlayers += (addmask & 1);
            }
            layermask |= cell.layermask();
            //printf(", ending mask %x (%d)\n", layermask, nlayers);
        }
    protected:
        unsigned int ieta_, iphi_;
        unsigned int nseedlayers_, seedlayermask_;
        unsigned int nlayers_, layermask_;
        unsigned int nmorelayers_, morelayermask_;
};

class HTHitMap {
    public:
        HTHitMap(unsigned int neta, unsigned int nphi) :
            neta_(neta), nphi_(nphi), cells_(neta*nphi), nclusters_(0) {}

        //HTCell & get(unsigned int ieta, unsigned int iphi) { 
        //   return cells_[ieta*nphi_+iphi];
        //}
        const HTCell & get(unsigned int ieta, unsigned int iphi) const { 
           assert(ieta < neta_);
           assert(iphi < nphi_);
           return cells_[ieta*nphi_+iphi];
        }
        void clear(unsigned int ieta, unsigned int iphi) {
           assert(ieta < neta_);
           assert(iphi < nphi_);
           return cells_[ieta*nphi_+iphi].clear();
        }
        void clearClusters() { clusters_.clear(); nclusters_ = 0; }
        void addHit(unsigned int ieta, unsigned int iphi, unsigned int hit, unsigned int layermask, unsigned int layerSeedCut) {
            HTCell &cell = cells_[ieta*nphi_+iphi];
            if (cell.addHit(hit,layermask) >= layerSeedCut && cell.icluster() == -1) {
                cell.setCluster(nclusters_);
                clusters_.push_back(HTCluster(ieta,iphi,cell));
                nclusters_++;
            }   
        }
        //void getNeighbours(unsigned int ieta, unsigned int iphi, std::array<HTCell *, 8> & cells, unsigned int &ncells);
        void getNeighbours(unsigned int ieta, unsigned int iphi, std::array<const HTCell *, 8> & cells, unsigned int &ncells) const ;
        unsigned int etabins() const { return neta_; }
        unsigned int phibins() const { return nphi_; }
        const std::vector<HTCluster> & clusters() const { return clusters_; }
        std::vector<HTCluster> & clusters() { return clusters_; }
        void clusterize() ;
    protected:
        unsigned int neta_, nphi_;
        std::vector<HTCell> cells_; 
        unsigned int nclusters_;
        std::vector<HTCluster> clusters_;
};
#endif
