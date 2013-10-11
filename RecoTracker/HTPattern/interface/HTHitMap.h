#ifndef RecoTracker_HTPattern_HTHitMap 
#define RecoTracker_HTPattern_HTHitMap

#include "HTHits.h"

class HTCell {   
    public:
        HTCell() : nlayers_(0), layermask_(0), nhits_(0) {}
        void addHit(unsigned int hit, unsigned int layermask) {
            hits_.push_back(hit); nhits_++;
            if ((layermask & layermask_) == 0) {
                layermask_ += layermask;
                nlayers_++;
            }
        }
        void clear() { if (nhits_) { nlayers_ = 0; layermask_ = 0; nhits_ = 0; hits_.clear(); } }
        unsigned int nhits() const { return nhits_; }
        unsigned int nlayers() const { return nlayers_; }
    protected:
        unsigned int nlayers_, layermask_, nhits_;
        std::vector<unsigned int> hits_; 
};

class HTHitMap {
    public:
        HTHitMap(unsigned int neta, unsigned int nphi) :
            neta_(neta), nphi_(nphi), cells_(neta*nphi) {}

        HTCell & get(unsigned int ieta, unsigned int iphi) { 
           return cells_[ieta*neta_+iphi];
        }
        const HTCell & get(unsigned int ieta, unsigned int iphi) const { 
           return cells_[ieta*neta_+iphi];
        }
        unsigned int etabins() const { return neta_; }
        unsigned int phibins() const { return nphi_; }
    protected:
        unsigned int neta_, nphi_;
        std::vector<HTCell> cells_; 
};
#endif
