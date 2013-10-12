#ifndef RecoTracker_HTPattern_HTHits
#define RecoTracker_HTPattern_HTHits

#include <vector>
#include <cmath>
class TrackingRecHit;
class HTHitsSpher;

class HTHits3D {
    public:
        HTHits3D() : size_(0) {}
        void reserve(int size) {
            hit_.reserve(size);
            rho_.reserve(size);
            phi_.reserve(size);
            z_.reserve(size);
            layermask_.reserve(size);
        }
        void push_back(const TrackingRecHit &hit, double rho, double phi, double z, int layermask) {
            size_++;
            hit_.push_back(&hit);
            rho_.push_back(rho);
            phi_.push_back(phi < 0 ? phi + float(2*M_PI) : phi);
            z_.push_back(z);
            layermask_.push_back(layermask);
        }            
        float rho(int i) const { return rho_[i]; }
        float z(int i)   const { return z_[i]; }
        float phi(int i) const { return phi_[i]; }
        const TrackingRecHit * hit(int i) const { return hit_[i]; }
        int   layermask(int i) const { return layermask_[i]; }
        unsigned int size() const { return size_; }
        friend class HTHitsSpher;
    protected:
        unsigned int size_;
        std::vector<float> rho_, phi_, z_;
        std::vector<const TrackingRecHit *> hit_;
        std::vector<int> layermask_;
};

class HTHitsSpher {
    public:
        HTHitsSpher(const HTHits3D &hits3d, float z0, unsigned int etabins) ;
        float eta(int i) const { return eta_[i]; }
        int   ieta(int i) const { return ieta_[i]; }
        float rho(int i) const { return (*rho_)[i]; }
        float phi(int i) const { return (*phi_)[i]; }
        float z(int i) const { return (*z_)[i]; }
        int   iphi(int i) const { return iphi_[i]; }
        const TrackingRecHit * hit(int i) const { return (*hit_)[i]; }
        int   layermask(int i) const { return (*layermask_)[i]; }
        void  filliphi(float alpha, unsigned int phibins) ;
        unsigned int size() const { return size_; }
        unsigned int etabins() const { return etabins_; }
        unsigned int phibins() const { return phibins_; }
        float alpha() const { return alpha_; }
        float z0() const { return z0_; }
    protected:
        unsigned int size_, etabins_, phibins_;
        float etascale_, phiscale_;
        float z0_, alpha_;
        const std::vector<float> *rho_, *phi_, *z_;
        std::vector<float> eta_;
        std::vector<int>   ieta_;
        std::vector<int>   iphi_;
        const std::vector<const TrackingRecHit *> *hit_;
        const std::vector<int> *layermask_;
};

#endif
