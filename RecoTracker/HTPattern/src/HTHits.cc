#include "RecoTracker/HTPattern/interface/HTHits.h"
#include <cmath>

HTHitsSpher::HTHitsSpher(const HTHits3D &hits3d, float z0, unsigned int etabins) :
    size_(hits3d.size_),
    etabins_(etabins),
    etascale_(etabins/6.0),
    x0_(hits3d.x0()),
    y0_(hits3d.y0()),
    z0_(z0),
    rho_(&hits3d.rho_),
    phi_(&hits3d.phi_),
    z_(&hits3d.z_),
    hit_(&hits3d.hit_),
    layermask_(&hits3d.layermask_)
{
    eta_.resize(hits3d.size_);
    ieta_.resize(hits3d.size_);
    for (unsigned int i = 0; i < size_; ++i) {
        eta_[i] = std::asinh((hits3d.z(i)-z0)/hits3d.rho(i));
        ieta_[i] = std::min<unsigned int>(std::round( std::max(0.f, (3.0f + eta_[i]) * etascale_ ) ), etabins-1); 
    }
}

void
HTHitsSpher::filliphi(float alpha, unsigned int phibins)
{
    alpha_ = alpha;
    phibins_ = phibins;
    phiscale_ = phibins/(2*M_PI);
    iphi_.resize(size_);
    for (unsigned int i = 0; i < size_; ++i) {
        iphi_[i] = std::floor( phi(i) * phiscale_ );
    }
}


