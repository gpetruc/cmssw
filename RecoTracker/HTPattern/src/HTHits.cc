#include "RecoTracker/HTPattern/interface/HTHits.h"
#include <cmath>
#include <cassert>
#include <cstdio>
#include "vdt/vdtMath.h"

namespace {
    void vector_iphi(const uint32_t size, float const * __restrict__ phi0, float const * __restrict__ rho, float alpha, float phiscale, unsigned int phimask, unsigned int * __restrict__ iphi) {
        for (unsigned int i = 0; i < size; ++i) {
            iphi[i] = (phi0[i] - alpha*rho[i]) * phiscale;
        //}
        //for (unsigned int i = 0; i < size; ++i) {
            iphi[i] = iphi[i] & phimask;
        }
    }
    void vector_eta(const uint32_t size,float const * __restrict__ z, float const * __restrict__ rho, float z0, float etascale, unsigned int etamask, float * __restrict__ eta, unsigned int * __restrict__ ieta) {
        using vdt::fast_invf;
        for (unsigned int i = 0; i < size; ++i) {
            eta[i]  = HTHitsSpher::pseudoEta((z[i]-z0)*fast_invf(rho[i]));
        //} 
        //for (unsigned int i = 0; i < size; ++i) {
            ieta[i] = (eta[i]+3.0f)*etascale;
        //} 
        //for (unsigned int i = 0; i < size; ++i) {
            ieta[i] = ieta[i] & etamask;
        } 
    }
}

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
    using vdt::details::fpfloor;
    eta_.resize(hits3d.size_);
    ieta_.resize(hits3d.size_);
#if 1
    unsigned int etamask = etabins - 1;
    ::vector_eta(size_, &(*z_)[0], &(*rho_)[0], z0, etascale_, etamask, &eta_[0], &ieta_[0]);
#else
    for (unsigned int i = 0; i < size_; ++i) {
        //eta_[i] = std::asinh((hits3d.z(i)-z0)/hits3d.rho(i));
        eta_[i] = pseudoEta((hits3d.z(i)-z0)/hits3d.rho(i));
        ieta_[i] = std::min<unsigned int>(fpfloor( std::max(0.f, (3.0f + eta_[i]) * etascale_ ) ), etabins-1); 
    }
#endif
}

void
HTHitsSpher::filliphi(float alpha, unsigned int phibins)
{
    using vdt::details::fpfloor;
    alpha_ = alpha;
    phibins_ = phibins;
    phiscale_ = phibins/(2*M_PI);
    iphi_.resize(size_);
#if 1
    unsigned int phimask_ = phibins-1;
    ::vector_iphi(size_, &(*phi_)[0], &(*rho_)[0], alpha, phiscale_, phimask_, &iphi_[0]);
#else
    unsigned int phimask_ = phibins-1;
    for (unsigned int i = 0; i < size_; ++i) {
        iphi_[i] = unsigned(fpfloor( phiNoLoop(i) * phiscale_ )) & phimask_;
        if (iphi_[i] >= phibins_) {
            printf("phi in floating point, scaled: %.8f\n", phiNoLoop(i) * phiscale_);
            printf("floor of phi in floating point, scaled: %.8f\n", fpfloor(phiNoLoop(i) * phiscale_));
            printf("unsigned int corresponding to phi in floating point, scaled: %u\n",  unsigned(fpfloor( phiNoLoop(i) * phiscale_ )));
            printf("masked unsigned int corresponding to phi in floating point, scaled: %u\n", iphi_[i]);
        }
        assert(iphi_[i] < phibins_);
    }
#endif
}


