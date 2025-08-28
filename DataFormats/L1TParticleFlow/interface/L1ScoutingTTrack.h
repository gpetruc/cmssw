#ifndef DataFormats_L1TParticleFlow_TTrack_h
#define DataFormats_L1TParticleFlow_TTrack_h

#include <algorithm>
#include <array>
#include <bitset>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace l1Scouting {
  class TTrack {
  public:
    TTrack() {}
    TTrack(unsigned int valid,
      unsigned int rInv,
      unsigned int phi0, 
      unsigned int tanl,
      unsigned int z0,
      unsigned int d0,
      unsigned int chi2RPhi,
      unsigned int chi2RZ,
      unsigned int bendChi2,
      unsigned int hitPattern,
      unsigned int mvaQuality,
      unsigned int mvaOther) : 
        valid_(valid),
        rInv_(rInv),
        phi0_(phi0), 
        tanl_(tanl),
        z0_(z0),
        d0_(d0),
        chi2RPhi_(chi2RPhi),
        chi2RZ_(chi2RZ),
        bendChi2_(bendChi2),
        hitPattern_(hitPattern),
        mvaQuality_(mvaQuality),
        mvaOther_(mvaOther) {}

  unsigned int valid() const { return valid_; }
  unsigned int rInv() const { return rInv_; }
  unsigned int phi0() const { return phi0_; }
  unsigned int tanl() const { return tanl_; }
  unsigned int z0() const { return z0_; }
  unsigned int d0() const { return d0_; }
  unsigned int chi2RPhi() const { return chi2RPhi_; }
  unsigned int chi2RZ() const { return chi2RZ_; }
  unsigned int bendChi2() const { return bendChi2_; }
  unsigned int hitPattern() const { return hitPattern_; }
  unsigned int mvaQuality() const { return mvaQuality_; }
  unsigned int mvaOther() const { return mvaOther_; }

  private:
    unsigned int valid_, rInv_, phi0_, tanl_, z0_, d0_;
    unsigned int chi2RPhi_, chi2RZ_, bendChi2_, hitPattern_, mvaQuality_, mvaOther_;
  };

  struct TTrackSOA {
    std::vector<uint16_t> bx;
    std::vector<uint32_t> offsets;
    std::vector<unsigned int> valid, rInv, phi0, tanl, z0, d0;
    std::vector<unsigned int> chi2RPhi, chi2RZ, bendChi2, hitPattern, mvaQuality, mvaOther; //quality?
    TTrackSOA() : bx(), offsets(), valid(), rInv(), phi0(), tanl(), z0(), d0(), chi2RPhi(), chi2RZ(), bendChi2(), hitPattern() {}
    TTrackSOA(const TTrackSOA& other) = default;
    TTrackSOA(TTrackSOA&& other) = default;
    TTrackSOA& operator=(const TTrackSOA& other) = default;
    TTrackSOA& operator=(TTrackSOA&& other) = default;
  };
};
#endif
