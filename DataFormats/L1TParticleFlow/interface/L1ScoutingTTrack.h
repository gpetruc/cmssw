#ifndef DataFormats_L1TParticleFlow_TTrack_h
#define DataFormats_L1TParticleFlow_TTrack_h

#include <algorithm>
#include <array>
#include <bitset>
#include <cmath>
#include <limits>
#include <string>
#include <vector>


// template <typename T>
// TTTrack<T>::TTTrack(double aRinv,
//                     double aphi0,
//                     double aTanlambda,
//                     double az0,
//                     double ad0,
//                     double aChi2,
//                     double trkMVA1,
//                     double trkMVA2,
//                     double trkMVA3,
//                     double aHitPattern,
//                     double nPar,
//                     double aBfield) {
//   theStubRefs.clear();
//   
//   theMomentum_ = GlobalVector(GlobalVector::Cylindrical(thePT, aphi0, thePT * aTanlambda));
//   theRInv_ = aRinv;
//   thePOCA_ = GlobalPoint(ad0 * sin(aphi0), -ad0 * cos(aphi0), az0);
//   theD0_ = ad0;
//   theZ0_ = az0;
//   thePhi_ = aphi0;
//   theTanL_ = aTanlambda;
//   thePhiSector_ = 0;      // must be set externally
//   theEtaSector_ = 0;      // must be set externally
//   theTrackSeedType_ = 0;  // must be set externally
//   theChi2_ = aChi2;
//   theTrkMVA1_ = trkMVA1;
//   theTrkMVA2_ = trkMVA2;
//   theTrkMVA3_ = trkMVA3;
//   theStubPtConsistency_ = 0.0;  // must be set externally
//   theNumFitPars_ = nPar;
//   theHitPattern_ = aHitPattern;
//   theBField_ = aBfield;
//   theChi2_XY_ = -999.;
//   theChi2_Z_ = -999.;
// }

//TODO - check types, is it all supposed to be double?                                            

namespace l1Scouting {
  class TTrack {
  public:
    TTrack() {}
    TTrack(
      double rInv,
      double phi0, 
      double chi2RPhi,
      double tanl,
      double z0,
      double chi2RZ,
      double d0,
      double bendChi2,
      unsigned int hitPattern,
      double mvaQuality,
      double mvaOther) : 
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
        mvaOther_(mvaOther) {

          pt_ = 0; //TODO - calculate the pt

        }




  double rInv() const { return rInv_; }
  double phi0() const { return phi0_; }
  double tanl() const { return tanl_; }
  double z0() const { return z0_; }
  double d0() const { return d0_; }
  double chi2RPhi() const { return chi2RPhi_; }
  double chi2RZ() const { return chi2RZ_; }
  double bendChi2() const { return bendChi2_; }
  double hitPattern() const { return hitPattern_; }
  double mvaQuality() const { return mvaQuality_; }
  double mvaOther() const { return mvaOther_; }
  double pt() const { return pt_; }

  private:
    double rInv_, phi0_, tanl_, z0_, d0_, pt_;
    double chi2RPhi_, chi2RZ_, bendChi2_, hitPattern_, mvaQuality_, mvaOther_;
  };

  struct TTrackSOA {
    std::vector<uint16_t> bx;
    std::vector<uint32_t> offsets;
    std::vector<double> rInv, phi0, tanl, z0, d0, pt;
    std::vector<double> chi2RPhi, chi2RZ, bendChi2, hitPattern, mvaQuality, mvaOther;
    TTrackSOA() : bx(), offsets(), rInv(), phi0(), tanl(), z0(), d0(), chi2RPhi(), chi2RZ(), bendChi2(), hitPattern() {}
    TTrackSOA(const TTrackSOA& other) = default;
    TTrackSOA(TTrackSOA&& other) = default;
    TTrackSOA& operator=(const TTrackSOA& other) = default;
    TTrackSOA& operator=(TTrackSOA&& other) = default;
  };
};
#endif
