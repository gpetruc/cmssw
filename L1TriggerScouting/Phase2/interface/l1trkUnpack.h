#ifndef L1TriggerScouting_Phase2_l1trkUnpack_h
#define L1TriggerScouting_Phase2_l1trkUnpack_h
#include <cstdint>
#include <cmath>

#include "CLHEP/Units/GlobalPhysicalConstants.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack_TrackWord.h"

// 0 valid - double 1
// 15-1 Rinv - double 15
// 27-16 Phi - double 12
// 31-28 Chi2RPhi - double 4
// 47-32 Tanl - double 16
// 59-48 Z0 - double 12
// 63-60 Chi2RZ - double 4
// 76-64 D0 - double 13
// 79-77 BendChi2 - double 3
// 86-80 HitPattern - unsigned int 7
// 89-87 MVAQuality - double 3
// 95-90 MVAOther - double 6

namespace l1trkUnpack {
  // constant is 0.299792458; who knew c_light was in mm/ns?
  static constexpr float MagConstant = CLHEP::c_light / 1.0E3;
  static constexpr float BField = 3.81120228767395; // in T

  inline void read(const uint64_t datalow,
                   const uint32_t datahigh,
                   uint32_t &rInv,
                   uint32_t &phi0,
                   uint32_t &chi2RPhi,
                   uint32_t &tanl,
                   uint32_t &z0,
                   uint32_t &chi2RZ,
                   uint32_t &d0,
                   uint32_t &bendChi2,
                   uint32_t &hitPattern,
                   uint32_t &mvaQuality,
                   uint32_t &mvaOther) {
    // rInv = ((datalow >> 15) & 1) ? ((datalow >> 1) | (-0x4000)) : ((datalow >> 1) & (0x7FFF));  // 15 bits
    // phi0 = ((datalow >> 27) & 1) ? ((datalow >> 16) | (-0x800)) : ((datalow >> 16) & (0xFFF));  // 12 bits
    // chi2RPhi = ((datalow >> 31) & 1) ? ((datalow >> 28) | (-0x8)) : ((datalow >> 28) & (0x8));  // 4 bits
    // tanl = ((datalow >> 47) & 1) ? ((datalow >> 32) | (-0x8000)) : ((datalow >> 32) & (0xFFFF));  // 16 bits
    // z0 = ((datalow >> 59) & 1) ? ((datalow >> 48) | (-0x800)) : ((datalow >> 48) & (0xFFF));  // 12 bits
    // chi2RZ = ((datalow >> 63) & 1) ? ((datalow >> 60) | (-0x8)) : ((datalow >> 60) & (0x8));  // 4 bits

    // d0 = ((datahigh >> 12) & 1) ? ((datahigh >> 0) | (-0x1000)) : ((datahigh >> 0) & (0x1FFF));  // 13 bits
    // bendChi2 = ((datahigh >> 15) & 1) ? ((datahigh >> 13) | (-0x4)) : ((datahigh >> 13) & (0x4));  // 3 bits
    // hitPattern = ((datahigh >> 22) & 1) ? ((datahigh >> 16) | (-0x7F)) : ((datahigh >> 16) & (0x7F));  // 7 bits
    // mvaQuality = ((datahigh >> 25) & 1) ? ((datahigh >> 23) | (-0x7)) : ((datahigh >> 23) & (0x7));  // 3 bits
    // mvaOther = ((datahigh >> 31) & 1) ? ((datahigh >> 26) | (-0x7)) : ((datahigh >> 26) & (0x7));  // 6 bits
    rInv = (datalow >> 1) & (0x7FFF);  // 15 bits
    phi0 = (datalow >> 16) & (0xFFF);  // 12 bits
    chi2RPhi = (datalow >> 28) & (0x8);  // 4 bits
    tanl = (datalow >> 32) & (0xFFFF);  // 16 bits
    z0 = (datalow >> 48) & (0xFFF);  // 12 bits
    chi2RZ = (datalow >> 60) & (0x8);  // 4 bits

    d0 = (datahigh >> 0) & (0x1FFF);  // 13 bits
    bendChi2 = (datahigh >> 13) & (0x4);  // 3 bits
    hitPattern = (datahigh >> 16) & (0x7F);  // 7 bits
    mvaQuality = (datahigh >> 23) & (0x7);  // 3 bits
    mvaOther = (datahigh >> 26) & (0x7);  // 6 bits
  }

  inline unsigned int countSetBits(unsigned int n) {
    unsigned int count = 0;
    while (n) {
      n &= (n - 1);
      count++;
    }
    return count;
  }

  inline float undigitizeSignedValue(unsigned int twosValue, unsigned int nBits, double lsb, double offset = 0.5) {
    // Check that none of the bits above the nBits-1 bit, in a range of [0, nBits-1], are set.
    // This makes sure that it isn't possible for the value represented by `twosValue` to be
    //  any bigger than ((1 << nBits) - 1).
    assert((twosValue >> nBits) == 0);

    // Convert from twos complement to C++ signed integer (normal digitized value)
    int digitizedValue = twosValue;
    if (twosValue & (1 << (nBits - 1))) {  // check if the twosValue is negative
      digitizedValue -= (1 << nBits);
    }

    // Convert to floating point value
    return (float(digitizedValue) + offset) * lsb;
  }

  inline float getRinv(uint32_t rInvInt) {
    return undigitizeSignedValue(rInvInt, TTTrack_TrackWord::TrackBitWidths::kRinvSize, TTTrack_TrackWord::stepRinv);
  }

  inline float getPhi0(uint32_t phi0Int) {
    return undigitizeSignedValue(phi0Int, TTTrack_TrackWord::TrackBitWidths::kPhiSize, TTTrack_TrackWord::stepPhi0);
  }

  inline float getTanl(uint32_t tanlInt) {
    return undigitizeSignedValue(tanlInt, TTTrack_TrackWord::TrackBitWidths::kTanlSize, TTTrack_TrackWord::stepTanL);
  }

  inline float getZ0(uint32_t z0Int) {
    return undigitizeSignedValue(z0Int, TTTrack_TrackWord::TrackBitWidths::kZ0Size, TTTrack_TrackWord::stepZ0);
  }

  inline float getD0(uint32_t d0Int) {
    return undigitizeSignedValue(d0Int, TTTrack_TrackWord::TrackBitWidths::kD0Size, TTTrack_TrackWord::stepD0);
  }

  inline float getChi2RPhi(uint32_t chi2RPhiInt) {
    return TTTrack_TrackWord::chi2RPhiBins[chi2RPhiInt];
  }

  inline float getChi2RZ(uint32_t chi2RZInt) {
    return TTTrack_TrackWord::chi2RZBins[chi2RZInt];
  }

  inline float getBendChi2(uint32_t bendChi2Int) {
    return TTTrack_TrackWord::bendChi2Bins[bendChi2Int];
  }

  inline unsigned int getNStubs(uint32_t hitPattern) {
    return countSetBits(hitPattern);
  }

  inline float getMVAQuality(uint32_t mvaQuality) {
    return TTTrack_TrackWord::tqMVABins[mvaQuality];
  }

  inline float getPt(uint32_t rInvInt) {
    return std::abs(MagConstant / getRinv(rInvInt) * MagConstant / 100.0);  // Rinv is in cm-1
  }

  inline GlobalVector getMomentum(float pt, float phi0, float tanl) {
    return GlobalVector(GlobalVector::Cylindrical(pt, phi0, pt * tanl));
  }

  inline GlobalPoint getPOCA(float d0, float phi0, float z0) {
    return GlobalPoint(d0 * sin(phi0), -d0 * cos(phi0), z0);
  }
}

#endif
