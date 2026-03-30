#ifndef L1TriggerScouting_Phase2_l1trkUnpack_h
#define L1TriggerScouting_Phase2_l1trkUnpack_h
#include <cstdint>
#include <cmath>

#include "DataFormats/L1TrackTrigger/interface/TTTrack_TrackWord.h"
#include "DataFormats/Math/interface/deltaPhi.h"

namespace l1trkUnpack {
  static constexpr float MagConstant = 0.299792458;
  static constexpr float BField = 3.81120228767395;  // in T

  inline void read(const uint64_t datalow,
                   const uint32_t datahigh,
                   uint32_t &valid,
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
                   uint32_t &phiSector) {
    phiSector = datalow & 0x3F;          // 6 bits   (0)
    mvaQuality = (datalow >> 6) & 0x7;   // 3 bits   (6)
    hitPattern = (datalow >> 9) & 0x7F;  // 7 bits   (9)
    bendChi2 = (datalow >> 16) & 0x7;    // 3 bits   (16)
    d0 = (datalow >> 19) & 0x1FFF;       // 13 bits  (19)
    chi2RZ = (datalow >> 32) & 0xF;      // 4 bits   (32)
    z0 = (datalow >> 36) & 0xFFF;        // 12 bits  (36)
    tanl = (datalow >> 48) & 0xFFFF;     // 16 bits  (48)

    chi2RPhi = (datahigh >> 0) & 0xF;  // 4 bits   (64)
    phi0 = (datahigh >> 4) & 0xFFF;    // 12 bits  (68)
    rInv = (datahigh >> 16) & 0x7FFF;  // 15 bits  (80)
    valid = (datahigh >> 31) & 0x1;    // 1 bit    (95)
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
    // any bigger than ((1 << nBits) - 1).
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

  inline float getChi2RPhi(uint32_t chi2RPhiInt) { return TTTrack_TrackWord::chi2RPhiBins[chi2RPhiInt]; }

  inline float getChi2RZ(uint32_t chi2RZInt) { return TTTrack_TrackWord::chi2RZBins[chi2RZInt]; }

  inline float getBendChi2(uint32_t bendChi2Int) { return TTTrack_TrackWord::bendChi2Bins[bendChi2Int]; }

  inline unsigned int getNStubs(uint32_t hitPattern) { return countSetBits(hitPattern); }

  inline float getMVAQuality(uint32_t mvaQuality) { return TTTrack_TrackWord::tqMVABins[mvaQuality]; }

  inline float getPt(uint32_t rInvInt) {
    return std::abs(MagConstant / getRinv(rInvInt) * BField / 100.0);  // Rinv is in cm-1
  }

  inline float getGlobalPhi(float localPhi, unsigned int sector) {
    float phiCenter = sector * TTTrack_TrackWord::sectorWidth;
    return reco::reducePhiRange(localPhi + phiCenter);
  }

  inline GlobalVector getMomentum(float pt, float phi0, float tanl) {
    return GlobalVector(GlobalVector::Cylindrical(pt, phi0, pt * tanl));
  }

  inline GlobalPoint getPOCA(float d0, float phi0, float z0) {
    return GlobalPoint(d0 * sin(phi0), -d0 * cos(phi0), z0);
  }
}  // namespace l1trkUnpack

#endif
