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
    TTrack(float pt, float eta, float phi, float z0, float dxy, float mvaQuality, uint8_t nStub, uint8_t quality, int8_t charge)
         : pt_(pt), eta_(eta), phi_(phi), z0_(z0), dxy_(dxy), mvaQuality_(mvaQuality), nStub_(nStub), quality_(quality), charge_(charge) {}

    float pt() const { return pt_; }
    float eta() const { return eta_; }
    float phi() const { return phi_; }
    float z0() const { return z0_; }
    float dxy() const { return dxy_; }
    float mvaQuality() const { return mvaQuality_; }
    uint8_t nStub() const { return nStub_; }
    uint8_t quality() const { return quality_; }
    int8_t charge() const { return charge_; }

  private:
    float pt_, eta_, phi_, z0_, dxy_, mvaQuality_;
    uint8_t nStub_, quality_;
    int8_t charge_;
  };

  struct TTrackSOA {
    std::vector<uint16_t> bx;
    std::vector<uint32_t> offsets;
    std::vector<float> pt, eta, phi, z0, dxy, mvaQuality;
    std::vector<uint8_t> nStub, quality;
    std::vector<int8_t> charge;
    TTrackSOA() : bx(), offsets(), pt(), eta(), phi(), z0(), dxy(), mvaQuality(), nStub(), quality(), charge() {}
    TTrackSOA(const TTrackSOA& other) = default;
    TTrackSOA(TTrackSOA&& other) = default;
    TTrackSOA& operator=(const TTrackSOA& other) = default;
    TTrackSOA& operator=(TTrackSOA&& other) = default;
  };
};
#endif
