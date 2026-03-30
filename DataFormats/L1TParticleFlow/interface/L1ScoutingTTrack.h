#ifndef DataFormats_L1TParticleFlow_L1ScoutingTTrack_h
#define DataFormats_L1TParticleFlow_L1ScoutingTTrack_h

#include <cstdint>

namespace l1Scouting {
  class TTrack {
  public:
    TTrack() {}
    TTrack(float pt,
           float eta,
           float phi,
           float z0,
           float dxy,
           float mvaQuality,
           uint8_t nStub,
           uint8_t quality,
           int8_t charge)
        : pt_(pt),
          eta_(eta),
          phi_(phi),
          z0_(z0),
          dxy_(dxy),
          mvaQuality_(mvaQuality),
          nStub_(nStub),
          quality_(quality),
          charge_(charge) {}

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
};  // namespace l1Scouting
#endif
