#ifndef DataFormats_L1TParticleFlow_L1ScoutingTkEm_h
#define DataFormats_L1TParticleFlow_L1ScoutingTkEm_h

#include <vector>
#include <utility>
#include <cstdint>
#include <Math/Vector4D.h>

namespace l1Scouting {
  class TkEm {
  public:
    TkEm() {}
    TkEm(float pt, float eta, float phi, bool valid, uint8_t quality, float isolation) 
        : pt_(pt), eta_(eta), phi_(phi), valid_(valid), quality_(quality), isolation_(isolation) {}

    float pt() const { return pt_; }
    float eta() const { return eta_; }
    float phi() const { return phi_; }
    bool valid() const { return valid_; }
    uint8_t quality() const { return quality_; }
    float isolation() const { return isolation_; }

    void setPt(float pt) { pt_ = pt; }
    void setEta(float eta) { eta_ = eta; }
    void setPhi(float phi) { phi_ = phi; }
    void setValid(bool valid) { valid_ = valid; }
    void setQuality(uint8_t quality) { quality_ = quality; }
    void setIsolation(float isolation) { isolation_ = isolation; }

    ROOT::Math::PtEtaPhiMVector p4() const { return ROOT::Math::PtEtaPhiMVector(pt_, eta_, phi_, 0.0); }

  private:
    float pt_, eta_, phi_;
    bool valid_;
    uint8_t quality_;
    float isolation_;

  };

}  // namespace l1Scouting
#endif
