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
    TkEm(float pt, float eta, float phi, uint8_t quality, float isolation)
        : pt_(pt), eta_(eta), phi_(phi), quality_(quality), isolation_(isolation), id_(0) {}
    TkEm(float pt, float eta, float phi, uint8_t quality, float isolation, int8_t id)
        : pt_(pt), eta_(eta), phi_(phi), quality_(quality), isolation_(isolation), id_(id) {}

    float pt() const { return pt_; }
    float eta() const { return eta_; }
    float phi() const { return phi_; }
    uint8_t quality() const { return quality_; }
    float isolation() const { return isolation_; }
    int8_t id() const { return id_; }

    void setId(int8_t id) { id_ = id; }
    void setPt(float pt) { pt_ = pt; }
    void setEta(float eta) { eta_ = eta; }
    void setPhi(float phi) { phi_ = phi; }
    void setQuality(uint8_t quality) { quality_ = quality; }
    void setIsolation(float isolation) { isolation_ = isolation; }

    ROOT::Math::PtEtaPhiMVector p4() const { return ROOT::Math::PtEtaPhiMVector(pt_, eta_, phi_, 0.0); }

  private:
    float pt_, eta_, phi_;
    uint8_t quality_;
    float isolation_;
    int8_t id_;
  };

  class TkEle : public TkEm {
  public:
    TkEle() {}
    TkEle(float pt, float eta, float phi, uint8_t quality, float isolation, int8_t charge, float z0)
        : TkEm(pt, eta, phi, quality, isolation), charge_(charge), z0_(z0) {}

    int8_t charge() const { return charge_; }
    float z0() const { return z0_; }

    void setZ0(float z0) { z0_ = z0; }
    void setCharge(int8_t charge) { charge_ = charge; }

  private:
    int8_t charge_;
    float z0_;
  };
}  // namespace l1Scouting
#endif
