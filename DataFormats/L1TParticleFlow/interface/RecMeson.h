#ifndef DataFormats_L1TParticleFlow_RecMeson_h
#define DataFormats_L1TParticleFlow_RecMeson_h

#include <vector>
#include <utility>
#include <cstdint>
#include <Math/Vector4D.h>

namespace l1Scouting {
  class RecMeson {
  public:
    RecMeson() {}
    RecMeson(float pt, float eta, float phi, uint8_t id1, uint8_t id2)
      : pt_(pt), eta_(eta), phi_(phi), id1_(id1), id2_(id2) {}

    float pt() const { return pt_; }
    float eta() const { return eta_; }
    float phi() const { return phi_; }
    int id1() const { return id1_; }
    int id2() const { return id2_; }

    void setPt(float pt) { pt_ = pt; }
    void setEta(float eta) { eta_ = eta; }
    void setPhi(float phi) { phi_ = phi; }
    void id1(int id1) { id1_ = id1; }
    void id2(int id2) { id2_ = id2; }

  private:
    float pt_;
    float eta_;
    float phi_;
    int id1_;
    int id2_;
  };
}  // namespace l1Scouting
#endif
