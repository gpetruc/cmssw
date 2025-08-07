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
    RecMeson(float pt, float eta, float phi, float mass, int charge, float dmass1, float dmass2, int pdgId, int id1, int id2, float isoDR0p25)
      : pt_(pt), eta_(eta), phi_(phi), mass_(mass), charge_(charge), dmass1_(dmass1), dmass2_(dmass2), pdgId_(pdgId), id1_(id1), id2_(id2), isoDR0p25_(isoDR0p25) {}

    float pt() const { return pt_; }
    float eta() const { return eta_; }
    float phi() const { return phi_; }
    float mass() const { return mass_; }
    int charge() const { return charge_; }
    float dmass1() const { return dmass1_; }
    float dmass2() const { return dmass2_; }
    int pdgId() const { return pdgId_; }
    int id1() const { return id1_; }
    int id2() const { return id2_; }
    float isoDR0p25() const { return isoDR0p25_; }
    ROOT::Math::PtEtaPhiMVector p4() const { return ROOT::Math::PtEtaPhiMVector(pt_, eta_, phi_, mass()); }

    void setPt(float pt) { pt_ = pt; }
    void setEta(float eta) { eta_ = eta; }
    void setPhi(float phi) { phi_ = phi; }
    void setMass(float mass) { mass_ = mass; }
    void setCharge(int charge) { charge_ = charge; }
    void setDmass1(float dmass1) { dmass1_ = dmass1; }
    void setDmass2(float dmass2) { dmass2_ = dmass2; }
    void setPdgId(int pdgId) { pdgId_ = pdgId; }
    void setId1(int id1) { id1_ = id1; }
    void setId2(int id2) { id2_ = id2; }
    void setIsoDR0p25(float isoDR0p25) { isoDR0p25_ = isoDR0p25; }


  private:
    float pt_, eta_, phi_, mass_;
    int charge_;
    float dmass1_, dmass2_;
    int pdgId_;
    int id1_, id2_;
    float isoDR0p25_;
  };
}  // namespace l1Scouting
#endif
