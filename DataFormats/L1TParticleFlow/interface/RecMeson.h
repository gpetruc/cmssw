#ifndef DataFormats_L1TParticleFlow_RecMeson_h
#define DataFormats_L1TParticleFlow_RecMeson_h

#include <vector>
#include <array>
#include <utility>
#include <cstdint>
#include <Math/Vector4D.h>

namespace l1Scouting {
  class RecMesonBase {
  public:
    RecMesonBase() {}
    RecMesonBase(float pt, float eta, float phi, float mass, int charge, int pdgId, float isoDR0p25)
        : pt_(pt), eta_(eta), phi_(phi), mass_(mass), charge_(charge), pdgId_(pdgId), isoDR0p25_(isoDR0p25) {}

    float pt() const { return pt_; }
    float eta() const { return eta_; }
    float phi() const { return phi_; }
    float mass() const { return mass_; }
    int charge() const { return charge_; }
    int pdgId() const { return pdgId_; }
    float isoDR0p25() const { return isoDR0p25_; }
    ROOT::Math::PtEtaPhiMVector p4() const { return ROOT::Math::PtEtaPhiMVector(pt_, eta_, phi_, mass_); }

    void setPt(float pt) { pt_ = pt; }
    void setEta(float eta) { eta_ = eta; }
    void setPhi(float phi) { phi_ = phi; }
    void setMass(float mass) { mass_ = mass; }
    void setCharge(int charge) { charge_ = charge; }
    void setPdgId(int pdgId) { pdgId_ = pdgId; }
    void setIsoDR0p25(float isoDR0p25) { isoDR0p25_ = isoDR0p25; }

  private:
    float pt_, eta_, phi_, mass_;
    int charge_;
    int pdgId_;
    float isoDR0p25_;
  };

  template <std::size_t N>
  class RecMeson : public RecMesonBase {
  public:
    RecMeson() {}
    RecMeson(float pt,
             float eta,
             float phi,
             float mass,
             int charge,
             int pdgId,
             float isoDR0p25,
             const std::array<double, N>& daughterMasses,
             const std::array<unsigned int, N>& daughterIds)
        : RecMesonBase(pt, eta, phi, mass, charge, pdgId, isoDR0p25),
          daughterMasses_(daughterMasses),
          daughterIds_(daughterIds) {}

    const std::array<double, N>& daughterMasses() const { return daughterMasses_; }
    const std::array<unsigned int, N>& daughterIds() const { return daughterIds_; }

    double daughterMass(size_t i) const { return daughterMasses_.at(i); }
    unsigned int daughterIds(size_t i) const { return daughterIds_.at(i); }

    void setDaughterMasses(const std::array<double, N>& m) { daughterMasses_ = m; }
    void setDaughterIds(const std::array<unsigned int, N>& ids) { daughterIds_ = ids; }

  private:
    std::array<double, N> daughterMasses_;
    std::array<unsigned int, N> daughterIds_;
  };
}  // namespace l1Scouting

#endif
