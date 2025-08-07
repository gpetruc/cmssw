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
    RecMeson(float pt, float eta, float phi, uint8_t pid, uint8_t id1, uint8_t id2)
      : pt_(pt), eta_(eta), phi_(phi), pid_(pid), id1_(id1), id2_(id2) {}

    float pt() const { return pt_; }
    float eta() const { return eta_; }
    float phi() const { return phi_; }
    uint8_t pid() const { return pid_; }
    int16_t pdgId() const { return PDGID_[pid_]; }
    int id1() const { return id1_; }
    int id2() const { return id2_; }
    float mass() const { return MASS_[pid_]; }
    int charge() const { return (pid_ < 2) ? 0 : (2 * (pid_ & 1) - 1); }

    void setPt(float pt) { pt_ = pt; }
    void setEta(float eta) { eta_ = eta; }
    void setPhi(float phi) { phi_ = phi; }
    void setPid(int8_t pid) { pid_ = pid; }
    void id1(int id1) { id1_ = id1; }
    void id2(int id2) { id2_ = id2; }

    ROOT::Math::PtEtaPhiMVector p4() const { return ROOT::Math::PtEtaPhiMVector(pt_, eta_, phi_, mass()); }
  
    enum PIDs {
      HadZero = 0,
      Gamma = 1,
      HadMinus = 2,
      HadPlus = 3,
      EleMinus = 4,
      ElePlus = 5,
      MuMinus = 6,
      MuPlus = 7,
      nPIDs = 8
    };

  private:
    float pt_, eta_, phi_;
    uint8_t pid_;
    int id1_, id2_;

    static constexpr int16_t PDGID_[nPIDs] = {130, 22, -211, 211, 11, -11, 13, -13};
    static constexpr float MASS_[nPIDs] = {0.5, 0.0, 0.13, 0.13, 0.0005, 0.0005, 0.105, 0.105};
  };
}  // namespace l1Scouting
#endif
