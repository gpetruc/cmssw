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
    TkEm(float pt, float eta, float phi, bool valid, uint8_t quality, float isolation) //, uint8_t pid, float z0, float dxy, float puppiw, uint8_t quality)
        : pt_(pt), eta_(eta), phi_(phi), valid_(valid), quality_(quality), isolation_(isolation) {} //, z0_(z0), dxy_(dxy), puppiw_(puppiw), pid_(pid), quality_(quality) {}
    //TkEm(float pt, float eta, float phi) //, uint8_t pid, float z0, float dxy, uint8_t quality)
    //    : pt_(pt), eta_(eta), phi_(phi) {}//, z0_(z0), dxy_(dxy), puppiw_(1.0f), pid_(pid), quality_(quality) {}
    //TkEm(float pt, float eta, float phi) //, uint8_t pid, float puppiw, uint8_t quality)
    //    : pt_(pt), eta_(eta), phi_(phi) {}//, z0_(0.0f), dxy_(0.0f), puppiw_(puppiw), pid_(pid), quality_(quality) {}

    float pt() const { return pt_; }
    float eta() const { return eta_; }
    float phi() const { return phi_; }
    bool valid() const { return valid_; }
    uint8_t quality() const { return quality_; }
    float isolation() const { return isolation_; }
    /*float z0() const { return z0_; }
    float dxy() const { return dxy_; }
    float puppiw() const { return puppiw_; }
    uint8_t pid() const { return pid_; }
    int16_t pdgId() const { return PDGID_[pid_]; }
    uint8_t quality() const { return quality_; }
    float mass() const { return MASS_[pid_]; }
    int charge() const { return (pid_ < 2) ? 0 : (2 * (pid_ & 1) - 1); }*/

    void setPt(float pt) { pt_ = pt; }
    void setEta(float eta) { eta_ = eta; }
    void setPhi(float phi) { phi_ = phi; }
    void setValid(bool valid) { valid_ = valid; }
    void setQuality(uint8_t quality) { quality_ = quality; }
    void setIsolation(float isolation) { isolation_ = isolation; }
    /*void setZ0(float z0) { z0_ = z0; }
    void setDxy(float dxy) { dxy_ = dxy; }
    void setPuppiw(float puppiw) { puppiw_ = puppiw; }
    void setPid(int8_t pid) { pid_ = pid; }
    void setQuality(uint8_t quality) { quality_ = quality; }*/

    //ROOT::Math::PtEtaPhiMVector p4() const { return ROOT::Math::PtEtaPhiMVector(pt_, eta_, phi_, mass()); }
    ROOT::Math::PtEtaPhiMVector p4() const { return ROOT::Math::PtEtaPhiMVector(pt_, eta_, phi_, 0.0); }

    /*enum PIDs {
      HadZero = 0,
      Gamma = 1,
      HadMinus = 2,
      HadPlus = 3,
      EleMinus = 4,
      ElePlus = 5,
      MuMinus = 6,
      MuPlus = 7,
      nPIDs = 8
    };*/

  private:
    float pt_, eta_, phi_;//, z0_, dxy_, puppiw_;
    bool valid_;
    uint8_t quality_;
    float isolation_;
    //uint8_t pid_, quality_;

    //static constexpr int16_t PDGID_[nPIDs] = {130, 22, -211, 211, 11, -11, 13, -13};
    //static constexpr float MASS_[nPIDs] = {0.5, 0.0, 0.13, 0.13, 0.0005, 0.0005, 0.105, 0.105};
  };

}  // namespace l1Scouting
#endif
