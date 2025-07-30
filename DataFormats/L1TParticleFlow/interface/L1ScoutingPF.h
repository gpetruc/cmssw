#ifndef DataFormats_L1TParticleFlow_L1ScoutingPF_h
#define DataFormats_L1TParticleFlow_L1ScoutingPF_h

#include <vector>
#include <utility>
#include <cstdint>
#include <Math/Vector4D.h>

namespace l1Scouting {
  class PF {
  public:
    PF() {}
    PF(float pt, float eta, float phi, uint8_t pid, float z0, float dxy, float pfw, uint8_t quality)
        : pt_(pt), eta_(eta), phi_(phi), z0_(z0), dxy_(dxy), pfw_(pfw), pid_(pid), quality_(quality) {}
    PF(float pt, float eta, float phi, uint8_t pid, float z0, float dxy, uint8_t quality)
        : pt_(pt), eta_(eta), phi_(phi), z0_(z0), dxy_(dxy), pfw_(1.0f), pid_(pid), quality_(quality) {}
    PF(float pt, float eta, float phi, uint8_t pid, float pfw, uint8_t quality)
        : pt_(pt), eta_(eta), phi_(phi), z0_(0.0f), dxy_(0.0f), pfw_(pfw), pid_(pid), quality_(quality) {}

    float pt() const { return pt_; }
    float eta() const { return eta_; }
    float phi() const { return phi_; }
    float z0() const { return z0_; }
    float dxy() const { return dxy_; }
    float pfw() const { return pfw_; }
    uint8_t pid() const { return pid_; }
    int16_t pdgId() const { return PDGID_[pid_]; }
    uint8_t quality() const { return quality_; }
    float mass() const { return MASS_[pid_]; }
    int charge() const { return (pid_ < 2) ? 0 : (2 * (pid_ & 1) - 1); }

    void setPt(float pt) { pt_ = pt; }
    void setEta(float eta) { eta_ = eta; }
    void setPhi(float phi) { phi_ = phi; }
    void setZ0(float z0) { z0_ = z0; }
    void setDxy(float dxy) { dxy_ = dxy; }
    void setPFw(float pfw) { pfw_ = pfw; }
    void setPid(int8_t pid) { pid_ = pid; }
    void setQuality(uint8_t quality) { quality_ = quality; }

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
    float pt_, eta_, phi_, z0_, dxy_, pfw_;
    uint8_t pid_, quality_;

    static constexpr int16_t PDGID_[nPIDs] = {130, 22, -211, 211, 11, -11, 13, -13};
    static constexpr float MASS_[nPIDs] = {0.5, 0.0, 0.13, 0.13, 0.0005, 0.0005, 0.105, 0.105};
  };

  struct PFSOA {
    std::vector<uint16_t> bx;
    std::vector<uint32_t> offsets;
    std::vector<float> pt, eta, phi, z0, dxy, pfw;
    std::vector<int16_t> pdgId;
    std::vector<uint8_t> quality;
    PFSOA() : bx(), offsets(), pt(), eta(), phi(), z0(), dxy(), pfw(), pdgId(), quality() {}
    PFSOA(const PFSOA& other) = default;
    PFSOA(PFSOA&& other) = default;
    PFSOA& operator=(const PFSOA& other) = default;
    PFSOA& operator=(PFSOA&& other) = default;
    void swap(PFSOA& other) {
      using std::swap;
      swap(bx, other.bx);
      swap(offsets, other.offsets);
      swap(pt, other.pt);
      swap(eta, other.eta);
      swap(phi, other.phi);
      swap(z0, other.z0);
      swap(dxy, other.dxy);
      swap(pfw, other.pfw);
      swap(pdgId, other.pdgId);
      swap(quality, other.quality);
    }
  };
  inline void swap(PFSOA& a, PFSOA& b) { a.swap(b); }
}  // namespace l1Scouting
#endif
