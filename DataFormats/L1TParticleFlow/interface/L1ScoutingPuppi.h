#ifndef DataFormats_L1TParticleFlow_L1ScoutingPuppi_h
#define DataFormats_L1TParticleFlow_L1ScoutingPuppi_h

#include <vector>
#include <utility>
#include <cstdint>

namespace l1Scouting {
  struct Puppi {
    float pt, eta, phi, z0, dxy, puppiw;
    int16_t pdgId;
    uint8_t quality;
  };
  struct PuppiSOA {
    std::vector<uint16_t> bx;
    std::vector<uint32_t> offsets;
    std::vector<float> pt, eta, phi, z0, dxy, puppiw;
    std::vector<int16_t> pdgId;
    std::vector<uint8_t> quality;
    PuppiSOA() : bx(), offsets(), pt(), eta(), phi(), z0(), dxy(), puppiw(), pdgId(), quality() {}
    PuppiSOA(const PuppiSOA& other) = default;
    PuppiSOA(PuppiSOA&& other) = default;
    PuppiSOA& operator=(const PuppiSOA& other) = default;
    PuppiSOA& operator=(PuppiSOA&& other) = default;
    void swap(PuppiSOA& other) {
      using std::swap;
      swap(bx, other.bx);
      swap(offsets, other.offsets);
      swap(pt, other.pt);
      swap(eta, other.eta);
      swap(phi, other.phi);
      swap(z0, other.z0);
      swap(dxy, other.dxy);
      swap(puppiw, other.puppiw);
      swap(pdgId, other.pdgId);
      swap(quality, other.quality);
    }
  };
  inline void swap(PuppiSOA& a, PuppiSOA& b) { a.swap(b); }
}
#endif
