#ifndef L1TriggerScouting_Phase2_l1puppiUnpack_h
#define L1TriggerScouting_Phase2_l1puppiUnpack_h
#include "DataFormats/L1TParticleFlow/interface/PFCandidate.h"
#include <cstdint>
#include <cmath>

namespace l1puppiUnpack {
  template <typename U>
  inline void parseHeader(const uint64_t &header, uint16_t &run, uint16_t &bx, uint32_t &orbit, bool &good, U &npuppi) {
    npuppi = header & 0xFFF;
    bx = (header >> 12) & 0xFFF;
    orbit = (header >> 24) & 0X3FFFFFFF;
    run = (header >> 54);
    good = !(header & (1llu << 61));
  }

  inline void assignpdgid(uint8_t pid, short int &pdgid) {
    static constexpr int16_t PDGIDS[8] = {130, 22, -211, 211, 11, -11, 13, -13};
    pdgid = PDGIDS[pid];
  }
  inline void vassignpdgid(uint8_t pid, short int &pdgid) {
    // vectorizable version, ugly as it is...
    short int pdgId = pid ? 22 : 130;
    if (pid > 1) {  // charged
      if (pid / 2 == 1)
        pdgId = -211;
      else if (pid / 2 == 2)
        pdgId = 11;
      else
        pdgId = 13;
      if (pid & 1)
        pdgId = -pdgId;
    }
    pdgid = pdgId;
  }
  inline void assignCMSSWPFCandidateId(uint8_t pid, l1t::PFCandidate::ParticleType &id) {
    static constexpr l1t::PFCandidate::ParticleType PFIDS[8] = {l1t::PFCandidate::NeutralHadron,
                                                                l1t::PFCandidate::Photon,
                                                                l1t::PFCandidate::ChargedHadron,
                                                                l1t::PFCandidate::ChargedHadron,
                                                                l1t::PFCandidate::Electron,
                                                                l1t::PFCandidate::Electron,
                                                                l1t::PFCandidate::Muon,
                                                                l1t::PFCandidate::Muon};
    id = PFIDS[pid];
  }
  inline void assignmass(uint8_t pid, float &mass) {
    static constexpr float MASSES[8] = {0.5, 0.0, 0.13, 0.13, 0.0005, 0.0005, 0.105, 0.105};
    mass = MASSES[pid];
  }
  inline void assigncharge(uint8_t pid, int &charge) { charge = (pid > 1 ? (pid & 1 ? +1 : -1) : 0); }
  inline bool readpid(const uint64_t data, short int &pdgid) {
    uint8_t pid = (data >> 37) & 0x7;
    assignpdgid(pid, pdgid);
    return (pid > 1);
  }
  inline void readshared(const uint64_t data, uint16_t &pt, int16_t &eta, int16_t &phi) {  //int
    pt = data & 0x3FFF;
    eta = ((data >> 25) & 1) ? ((data >> 14) | (-0x800)) : ((data >> 14) & (0xFFF));
    phi = ((data >> 36) & 1) ? ((data >> 26) | (-0x400)) : ((data >> 26) & (0x7FF));
  }
  inline void readshared(const uint64_t data, float &pt, float &eta, float &phi) {  //float
    uint16_t ptint = data & 0x3FFF;
    pt = ptint * 0.25f;

    int etaint = ((data >> 25) & 1) ? ((data >> 14) | (-0x800)) : ((data >> 14) & (0xFFF));
    eta = etaint * float(M_PI / 720.);

    int phiint = ((data >> 36) & 1) ? ((data >> 26) | (-0x400)) : ((data >> 26) & (0x7FF));
    phi = phiint * float(M_PI / 720.);
  }
  inline void readcharged(const uint64_t data, int16_t &z0, int8_t &dxy, uint8_t &quality) {  //int
    z0 = ((data >> 49) & 1) ? ((data >> 40) | (-0x200)) : ((data >> 40) & 0x3FF);
    dxy = ((data >> 57) & 1) ? ((data >> 50) | (-0x100)) : ((data >> 50) & 0xFF);
    quality = (data >> 58) & 0x7;  //3 bits
  }
  inline void readcharged(const uint64_t data, float &z0, float &dxy, uint8_t &quality) {  //float
    int z0int = ((data >> 49) & 1) ? ((data >> 40) | (-0x200)) : ((data >> 40) & 0x3FF);
    z0 = z0int * .05f;  //conver to centimeters

    int dxyint = ((data >> 57) & 1) ? ((data >> 50) | (-0x100)) : ((data >> 50) & 0xFF);
    dxy = dxyint * 0.05f;          // PLACEHOLDER
    quality = (data >> 58) & 0x7;  //3 bits
  }
  inline void readcharged(const uint64_t data, uint8_t pid, float &z0, float &dxy) {  //float
    int z0int = ((data >> 49) & 1) ? ((data >> 40) | (-0x200)) : ((data >> 40) & 0x3FF);
    z0 = (pid > 1) * z0int * .05f;  //conver to centimeters
    int dxyint = ((data >> 57) & 1) ? ((data >> 50) | (-0x100)) : ((data >> 50) & 0xFF);
    dxy = (pid > 1) * dxyint * 0.05f;  // PLACEHOLDER
  }
  inline void readneutral(const uint64_t data, uint16_t &wpuppi, uint8_t &id) {
    wpuppi = (data >> 40) & 0x3FF;
    id = (data >> 50) & 0x3F;
  }
  inline void readneutral(const uint64_t data, float &wpuppi, uint8_t &id) {
    int wpuppiint = (data >> 40) & 0x3FF;
    wpuppi = wpuppiint * (1 / 256.f);
    id = (data >> 50) & 0x3F;
  }
  inline void readneutral(const uint64_t data, uint8_t pid, float &wpuppi) {
    int wpuppiint = (data >> 40) & 0x3FF;
    wpuppi = pid > 1 ? wpuppiint * float(1 / 256.f) : 1.0f;
  }
  inline void readquality(const uint64_t data, uint8_t pid, uint8_t &quality) {
    quality = pid > 1 ? (data >> 58) & 0x7 : (data >> 50) & 0x3F;
  }


  template <unsigned int start, unsigned int bits = 16, typename T>
     inline uint16_t extractBitsFromW(const T word) {
    return (word >> start) & ((1 << bits) - 1);
  }
  template <unsigned int start, unsigned int bits = 16, typename T>
    inline int16_t extractSignedBitsFromW(const T word) {
    uint16_t raw = extractBitsFromW<start, bits>(word);
    if ((bits < 16) && (raw & (1 << (bits - 1)))) {
      constexpr uint16_t ormask = (1 << 16) - (1 << bits);
      raw |= ormask;
    }
    return raw;
  }



}  // namespace l1puppiUnpack

#endif
