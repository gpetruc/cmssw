#ifndef L1TriggerScouting_Phase2_l1tkemUnpack_h
#define L1TriggerScouting_Phase2_l1tkemUnpack_h
#include <cstdint>
#include <cmath>

// Section 6.3.3.2 https://www.overleaf.com/project/5fbecbed9fa10637655869e2
// 0 valid
// 16-1 pt
// 29-17 phi
// 43-30 eta
// 47-44 quality
// 58-48 isolation
// 95-59 unassigned

namespace l1tkemUnpack {
  template <typename U>
  inline void parseHeader(const uint64_t &header, uint16_t &run, uint16_t &bx, uint32_t &orbit, bool &good, U &npackets) {
    npackets = header & 0xFFF; // lenght of the packets in 64bit words // FIXME is 12 bits correct?
    bx = (header >> 12) & 0xFFF;
    orbit = (header >> 24) & 0X3FFFFFFF;
    run = (header >> 54);
    good = !(header & (1llu << 61));
  } 

  inline void readshared(const uint64_t data, uint16_t &pt, int16_t &eta, int16_t &phi, bool &valid, uint8_t &quality, uint16_t isolation) {  //int
    valid = data & 0x001;                                                                       // 1 bit
    pt = ((data >> 16) & 1) ? ((data >> 1) | (-0x8000)) : ((data >> 1) & (0xFFFF));             // 16 bits
    phi = ((data >> 29) & 1) ? ((data >> 17) | (-0x1000)) : ((data >> 17) & (0x1FFF));          // 13 bits
    eta = ((data >> 43) & 1) ? ((data >> 30) | (-0x2000)) : ((data >> 30) & (0x3FFF));          // 14 bits
    quality = ((data >> 47) & 1) ? ((data >> 44) | (-0x8)) : ((data >> 44) & (0xF));            // 4 bits
    isolation = ((data >> 58) & 1) ? ((data >> 48) | (-0x400)) : ((data >> 48) & (0x7FF));     // 11 bits

  }
  inline void readshared(const uint64_t data, float &pt, float &eta, float &phi, bool &valid, uint8_t &quality, float isolation) {  //float

    valid = data & 0x001;

    uint16_t ptint = ((data >> 16) & 1) ? ((data >> 1) | (-0x8000)) : ((data >> 1) & (0xFFFF));
    pt = ptint * 0.03125f;

    int phiint = ((data >> 29) & 1) ? ((data >> 17) | (-0x1000)) : ((data >> 17) & (0x1FFF));
    phi = phiint * float(M_PI / 4096.);

    int etaint = ((data >> 43) & 1) ? ((data >> 30) | (-0x2000)) : ((data >> 30) & (0x3FFF));
    eta = etaint * float(M_PI / 4096.);

    quality = ((data >> 47) & 1) ? ((data >> 44) | (-0x8)) : ((data >> 44) & (0xF));

    int isolationint = ((data >> 58) & 1) ? ((data >> 48) | (-0x400)) : ((data >> 48) & (0x7FF));
    isolation = isolationint * 0.25f;

  }

}  // namespace l1tkemUnpack

#endif
