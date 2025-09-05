#ifndef L1TriggerScouting_Phase2_l1tkemUnpack_h
#define L1TriggerScouting_Phase2_l1tkemUnpack_h
#include <cstdint>
#include <cmath>

// 0 valid - double 1
// 15-1 Rinv - double 15
// 27-16 Phi - double 12
// 31-28 Chi2RPhi - double 4
// 47-32 Tanl - double 16
// 59-48 Z0 - double 12
// 63-60 Chi2RZ - double 4
// 76-64 D0 - double 13
// 79-77 BendChi2 - double 3
// 86-80 HitPattern - unsigned int 7
// 89-87 MVAQuality - double 3
// 95-90 MVAOther - double 6

//TODO - check types, is it all supposed to be double?                                            
namespace l1tkemUnpack {
  inline void read(const uint64_t datalow,
                         const uint32_t datahigh,
                         double &rinv,
                         double &phi,
                         double &chi2RPhi,
                         double &tanl,
                         double &z0,
                         double &chi2Rz,
                         double &d0,
                         double &bendChi2,
                         int16_t &hitPattern,
                         double &mvaQuality,
                         double &MVAOther) {

    //LOGIC HERE
    // 1-  just go to the most significant bit of each property
    // 2- IF STATMENT is it positive or negative? 
    // 2.1- If positive - take the n number of bits of that property as a positive number
    // 2.2- If negative - propagate? the sign. 
    // 3- scaling!

  }
} 

#endif
