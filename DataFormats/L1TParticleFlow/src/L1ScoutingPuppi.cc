#include "DataFormats/L1TParticleFlow/interface/L1ScoutingPuppi.h"

float l1Scouting::Puppi::mass() const {
    switch(std::abs(pdgId_)) {
        case 22: return 0.f;
        case 11: return 0.0005f;
        case 13: return 0.105f;
        case 130: return 0.5f;
        case 211: return 0.13f;
        default:
            return 0.f;
    }
}


int l1Scouting::Puppi::charge() const {
    switch(pdgId_) {
        case +11: return -1;
        case +13: return -1;
        case -11: return +1;
        case -13: return +1;
        case +211: return +1;
        case -211: return -1;
        default:
            return 0;
    }
}