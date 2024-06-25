#include "DataFormats/L1Scouting/interface/OrbitFlatTable.h"

l1ScoutingRun3::OrbitFlatTable::OrbitFlatTable(std::vector<unsigned> bxOffsets,
                                               const std::string &name,
                                               bool singleton,
                                               bool extension)
    : nanoaod::FlatTable(bxOffsets.back(), name, singleton, extension), bxOffsets_(bxOffsets) {
  if (bxOffsets.size() != orbitBufferSize_ + 1) {
    throw cms::Exception("LogicError") <<
                 "Mismatch between bxOffsets.size() " << bxOffsets.size() << 
                 " and orbitBufferSize_ + 1" << (orbitBufferSize_ + 1);
  }
}

void l1ScoutingRun3::OrbitFlatTable::throwBadBx(unsigned bx) const {
  throw cms::Exception("OrbitFlatTable")
            << "Trying to access bad bx " << bx;
}
