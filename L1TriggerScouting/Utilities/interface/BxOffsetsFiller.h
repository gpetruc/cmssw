#ifndef L1TriggerScouting_Utilities_BxOffsetsFiller_h
#define L1TriggerScouting_Utilities_BxOffsetsFiller_h

#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/L1Scouting/interface/OrbitFlatTable.h"
#include <vector>
#include <memory>

namespace l1ScoutingRun3 {
  class BxOffsetsFillter {
  public:
    static constexpr unsigned NBX = OrbitFlatTable::NBX;

    BxOffsetsFillter() : lastBx_(0), bxOffsets_(1, 0) {}

    void start() {
      lastBx_ = 0;
      bxOffsets_.resize(1, 0);
      bxOffsets_[0] = 0;
      bxOffsets_.reserve(NBX + 2);
    }

    void addBx(unsigned bx, unsigned size) {
      if (bx == 0 || bx < lastBx_ || bx > NBX)
        throw cms::Exception("LogicError") << "Bad BX or filled out of order, bx " << bx << ", last bx " << lastBx_;
      unsigned current = bxOffsets_.back(), end = current + size;
      bxOffsets_.resize(bx + 1, current);
      bxOffsets_.emplace_back(end);
    }

    std::vector<unsigned>&& done() {
      bxOffsets_.resize(NBX + 2, bxOffsets_.back());
      return std::move(bxOffsets_);
    }

  private:
    unsigned lastBx_;
    std::vector<unsigned> bxOffsets_;
  };
}  // namespace l1ScoutingRun3

#endif