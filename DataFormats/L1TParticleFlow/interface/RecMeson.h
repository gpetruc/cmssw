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
    RecMeson(float pt, int id)
        : pt_(pt), id_(id) {}

    float pt() const { return pt_; }
    int id() const { return id_; }

    void setPt(float pt) { pt_ = pt; }
    void setId(int id) { id_ = id; }

  private:
    float pt_;
    int id_;
  };
}  // namespace l1Scouting
#endif
