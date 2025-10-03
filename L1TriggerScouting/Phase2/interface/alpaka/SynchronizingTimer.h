#ifndef L1TriggerScouting_Phase2_interface_alpaka_SynchronizingTimer_h
#define L1TriggerScouting_Phase2_interface_alpaka_SynchronizingTimer_h

#include <alpaka/alpaka.hpp>
#include <chrono>
#include <fmt/format.h>

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "L1TriggerScouting/Phase2/interface/L1TScPhase2Common.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc {

  using namespace ::l1sc;
  using Clock_t = std::chrono::steady_clock;

  class SynchronizingTimer {
  public:
    SynchronizingTimer(const std::string &label, const Environment environment)
        : label_{label}, environment_{environment} {
      msg_ += fmt::format("[DEBUG] OK - {} ", label);
    }

    void start(Queue &queue) {
      if (environment_ < Environment::kTest)
        return;
      timestamp_ = Clock_t::now();
    }

    void sync(Queue &queue) {
      if (environment_ < Environment::kTest)
        return;

      alpaka::wait(queue);
      auto end = Clock_t::now();
      auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - timestamp_).count();
      auto log = msg_ + fmt::format("({} us)\n", duration);
      fmt::print("{}", log);
    }

  private:
    const std::string label_;
    const Environment environment_;
    Clock_t::time_point timestamp_;
    std::string msg_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc

#endif  // L1TriggerScouting_Phase2_interface_alpaka_SynchronizingTimer_h