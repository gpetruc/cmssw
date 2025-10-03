#include <cassert>

#include <fmt/format.h>

#include "DataFormats/L1ScoutingSoA/interface/BxLookupHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/PuppiHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/SelectedBxHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/W3PiHostTable.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/stream/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "L1TriggerScouting/Phase2/interface/L1TScPhase2Common.h"

namespace l1sc {

  using namespace cms::alpakatools;

  inline edm::InputTag getBackendTag(edm::InputTag const& tag) {
    return edm::InputTag(tag.label(), "backend", tag.process());
  }

  class L1TScPhase2AlpakaAnalyzer : public edm::stream::EDAnalyzer<> {
  public:
    L1TScPhase2AlpakaAnalyzer(const edm::ParameterSet& params)
        : puppi_token_{consumes(params.getUntrackedParameter<edm::InputTag>("puppi"))},
          bx_lookup_token_{consumes(params.getUntrackedParameter<edm::InputTag>("bx_lookup"))},
          selected_bxs_token_{consumes(params.getUntrackedParameter<edm::InputTag>("selected_bxs"))},
          w3pi_table_token_{consumes(params.getUntrackedParameter<edm::InputTag>("w3pi_table"))},
          puppi_backend_{consumes(getBackendTag(params.getUntrackedParameter<edm::InputTag>("puppi")))},
          bx_lookup_backend_{consumes(getBackendTag(params.getUntrackedParameter<edm::InputTag>("bx_lookup")))},
          selected_bxs_backend_{consumes(getBackendTag(params.getUntrackedParameter<edm::InputTag>("selected_bxs")))},
          w3pi_table_backend_{consumes(getBackendTag(params.getUntrackedParameter<edm::InputTag>("w3pi_table")))},
          environment_{static_cast<Environment>(params.getUntrackedParameter<int>("environment"))},
          fast_path_{params.getParameter<bool>("fast_path")} {}

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      desc.addUntracked<int>("environment", static_cast<int>(Environment::kDevelopment));
      desc.addUntracked<edm::InputTag>("puppi");
      desc.addUntracked<edm::InputTag>("bx_lookup");
      desc.addUntracked<edm::InputTag>("selected_bxs");
      desc.addUntracked<edm::InputTag>("w3pi_table");
      desc.add<bool>("fast_path", false);
      descriptions.addWithDefaultLabel(desc);
    }

    void analyze(edm::Event const& event, edm::EventSetup const&) override {
      if (environment_ >= Environment::kDevelopment) {
        constexpr int total_len = 72;
        auto label = fmt::format("EVENT: {}", event.id().event());
        int pad = total_len - static_cast<int>(label.size());
        int pad_left = pad / 2;
        int pad_right = pad - pad_left - 1;
        fmt::print("\n{0} {1} {2}\n\n", std::string(pad_left, '-'), label, std::string(pad_right, '-'));
      }

      // unpacking/decoding step debugging
      const auto puppi_handle = event.getHandle(puppi_token_);
      const auto bx_lookup_handle = event.getHandle(bx_lookup_token_);
      // w -> 3pi step debugging
      const auto selected_bxs_handle = event.getHandle(selected_bxs_token_);
      const auto w3pi_table_handle = event.getHandle(w3pi_table_token_);

      if (puppi_handle.isValid() && bx_lookup_handle.isValid() && selected_bxs_handle.isValid() &&
          w3pi_table_handle.isValid()) {
        // puppi & bx_lookup
        auto const& puppi = *puppi_handle;
        auto const& bx_lookup = *bx_lookup_handle;
        // w3pi selected bxs
        auto const& selected_bxs = *selected_bxs_handle;
        auto const& w3pi_table = *w3pi_table_handle;
        // update processed bx count
        bx_count_ += bx_lookup.const_view<BxIndexSoA>().metadata().size();
        selected_bx_count_ += countBxMask(selected_bxs.const_view());
        // backends
        auto const puppi_backend = static_cast<Backend>(event.get(puppi_backend_));
        auto const bx_lookup_backend = static_cast<Backend>(event.get(bx_lookup_backend_));
        // auto const selected_bxs_backend = static_cast<Backend>(event.get(selected_bxs_backend_));
        auto const w3pi_table_backend = static_cast<Backend>(event.get(w3pi_table_backend_));
        if (environment_ >= Environment::kDevelopment) {
          logDebug(puppi.const_view(),
                   bx_lookup.const_view<BxIndexSoA>(),
                   bx_lookup.const_view<OffsetsSoA>(),
                   toString(puppi_backend),
                   toString(bx_lookup_backend));
          logDebug(w3pi_table.const_view(),
                   puppi.const_view(),
                   bx_lookup.const_view<BxIndexSoA>(),
                   bx_lookup.const_view<OffsetsSoA>(),
                   toString(w3pi_table_backend));
          logDebug(selected_bxs.const_view(), bx_lookup.const_view<BxIndexSoA>().metadata().size());
        }
      }
    }

    void beginStream(edm::StreamID) override {
      if (environment_ >= Environment::kDevelopment) {
        fmt::print("=========================================================================\n");
      }
    }

    void endStream() override {
      fmt::print("=========================================================================\n");
      fmt::print("[INFO] W3Pi {} -> {}\n", bx_count_, selected_bx_count_);
      fmt::print("=========================================================================\n");
    }

    int countBxMask(const SelectedBxHostCollection::ConstView selected_bxs) {
      if (!fast_path_)
        return selected_bxs.metadata().size();
      int count = 0;
      for (int i = 0; i < selected_bxs.metadata().size(); ++i) {
        if (selected_bxs.bx()[i] > 0)
          count++;
      }
      return count;
    }

    void logDebug(const SelectedBxHostCollection::ConstView selected_bxs, const int nbxs) {
      auto n_selected = countBxMask(selected_bxs);
      if (n_selected == 0)
        return;
      std::vector<int> selected_bxs_indices;
      for (int i = 0; i < selected_bxs.metadata().size(); ++i) {
        if (fast_path_) {
          if (selected_bxs.bx()[i] == 0)
            continue;
          selected_bxs_indices.push_back(i);
        } else {
          selected_bxs_indices.push_back(selected_bxs.bx()[i]);
        }
      }
      fmt::print("[DEBUG] SelectedBxCollection[{}] {} -> {}: {}\n",
                 n_selected,
                 nbxs,
                 n_selected,
                 fmt::join(selected_bxs_indices, ", "));
    }

    void logDebug(const W3PiHostTable::ConstView& table,
                  const PuppiHostCollection::ConstView& puppi,
                  const BxIndexSoA::ConstView& bx_index,
                  const OffsetsSoA::ConstView& offsets,
                  const std::string_view table_backend) {
      const auto size = table.metadata().size();
      if (size == 0)
        return;
      fmt::print("[DEBUG] W3PiTable[{}] ({}):\n", size, table_backend);

      // table header
      const std::string separator = "+---------+---------+---------+---------+-------+";
      fmt::print("{}\n", separator);
      fmt::print("| {:>7} | {:>7} | {:>7} | {:>7} | {:>5} |\n", "idx", "pt", "eta", "phi", "bx");
      fmt::print("{}\n", separator);
      std::vector<int> bxs;
      for (int32_t idx = 0; idx < size; ++idx) {
        const auto& triplet = table[idx];
        auto bx = 0;
        for (int32_t bx_idx = 0; bx_idx < bx_index.metadata().size(); bx_idx++) {
          if (triplet.i() < offsets.offsets()[bx_idx + 1] && triplet.i() >= offsets.offsets()[bx_idx]) {
            bx = bx_idx;
            break;
          }
        }
        bxs.push_back(bx);
        fmt::print("| {:7d} | {:7.2f} | {:7.2f} | {:7.2f} | {:>5} |\n",
                   static_cast<int>(triplet.i()),
                   puppi.pt()[triplet.i()],
                   puppi.eta()[triplet.i()],
                   puppi.phi()[triplet.i()],
                   bx);
        fmt::print("| {:7d} | {:7.2f} | {:7.2f} | {:7.2f} | {:>5} |\n",
                   static_cast<int>(triplet.j()),
                   puppi.pt()[triplet.j()],
                   puppi.eta()[triplet.j()],
                   puppi.phi()[triplet.j()],
                   bx);
        fmt::print("| {:7d} | {:7.2f} | {:7.2f} | {:7.2f} | {:>5} |\n",
                   static_cast<int>(triplet.k()),
                   puppi.pt()[triplet.k()],
                   puppi.eta()[triplet.k()],
                   puppi.phi()[triplet.k()],
                   bx);
        fmt::print("{}\n", separator);
      }
    }

    void logDebug(const PuppiHostCollection::ConstView& puppi,
                  const BxIndexSoA::ConstView& bx_index,
                  const OffsetsSoA::ConstView& offsets,
                  const std::string_view puppi_backend,
                  const std::string_view bx_lookup_backend) {
      const auto size = puppi.metadata().size();
      if (size == 0)
        return;

      // Header
      fmt::print("[DEBUG] PuppiCollection[{}] ({}):\n", size, puppi_backend);

      constexpr auto sep = "+-------+-------+-------+---------+---------+---------+---------+-------+";
      auto printHeader = [&] {
        fmt::print("{}\n", sep);
        fmt::print("| {:>5} | {:>5} | {:>5} | {:>7} | {:>7} | {:>7} | {:>7} | {:>5} |\n",
                   "bx",
                   "range",
                   "local",
                   "global",
                   "pt",
                   "eta",
                   "phi",
                   "pdgid");
        fmt::print("{}\n", sep);
      };

      auto printRow = [&](int bx, int range, int local, int global, const auto& view) {
        fmt::print("| {:5d} | {:5d} | {:5d} | {:7d} | {:7.2f} | {:7.2f} | {:7.2f} | {:5d} |\n",
                   bx,
                   range,
                   local,
                   global,
                   view.pt(),
                   view.eta(),
                   view.phi(),
                   view.pdgid());
      };

      printHeader();

      const int max_entries = (environment_ > Environment::kTest) ? size : 5;
      int printed = 0;

      // Print first 5 entries until max_entries
      for (int i = 0; i < bx_index.metadata().size() && printed < max_entries; ++i) {
        const int bx = bx_index.bx()[i];
        const int start = offsets.offsets()[i];
        const int end = offsets.offsets()[i + 1];
        const int range = end - start;

        for (int j = 0; j < range && printed < max_entries; ++j) {
          const int globalIdx = start + j;
          if (globalIdx >= size)
            break;
          printRow(bx, range, j, globalIdx, puppi[globalIdx]);
          ++printed;
        }
      }

      // Ellipsis row if not all printed
      if (printed < size) {
        fmt::print("| {:>5} | {:>5} | {:>5} | {:>7} | {:>7} | {:>7} | {:>7} | {:>5} |\n",
                   "...",
                   "...",
                   "...",
                   "...",
                   "...",
                   "...",
                   "...",
                   "...");
        // Print last 5 entries of the last BX only
        const int lastBxIdx = bx_index.metadata().size() - 1;
        const int start = offsets.offsets()[lastBxIdx];
        const int end = offsets.offsets()[lastBxIdx + 1];
        const int range = end - start;
        const int n_last = std::min(5, range);

        for (int j = range - n_last; j < range; ++j) {
          const int globalIdx = start + j;
          printRow(bx_index.bx()[lastBxIdx], range, j, globalIdx, puppi[globalIdx]);
        }
      }

      fmt::print("{}\n", sep);
    }

  private:
    // get products
    const edm::EDGetTokenT<PuppiHostCollection> puppi_token_;
    const edm::EDGetTokenT<BxLookupHostCollection> bx_lookup_token_;
    const edm::EDGetTokenT<SelectedBxHostCollection> selected_bxs_token_;
    const edm::EDGetTokenT<W3PiHostTable> w3pi_table_token_;

    // backend query
    const edm::EDGetTokenT<unsigned short> puppi_backend_;
    const edm::EDGetTokenT<unsigned short> bx_lookup_backend_;
    const edm::EDGetTokenT<unsigned short> selected_bxs_backend_;
    const edm::EDGetTokenT<unsigned short> w3pi_table_backend_;

    // debug
    const Environment environment_;
    const bool fast_path_;
    int selected_bx_count_ = 0;
    int bx_count_ = 0;
  };

}  // namespace l1sc

DEFINE_FWK_MODULE(l1sc::L1TScPhase2AlpakaAnalyzer);