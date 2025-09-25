#include <cassert>

#include <fmt/format.h>

#include "DataFormats/L1ScoutingSoA/interface/HostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/HostObject.h"
#include "DataFormats/L1ScoutingSoA/interface/CounterHost.h"
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

namespace l1sc {

  inline edm::InputTag getBackendTag(edm::InputTag const& tag) {
    return edm::InputTag(tag.label(), "backend", tag.process());
  }

  class L1TScPhase2W3PiAnalyzer : public edm::stream::EDAnalyzer<> {
  public:
    L1TScPhase2W3PiAnalyzer(const edm::ParameterSet& params)
        : verbose_{params.getUntrackedParameter<bool>("verbose")},
          verbose_level_(params.getUntrackedParameter<int>("verboseLevel")),
          puppi_token_{consumes(params.getUntrackedParameter<edm::InputTag>("puppi"))},
          puppi_backend_{consumes(getBackendTag(params.getUntrackedParameter<edm::InputTag>("puppi")))},
          nbx_map_token_{consumes(params.getUntrackedParameter<edm::InputTag>("nbx_map"))},
          nbx_map_backend_{consumes(getBackendTag(params.getUntrackedParameter<edm::InputTag>("nbx_map")))},
          table_token_{consumes(params.getUntrackedParameter<edm::InputTag>("table"))},
          table_backend_{consumes(getBackendTag(params.getUntrackedParameter<edm::InputTag>("table")))},
          bx_ct_token_{consumes(params.getUntrackedParameter<edm::InputTag>("bx_ct"))} {}

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      desc.addUntracked<bool>("verbose", false);
      desc.addUntracked<int>("verboseLevel", 0);
      desc.addUntracked<edm::InputTag>("puppi");
      desc.addUntracked<edm::InputTag>("nbx_map");
      desc.addUntracked<edm::InputTag>("table");
      desc.addUntracked<edm::InputTag>("bx_ct");
      descriptions.addWithDefaultLabel(desc);
    }

    void analyze(edm::Event const& event, edm::EventSetup const&) override {
      auto const& puppi_collection = event.get(puppi_token_);
      auto const puppi_collection_backend = static_cast<cms::alpakatools::Backend>(event.get(puppi_backend_));

      auto const& nbx_map_collection = event.get(nbx_map_token_);
      auto const nbx_map_collection_backend = static_cast<cms::alpakatools::Backend>(event.get(nbx_map_backend_));

      auto const& table_collection = event.get(table_token_);
      auto const table_collection_backend = static_cast<cms::alpakatools::Backend>(event.get(table_backend_));

      auto const &bx_ct = event.get(bx_ct_token_);

      if (verbose_) {
        nbx_processed_ += nbx_map_collection.view<NbxSoA>().metadata().size();
        debugLog(event.id().event(), puppi_collection, puppi_collection_backend, nbx_map_collection, nbx_map_collection_backend, table_collection, table_collection_backend, bx_ct);
      }
    }

    void beginStream(edm::StreamID) override {
      if (verbose_) {
        fmt::print("============================================================\n");
        start_ = std::chrono::steady_clock::now();
      }
    }

    void endStream() override {
      if (verbose_) {
        end_ = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_ - start_);
        fmt::print("[DEBUG] l1sc::W3PiAnalysisResults: {} -> {} ({} ms)\n", nbx_processed_, nbx_selected_, duration.count());
        fmt::print("============================================================\n");
      }
    }

  private:
    const bool verbose_;
    const int verbose_level_; 
    const edm::EDGetTokenT<PuppiHostCollection> puppi_token_;
    const edm::EDGetTokenT<unsigned short> puppi_backend_;
    const edm::EDGetTokenT<NbxMapHostCollection> nbx_map_token_;
    const edm::EDGetTokenT<unsigned short> nbx_map_backend_;
    const edm::EDGetTokenT<W3PiPuppiTableHostCollection> table_token_;
    const edm::EDGetTokenT<unsigned short> table_backend_;  
    const edm::EDGetTokenT<CounterHost> bx_ct_token_;
    std::chrono::steady_clock::time_point start_, end_;
    uint32_t nbx_selected_ = 0;
    uint32_t nbx_processed_ = 0;

    const int32_t kMaxView = 5;

    void debugLog(int event_id, 
                              const PuppiHostCollection &puppi_host, 
                              const cms::alpakatools::Backend puppi_backend, 
                              const NbxMapHostCollection &nbx_map_host,
                              const cms::alpakatools::Backend nbx_map_backend,
                              const W3PiPuppiTableHostCollection &table_host,
                              const cms::alpakatools::Backend table_backend,
                              const CounterHost &bx_ct) {
      const auto size = table_host.const_view().metadata().size();
      const auto size_puppi = puppi_host.const_view().metadata().size();

      if (size == 0) return;
      fmt::print("[DEBUG] PuppiCollection[{}] (event: {}, backend: {}):\n", size_puppi, event_id, cms::alpakatools::toString(puppi_backend));

      // table header
      const std::string pseparator = "+-------+-------+-------+---------+---------+---------+---------+-------+";
      fmt::print("{}\n", pseparator);
      fmt::print("| {:>5} | {:>5} | {:>5} | {:>7} | {:>7} | {:>7} | {:>7} | {:>5} |\n",
                "bx", "range", "local", "idx", "pt", "eta", "phi", "pdgid");
      fmt::print("{}\n", pseparator);

      const auto &offsets = nbx_map_host.const_view<OffsetsSoA>().offsets();
      const auto &bx_array = nbx_map_host.const_view<NbxSoA>().bx();

      int max_entries = (verbose_level_ == 1) ? size_puppi : kMaxView; // suppress tail unless verbose

      int printed = 0;
      for (int32_t i = 0; i < nbx_map_host.const_view().metadata().size(); ++i) {
          int bx = bx_array[i];
          int start = offsets[i];
          int end = offsets[i + 1];
          int range = end - start;

          for (int j = 0; j < range; ++j) {
              if (printed >= max_entries) break;
              int localidx = j;
              int idx = start + j;
              if (idx >= size_puppi) break;

              const auto &view = puppi_host.const_view()[idx];
              fmt::print("| {:5d} | {:5d} | {:5d} | {:7d} | {:7.2f} | {:7.2f} | {:7.2f} | {:5d} |\n",
                        bx,
                        range,
                        localidx,
                        idx,
                        view.pt(),
                        view.eta(),
                        view.phi(),
                        view.pdgid());
              ++printed;
          }
          if (printed >= max_entries) break;
      }

      // log tail if suppressed
      if (printed < size_puppi) {
          fmt::print("| {:>5} | {:>5} | {:>5} | {:>7} | {:>7} | {:>7} | {:>7} | {:>5} |\n",
                    "...", "...", "...", "...", "...", "...", "...", "...");
      }

      printed = 0;
      for (int32_t i = nbx_map_host.const_view().metadata().size(); i > 0; i--) {
          int bx = bx_array[i];
          int start = offsets[i];
          int end = offsets[i+1];
          int range = end - start;

          for (int j = start; j < end; ++j) {
              if (printed >= max_entries) break;
              int localidx = j - start;
              int idx = j;
              if (idx >= size_puppi) break;

              const auto &view = puppi_host.const_view()[idx];
              fmt::print("| {:5d} | {:5d} | {:5d} | {:7d} | {:7.2f} | {:7.2f} | {:7.2f} | {:5d} |\n",
                        bx,
                        range,
                        localidx,
                        idx,
                        view.pt(),
                        view.eta(),
                        view.phi(),
                        view.pdgid());
              ++printed;
          }
          if (printed >= max_entries) break;
      }

      fmt::print("{}\n", pseparator);

      fmt::print("[DEBUG] TableCollection[{}] (event: {}, backend: {}):\n", size, event_id, cms::alpakatools::toString(puppi_backend));

      // table header
      const std::string separator = "+---------+---------+---------+---------+-------+";
      fmt::print("{}\n", separator);
      fmt::print("| {:>7} | {:>7} | {:>7} | {:>7} | {:>5} |\n",
                "idx", "pt", "eta", "phi", "bx");
      fmt::print("{}\n", separator);
      std::vector<int> bxs;
      for (int32_t idx = 0; idx < size; ++idx) {
        const auto &triplet = table_host.const_view()[idx];
        auto bx = 0;
        for (int32_t idx = 0; idx < nbx_map_host.view<NbxSoA>().metadata().size(); idx++) {
          if (triplet.i() < nbx_map_host.view<OffsetsSoA>().offsets()[idx+1] && triplet.i() >= nbx_map_host.view<OffsetsSoA>().offsets()[idx]) {
            bx = idx;
            break;
          }
        }
        bxs.push_back(bx);
        fmt::print("| {:7d} | {:7.2f} | {:7.2f} | {:7.2f} | {:>5} |\n",
                  static_cast<int>(triplet.i()),
                  puppi_host.const_view().pt()[triplet.i()],
                  puppi_host.const_view().eta()[triplet.i()],
                  puppi_host.const_view().phi()[triplet.i()],
                  bx);
        fmt::print("| {:7d} | {:7.2f} | {:7.2f} | {:7.2f} | {:>5} |\n",
                  static_cast<int>(triplet.j()),
                  puppi_host.const_view().pt()[triplet.j()],
                  puppi_host.const_view().eta()[triplet.j()],
                  puppi_host.const_view().phi()[triplet.j()],
                  bx);
        fmt::print("| {:7d} | {:7.2f} | {:7.2f} | {:7.2f} | {:>5} |\n",
                  static_cast<int>(triplet.k()),
                  puppi_host.const_view().pt()[triplet.k()],
                  puppi_host.const_view().eta()[triplet.k()],
                  puppi_host.const_view().phi()[triplet.k()],
                  bx);
        fmt::print("{}\n", separator);
      }
      nbx_selected_ += bxs.size();
      fmt::print("[DEBUG] Event {} -> BXs[{} out of {}]: {}\n", event_id, bxs.size(), bx_ct.value(), fmt::join(bxs, ", "));
      fmt::print("============================================================\n");
    }
  };

}  // namespace l1sc

DEFINE_FWK_MODULE(l1sc::L1TScPhase2W3PiAnalyzer);