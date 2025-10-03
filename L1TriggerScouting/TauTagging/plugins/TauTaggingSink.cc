#include <cassert>

#include <fmt/format.h>

#include "DataFormats/L1ScoutingSoA/interface/BxLookupHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/ClustersHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/PFCandidateHostCollection.h"
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

  class TauTaggingSink : public edm::stream::EDAnalyzer<> {
  public:
    TauTaggingSink(const edm::ParameterSet& params)
        : pf_token_{consumes(params.getUntrackedParameter<edm::InputTag>("src"))},
          bx_lookup_token_{consumes(params.getUntrackedParameter<edm::InputTag>("src"))},
          clusters_token_{consumes(params.getUntrackedParameter<edm::InputTag>("clusters"))},
          pf_backend_{consumes(getBackendTag(params.getUntrackedParameter<edm::InputTag>("src")))},
          bx_lookup_backend_{consumes(getBackendTag(params.getUntrackedParameter<edm::InputTag>("src")))},
          clusters_backend_{consumes(getBackendTag(params.getUntrackedParameter<edm::InputTag>("clusters")))},
          environment_{static_cast<Environment>(params.getUntrackedParameter<int>("environment"))},
          run_scout_{params.getParameter<bool>("run_scout")} {}

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<bool>("run_scout");
      desc.addUntracked<int>("environment", static_cast<int>(Environment::kDevelopment));
      desc.addUntracked<edm::InputTag>("src");
      desc.addUntracked<edm::InputTag>("clusters");
      descriptions.addWithDefaultLabel(desc);
    }

    void analyze(edm::Event const& event, edm::EventSetup const&) override {
      if (environment_ >= Environment::kDevelopment) {
        constexpr int total_len = 100;
        auto label = fmt::format("EVENT: {}", event.id().event());
        int pad = total_len - static_cast<int>(label.size());
        int pad_left = pad / 2;
        int pad_right = pad - pad_left - 1;
        fmt::print("\n{0} {1} {2}\n\n", std::string(pad_left, '-'), label, std::string(pad_right, '-'));

        // unpacking/decoding step debugging
        const auto pf_handle = event.getHandle(pf_token_);
        const auto bx_lookup_handle = event.getHandle(bx_lookup_token_);
        const auto clusters_handle = event.getHandle(clusters_token_);

        if (pf_handle.isValid()) {
          // pf
          auto const& pf = *pf_handle;
          auto const pf_backend = static_cast<Backend>(event.get(pf_backend_));

          if (clusters_handle.isValid()) {
            // clusters
            auto const& clusters = *clusters_handle;
            auto const clusters_backend = static_cast<Backend>(event.get(clusters_backend_));
            if (run_scout_ && bx_lookup_handle.isValid()) {
              // bx lookup
              auto const& bx_lookup = *bx_lookup_handle;
              auto const bx_lookup_backend = static_cast<Backend>(event.get(bx_lookup_backend_));
              assert(pf_backend == bx_lookup_backend);
              assert(pf_backend == clusters_backend);
              print(pf.const_view(),
                    bx_lookup.const_view<BxIndexSoA>(),
                    bx_lookup.const_view<OffsetsSoA>(),
                    clusters.const_view(),
                    toString(bx_lookup_backend));
            } else {
              assert(pf_backend == clusters_backend);
              print(pf.const_view(), clusters.const_view(), toString(clusters_backend));
            }
          } else {
            // debug unpacker only
            print(pf.const_view(), toString(pf_backend));
          }
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
      fmt::print("[INFO] OK - TauTagging\n");
      fmt::print("=========================================================================\n");
    }

    void print(const PFCandidateHostCollection::ConstView& pf,
               const ClustersHostCollection::ConstView& clusters,
               const std::string_view clusters_backend) {
      const auto size = pf.metadata().size();
      if (size == 0)
        return;

      // Header
      int clusters_num = 0;
      for (int i = 0; i < clusters.metadata().size(); ++i) {
        if (clusters.cluster()[i] > clusters_num)
          clusters_num = clusters.cluster()[i];
      }
      fmt::print("[DEBUG] CLUETaus[{}] ({}) found {} clusters:\n", size, clusters_backend, clusters_num);

      constexpr auto sep =
          "+---------+---------+---------+---------+---------+---------+---------+---------+---------+---------+-------"
          "--+---------+---------+";
      auto printHeader = [&] {
        fmt::print("{}\n", sep);
        fmt::print(
            "| {:>7} | {:>7} | {:>7} | {:>7} | {:>7} | {:>7} | {:>7} | {:>7} | {:>7} | {:>7} | {:>7} | {:>7} | {:>7} "
            "|\n",
            "index",
            "is_seed",
            "cluster",
            "pt",
            "eta",
            "phi",
            "mass",
            "z0",
            "dxy",
            "puppiw",
            "charge",
            "type",
            "pdgid");
        fmt::print("{}\n", sep);
      };

      auto printRow = [&](int index, const auto& pf_view, const auto& clusters_view) {
        fmt::print(
            "| {:>7d} | {:>7d} | {:>7d} | {:>7.2f} | {:>7.2f} | {:>7.2f} | {:>7.2f} | {:>7.2f} | {:>7.2f} | {:>7.2f} | "
            "{:>7d} | {:>7d} | {:>7d} |\n",
            index,
            clusters_view.is_seed(),
            clusters_view.cluster(),
            pf_view.pt(),
            pf_view.eta(),
            pf_view.phi(),
            pf_view.mass(),
            pf_view.z0(),
            pf_view.dxy(),
            pf_view.puppiw(),
            pf_view.charge(),
            pf_view.type(),
            pf_view.pdgid());
      };

      printHeader();

      const int max_entries = (environment_ > Environment::kTest) ? size : 5;

      for (int i = 0; i < pf.metadata().size() && i < max_entries; ++i) {
        if (i >= size)
          break;
        printRow(i, pf[i], clusters[i]);
      }

      if (max_entries < size) {
        fmt::print(
            "| {:>7} | {:>7} | {:>7} | {:>7} | {:>7} | {:>7} | {:>7} | {:>7} | {:>7} | {:>7} | {:>7} | {:>7} | {:>7} "
            "|\n",
            "...",
            "...",
            "...",
            "...",
            "...",
            "...",
            "...",
            "...",
            "...",
            "...",
            "...",
            "...",
            "...");
        for (int j = pf.metadata().size() - 5; j < pf.metadata().size(); ++j) {
          printRow(j, pf[j], clusters[j]);
        }
      }

      fmt::print("{}\n", sep);
    }

    void print(const PFCandidateHostCollection::ConstView& pf,
               const BxIndexSoA::ConstView& bx_index,
               const OffsetsSoA::ConstView& offsets,
               const ClustersHostCollection::ConstView& clusters,
               const std::string_view pf_backend) {
      const auto size = pf.metadata().size();
      if (size == 0)
        return;

      int all_clusters_num = 0;
      for (int i = 0; i < bx_index.metadata().size(); ++i) {
        auto bx_idx = bx_index.bx()[i];
        auto begin = offsets.offsets()[bx_idx];
        auto end = offsets.offsets()[bx_idx + 1];
        int clusters_num = 0;
        for (uint32_t j = begin; j < end; ++j)
          if (clusters.cluster()[i] > clusters_num)
            clusters_num = clusters.cluster()[i];
        all_clusters_num += clusters_num;
      }
      // Header
      fmt::print("[DEBUG] PFCandidateCollection[{}] ({}) found {} clusters per BX (avg. across {} BXs):\n",
                 size,
                 pf_backend,
                 static_cast<int>(all_clusters_num / bx_index.metadata().size()),
                 bx_index.metadata().size());

      constexpr auto sep =
          "+-------+---------+-------+---------+---------+---------+---------+---------+---------+---------+---------+-"
          "--------+---------+---------+---------+";
      auto printHeader = [&] {
        fmt::print("{}\n", sep);
        fmt::print(
            "| {:>5} | {:>7} | {:>5} | {:>7} | {:>7} | {:>7} | {:>7} | {:>7} | {:>7} | {:>7} | {:>7} | {:>7} | {:>7} | "
            "{:>7} | {:>7} |\n",
            "bx",
            "index",
            "local",
            "is_seed",
            "cluster",
            "pt",
            "eta",
            "phi",
            "mass",
            "z0",
            "dxy",
            "puppiw",
            "charge",
            "type",
            "pdgid");
        fmt::print("{}\n", sep);
      };

      auto printRow = [&](int bx, int global, int local, const auto& clusters_view, const auto& pf_view) {
        fmt::print(
            "| {:5d} | {:7d} | {:5d} | {:7d} | {:7d} | {:>7.2f} | {:>7.2f} | {:>7.2f} | {:>7.2f} | {:>7.2f} | {:>7.2f} "
            "| {:>7.2f} | {:>7d} | {:>7d} | {:>7d} |\n",
            bx,
            global,
            local,
            clusters_view.is_seed(),
            clusters_view.cluster(),
            pf_view.pt(),
            pf_view.eta(),
            pf_view.phi(),
            pf_view.mass(),
            pf_view.z0(),
            pf_view.dxy(),
            pf_view.puppiw(),
            pf_view.charge(),
            pf_view.type(),
            pf_view.pdgid());
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
          printRow(bx, globalIdx, j, clusters[globalIdx], pf[globalIdx]);
          ++printed;
        }
      }

      // Ellipsis row if not all printed
      if (printed < size) {
        fmt::print(
            "| {:>5} | {:>7} | {:>5} | {:>7} | {:>7} | {:>7} | {:>7} | {:>7} | {:>7} | {:>7} | {:>7} | {:>7} | {:>7} | "
            "{:>7} | {:>7} |\n",
            "...",
            "...",
            "...",
            "...",
            "...",
            "...",
            "...",
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
          printRow(bx_index.bx()[lastBxIdx], globalIdx, j, clusters[globalIdx], pf[globalIdx]);
        }
      }

      fmt::print("{}\n", sep);
    }

    void print(const PFCandidateHostCollection::ConstView& pf, const std::string_view pf_backend) {
      const auto size = pf.metadata().size();
      if (size == 0)
        return;

      // Header
      fmt::print("[DEBUG] PFCandidateCollection[{}] ({}):\n", size, pf_backend);

      constexpr auto sep =
          "+---------+---------+---------+---------+---------+---------+---------+---------+---------+---------+-------"
          "--+";
      auto printHeader = [&] {
        fmt::print("{}\n", sep);
        fmt::print("| {:>7} | {:>7} | {:>7} | {:>7} | {:>7} | {:>7} | {:>7} | {:>7} | {:>7} | {:>7} | {:>7} |\n",
                   "index",
                   "pt",
                   "eta",
                   "phi",
                   "mass",
                   "z0",
                   "dxy",
                   "puppiw",
                   "charge",
                   "type",
                   "pdgid");
        fmt::print("{}\n", sep);
      };

      auto printRow = [&](int index, const auto& pf_view) {
        fmt::print(
            "| {:7d} | {:7.2f} | {:7.2f} | {:7.2f} | {:7.2f} | {:7.2f} | {:7.2f} | {:7.2f} | {:7d} | {:7d} | {:7d} |\n",
            index,
            pf_view.pt(),
            pf_view.eta(),
            pf_view.phi(),
            pf_view.mass(),
            pf_view.z0(),
            pf_view.dxy(),
            pf_view.puppiw(),
            pf_view.charge(),
            pf_view.type(),
            pf_view.pdgid());
      };

      printHeader();

      const int max_entries = (environment_ > Environment::kTest) ? size : 5;

      for (int i = 0; i < pf.metadata().size() && i < max_entries; ++i) {
        if (i >= size)
          break;
        printRow(i, pf[i]);
      }

      if (max_entries < size) {
        fmt::print("| {:>7} | {:>7} | {:>7} | {:>7} | {:>7} | {:>7} | {:>7} | {:>7} | {:>7} | {:>7} |\n",
                   "...",
                   "...",
                   "...",
                   "...",
                   "...",
                   "...",
                   "...",
                   "...",
                   "...",
                   "...",
                   "...");
        for (int j = pf.metadata().size() - 5; j < pf.metadata().size(); ++j) {
          printRow(j, pf[j]);
        }
      }

      fmt::print("{}\n", sep);
    }

  private:
    // get products
    const edm::EDGetTokenT<PFCandidateHostCollection> pf_token_;
    const edm::EDGetTokenT<BxLookupHostCollection> bx_lookup_token_;
    const edm::EDGetTokenT<ClustersHostCollection> clusters_token_;
    // backend query
    const edm::EDGetTokenT<unsigned short> pf_backend_;
    const edm::EDGetTokenT<unsigned short> bx_lookup_backend_;
    const edm::EDGetTokenT<unsigned short> clusters_backend_;
    // debug
    const Environment environment_;
    const bool run_scout_;
  };

}  // namespace l1sc

DEFINE_FWK_MODULE(l1sc::TauTaggingSink);