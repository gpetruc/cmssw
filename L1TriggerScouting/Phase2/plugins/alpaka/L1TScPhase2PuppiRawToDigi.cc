#include "DataFormats/L1ScoutingSoA/interface/alpaka/OrbitEventIndexMapDeviceCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/PuppiDeviceCollection.h"
#include "DataFormats/L1ScoutingRawData/interface/SDSRawDataCollection.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/Event.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EventSetup.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/EDProducer.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "L1TriggerScouting/Phase2/interface/L1TScPhase2Common.h"
#include "L1TriggerScouting/Phase2/interface/alpaka/L1TScPhase2PuppiRawToDigiKernels.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc {

  using namespace ::l1sc;

  /**
   * @class L1TScPhase2PuppiRawToDigi
   * @brief Produces PuppiDeviceCollection (PortableCollection)
   */
  class L1TScPhase2PuppiRawToDigi : public stream::EDProducer<> {
  public:
    L1TScPhase2PuppiRawToDigi(const edm::ParameterSet &params);

    void produce(device::Event &event, const device::EventSetup &event_setup) override;
    static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

  private:
    const edm::EDGetTokenT<SDSRawDataCollection> raw_data_token_;                               // raw data
    const device::EDPutToken<PuppiDeviceCollection> puppi_token_;                               // PUPPI candidates
    const device::EDPutToken<OrbitEventIndexMapDeviceCollection> orbit_association_map_token_;  // orbit association map
    const std::vector<uint32_t> links_ids_;                                           // front-end devices stream links
    std::array<data_t, kOrbitSize> h_data_{};                                         // headers 64-bit words
    std::vector<data_t> p_data_{};                                                    // payload 64-bit words
    std::unique_ptr<kernels::L1TScPhase2PuppiRawToDigiKernels> raw_to_digi_kernels_;  // kernels for decoding
    const bool verbose_;                                                              // verbose output
    const int verbose_level_;                                                         // verbose level

    void collectBuffers(const SDSRawDataCollection &raw_data);
    void logDebugMessage(int event_id, const PuppiHostCollection &puppi_host) const;
  };

  // __________________________________________________________________________________________________________________
  // IMPLEMENTATION

  L1TScPhase2PuppiRawToDigi::L1TScPhase2PuppiRawToDigi(const edm::ParameterSet &params)
      : EDProducer<>(params),
        raw_data_token_{consumes(params.getParameter<edm::InputTag>("src"))},
        puppi_token_{produces()},
        orbit_association_map_token_{produces()},
        links_ids_(params.getParameter<std::vector<uint32_t>>("linksIds")),
        raw_to_digi_kernels_(std::make_unique<kernels::L1TScPhase2PuppiRawToDigiKernels>()),
        verbose_(params.getUntrackedParameter<bool>("verbose")),
        verbose_level_(params.getUntrackedParameter<int>("verboseLevel")) {}

  void L1TScPhase2PuppiRawToDigi::produce(device::Event &event, const device::EventSetup &event_setup) {
    // intialize device constant memory -> called only once
    raw_to_digi_kernels_->initialize(event.queue());

    // get raw data input
    auto raw_data = event.getHandle(raw_data_token_);

    // preprocess header -> payload
    collectBuffers(*raw_data);

    // orbit event index association map
    auto map_size = links_ids_.size() * kOrbitSize + 1;
    auto orbit_association_map = OrbitEventIndexMapDeviceCollection(map_size, event.queue());
    kernels::associateOrbitEventIndex(event.queue(), h_data_.data(), orbit_association_map);

    // pf candidates data
    auto puppi = PuppiDeviceCollection(p_data_.size(), event.queue());
    kernels::rawToDigi(event.queue(), p_data_.data(), puppi);

    // debug log to stdout
    if (verbose_) {
      auto puppi_host = PuppiHostCollection(puppi.view().metadata().size(), event.queue());
      alpaka::memcpy(event.queue(), puppi_host.buffer(), puppi.buffer());
      alpaka::wait(event.queue());  // explicitly synchronize the device
      logDebugMessage(event.id().event(), puppi_host);
    }

    // store data in the event
    event.emplace(orbit_association_map_token_, std::move(orbit_association_map));
    event.emplace(puppi_token_, std::move(puppi));
  }

  void L1TScPhase2PuppiRawToDigi::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<std::vector<uint32_t>>("linksIds");
    desc.add<edm::InputTag>("src");
    desc.addUntracked<bool>("verbose", false);
    desc.addUntracked<int>("verboseLevel", 0);
    descriptions.addWithDefaultLabel(desc);
  }

  void L1TScPhase2PuppiRawToDigi::collectBuffers(const SDSRawDataCollection &raw_data) {
    p_data_.clear();  // reset payload buffer

    size_t h_idx = 0;
    for (auto &link_id : links_ids_) {
      const auto &link = raw_data.FEDData(link_id);
      const auto chunk_begin = reinterpret_cast<const data_t *>(link.data());
      const auto chunk_end = reinterpret_cast<const data_t *>(link.data() + link.size());

      for (auto ptr = chunk_begin; ptr < chunk_end;) {
        // skip empty words
        if (*ptr == 0) {
          ++ptr;
          continue;
        }

        h_data_[h_idx++] = *ptr;           // store header
        auto chunk_size = (*ptr) & 0xFFF;  // unpack chunk size
        ++ptr;                             // move to the next word

        const size_t payload = chunk_end - ptr;                           // calculate payload size
        const size_t copy_count = std::min<size_t>(chunk_size, payload);  // block size

        // skip if no trailing payload
        if (copy_count == 0)
          continue;

        p_data_.insert(p_data_.end(), ptr, ptr + copy_count);  // copy payload
        ptr += copy_count;                                     // move to the next word
      }
    }
  }

  /**
   * Log converstion results to stdout
   */
  void L1TScPhase2PuppiRawToDigi::logDebugMessage(int event_id, const PuppiHostCollection &puppi_host) const {
    const auto size = puppi_host.const_view().metadata().size();
    fmt::print("[DEBUG] l1sc::L1TScPhase2PuppiRawToDigi: PuppiDeviceCollection[{}]:\n", size);

    // table header
    const std::string separator = "+-------+---------+---------+---------+---------+-------+";
    fmt::print("{}\n", separator);
    fmt::print("| {:>5} | {:>7} | {:>7} | {:>7} | {:>7} | {:>5} |\n", "index", "pt", "eta", "phi", "z0", "pdgid");
    fmt::print("{}\n", separator);

    // log head of collection
    auto span = (size > 10) ? 10 : size;
    if (verbose_level_ == 1)
      span = size;
    for (int32_t idx = 0; idx < span; ++idx) {
      const auto &view = puppi_host.const_view()[idx];
      fmt::print("| {:5d} | {:7.2f} | {:7.2f} | {:7.2f} | {:7.2f} | {:5d} |\n",
                 idx,
                 view.pt(),
                 view.eta(),
                 view.phi(),
                 view.z0(),
                 view.pdgid());
    }

    // log tail
    if (span < size) {
      fmt::print("| {:>5} | {:>7} | {:>7} | {:>7} | {:>7} | {:>5} |\n", "...", "...", "...", "...", "...", "...");
    }

    fmt::print("{}\n", separator);
    fmt::print("[DEBUG] l1sc::L1TScPhase2PuppiRawToDigi: OK ({})\n", event_id);
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc

DEFINE_FWK_ALPAKA_MODULE(l1sc::L1TScPhase2PuppiRawToDigi);