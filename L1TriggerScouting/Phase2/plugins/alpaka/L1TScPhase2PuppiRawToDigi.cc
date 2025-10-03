#include "DataFormats/L1ScoutingRawData/interface/SDSRawDataCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/CounterHost.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/BxLookupDeviceCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/PuppiDeviceCollection.h"
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
#include "L1TriggerScouting/Phase2/interface/alpaka/SynchronizingTimer.h"
#include "L1TriggerScouting/Phase2/interface/alpaka/L1TScPhase2PuppiRawToDigiKernels.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc {

  struct BxData {
    unsigned int bx;
    const data_t *header_ptr;
    const data_t *data_ptr;
    size_t data_size;

    // comparison operator for priority queue
    bool operator>(const BxData &other) const { return bx > other.bx; }
  };

  using MinHeap = std::priority_queue<BxData, std::vector<BxData>, std::greater<>>;

  using namespace ::l1sc;

  class L1TScPhase2PuppiRawToDigi : public stream::EDProducer<> {
  public:
    L1TScPhase2PuppiRawToDigi(const edm::ParameterSet &params)
        : EDProducer<>(params),
          raw_data_token_{consumes(params.getParameter<edm::InputTag>("src"))},
          puppi_token_{produces()},
          bx_lookup_token_{produces()},
          nbx_token_{produces("nbx")},
          streams_(params.getParameter<std::vector<uint32_t>>("streams")),
          splitFactor_(params.getParameter<unsigned int>("splitFactor")),
          environment_{static_cast<Environment>(params.getUntrackedParameter<int>("environment"))},
          sync_timer_(std::in_place, "L1TScPhase2PuppiRawToDigi", environment_) {}

    void produce(device::Event &event, const device::EventSetup &event_setup) override {
      // debug / test
      sync_timer_.value().start(event.queue());

      // get raw data input
      const auto &raw_data = event.get(raw_data_token_);

      // normalize header & payload
      auto ngoodbx = normalize(raw_data);
      const auto nbx = static_cast<int32_t>(h_data_.size());

      // allocate memory buffers
      auto bx_lookup = BxLookupDeviceCollection({{nbx, nbx + 1}}, event.queue());
      auto puppi = PuppiDeviceCollection(p_data_.size(), event.queue());

      // initialize device constant memory (called once)
      rtd_kernels_.initialize(event.queue());

      // decode raw data
      kernels::decode(event.queue(), h_data_.data(), bx_lookup);
      kernels::decode(event.queue(), p_data_.data(), puppi);

      // store data in the event (device-side products)
      event.emplace(bx_lookup_token_, std::move(bx_lookup));
      event.emplace(puppi_token_, std::move(puppi));

      // store nbx
      auto nbx_portable = CounterHost(event.queue(), static_cast<unsigned int>(ngoodbx));
      event.emplace(nbx_token_, std::move(nbx_portable));

      // debug / test end
      sync_timer_.value().sync(event.queue());
    };

    static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<std::vector<uint32_t>>("streams");
      desc.add<unsigned int>("splitFactor", 1)->setComment("Number of streams per BX");
      desc.add<edm::InputTag>("src");
      desc.addUntracked<int>("environment", static_cast<int>(Environment::kProduction));
      descriptions.addWithDefaultLabel(desc);
    };

    unsigned int normalize(const SDSRawDataCollection &raw_data) {
      p_data_.clear();
      h_data_.clear();

      MinHeap min_heap;
      // readout data from links breadth-first (order not guaranteed)
      for (auto stream_id : streams_) {
        const auto &stream = raw_data.FEDData(stream_id);
        const auto chunk_begin = reinterpret_cast<const data_t *>(stream.data());
        const auto chunk_end = reinterpret_cast<const data_t *>(stream.data() + stream.size());

        for (auto ptr = chunk_begin; ptr < chunk_end;) {
          if (*ptr == 0) {
            ++ptr;
            continue;
          }  // skip empty words

          unsigned int bx = ((*ptr) >> 12) & 0xFFF;  // unpack bx number
          auto chunk_size = (*ptr) & 0xFFF;          // unpack chunk size
          ++ptr;                                     // move past header

          const size_t payload = chunk_end - ptr;                           // calculate payload size
          const size_t copy_count = std::min<size_t>(chunk_size, payload);  // block size

          min_heap.push({bx, ptr - 1, ptr, copy_count});
          ptr += copy_count;  // move to the next word
        }
      }

      // restore bx order (size of heap at max 3564 not that heavy to track)
      // pop the first element outside the loop
      if (min_heap.empty())
        return 0;

      unsigned int ngoodbx = 0;
      BxData bx_data = min_heap.top();
      h_data_.push_back(*(bx_data.header_ptr));                                               // store header
      p_data_.insert(p_data_.end(), bx_data.data_ptr, bx_data.data_ptr + bx_data.data_size);  // copy payload
      min_heap.pop();
      unsigned int nslices = 1;
      while (!min_heap.empty()) {
        const auto &bx_data2 = min_heap.top();
        if (bx_data2.bx != bx_data.bx) {  // new bx
          if (nslices == splitFactor_)
            ngoodbx++;
          nslices = 1;
          h_data_.push_back(*(bx_data2.header_ptr));           // store header
        } else {                                               // same bx, skip header
          h_data_.back() += (*(bx_data2.header_ptr)) & 0xFFF;  // merge size in header
          nslices++;
        }
        bx_data = bx_data2;
        p_data_.insert(p_data_.end(), bx_data.data_ptr, bx_data.data_ptr + bx_data.data_size);  // copy payload
        min_heap.pop();
      }
      if (nslices == splitFactor_)
        ngoodbx++;

      return ngoodbx;
    }

  private:
    // consume host side input data
    const edm::EDGetTokenT<SDSRawDataCollection> raw_data_token_;

    // produce device-side products
    const device::EDPutToken<PuppiDeviceCollection> puppi_token_;
    const device::EDPutToken<BxLookupDeviceCollection> bx_lookup_token_;

    // produce host-side products
    const edm::EDPutTokenT<CounterHost> nbx_token_;

    // utility members
    const std::vector<uint32_t> streams_;
    const unsigned int splitFactor_;  // number of streams per BX
    const Environment environment_;

    // temporary storage
    std::vector<uint64_t> h_data_;
    std::vector<uint64_t> p_data_;

    // kernel
    kernels::L1TScPhase2PuppiRawToDigiKernels rtd_kernels_;

    // debug / test stats
    std::optional<SynchronizingTimer> sync_timer_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc

DEFINE_FWK_ALPAKA_MODULE(l1sc::L1TScPhase2PuppiRawToDigi);