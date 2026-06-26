#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/FEDRawData/interface/FEDRawData.h"
#include "DataFormats/L1ScoutingRawData/interface/SDSNumbering.h"
#include "DataFormats/L1ScoutingRawData/interface/SDSRawDataCollection.h"
#include "DataFormats/L1Scouting/interface/OrbitCollection.h"
#include "DataFormats/L1Trigger/interface/VertexWord.h"

class ScPhase2VertexPacker : public edm::stream::EDProducer<> {
public:
  explicit ScPhase2VertexPacker(const edm::ParameterSet &);
  ~ScPhase2VertexPacker() override;
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  void produce(edm::Event &, const edm::EventSetup &) override;
  uint64_t packVertex(l1t::VertexWord);

  edm::EDGetTokenT<std::vector<l1t::VertexWord>> src_;
  std::vector<unsigned int> fedIDs_;
  unsigned int splitFactor_;
  bool scoutingHeader_;
};

ScPhase2VertexPacker::ScPhase2VertexPacker(const edm::ParameterSet &iConfig)
    : src_(consumes<std::vector<l1t::VertexWord>>(iConfig.getParameter<edm::InputTag>("src"))),
      fedIDs_(iConfig.getParameter<std::vector<unsigned int>>("fedIDs")),
      splitFactor_(iConfig.getParameter<unsigned int>("splitFactor")),
      scoutingHeader_(iConfig.getParameter<bool>("scoutingHeader")) {
  produces<SDSRawDataCollection>();
  if (splitFactor_ != fedIDs_.size()) {
    throw cms::Exception("Configuration") << "splitFactor must match number of fedIDs";
  }
}

ScPhase2VertexPacker::~ScPhase2VertexPacker() {};

void ScPhase2VertexPacker::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  edm::Handle<std::vector<l1t::VertexWord>> vertices;
  iEvent.getByToken(src_, vertices);
  // Packing logic would go here, creating an SDSRawDataCollection
  auto outputCollection = std::make_unique<SDSRawDataCollection>();
  unsigned int ncands = vertices->size();
  unsigned int candsPerFragment = (ncands + splitFactor_ - 1) / splitFactor_;
  for (unsigned int i = 0; i < splitFactor_; ++i) {
    unsigned int startIdx = i * candsPerFragment;
    unsigned int endIdx = std::min(startIdx + candsPerFragment, ncands);
    unsigned int nCandsInFragment = endIdx - startIdx;
    auto &fedData = outputCollection->FEDData(fedIDs_[i]);
    if (startIdx >= endIdx) {
      // No candidates for this fragment, just add header if needed
      if (scoutingHeader_) {
        uint64_t header = (0b10llu << 62) | (uint64_t(iEvent.id().event()) << 24);
        fedData.resize(sizeof(header));
        *reinterpret_cast<uint64_t *>(fedData.data()) = header;
      }
      continue;
    } else {
      // Allocate space for candidates and optional header
      auto nwords = (nCandsInFragment + (scoutingHeader_ ? 1 : 0));
      fedData.resize(nwords * sizeof(uint64_t));
      uint64_t *output = reinterpret_cast<uint64_t *>(fedData.data());
      // Add header if needed
      if (scoutingHeader_) {
        uint64_t header = (0b10llu << 62) | (uint64_t(iEvent.id().event()) << 24) | (nCandsInFragment & 0xFFF);
        *output++ = header;
      }
      // Pack candidates from startIdx to endIdx into a fragment
      for (unsigned int j = startIdx; j < endIdx; ++j) {
        *output++ = packVertex((*vertices)[j]);
      }
    }
  }
  iEvent.put(std::move(outputCollection));
}

uint64_t ScPhase2VertexPacker::packVertex(l1t::VertexWord vertex) {
  return vertex.validBits() | ((uint64_t)vertex.z0Bits() << l1t::VertexWord::VertexBitLocations::kZ0LSB) |
         ((uint64_t)vertex.multiplicityBits() << l1t::VertexWord::VertexBitLocations::kNTrackInPVLSB) |
         ((uint64_t)vertex.ptBits() << l1t::VertexWord::VertexBitLocations::kSumPtLSB) |
         ((uint64_t)vertex.qualityBits() << l1t::VertexWord::VertexBitLocations::kQualityLSB) |
         ((uint64_t)vertex.inverseMultiplicityBits() << l1t::VertexWord::VertexBitLocations::kNTrackOutPVLSB) |
         ((uint64_t)vertex.unassignedBits() << l1t::VertexWord::VertexBitLocations::kUnassignedLSB);
}

void ScPhase2VertexPacker::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("rawDataCollector"));
  desc.add<std::vector<unsigned int>>("fedIDs")->setComment("List of FED IDs to use");
  desc.add<unsigned int>("splitFactor", 1)->setComment("Number of fragments per BX");
  desc.add<bool>("scoutingHeader", false)->setComment("Include scouting header in output");
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScPhase2VertexPacker);
