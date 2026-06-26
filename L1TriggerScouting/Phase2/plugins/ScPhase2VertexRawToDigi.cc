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

#include <ap_int.h>

class ScPhase2VertexRawToDigi : public edm::stream::EDProducer<> {
public:
  explicit ScPhase2VertexRawToDigi(const edm::ParameterSet &);
  ~ScPhase2VertexRawToDigi() override;
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  void produce(edm::Event &, const edm::EventSetup &) override;

  template <typename T>
  std::unique_ptr<OrbitCollection<T>> unpackObj(unsigned int orbit,
                                                const SDSRawDataCollection &feds,
                                                std::vector<std::vector<T>> &buffer);

  edm::EDGetTokenT<SDSRawDataCollection> rawToken_;
  std::vector<unsigned int> fedIDs_;
  uint8_t splitFactor_;  // number of fragments per BX

  // temporary storage
  std::vector<std::vector<l1t::VertexWord>> candBuffer_;
  unsigned int nbx_;

  void unpackFromRaw(uint64_t data, std::vector<l1t::VertexWord> &outBuffer);
};

ScPhase2VertexRawToDigi::ScPhase2VertexRawToDigi(const edm::ParameterSet &iConfig)
    : rawToken_(consumes<SDSRawDataCollection>(iConfig.getParameter<edm::InputTag>("src"))),
      fedIDs_(iConfig.getParameter<std::vector<unsigned int>>("fedIDs")),
      splitFactor_(iConfig.getParameter<unsigned int>("splitFactor")) {
  candBuffer_.resize(OrbitCollection<l1t::VertexWord>::NBX + 1);
  produces<OrbitCollection<l1t::VertexWord>>();
  produces<unsigned int>("nbx");
}

ScPhase2VertexRawToDigi::~ScPhase2VertexRawToDigi() {};

void ScPhase2VertexRawToDigi::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  edm::Handle<SDSRawDataCollection> scoutingRawDataCollection;
  iEvent.getByToken(rawToken_, scoutingRawDataCollection);
  iEvent.put(unpackObj(iEvent.id().event(), *scoutingRawDataCollection, candBuffer_));
  iEvent.put(std::make_unique<unsigned int>(nbx_), "nbx");
}

template <typename T>
std::unique_ptr<OrbitCollection<T>> ScPhase2VertexRawToDigi::unpackObj(unsigned int orbit,
                                                                       const SDSRawDataCollection &feds,
                                                                       std::vector<std::vector<T>> &buffer) {
  unsigned int ntot = 0;
  nbx_ = 0;
  std::array<uint8_t, OrbitCollection<T>::NBX> bxcount;
  std::fill(bxcount.begin(), bxcount.end(), 0);
  for (auto &fedId : fedIDs_) {
    const FEDRawData &src = feds.FEDData(fedId);
    const uint64_t *begin = reinterpret_cast<const uint64_t *>(src.data());
    const uint64_t *end = reinterpret_cast<const uint64_t *>(src.data() + src.size());
    for (auto p = begin; p != end;) {
      if ((*p) == 0) {
        ++p;
        continue;
      }
      unsigned int bx = ((*p) >> 12) & 0xFFF;
      unsigned int orbitno = ((*p) >> 24) & 0xFFFFFFFFFlu;
      unsigned int nwords = (*p) & 0xFFF;
      if (orbitno != orbit) {
        throw cms::Exception("CorruptData") << "Data for orbit " << orbit << ", fedId " << fedId
                                            << " has header with mismatching orbit number " << orbitno << std::endl;
      }
      assert(bx < OrbitCollection<T>::NBX);
      auto nfound = ++bxcount[bx];
      if (nfound > splitFactor_) {
        throw cms::Exception("CorruptData") << "Data for orbit " << orbit << " has " << nfound << " blocks for bx "
                                            << bx << ", expected " << splitFactor_ << std::endl;
      } else if (nfound == splitFactor_) {
        nbx_++;
      }
      ++p;
      std::vector<T> &outputBuffer = buffer[bx + 1];
      outputBuffer.reserve(nwords);
      for (unsigned int i = 0; i < nwords; ++i, ++p) {
        uint64_t data = *p;
        unpackFromRaw(data, outputBuffer);
        ntot++;
      }
    }
  }
  return std::make_unique<OrbitCollection<T>>(buffer, ntot);
}

void ScPhase2VertexRawToDigi::unpackFromRaw(uint64_t data, std::vector<l1t::VertexWord> &outBuffer) {
  l1t::VertexWord::vtxvalid_t valid = (data & (((uint64_t)1 << l1t::VertexWord::VertexBitWidths::kValidSize) - 1));

  ap_int<l1t::VertexWord::VertexBitWidths::kZ0Size> z0_tot =
      (data >> l1t::VertexWord::VertexBitLocations::kZ0LSB) &
      ((((uint64_t)1 << l1t::VertexWord::VertexBitWidths::kZ0Size) - 1));
  l1t::VertexWord::vtxz0_t z0;
  z0.range() = z0_tot;

  l1t::VertexWord::vtxmultiplicity_t multIn =
      (data >> l1t::VertexWord::VertexBitLocations::kNTrackInPVLSB) &
      ((((uint64_t)1 << l1t::VertexWord::VertexBitWidths::kNTrackInPVSize) - 1));

  auto sumpt_tot = (data >> l1t::VertexWord::VertexBitLocations::kSumPtLSB) &
                   ((((uint64_t)1 << l1t::VertexWord::VertexBitWidths::kSumPtSize) - 1));
  l1t::VertexWord::vtxsumpt_t sumpt;
  sumpt.range() = sumpt_tot;

  l1t::VertexWord::vtxquality_t quality = (data >> l1t::VertexWord::VertexBitLocations::kQualityLSB) &
                                          ((((uint64_t)1 << l1t::VertexWord::VertexBitWidths::kQualitySize) - 1));

  l1t::VertexWord::vtxinversemult_t multOut =
      (data >> l1t::VertexWord::VertexBitLocations::kNTrackOutPVLSB) &
      ((((uint64_t)1 << l1t::VertexWord::VertexBitWidths::kNTrackOutPVSize) - 1));

  l1t::VertexWord::vtxunassigned_t unassigned =
      (data >> l1t::VertexWord::VertexBitLocations::kUnassignedLSB) &
      ((((uint64_t)1 << l1t::VertexWord::VertexBitWidths::kUnassignedSize) - 1));

  outBuffer.emplace_back(valid, z0, multIn, sumpt, quality, multOut, unassigned);
}

void ScPhase2VertexRawToDigi::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("rawDataCollector"));
  desc.add<std::vector<unsigned int>>("fedIDs");
  desc.add<unsigned int>("splitFactor", 1)->setComment("Number of fragments per BX");
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScPhase2VertexRawToDigi);
