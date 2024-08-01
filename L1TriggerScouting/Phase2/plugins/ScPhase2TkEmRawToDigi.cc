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
#include "DataFormats/L1TParticleFlow/interface/L1ScoutingTkEm.h"
#include "L1TriggerScouting/Phase2/interface/l1tkemUnpack.h"

class ScPhase2TkEmRawToDigi : public edm::stream::EDProducer<> {
public:
  explicit ScPhase2TkEmRawToDigi(const edm::ParameterSet &);
  ~ScPhase2TkEmRawToDigi() override;
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  //void beginStream(edm::StreamID) override;
  void produce(edm::Event &, const edm::EventSetup &) override;
  //void endStream() override;

  template <typename T>
  std::unique_ptr<OrbitCollection<T>> unpackObj(const SDSRawDataCollection &feds, std::vector<std::vector<T>> &buffer);

  edm::EDGetTokenT<SDSRawDataCollection> rawToken_;
  std::vector<unsigned int> fedIDs_;
  bool doStruct_;

  // temporary storage
  std::vector<std::vector<l1Scouting::TkEm>> structBuffer_;

  void unpackFromRaw(uint64_t data, std::vector<l1Scouting::TkEm> &outBuffer);
};

ScPhase2TkEmRawToDigi::ScPhase2TkEmRawToDigi(const edm::ParameterSet &iConfig)
    : rawToken_(consumes<SDSRawDataCollection>(iConfig.getParameter<edm::InputTag>("src"))),
      fedIDs_(iConfig.getParameter<std::vector<unsigned int>>("fedIDs")),
      doStruct_(iConfig.getParameter<bool>("runStructUnpacker")) {
  if (doStruct_) {
    structBuffer_.resize(OrbitCollection<l1Scouting::TkEm>::NBX + 1);  // FIXME magic number
    produces<OrbitCollection<l1Scouting::TkEm>>();
  }
}

ScPhase2TkEmRawToDigi::~ScPhase2TkEmRawToDigi(){};

void ScPhase2TkEmRawToDigi::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  edm::Handle<SDSRawDataCollection> scoutingRawDataCollection;
  iEvent.getByToken(rawToken_, scoutingRawDataCollection);

  if (doStruct_) {
    iEvent.put(unpackObj(*scoutingRawDataCollection, structBuffer_));
  }
}

template <typename T>
std::unique_ptr<OrbitCollection<T>> ScPhase2TkEmRawToDigi::unpackObj(const SDSRawDataCollection &feds,
                                                                      std::vector<std::vector<T>> &buffer) {
  unsigned int ntot = 0;
  for (auto &fedId : fedIDs_) {
    const FEDRawData &src = feds.FEDData(fedId);
    const uint64_t *begin = reinterpret_cast<const uint64_t *>(src.data());
    const uint64_t *end = reinterpret_cast<const uint64_t *>(src.data() + src.size());
    for (auto p = begin; p != end;) {
      if ((*p) == 0)
        continue;
      unsigned int bx = ((*p) >> 12) & 0xFFF;
      unsigned int nwords = (*p) & 0xFFF;
      ++p;
      assert(bx < OrbitCollection<T>::NBX);
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

void ScPhase2TkEmRawToDigi::unpackFromRaw(uint64_t data, std::vector<l1Scouting::TkEm> &outBuffer) {
  float pt, eta, phi, isolation =0;//, z0 = 0, dxy = 0, puppiw = 1;
  uint8_t quality;
  bool valid;
  l1tkemUnpack::readshared(data, pt, eta, phi, valid, quality, isolation);
  /*uint8_t pid = (data >> 37) & 0x7;
  if (pid > 1) {
    l1tkemUnpack::readcharged(data, z0, dxy, quality);
  } else {
    l1tkemUnpack::readneutral(data, puppiw, quality);
  }*/
  //outBuffer.emplace_back(pt, eta, phi, pid, z0, dxy, puppiw, quality);
  outBuffer.emplace_back(pt, eta, phi, valid, quality, isolation);
}

void ScPhase2TkEmRawToDigi::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScPhase2TkEmRawToDigi);
