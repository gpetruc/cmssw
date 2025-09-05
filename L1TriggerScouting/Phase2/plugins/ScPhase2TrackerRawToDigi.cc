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
#include "DataFormats/L1TParticleFlow/interface/L1ScoutingTTrack.h"
#include "L1TriggerScouting/Phase2/interface/l1trackerUnpack.h"

class ScPhase2TrackerRawToDigi : public edm::stream::EDProducer<> {
public:
  explicit ScPhase2TrackerRawToDigi(const edm::ParameterSet &);
  ~ScPhase2TrackerRawToDigi() override;
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  void produce(edm::Event &, const edm::EventSetup &) override;

  template <typename T>
  std::unique_ptr<OrbitCollection<T>> unpackObj(const SDSRawDataCollection &feds, std::vector<std::vector<T>> &buffer);

  edm::EDGetTokenT<SDSRawDataCollection> rawToken_;
  std::vector<unsigned int> fedIDs_;

  // temporary storage
  std::vector<std::vector<l1Scouting::TTrack>> structBuffer_;

  void unpackFromRaw(uint64_t datalow, uint32_t datahigh, std::vector<l1Scouting::TTrack> &outBuffer);
};

ScPhase2TrackerRawToDigi::ScPhase2TrackerRawToDigi(const edm::ParameterSet &iConfig)
    : rawToken_(consumes<SDSRawDataCollection>(iConfig.getParameter<edm::InputTag>("src"))),
      fedIDs_(iConfig.getParameter<std::vector<unsigned int>>("fedIDs")) {
  structBuffer_.resize(OrbitCollection<l1Scouting::TTrack>::NBX + 1);
  produces<OrbitCollection<l1Scouting::TTrack>>();
  produces<unsigned int>("nbx");
}

ScPhase2TrackerRawToDigi::~ScPhase2TrackerRawToDigi() {};

void ScPhase2TrackerRawToDigi::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  edm::Handle<SDSRawDataCollection> feds;
  iEvent.getByToken(rawToken_, feds);

  unsigned int ntot = 0, nbx = 0, reforbit = iEvent.id().event();
  for (auto &fedId : fedIDs_) {
    const FEDRawData &src = feds->FEDData(fedId);
    const uint64_t *begin = reinterpret_cast<const uint64_t *>(src.data());
    const uint64_t *end = reinterpret_cast<const uint64_t *>(src.data() + src.size());
    for (auto p = begin; p != end;) {
      if ((*p) == 0) {
        ++p;
        continue;
      }
      unsigned int bx = ((*p) >> 12) & 0xFFF;
      unsigned int nwords = (*p) & 0xFFF;
      unsigned int orbit = ((*p) >> 24) & 0xFFFFFFFFFlu;
      if (reforbit != orbit) {
        throw cms::Exception("CorruptData") << "Data for orbit " << reforbit << ", fedId " << fedId
                                            << " has header with mismatching orbit number " << orbit << std::endl;
      }
      nbx++;
      unsigned int nTrackers = (2 * nwords) / 3;  // to count for the 96-bit muon words
      ++p;

      assert(bx < OrbitCollection<l1Scouting::TTrack>::NBX);  // asser fail --> unpacked wrong !
      std::vector<l1Scouting::TTrack> &outputBuffer = structBuffer_[bx + 1];
      outputBuffer.reserve(nwords);

      uint64_t datalow;
      uint32_t datahigh;

      const uint32_t *ptr32 = reinterpret_cast<const uint32_t *>(p);

      for (unsigned int i = 0; i < nTrackers; ++i, ptr32 += 3 /* jumping 96bits*/) {
        if ((i & 1) == 1)  // ODD Trackers
        {
          datalow = *reinterpret_cast<const uint64_t *>(ptr32 + 1);
          datahigh = *ptr32;
        } else {
          datalow = *reinterpret_cast<const uint64_t *>(ptr32);
          datahigh = *(ptr32 + 2);
        }
        if ((datalow == 0) and (datahigh == 0))
          continue;
        unpackFromRaw(datalow, datahigh, outputBuffer);
        ntot++;
      }
      p += nwords;
    }
  }
  iEvent.put(std::make_unique<OrbitCollection<l1Scouting::TTrack>>(structBuffer_, ntot));
  iEvent.put(std::make_unique<unsigned int>(nbx), "nbx");
}

void ScPhase2TrackerRawToDigi::unpackFromRaw(uint64_t datalow,
                                            uint32_t datahigh,
                                            std::vector<l1Scouting::TTrack> &outBuffer) {
  
  //TODO - check types, is it all supposed to be double?                                            
  double rinv, chi2RPhi, tanl, z0, chi2Rz, d0, bendChi2, mvaQuality, MVAOther;
  int16_t hitPattern;

  l1tkemUnpack::read(datalow, datahigh, rinv, phi, chi2RPhi, tanl, z0, chi2Rz, d0, bendChi2, hitPattern, mvaQuality, MVAOther);
  
  outBuffer.emplace_back(rinv, phi, chi2RPhi, tanl, z0, chi2Rz, d0, bendChi2, hitPattern, mvaQuality, MVAOther);
}

void ScPhase2TrackerRawToDigi::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("rawDataCollector"));
  desc.add<std::vector<unsigned int>>("fedIDs");
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScPhase2TrackerRawToDigi);
