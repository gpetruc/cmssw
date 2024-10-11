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
#include "DataFormats/L1TMuonPhase2/interface/L1ScoutingTrackerMuon.h"
#include "L1TriggerScouting/Phase2/interface/l1puppiUnpack.h"

class ScPhase2TrackerMuonRawToDigi : public edm::stream::EDProducer<> {
public:
  explicit ScPhase2TrackerMuonRawToDigi(const edm::ParameterSet &);
  ~ScPhase2TrackerMuonRawToDigi() override;
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  void produce(edm::Event &, const edm::EventSetup &) override;

  template <typename T>
  std::unique_ptr<OrbitCollection<T>> unpackObj(const SDSRawDataCollection &feds, std::vector<std::vector<T>> &buffer);

  edm::EDGetTokenT<SDSRawDataCollection> rawToken_;
  std::vector<unsigned int> fedIDs_;

  // temporary storage
  std::vector<std::vector<l1Scouting::TrackerMuon>> structBuffer_;

  void unpackFromRaw(uint64_t wlo, uint32_t whi, std::vector<l1Scouting::TrackerMuon> &outBuffer);
};

ScPhase2TrackerMuonRawToDigi::ScPhase2TrackerMuonRawToDigi(const edm::ParameterSet &iConfig)
    : rawToken_(consumes<SDSRawDataCollection>(iConfig.getParameter<edm::InputTag>("src"))),
      fedIDs_(iConfig.getParameter<std::vector<unsigned int>>("fedIDs")) {
  structBuffer_.resize(OrbitCollection<l1Scouting::TrackerMuon>::NBX + 1);
  produces<OrbitCollection<l1Scouting::TrackerMuon>>();
  produces<unsigned int>("nbx");
}

ScPhase2TrackerMuonRawToDigi::~ScPhase2TrackerMuonRawToDigi(){};

void ScPhase2TrackerMuonRawToDigi::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
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
      unsigned int nTrackerMuons = (2 * nwords) / 3;  // to count for the 96-bit muon words
      ++p;

      assert(bx < OrbitCollection<l1Scouting::TrackerMuon>::NBX);  // asser fail --> unpacked wrong !
      std::vector<l1Scouting::TrackerMuon> &outputBuffer = structBuffer_[bx + 1];
      outputBuffer.reserve(nwords);

      uint64_t wlo;
      uint32_t whi;

      const uint32_t *pMu = reinterpret_cast<const uint32_t *>(p);
      for (unsigned int i = 0; i < nTrackerMuons; ++i, pMu += 3 /* jumping 96bits*/) {
        if ((i & 1) == 1)  // ODD TrackerMuons
        {
          wlo = *reinterpret_cast<const uint64_t *>(pMu + 1);
          whi = *pMu;
        } else {
          wlo = *reinterpret_cast<const uint64_t *>(pMu);
          whi = *(pMu + 2);
        }
        if ((wlo == 0) and (whi == 0))
          continue;
        unpackFromRaw(wlo, whi, outputBuffer);
        ntot++;
      }
      p += nwords;
    }
  }
  iEvent.put(std::make_unique<OrbitCollection<l1Scouting::TrackerMuon>>(structBuffer_, ntot));
  iEvent.put(std::make_unique<unsigned int>(nbx), "nbx");
}

void ScPhase2TrackerMuonRawToDigi::unpackFromRaw(uint64_t wlo,
                                                 uint32_t whi,
                                                 std::vector<l1Scouting::TrackerMuon> &outBuffer) {
  float pt, eta, phi, z0 = 0, d0 = 0, beta;
  int8_t charge;
  uint8_t quality, isolation;

  pt = l1puppiUnpack::extractBitsFromW<1, 16>(wlo) * 0.03125f;
  phi = l1puppiUnpack::extractSignedBitsFromW<17, 13>(wlo) * float(M_PI / (1 << 12));
  eta = l1puppiUnpack::extractSignedBitsFromW<30, 14>(wlo) * float(M_PI / (1 << 12));
  z0 = l1puppiUnpack::extractSignedBitsFromW<44, 10>(wlo) * 0.05f;
  d0 = l1puppiUnpack::extractSignedBitsFromW<54, 10>(wlo) * 0.03f;
  quality = l1puppiUnpack::extractBitsFromW<1, 8>(whi);
  isolation = l1puppiUnpack::extractBitsFromW<9, 4>(whi);
  beta = l1puppiUnpack::extractBitsFromW<13, 4>(whi) * 0.06f;
  charge = (whi & 1) ? -1 : +1;
  outBuffer.emplace_back(pt, eta, phi, z0, d0, charge, quality, beta, isolation);
}

void ScPhase2TrackerMuonRawToDigi::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("rawDataCollector"));
  desc.add<std::vector<unsigned int>>("fedIDs");
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScPhase2TrackerMuonRawToDigi);
