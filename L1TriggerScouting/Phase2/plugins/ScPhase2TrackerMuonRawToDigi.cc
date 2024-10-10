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
#include "DataFormats/L1TParticleFlow/interface/PFCandidate.h"
#include "DataFormats/L1TMuonPhase2/interface/L1ScoutingTrackerMuon.h"
#include "L1TriggerScouting/Phase2/interface/l1puppiUnpack.h"

class ScPhase2TrackerMuonRawToDigi : public edm::stream::EDProducer<> {
public:
  explicit ScPhase2TrackerMuonRawToDigi(const edm::ParameterSet &);
  ~ScPhase2TrackerMuonRawToDigi() override;
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  //void beginStream(edm::StreamID) override;
  void produce(edm::Event &, const edm::EventSetup &) override;
  //void endStream() override;

  template <typename T>
  std::unique_ptr<OrbitCollection<T>> unpackObj(const SDSRawDataCollection &feds, std::vector<std::vector<T>> &buffer);

  edm::EDGetTokenT<SDSRawDataCollection> rawToken_;
  std::vector<unsigned int> fedIDs_;
  bool doCandidate_, doStruct_;

  // temporary storage
  std::vector<std::vector<l1t::PFCandidate>> candBuffer_;
  std::vector<std::vector<l1Scouting::TrackerMuon>> structBuffer_;

  void unpackFromRaw(uint64_t wlo, uint32_t whi, std::vector<l1Scouting::TrackerMuon> &outBuffer);
};

ScPhase2TrackerMuonRawToDigi::ScPhase2TrackerMuonRawToDigi(const edm::ParameterSet &iConfig)
    : rawToken_(consumes<SDSRawDataCollection>(iConfig.getParameter<edm::InputTag>("src"))),
      fedIDs_(iConfig.getParameter<std::vector<unsigned int>>("fedIDs")),
      doCandidate_(iConfig.getParameter<bool>("runCandidateUnpacker")),
      doStruct_(iConfig.getParameter<bool>("runStructUnpacker")) {
  if (doCandidate_) {
    produces<OrbitCollection<l1t::PFCandidate>>();
    candBuffer_.resize(OrbitCollection<l1t::PFCandidate>::NBX + 1);  // FIXME magic number
  }
  if (doStruct_) {
    structBuffer_.resize(OrbitCollection<l1Scouting::TrackerMuon>::NBX + 1);  // FIXME magic number
    produces<OrbitCollection<l1Scouting::TrackerMuon>>();
  }
}

ScPhase2TrackerMuonRawToDigi::~ScPhase2TrackerMuonRawToDigi(){};

void ScPhase2TrackerMuonRawToDigi::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  edm::Handle<SDSRawDataCollection> scoutingRawDataCollection;
  iEvent.getByToken(rawToken_, scoutingRawDataCollection);

  if (doStruct_) {
    iEvent.put(unpackObj(*scoutingRawDataCollection, structBuffer_));
  }
}

template <typename T>
std::unique_ptr<OrbitCollection<T>> ScPhase2TrackerMuonRawToDigi::unpackObj(const SDSRawDataCollection &feds,
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
      unsigned int nTrackerMuons = 2 * nwords / 3;  // tocount for the 96-bit muon words
      ++p;

      assert(bx < OrbitCollection<T>::NBX);  // asser fail --> unpacked wrong !
      std::vector<T> &outputBuffer = buffer[bx + 1];
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
  return std::make_unique<OrbitCollection<T>>(buffer, ntot);
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
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScPhase2TrackerMuonRawToDigi);
