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
#include "L1TriggerScouting/Phase2/interface/l1trkUnpack.h"

class ScPhase2TrackerTrackRawToDigi : public edm::stream::EDProducer<> {
public:
  explicit ScPhase2TrackerTrackRawToDigi(const edm::ParameterSet &);
  ~ScPhase2TrackerTrackRawToDigi() override;
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  void produce(edm::Event &, const edm::EventSetup &) override;

  template <typename T>
  std::unique_ptr<OrbitCollection<T>> unpackObj(const SDSRawDataCollection &feds, std::vector<std::vector<T>> &buffer);

  edm::EDGetTokenT<SDSRawDataCollection> rawToken_;
  std::vector<unsigned int> fedIDs_;
  unsigned int nFitPars_;

  // temporary storage
  std::vector<std::vector<l1Scouting::TTrack>> structBuffer_;

  void unpackFromRaw(uint64_t datalow, uint32_t datahigh, std::vector<l1Scouting::TTrack> &outBuffer);
};

ScPhase2TrackerTrackRawToDigi::ScPhase2TrackerTrackRawToDigi(const edm::ParameterSet &iConfig)
    : rawToken_(consumes<SDSRawDataCollection>(iConfig.getParameter<edm::InputTag>("src"))),
      fedIDs_(iConfig.getParameter<std::vector<unsigned int>>("fedIDs")),
      nFitPars_(iConfig.getParameter<unsigned int>("nFitPars")) {
  structBuffer_.resize(OrbitCollection<l1Scouting::TTrack>::NBX + 1);
  produces<OrbitCollection<l1Scouting::TTrack>>();
  produces<unsigned int>("nbx");
}

ScPhase2TrackerTrackRawToDigi::~ScPhase2TrackerTrackRawToDigi() {};

void ScPhase2TrackerTrackRawToDigi::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
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
        if ((i & 1) == 1)  // ODD tracks
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

void ScPhase2TrackerTrackRawToDigi::unpackFromRaw(uint64_t datalow,
                                                  uint32_t datahigh,
                                                  std::vector<l1Scouting::TTrack> &outBuffer) {
  unsigned int valid, rInv, phi0, chi2RPhi, tanl, z0, chi2RZ, d0, bendChi2, hitPattern, mvaQuality, phiSector;
  l1trkUnpack::read(
      datalow, datahigh, valid, rInv, phi0, chi2RPhi, tanl, z0, chi2RZ, d0, bendChi2, hitPattern, mvaQuality, phiSector);

  if (valid) {
    float ptF = l1trkUnpack::getPt(rInv);
    float rInvF = l1trkUnpack::getRinv(rInv);
    float phi0F = l1trkUnpack::getPhi0(phi0);
    float phiF = l1trkUnpack::getGlobalPhi(phi0F, phiSector);
    float tanlF = l1trkUnpack::getTanl(tanl);
    float d0F = l1trkUnpack::getD0(d0);
    float z0F = l1trkUnpack::getZ0(z0);
    float chi2RPhiF = l1trkUnpack::getChi2RPhi(chi2RPhi);
    float chi2RZF = l1trkUnpack::getChi2RZ(chi2RZ);
    int8_t charge = rInvF > 0 ? +1 : -1;
    GlobalVector momentum = l1trkUnpack::getMomentum(ptF, phiF, tanlF);
    GlobalPoint poca = l1trkUnpack::getPOCA(d0F, phiF, z0F);
    float dxyF = poca.perp();
    uint8_t nStub = l1trkUnpack::getNStubs(hitPattern);
    float chi2 = chi2RPhiF + chi2RZF;  // TODO: not fully sure about the chi2 sum
    float chi2Red = chi2 / (2 * nStub - nFitPars_);
    float mvaQualityF = l1trkUnpack::getMVAQuality(mvaQuality);
    float etaF = momentum.eta();

    // compute quality bits
    uint8_t quality = 0;
    if (nFitPars_ == 4) {
      if (ptF > 2 && nStub >= 4 && chi2Red < 15)
        quality += (1 << 0);
      if (ptF > 2 && nStub >= 6 && chi2Red < 15 && chi2 < 50)
        quality += (1 << 1);
      if (ptF > 5 && nStub >= 4)
        quality += (1 << 2);
    } else if (nFitPars_ == 5) {
      bool pocaCond = poca.x() < 1.0 && poca.x() > -1.0 && poca.y() < 1.0 && poca.y() > -1.0;
      if (ptF > 2 && nStub >= 4 && chi2Red < 15 && pocaCond)
        quality += (1 << 0);
      if (ptF > 2 && nStub >= 6 && chi2Red < 15 && chi2 < 50 && pocaCond)
        quality += (1 << 1);
      if (ptF > 5 && nStub >= 4 && pocaCond)
        quality |= (1 << 2);
    }

    outBuffer.emplace_back(ptF, etaF, phiF, z0F, dxyF, mvaQualityF, nStub, quality, charge);
  }
}

void ScPhase2TrackerTrackRawToDigi::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("rawDataCollector"));
  desc.add<std::vector<unsigned int>>("fedIDs");
  desc.add<unsigned int>("nFitPars");
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScPhase2TrackerTrackRawToDigi);
