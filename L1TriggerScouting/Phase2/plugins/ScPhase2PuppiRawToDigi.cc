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
#include "DataFormats/L1TParticleFlow/interface/L1ScoutingPuppi.h"
#include "L1TriggerScouting/Phase2/interface/l1puppiUnpack.h"

class ScPhase2PuppiRawToDigi : public edm::stream::EDProducer<> {
public:
  explicit ScPhase2PuppiRawToDigi(const edm::ParameterSet &);
  ~ScPhase2PuppiRawToDigi() override;
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  void produce(edm::Event &, const edm::EventSetup &) override;

  template <typename T>
  std::unique_ptr<OrbitCollection<T>> unpackObj(unsigned int orbit,
                                                const SDSRawDataCollection &feds,
                                                std::vector<std::vector<T>> &buffer);

  edm::EDGetTokenT<SDSRawDataCollection> rawToken_;
  std::vector<unsigned int> fedIDs_;
  uint8_t splitFactor_; // number of fragments per BX
  bool doCandidate_, doStruct_;

  // temporary storage
  std::vector<std::vector<l1t::PFCandidate>> candBuffer_;
  std::vector<std::vector<l1Scouting::Puppi>> structBuffer_;
  unsigned int nbx_;

  void unpackFromRaw(uint64_t data, std::vector<l1t::PFCandidate> &outBuffer);
  void unpackFromRaw(uint64_t data, std::vector<l1Scouting::Puppi> &outBuffer);
};

ScPhase2PuppiRawToDigi::ScPhase2PuppiRawToDigi(const edm::ParameterSet &iConfig)
    : rawToken_(consumes<SDSRawDataCollection>(iConfig.getParameter<edm::InputTag>("src"))),
      fedIDs_(iConfig.getParameter<std::vector<unsigned int>>("fedIDs")),
      splitFactor_(iConfig.getParameter<unsigned int>("splitFactor")),
      doCandidate_(iConfig.getParameter<bool>("runCandidateUnpacker")),
      doStruct_(iConfig.getParameter<bool>("runStructUnpacker")) {
  if (doCandidate_) {
    candBuffer_.resize(OrbitCollection<l1t::PFCandidate>::NBX + 1);
    produces<OrbitCollection<l1t::PFCandidate>>();
  }
  if (doStruct_) {
    structBuffer_.resize(OrbitCollection<l1Scouting::Puppi>::NBX + 1);
    produces<OrbitCollection<l1Scouting::Puppi>>();
  }
  if (doCandidate_ || doStruct_) {
    produces<unsigned int>("nbx");
  }
}

ScPhase2PuppiRawToDigi::~ScPhase2PuppiRawToDigi() {};

void ScPhase2PuppiRawToDigi::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  edm::Handle<SDSRawDataCollection> scoutingRawDataCollection;
  iEvent.getByToken(rawToken_, scoutingRawDataCollection);
  if (doCandidate_) {
    iEvent.put(unpackObj(iEvent.id().event(), *scoutingRawDataCollection, candBuffer_));
  }
  if (doStruct_) {
    iEvent.put(unpackObj(iEvent.id().event(), *scoutingRawDataCollection, structBuffer_));
  }
  if (doCandidate_ || doStruct_) {
    iEvent.put(std::make_unique<unsigned int>(nbx_), "nbx");
  }
}

template <typename T>
std::unique_ptr<OrbitCollection<T>> ScPhase2PuppiRawToDigi::unpackObj(unsigned int orbit,
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
        throw cms::Exception("CorruptData") << "Data for orbit " << orbit 
                                            << " has " << nfound << " blocks for bx " << bx << ", expected " << splitFactor_ << std::endl;
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

void ScPhase2PuppiRawToDigi::unpackFromRaw(uint64_t data, std::vector<l1t::PFCandidate> &outBuffer) {
  float pt, eta, phi, mass, z0 = 0, dxy = 0, puppiw = 1;
  uint16_t hwPt, hwPuppiW = 1 << 8;
  int16_t pdgId, hwEta, hwPhi, hwZ0 = 0;
  int8_t hwDxy = 0;
  uint8_t pid, hwQuality;
  l1t::PFCandidate::ParticleType type;
  int charge;
  l1puppiUnpack::readshared(data, pt, eta, phi);
  l1puppiUnpack::readshared(data, hwPt, hwEta, hwPhi);
  pid = (data >> 37) & 0x7;
  l1puppiUnpack::assignpdgid(pid, pdgId);
  l1puppiUnpack::assignCMSSWPFCandidateId(pid, type);
  l1puppiUnpack::assignmass(pid, mass);
  l1puppiUnpack::assigncharge(pid, charge);
  reco::Particle::PolarLorentzVector p4(pt, eta, phi, mass);
  if (pid > 1) {
    l1puppiUnpack::readcharged(data, z0, dxy, hwQuality);
    l1puppiUnpack::readcharged(data, hwZ0, hwDxy, hwQuality);
  } else {
    l1puppiUnpack::readneutral(data, puppiw, hwQuality);
    l1puppiUnpack::readneutral(data, hwPuppiW, hwQuality);
  }
  outBuffer.emplace_back(type, charge, p4, puppiw, hwPt, hwEta, hwPhi);
  if (pid > 1) {
    outBuffer.back().setZ0(z0);
    outBuffer.back().setDxy(dxy);
    outBuffer.back().setHwZ0(hwZ0);
    outBuffer.back().setHwDxy(hwDxy);
    outBuffer.back().setHwTkQuality(hwQuality);
  } else {
    outBuffer.back().setHwPuppiWeight(hwPuppiW);
    outBuffer.back().setHwEmID(hwQuality);
  }
  outBuffer.back().setEncodedPuppi64(data);
}

void ScPhase2PuppiRawToDigi::unpackFromRaw(uint64_t data, std::vector<l1Scouting::Puppi> &outBuffer) {
  float pt, eta, phi, z0 = 0, dxy = 0, puppiw = 1;
  uint8_t quality;
  l1puppiUnpack::readshared(data, pt, eta, phi);
  uint8_t pid = (data >> 37) & 0x7;
  if (pid > 1) {
    l1puppiUnpack::readcharged(data, z0, dxy, quality);
  } else {
    l1puppiUnpack::readneutral(data, puppiw, quality);
  }
  outBuffer.emplace_back(pt, eta, phi, pid, z0, dxy, puppiw, quality);
}

void ScPhase2PuppiRawToDigi::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("rawDataCollector"));
  desc.add<std::vector<unsigned int>>("fedIDs");
  desc.add<unsigned int>("splitFactor", 1)->setComment("Number of fragments per BX");
  desc.add<bool>("runCandidateUnpacker", false);
  desc.add<bool>("runStructUnpacker", true);
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScPhase2PuppiRawToDigi);
