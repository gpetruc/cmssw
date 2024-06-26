#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/FEDRawData/interface/FEDRawData.h"
#include "DataFormats/L1Scouting/interface/SRawDataCollection.h"
#include "DataFormats/L1Scouting/interface/SDSNumbering.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "DataFormats/L1Trigger/interface/SparseBXVector.h"
#include "DataFormats/L1TParticleFlow/interface/PFCandidate.h"
#include "DataFormats/L1TParticleFlow/interface/puppi.h"
#include "DataFormats/L1TParticleFlow/interface/struct_types.h"
#include "phase2Utils.h"

class ScPhase2PuppiRawToDigi : public edm::stream::EDProducer<> {
public:
  explicit ScPhase2PuppiRawToDigi(const edm::ParameterSet &);
  ~ScPhase2PuppiRawToDigi() override;
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  //void beginStream(edm::StreamID) override;
  void produce(edm::Event &, const edm::EventSetup &) override;
  //void endStream() override;

  template <typename BXVS>  // BXVector or SparseBXVector of l1t::PFCandidate
  std::unique_ptr<BXVS> unpackSimple(const SRDCollection &feds);

  template <typename BXVS>  // BXVector or SparseBXVector of l1ct::structs::Puppie
  std::unique_ptr<BXVS> unpackStruct(const SRDCollection &feds);

  std::unique_ptr<l1ct::structs::PuppiSOA> unpackSOA(const SRDCollection &feds);
  std::unique_ptr<l1ct::structs::PuppiASOA> vunpackSOA(const SRDCollection &feds);

  template <typename T>
  void appendToBXV(BXVector<T> &to, const std::vector<T> &from, int bx) {
    if (bx == to.getLastBX())
      to.append(from.begin(), from.end());
    else
      to.addBX(from.begin(), from.end());
  }

  template <typename T>
  void appendToBXV(SparseBXVector<T> &to, const std::vector<T> &from, int bx) {
    to.addBX(bx, from.begin(), from.end());
  }

  edm::EDGetTokenT<SRDCollection> rawToken_;
  std::vector<unsigned int> fedIDs_;
  bool doSimple_, doStruct_, doSOA_, dovSOA_;
  bool noWrite_, sparseBXVector_;
};

ScPhase2PuppiRawToDigi::ScPhase2PuppiRawToDigi(const edm::ParameterSet &iConfig)
    : rawToken_(consumes<SRDCollection>(iConfig.getParameter<edm::InputTag>("src"))),
      fedIDs_(iConfig.getParameter<std::vector<unsigned int>>("fedIDs")),
      doSimple_(iConfig.getParameter<bool>("runSimpleUnpacker")),
      doStruct_(iConfig.getParameter<bool>("runStructUnpacker")),
      doSOA_(iConfig.getParameter<bool>("runSOAUnpacker")),
      dovSOA_(doSOA_ ? iConfig.getParameter<bool>("vectorizeSOAUnpacker") : false),
      noWrite_(iConfig.getParameter<bool>("noWrite")),
      sparseBXVector_(iConfig.getParameter<bool>("sparseBXVector")) {
  if (doSimple_) {
    if (sparseBXVector_)
      produces<SparseBXVector<l1t::PFCandidate>>();
    else
      produces<BXVector<l1t::PFCandidate>>();
  }
  if (doStruct_) {
    if (sparseBXVector_)
      produces<SparseBXVector<l1ct::structs::Puppi>>();
    else
      produces<BXVector<l1ct::structs::Puppi>>();
  }
  if (doSOA_) {
    if (dovSOA_) 
      produces<l1ct::structs::PuppiASOA>();
    else
      produces<l1ct::structs::PuppiSOA>();
  }
}

ScPhase2PuppiRawToDigi::~ScPhase2PuppiRawToDigi(){};

void ScPhase2PuppiRawToDigi::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  edm::Handle<SRDCollection> scoutingRawDataCollection;
  iEvent.getByToken(rawToken_, scoutingRawDataCollection);

  if (doSimple_) {
    if (sparseBXVector_)
      iEvent.put(unpackSimple<SparseBXVector<l1t::PFCandidate>>(*scoutingRawDataCollection));
    else
      iEvent.put(unpackSimple<BXVector<l1t::PFCandidate>>(*scoutingRawDataCollection));
  }
  if (doStruct_) {
    if (sparseBXVector_)
      iEvent.put(unpackStruct<SparseBXVector<l1ct::structs::Puppi>>(*scoutingRawDataCollection));
    else
      iEvent.put(unpackStruct<BXVector<l1ct::structs::Puppi>>(*scoutingRawDataCollection));
  }
  if (doSOA_) {
    if (dovSOA_) 
      iEvent.put(vunpackSOA(*scoutingRawDataCollection));
    else
      iEvent.put(unpackSOA(*scoutingRawDataCollection));
  }
}

template <typename BXVS>
std::unique_ptr<BXVS> ScPhase2PuppiRawToDigi::unpackSimple(const SRDCollection &feds) {
  std::vector<std::pair<const uint64_t *, const uint64_t *>> buffers;
  for (auto &fedId : fedIDs_) {
    const FEDRawData &src = feds.FEDData(fedId);
    buffers.emplace_back(reinterpret_cast<const uint64_t *>(src.data()),
                         reinterpret_cast<const uint64_t *>(src.data() + src.size()));
  }

  auto ret = std::make_unique<BXVS>();
  std::vector<l1t::PFCandidate> candBuffer;
  for (int ibuff = 0, nbuffs = buffers.size(), lbuff = nbuffs - 1; buffers[ibuff].first != buffers[ibuff].second;
       ibuff = (ibuff == lbuff ? 0 : ibuff + 1)) {
    auto &pa = buffers[ibuff];
    while (pa.first != pa.second && *pa.first == 0) {
      pa.first++;
    }
    if (!*pa.first)
      continue;
    unsigned int bx = ((*pa.first) >> 12) & 0xFFF;
    unsigned int nwords = (*pa.first) & 0xFFF;
    pa.first++;
    float pt, eta, phi, mass, z0, dxy, puppiw;
    uint16_t hwPt, hwPuppiW;
    int16_t pdgId, hwEta, hwPhi, hwZ0;
    int8_t hwDxy;
    uint8_t pid, hwQuality;
    l1t::PFCandidate::ParticleType type;
    int charge;
    for (unsigned int i = 0; i < nwords; ++i, pa.first++) {
      uint64_t data = *pa.first;
      phase2Utils::readshared(data, pt, eta, phi);
      phase2Utils::readshared(data, hwPt, hwEta, hwPhi);
      pid = (data >> 37) & 0x7;
      phase2Utils::assignpdgid(pid, pdgId);
      phase2Utils::assignCMSSWPFCandidateId(pid, type);
      phase2Utils::assignmass(pid, mass);
      phase2Utils::assigncharge(pid, charge);
      reco::Particle::PolarLorentzVector p4(pt, eta, phi, mass);
      if (pid > 1) {
        phase2Utils::readcharged(data, z0, dxy, hwQuality);
        phase2Utils::readcharged(data, hwZ0, hwDxy, hwQuality);
        puppiw = 1.0;
      } else {
        phase2Utils::readneutral(data, puppiw, hwQuality);
        phase2Utils::readneutral(data, hwPuppiW, hwQuality);
        dxy = 0;
        z0 = 0;
      }
      candBuffer.emplace_back(type, charge, p4, puppiw, hwPt, hwEta, hwPhi);
      if (pid > 1) {
        candBuffer.back().setZ0(z0);
        candBuffer.back().setDxy(dxy);
        candBuffer.back().setHwZ0(hwZ0);
        candBuffer.back().setHwDxy(hwDxy);
        candBuffer.back().setHwTkQuality(hwQuality);
      } else {
        candBuffer.back().setHwPuppiWeight(hwPuppiW);
        candBuffer.back().setHwEmID(hwQuality);
      }
      candBuffer.back().setEncodedPuppi64(data);
    }
    if (!noWrite_) {
      appendToBXV(*ret, candBuffer, bx);
    }
    candBuffer.clear();
  }
  return ret;
}

template <typename BXVS>
std::unique_ptr<BXVS> ScPhase2PuppiRawToDigi::unpackStruct(const SRDCollection &feds) {
  std::vector<std::pair<const uint64_t *, const uint64_t *>> buffers;
  for (auto &fedId : fedIDs_) {
    const FEDRawData &src = feds.FEDData(fedId);
    buffers.emplace_back(reinterpret_cast<const uint64_t *>(src.data()),
                         reinterpret_cast<const uint64_t *>(src.data() + src.size()));
  }

  auto ret = std::make_unique<BXVS>();
  std::vector<l1ct::structs::Puppi> candBuffer;
  for (int ibuff = 0, nbuffs = buffers.size(), lbuff = nbuffs - 1; buffers[ibuff].first != buffers[ibuff].second;
       ibuff = (ibuff == lbuff ? 0 : ibuff + 1)) {
    auto &pa = buffers[ibuff];
    while (pa.first != pa.second && *pa.first == 0) {
      pa.first++;
    }
    if (!*pa.first)
      continue;
    unsigned int bx = ((*pa.first) >> 12) & 0xFFF;
    unsigned int nwords = (*pa.first) & 0xFFF;
    pa.first++;
    candBuffer.resize(nwords);
    for (unsigned int i = 0; i < nwords; ++i, pa.first++) {
      uint64_t data = *pa.first;
      phase2Utils::readshared(data, candBuffer[i].pt, candBuffer[i].eta, candBuffer[i].phi);
      uint8_t pid = (data >> 37) & 0x7;
      phase2Utils::assignpdgid(pid, candBuffer[i].pdgId);
      if (pid > 1) {
        phase2Utils::readcharged(data, candBuffer[i].z0, candBuffer[i].dxy, candBuffer[i].quality);
        candBuffer[i].puppiw = 1.0;
      } else {
        phase2Utils::readneutral(data, candBuffer[i].puppiw, candBuffer[i].quality);
        candBuffer[i].dxy = 0;
        candBuffer[i].z0 = 0;
      }
    }
    if (!noWrite_) {
      appendToBXV(*ret, candBuffer, bx);
    }
    candBuffer.clear();
  }
  return ret;
}

std::unique_ptr<l1ct::structs::PuppiSOA> ScPhase2PuppiRawToDigi::unpackSOA(const SRDCollection &feds) {
  std::vector<std::pair<const uint64_t *, const uint64_t *>> buffers;
  unsigned int sizeguess = 0;
  for (auto &fedId : fedIDs_) {
    const FEDRawData &src = feds.FEDData(fedId);
    buffers.emplace_back(reinterpret_cast<const uint64_t *>(src.data()),
                         reinterpret_cast<const uint64_t *>(src.data() + src.size()));
    sizeguess += src.size();
  }
  l1ct::structs::PuppiSOA ret;
  ret.bx.reserve(3564);
  ret.offsets.reserve(3564 + 1);
  for (std::vector<float> *v : {&ret.pt, &ret.eta, &ret.phi, &ret.z0, &ret.dxy, &ret.puppiw}) {
    v->resize(sizeguess);
  }
  ret.pdgId.resize(sizeguess);
  ret.quality.resize(sizeguess);
  unsigned int i0 = 0;
  for (int ibuff = 0, nbuffs = buffers.size(), lbuff = nbuffs - 1; buffers[ibuff].first != buffers[ibuff].second;
       ibuff = (ibuff == lbuff ? 0 : ibuff + 1)) {
    auto &pa = buffers[ibuff];
    while (pa.first != pa.second && *pa.first == 0) {
      pa.first++;
    }
    if (!*pa.first)
      continue;
    unsigned int bx = ((*pa.first) >> 12) & 0xFFF;
    unsigned int nwords = (*pa.first) & 0xFFF;
    pa.first++;
    ret.bx.push_back(bx);
    ret.offsets.push_back(i0);
    for (unsigned int i = 0; i < nwords; ++i, ++pa.first, ++i0) {
      uint64_t data = *pa.first;
      phase2Utils::readshared(data, ret.pt[i0], ret.eta[i0], ret.phi[i0]);
      uint8_t pid = (data >> 37) & 0x7;
      phase2Utils::assignpdgid(pid, ret.pdgId[i0]);
      if (pid > 1) {
        phase2Utils::readcharged(data, ret.z0[i0], ret.dxy[i0], ret.quality[i0]);
        ret.puppiw[i0] = 1.0f;
      } else {
        phase2Utils::readneutral(data, ret.puppiw[i0], ret.quality[i0]);
        ret.dxy[i0] = 0.0f;
        ret.z0[i0] = 0.0f;
      }
    }
  }
  ret.offsets.push_back(i0);
  for (std::vector<float> *v : {&ret.pt, &ret.eta, &ret.phi, &ret.z0, &ret.dxy, &ret.puppiw}) {
    v->resize(i0);
  }
  ret.pdgId.resize(i0);
  ret.quality.resize(i0);
  auto retptr = std::make_unique<l1ct::structs::PuppiSOA>(std::move(ret));
  return retptr;
}

std::unique_ptr<l1ct::structs::PuppiASOA> ScPhase2PuppiRawToDigi::vunpackSOA(const SRDCollection &feds) {
  std::vector<std::pair<const uint64_t *, const uint64_t *>> buffers;
  unsigned int sizeguess = 0;
  for (auto &fedId : fedIDs_) {
    const FEDRawData &src = feds.FEDData(fedId);
    buffers.emplace_back(reinterpret_cast<const uint64_t *>(src.data()),
                         reinterpret_cast<const uint64_t *>(src.data() + src.size()));
    sizeguess += src.size();
  }
  l1ct::structs::PuppiASOA ret;
  ret.bx.reserve(3564);
  ret.offsets.reserve(3564 + 1);
  std::vector<uint64_t> packed(sizeguess);
  unsigned int i0 = 0;
  for (int ibuff = 0, nbuffs = buffers.size(), lbuff = nbuffs - 1; buffers[ibuff].first != buffers[ibuff].second;
       ibuff = (ibuff == lbuff ? 0 : ibuff + 1)) {
    auto &pa = buffers[ibuff];
    while (pa.first != pa.second && *pa.first == 0) {
      pa.first++;
    }
    if (!*pa.first)
      continue;
    unsigned int bx = ((*pa.first) >> 12) & 0xFFF;
    unsigned int nwords = (*pa.first) & 0xFFF;
    pa.first++;
    ret.bx.push_back(bx);
    ret.offsets.push_back(i0);
    packed.insert(packed.end(), pa.first, pa.first + nwords);
    pa.first += nwords;
  }
  ret.offsets.push_back(i0);
  for (l1ct::structs::PuppiASOA::avector<float> *v : {&ret.pt, &ret.eta, &ret.phi, &ret.z0, &ret.dxy, &ret.puppiw}) {
    v->resize(i0);
  }
  ret.pdgId.resize(i0);
  ret.quality.resize(i0);
  phase2Utils::readpuppi(i0, packed.data(), ret.pt.data(),  ret.eta.data(), ret.phi.data(), ret.pdgId.data(), ret.quality.data(), ret.z0.data(), ret.dxy.data(), ret.puppiw.data());
  auto retptr = std::make_unique<l1ct::structs::PuppiASOA>(std::move(ret));
  return retptr;
}
void ScPhase2PuppiRawToDigi::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScPhase2PuppiRawToDigi);
