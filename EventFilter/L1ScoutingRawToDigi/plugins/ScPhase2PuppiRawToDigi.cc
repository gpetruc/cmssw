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
#include "DataFormats/L1TParticleFlow/interface/PFCandidate.h"
#include "DataFormats/L1TParticleFlow/interface/puppi.h"


class ScPhase2PuppiRawToDigi : public edm::stream::EDProducer<> {
public:
  explicit ScPhase2PuppiRawToDigi(const edm::ParameterSet&);
  ~ScPhase2PuppiRawToDigi() override;
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  //void beginStream(edm::StreamID) override;
  void produce(edm::Event&, const edm::EventSetup&) override;
  //void endStream() override;
  std::unique_ptr<BXVector<l1t::PFCandidate>> unpackSimple(const SRDCollection & feds) ;

  edm::EDGetTokenT<SRDCollection> rawToken_;
  std::vector<unsigned int> fedIDs_;
  bool doSimple_;
};

ScPhase2PuppiRawToDigi::ScPhase2PuppiRawToDigi(const edm::ParameterSet& iConfig) :
  rawToken_(consumes<SRDCollection>(iConfig.getParameter<edm::InputTag>("src"))),
  fedIDs_(iConfig.getParameter<std::vector<unsigned int>>("fedIDs")),
  doSimple_(iConfig.getParameter<bool>("runSimpleUnpacker"))
{
  if (doSimple_)
    produces<BXVector<l1t::PFCandidate>>();
}

ScPhase2PuppiRawToDigi::~ScPhase2PuppiRawToDigi() {};

void ScPhase2PuppiRawToDigi::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<SRDCollection> scoutingRawDataCollection;
  iEvent.getByToken(rawToken_, scoutingRawDataCollection);

  if (doSimple_) {
    iEvent.put(unpackSimple(*scoutingRawDataCollection));
  }
}

std::unique_ptr<BXVector<l1t::PFCandidate>> ScPhase2PuppiRawToDigi::unpackSimple(const SRDCollection & feds) {
  std::vector<std::pair<const uint64_t *,const uint64_t *>> buffers;
  for (auto & fedId : fedIDs_) {
    const FEDRawData& src = feds.FEDData(fedId);
    buffers.emplace_back(reinterpret_cast<const uint64_t *>(src.data()), reinterpret_cast<const uint64_t *>(src.data() + src.size()));
  }

  auto ret = std::make_unique<BXVector<l1t::PFCandidate>>();
  std::vector<l1t::PFCandidate> candBuffer;
  for (int ibuff = 0, nbuffs = buffers.size(), lbuff = nbuffs - 1; buffers[ibuff].first != buffers[ibuff].second;
       ibuff = (ibuff == lbuff ? 0 : ibuff + 1)) {
    auto &pa = buffers[ibuff];
    while (pa.first != pa.second && *pa.first == 0) {
      pa.first++;
    }
    if (!*pa.first)
      continue;
    unsigned int nwords = (*pa.first) & 0xFFF;
    pa.first++;
    for (unsigned int i = 0; i < nwords; ++i, pa.first++) {
      auto l1puppi = l1ct::PuppiObj::unpack(ap_uint<64>(*pa.first));
            l1t::PFCandidate::ParticleType type;
      float mass = 0.13f;
      if (l1puppi.hwId.charged()) {
        if (l1puppi.hwId.isMuon()) {
          type = l1t::PFCandidate::Muon;
          mass = 0.105;
        } else if (l1puppi.hwId.isElectron()) {
          type = l1t::PFCandidate::Electron;
          mass = 0.005;
        } else
          type = l1t::PFCandidate::ChargedHadron;
      } else {
        type = l1puppi.hwId.isPhoton() ? l1t::PFCandidate::Photon : l1t::PFCandidate::NeutralHadron;
        mass = l1puppi.hwId.isPhoton() ? 0.0 : 0.5;
      }
      reco::Particle::PolarLorentzVector p4(l1puppi.floatPt(), l1puppi.floatEta(), l1puppi.floatPhi(), mass);
      candBuffer.emplace_back(type, l1puppi.intCharge(), p4, l1puppi.floatPuppiW(), l1puppi.intPt(), l1puppi.intEta(), l1puppi.intPhi());
      if (l1puppi.hwId.charged()) {
        candBuffer.back().setZ0(l1puppi.floatZ0());
        candBuffer.back().setDxy(l1puppi.floatDxy());
        candBuffer.back().setHwZ0(l1puppi.hwZ0());
        candBuffer.back().setHwDxy(l1puppi.hwDxy());
        candBuffer.back().setHwTkQuality(l1puppi.hwTkQuality());
      } else {
        candBuffer.back().setHwPuppiWeight(l1puppi.hwPuppiW());
        candBuffer.back().setHwEmID(l1puppi.hwEmID());
      }
      candBuffer.back().setEncodedPuppi64(l1puppi.pack().to_uint64());
    }
    ret->addBX(candBuffer.begin(), candBuffer.end());
    candBuffer.clear();
  }
  return ret;
}

void ScPhase2PuppiRawToDigi::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScPhase2PuppiRawToDigi);
