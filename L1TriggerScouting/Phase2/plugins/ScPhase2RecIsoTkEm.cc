#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/L1Scouting/interface/OrbitCollection.h"
#include "DataFormats/L1Scouting/interface/OrbitFlatTable.h"
#include "DataFormats/L1TParticleFlow/interface/L1ScoutingPuppi.h"
#include "DataFormats/L1TParticleFlow/interface/L1ScoutingTkEm.h"
#include "L1TriggerScouting/Utilities/interface/BxOffsetsFiller.h"

#include <ROOT/RVec.hxx>
#include <Math/Vector4D.h>
#include <Math/GenVector/LorentzVector.h>
#include <Math/GenVector/PtEtaPhiM4D.h>
#include <algorithm>
#include <array>
#include <iostream>

class ScPhase2RecIsoTkEm : public edm::stream::EDProducer<> {
public:
  explicit ScPhase2RecIsoTkEm(const edm::ParameterSet &);
  ~ScPhase2RecIsoTkEm() override;
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  void beginStream(edm::StreamID) override;
  void produce(edm::Event &, const edm::EventSetup &) override;
  void endStream() override;
  template <typename T, typename U>
  void runObj(const OrbitCollection<T> &src, const OrbitCollection<U> &srcTkEm, edm::Event &out);

  edm::EDGetTokenT<OrbitCollection<l1Scouting::Puppi>> structToken_;
  edm::EDGetTokenT<OrbitCollection<l1Scouting::TkEm>> structTkEmToken_;

  double minPtGamma_;
  double minDeltaR_;
  double maxDeltaR_;
  double maxIso_;

  template <typename T>
  bool isolationTkEm(float pt, float eta, float phi, const T *cands, unsigned int size) const;
};

ScPhase2RecIsoTkEm::ScPhase2RecIsoTkEm(const edm::ParameterSet &iConfig)
    : structToken_(consumes<OrbitCollection<l1Scouting::Puppi>>(iConfig.getParameter<edm::InputTag>("src"))),
      structTkEmToken_(consumes<OrbitCollection<l1Scouting::TkEm>>(iConfig.getParameter<edm::InputTag>("srcTkEm"))),
      minPtGamma_(iConfig.getParameter<double>("minPtGamma")),
      minDeltaR_(iConfig.getParameter<double>("minDeltaR")),
      maxDeltaR_(iConfig.getParameter<double>("maxDeltaR")),
      maxIso_(iConfig.getParameter<double>("maxIso")) {
  produces<OrbitCollection<l1Scouting::TkEm>>();
}

ScPhase2RecIsoTkEm::~ScPhase2RecIsoTkEm() {};

void ScPhase2RecIsoTkEm::beginStream(edm::StreamID) {}

void ScPhase2RecIsoTkEm::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  edm::Handle<OrbitCollection<l1Scouting::TkEm>> srcTkEm;
  iEvent.getByToken(structTkEmToken_, srcTkEm);

  edm::Handle<OrbitCollection<l1Scouting::Puppi>> src;
  iEvent.getByToken(structToken_, src);

  runObj(*src, *srcTkEm, iEvent);
}

void ScPhase2RecIsoTkEm::endStream() {}

template <typename T, typename U>
void ScPhase2RecIsoTkEm::runObj(const OrbitCollection<T> &src, const OrbitCollection<U> &srcTkEm, edm::Event &iEvent) {
  std::vector<std::vector<l1Scouting::TkEm>> photons_vec;

  ROOT::RVec<unsigned int> ig;
  unsigned int ntotIsoPhoton = 0;

  for (unsigned int bx = 0; bx <= OrbitCollection<T>::NBX; ++bx) {
    std::vector<l1Scouting::TkEm> photon_thisBx;

    auto rangeTkEm = srcTkEm.bxIterator(bx);
    auto sizeTkEm = rangeTkEm.size();
    const U *candsTkEm = (sizeTkEm > 0 ? &rangeTkEm.front() : nullptr);

    auto range = src.bxIterator(bx);
    auto size = range.size();
    const T *cands = (size > 0 ? &range.front() : nullptr);

    ig.clear();
    photon_thisBx.clear();
    for (unsigned int i = 0; i < sizeTkEm; ++i) {  // make list of all photons
      if (candsTkEm[i].pt() >= minPtGamma_) {
        // photon isolation
        bool isop = isolationTkEm(candsTkEm[i].pt(), candsTkEm[i].eta(), candsTkEm[i].phi(), cands, size);
        if (!isop)
          continue;

        l1Scouting::TkEm isolatedPhoton(candsTkEm[i].pt(),
                                        candsTkEm[i].eta(),
                                        candsTkEm[i].phi(),
                                        candsTkEm[i].quality(),
                                        candsTkEm[i].isolation(),
                                        i);

        ig.push_back(i);
        photon_thisBx.push_back(isolatedPhoton);
      }
    }

    photons_vec.push_back(photon_thisBx);
    ntotIsoPhoton += photon_thisBx.size();
  }

  auto outIsoPhoton = std::make_unique<OrbitCollection<l1Scouting::TkEm>>(photons_vec, ntotIsoPhoton);
  iEvent.put(std::move(outIsoPhoton));
}

template <typename T>
bool ScPhase2RecIsoTkEm::isolationTkEm(float pt, float eta, float phi, const T *cands, unsigned int size) const {
  bool passed = false;
  float psum = 0;
  for (unsigned int j = 0u; j < size; ++j) {  //loop over other particles
    float deta = eta - cands[j].eta(), dphi = ROOT::VecOps::DeltaPhi<float>(phi, cands[j].phi());
    float dr2 = deta * deta + dphi * dphi;
    if (dr2 >= minDeltaR_ && dr2 <= maxDeltaR_)
      psum += cands[j].pt();
  }
  if (psum <= maxIso_ * pt)
    passed = true;
  return passed;
}

void ScPhase2RecIsoTkEm::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src");
  desc.add<edm::InputTag>("srcTkEm");
  desc.add<double>("minPtGamma");
  desc.add<double>("maxIso");
  desc.add<double>("minDeltaR");
  desc.add<double>("maxDeltaR");
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScPhase2RecIsoTkEm);
