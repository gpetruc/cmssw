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

//CHANGES TO IMPLEMENT
//- RETURN THE FULL 4 particles ?

class ScPhase2PhotonIsolation : public edm::stream::EDProducer<> {
public:
  explicit ScPhase2PhotonIsolation(const edm::ParameterSet &);
  ~ScPhase2PhotonIsolation() override;
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  void beginStream(edm::StreamID) override;
  void produce(edm::Event &, const edm::EventSetup &) override;
  void endStream() override;
  template <typename T, typename U>
  void runObj(const OrbitCollection<T> &src,
              const OrbitCollection<U> &srcTkEm,
              edm::Event &out,
              const std::string &bxLabel);

  bool doStruct_;
  edm::EDGetTokenT<OrbitCollection<l1Scouting::Puppi>> structToken_;
  edm::EDGetTokenT<OrbitCollection<l1Scouting::TkEm>> structTkEmToken_;

  struct Cuts {
    float minpt3 = 1;
    float mindr2 = 1;
    float maxdr2 = 2;
    float maxdeltar2 = 1;
    float mindr2tkem = 1;
    float maxdr2tkem = 1;
    float maxiso = 1;
    float maxisotkem = 1;
  } cuts;

  template <typename T>
  bool isolationTkEm(float pt, float eta, float phi, const T *cands, unsigned int size) const;

  unsigned long countStruct_;
  unsigned long passStruct_;
};

ScPhase2PhotonIsolation::ScPhase2PhotonIsolation(const edm::ParameterSet &iConfig)
    : doStruct_(iConfig.getParameter<bool>("runStruct")) {
  if (doStruct_) {
    //PUPPI input being given here
    structToken_ = consumes<OrbitCollection<l1Scouting::Puppi>>(iConfig.getParameter<edm::InputTag>("src"));
    structTkEmToken_ = consumes<OrbitCollection<l1Scouting::TkEm>>(iConfig.getParameter<edm::InputTag>("srcTkEm"));

    produces<OrbitCollection<l1Scouting::TkEmIsolated>>();
    produces<std::vector<unsigned>>("selectedBx");
    produces<unsigned int>("nbx");
  }
}

ScPhase2PhotonIsolation::~ScPhase2PhotonIsolation() {};

void ScPhase2PhotonIsolation::beginStream(edm::StreamID) {}

void ScPhase2PhotonIsolation::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  if (doStruct_) {
    edm::Handle<OrbitCollection<l1Scouting::TkEm>> srcTkEm;
    iEvent.getByToken(structTkEmToken_, srcTkEm);

    edm::Handle<OrbitCollection<l1Scouting::Puppi>> src;
    iEvent.getByToken(structToken_, src);

    runObj(*src, *srcTkEm, iEvent, "");
  }
}

void ScPhase2PhotonIsolation::endStream() {}

template <typename T, typename U>
void ScPhase2PhotonIsolation::runObj(const OrbitCollection<T> &src,
                                    const OrbitCollection<U> &srcTkEm,
                                    edm::Event &iEvent,
                                    const std::string &label) {
  auto ret = std::make_unique<std::vector<unsigned>>();
  auto selectedBx = std::make_unique<std::vector<unsigned>>();

  std::vector<std::vector<l1Scouting::TkEmIsolated>> photons_vec;

  ROOT::RVec<unsigned int> ig;
  unsigned int nbx = 0, ntotIsoPhoton = 0;

  for (unsigned int bx = 0; bx <= OrbitCollection<T>::NBX; ++bx) {
    nbx++;
    std::vector<l1Scouting::TkEmIsolated> photon_thisBx;

    auto rangeTkEm = srcTkEm.bxIterator(bx);
    const U *candsTkEm = &rangeTkEm.front();
    auto sizeTkEm = rangeTkEm.size();

    auto range = src.bxIterator(bx);
    const T *cands = &range.front();
    auto size = range.size();

    ig.clear();
    photon_thisBx.clear();
    for (unsigned int i = 0; i < sizeTkEm; ++i) {  // make list of all photons
      if (candsTkEm[i].pt() >= cuts.minpt3) {

        // photon isolation
        bool isop = isolationTkEm(
          candsTkEm[i].pt(), candsTkEm[i].eta(), candsTkEm[i].phi(), cands, size);
        if (!isop)
          continue;    

        l1Scouting::TkEmIsolated isolatedPhoton(
          candsTkEm[i].pt(),
          candsTkEm[i].eta(),
          candsTkEm[i].phi(),
          candsTkEm[i].quality(),
          candsTkEm[i].isolation(),
          i
        );

        ig.push_back(i);
        photon_thisBx.push_back(isolatedPhoton);
      }
    }

    photons_vec.push_back(photon_thisBx);
    ntotIsoPhoton++;
  } 

  auto outIsoPhoton = std::make_unique<OrbitCollection<l1Scouting::TkEmIsolated>>(photons_vec, ntotIsoPhoton);
  iEvent.put(std::move(outIsoPhoton));
  iEvent.put(std::make_unique<unsigned int>(nbx), "nbx");
  iEvent.put(std::move(selectedBx), "selectedBx");
}

template <typename T>
bool ScPhase2PhotonIsolation::isolationTkEm(
    float pt, float eta, float phi, const T *cands, unsigned int size) const {
  bool passed = false;
  float psum = 0;
  for (unsigned int j = 0u; j < size; ++j) {  //loop over other particles
    float deta = eta - cands[j].eta(), dphi = ROOT::VecOps::DeltaPhi<float>(phi, cands[j].phi());
    float dr2 = deta * deta + dphi * dphi;
    if (dr2 >= cuts.mindr2tkem && dr2 <= cuts.maxdr2tkem)
      psum += cands[j].pt();
  }
  if (psum <= cuts.maxisotkem * pt)
    passed = true;
  return passed;
}

void ScPhase2PhotonIsolation::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src");
  desc.add<edm::InputTag>("srcTkEm");
  desc.add<bool>("runStruct", true);
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScPhase2PhotonIsolation);
