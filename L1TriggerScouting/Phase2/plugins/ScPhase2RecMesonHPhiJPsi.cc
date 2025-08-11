#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/L1Scouting/interface/OrbitCollection.h"
#include "DataFormats/L1Scouting/interface/OrbitFlatTable.h"
#include "DataFormats/L1TParticleFlow/interface/RecMeson.h"
#include "DataFormats/L1TParticleFlow/interface/L1ScoutingTkEm.h"
#include "L1TriggerScouting/Utilities/interface/BxOffsetsFiller.h"

#include <ROOT/RVec.hxx>
#include <Math/Vector4D.h>
#include <Math/GenVector/LorentzVector.h>
#include <Math/GenVector/PtEtaPhiM4D.h>
#include <algorithm>
#include <array>
#include <iostream>

class ScPhase2RecMesonHPhiJPsi : public edm::stream::EDProducer<> {
public:
  explicit ScPhase2RecMesonHPhiJPsi(const edm::ParameterSet &);
  ~ScPhase2RecMesonHPhiJPsi() override;
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  void beginStream(edm::StreamID) override;
  void produce(edm::Event &, const edm::EventSetup &) override;
  void endStream() override;
  template <typename T>
  void runObj(const OrbitCollection<T> &srcPhi,
              const OrbitCollection<T> &srcJPsi,
              edm::Event &out,
              unsigned long &nTry,
              unsigned long &nPass,
              const std::string &bxLabel);

  bool doStruct_;
  edm::EDGetTokenT<OrbitCollection<l1Scouting::RecMeson>> structToken_;

  struct Cuts {
    float minmassH = 100;
    float maxmassH = 150;
    float minptQ = 30;
    float maxiso = 0.25;
  } cuts;

  static float quadrimass(const l1Scouting::RecMeson *candsa, int a, 
                   const l1Scouting::RecMeson *candsb, int b);

  unsigned long countStruct_;
  unsigned long passStruct_;
};

ScPhase2RecMesonHPhiJPsi::ScPhase2RecMesonHPhiJPsi(const edm::ParameterSet &iConfig)
    : doStruct_(iConfig.getParameter<bool>("runStruct")) {
  if (doStruct_) {
    structToken_ = consumes<OrbitCollection<l1Scouting::RecMeson>>(iConfig.getParameter<edm::InputTag>("srcPhi"));
    structToken_ = consumes<OrbitCollection<l1Scouting::RecMeson>>(iConfig.getParameter<edm::InputTag>("srcJPsi"));
    produces<std::vector<unsigned>>("selectedBx");
    produces<l1ScoutingRun3::OrbitFlatTable>("recMesonHphijpsi");
  }
}

ScPhase2RecMesonHPhiJPsi::~ScPhase2RecMesonHPhiJPsi() {};

void ScPhase2RecMesonHPhiJPsi::beginStream(edm::StreamID) {
  countStruct_ = 0;
  passStruct_ = 0;
}

void ScPhase2RecMesonHPhiJPsi::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  if (doStruct_) {
    edm::Handle<OrbitCollection<l1Scouting::RecMeson>> srcPhi;
    edm::Handle<OrbitCollection<l1Scouting::RecMeson>> srcJPsi;
    iEvent.getByToken(structToken_, srcPhi);
    iEvent.getByToken(structToken_, srcJPsi);
    runObj(*srcPhi, *srcJPsi, iEvent, countStruct_, passStruct_, "");
  }
}

void ScPhase2RecMesonHPhiJPsi::endStream() {
  if (doStruct_)
    edm::LogImportant("ScPhase2AnalysisSummary") << "Rec Meson HPhiJPsi Struct analysis: " << countStruct_ << " -> " << passStruct_;
}

template <typename T>
void ScPhase2RecMesonHPhiJPsi::runObj(const OrbitCollection<T> &srcPhi,
                                    const OrbitCollection<T> &srcJPsi,
                                    edm::Event &iEvent,
                                    unsigned long &nTry,
                                    unsigned long &nPass,
                                    const std::string &label) {
  l1ScoutingRun3::BxOffsetsFillter bxOffsetsFiller;
  bxOffsetsFiller.start();
  auto ret = std::make_unique<std::vector<unsigned>>();
  std::vector<float> masses;
  std::vector<uint8_t> i0s, i1s, i2s, i3s;
  std::array<unsigned int, 2> bestMesonPair;
  float bestMesonPairScore;
  bool bestMesonPairFound;

  for (unsigned int bx = 1; bx <= OrbitCollection<T>::NBX; ++bx) {
    nTry++;

    auto rangePhi = srcPhi.bxIterator(bx);
    const T *candsPhi = &rangePhi.front();
    unsigned int nPhi = rangePhi.size();

    auto rangeJPsi = srcJPsi.bxIterator(bx);
    const T *candsJPsi = &rangeJPsi.front();
    unsigned int nJPsi = rangeJPsi.size();

    bestMesonPairScore = 0.;
    bestMesonPairFound = false;

    if (nPhi < 1 || nJPsi < 1) continue;

    for (unsigned int i1 = 0; i1 < nPhi; ++i1) {
      // minimum pt and isolation of Q1
      if ((candsPhi[i1].pt() < cuts.minptQ) || (candsPhi[i1].isoDR0p25() >= cuts.maxiso))
        continue;
      for (unsigned int i2 = i1 + 1; i2 < nJPsi; ++i2) {
        // minimum pt and isolation of Q2
        if ((candsJPsi[i2].pt() < cuts.minptQ) || (candsJPsi[i2].isoDR0p25() >= cuts.maxiso))
          continue;

        // Four different dauther particles
        if ((candsPhi[i1].id1() == candsJPsi[i2].id1()) || (candsPhi[i1].id1() == candsJPsi[i2].id2()))
          continue;
        if ((candsPhi[i1].id2() == candsJPsi[i2].id1()) || (candsPhi[i1].id2() == candsJPsi[i2].id2()))
          continue;

        // Choose best pair of mesons based on score (e.g. max pt)
        float ptsum = candsPhi[i1].pt() + candsJPsi[i2].pt();
        if (ptsum > bestMesonPairScore) {
          bestMesonPairScore = ptsum;
          bestMesonPair = {{i1, i2}};
          bestMesonPairFound = true;
        }
      }
    }

    if (!bestMesonPairFound)
      continue;

    // H mass
    auto mass = quadrimass(candsPhi, bestMesonPair[0], candsJPsi, bestMesonPair[1]);
    if (!(mass >= cuts.minmassH and mass <= cuts.maxmassH))
      continue;    

    ret->emplace_back(bx);
    nPass++;
    masses.push_back(mass);
    i0s.push_back(candsPhi[bestMesonPair[0]].id1());
    i1s.push_back(candsPhi[bestMesonPair[0]].id2());
    i2s.push_back(candsJPsi[bestMesonPair[1]].id1());
    i3s.push_back(candsJPsi[bestMesonPair[1]].id2());
    bxOffsetsFiller.addBx(bx, 1);
  }  // loop on BXs

  iEvent.put(std::move(ret), "selectedBx" + label);
  // now we make the table
  auto bxOffsets = bxOffsetsFiller.done();
  auto tab = std::make_unique<l1ScoutingRun3::OrbitFlatTable>(bxOffsets, "recMesonHphijpsi" + label, true);
  tab->addColumn<float>("mass", masses, "2 kaons + 2 muons invariant mass");
  tab->addColumn<uint8_t>("i0", i0s, "1st kaon (phi)");
  tab->addColumn<uint8_t>("i1", i1s, "2nd kaon (phi)");
  tab->addColumn<uint8_t>("i2", i2s, "1st muon (jpsi)");
  tab->addColumn<uint8_t>("i3", i3s, "2nd muon (jpsi)");
  iEvent.put(std::move(tab), "recMesonHphijpsi" + label);
}

float ScPhase2RecMesonHPhiJPsi::quadrimass(const l1Scouting::RecMeson *candsa, int a, const l1Scouting::RecMeson *candsb, int b) {
  ROOT::Math::PtEtaPhiMVector p1(candsa[a].pt(), candsa[a].eta(), candsa[a].phi(), candsa[a].mass());
  ROOT::Math::PtEtaPhiMVector p2(candsb[b].pt(), candsb[b].eta(), candsb[b].phi(), candsb[b].mass());
  float mass = (p1 + p2).M();
  return mass;
}

void ScPhase2RecMesonHPhiJPsi::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("srcPhi");
  desc.add<edm::InputTag>("srcJPsi");
  desc.add<bool>("runStruct", true);
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScPhase2RecMesonHPhiJPsi);