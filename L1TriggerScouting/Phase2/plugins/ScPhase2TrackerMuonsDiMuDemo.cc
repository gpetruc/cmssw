#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/L1Scouting/interface/OrbitCollection.h"
#include "DataFormats/L1Scouting/interface/OrbitFlatTable.h"
#include "DataFormats/L1TParticleFlow/interface/PFCandidate.h"
#include "DataFormats/L1TMuonPhase2/interface/L1ScoutingTrackerMuon.h"
#include "L1TriggerScouting/Utilities/interface/BxOffsetsFiller.h"

#include <ROOT/RVec.hxx>
#include <Math/Vector4D.h>
#include <Math/GenVector/LorentzVector.h>
#include <Math/GenVector/PtEtaPhiM4D.h>
#include <algorithm>
#include <array>
#include <iostream>

class ScPhase2TrackerMuonDiMuDemo : public edm::stream::EDProducer<> {
public:
  explicit ScPhase2TrackerMuonDiMuDemo(const edm::ParameterSet &);
  ~ScPhase2TrackerMuonDiMuDemo() override;
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  void beginStream(edm::StreamID) override;
  void produce(edm::Event &, const edm::EventSetup &) override;
  void endStream() override;
  template <typename T>
  void runObj(const OrbitCollection<T> &src,
              edm::Event &out,
              unsigned long &nTry,
              unsigned long &nPass,
              const std::string &bxLabel);

  bool doCandidate_, doStruct_;
  edm::EDGetTokenT<OrbitCollection<l1Scouting::TrackerMuon>> structToken_;

  struct Cuts {
    float minptOverMass = 0.25;
    float qualityMin = 75;
    float massMin = 7.0;  //  after b's
    bool doOppsiteCharge = true;
    float etaMax = 2.0;
    float ptMin = 2.0;
    float z0Max = 1.0;
  } cuts;

  unsigned long countStruct_;
  unsigned long passStruct_;
};

ScPhase2TrackerMuonDiMuDemo::ScPhase2TrackerMuonDiMuDemo(const edm::ParameterSet &iConfig)
    : doStruct_(iConfig.getParameter<bool>("runStruct")) {
  if (doStruct_) {
    structToken_ = consumes<OrbitCollection<l1Scouting::TrackerMuon>>(iConfig.getParameter<edm::InputTag>("src"));
    produces<std::vector<unsigned>>("selectedBx");
    produces<l1ScoutingRun3::OrbitFlatTable>("DiMu");
  }
}

ScPhase2TrackerMuonDiMuDemo::~ScPhase2TrackerMuonDiMuDemo(){};

void ScPhase2TrackerMuonDiMuDemo::beginStream(edm::StreamID) {
  countStruct_ = 0;
  passStruct_ = 0;
}

void ScPhase2TrackerMuonDiMuDemo::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  if (doStruct_) {
    edm::Handle<OrbitCollection<l1Scouting::TrackerMuon>> src;
    iEvent.getByToken(structToken_, src);
    runObj(*src, iEvent, countStruct_, passStruct_, "");
  }
}

void ScPhase2TrackerMuonDiMuDemo::endStream() {
  if (doStruct_)
    std::cout << "DiTrackerMuon | Struct analysis : " << countStruct_ << " -> " << passStruct_ << std::endl;
}

template <typename T>
void ScPhase2TrackerMuonDiMuDemo::runObj(const OrbitCollection<T> &src,
                                         edm::Event &iEvent,
                                         unsigned long &nTry,
                                         unsigned long &nPass,
                                         const std::string &label) {
  l1ScoutingRun3::BxOffsetsFillter bxOffsetsFiller;
  bxOffsetsFiller.start();
  auto selectedBx_idx = std::make_unique<std::vector<unsigned>>();
  float mass;
  std::vector<float> masses;
  std::vector<uint8_t> i0s, i1s;
  for (unsigned int bx = 1; bx <= OrbitCollection<T>::NBX; ++bx) {
    nTry++;
    auto range = src.bxIterator(bx);
    const T *cands = &range.front();
    auto size = range.size();
    bool hasCandidate = false;
    for (unsigned int i = 0; i < size; ++i) {
      if (cands[i].pt() < cuts.ptMin)
        break;  // assumes pt ordering
      if (cands[i].z0() > cuts.z0Max)
        continue;
      if (abs(cands[i].eta()) > cuts.etaMax)
        continue;
      if (cands[i].quality() < cuts.qualityMin)
        continue;
      // if(    cands[i].isolation()   > XX  ) continue;

      for (unsigned int j = i + 1; j < size; ++j) {
        if (cands[j].pt() < cuts.ptMin)
          break;  // assumes pt ordering
        if (cuts.doOppsiteCharge and (cands[i].charge() * cands[j].charge() > 0))
          continue;
        if (abs(cands[j].z0()) > cuts.z0Max)
          continue;
        if (abs(cands[j].eta()) > cuts.etaMax)
          continue;
        if (cands[j].quality() < cuts.qualityMin)
          continue;
        //      if(    cands[j].isolation()   > YY  ) continue;

        mass = (cands[i].p4() + cands[j].p4()).mass();
        if (mass < cuts.massMin)
          continue;
        if (cands[j].pt() < cuts.minptOverMass * mass)
          continue;
        if (cands[i].pt() < cuts.minptOverMass * mass)
          continue;
        selectedBx_idx->emplace_back(bx);
        masses.push_back(mass);
        i0s.push_back(i);
        i1s.push_back(j);
        bxOffsetsFiller.addBx(bx, 1);
        hasCandidate = true;
      }
    }

    if (hasCandidate)
      nPass++;
  }  // loop on BXs

  iEvent.put(std::move(selectedBx_idx), "selectedBx" + label);
  // now we make the table
  auto bxOffsets = bxOffsetsFiller.done();
  auto tab = std::make_unique<l1ScoutingRun3::OrbitFlatTable>(bxOffsets, "DiMu" + label, true);
  tab->addColumn<float>("mass", masses, "Dimuon invariant mass");
  tab->addColumn<uint8_t>("i0", i0s, "leading muon");
  tab->addColumn<uint8_t>("i1", i1s, "subleading muon");
  iEvent.put(std::move(tab), "DiMu" + label);
}

void ScPhase2TrackerMuonDiMuDemo::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScPhase2TrackerMuonDiMuDemo);
