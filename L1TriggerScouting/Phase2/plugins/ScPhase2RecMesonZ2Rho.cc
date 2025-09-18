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

class ScPhase2RecMesonZ2Rho : public edm::stream::EDProducer<> {
public:
  explicit ScPhase2RecMesonZ2Rho(const edm::ParameterSet &);
  ~ScPhase2RecMesonZ2Rho() override;
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

  bool doStruct_;
  edm::EDGetTokenT<OrbitCollection<l1Scouting::RecMeson>> structToken_;

  struct Cuts {
    float minmassZ = 60;
    float maxmassZ = 120;
    float minptQ = 3;
    float maxiso = 0.25;
  } cuts;

  template <typename T>
  bool isolationQ(unsigned int pidex1, unsigned int pidex2, const T *cands, unsigned int size) const;

  std::tuple<bool, float> deltar(float eta1, float eta2, float phi1, float phi2) const;

  static float pairmass(const l1Scouting::RecMeson *cands, int a, int b);

  unsigned long countStruct_;
  unsigned long passStruct_;
};

ScPhase2RecMesonZ2Rho::ScPhase2RecMesonZ2Rho(const edm::ParameterSet &iConfig)
    : doStruct_(iConfig.getParameter<bool>("runStruct")) {
  if (doStruct_) {
    structToken_ = consumes<OrbitCollection<l1Scouting::RecMeson>>(iConfig.getParameter<edm::InputTag>("src"));
    produces<std::vector<unsigned>>("selectedBx");
    produces<l1ScoutingRun3::OrbitFlatTable>("recMesonZ2rho");
  }
}

ScPhase2RecMesonZ2Rho::~ScPhase2RecMesonZ2Rho() {};

void ScPhase2RecMesonZ2Rho::beginStream(edm::StreamID) {
  countStruct_ = 0;
  passStruct_ = 0;
}

void ScPhase2RecMesonZ2Rho::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  if (doStruct_) {
    edm::Handle<OrbitCollection<l1Scouting::RecMeson>> src;
    iEvent.getByToken(structToken_, src);

    runObj(*src, iEvent, countStruct_, passStruct_, "");
  }
}

void ScPhase2RecMesonZ2Rho::endStream() {
  if (doStruct_)
    edm::LogImportant("ScPhase2AnalysisSummary") << "Rec Meson Z2Rho Struct analysis: " << countStruct_ << " -> " << passStruct_;
}

template <typename T>
void ScPhase2RecMesonZ2Rho::runObj(const OrbitCollection<T> &src,
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
    auto range = src.bxIterator(bx);
    const T *cands = &range.front();
    auto size = range.size();
    unsigned int ndaus = size;

    bestMesonPairScore = 0.;
    bestMesonPairFound = false;

    if (ndaus < 2)
      continue;

    for (unsigned int i1 = 0; i1 < ndaus; ++i1) {
      // minimum pt and isolation of Q1
      if ((cands[i1].pt() < cuts.minptQ) || (cands[i1].isoDR0p25() >= cuts.maxiso))
        continue;
      for (unsigned int i2 = i1 + 1; i2 < ndaus; ++i2) {
        // minimum pt and isolation of Q2
        if ((cands[i2].pt() < cuts.minptQ) || (cands[i2].isoDR0p25() >= cuts.maxiso))
          continue;

        // Four different dauther particles
        if ((cands[i1].id1() == cands[i2].id1()) || (cands[i1].id1() == cands[i2].id2()))
          continue;
        if ((cands[i1].id2() == cands[i2].id1()) || (cands[i1].id2() == cands[i2].id2()))
          continue;

        auto p4_1 = ROOT::Math::PtEtaPhiMVector(cands[i1].pt(), cands[i1].eta(), cands[i1].phi(), cands[i1].mass());
        auto p4_2 = ROOT::Math::PtEtaPhiMVector(cands[i2].pt(), cands[i2].eta(), cands[i2].phi(), cands[i2].mass());

        // Choose best pair of mesons based on score (e.g. max pt)
        float ptsum = (p4_1 + p4_2).pt();

        if (ptsum > bestMesonPairScore) {
          bestMesonPairScore = ptsum;
          bestMesonPair = {{i1, i2}};
          bestMesonPairFound = true;
        }
      }
    }

    if (!bestMesonPairFound)
      continue;

    // Z mass
    auto mass = pairmass(cands, bestMesonPair[0], bestMesonPair[1]);
    if (!(mass >= cuts.minmassZ and mass <= cuts.maxmassZ))
      continue;

    ret->emplace_back(bx);

    nPass++;
    masses.push_back(mass);
    i0s.push_back(cands[bestMesonPair[0]].id1());
    i1s.push_back(cands[bestMesonPair[0]].id2());
    i2s.push_back(cands[bestMesonPair[1]].id1());
    i3s.push_back(cands[bestMesonPair[1]].id2());
    bxOffsetsFiller.addBx(bx, 1);
  }  // loop on BXs

  iEvent.put(std::move(ret), "selectedBx" + label);
  // now we make the table
  auto bxOffsets = bxOffsetsFiller.done();
  auto tab = std::make_unique<l1ScoutingRun3::OrbitFlatTable>(bxOffsets, "recMesonZ2rho" + label, true);
  tab->addColumn<float>("mass", masses, "4 pions invariant mass");
  tab->addColumn<uint8_t>("i0", i0s, "1st pion (rho1)");
  tab->addColumn<uint8_t>("i1", i1s, "2nd pion (rho1)");
  tab->addColumn<uint8_t>("i2", i2s, "1st pion (rho2)");
  tab->addColumn<uint8_t>("i3", i3s, "2nd pion (rho2)");
  iEvent.put(std::move(tab), "recMesonZ2rho" + label);
}

float ScPhase2RecMesonZ2Rho::pairmass(const l1Scouting::RecMeson *cands, int a, int b) {
  ROOT::Math::PtEtaPhiMVector p1(cands[a].pt(), cands[a].eta(), cands[a].phi(), cands[a].mass());
  ROOT::Math::PtEtaPhiMVector p2(cands[b].pt(), cands[b].eta(), cands[b].phi(), cands[b].mass());
  float mass = (p1 + p2).M();
  return mass;
}

void ScPhase2RecMesonZ2Rho::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src");
  desc.add<bool>("runStruct", true);
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScPhase2RecMesonZ2Rho);
