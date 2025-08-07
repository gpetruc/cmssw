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

class ScPhase2RecMesonH2Rho : public edm::stream::EDProducer<> {
public:
  explicit ScPhase2RecMesonH2Rho(const edm::ParameterSet &);
  ~ScPhase2RecMesonH2Rho() override;
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
    float minmassH = 100;
    float maxmassH = 150;
  } cuts;

  template <typename T>
  bool isolationQ(unsigned int pidex1, unsigned int pidex2, const T *cands, unsigned int size) const;

  std::tuple<bool, float> deltar(float eta1, float eta2, float phi1, float phi2) const;

  static float quadrupletmass(const l1Scouting::RecMeson *cands);

  unsigned long countStruct_;
  unsigned long passStruct_;
};

ScPhase2RecMesonH2Rho::ScPhase2RecMesonH2Rho(const edm::ParameterSet &iConfig)
    : doStruct_(iConfig.getParameter<bool>("runStruct")) {
  if (doStruct_) {
    structToken_ = consumes<OrbitCollection<l1Scouting::RecMeson>>(iConfig.getParameter<edm::InputTag>("src"));
    produces<std::vector<unsigned>>("selectedBx");
    produces<l1ScoutingRun3::OrbitFlatTable>("recMesonH2rho");
  }
}

ScPhase2RecMesonH2Rho::~ScPhase2RecMesonH2Rho() {};

void ScPhase2RecMesonH2Rho::beginStream(edm::StreamID) {
  countStruct_ = 0;
  passStruct_ = 0;
}

void ScPhase2RecMesonH2Rho::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  if (doStruct_) {
    edm::Handle<OrbitCollection<l1Scouting::RecMeson>> src;
    iEvent.getByToken(structToken_, src);

    runObj(*src, iEvent, countStruct_, passStruct_, "");
  }
}

void ScPhase2RecMesonH2Rho::endStream() {
  if (doStruct_)
    edm::LogImportant("ScPhase2AnalysisSummary") << "Rec Meson H2Rho Struct analysis: " << countStruct_ << " -> " << passStruct_;
}

template <typename T>
void ScPhase2RecMesonH2Rho::runObj(const OrbitCollection<T> &src,
                                    edm::Event &iEvent,
                                    unsigned long &nTry,
                                    unsigned long &nPass,
                                    const std::string &label) {
  l1ScoutingRun3::BxOffsetsFillter bxOffsetsFiller;
  bxOffsetsFiller.start();
  auto ret = std::make_unique<std::vector<unsigned>>();
  std::vector<float> masses;
  std::vector<uint8_t> i0s, i1s, i2s, i3s;
  std::array<unsigned int, 2> bestPair1, bestPair2;
  bool bestPair1Found, bestPair2Found;
  float bestPair1Score, bestPair2Score;

  for (unsigned int bx = 1; bx <= OrbitCollection<T>::NBX; ++bx) {
    nTry++;
    auto range = src.bxIterator(bx);
    const T *cands = &range.front();
    auto size = range.size();
    unsigned int ndaus = size;

    if ( ndaus >= 2) {
      // std::cout << "NEW" << std::endl;
      // std::cout << "BX = " << bx << " ; number of mesons = " << ndaus << std::endl;
      // std::cout << "cand 0 - ids = " << cands[0].id1() << " and " << cands[0].id2() << std::endl ;
      // std::cout << "cand 1 - ids = " << cands[1].id1() << " and " << cands[1].id2() << std::endl << std::endl;
    }

    if (ndaus < 2)
      continue;

    // H mass
    auto mass = quadrupletmass(cands);
    if (!(mass >= cuts.minmassH and mass <= cuts.maxmassH))
      continue;

    ret->emplace_back(bx);

    nPass++;
    masses.push_back(mass);
    i0s.push_back(cands[0].id1());
    i1s.push_back(cands[0].id2());
    i2s.push_back(cands[1].id1());
    i3s.push_back(cands[1].id2());
    bxOffsetsFiller.addBx(bx, 1);
  }  // loop on BXs

  iEvent.put(std::move(ret), "selectedBx" + label);
  // now we make the table
  auto bxOffsets = bxOffsetsFiller.done();
  auto tab = std::make_unique<l1ScoutingRun3::OrbitFlatTable>(bxOffsets, "recMesonH2rho" + label, true);
  tab->addColumn<float>("mass", masses, "4 kaons invariant mass");
  tab->addColumn<uint8_t>("i0", i0s, "1st kaon (rho1)");
  tab->addColumn<uint8_t>("i1", i1s, "2nd kaon (rho1)");
  tab->addColumn<uint8_t>("i2", i2s, "1st kaon (rho2)");
  tab->addColumn<uint8_t>("i3", i3s, "2nd kaon (rho2)");
  iEvent.put(std::move(tab), "recMesonH2rho" + label);
}

float ScPhase2RecMesonH2Rho::quadrupletmass(const l1Scouting::RecMeson *cands) {
  ROOT::Math::PtEtaPhiMVector p1(cands[0].pt(), cands[0].eta(), cands[0].phi(), cands[0].mass());
  ROOT::Math::PtEtaPhiMVector p2(cands[1].pt(), cands[1].eta(), cands[1].phi(), cands[1].mass());
  ROOT::Math::PtEtaPhiMVector p3(cands[2].pt(), cands[2].eta(), cands[2].phi(), cands[2].mass());
  ROOT::Math::PtEtaPhiMVector p4(cands[3].pt(), cands[3].eta(), cands[3].phi(), cands[3].mass());
  float mass = (p1 + p2 + p3 + p4).M();
  return mass;
}

void ScPhase2RecMesonH2Rho::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src");
  desc.add<bool>("runStruct", true);
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScPhase2RecMesonH2Rho);
