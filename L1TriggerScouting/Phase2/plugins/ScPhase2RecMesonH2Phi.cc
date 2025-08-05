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

class ScPhase2RecMesonH2Phi : public edm::stream::EDProducer<> {
public:
  explicit ScPhase2RecMesonH2Phi(const edm::ParameterSet &);
  ~ScPhase2RecMesonH2Phi() override;
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

  template <typename T>
  static float pairmass(const std::array<unsigned int, 2> &t, const T *cands, const std::array<float, 2> &massD);

  template <typename T>
  static float quadrupletmass(const std::array<unsigned int, 4> &t, const T *cands, const std::array<float, 4> &massD);

  unsigned long countStruct_;
  unsigned long passStruct_;
};

ScPhase2RecMesonH2Phi::ScPhase2RecMesonH2Phi(const edm::ParameterSet &iConfig)
    : doStruct_(iConfig.getParameter<bool>("runStruct")) {
  if (doStruct_) {
    structToken_ = consumes<OrbitCollection<l1Scouting::RecMeson>>(iConfig.getParameter<edm::InputTag>("src"));
    produces<std::vector<unsigned>>("selectedBx");
    produces<l1ScoutingRun3::OrbitFlatTable>("recMesonH2phi");
  }
}

ScPhase2RecMesonH2Phi::~ScPhase2RecMesonH2Phi() {};

void ScPhase2RecMesonH2Phi::beginStream(edm::StreamID) {
  countStruct_ = 0;
  passStruct_ = 0;
}

void ScPhase2RecMesonH2Phi::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  if (doStruct_) {
    edm::Handle<OrbitCollection<l1Scouting::RecMeson>> src;
    iEvent.getByToken(structToken_, src);

    runObj(*src, iEvent, countStruct_, passStruct_, "");
  }
}

void ScPhase2RecMesonH2Phi::endStream() {
  if (doStruct_)
    edm::LogImportant("ScPhase2AnalysisSummary") << "Rec Meson H2Phi Struct analysis: " << countStruct_ << " -> " << passStruct_;
}

template <typename T>
void ScPhase2RecMesonH2Phi::runObj(const OrbitCollection<T> &src,
                                    edm::Event &iEvent,
                                    unsigned long &nTry,
                                    unsigned long &nPass,
                                    const std::string &label) {
  l1ScoutingRun3::BxOffsetsFillter bxOffsetsFiller;
  bxOffsetsFiller.start();
  auto ret = std::make_unique<std::vector<unsigned>>();
  std::vector<float> masses;
  std::vector<uint8_t> i0s, i1s, i2s, i3s;
  ROOT::RVec<unsigned int> ix;  //
  std::array<unsigned int, 2> bestPair1, bestPair2;
  bool bestPair1Found, bestPair2Found;
  float bestPair1Score, bestPair2Score;

  for (unsigned int bx = 1; bx <= OrbitCollection<T>::NBX; ++bx) {
    nTry++;
    auto range = src.bxIterator(bx);
    const T *cands = &range.front();
    auto size = range.size();

    ix.clear();
    unsigned int ndaus = ix.size();
    if (ndaus < 4)
      continue;

    nPass++;
    masses.push_back(1.);
    i0s.push_back(1.); //bestQuadruplet[0]);
    i1s.push_back(1.); //bestQuadruplet[1]);
    i2s.push_back(1.); //bestQuadruplet[2]);
    i3s.push_back(1.); //bestQuadruplet[3]);
    bxOffsetsFiller.addBx(bx, 1);
  }  // loop on BXs

  iEvent.put(std::move(ret), "selectedBx" + label);
  // now we make the table
  auto bxOffsets = bxOffsetsFiller.done();
  auto tab = std::make_unique<l1ScoutingRun3::OrbitFlatTable>(bxOffsets, "recMesonH2phi" + label, true);
  tab->addColumn<float>("mass", masses, "4 kaons invariant mass");
  tab->addColumn<uint8_t>("i0", i0s, "1st kaon (phi1)");
  tab->addColumn<uint8_t>("i1", i1s, "2nd kaon (phi1)");
  tab->addColumn<uint8_t>("i2", i2s, "1st kaon (phi2)");
  tab->addColumn<uint8_t>("i3", i3s, "2nd kaon (phi2)");
  iEvent.put(std::move(tab), "recMesonH2phi" + label);
}

template <typename T>
float ScPhase2RecMesonH2Phi::pairmass(const std::array<unsigned int, 2> &t,
                                       const T *cands,
                                       const std::array<float, 2> &massD) {
  ROOT::Math::PtEtaPhiMVector p1(cands[t[0]].pt(), cands[t[0]].eta(), cands[t[0]].phi(), massD[0]);
  ROOT::Math::PtEtaPhiMVector p2(cands[t[1]].pt(), cands[t[1]].eta(), cands[t[1]].phi(), massD[1]);
  float mass = (p1 + p2).M();
  return mass;
}

template <typename T>
float ScPhase2RecMesonH2Phi::quadrupletmass(const std::array<unsigned int, 4> &t,
                                             const T *cands,
                                             const std::array<float, 4> &massD) {
  ROOT::Math::PtEtaPhiMVector p1(cands[t[0]].pt(), cands[t[0]].eta(), cands[t[0]].phi(), massD[0]);
  ROOT::Math::PtEtaPhiMVector p2(cands[t[1]].pt(), cands[t[1]].eta(), cands[t[1]].phi(), massD[1]);
  ROOT::Math::PtEtaPhiMVector p3(cands[t[2]].pt(), cands[t[2]].eta(), cands[t[2]].phi(), massD[2]);
  ROOT::Math::PtEtaPhiMVector p4(cands[t[3]].pt(), cands[t[3]].eta(), cands[t[3]].phi(), massD[3]);
  float mass = (p1 + p2 + p3 + p4).M();
  return mass;
}

void ScPhase2RecMesonH2Phi::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src");
  desc.add<bool>("runStruct", true);
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScPhase2RecMesonH2Phi);
