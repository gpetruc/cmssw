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

#include "DataFormats/Math/interface/deltaR.h"
#include <ROOT/RVec.hxx>
#include <Math/Vector4D.h>
#include <Math/GenVector/LorentzVector.h>
#include <Math/GenVector/PtEtaPhiM4D.h>
#include <algorithm>
#include <array>
#include <iostream>

class ScPhase2PuppiH2RhoDemo : public edm::stream::EDProducer<> {
public:
  explicit ScPhase2PuppiH2RhoDemo(const edm::ParameterSet &);
  ~ScPhase2PuppiH2RhoDemo() override;
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
  edm::EDGetTokenT<OrbitCollection<l1Scouting::Puppi>> structToken_;

  struct Cuts {
    float minptD = 5;
    float minptQ = 3;
    float maxdeltarD2 = 0.40 * 0.40;
    float minmassH = 100;
    float maxmassH = 150;
    float minmassQ = 0.40;
    float maxmassQ = 1.30;
    float mindr2 = 0.05 * 0.05;
    float maxdr2 = 0.25 * 0.25;
    float maxiso = 0.25;
  } cuts;

  template <typename T>
  std::tuple<bool, float> isolationQ(
      float eta, float phi, unsigned int pidex1, unsigned int pidex2, const T *cands, unsigned int size) const;

  std::tuple<bool, float> deltar(float eta1, float eta2, float phi1, float phi2) const;

  template <typename T>
  static ROOT::Math::PtEtaPhiMVector pairP4(const std::array<unsigned int, 2> &t,
                                            const T *cands,
                                            const std::array<float, 2> &massD);

  unsigned long countStruct_;
  unsigned long passStruct_;
};

ScPhase2PuppiH2RhoDemo::ScPhase2PuppiH2RhoDemo(const edm::ParameterSet &iConfig)
    : doStruct_(iConfig.getParameter<bool>("runStruct")) {
  if (doStruct_) {
    structToken_ = consumes<OrbitCollection<l1Scouting::Puppi>>(iConfig.getParameter<edm::InputTag>("src"));
    produces<std::vector<unsigned>>("selectedBx");
    produces<l1ScoutingRun3::OrbitFlatTable>("h2rho");
  }
}

ScPhase2PuppiH2RhoDemo::~ScPhase2PuppiH2RhoDemo() {};

void ScPhase2PuppiH2RhoDemo::beginStream(edm::StreamID) {
  countStruct_ = 0;
  passStruct_ = 0;
}

void ScPhase2PuppiH2RhoDemo::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  if (doStruct_) {
    edm::Handle<OrbitCollection<l1Scouting::Puppi>> src;
    iEvent.getByToken(structToken_, src);

    runObj(*src, iEvent, countStruct_, passStruct_, "");
  }
}

void ScPhase2PuppiH2RhoDemo::endStream() {
  if (doStruct_)
    edm::LogImportant("ScPhase2AnalysisSummary") << "H2Rho Struct analysis: " << countStruct_ << " -> " << passStruct_;
}

template <typename T>
void ScPhase2PuppiH2RhoDemo::runObj(const OrbitCollection<T> &src,
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
  ROOT::Math::PtEtaPhiMVector p4Q1, p4Q2;
  for (unsigned int bx = 1; bx <= OrbitCollection<T>::NBX; ++bx) {
    nTry++;
    auto range = src.bxIterator(bx);
    auto size = range.size();
    const T *cands = (size > 0) ? &range.front() : nullptr;

    ix.clear();
    for (unsigned int i = 0; i < size; ++i) {  //make list of all hadrons
      if ((std::abs(cands[i].pdgId()) == 211 or std::abs(cands[i].pdgId()) == 11)) {
        if (cands[i].pt() >= cuts.minptD)
          ix.push_back(i);
      }
    }
    unsigned int ndaus = ix.size();
    if (ndaus < 4)
      continue;

    // Q1 candidate from highest pT OS pair with mass compatible with mQ
    bestPair1Found = false;
    bestPair1Score = cuts.minptQ;
    for (unsigned int i1 = 0; i1 < ndaus; ++i1) {
      for (unsigned int i2 = i1 + 1; i2 < ndaus; ++i2) {
        if (!(cands[ix[i1]].charge() * cands[ix[i2]].charge() < 0))
          continue;

        auto p4 = pairP4({{ix[i1], ix[i2]}}, cands, {{0.1396, 0.1396}});
        auto mass2 = p4.M();
        if (!(mass2 >= cuts.minmassQ and mass2 <= cuts.maxmassQ))
          continue;

        auto [drcond, dr2Q] =
            deltar(cands[ix[i1]].eta(), cands[ix[i2]].eta(), cands[ix[i1]].phi(), cands[ix[i2]].phi());
        if (!drcond)
          continue;  // angular sep of top 2 tracks

        auto [isolated, psum] = isolationQ(p4.eta(), p4.phi(), ix[i1], ix[i2], cands, size);
        if (!isolated)
          continue;  // Q isolation

        std::array<unsigned int, 2> pair{{ix[i1], ix[i2]}};  // pair of indices
        if (p4.pt() > bestPair1Score) {
          std::copy_n(pair.begin(), 2, bestPair1.begin());
          bestPair1Found = true;
          bestPair1Score = p4.pt();
          p4Q1 = p4;
        }
      }
    }
    if (!bestPair1Found)
      continue;  // pair was found

    // Q2 candidate from highest pT OS pair with mass compatible with mQ
    bestPair2Found = false;
    bestPair2Score = cuts.minptQ;
    for (unsigned int i3 = 0; i3 < ndaus; ++i3) {
      if (ix[i3] == bestPair1[0] or ix[i3] == bestPair1[1])
        continue;  // don't reuse candidates from previous pair
      for (unsigned int i4 = i3 + 1; i4 < ndaus; ++i4) {
        if (ix[i4] == bestPair1[0] or ix[i4] == bestPair1[1])
          continue;  // don't reuse candidates from previous pair
        if (!(cands[ix[i3]].charge() * cands[ix[i4]].charge() < 0))
          continue;  // OS pair
        auto p4 = pairP4({{ix[i3], ix[i4]}}, cands, {{0.1396, 0.1396}});
        auto mass2 = p4.M();
        if (!(mass2 >= cuts.minmassQ and mass2 <= cuts.maxmassQ))
          continue;  // Q mass
        auto [drcond, dr2Q] =
            deltar(cands[ix[i3]].eta(), cands[ix[i4]].eta(), cands[ix[i3]].phi(), cands[ix[i4]].phi());
        if (!drcond)
          continue;  // angular sep of top 2 tracks

        auto [isolated, psum] = isolationQ(p4.eta(), p4.phi(), ix[i3], ix[i4], cands, size);
        if (!isolated)
          continue;  // Q isolation

        std::array<unsigned int, 2> pair{{ix[i3], ix[i4]}};  // pair of indices
        if (p4.pt() > bestPair2Score) {
          std::copy_n(pair.begin(), 2, bestPair2.begin());
          bestPair2Found = true;
          bestPair2Score = p4.pt();
          p4Q2 = p4;
        }
      }
    }
    if (!bestPair2Found)
      continue;  // pair was found

    std::array<unsigned int, 4> bestQuadruplet{{bestPair1[0], bestPair1[1], bestPair2[0], bestPair2[1]}};
    // H mass
    auto mass = (p4Q1 + p4Q2).M();
    if (!(mass >= cuts.minmassH and mass <= cuts.maxmassH))
      continue;

    ret->emplace_back(bx);
    nPass++;
    masses.push_back(mass);
    i0s.push_back(bestQuadruplet[0]);
    i1s.push_back(bestQuadruplet[1]);
    i2s.push_back(bestQuadruplet[2]);
    i3s.push_back(bestQuadruplet[3]);
    bxOffsetsFiller.addBx(bx, 1);
  }  // loop on BXs

  iEvent.put(std::move(ret), "selectedBx" + label);
  // now we make the table
  auto bxOffsets = bxOffsetsFiller.done();
  auto tab = std::make_unique<l1ScoutingRun3::OrbitFlatTable>(bxOffsets, "H2Rho" + label, true);
  tab->addColumn<float>("mass", masses, "4 pions invariant mass");
  tab->addColumn<uint8_t>("i0", i0s, "1st pion (rho1)");
  tab->addColumn<uint8_t>("i1", i1s, "2nd pion (rho1)");
  tab->addColumn<uint8_t>("i2", i2s, "1st pion (rho2)");
  tab->addColumn<uint8_t>("i3", i3s, "2nd pion (rho2)");
  iEvent.put(std::move(tab), "h2rho" + label);
}

//TEST functions
template <typename T>
std::tuple<bool, float> ScPhase2PuppiH2RhoDemo::isolationQ(
    float eta, float phi, unsigned int pidex1, unsigned int pidex2, const T *cands, unsigned int size) const {
  float psum = 0;
  for (unsigned int j = 0u; j < size; ++j) {  //loop over other particles
    if (pidex1 == j or pidex2 == j)
      continue;
    float dr2 = reco::deltaR2(eta, phi, cands[j].eta(), cands[j].phi());
    if (dr2 >= cuts.mindr2 && dr2 <= cuts.maxdr2)
      psum += cands[j].pt();
  }
  bool passed = (psum <= cuts.maxiso * (cands[pidex1].pt() + cands[pidex2].pt()));
  return std::tuple(passed, psum);
}

std::tuple<bool, float> ScPhase2PuppiH2RhoDemo::deltar(float eta1, float eta2, float phi1, float phi2) const {
  float dr2 = reco::deltaR2(eta1, phi1, eta2, phi2);
  bool passed = (dr2 <= cuts.maxdeltarD2);
  return std::tuple(passed, dr2);
}

template <typename T>
ROOT::Math::PtEtaPhiMVector ScPhase2PuppiH2RhoDemo::pairP4(const std::array<unsigned int, 2> &t,
                                                           const T *cands,
                                                           const std::array<float, 2> &massD) {
  ROOT::Math::PtEtaPhiMVector p1(cands[t[0]].pt(), cands[t[0]].eta(), cands[t[0]].phi(), massD[0]);
  ROOT::Math::PtEtaPhiMVector p2(cands[t[1]].pt(), cands[t[1]].eta(), cands[t[1]].phi(), massD[1]);
  return p1 + p2;
}

void ScPhase2PuppiH2RhoDemo::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src");
  desc.add<bool>("runStruct", true);
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScPhase2PuppiH2RhoDemo);
