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

class ScPhase2PuppiH2PhiDemo : public edm::stream::EDProducer<> {
public:
  explicit ScPhase2PuppiH2PhiDemo(const edm::ParameterSet &);
  ~ScPhase2PuppiH2PhiDemo() override;
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
    float minptD = 10;
    float minptQ = 30;
    float maxdeltarD2 = 0.40 * 0.40;
    float minmassH = 100;
    float maxmassH = 150;
    float minmassQ = 0.95;
    float maxmassQ = 1.25;
    float mindr2 = 0.05 * 0.05;
    float maxdr2 = 0.25 * 0.25;
    float maxiso = 0.25;
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

ScPhase2PuppiH2PhiDemo::ScPhase2PuppiH2PhiDemo(const edm::ParameterSet &iConfig)
    : doStruct_(iConfig.getParameter<bool>("runStruct")) {
  if (doStruct_) {
    structToken_ = consumes<OrbitCollection<l1Scouting::Puppi>>(iConfig.getParameter<edm::InputTag>("src"));
    produces<std::vector<unsigned>>("selectedBx");
    produces<l1ScoutingRun3::OrbitFlatTable>("h2phi");
  }
}

ScPhase2PuppiH2PhiDemo::~ScPhase2PuppiH2PhiDemo() {};

void ScPhase2PuppiH2PhiDemo::beginStream(edm::StreamID) {
  countStruct_ = 0;
  passStruct_ = 0;
}

void ScPhase2PuppiH2PhiDemo::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  if (doStruct_) {
    edm::Handle<OrbitCollection<l1Scouting::Puppi>> src;
    iEvent.getByToken(structToken_, src);

    runObj(*src, iEvent, countStruct_, passStruct_, "");
  }
}

void ScPhase2PuppiH2PhiDemo::endStream() {
  if (doStruct_)
    edm::LogImportant("ScPhase2AnalysisSummary") << "H2Phi Struct analysis: " << countStruct_ << " -> " << passStruct_;
}

template <typename T>
void ScPhase2PuppiH2PhiDemo::runObj(const OrbitCollection<T> &src,
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
    for (unsigned int i = 0; i < size; ++i) {  //make list of all hadrons
      if ((std::abs(cands[i].pdgId()) == 211 or std::abs(cands[i].pdgId()) == 11)) {
        if (cands[i].pt() >= cuts.minptD)
          ix.push_back(i);
      }
    }
    unsigned int ndaus = ix.size();
    if (ndaus < 4)
      continue;

    // Q1 candidate from closest OS pair with mass compatible with mQ
    bestPair1Found = false;
    bestPair1Score = 999;
    for (unsigned int i1 = 0; i1 < ndaus; ++i1) {
      if (cands[ix[i1]].pt() < cuts.minptD)
        continue;  // D1 pt cut
      for (unsigned int i2 = 0; i2 < ndaus; ++i2) {
        if (i2 == i1 || cands[ix[i2]].pt() < cuts.minptD)
          continue;  // D2 pt cut

        if (!(cands[ix[i1]].charge() * cands[ix[i2]].charge() < 0))
          continue;

        auto mass2 = pairmass({{ix[i1], ix[i2]}}, cands, {{0.4937, 0.4937}});
        if (mass2 >= cuts.minmassQ and mass2 <= cuts.maxmassQ)
          continue;

        auto [drcond, drQ] = deltar(cands[ix[i1]].eta(), cands[ix[i2]].eta(), cands[ix[i1]].phi(), cands[ix[i2]].phi());
        if (!drcond)
          continue;  // angular sep of top 2 tracks

        std::array<unsigned int, 2> pair{{ix[i1], ix[i2]}};  // pair of indices
        if (drQ < bestPair1Score) {
          std::copy_n(pair.begin(), 2, bestPair1.begin());
          bestPair1Score = drQ;
          if (bestPair1Score * bestPair1Score < cuts.maxdeltarD2)
            bestPair1Found = true;
        }
      }
    }
    if (!bestPair1Found)
      continue;  // pair was found
    auto ptQ = (cands[bestPair1[0]].p4() + cands[bestPair1[1]].p4()).pt();
    if (ptQ < cuts.minptQ)
      continue;  // Q pt
    if (!isolationQ(bestPair1[0], bestPair1[1], cands, size))
      continue;  // Q isolation

    // Q2 candidate from closest OS pair with mass compatible with mQ
    bestPair2Found = false;
    bestPair2Score = 999;
    for (unsigned int i3 = 0; i3 < ndaus; ++i3) {
      if (cands[ix[i3]].pt() < cuts.minptD)
        continue;  // D1 pt cut
      if (ix[i3] == bestPair1[0] or ix[i3] == bestPair1[1])
        continue;  // don't reuse candidates from previous pair
      for (unsigned int i4 = 0; i4 < ndaus; ++i4) {
        if (i4 == i3 || cands[ix[i4]].pt() < cuts.minptD)
          continue;  // D2 pt cut
        if (ix[i4] == bestPair1[0] or ix[i4] == bestPair1[1])
          continue;  // don't reuse candidates from previous pair
        if (!(cands[ix[i3]].charge() * cands[ix[i4]].charge() < 0))
          continue;  // OS pair
        auto mass2 = pairmass(
            {{ix[i3], ix[i4]}}, cands, {{0.4937, 0.4937}});  // (cands[ix[i3]].p4() + cands[ix[i4]].p4()).mass();
        if (mass2 >= cuts.minmassQ and mass2 <= cuts.maxmassQ)
          continue;  // Q mass
        auto [drcond, drQ] = deltar(cands[ix[i3]].eta(), cands[ix[i4]].eta(), cands[ix[i3]].phi(), cands[ix[i4]].phi());
        if (!drcond)
          continue;  // angular sep of top 2 tracks

        std::array<unsigned int, 2> pair{{ix[i3], ix[i4]}};  // pair of indices
        if (drQ < bestPair2Score) {
          std::copy_n(pair.begin(), 2, bestPair2.begin());
          bestPair2Score = drQ;
          if (bestPair2Score * bestPair2Score < cuts.maxdeltarD2)
            bestPair2Found = true;
        }
      }
    }
    if (!bestPair2Found)
      continue;  // pair was found
    ptQ = (cands[bestPair2[0]].p4() + cands[bestPair2[1]].p4()).pt();
    if (ptQ < cuts.minptQ)
      continue;  // Q pt
    if (!isolationQ(bestPair2[0], bestPair2[1], cands, size))
      continue;  // Q isolation

    std::array<unsigned int, 4> bestQuadruplet{{bestPair1[0], bestPair1[1], bestPair2[0], bestPair2[1]}};
    // H mass
    auto mass = quadrupletmass(bestQuadruplet, cands, {{0.4937, 0.4937, 0.4937, 0.4937}});
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
  auto tab = std::make_unique<l1ScoutingRun3::OrbitFlatTable>(bxOffsets, "H2Phi" + label, true);
  tab->addColumn<float>("mass", masses, "4 kaons invariant mass");
  tab->addColumn<uint8_t>("i0", i0s, "1st kaon (phi1)");
  tab->addColumn<uint8_t>("i1", i1s, "2nd kaon (phi1)");
  tab->addColumn<uint8_t>("i2", i2s, "1st kaon (phi2)");
  tab->addColumn<uint8_t>("i3", i3s, "2nd kaon (phi2)");
  iEvent.put(std::move(tab), "h2phi" + label);
}

//TEST functions
template <typename T>
bool ScPhase2PuppiH2PhiDemo::isolationQ(unsigned int pidex1,
                                        unsigned int pidex2,
                                        const T *cands,
                                        unsigned int size) const {
  bool passed = false;
  float psum = 0;
  float eta = cands[pidex1].eta();  //center cone around leading track
  float phi = cands[pidex1].phi();
  for (unsigned int j = 0u; j < size; ++j) {  //loop over other particles
    if (pidex1 == j or pidex2 == j)
      continue;
    float deta = eta - cands[j].eta(), dphi = ROOT::VecOps::DeltaPhi<float>(phi, cands[j].phi());
    float dr2 = deta * deta + dphi * dphi;
    if (dr2 >= cuts.mindr2 && dr2 <= cuts.maxdr2)
      psum += cands[j].pt();
  }
  if (psum <= cuts.maxiso * (cands[pidex1].pt() + cands[pidex2].pt()))
    passed = true;
  return passed;
}

std::tuple<bool, float> ScPhase2PuppiH2PhiDemo::deltar(float eta1, float eta2, float phi1, float phi2) const {
  bool passed = true;
  float deta = eta1 - eta2;
  float dphi = ROOT::VecOps::DeltaPhi<float>(phi1, phi2);
  float dr2 = deta * deta + dphi * dphi;
  if (dr2 > cuts.maxdeltarD2) {
    passed = false;
    return std::tuple(passed, dr2);
  }
  return std::tuple(passed, dr2);
}

template <typename T>
float ScPhase2PuppiH2PhiDemo::pairmass(const std::array<unsigned int, 2> &t,
                                       const T *cands,
                                       const std::array<float, 2> &massD) {
  ROOT::Math::PtEtaPhiMVector p1(cands[t[0]].pt(), cands[t[0]].eta(), cands[t[0]].phi(), massD[0]);
  ROOT::Math::PtEtaPhiMVector p2(cands[t[1]].pt(), cands[t[1]].eta(), cands[t[1]].phi(), massD[1]);
  float mass = (p1 + p2).M();
  return mass;
}

template <typename T>
float ScPhase2PuppiH2PhiDemo::quadrupletmass(const std::array<unsigned int, 4> &t,
                                             const T *cands,
                                             const std::array<float, 4> &massD) {
  ROOT::Math::PtEtaPhiMVector p1(cands[t[0]].pt(), cands[t[0]].eta(), cands[t[0]].phi(), massD[0]);
  ROOT::Math::PtEtaPhiMVector p2(cands[t[1]].pt(), cands[t[1]].eta(), cands[t[1]].phi(), massD[1]);
  ROOT::Math::PtEtaPhiMVector p3(cands[t[2]].pt(), cands[t[2]].eta(), cands[t[2]].phi(), massD[2]);
  ROOT::Math::PtEtaPhiMVector p4(cands[t[3]].pt(), cands[t[3]].eta(), cands[t[3]].phi(), massD[3]);
  float mass = (p1 + p2 + p3 + p4).M();
  return mass;
}

void ScPhase2PuppiH2PhiDemo::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src");
  desc.add<bool>("runStruct", true);
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScPhase2PuppiH2PhiDemo);
