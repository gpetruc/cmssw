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

class ScPhase2PuppiHRhoGammaDemo : public edm::stream::EDProducer<> {
public:
  explicit ScPhase2PuppiHRhoGammaDemo(const edm::ParameterSet &);
  ~ScPhase2PuppiHRhoGammaDemo() override;
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  void beginStream(edm::StreamID) override;
  void produce(edm::Event &, const edm::EventSetup &) override;
  void endStream() override;
  template <typename T, typename U>
  void runObj(const OrbitCollection<T> &src,
              const OrbitCollection<U> &src2,
              edm::Event &out,
              unsigned long &nTry,
              unsigned long &nPass,
              const std::string &bxLabel);

  bool doStruct_;
  edm::EDGetTokenT<OrbitCollection<l1Scouting::Puppi>> structPuppiToken_;
  edm::EDGetTokenT<OrbitCollection<l1Scouting::TkEm>> structTkEmToken_;

  struct Cuts {
    float minpt1 = 10;
    float minpt2 = 10;
    float minpt3 = 30;
    float minptQ = 30;
    float maxdeltar2 = 0.40 * 0.40;
    float minmass = 100;
    float maxmass = 150;
    float minmass2 = 0.40;
    float maxmass2 = 1.30;
    float mindr2 = 0.05 * 0.05;
    float maxdr2 = 0.25 * 0.25;
    float maxiso = 0.25;
    float mindr2tkem = 0.05 * 0.05;
    float maxdr2tkem = 0.25 * 0.25;
    float maxisotkem = 0.25;
  } cuts;

  template <typename T>
  bool isolationQ(unsigned int pidex1, unsigned int pidex2, const T *cands, unsigned int size) const;

  template <typename T>
  bool isolationTkEm(float pt, float eta, float phi, const T *cands, unsigned int size) const;

  std::tuple<bool, float> deltar(float eta1, float eta2, float phi1, float phi2) const;

  template <typename T>
  static float pairmass(const std::array<unsigned int, 2> &t, const T *cands, const std::array<float, 2> &massD);

  template <typename T, typename U>
  float tripletmass(const std::array<unsigned int, 3> &t,
                    const T *cands,
                    const U *cands2,
                    const std::array<float, 3> &masses);

  unsigned long countStruct_;
  unsigned long passStruct_;
};

ScPhase2PuppiHRhoGammaDemo::ScPhase2PuppiHRhoGammaDemo(const edm::ParameterSet &iConfig)
    : doStruct_(iConfig.getParameter<bool>("runStruct")) {
  if (doStruct_) {
    structPuppiToken_ = consumes<OrbitCollection<l1Scouting::Puppi>>(iConfig.getParameter<edm::InputTag>("srcPuppi"));
    structTkEmToken_ = consumes<OrbitCollection<l1Scouting::TkEm>>(iConfig.getParameter<edm::InputTag>("srcTkEm"));
    produces<std::vector<unsigned>>("selectedBx");
    produces<l1ScoutingRun3::OrbitFlatTable>("hrhogamma");
  }
}

ScPhase2PuppiHRhoGammaDemo::~ScPhase2PuppiHRhoGammaDemo() {};

void ScPhase2PuppiHRhoGammaDemo::beginStream(edm::StreamID) {
  countStruct_ = 0;
  passStruct_ = 0;
}

void ScPhase2PuppiHRhoGammaDemo::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  if (doStruct_) {
    edm::Handle<OrbitCollection<l1Scouting::Puppi>> srcPuppi;
    iEvent.getByToken(structPuppiToken_, srcPuppi);

    edm::Handle<OrbitCollection<l1Scouting::TkEm>> srcTkEm;
    iEvent.getByToken(structTkEmToken_, srcTkEm);

    runObj(*srcPuppi, *srcTkEm, iEvent, countStruct_, passStruct_, "");
  }
}

void ScPhase2PuppiHRhoGammaDemo::endStream() {
  if (doStruct_)
    edm::LogImportant("ScPhase2AnalysisSummary")
        << "HRhoGamma Struct analysis: " << countStruct_ << " -> " << passStruct_;
}

template <typename T, typename U>
void ScPhase2PuppiHRhoGammaDemo::runObj(const OrbitCollection<T> &srcPuppi,
                                        const OrbitCollection<U> &srcTkEm,
                                        edm::Event &iEvent,
                                        unsigned long &nTry,
                                        unsigned long &nPass,
                                        const std::string &label) {
  l1ScoutingRun3::BxOffsetsFillter bxOffsetsFiller;
  bxOffsetsFiller.start();
  auto ret = std::make_unique<std::vector<unsigned>>();
  std::vector<float> masses;
  std::vector<uint8_t> i0s, i1s, i2s;  // i2s is the photon
  ROOT::RVec<unsigned int> ix;         //
  ROOT::RVec<unsigned int> ig;         // photons
  std::array<unsigned int, 2> bestPair;
  std::array<unsigned int, 3> bestTriplet;
  bool bestPairFound;
  float bestPairScore;
  float bestPhoScore;
  for (unsigned int bx = 1; bx <= OrbitCollection<T>::NBX; ++bx) {
    nTry++;
    auto range = srcPuppi.bxIterator(bx);
    const T *cands = &range.front();
    auto size = range.size();

    auto range2 = srcTkEm.bxIterator(bx);
    const U *cands2 = &range2.front();
    auto size2 = range2.size();

    ix.clear();
    int highcut = 0;
    for (unsigned int i = 0; i < size; ++i) {  //make list of all hadrons
      if ((std::abs(cands[i].pdgId()) == 211 or std::abs(cands[i].pdgId()) == 11)) {
        if (cands[i].pt() >= cuts.minpt1) {
          ix.push_back(i);
          if (cands[i].pt() >= cuts.minpt2)
            highcut++;
        }
      }
    }
    unsigned int ndaus = ix.size();
    if (highcut < 1 || ndaus < 2)
      continue;

    ig.clear();
    for (unsigned int i = 0; i < size2; ++i) {  // make list of all photons
      if (cands2[i].pt() >= cuts.minpt3) {
        ig.push_back(i);
      }
    }
    unsigned int ngammas = ig.size();
    if (ngammas < 1)
      continue;

    // Q candidate from closest OS pair with mass compatible with mQ
    bestPairFound = false;
    bestPairScore = 999;
    for (unsigned int i1 = 0; i1 < ndaus; ++i1) {
      if (cands[ix[i1]].pt() < cuts.minpt2)
        continue;  // high pt cut
      for (unsigned int i2 = 0; i2 < ndaus; ++i2) {
        if (i2 == i1 || cands[ix[i2]].pt() < cuts.minpt1)
          continue;

        if (!(cands[ix[i1]].charge() * cands[ix[i2]].charge() < 0))
          continue;

        auto mass2 = pairmass({{ix[i1], ix[i2]}}, cands, {{0.1396, 0.1396}});
        if (!(mass2 >= cuts.minmass2 and mass2 <= cuts.maxmass2))
          continue;

        auto [drcond, drQ] = deltar(cands[ix[i1]].eta(), cands[ix[i2]].eta(), cands[ix[i1]].phi(), cands[ix[i2]].phi());
        if (!drcond)
          continue;  //angular sep of top 2 tracks

        std::array<unsigned int, 2> pair{{ix[i1], ix[i2]}};  // pair of indices
        if (drQ < bestPairScore) {
          std::copy_n(pair.begin(), 2, bestPair.begin());
          bestPairScore = drQ;
          if (bestPairScore * bestPairScore < cuts.maxdeltar2)
            bestPairFound = true;
        }
      }
    }
    if (!bestPairFound)
      continue;

    // Q pt
    auto ptQ = (cands[bestPair[0]].p4() + cands[bestPair[1]].p4()).pt();
    if (ptQ < cuts.minptQ)
      continue;

    // Q isolation
    if (!isolationQ(bestPair[0], bestPair[1], cands, size))
      continue;

    // photon
    bestPhoScore = 0;
    for (unsigned int i3 = 0; i3 < ngammas; ++i3) {
      if (cands2[ig[i3]].pt() < cuts.minpt3)
        continue;  // photon pt cut

      std::array<unsigned int, 3> tr{{bestPair[0], bestPair[1], ig[i3]}};  // triplet of indices

      if (cands2[ig[i3]].pt() > bestPhoScore) {
        std::copy_n(tr.begin(), 3, bestTriplet.begin());
        bestPhoScore = cands2[ig[i3]].pt();
      }
    }
    if (bestPhoScore < 0)
      continue;

    // photon isolation
    bool isop = isolationTkEm(
        cands2[bestTriplet[2]].pt(), cands2[bestTriplet[2]].eta(), cands2[bestTriplet[2]].phi(), cands, size);
    if (!isop)
      continue;

    // H mass
    auto mass = tripletmass(bestTriplet, cands, cands2, {{0.1396, 0.1396, 0.}});
    if (!(mass >= cuts.minmass and mass <= cuts.maxmass))
      continue;

    ret->emplace_back(bx);
    nPass++;
    masses.push_back(mass);
    i0s.push_back(bestTriplet[0]);
    i1s.push_back(bestTriplet[1]);
    i2s.push_back(bestTriplet[2]);
    bxOffsetsFiller.addBx(bx, 1);
  }  // loop on BXs

  iEvent.put(std::move(ret), "selectedBx" + label);
  // now we make the table
  auto bxOffsets = bxOffsetsFiller.done();
  auto tab = std::make_unique<l1ScoutingRun3::OrbitFlatTable>(bxOffsets, "HRhoGamma" + label, true);
  tab->addColumn<float>("mass", masses, "2 pions plus photon invariant mass");
  tab->addColumn<uint8_t>("i0", i0s, "leading pion");
  tab->addColumn<uint8_t>("i1", i1s, "subleading pion");
  tab->addColumn<uint8_t>("i2", i2s, "photon");
  iEvent.put(std::move(tab), "hrhogamma" + label);
}

//TEST functions
template <typename T>
bool ScPhase2PuppiHRhoGammaDemo::isolationQ(unsigned int pidex1,
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

template <typename T>
bool ScPhase2PuppiHRhoGammaDemo::isolationTkEm(float pt, float eta, float phi, const T *cands, unsigned int size) const {
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

std::tuple<bool, float> ScPhase2PuppiHRhoGammaDemo::deltar(float eta1, float eta2, float phi1, float phi2) const {
  bool passed = true;
  float deta = eta1 - eta2;
  float dphi = ROOT::VecOps::DeltaPhi<float>(phi1, phi2);
  float dr2 = deta * deta + dphi * dphi;
  if (dr2 > cuts.maxdeltar2) {
    passed = false;
    return std::tuple(passed, dr2);
  }
  return std::tuple(passed, dr2);
}

template <typename T>
float ScPhase2PuppiHRhoGammaDemo::pairmass(const std::array<unsigned int, 2> &t,
                                           const T *cands,
                                           const std::array<float, 2> &massD) {
  ROOT::Math::PtEtaPhiMVector p1(cands[t[0]].pt(), cands[t[0]].eta(), cands[t[0]].phi(), massD[0]);
  ROOT::Math::PtEtaPhiMVector p2(cands[t[1]].pt(), cands[t[1]].eta(), cands[t[1]].phi(), massD[1]);
  float mass = (p1 + p2).M();
  return mass;
}

template <typename T, typename U>
float ScPhase2PuppiHRhoGammaDemo::tripletmass(const std::array<unsigned int, 3> &t,
                                              const T *cands,
                                              const U *cands2,
                                              const std::array<float, 3> &masses) {
  ROOT::Math::PtEtaPhiMVector p1(cands[t[0]].pt(), cands[t[0]].eta(), cands[t[0]].phi(), masses[0]);
  ROOT::Math::PtEtaPhiMVector p2(cands[t[1]].pt(), cands[t[1]].eta(), cands[t[1]].phi(), masses[1]);
  ROOT::Math::PtEtaPhiMVector p3(cands2[t[2]].pt(), cands2[t[2]].eta(), cands2[t[2]].phi(), masses[2]);
  float mass = (p1 + p2 + p3).M();
  return mass;
}

void ScPhase2PuppiHRhoGammaDemo::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("srcPuppi");
  desc.add<edm::InputTag>("srcTkEm");
  desc.add<bool>("runStruct", true);
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScPhase2PuppiHRhoGammaDemo);
