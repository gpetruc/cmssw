#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/L1Scouting/interface/OrbitCollection.h"
#include "DataFormats/NanoAOD/interface/OrbitFlatTable.h"
#include "DataFormats/L1TParticleFlow/interface/PFCandidate.h"
#include "DataFormats/L1TParticleFlow/interface/L1ScoutingPuppi.h"
#include "L1TriggerScouting/Utilities/interface/BxOffsetsFiller.h"

#include <ROOT/RVec.hxx>
#include <Math/Vector4D.h>
#include <Math/GenVector/LorentzVector.h>
#include <Math/GenVector/PtEtaPhiM4D.h>
#include <algorithm>
#include <array>
#include <iostream>

class ScPhase2PuppiW3PiDemo : public edm::stream::EDProducer<> {
public:
  explicit ScPhase2PuppiW3PiDemo(const edm::ParameterSet &);
  ~ScPhase2PuppiW3PiDemo() override;
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
  edm::EDGetTokenT<OrbitCollection<l1t::PFCandidate>> candidateToken_;
  edm::EDGetTokenT<OrbitCollection<l1Scouting::Puppi>> structToken_;

  struct Cuts {
    float minpt1 = 7;   // 9
    float minpt2 = 12;  // 15
    float minpt3 = 15;  // 20
    float mindeltar2 = 0.5 * 0.5;
    float minmass = 40;   // 60
    float maxmass = 150;  // 100
    float mindr2 = 0.01 * 0.01;
    float maxdr2 = 0.25 * 0.25;
    float maxiso = 2.0;  //0.4
  } cuts;

  template <typename T>
  bool isolation(unsigned int pidex, const T *cands, unsigned int size) const;

  template <typename T>
  bool isolation(unsigned int pidex, const T *cands, unsigned int size, unsigned int &cache) const {
    if (cache == 0)
      cache = isolation(pidex, cands, size) ? 1 : 2;
    return (cache == 1);
  }

  bool deltar(float eta1, float eta2, float phi1, float phi2) const;

  unsigned long countCandidate_, countStruct_;
  unsigned long passCandidate_, passStruct_;
};

ScPhase2PuppiW3PiDemo::ScPhase2PuppiW3PiDemo(const edm::ParameterSet &iConfig)
    : doCandidate_(iConfig.getParameter<bool>("runCandidate")), doStruct_(iConfig.getParameter<bool>("runStruct")) {
  if (doCandidate_) {
    candidateToken_ = consumes<OrbitCollection<l1t::PFCandidate>>(iConfig.getParameter<edm::InputTag>("src"));
    produces<std::vector<unsigned>>("selectedBxCandidate");
    produces<l1ScoutingRun3::OrbitFlatTable>("w3piCandidate");
  }
  if (doStruct_) {
    structToken_ = consumes<OrbitCollection<l1Scouting::Puppi>>(iConfig.getParameter<edm::InputTag>("src"));
    produces<std::vector<unsigned>>("selectedBx");
    produces<l1ScoutingRun3::OrbitFlatTable>("w3pi");
  }
}

ScPhase2PuppiW3PiDemo::~ScPhase2PuppiW3PiDemo() {};

void ScPhase2PuppiW3PiDemo::beginStream(edm::StreamID) {
  countCandidate_ = 0;
  countStruct_ = 0;
  passCandidate_ = 0;
  passStruct_ = 0;
}

void ScPhase2PuppiW3PiDemo::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  if (doCandidate_) {
    edm::Handle<OrbitCollection<l1t::PFCandidate>> src;
    iEvent.getByToken(candidateToken_, src);
    runObj(*src, iEvent, countCandidate_, passCandidate_, "Candidate");
  }
  if (doStruct_) {
    edm::Handle<OrbitCollection<l1Scouting::Puppi>> src;
    iEvent.getByToken(structToken_, src);
    runObj(*src, iEvent, countStruct_, passStruct_, "");
  }
}

void ScPhase2PuppiW3PiDemo::endStream() {
  if (doCandidate_)
    edm::LogImportant("ScPhase2AnalysisSummary")
        << "W3Pi Candidate analysis: " << countCandidate_ << " -> " << passCandidate_;
  if (doStruct_)
    edm::LogImportant("ScPhase2AnalysisSummary") << "W3Pi Struct analysis: " << countStruct_ << " -> " << passStruct_;
}

template <typename T>
void ScPhase2PuppiW3PiDemo::runObj(const OrbitCollection<T> &src,
                                   edm::Event &iEvent,
                                   unsigned long &nTry,
                                   unsigned long &nPass,
                                   const std::string &label) {
  l1ScoutingRun3::BxOffsetsFillter bxOffsetsFiller;
  bxOffsetsFiller.start();
  auto ret = std::make_unique<std::vector<unsigned>>();
  std::vector<float> masses;
  std::vector<uint8_t> i0s, i1s, i2s;
  ROOT::RVec<unsigned int> ix;   // pions
  ROOT::RVec<unsigned int> iso;  //stores whether a particle passes isolation test so we don't calculate reliso twice
  std::array<unsigned int, 3> bestTriplet;  // best triplet
  float bestTripletScore, bestTripletMass;
  for (unsigned int bx = 1; bx <= OrbitCollection<T>::NBX; ++bx) {
    nTry++;
    auto range = src.bxIterator(bx);
    const T *cands = &range.front();
    auto size = range.size();
    ix.clear();
    int intermediatecut = 0;
    int highcut = 0;
    for (unsigned int i = 0; i < size; ++i) {  //make list of all hadrons
      if ((std::abs(cands[i].pdgId()) == 211 or std::abs(cands[i].pdgId()) == 11)) {
        if (cands[i].pt() >= cuts.minpt1) {
          ix.push_back(i);
          if (cands[i].pt() >= cuts.minpt2)
            intermediatecut++;
          if (cands[i].pt() >= cuts.minpt3)
            highcut++;
        }
      }
    }
    unsigned int npions = ix.size();
    if (highcut < 1 || intermediatecut < 2 || npions < 3)
      continue;
    iso.resize(npions);
    std::fill(iso.begin(), iso.end(), 0);
    bestTripletScore = 0;

    for (unsigned int i1 = 0; i1 < npions; ++i1) {
      if (cands[ix[i1]].pt() < cuts.minpt3)
        continue;  //high pt cut
      if (isolation(ix[i1], cands, size, iso[i1]) == 0)
        continue;  //check iso of high pt pion
      for (unsigned int i2 = 0; i2 < npions; ++i2) {
        if (i2 == i1 || cands[ix[i2]].pt() < cuts.minpt2)
          continue;
        if (cands[ix[i2]].pt() > cands[ix[i1]].pt() || (cands[ix[i2]].pt() == cands[ix[i1]].pt() and i2 < i1))
          continue;  //intermediate pt cut
        if (!deltar(cands[ix[i1]].eta(), cands[ix[i2]].eta(), cands[ix[i1]].phi(), cands[ix[i2]].phi()))
          continue;  //angular sep of top 2 pions
        for (unsigned int i3 = 0; i3 < npions; ++i3) {
          if (i3 == i1 or i3 == i2)
            continue;
          if (cands[ix[i2]].pt() < cuts.minpt1)
            continue;  //low pt cut
          if (cands[ix[i3]].pt() > cands[ix[i1]].pt() || (cands[ix[i3]].pt() == cands[ix[i1]].pt() and i3 < i1))
            continue;
          if (cands[ix[i3]].pt() > cands[ix[i2]].pt() || (cands[ix[i3]].pt() == cands[ix[i2]].pt() and i3 < i2))
            continue;
          std::array<unsigned int, 3> tr{{ix[i1], ix[i2], ix[i3]}};  //triplet of indeces

          if (std::abs(cands[ix[i1]].charge() + cands[ix[i2]].charge() + cands[ix[i3]].charge()) == 1) {
            //make Lorentz vectors for each triplet
            auto mass = (cands[ix[i1]].p4() + cands[ix[i2]].p4() + cands[ix[i3]].p4()).mass();
            if (mass >= cuts.minmass and mass <= cuts.maxmass) {  //MASS test
              if (deltar(cands[ix[i1]].eta(), cands[ix[i3]].eta(), cands[ix[i1]].phi(), cands[ix[i3]].phi()) and
                  deltar(cands[ix[i2]].eta(), cands[ix[i3]].eta(), cands[ix[i2]].phi(), cands[ix[i3]].phi())) {
                //ISOLATION test for lower 4 pions
                bool isop = isolation(ix[i2], cands, size, iso[i2]) && isolation(ix[i3], cands, size, iso[i3]);
                if (isop == true) {
                  float ptsum = cands[ix[i1]].pt() + cands[ix[i2]].pt() + cands[ix[i3]].pt();
                  if (ptsum > bestTripletScore) {
                    std::copy_n(tr.begin(), 3, bestTriplet.begin());
                    bestTripletScore = ptsum;
                    bestTripletMass = mass;
                  }
                }  // iso
              }  // delta R
            }  // mass
          }  //charge
        }  //low pt cut
      }  //intermediate pt cut
    }  //high pt cut

    if (bestTripletScore > 0) {
      ret->emplace_back(bx);
      nPass++;
      masses.push_back(bestTripletMass);
      i0s.push_back(bestTriplet[0]);
      i1s.push_back(bestTriplet[1]);
      i2s.push_back(bestTriplet[2]);
      bxOffsetsFiller.addBx(bx, 1);
    }
  }  // loop on BXs

  iEvent.put(std::move(ret), "selectedBx" + label);

  // now we make the table
  auto bxOffsets = bxOffsetsFiller.done();
  auto tab = std::make_unique<l1ScoutingRun3::OrbitFlatTable>(bxOffsets, "W3Pi" + label, true);
  tab->addColumn<float>("mass", masses, "3-pion invariant mass");
  tab->addColumn<uint8_t>("i0", i0s, "leading pion");
  tab->addColumn<uint8_t>("i1", i1s, "subleading pion");
  tab->addColumn<uint8_t>("i2", i2s, "trailing pion");
  iEvent.put(std::move(tab), "w3pi" + label);
}

//TEST functions
template <typename T>
bool ScPhase2PuppiW3PiDemo::isolation(unsigned int pidex, const T *cands, unsigned int size) const {
  bool passed = false;
  float psum = 0;
  float eta = cands[pidex].eta();
  float phi = cands[pidex].phi();
  for (unsigned int j = 0u; j < size; ++j) {  //loop over other particles
    if (pidex == j)
      continue;
    float deta = eta - cands[j].eta(), dphi = ROOT::VecOps::DeltaPhi<float>(phi, cands[j].phi());
    float dr2 = deta * deta + dphi * dphi;
    if (dr2 >= cuts.mindr2 && dr2 <= cuts.maxdr2)
      psum += cands[j].pt();
  }
  if (psum <= cuts.maxiso * cands[pidex].pt())
    passed = true;
  return passed;
}
bool ScPhase2PuppiW3PiDemo::deltar(float eta1, float eta2, float phi1, float phi2) const {
  bool passed = true;
  float deta = eta1 - eta2;
  float dphi = ROOT::VecOps::DeltaPhi<float>(phi1, phi2);
  float dr2 = deta * deta + dphi * dphi;
  if (dr2 < cuts.mindeltar2) {
    passed = false;
    return passed;
  }
  return passed;
}

void ScPhase2PuppiW3PiDemo::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src");
  desc.add<bool>("runStruct", true);
  desc.add<bool>("runCandidate", false);
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScPhase2PuppiW3PiDemo);
