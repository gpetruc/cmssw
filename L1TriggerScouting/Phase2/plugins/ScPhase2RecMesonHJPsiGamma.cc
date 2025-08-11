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

class ScPhase2RecMesonHJPsiGamma : public edm::stream::EDProducer<> {
public:
  explicit ScPhase2RecMesonHJPsiGamma(const edm::ParameterSet &);
  ~ScPhase2RecMesonHJPsiGamma() override;
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  void beginStream(edm::StreamID) override;
  void produce(edm::Event &, const edm::EventSetup &) override;
  void endStream() override;
  template <typename T>
  void runObj(const OrbitCollection<T> &src,
              const OrbitCollection<T> &srcMeson,
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
    float minptGamma = 30;
    float maxiso = 0.25;
    float mindr2tkem = 0.05 * 0.05;
    float maxdr2tkem = 0.25 * 0.25;
    float maxisotkem = 0.25;
  } cuts;

  template <typename T>
  bool isolationQ(unsigned int pidex1, unsigned int pidex2, const T *cands, unsigned int size) const;

  std::tuple<bool, float> deltar(float eta1, float eta2, float phi1, float phi2) const;

  template <typename T, typename U>
  float tripletmass(const std::array<unsigned int, 2> &t,
                    const T *candsMeson,
                    const U *cands);

  template <typename T>
  bool isolationTkEm(float pt, float eta, float phi, const T *cands, unsigned int size) const;

  template <typename T, typename U>
  float tripletpt(const std::array<unsigned int, 2> &t,
                    const T *candsMeson,
                    const U *cands);

  unsigned long countStruct_;
  unsigned long passStruct_;
};

ScPhase2RecMesonHJPsiGamma::ScPhase2RecMesonHJPsiGamma(const edm::ParameterSet &iConfig)
    : doStruct_(iConfig.getParameter<bool>("runStruct")) {
  if (doStruct_) {
    structToken_ = consumes<OrbitCollection<l1Scouting::RecMeson>>(iConfig.getParameter<edm::InputTag>("src"));
    structToken_ = consumes<OrbitCollection<l1Scouting::RecMeson>>(iConfig.getParameter<edm::InputTag>("srcMeson"));
    produces<std::vector<unsigned>>("selectedBx");
    produces<l1ScoutingRun3::OrbitFlatTable>("recMesonHjpsigamma");
  }
}

ScPhase2RecMesonHJPsiGamma::~ScPhase2RecMesonHJPsiGamma() {};

void ScPhase2RecMesonHJPsiGamma::beginStream(edm::StreamID) {
  countStruct_ = 0;
  passStruct_ = 0;
}

void ScPhase2RecMesonHJPsiGamma::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  if (doStruct_) {
    edm::Handle<OrbitCollection<l1Scouting::RecMeson>> src;
    edm::Handle<OrbitCollection<l1Scouting::RecMeson>> srcMeson;
    iEvent.getByToken(structToken_, src);
    iEvent.getByToken(structToken_, srcMeson);
    runObj(*src, *srcMeson, iEvent, countStruct_, passStruct_, "");
  }
}

void ScPhase2RecMesonHJPsiGamma::endStream() {
  if (doStruct_)
    edm::LogImportant("ScPhase2AnalysisSummary") << "Rec Meson HJPsiGamma Struct analysis: " << countStruct_ << " -> " << passStruct_;
}

template <typename T>
void ScPhase2RecMesonHJPsiGamma::runObj(const OrbitCollection<T> &src,
                                    const OrbitCollection<T> &srcMeson,
                                    edm::Event &iEvent,
                                    unsigned long &nTry,
                                    unsigned long &nPass,
                                    const std::string &label) {
  l1ScoutingRun3::BxOffsetsFillter bxOffsetsFiller;
  bxOffsetsFiller.start();
  auto ret = std::make_unique<std::vector<unsigned>>();
  std::vector<float> masses;
  std::vector<uint8_t> i0s, i1s, i2s, i3s;
  ROOT::RVec<unsigned int> ig;
  std::array<unsigned int, 2> bestMesonPair;
  float bestTripletScore;
  bool bestTripletFound;

  for (unsigned int bx = 1; bx <= OrbitCollection<T>::NBX; ++bx) {
    nTry++;

    auto range = src.bxIterator(bx);
    const T *cands = &range.front();
    auto size = range.size();

    auto rangeMesons = srcMeson.bxIterator(bx);
    const T *candsMeson = &rangeMesons.front();
    unsigned int nMesons = rangeMesons.size();

    if (size < 1 || nMesons < 1) continue;

    ig.clear();
    for (unsigned int i = 0; i < size; ++i) {  // make list of all photons
      
      // photon isolation
      bool isop = isolationTkEm(
        cands[i].pt(), cands[i].eta(), cands[i].phi(), cands, size);
      if (!isop)
        continue;

      if (cands[i].pt() >= cuts.minptGamma) {
        ig.push_back(i);
      }
    }
    unsigned int ngammas = ig.size();
    if (ngammas < 1)
      continue;

    bestTripletScore = 0.;
    bestTripletFound = false;
    std::array<unsigned int, 2> bestTriplet{{0, 0}};
    auto mass = 0.;
    auto pt = 0.;

    for (unsigned int i1 = 0; i1 < nMesons; ++i1) {

      if (candsMeson[i1].pt() < cuts.minptQ)
        continue;
                
      for (unsigned int i2 = 0; i2 < ngammas; ++i2) {

        std::array<unsigned int, 2> pair{{i1, ig[i2]}};
        mass = tripletmass(pair, candsMeson, cands);
        pt = tripletpt(pair, candsMeson, cands);

        if (!(mass >= cuts.minmassH and mass <= cuts.maxmassH))
          continue;

        if (pt > bestTripletScore){
          bestTripletFound = true;
          bestTriplet = pair;
          bestTripletScore = pt;
        }
      }
    }
    
    if (!bestTripletFound)
      continue;

    ret->emplace_back(bx);
    nPass++;
    masses.push_back(mass);
    i0s.push_back(candsMeson[bestTriplet[0]].id1());
    i1s.push_back(candsMeson[bestTriplet[0]].id2());
    i2s.push_back(bestTriplet[1]);
    bxOffsetsFiller.addBx(bx, 1);
  }  // loop on BXs

  iEvent.put(std::move(ret), "selectedBx" + label);
  // now we make the table
  auto bxOffsets = bxOffsetsFiller.done();
  auto tab = std::make_unique<l1ScoutingRun3::OrbitFlatTable>(bxOffsets, "recMesonHjpsigamma" + label, true);
  tab->addColumn<float>("mass", masses, "2 muons plus photon invariant mass");
  tab->addColumn<uint8_t>("i0", i0s, "leading muon");
  tab->addColumn<uint8_t>("i1", i1s, "subleading muon");
  tab->addColumn<uint8_t>("i2", i2s, "photon");
  iEvent.put(std::move(tab), "recMesonHjpsigamma" + label);
}

template <typename T>
bool ScPhase2RecMesonHJPsiGamma::isolationTkEm(float pt, float eta, float phi, const T *cands, unsigned int size) const {
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

template <typename T, typename U>
float ScPhase2RecMesonHJPsiGamma::tripletmass(const std::array<unsigned int, 2> &t,
                                              const T *candsMeson,
                                              const U *cands) {
  ROOT::Math::PtEtaPhiMVector p1(candsMeson[t[0]].pt(), candsMeson[t[0]].eta(), candsMeson[t[0]].phi(), candsMeson[t[0]].mass());
  ROOT::Math::PtEtaPhiMVector p2(cands[t[1]].pt(), cands[t[1]].eta(), cands[t[1]].phi(), 0);
  float mass = (p1 + p2).M();
  return mass;
}

template <typename T, typename U>
float ScPhase2RecMesonHJPsiGamma::tripletpt(const std::array<unsigned int, 2> &t,
                                              const T *candsMeson,
                                              const U *cands) {
  ROOT::Math::PtEtaPhiMVector p1(candsMeson[t[0]].pt(), candsMeson[t[0]].eta(), candsMeson[t[0]].phi(), candsMeson[t[0]].mass());
  ROOT::Math::PtEtaPhiMVector p2(cands[t[1]].pt(), cands[t[1]].eta(), cands[t[1]].phi(), 0);
  float pt = (p1 + p2).Pt();
  return pt;
}

void ScPhase2RecMesonHJPsiGamma::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src");
  desc.add<edm::InputTag>("srcMeson");
  desc.add<bool>("runStruct", true);
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScPhase2RecMesonHJPsiGamma);