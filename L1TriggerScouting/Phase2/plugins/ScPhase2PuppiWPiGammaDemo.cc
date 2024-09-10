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

class ScPhase2PuppiWPiGammaDemo : public edm::stream::EDProducer<> {
public:
  explicit ScPhase2PuppiWPiGammaDemo(const edm::ParameterSet &);
  ~ScPhase2PuppiWPiGammaDemo() override;
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
  edm::EDGetTokenT<OrbitCollection<l1Scouting::Puppi>> structToken_;
  edm::EDGetTokenT<OrbitCollection<l1Scouting::TkEm>> struct2Token_;

  struct Cuts {
    float minpt_pi = 25;
    float minpt_tkem = 20;
    float minmass = 60;
    float maxmass = 100;
    float maxiso_pi = 0.30;
    float maxiso_tkem = 0.30;
    float mindeltar2 = 0.50 * 0.50;
    float mindr2 = 0.00 * 0.00;
    float maxdr2 = 0.50 * 0.50;
    float mindr2tkem = 0.02 * 0.02;
    float maxdr2tkem = 0.50 * 0.50;
  } cuts;

  template <typename T>
  bool isolationPi(
      unsigned int pidex1, const T *cands, unsigned int size) const;

  template <typename T>
  bool isolationPi(unsigned int pidex1,
                   const T *cands,
                   unsigned int size,
                   unsigned int &cache) const {
    if (cache == 0)
      cache = isolationPi(pidex1, cands, size) ? 1 : 2;
    return (cache == 1);
  }

  template <typename T>
  bool isolationTkEm(float pt, float eta, float phi, const T *cands, unsigned int size) const;

  template <typename T>
  bool isolationTkEm(float pt, float eta, float phi, const T *cands, unsigned int size, unsigned int &cache) const {
    if (cache == 0)
      cache = isolationTkEm(pt, eta, phi, cands, size) ? 1 : 2;
    return (cache == 1);
  }

  bool deltar(float eta1, float eta2, float phi1, float phi2) const;
  bool deltarmin(float eta1, float eta2, float phi1, float phi2) const;
  bool deltaphi(float phi1, float phi2) const;
  static float doubletmass(const std::array<unsigned int, 2> &t, const float *pts, const float *etas, const float *phis);

  unsigned long countStruct_;
  unsigned long passStruct_;
};

ScPhase2PuppiWPiGammaDemo::ScPhase2PuppiWPiGammaDemo(const edm::ParameterSet &iConfig)
    : doStruct_(iConfig.getParameter<bool>("runStruct")) {
  if (doStruct_) {
    structToken_ = consumes<OrbitCollection<l1Scouting::Puppi>>(iConfig.getParameter<edm::InputTag>("src"));
    struct2Token_ = consumes<OrbitCollection<l1Scouting::TkEm>>(iConfig.getParameter<edm::InputTag>("src2"));
    produces<std::vector<unsigned>>("selectedBx");
    produces<l1ScoutingRun3::OrbitFlatTable>("wdsgamma");
  }
}

ScPhase2PuppiWPiGammaDemo::~ScPhase2PuppiWPiGammaDemo(){};

void ScPhase2PuppiWPiGammaDemo::beginStream(edm::StreamID) {
  countStruct_ = 0;
  passStruct_ = 0;
}

void ScPhase2PuppiWPiGammaDemo::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  if (doStruct_) {
    edm::Handle<OrbitCollection<l1Scouting::Puppi>> src;
    iEvent.getByToken(structToken_, src);

    edm::Handle<OrbitCollection<l1Scouting::TkEm>> src2;
    iEvent.getByToken(struct2Token_, src2);

    runObj(*src, *src2, iEvent, countStruct_, passStruct_, "");
  }
}

void ScPhase2PuppiWPiGammaDemo::endStream() {
  if (doStruct_)
    std::cout << "Struct analysis: " << countStruct_ << " -> " << passStruct_ << std::endl;
}

template <typename T, typename U>
void ScPhase2PuppiWPiGammaDemo::runObj(const OrbitCollection<T> &src,
                                       const OrbitCollection<U> &src2,
                                       edm::Event &iEvent,
                                       unsigned long &nTry,
                                       unsigned long &nPass,
                                       const std::string &label) {
  l1ScoutingRun3::BxOffsetsFillter bxOffsetsFiller;
  bxOffsetsFiller.start();
  auto ret = std::make_unique<std::vector<unsigned>>();
  std::vector<float> masses;
  std::vector<uint8_t> i0s, i1s;            // i1 is the photon
  ROOT::RVec<unsigned int> ix;              // pions
  ROOT::RVec<unsigned int> ig;              // photons
  ROOT::RVec<unsigned int>
      iso;  //stores whether the Pi or photon passes isolation test so we don't calculate reliso twice
  std::array<unsigned int, 4> bestDoublet;
  float bestDoubletScore, bestDoubletMass;
  for (unsigned int bx = 1; bx <= OrbitCollection<T>::NBX; ++bx) {
    nTry++;
    auto range = src.bxIterator(bx);
    const T *cands = &range.front();
    auto size = range.size();

    auto range2 = src2.bxIterator(bx);
    const U *cands2 = &range2.front();
    auto size2 = range2.size();

    ix.clear();
    for (unsigned int i = 0; i < size; ++i) {  //make list of all hadrons
      if ((std::abs(cands[i].pdgId()) == 211 or std::abs(cands[i].pdgId()) == 11)) {
        if (cands[i].pt() >= cuts.minpt_pi) {
          ix.push_back(i);
	}
      }
    }
    unsigned int npions = ix.size();
    if (npions < 1)
      continue;

    ig.clear();
    for (unsigned int i = 0; i < size2; ++i) {  //make list of all photons
      if (cands2[i].pt() >= cuts.minpt_tkem) {
        ig.push_back(i);
      }
    }
    unsigned int ngammas = ig.size();
    if (ngammas < 1)
      continue;

    iso.resize(2);  //gamma and Pi isolations
    std::fill(iso.begin(), iso.end(), 0);
    bestDoubletScore = 0;

    for (unsigned int i1 = 0; i1 < npions; ++i1) {
      if (cands[ix[i1]].pt() < cuts.minpt_pi)
        continue;  //pion pt cut

      if (!isolationPi(ix[i1], cands, size, iso[0]))
        continue; //ISOLATION test for pion

      for (unsigned int i2 = 0; i2 < ngammas; ++i2) {
        if (cands2[ig[i2]].pt() < cuts.minpt_tkem)
          continue;  //photon pt cut

        std::array<unsigned int, 2> tr{{ix[i1], ig[i2]}};  //doublet of indices

        auto mass = (cands[ix[i1]].p4() + cands2[ig[i2]].p4()).mass();
        if (mass >= cuts.minmass and mass <= cuts.maxmass) {  //MASS test
          bool isop = isolationTkEm(
                      cands2[ig[i2]].pt(), cands2[ig[i2]].eta(), cands2[ig[i2]].phi(), cands, size, iso[1]);
          bool pass_deltar = deltar(cands[ix[i1]].eta(), cands[ig[i2]].eta(), cands[ix[i1]].phi(), cands[ig[i2]].phi());
          if (isop == true and pass_deltar == true) { //ISOLATION and DR tests
            float ptsum = cands[ix[i1]].pt() + cands2[ig[i2]].pt();
            if (ptsum > bestDoubletScore) {
              std::copy_n(tr.begin(), 2, bestDoublet.begin());
              bestDoubletScore = ptsum;
              bestDoubletMass = mass;
            }
          }  //iso and dr tests
        }    //mass test
      }      //photon loop
    }        //pion loop

    if (bestDoubletScore > 0) {
      ret->emplace_back(bx);
      nPass++;
      masses.push_back(bestDoubletMass);
      i0s.push_back(bestDoublet[0]);
      i1s.push_back(bestDoublet[1]);
      bxOffsetsFiller.addBx(bx, 1);
    }
  }  // loop on BXs

  iEvent.put(std::move(ret), "selectedBx" + label);
  // now we make the table
  auto bxOffsets = bxOffsetsFiller.done();
  auto tab = std::make_unique<l1ScoutingRun3::OrbitFlatTable>(bxOffsets, "WPiGamma" + label, true);
  tab->addColumn<float>("mass", masses, "pi-gamma invariant mass");
  tab->addColumn<uint8_t>("i0", i0s, "Pion");
  tab->addColumn<uint8_t>("i1", i1s, "Photon");
  iEvent.put(std::move(tab), "wpigamma" + label);
}

//TEST functions
template <typename T>
bool ScPhase2PuppiWPiGammaDemo::isolationPi(
    unsigned int pidex1, const T *cands, unsigned int size) const {
  bool passed = false;
  float psum = 0;
  float eta = cands[pidex1].eta();  //center cone around leading track
  float phi = cands[pidex1].phi();
  for (unsigned int j = 0u; j < size; ++j) {  //loop over other particles
    if (pidex1 == j)
      continue;
    float deta = eta - cands[j].eta(), dphi = ROOT::VecOps::DeltaPhi<float>(phi, cands[j].phi());
    float dr2 = deta * deta + dphi * dphi;
    if (dr2 >= cuts.mindr2 && dr2 <= cuts.maxdr2)
      psum += cands[j].pt();
  }
  if (psum <= cuts.maxiso_pi * cands[pidex1].pt())
    passed = true;
  return passed;
}

template <typename T>
bool ScPhase2PuppiWPiGammaDemo::isolationTkEm(float pt, float eta, float phi, const T *cands, unsigned int size) const {
  bool passed = false;
  float psum = 0;
  for (unsigned int j = 0u; j < size; ++j) {  //loop over other particles
    float deta = eta - cands[j].eta(), dphi = ROOT::VecOps::DeltaPhi<float>(phi, cands[j].phi());
    float dr2 = deta * deta + dphi * dphi;
    if (dr2 >= cuts.mindr2tkem && dr2 <= cuts.maxdr2tkem)
      psum += cands[j].pt();
  }
  if (psum <= cuts.maxiso_tkem * pt)
    passed = true;
  return passed;
}

bool ScPhase2PuppiWPiGammaDemo::deltar(float eta1, float eta2, float phi1, float phi2) const {
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

float ScPhase2PuppiWPiGammaDemo::doubletmass(const std::array<unsigned int, 2> &t,
                                                const float *pts,
                                                const float *etas,
                                                const float *phis) {
  ROOT::Math::PtEtaPhiMVector p1(pts[t[0]], etas[t[0]], phis[t[0]], 0.2);
  ROOT::Math::PtEtaPhiMVector p2(pts[t[1]], etas[t[1]], phis[t[1]], 0.00);
  float mass = (p1 + p2).M();
  return mass;
}

void ScPhase2PuppiWPiGammaDemo::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScPhase2PuppiWPiGammaDemo);
