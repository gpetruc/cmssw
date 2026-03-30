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

class ScPhase2PuppiHJPsiGammaDemo : public edm::stream::EDProducer<> {
public:
  explicit ScPhase2PuppiHJPsiGammaDemo(const edm::ParameterSet &);
  ~ScPhase2PuppiHJPsiGammaDemo() override;
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
    float minpt1 = 5;
    float minpt2 = 10;
    float minpt3 = 30;
    float minptQ = 30;
    float maxdeltar2 = 0.40 * 0.40;
    float minmass = 100;
    float maxmass = 150;
    float minmass2 = 2.50;
    float maxmass2 = 3.50;
    float mindr2 = 0.05 * 0.05;
    float maxdr2 = 0.25 * 0.25;
    float maxiso = 0.25;
    float mindr2tkem = 0.05 * 0.05;
    float maxdr2tkem = 0.25 * 0.25;
    float maxisotkem = 0.25;
  } cuts;

  template <typename T>
  std::tuple<bool, float> isolationQ(
      float eta, float phi, unsigned int pidex1, unsigned int pidex2, const T *cands, unsigned int size) const;

  template <typename T>
  std::tuple<bool, float> isolationTkEm(float pt, float eta, float phi, const T *cands, unsigned int size) const;

  std::tuple<bool, float> deltar(float eta1, float eta2, float phi1, float phi2) const;

  template <typename T>
  static ROOT::Math::PtEtaPhiMVector pairP4(const std::array<unsigned int, 2> &t,
                                            const T *cands,
                                            const std::array<float, 2> &massD);

  template <typename T, typename U>
  float tripletmass(const std::array<unsigned int, 3> &t,
                    const T *cands,
                    const U *cands2,
                    const std::array<float, 3> &masses);

  unsigned long countStruct_;
  unsigned long passStruct_;
};

ScPhase2PuppiHJPsiGammaDemo::ScPhase2PuppiHJPsiGammaDemo(const edm::ParameterSet &iConfig)
    : doStruct_(iConfig.getParameter<bool>("runStruct")) {
  if (doStruct_) {
    structPuppiToken_ = consumes<OrbitCollection<l1Scouting::Puppi>>(iConfig.getParameter<edm::InputTag>("srcPuppi"));
    structTkEmToken_ = consumes<OrbitCollection<l1Scouting::TkEm>>(iConfig.getParameter<edm::InputTag>("srcTkEm"));
    produces<std::vector<unsigned>>("selectedBx");
    produces<l1ScoutingRun3::OrbitFlatTable>("hjpsigamma");
  }
}

ScPhase2PuppiHJPsiGammaDemo::~ScPhase2PuppiHJPsiGammaDemo() {};

void ScPhase2PuppiHJPsiGammaDemo::beginStream(edm::StreamID) {
  countStruct_ = 0;
  passStruct_ = 0;
}

void ScPhase2PuppiHJPsiGammaDemo::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  if (doStruct_) {
    edm::Handle<OrbitCollection<l1Scouting::Puppi>> srcPuppi;
    iEvent.getByToken(structPuppiToken_, srcPuppi);

    edm::Handle<OrbitCollection<l1Scouting::TkEm>> srcTkEm;
    iEvent.getByToken(structTkEmToken_, srcTkEm);

    runObj(*srcPuppi, *srcTkEm, iEvent, countStruct_, passStruct_, "");
  }
}

void ScPhase2PuppiHJPsiGammaDemo::endStream() {
  if (doStruct_)
    edm::LogImportant("ScPhase2AnalysisSummary")
        << "HJPsiGamma Struct analysis: " << countStruct_ << " -> " << passStruct_;
}

template <typename T, typename U>
void ScPhase2PuppiHJPsiGammaDemo::runObj(const OrbitCollection<T> &srcPuppi,
                                         const OrbitCollection<U> &srcTkEm,
                                         edm::Event &iEvent,
                                         unsigned long &nTry,
                                         unsigned long &nPass,
                                         const std::string &label) {
  l1ScoutingRun3::BxOffsetsFillter bxOffsetsFiller;
  bxOffsetsFiller.start();
  auto ret = std::make_unique<std::vector<unsigned>>();
  std::vector<float> masses;
  std::vector<uint8_t> i0s, i1s, i2s;                  // i2s is the photon
  ROOT::RVec<unsigned int> ix;                         //
  ROOT::RVec<unsigned int> ig;                         // photons
  std::array<unsigned int, 3> bestTriplet{{0, 0, 0}};  // indices of the best triplet (2 hadrons + 1 photon)
  bool bestPairFound;
  float bestPairScore;
  float bestPhoScore;
  ROOT::Math::PtEtaPhiMVector p4Pair, p4Pho;
  for (unsigned int bx = 1; bx <= OrbitCollection<T>::NBX; ++bx) {
    nTry++;
    auto range = srcPuppi.bxIterator(bx);
    auto size = range.size();
    const T *cands = (size > 0) ? &range.front() : nullptr;

    auto range2 = srcTkEm.bxIterator(bx);
    auto size2 = range2.size();
    const U *cands2 = (size2 > 0) ? &range2.front() : nullptr;

    ix.clear();
    int highcut = 0;
    for (unsigned int i = 0; i < size; ++i) {  // make list of puppi with a certain threshold
      if (cands[i].pt() >= cuts.minpt1) {
        ix.push_back(i);
        if (cands[i].pt() >= cuts.minpt2)
          highcut++;
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
    std::stable_sort(
        ix.begin(), ix.end(), [&](unsigned int i1, unsigned int i2) { return cands[i1].pt() > cands[i2].pt(); });
    bestPairFound = false;
    bestPairScore = cuts.minptQ;
    for (unsigned int i1 = 0; i1 < ndaus; ++i1) {
      if (cands[ix[i1]].pt() < cuts.minpt2)
        break;  // high pt cut
      for (unsigned int i2 = i1 + 1; i2 < ndaus; ++i2) {
        if (!(cands[ix[i1]].charge() * cands[ix[i2]].charge() < 0))
          continue;

        auto p4 = pairP4({{ix[i1], ix[i2]}}, cands, {{0.1057, 0.1057}});
        auto mass2 = p4.M();
        if (!(mass2 >= cuts.minmass2 and mass2 <= cuts.maxmass2))
          continue;

        auto [drcond, dr2Q] =
            deltar(cands[ix[i1]].eta(), cands[ix[i2]].eta(), cands[ix[i1]].phi(), cands[ix[i2]].phi());
        if (!drcond)
          continue;  // angular sep of top 2 tracks

        auto [isolated, psum] = isolationQ(p4.eta(), p4.phi(), ix[i1], ix[i2], cands, size);
        if (!isolated)
          continue;  // Q isolation

        std::array<unsigned int, 2> pair{{ix[i1], ix[i2]}};  // pair of indices
        if (p4.pt() > bestPairScore) {
          std::copy_n(pair.begin(), 2, bestTriplet.begin());
          bestPairFound = true;
          bestPairScore = p4.pt();
          p4Pair = p4;
        }
      }
    }
    if (!bestPairFound)
      continue;

    // photon
    bestPhoScore = 0;
    for (unsigned int i3 = 0; i3 < ngammas; ++i3) {
      auto [isop, psum] = isolationTkEm(cands2[ig[i3]].pt(), cands2[ig[i3]].eta(), cands2[ig[i3]].phi(), cands, size);
      if (!isop)
        continue;

      if (cands2[ig[i3]].pt() > bestPhoScore) {
        bestPhoScore = cands2[ig[i3]].pt();
        bestTriplet[2] = ig[i3];
        p4Pho = ROOT::Math::PtEtaPhiMVector(cands2[ig[i3]].pt(), cands2[ig[i3]].eta(), cands2[ig[i3]].phi(), 0.);
      }
    }
    if (bestPhoScore <= 0)
      continue;

    // H mass
    auto mass = (p4Pair + p4Pho).M();
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
  auto tab = std::make_unique<l1ScoutingRun3::OrbitFlatTable>(bxOffsets, "HJPsiGamma" + label, true);
  tab->addColumn<float>("mass", masses, "2 muons plus photon invariant mass");
  tab->addColumn<uint8_t>("i0", i0s, "leading muon");
  tab->addColumn<uint8_t>("i1", i1s, "subleading muon");
  tab->addColumn<uint8_t>("i2", i2s, "photon");
  iEvent.put(std::move(tab), "hjpsigamma" + label);
}

//TEST functions
template <typename T>
std::tuple<bool, float> ScPhase2PuppiHJPsiGammaDemo::isolationQ(
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

template <typename T>
std::tuple<bool, float> ScPhase2PuppiHJPsiGammaDemo::isolationTkEm(
    float pt, float eta, float phi, const T *cands, unsigned int size) const {
  float psum = 0;
  for (unsigned int j = 0u; j < size; ++j) {  //loop over other particles
    float dr2 = reco::deltaR2(eta, phi, cands[j].eta(), cands[j].phi());
    if (dr2 >= cuts.mindr2tkem && dr2 <= cuts.maxdr2tkem)
      psum += cands[j].pt();
  }
  bool passed = (psum <= cuts.maxisotkem * pt);
  return std::tuple(passed, psum);
}

std::tuple<bool, float> ScPhase2PuppiHJPsiGammaDemo::deltar(float eta1, float eta2, float phi1, float phi2) const {
  float dr2 = reco::deltaR2(eta1, phi1, eta2, phi2);
  bool passed = (dr2 <= cuts.maxdeltar2);
  return std::tuple(passed, dr2);
}

template <typename T>
ROOT::Math::PtEtaPhiMVector ScPhase2PuppiHJPsiGammaDemo::pairP4(const std::array<unsigned int, 2> &t,
                                                                const T *cands,
                                                                const std::array<float, 2> &massD) {
  ROOT::Math::PtEtaPhiMVector p1(cands[t[0]].pt(), cands[t[0]].eta(), cands[t[0]].phi(), massD[0]);
  ROOT::Math::PtEtaPhiMVector p2(cands[t[1]].pt(), cands[t[1]].eta(), cands[t[1]].phi(), massD[1]);
  return p1 + p2;
}

void ScPhase2PuppiHJPsiGammaDemo::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("srcPuppi");
  desc.add<edm::InputTag>("srcTkEm");
  desc.add<bool>("runStruct", true);
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScPhase2PuppiHJPsiGammaDemo);
