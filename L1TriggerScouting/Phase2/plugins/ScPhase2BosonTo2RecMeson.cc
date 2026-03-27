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

class ScPhase2BosonTo2RecMeson : public edm::stream::EDProducer<> {
public:
  explicit ScPhase2BosonTo2RecMeson(const edm::ParameterSet &);
  ~ScPhase2BosonTo2RecMeson() override;
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  void beginStream(edm::StreamID) override;
  void produce(edm::Event &, const edm::EventSetup &) override;
  void endStream() override;
  template <typename T>
  void runObj(const OrbitCollection<T> &srcMeson1,
              const OrbitCollection<T> &srcMeson2,
              edm::Event &out,
              unsigned long &nTry,
              unsigned long &nPass,
              const std::string &bxLabel);
  template <typename T>
  void runObj(const OrbitCollection<T> &srcMeson,
              edm::Event &out,
              unsigned long &nTry,
              unsigned long &nPass,
              const std::string &bxLabel);

  edm::EDGetTokenT<OrbitCollection<l1Scouting::RecMeson<2>>> structToken1_;
  edm::EDGetTokenT<OrbitCollection<l1Scouting::RecMeson<2>>> structToken2_;

  double minMassBoson_;
  double maxMassBoson_;
  double minPtQ_;
  double maxIso_;
  std::string analysisName_;

  static float quadrimass(const l1Scouting::RecMeson<2> *candsa, int a, const l1Scouting::RecMeson<2> *candsb, int b);

  unsigned long countStruct_;
  unsigned long passStruct_;

  // If the decay is to Q+Q, the meson candidate collections 1 & 2 will be the same
  // If the deacy is to Q_1+Q_2, then no
  bool sameMesonCollection_;
};

ScPhase2BosonTo2RecMeson::ScPhase2BosonTo2RecMeson(const edm::ParameterSet &iConfig)
    : minMassBoson_(iConfig.getParameter<double>("minMassBoson")),
      maxMassBoson_(iConfig.getParameter<double>("maxMassBoson")),
      minPtQ_(iConfig.getParameter<double>("minPtQ")),
      maxIso_(iConfig.getParameter<double>("maxIso")),
      analysisName_(iConfig.getParameter<std::string>("analysisName")) {
  structToken1_ = consumes<OrbitCollection<l1Scouting::RecMeson<2>>>(iConfig.getParameter<edm::InputTag>("srcMeson1"));
  structToken2_ = consumes<OrbitCollection<l1Scouting::RecMeson<2>>>(iConfig.getParameter<edm::InputTag>("srcMeson2"));
  sameMesonCollection_ =
      (iConfig.getParameter<edm::InputTag>("srcMeson1") == iConfig.getParameter<edm::InputTag>("srcMeson2"));
  produces<std::vector<unsigned>>("selectedBx");
  produces<l1ScoutingRun3::OrbitFlatTable>(analysisName_);
}

ScPhase2BosonTo2RecMeson::~ScPhase2BosonTo2RecMeson() {};

void ScPhase2BosonTo2RecMeson::beginStream(edm::StreamID) {
  countStruct_ = 0;
  passStruct_ = 0;
}

void ScPhase2BosonTo2RecMeson::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  edm::Handle<OrbitCollection<l1Scouting::RecMeson<2>>> srcMeson1;
  edm::Handle<OrbitCollection<l1Scouting::RecMeson<2>>> srcMeson2;
  iEvent.getByToken(structToken1_, srcMeson1);
  iEvent.getByToken(structToken2_, srcMeson2);
  runObj(*srcMeson1, *srcMeson2, iEvent, countStruct_, passStruct_, "");
}

void ScPhase2BosonTo2RecMeson::endStream() {
  edm::LogImportant("ScPhase2AnalysisSummary")
      << "Rec Meson " << analysisName_ << " Struct analysis: " << countStruct_ << " -> " << passStruct_;
}

template <typename T>
void ScPhase2BosonTo2RecMeson::runObj(const OrbitCollection<T> &srcMeson1,
                                      const OrbitCollection<T> &srcMeson2,
                                      edm::Event &iEvent,
                                      unsigned long &nTry,
                                      unsigned long &nPass,
                                      const std::string &label) {
  l1ScoutingRun3::BxOffsetsFillter bxOffsetsFiller;
  bxOffsetsFiller.start();
  auto ret = std::make_unique<std::vector<unsigned>>();
  std::vector<float> masses;
  std::vector<uint8_t> i0s, i1s, i2s, i3s;
  std::array<unsigned int, 2> bestMesonPair{{0, 0}};
  float bestMesonPairScore;
  bool bestMesonPairFound;

  for (unsigned int bx = 1; bx <= OrbitCollection<T>::NBX; ++bx) {
    nTry++;

    auto rangeMeson1 = srcMeson1.bxIterator(bx);
    unsigned int nMeson1 = rangeMeson1.size();
    const T *candsMeson1 = (nMeson1 > 0) ? &rangeMeson1.front() : nullptr;

    auto rangeMeson2 = srcMeson2.bxIterator(bx);
    unsigned int nMeson2 = rangeMeson2.size();
    const T *candsMeson2 = (nMeson2 > 0) ? &rangeMeson2.front() : nullptr;

    bestMesonPairScore = 0.;
    bestMesonPairFound = false;

    if (nMeson1 < 1 || nMeson2 < 1)
      continue;

    for (unsigned int i1 = 0; i1 < nMeson1; ++i1) {
      // minimum pt and isolation of Q1
      if ((candsMeson1[i1].pt() < minPtQ_) || (candsMeson1[i1].isoDR0p25() >= maxIso_))
        continue;
      unsigned int start = sameMesonCollection_ ? i1 + 1 : 0;
      for (unsigned int i2 = start; i2 < nMeson2; ++i2) {
        // minimum pt and isolation of Q2
        if ((candsMeson2[i2].pt() < minPtQ_) || (candsMeson2[i2].isoDR0p25() >= maxIso_))
          continue;

        // Four different dauther particles
        if ((candsMeson1[i1].daughterIds(0) == candsMeson2[i2].daughterIds(0)) ||
            (candsMeson1[i1].daughterIds(0) == candsMeson2[i2].daughterIds(1)))
          continue;
        if ((candsMeson1[i1].daughterIds(1) == candsMeson2[i2].daughterIds(0)) ||
            (candsMeson1[i1].daughterIds(1) == candsMeson2[i2].daughterIds(1)))
          continue;

        // Choose best pair of mesons based on score (e.g. max pt)
        float ptsum = candsMeson1[i1].pt() + candsMeson2[i2].pt();
        if (ptsum > bestMesonPairScore) {
          bestMesonPairScore = ptsum;
          bestMesonPair = {{i1, i2}};
          bestMesonPairFound = true;
        }
      }
    }

    if (!bestMesonPairFound)
      continue;

    // Boson mass
    auto mass = quadrimass(candsMeson1, bestMesonPair[0], candsMeson2, bestMesonPair[1]);
    if (!(mass >= minMassBoson_ and mass <= maxMassBoson_))
      continue;

    ret->emplace_back(bx);
    nPass++;
    masses.push_back(mass);
    i0s.push_back(candsMeson1[bestMesonPair[0]].daughterIds(0));
    i1s.push_back(candsMeson1[bestMesonPair[0]].daughterIds(1));
    i2s.push_back(candsMeson2[bestMesonPair[1]].daughterIds(0));
    i3s.push_back(candsMeson2[bestMesonPair[1]].daughterIds(1));
    bxOffsetsFiller.addBx(bx, 1);
  }  // loop on BXs

  iEvent.put(std::move(ret), "selectedBx" + label);
  // now we make the table
  auto bxOffsets = bxOffsetsFiller.done();
  auto tab = std::make_unique<l1ScoutingRun3::OrbitFlatTable>(bxOffsets, analysisName_ + label, true);
  tab->addColumn<float>("mass", masses, "4 daughters invariant mass");
  tab->addColumn<uint8_t>("i0", i0s, "1st daughter (1st meson)");
  tab->addColumn<uint8_t>("i1", i1s, "2nd daughter (1st meson)");
  tab->addColumn<uint8_t>("i2", i2s, "1st daughter (2nd meson)");
  tab->addColumn<uint8_t>("i3", i3s, "2nd daughter (2nd meson)");
  iEvent.put(std::move(tab), analysisName_ + label);
}

float ScPhase2BosonTo2RecMeson::quadrimass(const l1Scouting::RecMeson<2> *candsa,
                                           int a,
                                           const l1Scouting::RecMeson<2> *candsb,
                                           int b) {
  ROOT::Math::PtEtaPhiMVector p1(candsa[a].pt(), candsa[a].eta(), candsa[a].phi(), candsa[a].mass());
  ROOT::Math::PtEtaPhiMVector p2(candsb[b].pt(), candsb[b].eta(), candsb[b].phi(), candsb[b].mass());
  float mass = (p1 + p2).M();
  return mass;
}

void ScPhase2BosonTo2RecMeson::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("srcMeson1");
  desc.add<edm::InputTag>("srcMeson2");
  desc.add<double>("minMassBoson");
  desc.add<double>("maxMassBoson");
  desc.add<double>("minPtQ");
  desc.add<double>("maxIso");
  desc.add<std::string>("analysisName");
  desc.add<bool>("runStruct", true);
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScPhase2BosonTo2RecMeson);
