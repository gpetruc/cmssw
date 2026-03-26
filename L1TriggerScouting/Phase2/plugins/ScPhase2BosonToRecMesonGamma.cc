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

class ScPhase2BosonToRecMesonGamma : public edm::stream::EDProducer<> {
public:
  explicit ScPhase2BosonToRecMesonGamma(const edm::ParameterSet &);
  ~ScPhase2BosonToRecMesonGamma() override;
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  void beginStream(edm::StreamID) override;
  void produce(edm::Event &, const edm::EventSetup &) override;
  void endStream() override;
  template <typename T, typename U>
  void runObj(const OrbitCollection<T> &srcGamma,
              const OrbitCollection<U> &srcMeson,
              edm::Event &out,
              unsigned long &nTry,
              unsigned long &nPass,
              const std::string &bxLabel);

  edm::EDGetTokenT<OrbitCollection<l1Scouting::RecMeson<2>>> structMesonToken_;
  edm::EDGetTokenT<OrbitCollection<l1Scouting::TkEm>> structGammaToken_;

  double minMassBoson_;
  double maxMassBoson_;
  double minPtQ_;
  double minPtGamma_;
  std::string analysisName_;

  std::tuple<bool, float> deltar(float eta1, float eta2, float phi1, float phi2) const;

  template <typename T, typename U>
  ROOT::Math::PtEtaPhiMVector tripletP4(const std::array<unsigned int, 2> &t, const T *candsGamma, const U *candsMeson);

  unsigned long countStruct_;
  unsigned long passStruct_;
};

ScPhase2BosonToRecMesonGamma::ScPhase2BosonToRecMesonGamma(const edm::ParameterSet &iConfig)
    : minMassBoson_(iConfig.getParameter<double>("minMassBoson")),
      maxMassBoson_(iConfig.getParameter<double>("maxMassBoson")),
      minPtQ_(iConfig.getParameter<double>("minPtQ")),
      minPtGamma_(iConfig.getParameter<double>("minPtGamma")),
      analysisName_(iConfig.getParameter<std::string>("analysisName")) {
  structGammaToken_ = consumes<OrbitCollection<l1Scouting::TkEm>>(iConfig.getParameter<edm::InputTag>("srcGamma"));
  structMesonToken_ =
      consumes<OrbitCollection<l1Scouting::RecMeson<2>>>(iConfig.getParameter<edm::InputTag>("srcMeson"));
  produces<std::vector<unsigned>>("selectedBx");
  produces<l1ScoutingRun3::OrbitFlatTable>(analysisName_);
}

ScPhase2BosonToRecMesonGamma::~ScPhase2BosonToRecMesonGamma() {};

void ScPhase2BosonToRecMesonGamma::beginStream(edm::StreamID) {
  countStruct_ = 0;
  passStruct_ = 0;
}

void ScPhase2BosonToRecMesonGamma::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  edm::Handle<OrbitCollection<l1Scouting::TkEm>> srcGamma;
  edm::Handle<OrbitCollection<l1Scouting::RecMeson<2>>> srcMeson;
  iEvent.getByToken(structGammaToken_, srcGamma);
  iEvent.getByToken(structMesonToken_, srcMeson);
  runObj(*srcGamma, *srcMeson, iEvent, countStruct_, passStruct_, "");
}

void ScPhase2BosonToRecMesonGamma::endStream() {
  edm::LogImportant("ScPhase2AnalysisSummary")
      << "Rec Meson " << analysisName_ << " Struct analysis: " << countStruct_ << " -> " << passStruct_;
}

template <typename T, typename U>
void ScPhase2BosonToRecMesonGamma::runObj(const OrbitCollection<T> &srcGamma,
                                          const OrbitCollection<U> &srcMeson,
                                          edm::Event &iEvent,
                                          unsigned long &nTry,
                                          unsigned long &nPass,
                                          const std::string &label) {
  l1ScoutingRun3::BxOffsetsFillter bxOffsetsFiller;
  bxOffsetsFiller.start();
  auto ret = std::make_unique<std::vector<unsigned>>();
  std::vector<float> masses;
  std::vector<uint8_t> i0s, i1s, i2s;
  std::array<unsigned int, 2> bestTriplet;
  float mass = 0.0f;
  float bestTripletScore = 0.;
  float bestTripletMass = 0.;
  bool bestTripletFound;

  for (unsigned int bx = 1; bx <= OrbitCollection<T>::NBX; ++bx) {
    nTry++;
    bestTripletFound = false;
    bestTripletScore = 0.;

    auto range = srcGamma.bxIterator(bx);
    auto nGamma = range.size();

    auto rangeMesons = srcMeson.bxIterator(bx);
    unsigned int nMesons = rangeMesons.size();

    if (nGamma < 1 || nMesons < 1)
      continue;

    const T *candsGamma = &range.front();
    const U *candsMeson = &rangeMesons.front();
    for (unsigned int i1 = 0; i1 < nMesons; ++i1) {
      if (candsMeson[i1].pt() < minPtQ_)
        continue;

      for (unsigned int i2 = 0; i2 < nGamma; ++i2) {
        if (candsGamma[i2].pt() < minPtGamma_)
          continue;

        std::array<unsigned int, 2> pair{{i1, i2}};
        auto p4 = tripletP4(pair, candsGamma, candsMeson);
        float mass = p4.M(), pt = p4.Pt();

        if (!(mass >= minMassBoson_ and mass <= maxMassBoson_))
          continue;

        if (pt > bestTripletScore) {
          bestTripletFound = true;
          bestTriplet = pair;
          bestTripletScore = pt;
          bestTripletMass = mass;
        }
      }
    }

    if (!bestTripletFound)
      continue;

    ret->emplace_back(bx);
    nPass++;
    masses.push_back(bestTripletMass);
    i0s.push_back(candsMeson[bestTriplet[0]].daughterIds(0));
    i1s.push_back(candsMeson[bestTriplet[0]].daughterIds(1));
    i2s.push_back(bestTriplet[1]);
    bxOffsetsFiller.addBx(bx, 1);
  }  // loop on BXs

  iEvent.put(std::move(ret), "selectedBx" + label);
  // now we make the table
  auto bxOffsets = bxOffsetsFiller.done();
  auto tab = std::make_unique<l1ScoutingRun3::OrbitFlatTable>(bxOffsets, analysisName_ + label, true);
  tab->addColumn<float>("mass", masses, "2 daughters plus photon invariant mass");
  tab->addColumn<uint8_t>("i0", i0s, "leading daughter (from meson)");
  tab->addColumn<uint8_t>("i1", i1s, "subleading daughter (from meson)");
  tab->addColumn<uint8_t>("i2", i2s, "photon");
  iEvent.put(std::move(tab), analysisName_ + label);
}

template <typename T, typename U>
ROOT::Math::PtEtaPhiMVector ScPhase2BosonToRecMesonGamma::tripletP4(const std::array<unsigned int, 2> &t,
                                                                    const T *candsGamma,
                                                                    const U *candsMeson) {
  ROOT::Math::PtEtaPhiMVector p1(
      candsMeson[t[0]].pt(), candsMeson[t[0]].eta(), candsMeson[t[0]].phi(), candsMeson[t[0]].mass());
  ROOT::Math::PtEtaPhiMVector p2(candsGamma[t[1]].pt(), candsGamma[t[1]].eta(), candsGamma[t[1]].phi(), 0);
  return (p1 + p2);
}

void ScPhase2BosonToRecMesonGamma::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("srcGamma");
  desc.add<edm::InputTag>("srcMeson");
  desc.add<bool>("runStruct", true);
  desc.add<double>("minMassBoson");
  desc.add<double>("maxMassBoson");
  desc.add<double>("minPtQ");
  desc.add<double>("minPtGamma");
  desc.add<std::string>("analysisName");
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScPhase2BosonToRecMesonGamma);
