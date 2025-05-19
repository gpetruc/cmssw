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

class ScPhase2TkEmDarkPhotonDiEle : public edm::stream::EDProducer<> {
public:
  explicit ScPhase2TkEmDarkPhotonDiEle(const edm::ParameterSet &);
  ~ScPhase2TkEmDarkPhotonDiEle() override;
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
  edm::EDGetTokenT<OrbitCollection<l1Scouting::TkEle>> structTkEmToken_;

  struct Cuts {
    float minpt = 1;
    float maxeta = 1.479;
    float maxdz = 1;
  } cuts;

  template <typename T>
  static float pairmass(const std::array<unsigned int, 2> &t, const T *cands, const std::array<float, 2> &massD);

  unsigned long countStruct_;
  unsigned long passStruct_;
};

ScPhase2TkEmDarkPhotonDiEle::ScPhase2TkEmDarkPhotonDiEle(const edm::ParameterSet &iConfig)
    : doStruct_(iConfig.getParameter<bool>("runStruct")) {
  if (doStruct_) {
    structTkEmToken_ = consumes<OrbitCollection<l1Scouting::TkEle>>(iConfig.getParameter<edm::InputTag>("srcTkEm"));
    produces<std::vector<unsigned>>("selectedBx");
    produces<l1ScoutingRun3::OrbitFlatTable>("zdee");
  }
}

ScPhase2TkEmDarkPhotonDiEle::~ScPhase2TkEmDarkPhotonDiEle(){};

void ScPhase2TkEmDarkPhotonDiEle::beginStream(edm::StreamID) {
  countStruct_ = 0;
  passStruct_ = 0;
}

void ScPhase2TkEmDarkPhotonDiEle::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  if (doStruct_) {
    //edm::Handle<OrbitCollection<l1Scouting::TkEm>> srcTkEm;
    edm::Handle<OrbitCollection<l1Scouting::TkEle>> srcTkEm;
    iEvent.getByToken(structTkEmToken_, srcTkEm);

    runObj(*srcTkEm, iEvent, countStruct_, passStruct_, "");
  }
}

void ScPhase2TkEmDarkPhotonDiEle::endStream() {
  if (doStruct_)
    edm::LogImportant("ScPhase2AnalysisSummary")
        << "zdee Struct analysis: " << countStruct_ << " -> " << passStruct_;
}

template <typename T>
void ScPhase2TkEmDarkPhotonDiEle::runObj(const OrbitCollection<T> &srcTkEm,
                                        edm::Event &iEvent,
                                        unsigned long &nTry,
                                        unsigned long &nPass,
                                        const std::string &label) {
  l1ScoutingRun3::BxOffsetsFillter bxOffsetsFiller;
  bxOffsetsFiller.start();
  auto ret = std::make_unique<std::vector<unsigned>>();
  
  std::vector<float> masses;
  ROOT::RVec<unsigned int> iEle;
  std::array<unsigned int, 2> bestPair;
  
  bool bestPairFound;
  float maxDeltaPhi;
  
  for (unsigned int bx = 1; bx <= OrbitCollection<T>::NBX; ++bx) {
    nTry++;
    auto range = srcTkEm.bxIterator(bx);
    const T *cands = &range.front();
    auto size = range.size();

    // Select events with two or more electrons with pT > 5 GeV and in barrel
    iEle.clear();
    for (unsigned int i = 0; i < size; ++i) {  //make list of all electrons
      if ((cands[i].pt() >= cuts.minpt) && (std::abs(cands[i].eta()) <= cuts.maxeta)) {
        iEle.push_back(i);
      }
    }
    if (iEle.size() < 2) continue;

    // Loop over possible ee pairs; get the best pair
    bestPairFound = false;
    maxDeltaPhi = -999;
    for (unsigned int i1 = 0; i1 < iEle.size(); ++i1) {
      if ((cands[iEle[i1]].pt() < cuts.minpt) || (std::abs(cands[iEle[i1]].eta()) > cuts.maxeta)) continue;
      
      for (unsigned int i2 = 0; i2 < iEle.size(); ++i2) {
        if ((cands[iEle[i2]].pt() < cuts.minpt) || (std::abs(cands[iEle[i2]].eta()) > cuts.maxeta)) continue;

	// OS requirement
        if (!(cands[iEle[i1]].charge() * cands[iEle[i2]].charge() < 0)) continue;

	// dz requirement
        if (std::abs(cands[iEle[i1]].z0() - cands[iEle[i2]].z0()) > 1) continue;

	// Find the one with the max dPhi
	auto dPhi = ROOT::VecOps::DeltaPhi<float>(cands[iEle[i1]].phi(), cands[iEle[i2]].phi());

        std::array<unsigned int, 2> pair{{iEle[i1], iEle[i2]}};  // pair of indices
        if (dPhi > maxDeltaPhi) {
          std::copy_n(pair.begin(), 2, bestPair.begin());
          maxDeltaPhi = dPhi;
          bestPairFound = true;
        }
      }
    }
    if (!bestPairFound)
      continue;

    // Best ee pair mass
    auto mass = pairmass({{bestPair[0], bestPair[1]}}, cands, {{0.51*1e-3, 0.51*1e-3}});

    ret->emplace_back(bx);
    nPass++;
    
    masses.push_back(mass);
    bxOffsetsFiller.addBx(bx, 1);
  }  // loop on BXs

  iEvent.put(std::move(ret), "selectedBx" + label);
  // now we make the table
  auto bxOffsets = bxOffsetsFiller.done();
  auto tab = std::make_unique<l1ScoutingRun3::OrbitFlatTable>(bxOffsets, "Zdee" + label, true);
  
  tab->addColumn<float>("mass", masses, "2 kaons plus photon invariant mass");
  
  iEvent.put(std::move(tab), "zdee" + label);
}


template <typename T>
float ScPhase2TkEmDarkPhotonDiEle::pairmass(const std::array<unsigned int, 2> &t,
                                           const T *cands,
                                           const std::array<float, 2> &massD) {
  ROOT::Math::PtEtaPhiMVector p1(cands[t[0]].pt(), cands[t[0]].eta(), cands[t[0]].phi(), massD[0]);
  ROOT::Math::PtEtaPhiMVector p2(cands[t[1]].pt(), cands[t[1]].eta(), cands[t[1]].phi(), massD[1]);
  float mass = (p1 + p2).M();
  return mass;
}

void ScPhase2TkEmDarkPhotonDiEle::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("srcTkEm");
  desc.add<bool>("runStruct", true);
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScPhase2TkEmDarkPhotonDiEle);
