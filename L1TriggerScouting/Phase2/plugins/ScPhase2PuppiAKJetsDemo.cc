#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/L1Scouting/interface/OrbitCollection.h"
#include "DataFormats/NanoAOD/interface/OrbitFlatTable.h"
#include "DataFormats/L1TParticleFlow/interface/L1ScoutingPuppi.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "L1TriggerScouting/Utilities/interface/BxOffsetsFiller.h"
#include "fastjet/ClusterSequence.hh"

#include <ROOT/RVec.hxx>
#include <Math/Vector4D.h>
#include <Math/GenVector/LorentzVector.h>
#include <Math/GenVector/PtEtaPhiM4D.h>
#include <algorithm>
#include <array>
#include <iostream>

class ScPhase2PuppiAKJetsDemo : public edm::stream::EDProducer<> {
public:
  explicit ScPhase2PuppiAKJetsDemo(const edm::ParameterSet &);
  ~ScPhase2PuppiAKJetsDemo() override;
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  void beginStream(edm::StreamID) override;
  void produce(edm::Event &, const edm::EventSetup &) override;
  void endStream() override;

  edm::EDGetTokenT<OrbitCollection<l1Scouting::Puppi>> src_;
  double R_;
};

ScPhase2PuppiAKJetsDemo::ScPhase2PuppiAKJetsDemo(const edm::ParameterSet &iConfig)
    : src_(consumes<OrbitCollection<l1Scouting::Puppi>>(iConfig.getParameter<edm::InputTag>("src"))),
      R_(iConfig.getParameter<double>("rParam")) {
  produces<l1ScoutingRun3::OrbitFlatTable>("jets");
  produces<l1ScoutingRun3::OrbitFlatTable>("clusters");
}

ScPhase2PuppiAKJetsDemo::~ScPhase2PuppiAKJetsDemo() {};

void ScPhase2PuppiAKJetsDemo::beginStream(edm::StreamID) {}

void ScPhase2PuppiAKJetsDemo::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  using namespace fastjet;

  l1ScoutingRun3::BxOffsetsFillter bxOffsetsFiller;
  bxOffsetsFiller.start();

  edm::Handle<OrbitCollection<l1Scouting::Puppi>> src;
  iEvent.getByToken(src_, src);

  // containers for output products
  // jet 3 momenta
  std::vector<float> pt;
  std::vector<float> eta;
  std::vector<float> phi;
  // clustering information
  std::vector<int> cluster(src->size(), -1);

  // choose a jet definition
  JetDefinition jet_def(antikt_algorithm, R_);
  std::vector<PseudoJet> particles;
  for (unsigned int bx = 1; bx <= OrbitCollection<l1Scouting::Puppi>::NBX; ++bx) {
    auto span = src->bxIterator(bx);
    unsigned int i0 = std::distance(&*src->begin(), &*span.begin()); // global index

    particles.clear();
    for (unsigned i = 0; i < span.size(); i++) {
      auto p4 = span[i].p4();
      particles.emplace_back(p4.px(), p4.py(), p4.pz(), p4.energy());
      particles.back().set_user_index(i);
    }

    // run the clustering, extract the jets
    ClusterSequence cs(particles, jet_def);
    auto jets = sorted_by_pt(cs.inclusive_jets());

    unsigned int icluster = 0;
    for (auto j : jets) {
      pt.push_back(j.pt());
      eta.push_back(j.eta());
      phi.push_back(reco::reducePhiRange(j.phi()));
      for (const auto & dau : j.constituents()) {
        unsigned int idx = dau.user_index() + i0;
        cluster[idx] = icluster;
      } 
      icluster++;
    }

    // add the number of jets to the BX offset
    bxOffsetsFiller.addBx(bx, jets.size());
  }  // loop over BX

  auto bxOffsets = bxOffsetsFiller.done();
  auto tab = std::make_unique<l1ScoutingRun3::OrbitFlatTable>(bxOffsets, "AK4Jets");
  tab->addColumn<float>("pt", pt, "Jet pt");
  tab->addColumn<float>("eta", eta, "Jet eta");
  tab->addColumn<float>("phi", phi, "Jet phi");
  iEvent.put(std::move(tab), "jets");

  auto cltab = std::make_unique<l1ScoutingRun3::OrbitFlatTable>(src->bxOffsets(), "AK4Clusters");
  cltab->addColumn<int>("cluster", cluster, "cluster index (-1 if unclustered)");
  iEvent.put(std::move(cltab), "clusters");
}

void ScPhase2PuppiAKJetsDemo::endStream() {}

void ScPhase2PuppiAKJetsDemo::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src");
  desc.add<double>("rParam", 0.4);
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScPhase2PuppiAKJetsDemo);