#include "FWCore/Framework/interface/MakerMacros.h"

#include <fstream>
#include <iomanip>
#include <memory>
#include <string>
#include <cmath>
#include <ROOT/RVec.hxx>
#include <Math/Vector4D.h>
#include <Math/GenVector/LorentzVector.h>
#include <Math/GenVector/PtEtaPhiM4D.h>

#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/MessageLogger/interface/MessageDrop.h"

#include "DataFormats/L1Scouting/interface/OrbitCollection.h"
#include "DataFormats/L1TParticleFlow/interface/L1ScoutingPuppi.h"
#include "DataFormats/L1Scouting/interface/OrbitFlatTable.h"
#include "DataFormats/L1TParticleFlow/interface/RecMeson.h"

class ScRecMesonToOrbitFlatTable : public edm::global::EDProducer<> {
public:
  // constructor and destructor
  explicit ScRecMesonToOrbitFlatTable(const edm::ParameterSet&);
  ~ScRecMesonToOrbitFlatTable() override {};

  void produce(edm::StreamID, edm::Event&, edm::EventSetup const&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  // the tokens to access the data
  edm::EDGetTokenT<OrbitCollection<l1Scouting::Puppi>> src_;

  std::string name_, doc_;

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
  static float pairmass(const std::array<unsigned int, 2>& t, const std::vector<T>& cands, const std::array<float, 2>& massD);
  
  std::tuple<bool, float> deltar(float eta1, float eta2, float phi1, float phi2) const;

};
// -----------------------------------------------------------------------------

// -------------------------------- constructor  -------------------------------

ScRecMesonToOrbitFlatTable::ScRecMesonToOrbitFlatTable(const edm::ParameterSet& iConfig)
    : src_(consumes<OrbitCollection<l1Scouting::Puppi>>(iConfig.getParameter<edm::InputTag>("src"))),
      name_(iConfig.getParameter<std::string>("name")),
      doc_(iConfig.getParameter<std::string>("doc")) {
  produces<l1ScoutingRun3::OrbitFlatTable>();
}
// -----------------------------------------------------------------------------

// ----------------------- method called for each orbit  -----------------------
void ScRecMesonToOrbitFlatTable::produce(edm::StreamID, edm::Event& iEvent, edm::EventSetup const&) const {
  
  edm::Handle<OrbitCollection<l1Scouting::Puppi>> src;
  iEvent.getByToken(src_, src);
  auto out = std::make_unique<l1ScoutingRun3::OrbitFlatTable>(src->bxOffsets(), name_);
  out->setDoc(doc_);
  std::vector<l1Scouting::RecMeson> phiRecMesons, rhoRecMesons;
  std::vector<float> phiRecMesons_pt, phiRecMesons_eta, phiRecMesons_phi;
  std::vector<unsigned int> phiRecMesons_i1, phiRecMesons_i2;

  unsigned int i = 0;
  std::vector<int> charged_daus_index;
  std::vector<l1Scouting::Puppi> charged_daus;

  for (const l1Scouting::Puppi& puppi : *src) {
    if ((std::abs(puppi.pdgId()) == 211 or std::abs(puppi.pdgId()) == 11)) {
      if (puppi.pt() >= cuts.minptD){
        charged_daus_index.push_back(i);
        charged_daus.push_back(puppi);
      }
    }
    ++i;
  }

  unsigned int ndaus = charged_daus.size();

  for (unsigned int i1 = 0; i1 < ndaus; ++i1) {
    if (charged_daus[i1].pt() < cuts.minptD)
      continue;  // D1 pt cut
    for (unsigned int i2 = 0; i2 < ndaus; ++i2) {
      if (i2 == i1 || charged_daus[i2].pt() < cuts.minptD)
        continue;  // D2 pt cut

      if (!(charged_daus[i1].charge() * charged_daus[i2].charge() < 0))
        continue;

      auto mass2 = pairmass({{i1, i2}}, charged_daus, {{0.4937, 0.4937}});
      if (mass2 >= cuts.minmassQ and mass2 <= cuts.maxmassQ)
        continue;

      auto [drcond, drQ] = deltar(charged_daus[i1].eta(), charged_daus[i2].eta(), charged_daus[i1].phi(), charged_daus[i2].phi());
      if (!drcond)
        continue;  // angular sep of top 2 tracks

      auto phiMeson = l1Scouting::RecMeson(1.1, 2.2, 3.3, i1, i2);
      phiRecMesons.push_back(phiMeson);
    }
  }

  // out->addColumn<float>("pt", phiRecMesons.pt(), "pt (GeV)");
  // out->addColumn<float>("eta", phiRecMesons.eta(), "eta (natural units)");
  // out->addColumn<float>("phi", phiRecMesons.phi(), "phi (natural units)");
  // out->addColumn<float>("i1", phiRecMesons.i1(), "index of 1st daughter");
  // out->addColumn<float>("i2", phiRecMesons.i2(), "index of 2nd daughter");

  for (const auto& meson : phiRecMesons) {
    // Assuming RecMeson has methods to get these properties:
    phiRecMesons_pt.push_back(meson.pt());
    phiRecMesons_eta.push_back(meson.eta());
    phiRecMesons_phi.push_back(meson.phi());
    phiRecMesons_i1.push_back(meson.id1());
    phiRecMesons_i2.push_back(meson.id2());
  }

  out->addColumn<float>("pt", phiRecMesons_pt, "pt (GeV)");
  out->addColumn<float>("eta", phiRecMesons_eta, "eta (natural units)");
  out->addColumn<float>("phi", phiRecMesons_phi, "phi (natural units)");
  out->addColumn<unsigned int>("id1", phiRecMesons_i1, "index of 1st daughter");
  out->addColumn<unsigned int>("id2", phiRecMesons_i2, "index of 2nd daughter");

  iEvent.put(std::move(out));
}

void ScRecMesonToOrbitFlatTable::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("src");
  desc.add<std::string>("name");
  desc.add<std::string>("doc");

  descriptions.addDefault(desc);
}

std::tuple<bool, float> ScRecMesonToOrbitFlatTable::deltar(float eta1, float eta2, float phi1, float phi2) const {
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
float ScRecMesonToOrbitFlatTable::pairmass(const std::array<unsigned int, 2>& t, const std::vector<T>& cands, const std::array<float, 2>& massD) {
  ROOT::Math::PtEtaPhiMVector p1(cands[t[0]].pt(), cands[t[0]].eta(), cands[t[0]].phi(), massD[0]);
  ROOT::Math::PtEtaPhiMVector p2(cands[t[1]].pt(), cands[t[1]].eta(), cands[t[1]].phi(), massD[1]);
  return (p1 + p2).M();
}

DEFINE_FWK_MODULE(ScRecMesonToOrbitFlatTable);
