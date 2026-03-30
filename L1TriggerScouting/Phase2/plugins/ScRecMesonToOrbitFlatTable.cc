#include "FWCore/Framework/interface/MakerMacros.h"

#include <fstream>
#include <iomanip>
#include <memory>
#include <string>
#include <cmath>

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
  edm::EDGetTokenT<OrbitCollection<l1Scouting::RecMeson<2>>> src_;

  std::string name_, doc_;
};
// -----------------------------------------------------------------------------

// -------------------------------- constructor  -------------------------------

ScRecMesonToOrbitFlatTable::ScRecMesonToOrbitFlatTable(const edm::ParameterSet& iConfig)
    : src_(consumes<OrbitCollection<l1Scouting::RecMeson<2>>>(iConfig.getParameter<edm::InputTag>("src"))),
      name_(iConfig.getParameter<std::string>("name")),
      doc_(iConfig.getParameter<std::string>("doc")) {
  produces<l1ScoutingRun3::OrbitFlatTable>();
}
// -----------------------------------------------------------------------------

// ----------------------- method called for each orbit  -----------------------
void ScRecMesonToOrbitFlatTable::produce(edm::StreamID, edm::Event& iEvent, edm::EventSetup const&) const {
  edm::Handle<OrbitCollection<l1Scouting::RecMeson<2>>> src;
  iEvent.getByToken(src_, src);
  auto out = std::make_unique<l1ScoutingRun3::OrbitFlatTable>(src->bxOffsets(), name_);
  out->setDoc(doc_);
  std::vector<float> pt(out->size()), eta(out->size()), phi(out->size()), mass(out->size());
  std::vector<int16_t> id1(out->size()), id2(out->size());
  unsigned int i = 0;
  for (const l1Scouting::RecMeson<2>& RecMeson : *src) {
    pt[i] = RecMeson.pt();
    eta[i] = RecMeson.eta();
    phi[i] = RecMeson.phi();
    id1[i] = RecMeson.daughterIds(0);
    id2[i] = RecMeson.daughterIds(1);
    mass[i] = RecMeson.mass();

    ++i;
  }
  out->addColumn<float>("pt", pt, "pt (GeV)");
  out->addColumn<float>("eta", eta, "eta (natural units)");
  out->addColumn<float>("phi", phi, "phi (natural units)");
  out->addColumn<int16_t>("id1", id1, "index of 1st daughter");
  out->addColumn<int16_t>("id2", id2, "index of 2nd daughter");
  out->addColumn<float>("mass", mass, "mass (natural units)");
  iEvent.put(std::move(out));
}

void ScRecMesonToOrbitFlatTable::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("src");
  desc.add<std::string>("name");
  desc.add<std::string>("doc");

  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScRecMesonToOrbitFlatTable);
