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
#include "DataFormats/L1TParticleFlow/interface/L1ScoutingTkEm.h"
#include "DataFormats/L1Scouting/interface/OrbitFlatTable.h"

class ScTkEleToOrbitFlatTable : public edm::global::EDProducer<> {
public:
  // constructor and destructor
  explicit ScTkEleToOrbitFlatTable(const edm::ParameterSet&);
  ~ScTkEleToOrbitFlatTable() override{};

  void produce(edm::StreamID, edm::Event&, edm::EventSetup const&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  // the tokens to access the data
  edm::EDGetTokenT<OrbitCollection<l1Scouting::TkEle>> src_;

  std::string name_, doc_;
};
// -----------------------------------------------------------------------------

// -------------------------------- constructor  -------------------------------

ScTkEleToOrbitFlatTable::ScTkEleToOrbitFlatTable(const edm::ParameterSet& iConfig)
    : src_(consumes<OrbitCollection<l1Scouting::TkEle>>(iConfig.getParameter<edm::InputTag>("src"))),
      name_(iConfig.getParameter<std::string>("name")),
      doc_(iConfig.getParameter<std::string>("doc")) {
  produces<l1ScoutingRun3::OrbitFlatTable>();
}
// -----------------------------------------------------------------------------

// ----------------------- method called for each orbit  -----------------------
void ScTkEleToOrbitFlatTable::produce(edm::StreamID, edm::Event& iEvent, edm::EventSetup const&) const {
  edm::Handle<OrbitCollection<l1Scouting::TkEle>> src;
  iEvent.getByToken(src_, src);
  auto out = std::make_unique<l1ScoutingRun3::OrbitFlatTable>(src->bxOffsets(), name_);
  out->setDoc(doc_);
  std::vector<float> pt(out->size()), eta(out->size()), phi(out->size()), isolation(out->size()), z0(out->size());
  std::vector<uint8_t> quality(out->size());
  std::vector<int16_t> charge(out->size());
  unsigned int i = 0;
  for (const l1Scouting::TkEle& tkem : *src) {
    pt[i] = tkem.pt();
    eta[i] = tkem.eta();
    phi[i] = tkem.phi();
    quality[i] = tkem.quality();
    isolation[i] = tkem.isolation();
    charge[i] = tkem.charge();
    z0[i] = tkem.z0();
    ++i;
  }
  out->addColumn<float>("pt", pt, "pt (GeV)");
  out->addColumn<float>("eta", eta, "eta (natural units)");
  out->addColumn<float>("phi", phi, "phi (natural units)");
  out->addColumn<uint8_t>("quality", quality, "quality (8 bits");
  out->addColumn<float>("isolation", isolation, "isolation (natural units)");
  out->addColumn<int16_t>("charge", charge, "charge");
  out->addColumn<float>("z0", z0, "z0 (cm)");
  iEvent.put(std::move(out));
}

void ScTkEleToOrbitFlatTable::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("src");
  desc.add<std::string>("name");
  desc.add<std::string>("doc");

  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScTkEleToOrbitFlatTable);
