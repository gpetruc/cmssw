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
#include "DataFormats/L1TMuonPhase2/interface/L1ScoutingTrackerMuon.h"
#include "DataFormats/L1Scouting/interface/OrbitFlatTable.h"

class ScTrackerMuonToOrbitFlatTable : public edm::global::EDProducer<> {
public:
  // constructor and destructor
  explicit ScTrackerMuonToOrbitFlatTable(const edm::ParameterSet&);
  ~ScTrackerMuonToOrbitFlatTable() override{};

  void produce(edm::StreamID, edm::Event&, edm::EventSetup const&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  // the tokens to access the data
  edm::EDGetTokenT<OrbitCollection<l1Scouting::TrackerMuon>> src_;

  std::string name_, doc_;
};
// -----------------------------------------------------------------------------

// -------------------------------- constructor  -------------------------------

ScTrackerMuonToOrbitFlatTable::ScTrackerMuonToOrbitFlatTable(const edm::ParameterSet& iConfig)
    : src_(consumes<OrbitCollection<l1Scouting::TrackerMuon>>(iConfig.getParameter<edm::InputTag>("src"))),
      name_(iConfig.getParameter<std::string>("name")),
      doc_(iConfig.getParameter<std::string>("doc")) {
  produces<l1ScoutingRun3::OrbitFlatTable>();
}
// -----------------------------------------------------------------------------

// ----------------------- method called for each orbit  -----------------------
void ScTrackerMuonToOrbitFlatTable::produce(edm::StreamID, edm::Event& iEvent, edm::EventSetup const&) const {
  
  edm::Handle<OrbitCollection<l1Scouting::TrackerMuon>> src;
  iEvent.getByToken(src_, src);
  
  auto out = std::make_unique<l1ScoutingRun3::OrbitFlatTable>(src->bxOffsets(), name_);
  out->setDoc(doc_);
  std::vector<float> pt(out->size()), eta(out->size()), phi(out->size()), z0(out->size()), d0(out->size()),beta(out->size());
  std::vector<uint8_t> quality(out->size()),isolation(out->size());
  std::vector<int32_t>  charge(out->size());
  unsigned int i = 0;
  for (const l1Scouting::TrackerMuon& muon : *src) {
    pt[i]      = muon.pt();
    eta[i]     = muon.eta();
    phi[i]     = muon.phi();
    z0[i]      = muon.z0();
    d0[i]      = muon.d0();
    quality[i] = muon.quality();
    isolation[i] = muon.isolation();
    charge[i]  = muon.charge();
    ++i;
  }

  out->addColumn<float>("pt", pt, "pt (GeV)");
  out->addColumn<float>("eta", eta, "eta (natural units)");
  out->addColumn<float>("phi", phi, "phi (natural units)");
  out->addColumn<float>("z0", z0, "z0 (cm)");
  out->addColumn<float>("d0", d0, "d0 (cm)");
  out->addColumn<int32_t>("charge", charge, "charge (+1/-1)");
  out->addColumn<uint8_t>("isolation", isolation, "isolation (8 bits");
  out->addColumn<uint8_t>("quality",quality, "quality (8 bits");
  iEvent.put(std::move(out));
}

void ScTrackerMuonToOrbitFlatTable::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("src");
  desc.add<std::string>("name");
  desc.add<std::string>("doc");

  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScTrackerMuonToOrbitFlatTable);
