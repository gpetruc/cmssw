#include "FWCore/Framework/interface/MakerMacros.h"

#include <fstream>
#include <memory>
#include <string>
#include <cmath>
#include <cstdint>

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
#include "DataFormats/L1Scouting/interface/L1ScoutingBMTFStub.h"
#include "DataFormats/L1Scouting/interface/OrbitFlatTable.h"

using namespace l1ScoutingRun3;

class ConvertScoutingStubsToOrbitFlatTable : public edm::global::EDProducer<> {
public:
  // constructor and destructor
  explicit ConvertScoutingStubsToOrbitFlatTable(const edm::ParameterSet&);
  ~ConvertScoutingStubsToOrbitFlatTable() override{};

  void produce(edm::StreamID, edm::Event&, edm::EventSetup const&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  // the tokens to access the data
  edm::EDGetTokenT<OrbitCollection<l1ScoutingRun3::BMTFStub>> src_;

  std::string name_, doc_;
};
// -----------------------------------------------------------------------------

// -------------------------------- constructor  -------------------------------

ConvertScoutingStubsToOrbitFlatTable::ConvertScoutingStubsToOrbitFlatTable(const edm::ParameterSet& iConfig)
    : src_(consumes<OrbitCollection<l1ScoutingRun3::BMTFStub>>(iConfig.getParameter<edm::InputTag>("src"))),
      name_(iConfig.getParameter<std::string>("name")),
      doc_(iConfig.getParameter<std::string>("doc")) {
  produces<OrbitFlatTable>();
}
// -----------------------------------------------------------------------------

// ----------------------- method called for each orbit  -----------------------
void ConvertScoutingStubsToOrbitFlatTable::produce(edm::StreamID, edm::Event& iEvent, edm::EventSetup const&) const {
  edm::Handle<OrbitCollection<l1ScoutingRun3::BMTFStub>> src;
  iEvent.getByToken(src_, src);
  auto out = std::make_unique<OrbitFlatTable>(src->bxOffsets(), name_);
  out->setDoc(doc_);

  std::vector<int16_t> hwPhi(out->size());
  std::vector<int16_t> hwPhiB(out->size());
  std::vector<int16_t> hwQual(out->size());
  std::vector<int16_t> hwEta(out->size());
  std::vector<int16_t> hwQEta(out->size());
  std::vector<int16_t> station(out->size());
  std::vector<int16_t> wheel(out->size());
  std::vector<int16_t> sector(out->size());
  std::vector<int16_t> tag(out->size());

  unsigned int i = 0;
  for (const l1ScoutingRun3::BMTFStub& stub : *src) {
    hwPhi[i] = stub.hwPhi();
    hwPhiB[i] = stub.hwPhiB();
    hwQual[i] = stub.hwQual();
    hwEta[i] = stub.hwEta();
    hwQEta[i] = stub.hwQEta();
    station[i] = stub.station();
    wheel[i] = stub.wheel();
    sector[i] = stub.sector();
    tag[i] = stub.tag();
    ++i;
  }

  out->addColumn<int16_t>("hwPhi", hwPhi, "hwPhi (raw L1T units)");
  out->addColumn<int16_t>("hwPhiB", hwPhiB, "hwPhiB (raw L1T units)");
  out->addColumn<int16_t>("hwQual", hwQual, "hwQual (raw L1T units)");
  out->addColumn<int16_t>("hwEta", hwEta, "hwEta (raw L1T units)");
  out->addColumn<int16_t>("hwQEta", hwQEta, "hwQEta (raw L1T units)");
  out->addColumn<int16_t>("station", station, "station (raw L1T units)");
  out->addColumn<int16_t>("wheel", wheel, "wheel (raw L1T units)");
  out->addColumn<int16_t>("sector", sector, "sector (raw L1T units)");
  out->addColumn<int16_t>("tag", tag, "tag (raw L1T units)");

  iEvent.put(std::move(out));
}

void ConvertScoutingStubsToOrbitFlatTable::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("src");
  desc.add<std::string>("name");
  desc.add<std::string>("doc");

  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ConvertScoutingStubsToOrbitFlatTable);