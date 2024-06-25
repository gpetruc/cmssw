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
#include "DataFormats/L1Scouting/interface/L1ScoutingMuon.h"
#include "DataFormats/L1Scouting/interface/L1ScoutingCalo.h"
#include "L1TriggerScouting/Utilities/interface/convertToL1TFormat.h"
#include "DataFormats/L1Scouting/interface/OrbitFlatTable.h"

using namespace l1ScoutingRun3;

class ConvertScoutingSumsToOrbitFlatTable : public edm::global::EDProducer<> {
public:
  // constructor and destructor
  explicit ConvertScoutingSumsToOrbitFlatTable(const edm::ParameterSet&);
  ~ConvertScoutingSumsToOrbitFlatTable() override{};

  void produce(edm::StreamID, edm::Event&, edm::EventSetup const&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  // the tokens to access the data
  edm::EDGetTokenT<OrbitCollection<l1ScoutingRun3::BxSums>> src_;

  std::string name_, doc_;
  bool writeHF_;
  bool writeMinBias_;
  bool writeCentrality_;
  bool writeAsym_;
};
// -----------------------------------------------------------------------------

// -------------------------------- constructor  -------------------------------

ConvertScoutingSumsToOrbitFlatTable::ConvertScoutingSumsToOrbitFlatTable(const edm::ParameterSet& iConfig)
    : src_(consumes<OrbitCollection<l1ScoutingRun3::BxSums>>(iConfig.getParameter<edm::InputTag>("src"))),
      name_(iConfig.getParameter<std::string>("name")),
      doc_(iConfig.getParameter<std::string>("doc")),
      writeHF_(iConfig.getParameter<bool>("writeHF")),
      writeMinBias_(iConfig.getParameter<bool>("writeMinBias")),
      writeCentrality_(iConfig.getParameter<bool>("writeCentrality")),
      writeAsym_(iConfig.getParameter<bool>("writeAsym")) {
  produces<OrbitFlatTable>();
}
// -----------------------------------------------------------------------------

// ----------------------- method called for each orbit  -----------------------
void ConvertScoutingSumsToOrbitFlatTable::produce(edm::StreamID, edm::Event& iEvent, edm::EventSetup const&) const {
  edm::Handle<OrbitCollection<l1ScoutingRun3::BxSums>> src;
  iEvent.getByToken(src_, src);
  auto out = std::make_unique<OrbitFlatTable>(src->bxOffsets(), name_, /*singleton=*/true);
  out->setDoc(doc_);

  std::vector<float> totalEt(out->size());
  std::vector<float> totalEtEm(out->size());
  std::vector<float> missEt(out->size());
  std::vector<float> missEtPhi(out->size());
  std::vector<float> missEtHF(out->size());
  std::vector<float> missEtHFPhi(out->size());
  std::vector<float> totalHt(out->size());
  std::vector<float> missHt(out->size());
  std::vector<float> missHtPhi(out->size());
  std::vector<float> missHtHF(out->size());
  std::vector<float> missHtHFPhi(out->size());
  std::vector<float> asymEt(out->size());
  std::vector<float> asymHt(out->size());
  std::vector<float> asymEtHF(out->size());
  std::vector<float> asymHtHF(out->size());
  std::vector<int> minBiasHFP0(out->size());
  std::vector<int> minBiasHFM0(out->size());
  std::vector<int> minBiasHFP1(out->size());
  std::vector<int> minBiasHFM1(out->size());
  std::vector<int> centrality(out->size());
  std::vector<int> towerCount(out->size());

  unsigned int i = 0;
  for (const l1ScoutingRun3::BxSums& sums : *src) {
    totalEt[i] = demux::fEt(sums.hwTotalEt());
    totalEtEm[i] = demux::fEt(sums.hwTotalEtEm());
    missEt[i] = demux::fEt(sums.hwMissEt());
    missEtPhi[i] = demux::fPhi(sums.hwMissEtPhi());
    missEtHF[i] = demux::fEt(sums.hwMissEtHF());
    missEtHFPhi[i] = demux::fPhi(sums.hwMissEtHFPhi());
    totalHt[i] = demux::fEt(sums.hwTotalHt());
    missHt[i] = demux::fEt(sums.hwMissHt());
    missHtPhi[i] = demux::fPhi(sums.hwMissHtPhi());
    missHtHF[i] = demux::fEt(sums.hwMissHtHF());
    missHtHFPhi[i] = demux::fPhi(sums.hwMissHtHFPhi());
    asymEt[i] = demux::fEt(sums.hwAsymEt());
    asymHt[i] = demux::fEt(sums.hwAsymHt());
    asymEtHF[i] = demux::fEt(sums.hwAsymEtHF());
    asymHtHF[i] = demux::fEt(sums.hwAsymHtHF());
    minBiasHFP0[i] = sums.minBiasHFP0();
    minBiasHFM0[i] = sums.minBiasHFM0();
    minBiasHFP1[i] = sums.minBiasHFP1();
    minBiasHFM1[i] = sums.minBiasHFM1();
    towerCount[i] = sums.towerCount();
    centrality[i] = sums.centrality();
    ++i;
  }

  out->template addColumn<float>("totalEt", totalEt, "totalEt from Calo (GeV)");
  out->template addColumn<float>("totalEtEm", totalEtEm, "totalEtEm from Calo (GeV)");
  out->template addColumn<float>("missEt", missEt, "missEt from Calo (GeV)");
  out->template addColumn<float>("missEtPhi", missEtPhi, "missEtPhi from Calo");
  out->template addColumn<float>("totalHt", totalHt, "totalHt from Calo (GeV)");
  out->template addColumn<float>("missHt", missHt, "missHt from Calo (GeV)");
  out->template addColumn<float>("missHtPhi", missHtPhi, "missHtPhi from Calo");
  if (writeHF_) {
    out->template addColumn<float>("missEtHF", missEtHF, "missEtHF from Calo (GeV)");
    out->template addColumn<float>("missEtHFPhi", missEtHFPhi, "missEtHFPhi from Calo");
    out->template addColumn<float>("missHtHF", missHtHF, "missHtHF from Calo (GeV)");
    out->template addColumn<float>("missHtHFPhi", missHtHFPhi, "missHtHFPhi from Calo");
  }
  if (writeAsym_) {
    out->template addColumn<float>("asymEt", asymEt, "asymEt from Calo (GeV)");
    out->template addColumn<float>("asymHt", asymHt, "asymHt from Calo (GeV)");
    if (writeHF_) {
      out->template addColumn<float>("asymEtHF", asymEtHF, "asymEtHF from Calo (GeV)");
      out->template addColumn<float>("asymHtHF", asymHtHF, "asymHtHF from Calo (GeV)");
    }
  }
  if (writeMinBias_) {
    out->template addColumn<int>("minBiasHFP0", minBiasHFP0, "minBiasHFP0 from Calo");
    out->template addColumn<int>("minBiasHFM0", minBiasHFM0, "minBiasHFM0 from Calo");
    out->template addColumn<int>("minBiasHFP1", minBiasHFP1, "minBiasHFP1 from Calo");
    out->template addColumn<int>("minBiasHFM1", minBiasHFM1, "minBiasHFM1 from Calo");
  }
  if (writeCentrality_) {
    out->template addColumn<int>("centrality", centrality, "centrality from Calo");
  }
  out->template addColumn<int>("towerCount", towerCount, "towerCount from Calo");

  iEvent.put(std::move(out));
}

void ConvertScoutingSumsToOrbitFlatTable::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("src");
  desc.add<std::string>("name");
  desc.add<std::string>("doc");
  desc.add<bool>("writeHF", true);
  desc.add<bool>("writeMinBias", true);
  desc.add<bool>("writeCentrality", true);
  desc.add<bool>("writeAsym", true);

  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ConvertScoutingSumsToOrbitFlatTable);