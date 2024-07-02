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

  std::unique_ptr<OrbitFlatTable> produceSingleObject(OrbitCollection<l1ScoutingRun3::BxSums> const& src) const;
  std::unique_ptr<OrbitFlatTable> produceManyObjects(OrbitCollection<l1ScoutingRun3::BxSums> const& src) const;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  // the tokens to access the data
  edm::EDGetTokenT<OrbitCollection<l1ScoutingRun3::BxSums>> src_;

  std::string name_, doc_;
  bool singleObject_;
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
      singleObject_(iConfig.getParameter<bool>("singleObject")),
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
  auto out = singleObject_ ? produceSingleObject(*src) : produceManyObjects(*src);
  iEvent.put(std::move(out));
}

std::unique_ptr<OrbitFlatTable> ConvertScoutingSumsToOrbitFlatTable::produceSingleObject(
    OrbitCollection<l1ScoutingRun3::BxSums> const& src) const {
  auto out = std::make_unique<OrbitFlatTable>(src.bxOffsets(), name_, /*singleton=*/true);
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
  for (const l1ScoutingRun3::BxSums& sums : src) {
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

  return out;
}

std::unique_ptr<OrbitFlatTable> ConvertScoutingSumsToOrbitFlatTable::produceManyObjects(
    OrbitCollection<l1ScoutingRun3::BxSums> const& src) const {
  unsigned int nitems = 6;  // totalEt, totalEtEm, missEt, totalHt, missHt, towerCount
  if (writeHF_)
    nitems += 2;
  if (writeAsym_)
    nitems += (writeHF_ ? 4 : 2);
  if (writeMinBias_)
    nitems += 4;
  if (writeCentrality_)
    nitems += 1;
  std::vector<unsigned> offsets(src.bxOffsets());
  for (auto& v : offsets)
    v *= nitems;
  auto out = std::make_unique<OrbitFlatTable>(offsets, name_, /*singleton=*/false);
  std::vector<float> pt(out->size()), phi(out->size(), 0);
  std::vector<int> sumType(out->size());
  unsigned int i = 0, n = out->size();
  for (const l1ScoutingRun3::BxSums& sums : src) {
    assert(i + nitems <= n && i % nitems == 0);
    pt[i] = demux::fEt(sums.hwTotalEt());
    sumType[i++] = l1t::EtSum::kTotalEt;
    pt[i] = demux::fEt(sums.hwTotalEtEm());
    sumType[i++] = l1t::EtSum::kTotalEtEm;
    pt[i] = demux::fEt(sums.hwMissEt());
    phi[i] = demux::fPhi(sums.hwMissEtPhi());
    sumType[i++] = l1t::EtSum::kMissingEt;
    pt[i] = demux::fEt(sums.hwTotalHt());
    sumType[i++] = l1t::EtSum::kTotalHt;
    pt[i] = demux::fEt(sums.hwMissHt());
    phi[i] = demux::fPhi(sums.hwMissHtPhi());
    sumType[i++] = l1t::EtSum::kMissingHt;
    if (writeHF_) {
      pt[i] = demux::fEt(sums.hwMissEtHF());
      phi[i] = demux::fPhi(sums.hwMissEtHFPhi());
      sumType[i++] = l1t::EtSum::kMissingEtHF;
      pt[i] = demux::fEt(sums.hwMissHtHF());
      phi[i] = demux::fPhi(sums.hwMissHtHFPhi());
      sumType[i++] = l1t::EtSum::kMissingHtHF;
    }
    if (writeAsym_) {
      pt[i] = demux::fEt(sums.hwAsymEt());
      sumType[i++] = l1t::EtSum::kAsymEt;
      pt[i] = demux::fEt(sums.hwAsymHt());
      sumType[i++] = l1t::EtSum::kAsymHt;
      if (writeHF_) {
        pt[i] = demux::fEt(sums.hwAsymEtHF());
        sumType[i++] = l1t::EtSum::kAsymEtHF;
        pt[i] = demux::fEt(sums.hwAsymHtHF());
        sumType[i++] = l1t::EtSum::kAsymHtHF;
      }
    }
    if (writeMinBias_) {
      pt[i] = sums.minBiasHFP0();
      sumType[i++] = l1t::EtSum::kMinBiasHFP0;
      pt[i] = sums.minBiasHFM0();
      sumType[i++] = l1t::EtSum::kMinBiasHFM0;
      pt[i] = sums.minBiasHFP1();
      sumType[i++] = l1t::EtSum::kMinBiasHFP1;
      pt[i] = sums.minBiasHFM1();
      sumType[i++] = l1t::EtSum::kMinBiasHFM1;
    }
    pt[i] = sums.towerCount();
    sumType[i++] = l1t::EtSum::kTowerCount;
    if (writeCentrality_) {
      pt[i] = sums.centrality();
      sumType[i++] = l1t::EtSum::kCentrality;
    }
  }
  out->addColumn<float>("pt", pt, "pt (GeV)", 8);
  out->addColumn<float>("phi", phi, "phi (rad)", 8);
  out->addColumn<int>("etSumType",
                      sumType,
                      "the type of the ET Sum "
                      "(https://github.com/cms-sw/cmssw/blob/master/DataFormats/L1Trigger/interface/EtSum.h#L27-L56)");
  out->setDoc(doc_);
  return out;
}

void ConvertScoutingSumsToOrbitFlatTable::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("src");
  desc.add<std::string>("name");
  desc.add<std::string>("doc");
  desc.add<bool>("singleObject", true);
  desc.add<bool>("writeHF", true);
  desc.add<bool>("writeMinBias", true);
  desc.add<bool>("writeCentrality", true);
  desc.add<bool>("writeAsym", true);

  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ConvertScoutingSumsToOrbitFlatTable);