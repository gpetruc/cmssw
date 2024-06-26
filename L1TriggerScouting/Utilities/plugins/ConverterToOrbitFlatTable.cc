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

template <typename T>
class ConverterToOrbitFlatTable : public edm::global::EDProducer<> {
public:
  // constructor and destructor
  explicit ConverterToOrbitFlatTable(const edm::ParameterSet&);
  ~ConverterToOrbitFlatTable() override{};

  void produce(edm::StreamID, edm::Event&, edm::EventSetup const&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  // the tokens to access the data
  edm::EDGetTokenT<OrbitCollection<T>> src_;

  std::string name_, doc_;
};
// -----------------------------------------------------------------------------

// -------------------------------- constructor  -------------------------------

template <typename T>
ConverterToOrbitFlatTable<T>::ConverterToOrbitFlatTable(const edm::ParameterSet& iConfig)
    : src_(consumes<OrbitCollection<T>>(iConfig.getParameter<edm::InputTag>("src"))),
      name_(iConfig.getParameter<std::string>("name")),
      doc_(iConfig.getParameter<std::string>("doc")) {
  produces<OrbitFlatTable>();
}
// -----------------------------------------------------------------------------

// ----------------------- method called for each orbit  -----------------------
template <typename T>
void ConverterToOrbitFlatTable<T>::produce(edm::StreamID, edm::Event& iEvent, edm::EventSetup const&) const {
  edm::Handle<OrbitCollection<T>> src;
  iEvent.getByToken(src_, src);
  auto out = std::make_unique<OrbitFlatTable>(src->bxOffsets(), name_);
  out->setDoc(doc_);
  std::vector<float> pt(out->size()), eta(out->size()), phi(out->size());
  if constexpr (std::is_same<T, l1ScoutingRun3::Muon>()) {
    std::vector<float> ptUnconstrained(out->size()), etaAtVtx(out->size()), phiAtVtx(out->size());
    std::vector<int> charge(out->size()), quality(out->size()), dxy(out->size()), index(out->size());
    unsigned int i = 0;
    for (const l1ScoutingRun3::Muon& muon : *src) {
      pt[i] = ugmt::fPt(muon.hwPt());
      eta[i] = ugmt::fEta(muon.hwEta());
      phi[i] = ugmt::fPhi(muon.hwPhi());
      charge[i] = muon.hwCharge();
      quality[i] = muon.hwQual();
      dxy[i] = muon.hwDXY();
      index[i] = muon.tfMuonIndex();
      ptUnconstrained[i] = ugmt::fPtUnconstrained(muon.hwPtUnconstrained());
      etaAtVtx[i] = ugmt::fEtaAtVtx(muon.hwEtaAtVtx());
      phiAtVtx[i] = ugmt::fPhiAtVtx(muon.hwPhiAtVtx());
      ++i;
    }
    out->template addColumn<int>("hwCharge", charge, "charge (0 = not valid)");
    out->template addColumn<int>("hwQual", quality, "quality");
    out->template addColumn<int>("hwDXY", dxy, "hwDXY");
    out->template addColumn<int>("tfMuonIndex",
                                 index,
                                 "Index of muon at the uGMT input. 3 indices per link/sector/wedge. EMTF+ are 0-17, "
                                 "OMTF+ are 18-35, BMTF are 36-71, OMTF- are 72-89, EMTF- are 90-107");
    out->template addColumn<float>("ptUnconstrained", ptUnconstrained, "Unconstrained p_{T} (GeV)");
    out->template addColumn<float>("etaAtVtx", etaAtVtx, "eta extrapolated at beam line (natural units)");
    out->template addColumn<float>("phiAtVtx", phiAtVtx, "phi extrapolated at beam line (natural units)");
  } else {
    std::vector<int> isolation(out->size());
    unsigned int i = 0;
    for (const T& calo : *src) {
      pt[i] = demux::fEt(calo.hwEt());
      eta[i] = demux::fEta(calo.hwEta());
      phi[i] = demux::fPhi(calo.hwPhi());
      isolation[i] = calo.hwIso();
      ++i;
    }
    out->template addColumn<int>("hwIso", isolation, "isolation (trigger units)");
  }
  out->template addColumn<float>("pt", pt, "pt (GeV)");
  out->template addColumn<float>("eta", eta, "eta (natural units)");
  out->template addColumn<float>("phi", phi, "phi (natural units)");

  iEvent.put(std::move(out));
}

template <typename T>
void ConverterToOrbitFlatTable<T>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("src");
  desc.add<std::string>("name");
  desc.add<std::string>("doc");

  descriptions.addDefault(desc);
}

typedef ConverterToOrbitFlatTable<l1ScoutingRun3::Muon> ConvertScoutingMuonsToOrbitFlatTable;
typedef ConverterToOrbitFlatTable<l1ScoutingRun3::Jet> ConvertScoutingJetsToOrbitFlatTable;
typedef ConverterToOrbitFlatTable<l1ScoutingRun3::EGamma> ConvertScoutingEGammasToOrbitFlatTable;
typedef ConverterToOrbitFlatTable<l1ScoutingRun3::Tau> ConvertScoutingTausToOrbitFlatTable;
DEFINE_FWK_MODULE(ConvertScoutingMuonsToOrbitFlatTable);
DEFINE_FWK_MODULE(ConvertScoutingJetsToOrbitFlatTable);
DEFINE_FWK_MODULE(ConvertScoutingEGammasToOrbitFlatTable);
DEFINE_FWK_MODULE(ConvertScoutingTausToOrbitFlatTable);
