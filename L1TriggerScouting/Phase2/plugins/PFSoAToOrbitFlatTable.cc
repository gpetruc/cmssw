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

#include "DataFormats/NanoAOD/interface/OrbitFlatTable.h"
#include "DataFormats/L1ScoutingSoA/interface/BxLookupHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/PFCandidateHostCollection.h"


class PFSoAToOrbitFlatTable : public edm::global::EDProducer<> {
public:
  // constructor and destructor
  explicit PFSoAToOrbitFlatTable(const edm::ParameterSet&);
  ~PFSoAToOrbitFlatTable() override {};

  void produce(edm::StreamID, edm::Event&, edm::EventSetup const&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  // the tokens to access the data
  edm::EDGetTokenT<l1sc::BxLookupHostCollection> srcBx_;
  edm::EDGetTokenT<l1sc::PFCandidateHostCollection> srcPF_;

  std::string name_, doc_;
};
// -----------------------------------------------------------------------------

// -------------------------------- constructor  -------------------------------

PFSoAToOrbitFlatTable::PFSoAToOrbitFlatTable(const edm::ParameterSet& iConfig)
    : srcBx_(consumes<l1sc::BxLookupHostCollection>(iConfig.getParameter<edm::InputTag>("srcBx"))),
      srcPF_(consumes<l1sc::PFCandidateHostCollection>(iConfig.getParameter<edm::InputTag>("srcPF"))),
      name_(iConfig.getParameter<std::string>("name")),
      doc_(iConfig.getParameter<std::string>("doc")) {
  produces<l1ScoutingRun3::OrbitFlatTable>();
}
// -----------------------------------------------------------------------------

// ----------------------- method called for each orbit  -----------------------
void PFSoAToOrbitFlatTable::produce(edm::StreamID, edm::Event& iEvent, edm::EventSetup const&) const {
  edm::Handle<l1sc::BxLookupHostCollection> srcBx;
  iEvent.getByToken(srcBx_, srcBx);
  edm::Handle<l1sc::PFCandidateHostCollection> srcPF;
  iEvent.getByToken(srcPF_, srcPF);

  const auto *offsets = srcBx->const_view<l1sc::OffsetsSoA>().offsets().data();
  const unsigned int noff = srcBx->const_view<l1sc::OffsetsSoA>().metadata().size();
  std::vector<unsigned int> bxOffsets{0u};
  bxOffsets.insert(bxOffsets.end(), offsets, offsets + noff);

  const auto *pt = srcPF->const_view().pt().data();
  const auto *eta = srcPF->const_view().eta().data();
  const auto *phi = srcPF->const_view().phi().data();
  const auto *z0 = srcPF->const_view().z0().data();
  const auto *dxy = srcPF->const_view().dxy().data();
  const auto *puppiw = srcPF->const_view().puppiw().data();
  const auto *quality = srcPF->const_view().quality().data();
  const auto *pdgid = srcPF->const_view().pdgid().data();

  const unsigned int npf = srcPF->const_view().metadata().size();

  std::vector<float> pf_pt{pt, pt + npf};
  std::vector<float> pf_eta{eta, eta + npf};
  std::vector<float> pf_phi{phi, phi + npf};
  std::vector<float> pf_z0{z0, z0 + npf};
  std::vector<float> pf_dxy{dxy, dxy + npf};
  std::vector<float> pf_puppiw{puppiw, puppiw + npf};
  std::vector<uint8_t> pf_quality{quality, quality + npf};
  std::vector<int16_t> pf_pdgid{pdgid, pdgid + npf};
  
  auto out = std::make_unique<l1ScoutingRun3::OrbitFlatTable>(bxOffsets, name_);
  out->setDoc(doc_);
  out->addColumn<float>("pt", pf_pt, "L1PF pt");
  out->addColumn<float>("eta", pf_eta, "L1PF eta");
  out->addColumn<float>("phi", pf_phi, "L1PF phi");
  out->addColumn<float>("z0", pf_z0, "L1PF z0");
  out->addColumn<float>("dxy", pf_dxy, "L1PF dxy");
  out->addColumn<float>("puppiw", pf_puppiw, "L1PF puppiw");
  out->addColumn<uint8_t>("quality", pf_quality, "L1PF quality");
  out->addColumn<int16_t>("pdgid", pf_pdgid, "L1PF pdgid");
  iEvent.put(std::move(out));
}

void PFSoAToOrbitFlatTable::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("srcBx");
  desc.add<edm::InputTag>("srcPF");
  desc.add<std::string>("name");
  desc.add<std::string>("doc", "");
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(PFSoAToOrbitFlatTable);
