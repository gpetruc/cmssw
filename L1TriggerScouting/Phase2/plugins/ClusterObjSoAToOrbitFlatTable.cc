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
#include "DataFormats/L1ScoutingSoA/interface/ClustersHostCollection.h"


class ClusterObjSoAToOrbitFlatTable : public edm::global::EDProducer<> {
public:
  // constructor and destructor
  explicit ClusterObjSoAToOrbitFlatTable(const edm::ParameterSet&);
  ~ClusterObjSoAToOrbitFlatTable() override {};

  void produce(edm::StreamID, edm::Event&, edm::EventSetup const&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  // the tokens to access the data
  edm::EDGetTokenT<l1sc::BxLookupHostCollection> srcBx_;
  edm::EDGetTokenT<l1sc::ClusterObjHostCollection> srcClusters_;

  std::string name_, doc_;
};
// -----------------------------------------------------------------------------

// -------------------------------- constructor  -------------------------------

ClusterObjSoAToOrbitFlatTable::ClusterObjSoAToOrbitFlatTable(const edm::ParameterSet& iConfig)
    : srcBx_(consumes<l1sc::BxLookupHostCollection>(iConfig.getParameter<edm::InputTag>("srcBx"))),
      srcClusters_(consumes<l1sc::ClusterObjHostCollection>(iConfig.getParameter<edm::InputTag>("srcClusters"))),
      name_(iConfig.getParameter<std::string>("name")),
      doc_(iConfig.getParameter<std::string>("doc")) {
  produces<l1ScoutingRun3::OrbitFlatTable>();
}
// -----------------------------------------------------------------------------

// ----------------------- method called for each orbit  -----------------------
void ClusterObjSoAToOrbitFlatTable::produce(edm::StreamID, edm::Event& iEvent, edm::EventSetup const&) const {
  edm::Handle<l1sc::BxLookupHostCollection> srcBx;
  iEvent.getByToken(srcBx_, srcBx);
  edm::Handle<l1sc::ClusterObjHostCollection> srcClusters;
  iEvent.getByToken(srcClusters_, srcClusters);

  const auto *bxs = srcBx->const_view<l1sc::OffsetsSoA>().offsets();
  const unsigned int nbx = srcBx->const_view<l1sc::OffsetsSoA>().metadata().size();

  const auto *pt = srcClusters->const_view().pt();
  const auto *eta = srcClusters->const_view().eta();
  const auto *phi = srcClusters->const_view().phi();
  const auto *cluster = srcClusters->const_view().cluster();

  std::vector<unsigned int> bxOffsets{1, 0};
  std::vector<float> pts, etas, phis;
  std::vector<int32_t> clusters;
  for (unsigned int bx = 1; bx < nbx; ++bx) {
    for (unsigned int i = bxs[bx - 1]; i < bxs[bx]; ++i) {
        if (pt[i] > 0) {
            pts.push_back(pt[i]);
            etas.push_back(eta[i]);
            phis.push_back(phi[i]);
            clusters.push_back(cluster[i]);
        }
    }
    bxOffsets.push_back(pts.size());
  }
  
  auto out = std::make_unique<l1ScoutingRun3::OrbitFlatTable>(bxOffsets, name_);
  out->setDoc(doc_);
  out->addColumn<float>("pt", pts, "cluster pt (GeV)");
  out->addColumn<float>("eta", etas, "cluster eta (GeV)");
  out->addColumn<float>("phi", phis, "cluster phi (GeV)");
  out->addColumn<int32_t>("cluster", clusters, "cluster index");
  iEvent.put(std::move(out));
}

void ClusterObjSoAToOrbitFlatTable::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("srcBx");
  desc.add<edm::InputTag>("srcClusters");
  desc.add<std::string>("name");
  desc.add<std::string>("doc");
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ClusterObjSoAToOrbitFlatTable);
