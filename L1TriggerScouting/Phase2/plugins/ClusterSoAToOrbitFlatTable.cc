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


class ClusterSoAToOrbitFlatTable : public edm::global::EDProducer<> {
public:
  // constructor and destructor
  explicit ClusterSoAToOrbitFlatTable(const edm::ParameterSet&);
  ~ClusterSoAToOrbitFlatTable() override {};

  void produce(edm::StreamID, edm::Event&, edm::EventSetup const&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  // the tokens to access the data
  edm::EDGetTokenT<l1sc::BxLookupHostCollection> srcBx_;
  edm::EDGetTokenT<l1sc::ClustersHostCollection> srcClusters_;

  std::string name_, doc_;
};
// -----------------------------------------------------------------------------

// -------------------------------- constructor  -------------------------------

ClusterSoAToOrbitFlatTable::ClusterSoAToOrbitFlatTable(const edm::ParameterSet& iConfig)
    : srcBx_(consumes<l1sc::BxLookupHostCollection>(iConfig.getParameter<edm::InputTag>("srcBx"))),
      srcClusters_(consumes<l1sc::ClustersHostCollection>(iConfig.getParameter<edm::InputTag>("srcClusters"))),
      name_(iConfig.getParameter<std::string>("name")),
      doc_(iConfig.getParameter<std::string>("doc")) {
  produces<l1ScoutingRun3::OrbitFlatTable>();
}
// -----------------------------------------------------------------------------

// ----------------------- method called for each orbit  -----------------------
void ClusterSoAToOrbitFlatTable::produce(edm::StreamID, edm::Event& iEvent, edm::EventSetup const&) const {
  edm::Handle<l1sc::BxLookupHostCollection> srcBx;
  iEvent.getByToken(srcBx_, srcBx);
  edm::Handle<l1sc::ClustersHostCollection> srcClusters;
  iEvent.getByToken(srcClusters_, srcClusters);

  const auto *bxs = srcBx->const_view<l1sc::OffsetsSoA>().offsets();
  const unsigned int nbx = srcBx->const_view<l1sc::OffsetsSoA>().metadata().size();
  std::vector<unsigned int> bxOffsets;
  bxOffsets.push_back(0);
  bxOffsets.insert(bxOffsets.end(), bxs, bxs + nbx);

  const auto *cluster = srcClusters->const_view().cluster();
  const auto *seed = srcClusters->const_view().is_seed();
  const unsigned int nclusters = srcClusters->const_view().metadata().size();
  std::vector<int32_t> clusters{cluster, cluster + nclusters};
  std::vector<int32_t> is_seed{seed, seed + nclusters};
  
  auto out = std::make_unique<l1ScoutingRun3::OrbitFlatTable>(bxOffsets, name_);
  out->setDoc(doc_);
  out->addColumn<int32_t>("cluster", clusters, "cluster index");
  out->addColumn<int32_t>("is_seed", is_seed, "whether the cell is a seed");
  iEvent.put(std::move(out));
}

void ClusterSoAToOrbitFlatTable::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("srcBx");
  desc.add<edm::InputTag>("srcClusters");
  desc.add<std::string>("name");
  desc.add<std::string>("doc", "");
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ClusterSoAToOrbitFlatTable);
