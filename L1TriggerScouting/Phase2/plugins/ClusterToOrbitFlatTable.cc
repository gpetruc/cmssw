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
#include "DataFormats/L1ScoutingSoA/interface/AssociationMapHost.h"
#include "DataFormats/L1ScoutingSoA/interface/BxLookupHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/PFCandidateHostCollection.h"


class ClusterToOrbitFlatTable : public edm::global::EDProducer<> {
public:
  // constructor and destructor
  explicit ClusterToOrbitFlatTable(const edm::ParameterSet&);
  ~ClusterToOrbitFlatTable() override {};

  void produce(edm::StreamID, edm::Event&, edm::EventSetup const&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  // the tokens to access the data
  edm::EDGetTokenT<l1sc::BxLookupHostCollection> srcBxClustersMap_;
  edm::EDGetTokenT<l1sc::AssociationMapHost> srcClusterCandsMap_;
  edm::EDGetTokenT<l1sc::PFCandidateHostCollection> srcCandidates_;

  std::string name_, doc_;
};
// -----------------------------------------------------------------------------

// -------------------------------- constructor  -------------------------------

ClusterToOrbitFlatTable::ClusterToOrbitFlatTable(const edm::ParameterSet& iConfig) :
      srcBxClustersMap_(consumes<l1sc::BxLookupHostCollection>(iConfig.getParameter<edm::InputTag>("srcClusters"))),
      srcClusterCandsMap_(consumes<l1sc::AssociationMapHost>(iConfig.getParameter<edm::InputTag>("srcClusters"))),
      srcCandidates_(consumes<l1sc::PFCandidateHostCollection>(iConfig.getParameter<edm::InputTag>("srcCandidates"))),
      name_(iConfig.getParameter<std::string>("name")),
      doc_(iConfig.getParameter<std::string>("doc")) {
  produces<l1ScoutingRun3::OrbitFlatTable>();
}
// -----------------------------------------------------------------------------

// ----------------------- method called for each orbit  -----------------------
void ClusterToOrbitFlatTable::produce(edm::StreamID, edm::Event& iEvent, edm::EventSetup const&) const {
  edm::Handle<l1sc::BxLookupHostCollection> srcBxClustersMap;
  iEvent.getByToken(srcBxClustersMap_, srcBxClustersMap);
  edm::Handle<l1sc::AssociationMapHost> srcClusterCandsMap;
  iEvent.getByToken(srcClusterCandsMap_, srcClusterCandsMap);
  edm::Handle<l1sc::PFCandidateHostCollection> srcCandidates;
  iEvent.getByToken(srcCandidates_, srcCandidates);

  std::cout << "Serializing clusters constituents for the current orbit" << std::endl;

  const unsigned int nbx = srcBxClustersMap->const_view<l1sc::OffsetsSoA>().metadata().size() - 1;

  const auto *bx_offsets = srcBxClustersMap->const_view<l1sc::OffsetsSoA>().offsets().data();
  const auto *cluster_idx = srcBxClustersMap->const_view<l1sc::BxIndexSoA>().bx().data();
  const auto *cluster_off = srcClusterCandsMap->const_view<l1sc::OffsetsSoA>().offsets().data();
  const auto *candidate_idx = srcClusterCandsMap->const_view<l1sc::IndexSoA>().indexes().data();

  const auto *pt = srcCandidates->const_view().pt().data();
  const auto *eta = srcCandidates->const_view().eta().data();
  const auto *phi = srcCandidates->const_view().phi().data();

  std::vector<float> pts, etas, phis;
  std::vector<uint32_t> clusters;

  // bxOffsets is needed in order to keep track of the offsets
  // for every bx with respect to the candidates that are within it.
  // But none of the offsets inside the CandsClusterBxHostCollection 
  // contain such information
  std::vector<uint32_t> bxOffsets{1u, 0u};

  // loop over the bxs in the orbit
  for (unsigned int bx_idx = 0; bx_idx < nbx; ++bx_idx) {
    auto bx_start = bx_offsets[bx_idx];
    auto bx_end = bx_offsets[bx_idx + 1];
    auto bx_size = bx_end - bx_start;
    
    // loop over the clusters indexes of the clusters in the current bx
    for (unsigned int ii = 0; ii < bx_size; ++ii) {
      auto cl_idx = cluster_idx[bx_start + ii];
      auto cluster_start = cluster_off[cl_idx];
      auto cluster_end = cluster_off[cl_idx + 1];
      auto cluster_size = cluster_end - cluster_start;

      // loop over the candidates indexes of the candidats in the current cluster
      for (unsigned int jj = 0; jj < cluster_size; ++jj) {
        auto cand_idx = candidate_idx[cluster_start + jj];
        pts.push_back(pt[cand_idx]);
        etas.push_back(eta[cand_idx]);
        phis.push_back(phi[cand_idx]);
        clusters.push_back(cl_idx);
      }
    }

    bxOffsets.push_back(pts.size());
  }

  auto out = std::make_unique<l1ScoutingRun3::OrbitFlatTable>(bxOffsets, name_);
  out->setDoc(doc_);
  out->addColumn<float>("pt", pts, "candidate pt (GeV)");
  out->addColumn<float>("eta", etas, "candidate eta (GeV)");
  out->addColumn<float>("phi", phis, "candidate phi (GeV)");
  out->addColumn<uint32_t>("cluster", clusters, "cluster index");
  iEvent.put(std::move(out));
}

void ClusterToOrbitFlatTable::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("srcClusters");
  desc.add<edm::InputTag>("srcCandidates");
  desc.add<std::string>("name");
  desc.add<std::string>("doc");
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ClusterToOrbitFlatTable);
