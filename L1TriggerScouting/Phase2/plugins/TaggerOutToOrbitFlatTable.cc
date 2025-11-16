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
#include "DataFormats/L1ScoutingSoA/interface/CandsClusterBxHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/PFCandidateHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/SoftTauHostTensor.h"

#include <ROOT/RVec.hxx>
#include <Math/Vector4D.h>
#include <Math/GenVector/LorentzVector.h>
#include <Math/GenVector/PtEtaPhiM4D.h>

using LorentzVector = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>>;

class TaggerOutToOrbitFlatTable : public edm::global::EDProducer<> {
public:
  // constructor and destructor
  explicit TaggerOutToOrbitFlatTable(const edm::ParameterSet&);
  ~TaggerOutToOrbitFlatTable() override {};

  void produce(edm::StreamID, edm::Event&, edm::EventSetup const&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  // the tokens to access the data
  edm::EDGetTokenT<l1sc::CandsClusterBxHostCollection> srcClusters_;
  edm::EDGetTokenT<l1sc::PFCandidateHostCollection> srcCandidates_;
  edm::EDGetTokenT<l1sc::SoftTauOutputHostTensor> srcOut_;

  std::string name_, doc_;
};
// -----------------------------------------------------------------------------

// -------------------------------- constructor  -------------------------------

TaggerOutToOrbitFlatTable::TaggerOutToOrbitFlatTable(const edm::ParameterSet& iConfig) :
      srcClusters_(consumes<l1sc::CandsClusterBxHostCollection>(iConfig.getParameter<edm::InputTag>("srcClusters"))),
      srcCandidates_(consumes<l1sc::PFCandidateHostCollection>(iConfig.getParameter<edm::InputTag>("srcCandidates"))),
      srcOut_(consumes<l1sc::SoftTauOutputHostTensor>(iConfig.getParameter<edm::InputTag>("srcOut"))),
      name_(iConfig.getParameter<std::string>("name")),
      doc_(iConfig.getParameter<std::string>("doc")) {
  produces<l1ScoutingRun3::OrbitFlatTable>();
}
// -----------------------------------------------------------------------------

// ----------------------- method called for each orbit  -----------------------
void TaggerOutToOrbitFlatTable::produce(edm::StreamID, edm::Event& iEvent, edm::EventSetup const&) const {
  edm::Handle<l1sc::CandsClusterBxHostCollection> srcClusters;
  iEvent.getByToken(srcClusters_, srcClusters);
  edm::Handle<l1sc::PFCandidateHostCollection> srcCandidates;
  iEvent.getByToken(srcCandidates_, srcCandidates);
  edm::Handle<l1sc::SoftTauOutputHostTensor> srcOut;
  iEvent.getByToken(srcOut_, srcOut);

  const unsigned int nbx = srcClusters->const_view<l1sc::LongOffsetsSoA>().metadata().size() - 1;

  const auto *bx_offsets = srcClusters->const_view<l1sc::LongOffsetsSoA>().offsets().data();
  const auto *cluster_idx = srcClusters->const_view<l1sc::ClusterIndexSoA>().indexes().data();
  const auto *cluster_off = srcClusters->const_view<l1sc::ClusterOffsetsSoA>().offsets().data();
  const auto *candidate_idx = srcClusters->const_view<l1sc::LongIndexSoA>().indexes().data();
  const auto num_clusters = srcClusters->const_view<l1sc::ClusterIndexSoA>().metadata().size();

  const auto *pt = srcCandidates->const_view().pt().data();
  const auto *eta = srcCandidates->const_view().eta().data();
  const auto *phi = srcCandidates->const_view().phi().data();

  const auto *cls = srcOut->const_view().cls().data();
  const auto *pt_reg = srcOut->const_view().vz().data();
  const auto *vz_reg = srcOut->const_view().vz().data();
  const auto *charge_reg = srcOut->const_view().charge().data();
  const auto num_outputs = srcOut->const_view().metadata().size(); 

  // num_outputs is expected to be the total number of clusters
  std::cout << "num_clusters = " << num_clusters << std::endl;
  std::cout << "num_outputs = " << num_outputs << std::endl;
  assert(num_outputs == num_clusters);

  // features that are going to be serialized for each 
  // cluster that has been digested by the tagger
  std::vector<float> pts, etas, phis;
  std::vector<int32_t> clusters;
  std::vector<float> clss{cls, cls + num_outputs};
  std::vector<float> pts_reg{pt_reg, pt_reg + num_outputs};
  std::vector<float> vzs_reg{vz_reg, vz_reg + num_outputs};
  std::vector<float> charges_reg{charge_reg, charge_reg + num_outputs};

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
    for (auto cl_idx = cluster_idx[bx_start]; cl_idx < cluster_idx[bx_start] + bx_size; ++cl_idx) {
      auto cluster_start = cluster_off[cl_idx];
      auto cluster_end = cluster_off[cl_idx + 1];
      auto cluster_size = cluster_end - cluster_start;

      LorentzVector sum(0., 0., 0., 0.);

      // loop over the candidates indexes of the candidats in the current cluster
      int32_t cand_counter = 0;
      for (auto cand_idx = candidate_idx[cluster_start]; cand_idx < candidate_idx[cluster_start] + cluster_size && cand_counter < 16; ++cand_idx) {
        LorentzVector cand_lv(pt[cand_idx], eta[cand_idx], phi[cand_idx], 0.13957f);
        sum += cand_lv;
        ++cand_counter;
      }

      pts.push_back(sum.Pt());
      etas.push_back(sum.Eta());
      phis.push_back(sum.Phi());
      clusters.push_back(cl_idx);
    }

    bxOffsets.push_back(pts.size());
  }
  
  auto out = std::make_unique<l1ScoutingRun3::OrbitFlatTable>(bxOffsets, name_);
  out->setDoc(doc_);
  out->addColumn<float>("pt", pts, "cluster pt (GeV)");
  out->addColumn<float>("eta", etas, "cluster eta (GeV)");
  out->addColumn<float>("phi", phis, "cluster phi (GeV)");
  out->addColumn<float>("cls", clss, "tagger classification");
  out->addColumn<float>("pt_reg", pts_reg, "tagger pt regression");
  out->addColumn<float>("vz_reg", vzs_reg, "tagger vz regression");
  out->addColumn<float>("charge_reg", charges_reg, "tagger charge regression");
  out->addColumn<int32_t>("cluster", clusters, "cluster index");
  iEvent.put(std::move(out));
}

void TaggerOutToOrbitFlatTable::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("srcClusters");
  desc.add<edm::InputTag>("srcCandidates");
  desc.add<edm::InputTag>("srcOut");
  desc.add<std::string>("name");
  desc.add<std::string>("doc");
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(TaggerOutToOrbitFlatTable);
