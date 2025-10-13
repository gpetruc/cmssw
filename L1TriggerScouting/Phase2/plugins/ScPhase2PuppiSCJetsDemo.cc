#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/L1Scouting/interface/OrbitCollection.h"
#include "DataFormats/NanoAOD/interface/OrbitFlatTable.h"
#include "DataFormats/L1TParticleFlow/interface/L1ScoutingPuppi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "L1TriggerScouting/Utilities/interface/BxOffsetsFiller.h"

#include <ROOT/RVec.hxx>
#include <Math/Vector4D.h>
#include <Math/GenVector/LorentzVector.h>
#include <Math/GenVector/PtEtaPhiM4D.h>
#include <algorithm>
#include <array>
#include <iostream>

class ScPhase2PuppiSCJetsDemo : public edm::stream::EDProducer<> {
public:
  explicit ScPhase2PuppiSCJetsDemo(const edm::ParameterSet &);
  ~ScPhase2PuppiSCJetsDemo() override;
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  void beginStream(edm::StreamID) override;
  void produce(edm::Event &, const edm::EventSetup &) override;
  void endStream() override;

  edm::EDGetTokenT<OrbitCollection<l1Scouting::Puppi>> src_;
  float R2_, minSeedPt_;
  unsigned int nJets_;

  struct Particle {
    float pt, eta, phi;
    unsigned int i;
  };
};

ScPhase2PuppiSCJetsDemo::ScPhase2PuppiSCJetsDemo(const edm::ParameterSet &iConfig)
    : src_(consumes<OrbitCollection<l1Scouting::Puppi>>(iConfig.getParameter<edm::InputTag>("src"))),
      R2_(std::pow(iConfig.getParameter<double>("rParam"), 2)),
      minSeedPt_(iConfig.getParameter<double>("minSeedPt")),
      nJets_(iConfig.getParameter<unsigned int>("nJets")) {
  produces<l1ScoutingRun3::OrbitFlatTable>("jets");
  produces<l1ScoutingRun3::OrbitFlatTable>("clusters");
}

ScPhase2PuppiSCJetsDemo::~ScPhase2PuppiSCJetsDemo() {};

void ScPhase2PuppiSCJetsDemo::beginStream(edm::StreamID) {}

void ScPhase2PuppiSCJetsDemo::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  l1ScoutingRun3::BxOffsetsFillter bxOffsetsFiller;
  bxOffsetsFiller.start();

  edm::Handle<OrbitCollection<l1Scouting::Puppi>> src;
  iEvent.getByToken(src_, src);

  // containers for output products
  // jet 3 momenta
  std::vector<float> pt;
  std::vector<float> eta;
  std::vector<float> phi;

  // clustering information
  std::vector<int> cluster(src->size(), -1);
  std::vector<bool> is_seed(src->size(), false);

  std::vector<Particle> work, work2;
  for (unsigned int bx = 1; bx <= OrbitCollection<l1Scouting::Puppi>::NBX; ++bx) {
    auto span = src->bxIterator(bx);
    unsigned int i0 = std::distance(&*src->begin(), &*span.begin());  // global index

    work.clear();
#ifdef SC_JETS_DEBUG    
    if (bx <= 2)
      printf("Clustering BX %u with %lu particles\n", bx, span.size());
#endif
    for (unsigned int i = 0, n = span.size(); i < n; ++i) {
      auto p = span[i];
      work.emplace_back(p.pt(), p.eta(), p.phi(), i);
#ifdef SC_JETS_DEBUG    
      if (bx <= 2)
        printf("  %4u: pt %7.2f eta %+6.3f phi %+6.3f\n", i, p.pt(), p.eta(), p.phi());
#endif
    }

    unsigned int njets = 0;
    std::sort(work.begin(), work.end(), [](const Particle &i, const Particle &j) { return (i.pt > j.pt); });
    // Get the particles within a coneSize of the seed
    while (!work.empty() && njets < nJets_) {
      // Take the first (highest pt) candidate as a seed
      Particle seed = work.front();
      if (seed.pt < minSeedPt_)
        break;
      cluster[i0 + seed.i] = njets;
      is_seed[i0 + seed.i] = true;
      float sumpt = seed.pt, sumeta = 0, sumphi = 0;
#ifdef SC_JETS_DEBUG    
      if (bx <= 2)
        printf("selected %u (pt %7.2f eta %+6.3f phi %+6.3f) as seed for iteration %u (%lu/%lu particles left)\n",
               seed.i,
               seed.pt,
               seed.eta,
               seed.phi,
               njets,
               work.size(),
               span.size());
#endif
      // loop on the others
      work2.clear();
      for (auto it = work.begin() + 1, ed = work.end(); it < ed; ++it) {
        float deta = it->eta - seed.eta;
        float dphi = reco::deltaPhi(it->phi, seed.phi);
        if (deta * deta + dphi * dphi < R2_) {
          sumpt += it->pt;
          sumeta += it->pt * deta;
          sumphi += it->pt * dphi;
          cluster[i0 + it->i] = njets;
#ifdef SC_JETS_DEBUG    
          if (bx <= 2)
            printf("  %4u: pt %7.2f eta %+6.3f phi %+6.3f <<= selected (dr %7.4f)\n",
                   it->i,
                   it->pt,
                   it->eta,
                   it->phi,
                   std::hypot(deta, dphi));
#endif
        } else {
          work2.push_back(*it);
        }
      }
      njets++;
      pt.push_back(sumpt);
      eta.push_back(seed.eta + sumeta / sumpt);
      phi.push_back(reco::reducePhiRange(seed.phi + sumphi / sumpt));
      std::swap(work2, work);
#ifdef SC_JETS_DEBUG    
      if (bx <= 2)
        printf("Jet pt %7.2f eta %+6.3f phi %+6.3f, seed %u\n\n", pt.back(), eta.back(), phi.back(), seed.i);
#endif
    }

    // add the number of jets to the BX offset
    bxOffsetsFiller.addBx(bx, njets);
  }  // loop over BX

  auto bxOffsets = bxOffsetsFiller.done();
  auto tab = std::make_unique<l1ScoutingRun3::OrbitFlatTable>(bxOffsets, "SC4Jets");
  tab->addColumn<float>("pt", pt, "Jet pt");
  tab->addColumn<float>("eta", eta, "Jet eta");
  tab->addColumn<float>("phi", phi, "Jet phi");
  iEvent.put(std::move(tab), "jets");

  auto cltab = std::make_unique<l1ScoutingRun3::OrbitFlatTable>(src->bxOffsets(), "SC4Clusters");
  cltab->addColumn<int>("cluster", cluster, "cluster index (-1 if unclustered)");
  cltab->addColumn<bool>("is_seed", is_seed, "whether the particle was a seed");
  iEvent.put(std::move(cltab), "clusters");
}

void ScPhase2PuppiSCJetsDemo::endStream() {}

void ScPhase2PuppiSCJetsDemo::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src");
  desc.add<double>("rParam", 0.4);
  desc.add<double>("minSeedPt", 0.0);
  desc.add<unsigned int>("nJets", 16);
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScPhase2PuppiSCJetsDemo);