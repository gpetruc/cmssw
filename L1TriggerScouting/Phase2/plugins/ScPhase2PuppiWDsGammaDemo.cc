#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/L1Scouting/interface/OrbitCollection.h"
#include "DataFormats/L1Scouting/interface/OrbitFlatTable.h"
#include "DataFormats/L1TParticleFlow/interface/L1ScoutingPuppi.h"
#include "DataFormats/L1TParticleFlow/interface/L1ScoutingTkEm.h"
#include "L1TriggerScouting/Utilities/interface/BxOffsetsFiller.h"

#include <ROOT/RVec.hxx>
#include <Math/Vector4D.h>
#include <Math/GenVector/LorentzVector.h>
#include <Math/GenVector/PtEtaPhiM4D.h>
#include <algorithm>
#include <array>
#include <iostream>

class ScPhase2PuppiWDsGammaDemo : public edm::stream::EDProducer<> {
public:
  explicit ScPhase2PuppiWDsGammaDemo(const edm::ParameterSet &);
  ~ScPhase2PuppiWDsGammaDemo() override;
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  void beginStream(edm::StreamID) override;
  void produce(edm::Event &, const edm::EventSetup &) override;
  void endStream() override;
  template <typename T, typename U>
  void runObj(const OrbitCollection<T> &src,
	      const OrbitCollection<U> &src2,
              edm::Event &out,
              unsigned long &nTry,
              unsigned long &nPass,
              const std::string &bxLabel);

  bool doStruct_;
  edm::EDGetTokenT<OrbitCollection<l1Scouting::Puppi>> structToken_;
  edm::EDGetTokenT<OrbitCollection<l1Scouting::TkEm>> struct2Token_;

  struct Cuts {
    float minpt1 = 3;  
    float minpt2 = 5;  
    float minpt3 = 10; 
    float minpt4 = 25; 
    float maxdeltar2 = 0.15*0.15;
    float mindeltarDsGamma2 = 3.5*3.5;
    float mindeltaphiDsGamma = 2.5;
    float minmass = 60;   
    float maxmass = 100;  
    float minmass3 = 1.75;
    float maxmass3 = 2.30;
    float mindr2 = 0.00 * 0.00;
    float maxdr2 = 0.50 * 0.50;
    float maxiso = 0.45;  
    float mindr2tkem = 0.02 * 0.02;
    float maxdr2tkem = 0.50 * 0.50;
    float maxisotkem = 0.25;  
  } cuts;

  template <typename T>
  bool isolationDs(unsigned int pidex1, unsigned int pidex2, unsigned int pidex3, const T *cands, unsigned int size) const;

  template <typename T>
  bool isolationDs(unsigned int pidex1, unsigned int pidex2, unsigned int pidex3, const T *cands, unsigned int size, unsigned int &cache) const {
    if (cache == 0)
      cache = isolationDs(pidex1, pidex2, pidex3, cands, size) ? 1 : 2;
    return (cache == 1);
  }

  template <typename T>
  bool isolationTkEm(float pt, float eta, float phi, const T *cands, unsigned int size) const;

  template <typename T>
  bool isolationTkEm(float pt, float eta, float phi, const T *cands, unsigned int size, unsigned int &cache) const {
    if (cache == 0)
      cache = isolationTkEm(pt, eta, phi, cands, size) ? 1 : 2;
    return (cache == 1);
  }


  bool deltar(float eta1, float eta2, float phi1, float phi2) const;
  bool deltarmin(float eta1, float eta2, float phi1, float phi2) const;
  bool deltaphi(float phi1, float phi2) const;
  static float tripletmass(const std::array<unsigned int, 3> &t, const float *pts, const float *etas, const float *phis);
  static float quadrupletmass(const std::array<unsigned int, 4> &t, const float *pts, const float *etas, const float *phis);

  unsigned long countStruct_;
  unsigned long passStruct_;
};

ScPhase2PuppiWDsGammaDemo::ScPhase2PuppiWDsGammaDemo(const edm::ParameterSet &iConfig)
    : doStruct_(iConfig.getParameter<bool>("runStruct")) {
  if (doStruct_) {
    structToken_ = consumes<OrbitCollection<l1Scouting::Puppi>>(iConfig.getParameter<edm::InputTag>("src"));
    struct2Token_ = consumes<OrbitCollection<l1Scouting::TkEm>>(iConfig.getParameter<edm::InputTag>("src2"));
    produces<std::vector<unsigned>>("selectedBx");
    produces<l1ScoutingRun3::OrbitFlatTable>("wdsgamma");
  }
}

ScPhase2PuppiWDsGammaDemo::~ScPhase2PuppiWDsGammaDemo(){};

void ScPhase2PuppiWDsGammaDemo::beginStream(edm::StreamID) {
  countStruct_ = 0;
  passStruct_ = 0;
}

void ScPhase2PuppiWDsGammaDemo::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  if (doStruct_) {
    edm::Handle<OrbitCollection<l1Scouting::Puppi>> src;
    iEvent.getByToken(structToken_, src);

    edm::Handle<OrbitCollection<l1Scouting::TkEm>> src2;
    iEvent.getByToken(struct2Token_, src2);

    runObj(*src, *src2, iEvent, countStruct_, passStruct_, "");
  }
}

void ScPhase2PuppiWDsGammaDemo::endStream() {
  if (doStruct_)
    std::cout << "Struct analysis: " << countStruct_ << " -> " << passStruct_ << std::endl;
}

template <typename T, typename U>
void ScPhase2PuppiWDsGammaDemo::runObj(const OrbitCollection<T> &src,
				   const OrbitCollection<U> &src2,
                                   edm::Event &iEvent,
                                   unsigned long &nTry,
                                   unsigned long &nPass,
                                   const std::string &label) {
  l1ScoutingRun3::BxOffsetsFillter bxOffsetsFiller;
  bxOffsetsFiller.start();
  auto ret = std::make_unique<std::vector<unsigned>>();
  std::vector<float> masses;
  std::vector<uint8_t> i0s, i1s, i2s, i3s; //i3s is the photon
  ROOT::RVec<unsigned int> ix;   // pions, kaons
  ROOT::RVec<unsigned int> ig;   // photons
  ROOT::RVec<unsigned int> iso;  //stores whether the Ds or photon passes isolation test so we don't calculate reliso twice
  std::array<unsigned int, 4> bestQuadruplet; 
  float bestQuadrupletScore, bestQuadrupletMass;
  for (unsigned int bx = 1; bx <= OrbitCollection<T>::NBX; ++bx) {
    nTry++;
    auto range = src.bxIterator(bx);
    const T *cands = &range.front();
    auto size = range.size();

    auto range2 = src2.bxIterator(bx);
    const U *cands2 = &range2.front();
    auto size2 = range2.size();

    ix.clear();
    int intermediatecut = 0;
    int highcut = 0;
    for (unsigned int i = 0; i < size; ++i) {  //make list of all hadrons
      if ((std::abs(cands[i].pdgId()) == 211 or std::abs(cands[i].pdgId()) == 11)) {
        if (cands[i].pt() >= cuts.minpt1) {
          ix.push_back(i);
          if (cands[i].pt() >= cuts.minpt2)
            intermediatecut++;
          if (cands[i].pt() >= cuts.minpt3)
            highcut++;
        }
      }
    }
    unsigned int npions = ix.size();
    if (highcut < 1 || intermediatecut < 2 || npions < 3)
      continue;

    ig.clear();
    for (unsigned int i = 0; i < size2; ++i) {  //make list of all hadrons
        if (cands2[i].pt() >= cuts.minpt4) {
          ig.push_back(i);
        }
    }
    unsigned int ngammas = ig.size();
    if (ngammas < 1)
      continue;

    iso.resize(2); // gamma and Ds isolations
    std::fill(iso.begin(), iso.end(), 0);
    bestQuadrupletScore = 0;

    for (unsigned int i1 = 0; i1 < npions; ++i1) {
      if (cands[ix[i1]].pt() < cuts.minpt3)
        continue;  //high pt cut
      //if (isolation(ix[i1], cands, size, iso[i1]) == 0)
      //  continue;  //check iso of high pt pion
      for (unsigned int i2 = 0; i2 < npions; ++i2) {
        if (i2 == i1 || cands[ix[i2]].pt() < cuts.minpt2)
          continue;
        if (cands[ix[i2]].pt() > cands[ix[i1]].pt() || (cands[ix[i2]].pt() == cands[ix[i1]].pt() and i2 < i1))
          continue;  //intermediate pt cut
        if (!deltar(cands[ix[i1]].eta(), cands[ix[i2]].eta(), cands[ix[i1]].phi(), cands[ix[i2]].phi()))
          continue;  //angular sep of top 2 tracks
        for (unsigned int i3 = 0; i3 < npions; ++i3) {
          if (i3 == i1 or i3 == i2)
            continue;
          if (cands[ix[i2]].pt() < cuts.minpt1)
            continue;  //low pt cut
          if (cands[ix[i3]].pt() > cands[ix[i1]].pt() || (cands[ix[i3]].pt() == cands[ix[i1]].pt() and i3 < i1))
            continue;
          if (cands[ix[i3]].pt() > cands[ix[i2]].pt() || (cands[ix[i3]].pt() == cands[ix[i2]].pt() and i3 < i2))
            continue;

	  if (!isolationDs(ix[i1], ix[i2], ix[i3], cands, size, iso[0]))
	    continue;

	  auto mass3 = (cands[ix[i1]].p4() + cands[ix[i2]].p4() + cands[ix[i3]].p4()).mass(); //FIXME switch to tripletmass function
          if (mass3 >= cuts.minmass3 and mass3 <= cuts.maxmass3)
            continue;

          for (unsigned int i4 = 0; i4 < ngammas; ++i4) {
	    if (cands2[ig[i4]].pt() < cuts.minpt4)
            continue;  //photon pt cut

            std::array<unsigned int, 4> tr{{ix[i1], ix[i2], ix[i3], ig[i4]}};  //quadruplet of indices

            if (std::abs(cands[ix[i1]].charge() + cands[ix[i2]].charge() + cands[ix[i3]].charge()) == 1) {
              //make Lorentz vectors for each quadruplet
              auto mass = (cands[ix[i1]].p4() + cands[ix[i2]].p4() + cands[ix[i3]].p4() + cands2[ig[i4]].p4()).mass();
              if (mass >= cuts.minmass and mass <= cuts.maxmass) {  //MASS test
                if (deltar(cands[ix[i1]].eta(), cands[ix[i3]].eta(), cands[ix[i1]].phi(), cands[ix[i3]].phi()) and
                    deltar(cands[ix[i2]].eta(), cands[ix[i3]].eta(), cands[ix[i2]].phi(), cands[ix[i3]].phi())) {
                  //ISOLATION test for photon
                  bool isop = isolationTkEm(cands2[ig[i4]].pt(), cands2[ig[i4]].eta(), cands2[ig[i4]].phi(), cands, size, iso[1]);
		  bool pass_deltaphi = deltaphi(cands[ix[i1]].phi(), cands2[ig[i4]].phi());
		  bool pass_deltar = deltarmin(cands[ix[i1]].eta(), cands2[ig[i4]].eta(), cands[ix[i1]].phi(), cands2[ig[i4]].phi());
                  if (isop == true and pass_deltaphi==true and pass_deltar==true) {
                    float ptsum = cands[ix[i1]].pt() + cands[ix[i2]].pt() + cands[ix[i3]].pt() + cands2[ig[i4]].pt();
                    if (ptsum > bestQuadrupletScore) {
                      std::copy_n(tr.begin(), 4, bestQuadruplet.begin());
                      bestQuadrupletScore = ptsum;
                      bestQuadrupletMass = mass;
                    }
                  }  // iso
                }    // delta R
              }      // mass
            }        //charge
	  }          // photon pt cut
        }          //low pt cut
      }            //intermediate pt cut
    }              //high pt cut

    if (bestQuadrupletScore > 0) {
      ret->emplace_back(bx);
      nPass++;
      masses.push_back(bestQuadrupletMass);
      i0s.push_back(bestQuadruplet[0]);
      i1s.push_back(bestQuadruplet[1]);
      i2s.push_back(bestQuadruplet[2]);
      i3s.push_back(bestQuadruplet[3]);
      bxOffsetsFiller.addBx(bx, 1);
    }
  }  // loop on BXs

  iEvent.put(std::move(ret), "selectedBx" + label);
  // now we make the table
  auto bxOffsets = bxOffsetsFiller.done();
  auto tab = std::make_unique<l1ScoutingRun3::OrbitFlatTable>(bxOffsets, "WDsGamma" + label, true);
  tab->addColumn<float>("mass", masses, "3-pion invariant mass");
  tab->addColumn<uint8_t>("i0", i0s, "leading pion");
  tab->addColumn<uint8_t>("i1", i1s, "subleading pion");
  tab->addColumn<uint8_t>("i2", i2s, "trailing pion");
  tab->addColumn<uint8_t>("i3", i3s, "photon");
  iEvent.put(std::move(tab), "wdsgamma" + label);
}

//TEST functions
template <typename T>
bool ScPhase2PuppiWDsGammaDemo::isolationDs(unsigned int pidex1, unsigned int pidex2, unsigned int pidex3, const T *cands, unsigned int size) const {
  bool passed = false;
  float psum = 0;
  float eta = cands[pidex1].eta(); //center cone around leading track
  float phi = cands[pidex1].phi();
  for (unsigned int j = 0u; j < size; ++j) {  //loop over other particles
    if (pidex1 == j or pidex2==j or pidex3==j)
      continue;
    float deta = eta - cands[j].eta(), dphi = ROOT::VecOps::DeltaPhi<float>(phi, cands[j].phi());
    float dr2 = deta * deta + dphi * dphi;
    if (dr2 >= cuts.mindr2 && dr2 <= cuts.maxdr2)
      psum += cands[j].pt();
  }
  if (psum <= cuts.maxiso * (cands[pidex1].pt() + cands[pidex2].pt() + cands[pidex3].pt()))
    passed = true;
  return passed;
}

template <typename T>
bool ScPhase2PuppiWDsGammaDemo::isolationTkEm(float pt, float eta, float phi, const T *cands, unsigned int size) const {
  bool passed = false;
  float psum = 0;
  for (unsigned int j = 0u; j < size; ++j) {  //loop over other particles
    float deta = eta - cands[j].eta(), dphi = ROOT::VecOps::DeltaPhi<float>(phi, cands[j].phi());
    float dr2 = deta * deta + dphi * dphi;
    if (dr2 >= cuts.mindr2tkem && dr2 <= cuts.maxdr2tkem)
      psum += cands[j].pt();
  }
  if (psum <= cuts.maxisotkem * pt)
    passed = true;
  return passed;
}

bool ScPhase2PuppiWDsGammaDemo::deltar(float eta1, float eta2, float phi1, float phi2) const {
  bool passed = true;
  float deta = eta1 - eta2;
  float dphi = ROOT::VecOps::DeltaPhi<float>(phi1, phi2);
  float dr2 = deta * deta + dphi * dphi;
  if (dr2 > cuts.maxdeltar2) {
    passed = false;
    return passed;
  }
  return passed;
}

bool ScPhase2PuppiWDsGammaDemo::deltarmin(float eta1, float eta2, float phi1, float phi2) const {
  bool passed = true;
  float deta = eta1 - eta2;
  float dphi = ROOT::VecOps::DeltaPhi<float>(phi1, phi2);
  float dr2 = deta * deta + dphi * dphi;
  if (dr2 < cuts.mindeltarDsGamma2) {
    passed = false;
    return passed;
  }
  return passed;
}

bool ScPhase2PuppiWDsGammaDemo::deltaphi(float phi1, float phi2) const {
  bool passed = true;
  float dphi = ROOT::VecOps::DeltaPhi<float>(phi1, phi2);
  if (fabs(dphi) < cuts.mindeltaphiDsGamma) {
    passed = false;
    return passed;
  }
  return passed;
}

float ScPhase2PuppiWDsGammaDemo::tripletmass(const std::array<unsigned int, 3> &t,
                                         const float *pts,
                                         const float *etas,
                                         const float *phis) {
  ROOT::Math::PtEtaPhiMVector p1(pts[t[0]], etas[t[0]], phis[t[0]], 0.1396);
  ROOT::Math::PtEtaPhiMVector p2(pts[t[1]], etas[t[1]], phis[t[1]], 0.1396);
  ROOT::Math::PtEtaPhiMVector p3(pts[t[2]], etas[t[2]], phis[t[2]], 0.1396);
  float mass = (p1 + p2 + p3).M();
  return mass;
}

float ScPhase2PuppiWDsGammaDemo::quadrupletmass(const std::array<unsigned int, 4> &t,
                                         const float *pts,
                                         const float *etas,
                                         const float *phis) {
  ROOT::Math::PtEtaPhiMVector p1(pts[t[0]], etas[t[0]], phis[t[0]], 0.1396);
  ROOT::Math::PtEtaPhiMVector p2(pts[t[1]], etas[t[1]], phis[t[1]], 0.1396);
  ROOT::Math::PtEtaPhiMVector p3(pts[t[2]], etas[t[2]], phis[t[2]], 0.1396);
  ROOT::Math::PtEtaPhiMVector p4(pts[t[3]], etas[t[3]], phis[t[3]], 0.0000);
  float mass = (p1 + p2 + p3 + p4).M();
  return mass;
}

void ScPhase2PuppiWDsGammaDemo::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScPhase2PuppiWDsGammaDemo);
