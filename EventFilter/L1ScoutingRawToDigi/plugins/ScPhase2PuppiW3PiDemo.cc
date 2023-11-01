#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/FEDRawData/interface/FEDRawData.h"
#include "DataFormats/L1Scouting/interface/SRawDataCollection.h"
#include "DataFormats/L1Scouting/interface/SDSNumbering.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "DataFormats/L1Trigger/interface/SparseBXVector.h"
#include "DataFormats/L1TParticleFlow/interface/PFCandidate.h"
#include "DataFormats/L1TParticleFlow/interface/puppi.h"
#include "DataFormats/L1TParticleFlow/interface/struct_types.h"
#include "phase2Utils.h"
#include <ROOT/RVec.hxx>
#include <Math/Vector4D.h>
#include <Math/GenVector/LorentzVector.h>
#include <Math/GenVector/PtEtaPhiM4D.h>
#include <algorithm>
#include <array>
#include <iostream>

class ScPhase2PuppiW3PiDemo : public edm::stream::EDProducer<> {
public:
  explicit ScPhase2PuppiW3PiDemo(const edm::ParameterSet &);
  ~ScPhase2PuppiW3PiDemo() override;
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  void beginStream(edm::StreamID) override;
  void produce(edm::Event &, const edm::EventSetup &) override;
  void endStream() override;
  std::unique_ptr<SparseBXVector<l1t::PFCandidate>> runSimple(const SparseBXVector<l1t::PFCandidate> &src);
  std::unique_ptr<SparseBXVector<l1ct::structs::Puppi>> runStruct(const SparseBXVector<l1ct::structs::Puppi> &src);
  std::unique_ptr<l1ct::structs::PuppiSOA> runSOA(const l1ct::structs::PuppiSOA &src);

  bool doSimple_, doStruct_, doSOA_, noWrite_;
  edm::EDGetTokenT<SparseBXVector<l1t::PFCandidate>> simpleToken_;
  edm::EDGetTokenT<SparseBXVector<l1ct::structs::Puppi>> structToken_;
  edm::EDGetTokenT<l1ct::structs::PuppiSOA> soaToken_;

  struct Cuts {
    float minpt1 = 7;   // 9
    float minpt2 = 12;  // 15
    float minpt3 = 15;  // 20
    float mindeltar2 = 0.5 * 0.5;
    float minmass = 40;   // 60
    float maxmass = 150;  // 100
    float mindr2 = 0.01 * 0.01;
    float maxdr2 = 0.25 * 0.25;
    float maxiso = 2.0;  //0.4
  } cuts;

  bool isolation(unsigned int pidex, const l1t::PFCandidate *cands, unsigned int size) const;
  bool isolation(unsigned int pidex, const l1ct::structs::Puppi *cands, unsigned int size) const;
  bool isolation(unsigned int pidex, const l1t::PFCandidate *cands, unsigned int size, unsigned int &cache) const {
    if (cache == 0)
      cache = isolation(pidex, cands, size) ? 1 : 2;
    return (cache == 1);
  }
  bool isolation(unsigned int pidex, const l1ct::structs::Puppi *cands, unsigned int size, unsigned int &cache) const {
    if (cache == 0)
      cache = isolation(pidex, cands, size) ? 1 : 2;
    return (cache == 1);
  }

  bool isolation(unsigned int pidex,
                 unsigned int npx,
                 const float *eta,
                 const float *phi,
                 const float *pt,
                 unsigned int &cache) const {
    if (cache == 0)
      cache = isolation(pidex, npx, eta, phi, pt) ? 1 : 2;
    return (cache == 1);
  }
  bool isolation(unsigned int pidex, unsigned int npx, const float *eta, const float *phi, const float *pt) const;
  bool deltar(float eta1, float eta2, float phi1, float phi2) const;
  static float tripletmass(const std::array<unsigned int, 3> &t, const l1ct::structs::Puppi *cands);
  static float tripletmass(const std::array<unsigned int, 3> &t, const float *pts, const float *etas, const float *phis);

  unsigned long countSimple_, countStruct_, countSOA_;
  unsigned long passSimple_, passStruct_, passSOA_;
};

ScPhase2PuppiW3PiDemo::ScPhase2PuppiW3PiDemo(const edm::ParameterSet &iConfig)
    : doSimple_(iConfig.getParameter<bool>("runSimple")),
      doStruct_(iConfig.getParameter<bool>("runStruct")),
      doSOA_(iConfig.getParameter<bool>("runSOA")),
      noWrite_(iConfig.getParameter<bool>("noWrite")) {
  if (doSimple_) {
    simpleToken_ = consumes<SparseBXVector<l1t::PFCandidate>>(iConfig.getParameter<edm::InputTag>("src"));
    produces<SparseBXVector<l1t::PFCandidate>>();
  }
  if (doStruct_) {
    structToken_ = consumes<SparseBXVector<l1ct::structs::Puppi>>(iConfig.getParameter<edm::InputTag>("src"));
    produces<SparseBXVector<l1ct::structs::Puppi>>();
  }
  if (doSOA_) {
    soaToken_ = consumes<l1ct::structs::PuppiSOA>(iConfig.getParameter<edm::InputTag>("src"));
    produces<l1ct::structs::PuppiSOA>();
  }
}

ScPhase2PuppiW3PiDemo::~ScPhase2PuppiW3PiDemo(){};

void ScPhase2PuppiW3PiDemo::beginStream(edm::StreamID) {
  countSimple_ = 0;
  countStruct_ = 0;
  countSOA_ = 0;
  passSimple_ = 0;
  passStruct_ = 0;
  passSOA_ = 0;
}

void ScPhase2PuppiW3PiDemo::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  if (doSimple_) {
    edm::Handle<SparseBXVector<l1t::PFCandidate>> src;
    iEvent.getByToken(simpleToken_, src);
    iEvent.put(runSimple(*src));
  }
  if (doStruct_) {
    edm::Handle<SparseBXVector<l1ct::structs::Puppi>> src;
    iEvent.getByToken(structToken_, src);
    iEvent.put(runStruct(*src));
  }
  if (doSOA_) {
    edm::Handle<l1ct::structs::PuppiSOA> src;
    iEvent.getByToken(soaToken_, src);
    iEvent.put(runSOA(*src));
  }
}

void ScPhase2PuppiW3PiDemo::endStream() {
  if (doSimple_)
    std::cout << "Simple analysis: " << countSimple_ << " -> " << passSimple_ << std::endl;
  if (doStruct_)
    std::cout << "Struct analysis: " << countStruct_ << " -> " << passStruct_ << std::endl;
  if (doSOA_)
    std::cout << "SOA analysis: " << countSOA_ << " -> " << passSOA_ << std::endl;
}

std::unique_ptr<SparseBXVector<l1t::PFCandidate>> ScPhase2PuppiW3PiDemo::runSimple(
    const SparseBXVector<l1t::PFCandidate> &src) {
  auto ret = std::make_unique<SparseBXVector<l1t::PFCandidate>>();
  ROOT::RVec<unsigned int> ix;   // pions
  ROOT::RVec<unsigned int> iso;  //stores whether a particle passes isolation test so we don't calculate reliso twice
  std::array<unsigned int, 3> bestTriplet;  // best triplet
  float bestTripletScore;
  for (auto iter = src.iter(); iter.good(); ++iter) {
    countSimple_++;
    const l1t::PFCandidate *cands = iter.data();
    auto size = iter.size();
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
    iso.resize(npions);
    std::fill(iso.begin(), iso.end(), 0);
    bestTripletScore = 0;

    for (unsigned int i1 = 0; i1 < npions; ++i1) {
      if (cands[ix[i1]].pt() < cuts.minpt3)
        continue;  //high pt cut
      if (isolation(ix[i1], cands, size, iso[i1]) == 0)
        continue;  //check iso of high pt pion
      for (unsigned int i2 = 0; i2 < npions; ++i2) {
        if (i2 == i1 || cands[ix[i2]].pt() < cuts.minpt2)
          continue;
        if (cands[ix[i2]].pt() > cands[ix[i1]].pt() || (cands[ix[i2]].pt() == cands[ix[i1]].pt() and i2 < i1))
          continue;  //intermediate pt cut
        if (!deltar(cands[ix[i1]].eta(), cands[ix[i2]].eta(), cands[ix[i1]].phi(), cands[ix[i2]].phi()))
          continue;  //angular sep of top 2 pions
        for (unsigned int i3 = 0; i3 < npions; ++i3) {
          if (i3 == i1 or i3 == i2)
            continue;
          if (cands[ix[i2]].pt() < cuts.minpt1)
            continue;  //low pt cut
          if (cands[ix[i3]].pt() > cands[ix[i1]].pt() || (cands[ix[i3]].pt() == cands[ix[i1]].pt() and i3 < i1))
            continue;
          if (cands[ix[i3]].pt() > cands[ix[i2]].pt() || (cands[ix[i3]].pt() == cands[ix[i2]].pt() and i3 < i2))
            continue;
          std::array<unsigned int, 3> tr{{ix[i1], ix[i2], ix[i3]}};  //triplet of indeces

          if (std::abs(cands[ix[i1]].charge() + cands[ix[i2]].charge() + cands[ix[i3]].charge()) == 1) {
            //make Lorentz vectors for each triplet
            auto mass = (cands[ix[i1]].p4() + cands[ix[i2]].p4() + cands[ix[i3]].p4()).mass();
            if (mass >= cuts.minmass and mass <= cuts.maxmass) {  //MASS test
              if (deltar(cands[ix[i1]].eta(), cands[ix[i3]].eta(), cands[ix[i1]].phi(), cands[ix[i3]].phi()) and
                  deltar(cands[ix[i2]].eta(), cands[ix[i3]].eta(), cands[ix[i2]].phi(), cands[ix[i3]].phi())) {
                //ISOLATION test for lower 4 pions
                bool isop = isolation(ix[i2], cands, size, iso[i2]) && isolation(ix[i3], cands, size, iso[i3]);
                if (isop == true) {
                  float ptsum = cands[ix[i1]].pt() + cands[ix[i2]].pt() + cands[ix[i3]].pt();
                  if (ptsum > bestTripletScore) {
                    std::copy_n(tr.begin(), 3, bestTriplet.begin());
                    bestTripletScore = ptsum;
                  }
                }  // iso
              }    // delta R
            }      // mass
          }        //charge
        }          //low pt cut
      }            //intermediate pt cut
    }              //high pt cut

    if (bestTripletScore > 0) {
      if (!noWrite_)
        ret->addBX(iter.bx(), cands, cands + size);
      passSimple_++;
    }
  }  // loop on BXs

  return ret;
}

std::unique_ptr<SparseBXVector<l1ct::structs::Puppi>> ScPhase2PuppiW3PiDemo::runStruct(
    const SparseBXVector<l1ct::structs::Puppi> &src) {
  auto ret = std::make_unique<SparseBXVector<l1ct::structs::Puppi>>();
  ROOT::RVec<unsigned int> ix;   // pions
  ROOT::RVec<unsigned int> iso;  //stores whether a particle passes isolation test so we don't calculate reliso twice
  ROOT::RVec<int> charge;        //stores whether a particle passes isolation test so we don't calculate reliso twice
  std::array<unsigned int, 3> bestTriplet;  // best triplet
  float bestTripletScore;
  for (auto iter = src.iter(); iter.good(); ++iter) {
    countStruct_++;
    const l1ct::structs::Puppi *cands = &*iter.begin();
    auto size = iter.size();
    ix.clear();
    charge.clear();
    int intermediatecut = 0;
    int highcut = 0;
    for (unsigned int i = 0; i < size; ++i) {  //make list of all hadrons
      if ((std::abs(cands[i].pdgId) == 211 or std::abs(cands[i].pdgId) == 11)) {
        if (cands[i].pt >= cuts.minpt1) {
          ix.push_back(i);
          charge.push_back(abs(cands[i].pdgId) == 11 ? (cands[i].pdgId > 0 ? -1 : +1) : (cands[i].pdgId > 0 ? +1 : -1));
          if (cands[i].pt >= cuts.minpt2)
            intermediatecut++;
          if (cands[i].pt >= cuts.minpt3)
            highcut++;
        }
      }
    }
    unsigned int npions = ix.size();
    if (highcut < 1 || intermediatecut < 2 || npions < 3)
      continue;
    iso.resize(npions);
    std::fill(iso.begin(), iso.end(), 0);
    bestTripletScore = 0;

    for (unsigned int i1 = 0; i1 < npions; ++i1) {
      if (cands[ix[i1]].pt < cuts.minpt3)
        continue;  //high pt cut
      if (isolation(ix[i1], cands, size, iso[i1]) == 0)
        continue;  //check iso of high pt pion
      for (unsigned int i2 = 0; i2 < npions; ++i2) {
        if (i2 == i1 || cands[ix[i2]].pt < cuts.minpt2)
          continue;
        if (cands[ix[i2]].pt > cands[ix[i1]].pt || (cands[ix[i2]].pt == cands[ix[i1]].pt and i2 < i1))
          continue;  //intermediate pt cut
        if (!deltar(cands[ix[i1]].eta, cands[ix[i2]].eta, cands[ix[i1]].phi, cands[ix[i2]].phi))
          continue;  //angular sep of top 2 pions
        for (unsigned int i3 = 0; i3 < npions; ++i3) {
          if (i3 == i1 or i3 == i2)
            continue;
          if (cands[ix[i2]].pt < cuts.minpt1)
            continue;  //low pt cut
          if (cands[ix[i3]].pt > cands[ix[i1]].pt || (cands[ix[i3]].pt == cands[ix[i1]].pt and i3 < i1))
            continue;
          if (cands[ix[i3]].pt > cands[ix[i2]].pt || (cands[ix[i3]].pt == cands[ix[i2]].pt and i3 < i2))
            continue;
          std::array<unsigned int, 3> tr{{ix[i1], ix[i2], ix[i3]}};  //triplet of indeces

          if (std::abs(charge[i1] + charge[i2] + charge[i3]) == 1) {
            //make Lorentz vectors for each triplet
            auto mass = tripletmass(tr, cands);
            if (mass >= cuts.minmass and mass <= cuts.maxmass) {  //MASS test
              if (deltar(cands[ix[i1]].eta, cands[ix[i3]].eta, cands[ix[i1]].phi, cands[ix[i3]].phi) and
                  deltar(cands[ix[i2]].eta, cands[ix[i3]].eta, cands[ix[i2]].phi, cands[ix[i3]].phi)) {
                //ISOLATION test for lower 4 pions
                bool isop = isolation(ix[i2], cands, size, iso[i2]) && isolation(ix[i3], cands, size, iso[i3]);
                if (isop == true) {
                  float ptsum = cands[ix[i1]].pt + cands[ix[i2]].pt + cands[ix[i3]].pt;
                  if (ptsum > bestTripletScore) {
                    std::copy_n(tr.begin(), 3, bestTriplet.begin());
                    bestTripletScore = ptsum;
                  }
                }  // iso
              }    // delta R
            }      // mass
          }        //charge
        }          //low pt cut
      }            //intermediate pt cut
    }              //high pt cut

    if (bestTripletScore > 0) {
      if (!noWrite_)
        ret->addBX(iter.bx(), cands, cands + size);
      passStruct_++;
    }
  }  // loop on BXs
  return ret;
}

std::unique_ptr<l1ct::structs::PuppiSOA> ScPhase2PuppiW3PiDemo::runSOA(const l1ct::structs::PuppiSOA &src) {
  auto ret = std::make_unique<l1ct::structs::PuppiSOA>();
  std::vector<uint32_t> &offsets = ret->offsets;
  offsets.push_back(0);
  ROOT::RVec<unsigned int> ix;   // pions
  ROOT::RVec<unsigned int> iso;  //stores whether a particle passes isolation test so we don't calculate reliso twice
  ROOT::RVec<int> charge;        //stores whether a particle passes isolation test so we don't calculate reliso twice
  std::array<unsigned int, 3> bestTriplet;  // best triplet
  float bestTripletScore;
  for (unsigned int ibx = 0, nbx = src.bx.size(); ibx < nbx; ++ibx) {
    countSOA_++;
    unsigned int offs = src.offsets[ibx];
    unsigned int size = src.offsets[ibx + 1] - offs;
    const float *pts = &src.pt[offs];
    const float *etas = &src.eta[offs];
    const float *phis = &src.phi[offs];
    const int16_t *pdgIds = &src.pdgId[offs];
    ix.clear();
    charge.clear();
    int intermediatecut = 0;
    int highcut = 0;
    for (unsigned int i = 0; i < size; ++i) {  //make list of all hadrons
      if ((std::abs(pdgIds[i]) == 211 or std::abs(pdgIds[i]) == 11)) {
        if (pts[i] >= cuts.minpt1) {
          ix.push_back(i);
          charge.push_back(abs(pdgIds[i]) == 11 ? (pdgIds[i] > 0 ? -1 : +1) : (pdgIds[i] > 0 ? +1 : -1));
          if (pts[i] >= cuts.minpt2)
            intermediatecut++;
          if (pts[i] >= cuts.minpt3)
            highcut++;
        }
      }
    }
    unsigned int npions = ix.size();
    if (highcut < 1 || intermediatecut < 2 || npions < 3)
      continue;
    iso.resize(npions);
    std::fill(iso.begin(), iso.end(), 0);
    bestTripletScore = 0;

    for (unsigned int i1 = 0; i1 < npions; ++i1) {
      if (pts[ix[i1]] < cuts.minpt3)
        continue;  //high pt cut
      if (isolation(ix[i1], size, etas, phis, pts, iso[i1]) == 0)
        continue;  //check iso of high pt pion
      for (unsigned int i2 = 0; i2 < npions; ++i2) {
        if (i2 == i1 || pts[ix[i2]] < cuts.minpt2)
          continue;
        if (pts[ix[i2]] > pts[ix[i1]] || (pts[ix[i2]] == pts[ix[i1]] and i2 < i1))
          continue;  //intermediate pt cut
        if (!deltar(etas[ix[i1]], etas[ix[i2]], phis[ix[i1]], phis[ix[i2]]))
          continue;  //angular sep of top 2 pions
        for (unsigned int i3 = 0; i3 < npions; ++i3) {
          if (i3 == i1 or i3 == i2)
            continue;
          if (pts[ix[i2]] < cuts.minpt1)
            continue;  //low pt cut
          if (pts[ix[i3]] > pts[ix[i1]] || (pts[ix[i3]] == pts[ix[i1]] and i3 < i1))
            continue;
          if (pts[ix[i3]] > pts[ix[i2]] || (pts[ix[i3]] == pts[ix[i2]] and i3 < i2))
            continue;
          std::array<unsigned int, 3> tr{{ix[i1], ix[i2], ix[i3]}};  //triplet of indeces

          if (std::abs(charge[i1] + charge[i2] + charge[i3]) == 1) {
            //make Lorentz vectors for each triplet
            auto mass = tripletmass(tr, pts, etas, phis);
            if (mass >= cuts.minmass and mass <= cuts.maxmass) {  //MASS test
              if (deltar(etas[ix[i1]], etas[ix[i3]], phis[ix[i1]], phis[ix[i3]]) and
                  deltar(etas[ix[i2]], etas[ix[i3]], phis[ix[i2]], phis[ix[i3]])) {
                //ISOLATION test for lower 4 pions
                bool isop = isolation(ix[i2], size, etas, phis, pts, iso[i2]) &&
                            isolation(ix[i3], size, etas, phis, pts, iso[i3]);
                if (isop == true) {
                  float ptsum = pts[ix[i1]] + pts[ix[i2]] + pts[ix[i3]];
                  if (ptsum > bestTripletScore) {
                    std::copy_n(tr.begin(), 3, bestTriplet.begin());
                    bestTripletScore = ptsum;
                  }
                }  // iso
              }    // delta R
            }      // mass
          }        //charge
        }          //low pt cut
      }            //intermediate pt cut
    }              //high pt cut

    if (bestTripletScore > 0) {
      if (!noWrite_) {
        offsets.push_back(offsets.back() + size);
        ret->bx.push_back(src.bx[ibx]);
        ret->pt.insert(ret->pt.end(), pts, pts + size);
        ret->eta.insert(ret->eta.end(), etas, etas + size);
        ret->phi.insert(ret->phi.end(), phis, phis + size);
        ret->pdgId.insert(ret->pdgId.end(), pdgIds, pdgIds + size);
        ret->z0.insert(ret->z0.end(), &src.z0[offs], &src.z0[offs + size]);
        ret->dxy.insert(ret->dxy.end(), &src.dxy[offs], &src.dxy[offs + size]);
        ret->puppiw.insert(ret->puppiw.end(), &src.puppiw[offs], &src.puppiw[offs + size]);
        ret->quality.insert(ret->quality.end(), &src.quality[offs], &src.quality[offs + size]);
      }
      passSOA_++;
    }
  }  // loop on BXs
  return ret;
}

//TEST functions
bool ScPhase2PuppiW3PiDemo::isolation(unsigned int pidex, const l1t::PFCandidate *cands, unsigned int size) const {
  bool passed = false;
  float psum = 0;
  float eta = cands[pidex].eta();
  float phi = cands[pidex].phi();
  for (unsigned int j = 0u; j < size; ++j) {  //loop over other particles
    if (pidex == j)
      continue;
    float deta = eta - cands[j].eta(), dphi = ROOT::VecOps::DeltaPhi<float>(phi, cands[j].phi());
    float dr2 = deta * deta + dphi * dphi;
    if (dr2 >= cuts.mindr2 && dr2 <= cuts.maxdr2)
      psum += cands[j].pt();
  }
  if (psum <= cuts.maxiso * cands[pidex].pt())
    passed = true;
  return passed;
}

bool ScPhase2PuppiW3PiDemo::isolation(unsigned int pidex, const l1ct::structs::Puppi *cands, unsigned int size) const {
  bool passed = false;
  float psum = 0;
  float eta = cands[pidex].eta;
  float phi = cands[pidex].phi;
  for (unsigned int j = 0u; j < size; ++j) {  //loop over other particles
    if (pidex == j)
      continue;
    float deta = eta - cands[j].eta, dphi = ROOT::VecOps::DeltaPhi<float>(phi, cands[j].phi);
    float dr2 = deta * deta + dphi * dphi;
    if (dr2 >= cuts.mindr2 && dr2 <= cuts.maxdr2)
      psum += cands[j].pt;
  }
  if (psum <= cuts.maxiso * cands[pidex].pt)
    passed = true;
  return passed;
}

bool ScPhase2PuppiW3PiDemo::isolation(
    unsigned int pidex, unsigned int npx, const float *eta, const float *phi, const float *pt) const {
  bool passed = false;
  float psum = 0;
  for (unsigned int j = 0u, n = npx; j < n; ++j) {  //loop over other particles
    if (pidex == j)
      continue;
    float deta = eta[pidex] - eta[j], dphi = ROOT::VecOps::DeltaPhi<float>(phi[pidex], phi[j]);
    float dr2 = deta * deta + dphi * dphi;
    if (dr2 >= cuts.mindr2 && dr2 <= cuts.maxdr2)
      psum += pt[j];
  }
  if (psum <= cuts.maxiso * pt[pidex])
    passed = true;
  return passed;
}

bool ScPhase2PuppiW3PiDemo::deltar(float eta1, float eta2, float phi1, float phi2) const {
  bool passed = true;
  float deta = eta1 - eta2;
  float dphi = ROOT::VecOps::DeltaPhi<float>(phi1, phi2);
  float dr2 = deta * deta + dphi * dphi;
  if (dr2 < cuts.mindeltar2) {
    passed = false;
    return passed;
  }
  return passed;
}

float ScPhase2PuppiW3PiDemo::tripletmass(const std::array<unsigned int, 3> &t, const l1ct::structs::Puppi *cands) {
  ROOT::Math::PtEtaPhiMVector p1(cands[t[0]].pt, cands[t[0]].eta, cands[t[0]].phi, 0.1396);
  ROOT::Math::PtEtaPhiMVector p2(cands[t[1]].pt, cands[t[1]].eta, cands[t[1]].phi, 0.1396);
  ROOT::Math::PtEtaPhiMVector p3(cands[t[2]].pt, cands[t[2]].eta, cands[t[2]].phi, 0.1396);
  float mass = (p1 + p2 + p3).M();
  return mass;
}

float ScPhase2PuppiW3PiDemo::tripletmass(const std::array<unsigned int, 3> &t,
                                         const float *pts,
                                         const float *etas,
                                         const float *phis) {
  ROOT::Math::PtEtaPhiMVector p1(pts[t[0]], etas[t[0]], phis[t[0]], 0.1396);
  ROOT::Math::PtEtaPhiMVector p2(pts[t[1]], etas[t[1]], phis[t[1]], 0.1396);
  ROOT::Math::PtEtaPhiMVector p3(pts[t[2]], etas[t[2]], phis[t[2]], 0.1396);
  float mass = (p1 + p2 + p3).M();
  return mass;
}

void ScPhase2PuppiW3PiDemo::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScPhase2PuppiW3PiDemo);
