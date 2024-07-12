#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/L1Scouting/interface/OrbitCollection.h"
#include "DataFormats/L1TParticleFlow/interface/PFCandidate.h"
#include "DataFormats/L1TParticleFlow/interface/L1ScoutingPuppi.h"
#include "L1TriggerScouting/Phase2/interface/phase2Utils.h"
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
  template <typename T>
  void runObj(const OrbitCollection<T> &src,
              edm::Event &out,
              unsigned long &nTry,
              unsigned long &nPass,
              const std::string &bxLabel);
  void runSOA(const l1Scouting::PuppiSOA &src, edm::Event &out);

  bool doCandidate_, doStruct_, doSOA_;
  edm::EDGetTokenT<OrbitCollection<l1t::PFCandidate>> candidateToken_;
  edm::EDGetTokenT<OrbitCollection<l1Scouting::Puppi>> structToken_;
  edm::EDGetTokenT<l1Scouting::PuppiSOA> soaToken_;

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

  template <typename T>
  bool isolation(unsigned int pidex, const T *cands, unsigned int size) const;

  template <typename T>
  bool isolation(unsigned int pidex, const T *cands, unsigned int size, unsigned int &cache) const {
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
  static float tripletmass(const std::array<unsigned int, 3> &t, const float *pts, const float *etas, const float *phis);

  unsigned long countCandidate_, countStruct_, countSOA_;
  unsigned long passCandidate_, passStruct_, passSOA_;
};

ScPhase2PuppiW3PiDemo::ScPhase2PuppiW3PiDemo(const edm::ParameterSet &iConfig)
    : doCandidate_(iConfig.getParameter<bool>("runCandidate")),
      doStruct_(iConfig.getParameter<bool>("runStruct")),
      doSOA_(iConfig.getParameter<bool>("runSOA")) {
  if (doCandidate_) {
    candidateToken_ = consumes<OrbitCollection<l1t::PFCandidate>>(iConfig.getParameter<edm::InputTag>("src"));
    produces<std::vector<unsigned>>("selectedBxCandidate");
  }
  if (doStruct_) {
    structToken_ = consumes<OrbitCollection<l1Scouting::Puppi>>(iConfig.getParameter<edm::InputTag>("src"));
    produces<std::vector<unsigned>>("selectedBx");
  }
  if (doSOA_) {
    soaToken_ = consumes<l1Scouting::PuppiSOA>(iConfig.getParameter<edm::InputTag>("src"));
    produces<l1Scouting::PuppiSOA>();
  }
}

ScPhase2PuppiW3PiDemo::~ScPhase2PuppiW3PiDemo(){};

void ScPhase2PuppiW3PiDemo::beginStream(edm::StreamID) {
  countCandidate_ = 0;
  countStruct_ = 0;
  countSOA_ = 0;
  passCandidate_ = 0;
  passStruct_ = 0;
  passSOA_ = 0;
}

void ScPhase2PuppiW3PiDemo::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  if (doCandidate_) {
    edm::Handle<OrbitCollection<l1t::PFCandidate>> src;
    iEvent.getByToken(candidateToken_, src);
    runObj(*src, iEvent, countCandidate_, passCandidate_, "Candidate");
  }
  if (doStruct_) {
    edm::Handle<OrbitCollection<l1Scouting::Puppi>> src;
    iEvent.getByToken(structToken_, src);
    runObj(*src, iEvent, countStruct_, passStruct_, "");
  }
  if (doSOA_) {
    edm::Handle<l1Scouting::PuppiSOA> src;
    iEvent.getByToken(soaToken_, src);
    runSOA(*src, iEvent);
  }
}

void ScPhase2PuppiW3PiDemo::endStream() {
  if (doCandidate_)
    std::cout << "Candidate analysis: " << countCandidate_ << " -> " << passCandidate_ << std::endl;
  if (doStruct_)
    std::cout << "Struct analysis: " << countStruct_ << " -> " << passStruct_ << std::endl;
  if (doSOA_)
    std::cout << "SOA analysis: " << countSOA_ << " -> " << passSOA_ << std::endl;
}

template <typename T>
void ScPhase2PuppiW3PiDemo::runObj(const OrbitCollection<T> &src,
                                   edm::Event &iEvent,
                                   unsigned long &nTry,
                                   unsigned long &nPass,
                                   const std::string &label) {
  auto ret = std::make_unique<std::vector<unsigned>>();
  ROOT::RVec<unsigned int> ix;   // pions
  ROOT::RVec<unsigned int> iso;  //stores whether a particle passes isolation test so we don't calculate reliso twice
  std::array<unsigned int, 3> bestTriplet;  // best triplet
  float bestTripletScore;
  for (unsigned int bx = 1; bx <= OrbitCollection<T>::NBX; ++bx) {
    nTry++;
    auto range = src.bxIterator(bx);
    const T *cands = &range.front();
    auto size = range.size();
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
      ret->emplace_back(bx);
      nPass++;
    }
  }  // loop on BXs

  iEvent.put(std::move(ret), "selectedBx" + label);
}

void ScPhase2PuppiW3PiDemo::runSOA(const l1Scouting::PuppiSOA &src, edm::Event &iEvent) {
  auto ret = std::make_unique<l1Scouting::PuppiSOA>();
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
      passSOA_++;
    }
  }  // loop on BXs

  iEvent.put(std::move(ret));
}

//TEST functions
template <typename T>
bool ScPhase2PuppiW3PiDemo::isolation(unsigned int pidex, const T *cands, unsigned int size) const {
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
