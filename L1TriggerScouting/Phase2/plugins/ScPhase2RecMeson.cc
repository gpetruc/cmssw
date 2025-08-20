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
#include "DataFormats/L1TParticleFlow/interface/RecMeson.h"
#include "L1TriggerScouting/Utilities/interface/BxOffsetsFiller.h"

#include <ROOT/RVec.hxx>
#include <Math/Vector4D.h>
#include <Math/GenVector/LorentzVector.h>
#include <Math/GenVector/PtEtaPhiM4D.h>
#include <algorithm>
#include <array>
#include <iostream>

//CHANGES TO IMPLEMENT
//- RETURN THE FULL 4 particles ?

class ScPhase2RecMeson : public edm::stream::EDProducer<> {
public:
  explicit ScPhase2RecMeson(const edm::ParameterSet &);
  ~ScPhase2RecMeson() override;
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  void beginStream(edm::StreamID) override;
  void produce(edm::Event &, const edm::EventSetup &) override;
  void endStream() override;
  template <typename T>
  void runObj(const OrbitCollection<T> &src,
              edm::Event &out,
              const std::string &bxLabel);

  bool doStruct_;
  edm::EDGetTokenT<OrbitCollection<l1Scouting::Puppi>> structToken_;
  std::string mesonType_;
  std::vector<float> mesonMassRange_ = {0.0f, 0.0f};
  float dmass1_ = 0;
  float dmass2_ = 0;

  double minPtDau_;
  double maxDeltaRDaus_;
  double minDeltaR_;
  double maxDeltaR_;

  template <typename T>
  float isolationQ(unsigned int pidex1, unsigned int pidex2, const T *cands, unsigned int size) const;

  std::tuple<bool, float> deltar(float eta1, float eta2, float phi1, float phi2) const;

  template <typename T>
  static float pairmass(const std::array<unsigned int, 2> &t, const T *cands, const std::array<float, 2> &massD);

  unsigned long countStruct_;
  unsigned long passStruct_;
};

ScPhase2RecMeson::ScPhase2RecMeson(const edm::ParameterSet &iConfig)
    : doStruct_(iConfig.getParameter<bool>("runStruct")),
      mesonType_(iConfig.getParameter<std::string>("mesonType")),
      minPtDau_(iConfig.getParameter<double>("minPtDau")),
      maxDeltaRDaus_(iConfig.getParameter<double>("maxDeltaRDaus")),
      minDeltaR_(iConfig.getParameter<double>("minDeltaR")),
      maxDeltaR_(iConfig.getParameter<double>("maxDeltaR"))
  {
  if (doStruct_) {
    //PUPPI input being given here
    structToken_ = consumes<OrbitCollection<l1Scouting::Puppi>>(iConfig.getParameter<edm::InputTag>("src"));

    if (mesonType_ == "phi") {
      mesonMassRange_ = {0.95, 1.25};
      dmass1_ = 0.4937;
      dmass2_ = 0.4937;
    } else if (mesonType_ == "rho") {
      mesonMassRange_ = {0.40, 1.30};
      dmass1_ = 0.1396;
      dmass2_ = 0.1396;
    } else if (mesonType_ == "jpsi") {
      mesonMassRange_ = {2.50, 3.50};
      dmass1_ = 0.1057;
      dmass2_ = 0.1057;
    }

    produces<OrbitCollection<l1Scouting::RecMeson>>();
    produces<std::vector<unsigned>>("selectedBx");
    produces<unsigned int>("nbx");
  }
}

ScPhase2RecMeson::~ScPhase2RecMeson() {};

void ScPhase2RecMeson::beginStream(edm::StreamID) {}

void ScPhase2RecMeson::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  if (doStruct_) {
    edm::Handle<OrbitCollection<l1Scouting::Puppi>> src;
    iEvent.getByToken(structToken_, src);

    runObj(*src, iEvent, "");
  }
}

void ScPhase2RecMeson::endStream() {}

template <typename T>
void ScPhase2RecMeson::runObj(const OrbitCollection<T> &src,
                                    edm::Event &iEvent,
                                    const std::string &label) {
  // l1ScoutingRun3::BxOffsetsFillter bxOffsetsFiller;
  // bxOffsetsFiller.start();
  auto ret = std::make_unique<std::vector<unsigned>>();
  auto selectedBx = std::make_unique<std::vector<unsigned>>();

  ROOT::RVec<unsigned int> ix;
  std::vector<std::vector<l1Scouting::RecMeson>> mesonVec;
  unsigned int ntotRecMeson = 0, nbx = 0;

  for (unsigned int bx = 0; bx <= OrbitCollection<T>::NBX; ++bx) {
    nbx++;
    std::vector<l1Scouting::RecMeson> mesonVec_thisBx;

    auto range = src.bxIterator(bx);
    const T *cands = &range.front();
    auto size = range.size();

    ix.clear();
    for (unsigned int i = 0; i < size; ++i) {  //make list of all hadrons
      if ((std::abs(cands[i].pdgId()) == 211 or std::abs(cands[i].pdgId()) == 11)) {
        if (cands[i].pt() >= minPtDau_)
          ix.push_back(i);
      }
    }
    unsigned int ndaus = ix.size();
    //std::cout << "BX = " << bx << " ; number of daugthers = " << ndaus << std::endl;

    std::set<unsigned int> usedIndices;

    for (unsigned int i1 = 0; i1 < ndaus; ++i1) {
      for (unsigned int i2 = i1 + 1; i2 < ndaus; ++i2) {
        if (!(cands[ix[i1]].charge() * cands[ix[i2]].charge() < 0))
          continue;

        auto mass2 = pairmass({{ix[i1], ix[i2]}}, cands, {{dmass1_, dmass2_}});
        if (!(mass2 >= mesonMassRange_[0] and mass2 <= mesonMassRange_[1]))
          continue;

        auto [drcond, drQ] = deltar(cands[ix[i1]].eta(), cands[ix[i2]].eta(), cands[ix[i1]].phi(), cands[ix[i2]].phi());
        if (!drcond)
          continue;  // angular sep of top 2 tracks

        //std::array<unsigned int, 2> pair{{ix[i1], ix[i2]}};  // pair of indices
        //std::cout << "found a meson!!" << std::endl;

        auto p4_1 = ROOT::Math::PtEtaPhiMVector(cands[ix[i1]].pt(), cands[ix[i1]].eta(), cands[ix[i1]].phi(), dmass1_);
        auto p4_2 = ROOT::Math::PtEtaPhiMVector(cands[ix[i2]].pt(), cands[ix[i2]].eta(), cands[ix[i2]].phi(), dmass2_);
        auto recMeson_quad = p4_1 + p4_2;
        // Do we want to put isolation computation here or outside?
        float isoDR0p25 = isolationQ(ix[i1], ix[i2], cands, size);

        // charge set to 0 because of opposite sign condition
        auto recMeson = l1Scouting::RecMeson(recMeson_quad.pt(), recMeson_quad.eta(),
                                             recMeson_quad.phi(), recMeson_quad.mass(),
                                             0, dmass1_, dmass2_, 211, ix[i1], ix[i2],
                                             isoDR0p25);
        mesonVec_thisBx.push_back(recMeson);

        ntotRecMeson++;
        break;
      }
    }

    // if ( mesonVec_thisBx.size() >= 2) {
    //   std::cout << "BX = " << bx << " ; number of mesons = " << mesonVec_thisBx.size() << std::endl;
    //   std::cout << "cand 0 - ids = " << mesonVec_thisBx[0].id1() << " and " << mesonVec_thisBx[0].id2() << std::endl ;
    //   std::cout << "cand 1 - ids = " << mesonVec_thisBx[1].id1() << " and " << mesonVec_thisBx[1].id2() << std::endl << std::endl;
    // }

    mesonVec.push_back(mesonVec_thisBx);
    // if(bx == 1) mesonVec.push_back(mesonVec_thisBx);
    // bxOffsetsFiller.addBx(bx, 1);
  }  // loop on BXs

  // std::cout << "Rec Meson - mesonVec.size() - " << mesonVec.size() << std::endl;
  // std::cout << "Rec Meson - nbx - " << nbx << std::endl;

  // auto bxOffsets = bxOffsetsFiller.done();

  // Put flat table into event
  auto outRecMeson = std::make_unique<OrbitCollection<l1Scouting::RecMeson>>(mesonVec, ntotRecMeson);
  iEvent.put(std::move(outRecMeson));
  iEvent.put(std::make_unique<unsigned int>(nbx), "nbx");
  iEvent.put(std::move(selectedBx), "selectedBx");
}

//TEST functions
template <typename T>
float ScPhase2RecMeson::isolationQ(unsigned int pidex1,
                                   unsigned int pidex2,
                                   const T *cands,
                                   unsigned int size) const {
  float psum = 0;
  float eta = cands[pidex1].eta();  //center cone around leading track
  float phi = cands[pidex1].phi();
  for (unsigned int j = 0u; j < size; ++j) {  //loop over other particles
    if (pidex1 == j or pidex2 == j)
      continue;
    float deta = eta - cands[j].eta(), dphi = ROOT::VecOps::DeltaPhi<float>(phi, cands[j].phi());
    float dr2 = deta * deta + dphi * dphi;
    if (dr2 >= minDeltaR_ && dr2 <= maxDeltaR_)
      psum += cands[j].pt();
  }
  // protect from 0 division?
  return psum / (cands[pidex1].pt() + cands[pidex2].pt());
}

std::tuple<bool, float> ScPhase2RecMeson::deltar(float eta1, float eta2, float phi1, float phi2) const {
  bool passed = true;
  float deta = eta1 - eta2;
  float dphi = ROOT::VecOps::DeltaPhi<float>(phi1, phi2);
  float dr2 = deta * deta + dphi * dphi;
  if (dr2 > maxDeltaRDaus_) {
    passed = false;
  }
  return std::tuple(passed, dr2);
}

template <typename T>
float ScPhase2RecMeson::pairmass(const std::array<unsigned int, 2> &t,
                                       const T *cands,
                                       const std::array<float, 2> &massD) {
  ROOT::Math::PtEtaPhiMVector p1(cands[t[0]].pt(), cands[t[0]].eta(), cands[t[0]].phi(), massD[0]);
  ROOT::Math::PtEtaPhiMVector p2(cands[t[1]].pt(), cands[t[1]].eta(), cands[t[1]].phi(), massD[1]);
  float mass = (p1 + p2).M();
  return mass;
}

void ScPhase2RecMeson::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src");
  desc.add<bool>("runStruct", true);
  desc.add<std::string>("mesonType");
  desc.add<double>("minPtDau");
  desc.add<double>("maxDeltaRDaus");
  desc.add<double>("minDeltaR");
  desc.add<double>("maxDeltaR");

  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScPhase2RecMeson);
