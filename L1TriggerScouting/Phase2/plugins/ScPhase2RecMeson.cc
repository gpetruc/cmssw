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
              unsigned long &nTry,
              unsigned long &nPass,
              const std::string &bxLabel);

  bool doStruct_;
  edm::EDGetTokenT<OrbitCollection<l1Scouting::Puppi>> structToken_;
  std::string mesonType_;
  std::vector<float> mesonMassRange_ = {0.0f, 0.0f};

  struct Cuts {
    float minptD = 1;
    float maxdeltarD2 = 0.40 * 0.40;
    float mindr2 = 0.05 * 0.05;
    float maxdr2 = 0.25 * 0.25;
    float maxiso = 0.25;
  } cuts;

  template <typename T>
  bool isolationQ(unsigned int pidex1, unsigned int pidex2, const T *cands, unsigned int size) const;

  std::tuple<bool, float> deltar(float eta1, float eta2, float phi1, float phi2) const;

  template <typename T>
  static float pairmass(const std::array<unsigned int, 2> &t, const T *cands, const std::array<float, 2> &massD);

  template <typename T>
  static float quadrupletmass(const std::array<unsigned int, 4> &t, const T *cands, const std::array<float, 4> &massD);

  unsigned long countStruct_;
  unsigned long passStruct_;
};

ScPhase2RecMeson::ScPhase2RecMeson(const edm::ParameterSet &iConfig)
    : doStruct_(iConfig.getParameter<bool>("runStruct")),
      mesonType_(iConfig.getParameter<std::string>("mesonType")) {
  if (doStruct_) {
    //PUPPI input being given here
    structToken_ = consumes<OrbitCollection<l1Scouting::Puppi>>(iConfig.getParameter<edm::InputTag>("src"));

    if (mesonType_ == "phi") {
      mesonMassRange_ = {0.95, 1.25};
    } else if (mesonType_ == "rho") {
      mesonMassRange_ = {0.40, 1.30};
    }

    produces<OrbitCollection<l1Scouting::RecMeson>>();
    produces<std::vector<unsigned>>("selectedBx");
    produces<unsigned int>("nbx");
  }
}

ScPhase2RecMeson::~ScPhase2RecMeson() {};

void ScPhase2RecMeson::beginStream(edm::StreamID) {
  countStruct_ = 0;
  passStruct_ = 0;
}

void ScPhase2RecMeson::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  if (doStruct_) {
    edm::Handle<OrbitCollection<l1Scouting::Puppi>> src;
    iEvent.getByToken(structToken_, src);

    runObj(*src, iEvent, countStruct_, passStruct_, "");
  }
}

void ScPhase2RecMeson::endStream() {
  if (doStruct_)
    edm::LogImportant("ScPhase2AnalysisSummary") << "RecMeson Struct analysis: " << countStruct_ << " -> " << passStruct_;
}

template <typename T>
void ScPhase2RecMeson::runObj(const OrbitCollection<T> &src,
                                    edm::Event &iEvent,
                                    unsigned long &nTry,
                                    unsigned long &nPass,
                                    const std::string &label) {
  l1ScoutingRun3::BxOffsetsFillter bxOffsetsFiller;
  bxOffsetsFiller.start();
  auto ret = std::make_unique<std::vector<unsigned>>();
  auto selectedBx = std::make_unique<std::vector<unsigned>>();

  ROOT::RVec<unsigned int> ix;
  std::vector<std::vector<l1Scouting::RecMeson>> mesonVec;
  unsigned int ntotRecMeson = 0, nbx = 0;

  for (unsigned int bx = 1; bx <= OrbitCollection<T>::NBX; ++bx) {
    nbx++;
    std::vector<l1Scouting::RecMeson> mesonVec_thisBx;

    nTry++;
    auto range = src.bxIterator(bx);
    const T *cands = &range.front();
    auto size = range.size();

    ix.clear();
    for (unsigned int i = 0; i < size; ++i) {  //make list of all hadrons
      if ((std::abs(cands[i].pdgId()) == 211 or std::abs(cands[i].pdgId()) == 11)) {
        if (cands[i].pt() >= cuts.minptD)
          ix.push_back(i);
      }
    }
    unsigned int ndaus = ix.size();
    //std::cout << "BX = " << bx << " ; number of daugthers = " << ndaus << std::endl;
    
    for (unsigned int i1 = 0; i1 < ndaus; ++i1) {
      if (cands[ix[i1]].pt() < cuts.minptD)
        continue;  // D1 pt cut
      for (unsigned int i2 = 0; i2 < ndaus; ++i2) {
        if (i2 == i1 || cands[ix[i2]].pt() < cuts.minptD)
          continue;  // D2 pt cut

        if (!(cands[ix[i1]].charge() * cands[ix[i2]].charge() < 0))
          continue;

        auto mass2 = pairmass({{ix[i1], ix[i2]}}, cands, {{0.4937, 0.4937}});
        if (!(mass2 >= mesonMassRange_[0] and mass2 <= mesonMassRange_[1]))
          continue;

        auto [drcond, drQ] = deltar(cands[ix[i1]].eta(), cands[ix[i2]].eta(), cands[ix[i1]].phi(), cands[ix[i2]].phi());
        if (!drcond)
          continue;  // angular sep of top 2 tracks

        //std::array<unsigned int, 2> pair{{ix[i1], ix[i2]}};  // pair of indices
        //std::cout << "found a meson!!" << std::endl;

        auto p4_1 = cands[ix[i1]].p4();
        auto p4_2 = cands[ix[i2]].p4();
        auto recMeson_quad = p4_1 + p4_2;  
        
        auto recMeson = l1Scouting::RecMeson(recMeson_quad.pt(), recMeson_quad.eta(), recMeson_quad.phi(), i1, i2);
        mesonVec_thisBx.push_back(recMeson);
        ntotRecMeson++;
      }
    }

    //std::cout << "BX = " << bx << " ; number of mesons = " << mesonVec_thisBx.size() << std::endl;
    mesonVec.push_back(mesonVec_thisBx);
    if(bx == 1) mesonVec.push_back(mesonVec_thisBx);
    bxOffsetsFiller.addBx(bx, 1);
  }  // loop on BXs

  // std::cout << "Rec Meson - mesonVec.size() - " << mesonVec.size() << std::endl;
  // std::cout << "Rec Meson - nbx - " << nbx << std::endl;

  auto bxOffsets = bxOffsetsFiller.done();

  // Put flat table into event
  auto outRecMeson = std::make_unique<OrbitCollection<l1Scouting::RecMeson>>(mesonVec, ntotRecMeson);
  iEvent.put(std::move(outRecMeson));
  iEvent.put(std::make_unique<unsigned int>(nbx), "nbx");
  iEvent.put(std::move(selectedBx), "selectedBx");
}

//TEST functions
template <typename T>
bool ScPhase2RecMeson::isolationQ(unsigned int pidex1,
                                        unsigned int pidex2,
                                        const T *cands,
                                        unsigned int size) const {
  bool passed = false;
  float psum = 0;
  float eta = cands[pidex1].eta();  //center cone around leading track
  float phi = cands[pidex1].phi();
  for (unsigned int j = 0u; j < size; ++j) {  //loop over other particles
    if (pidex1 == j or pidex2 == j)
      continue;
    float deta = eta - cands[j].eta(), dphi = ROOT::VecOps::DeltaPhi<float>(phi, cands[j].phi());
    float dr2 = deta * deta + dphi * dphi;
    if (dr2 >= cuts.mindr2 && dr2 <= cuts.maxdr2)
      psum += cands[j].pt();
  }
  if (psum <= cuts.maxiso * (cands[pidex1].pt() + cands[pidex2].pt()))
    passed = true;
  return passed;
}

std::tuple<bool, float> ScPhase2RecMeson::deltar(float eta1, float eta2, float phi1, float phi2) const {
  bool passed = true;
  float deta = eta1 - eta2;
  float dphi = ROOT::VecOps::DeltaPhi<float>(phi1, phi2);
  float dr2 = deta * deta + dphi * dphi;
  if (dr2 > cuts.maxdeltarD2) {
    passed = false;
    return std::tuple(passed, dr2);
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

template <typename T>
float ScPhase2RecMeson::quadrupletmass(const std::array<unsigned int, 4> &t,
                                             const T *cands,
                                             const std::array<float, 4> &massD) {
  ROOT::Math::PtEtaPhiMVector p1(cands[t[0]].pt(), cands[t[0]].eta(), cands[t[0]].phi(), massD[0]);
  ROOT::Math::PtEtaPhiMVector p2(cands[t[1]].pt(), cands[t[1]].eta(), cands[t[1]].phi(), massD[1]);
  ROOT::Math::PtEtaPhiMVector p3(cands[t[2]].pt(), cands[t[2]].eta(), cands[t[2]].phi(), massD[2]);
  ROOT::Math::PtEtaPhiMVector p4(cands[t[3]].pt(), cands[t[3]].eta(), cands[t[3]].phi(), massD[3]);
  float mass = (p1 + p2 + p3 + p4).M();
  return mass;
}

void ScPhase2RecMeson::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src");
  desc.add<bool>("runStruct", true);
  desc.add<std::string>("mesonType");
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScPhase2RecMeson);
