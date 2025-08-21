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

class ScPhase2RecMesonAll : public edm::stream::EDProducer<> {
public:
  explicit ScPhase2RecMesonAll(const edm::ParameterSet &);
  ~ScPhase2RecMesonAll() override;
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
  std::vector<std::string> mesonTypes_;

  std::vector<float> minDeltaR_;
  std::vector<float> maxDeltaR_;
  std::vector<float> maxDeltaRDaus_;
  float minPtDau_ = 5.0;
  std::vector<std::array<float, 2>> massRange_;
  std::vector<float> dmass1_;
  std::vector<float> dmass2_;

  template <typename T>
  float isolationQ(int itype, unsigned int pidex1, unsigned int pidex2, const T *cands, unsigned int size) const;

  float deltar(float eta1, float eta2, float phi1, float phi2) const;

  template <typename T>
  static float pairmass(const std::array<unsigned int, 2> &t, const T *cands, const std::array<float, 2> &massD);

  unsigned long countStruct_;
  unsigned long passStruct_;
};

ScPhase2RecMesonAll::ScPhase2RecMesonAll(const edm::ParameterSet &iConfig)
    : doStruct_(iConfig.getParameter<bool>("runStruct")),
      mesonTypes_(iConfig.getParameter<std::vector<std::string>>("mesonTypes"))
  {
  if (doStruct_) {
    //PUPPI input being given here
    structToken_ = consumes<OrbitCollection<l1Scouting::Puppi>>(iConfig.getParameter<edm::InputTag>("src"));
    
    for (const auto &mt : mesonTypes_) {
    
      if (mt == "phi") {
        minDeltaR_.push_back(0.05 * 0.05);
        maxDeltaR_.push_back(0.25 * 0.25);
        maxDeltaRDaus_.push_back(0.40 * 0.40);
        massRange_.push_back({{0.95, 1.25}});
        dmass1_.push_back(0.4937);
        dmass2_.push_back(0.4937);
      } else if (mt == "rho") {
        minDeltaR_.push_back(0.05 * 0.05);
        maxDeltaR_.push_back(0.25 * 0.25);
        maxDeltaRDaus_.push_back(0.40 * 0.40);
        massRange_.push_back({{0.40, 1.30}});
        dmass1_.push_back(0.1396);
        dmass2_.push_back(0.1396);
      } else if (mt == "jpsi") {
        minDeltaR_.push_back(0.05 * 0.05);
        maxDeltaR_.push_back(0.25 * 0.25);
        maxDeltaRDaus_.push_back(0.40 * 0.40);
        massRange_.push_back({{2.50, 3.50}});
        dmass1_.push_back(0.1057);
        dmass2_.push_back(0.1057);
      }
    
      produces<OrbitCollection<l1Scouting::RecMeson>>("recMeson" + mt);
    }

    produces<std::vector<unsigned>>("selectedBx");
    produces<unsigned int>("nbx");
  }
}

ScPhase2RecMesonAll::~ScPhase2RecMesonAll() {};

void ScPhase2RecMesonAll::beginStream(edm::StreamID) {}

void ScPhase2RecMesonAll::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  if (doStruct_) {
    edm::Handle<OrbitCollection<l1Scouting::Puppi>> src;
    iEvent.getByToken(structToken_, src);

    runObj(*src, iEvent, "");
  }
}

void ScPhase2RecMesonAll::endStream() {}

template <typename T>
void ScPhase2RecMesonAll::runObj(const OrbitCollection<T> &src,
                                    edm::Event &iEvent,
                                    const std::string &label) {
  // l1ScoutingRun3::BxOffsetsFillter bxOffsetsFiller;
  // bxOffsetsFiller.start();
  auto ret = std::make_unique<std::vector<unsigned>>();
  auto selectedBx = std::make_unique<std::vector<unsigned>>();

  ROOT::RVec<unsigned int> ix;

  std::vector<std::vector<std::vector<l1Scouting::RecMeson>>> all_mesonVec;
  std::vector<unsigned int> all_ntotRecMeson(mesonTypes_.size(), 0);;
  int nbx = 0;

  for (unsigned int bx = 0; bx <= OrbitCollection<T>::NBX; ++bx) {
    nbx++;
    std::vector<std::vector<l1Scouting::RecMeson>> all_mesonVec_thisBx(mesonTypes_.size());

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

        for (unsigned int itype = 0; itype < mesonTypes_.size(); ++itype) {

          auto mass2 = pairmass({{ix[i1], ix[i2]}}, cands, {{dmass1_[itype], dmass2_[itype]}});
          if (!(mass2 >= massRange_[itype][0] and mass2 <= massRange_[itype][1]))
            continue;

          float drQ = deltar(cands[ix[i1]].eta(), cands[ix[i2]].eta(), cands[ix[i1]].phi(), cands[ix[i2]].phi());
          if (drQ > maxDeltaRDaus_[itype]) 
            continue;
  
          float isoDR0p25 = isolationQ(itype, ix[i1], ix[i2], cands, size);
  
          auto p4_1 = ROOT::Math::PtEtaPhiMVector(cands[ix[i1]].pt(), cands[ix[i1]].eta(), cands[ix[i1]].phi(), dmass1_[itype]);
          auto p4_2 = ROOT::Math::PtEtaPhiMVector(cands[ix[i2]].pt(), cands[ix[i2]].eta(), cands[ix[i2]].phi(), dmass2_[itype]);
          auto recMeson_quad = p4_1 + p4_2;

          auto recMeson = l1Scouting::RecMeson(recMeson_quad.pt(), recMeson_quad.eta(),
                                              recMeson_quad.phi(), recMeson_quad.mass(),
                                              0, dmass1_[itype], dmass2_[itype], 211, ix[i1], ix[i2],
                                              isoDR0p25);

          all_mesonVec_thisBx[itype].push_back(recMeson);

          all_ntotRecMeson[itype]++;
          break;
        }
      }
    }
    all_mesonVec.push_back(all_mesonVec_thisBx);
  }  

  for (unsigned int itype = 0; itype < mesonTypes_.size(); ++itype) {
    std::vector<std::vector<l1Scouting::RecMeson>> mesonVec_perType;
    for (auto &bxVec : all_mesonVec) {
        mesonVec_perType.push_back(std::move(bxVec[itype]));
    }
    auto outRecMeson = std::make_unique<OrbitCollection<l1Scouting::RecMeson>>(
        mesonVec_perType, all_ntotRecMeson[itype]
    );
    iEvent.put(std::move(outRecMeson), "recMeson" + mesonTypes_[itype]);
  }
  iEvent.put(std::make_unique<unsigned int>(nbx), "nbx");
  iEvent.put(std::move(selectedBx), "selectedBx");
}

//TEST functions
template <typename T>
float ScPhase2RecMesonAll::isolationQ(int itype, unsigned int pidex1,
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
    if (dr2 >= minDeltaR_[0] && dr2 <= maxDeltaR_[0])
      psum += cands[j].pt();
  }
  // protect from 0 division?
  return psum / (cands[pidex1].pt() + cands[pidex2].pt());
}

float ScPhase2RecMesonAll::deltar(float eta1, float eta2, float phi1, float phi2) const {
  float deta = eta1 - eta2;
  float dphi = ROOT::VecOps::DeltaPhi<float>(phi1, phi2);
  float dr2 = deta * deta + dphi * dphi;

  return dr2;
}

template <typename T>
float ScPhase2RecMesonAll::pairmass(const std::array<unsigned int, 2> &t,
                                       const T *cands,
                                       const std::array<float, 2> &massD) {
  ROOT::Math::PtEtaPhiMVector p1(cands[t[0]].pt(), cands[t[0]].eta(), cands[t[0]].phi(), massD[0]);
  ROOT::Math::PtEtaPhiMVector p2(cands[t[1]].pt(), cands[t[1]].eta(), cands[t[1]].phi(), massD[1]);
  float mass = (p1 + p2).M();
  return mass;
}

void ScPhase2RecMesonAll::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src");
  desc.add<bool>("runStruct", true);
  desc.add<std::vector<std::string>>("mesonTypes");
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScPhase2RecMesonAll);


// process.scoutingMesons = cms.EDProducer("ScPhase2RecMesonAll",
//   src = cms.InputTag("l1tPuppi"),
//   runStruct = cms.bool(True),
//   mesonTypes = cms.vstring("phi", "rho", "jpsi")
// )