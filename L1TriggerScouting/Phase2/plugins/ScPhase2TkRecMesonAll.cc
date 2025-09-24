#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/L1Scouting/interface/OrbitCollection.h"
#include "DataFormats/L1Scouting/interface/OrbitFlatTable.h"
#include "DataFormats/L1TParticleFlow/interface/L1ScoutingTTrack.h"
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

class ScPhase2TkRecMesonAll : public edm::stream::EDProducer<> {
public:
  explicit ScPhase2TkRecMesonAll(const edm::ParameterSet &);
  ~ScPhase2TkRecMesonAll() override;
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
  edm::EDGetTokenT<OrbitCollection<l1Scouting::TTrack>> structToken_;

  struct mesonTypeStruct {
    std::string name;
    double minMesonMass;
    double maxMesonMass;
    double dmass1;
    double dmass2;
  };

  std::vector<mesonTypeStruct> mesonTypes_;

  double minDeltaR_;
  double maxDeltaR_;
  double maxDeltaRDaus_;
  double maxDeltaZ_;
  double minPtDau_;
  double maxZIsolation_;

  template <typename T>
  float isolationQ(int itype, unsigned int pidex1, unsigned int pidex2, const T *cands, unsigned int size) const;

  float deltar(float eta1, float eta2, float phi1, float phi2) const;

  template <typename T>
  static float pairmass(const std::array<unsigned int, 2> &t, const T *cands, const std::array<double, 2> &massD);

  unsigned long countStruct_;
  unsigned long passStruct_;
};

ScPhase2TkRecMesonAll::ScPhase2TkRecMesonAll(const edm::ParameterSet &iConfig)
    : doStruct_(iConfig.getParameter<bool>("runStruct")),
      minDeltaR_(iConfig.getParameter<double>("minDeltaR")),
      maxDeltaR_(iConfig.getParameter<double>("maxDeltaR")),
      maxDeltaRDaus_(iConfig.getParameter<double>("maxDeltaRDaus")),
      maxDeltaZ_(iConfig.getParameter<double>("maxDeltaZ")),
      minPtDau_(iConfig.getParameter<double>("minPtDau")),
      maxZIsolation_(iConfig.getParameter<double>("maxZIsolation"))
  {
  if (doStruct_) {
    //TTrack input being given here
    structToken_ = consumes<OrbitCollection<l1Scouting::TTrack>>(iConfig.getParameter<edm::InputTag>("src"));

    std::vector<edm::ParameterSet> mesonPsets =
        iConfig.getParameter<std::vector<edm::ParameterSet>>("mesonTypes");

    for (const auto& pset : mesonPsets) {
      mesonTypeStruct mt;
      mt.name = pset.getParameter<std::string>("name");
      mt.minMesonMass = pset.getParameter<double>("minMesonMass");
      mt.maxMesonMass = pset.getParameter<double>("maxMesonMass");
      mt.dmass1 = pset.getParameter<double>("dmass1");
      mt.dmass2 = pset.getParameter<double>("dmass2");
      mesonTypes_.push_back(mt);

      // Register output collections with instance labels
      produces<OrbitCollection<l1Scouting::RecMeson>>(mt.name);
    }
  }
}

ScPhase2TkRecMesonAll::~ScPhase2TkRecMesonAll() {};

void ScPhase2TkRecMesonAll::beginStream(edm::StreamID) {}

void ScPhase2TkRecMesonAll::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  if (doStruct_) {
    edm::Handle<OrbitCollection<l1Scouting::TTrack>> src;
    iEvent.getByToken(structToken_, src);

    runObj(*src, iEvent, "");
  }
}

void ScPhase2TkRecMesonAll::endStream() {}

template <typename T>
void ScPhase2TkRecMesonAll::runObj(const OrbitCollection<T> &src,
                                    edm::Event &iEvent,
                                    const std::string &label) {
  // l1ScoutingRun3::BxOffsetsFillter bxOffsetsFiller;
  // bxOffsetsFiller.start();
  auto ret = std::make_unique<std::vector<unsigned>>();

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
    for (unsigned int i = 0; i < size; ++i) {  //make list of all tracks
      if (cands[i].pt() >= minPtDau_)
        ix.push_back(i);
    }
    unsigned int ndaus = ix.size();
    
    for (unsigned int i1 = 0; i1 < ndaus; ++i1) {
      for (unsigned int i2 = i1 + 1; i2 < ndaus; ++i2) {
        if (!(cands[ix[i1]].charge() * cands[ix[i2]].charge() < 0))
          continue;

        float drQ = deltar(cands[ix[i1]].eta(), cands[ix[i2]].eta(), cands[ix[i1]].phi(), cands[ix[i2]].phi());
        if (drQ > maxDeltaRDaus_) 
          continue;

        float dZ = abs(cands[ix[i1]].z0() - cands[ix[i2]].z0());
        if (dZ > maxDeltaZ_) 
          continue;

        float isoDR0p25 = isolationQ(0, ix[i1], ix[i2], cands, size);

        for (unsigned int itype = 0; itype < mesonTypes_.size(); ++itype) {

          auto mass2 = pairmass({{ix[i1], ix[i2]}}, cands, {{mesonTypes_[itype].dmass1, mesonTypes_[itype].dmass2}});
          if (!(mass2 >= mesonTypes_[itype].minMesonMass and mass2 <= mesonTypes_[itype].maxMesonMass))
            continue;

          auto p4_1 = ROOT::Math::PtEtaPhiMVector(cands[ix[i1]].pt(), cands[ix[i1]].eta(), cands[ix[i1]].phi(), mesonTypes_[itype].dmass1);
          auto p4_2 = ROOT::Math::PtEtaPhiMVector(cands[ix[i2]].pt(), cands[ix[i2]].eta(), cands[ix[i2]].phi(), mesonTypes_[itype].dmass2);
          auto recMeson_quad = p4_1 + p4_2;

          auto recMeson = l1Scouting::RecMeson(recMeson_quad.pt(), recMeson_quad.eta(),
                                              recMeson_quad.phi(), recMeson_quad.mass(),
                                              0, mesonTypes_[itype].dmass1, mesonTypes_[itype].dmass2, 211, ix[i1], ix[i2],
                                              isoDR0p25);

          all_mesonVec_thisBx[itype].push_back(recMeson);

          all_ntotRecMeson[itype]++;
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
    iEvent.put(std::move(outRecMeson), mesonTypes_[itype].name);
  }
}

template <typename T>
float ScPhase2TkRecMesonAll::isolationQ(int itype, unsigned int pidex1,
                                   unsigned int pidex2,
                                   const T *cands,
                                   unsigned int size) const {
  float psum = 0;

  auto p4_1 = ROOT::Math::PtEtaPhiMVector(cands[pidex1].pt(), cands[pidex1].eta(), cands[pidex1].phi(), mesonTypes_[itype].dmass1);
  auto p4_2 = ROOT::Math::PtEtaPhiMVector(cands[pidex2].pt(), cands[pidex2].eta(), cands[pidex2].phi(), mesonTypes_[itype].dmass2);
  float ptQ = (p4_1 + p4_2).pt();
  float etaQ = (p4_1 + p4_2).eta();
  float phiQ = (p4_1 + p4_2).phi();
  for (unsigned int j = 0u; j < size; ++j) {  //loop over other particles
    if (pidex1 == j or pidex2 == j)
      continue;

    //only consider particles with a small distance in z from the candidates for the isolation calculation
    float z_boson = (cands[pidex1].z0()*cands[pidex1].pt() + cands[pidex2].z0()*cands[pidex2].pt())/(cands[pidex1].pt() + cands[pidex2].pt());
    if (abs(z_boson - cands[j].z0()) > maxZIsolation_) 
      continue;


    float deta = etaQ - cands[j].eta(), dphi = ROOT::VecOps::DeltaPhi<float>(phiQ, cands[j].phi());
    float dr2 = deta * deta + dphi * dphi;
    if (dr2 >= minDeltaR_ && dr2 <= maxDeltaR_)
      psum += cands[j].pt();
  }
  // protect from 0 division?
  return psum / ptQ;
}

float ScPhase2TkRecMesonAll::deltar(float eta1, float eta2, float phi1, float phi2) const {
  float deta = eta1 - eta2;
  float dphi = ROOT::VecOps::DeltaPhi<float>(phi1, phi2);
  float dr2 = deta * deta + dphi * dphi;

  return dr2;
}

template <typename T>
float ScPhase2TkRecMesonAll::pairmass(const std::array<unsigned int, 2> &t,
                                    const T *cands,
                                    const std::array<double, 2> &massD) {
  ROOT::Math::PtEtaPhiMVector p1(cands[t[0]].pt(), cands[t[0]].eta(), cands[t[0]].phi(), massD[0]);
  ROOT::Math::PtEtaPhiMVector p2(cands[t[1]].pt(), cands[t[1]].eta(), cands[t[1]].phi(), massD[1]);
  float mass = (p1 + p2).M();
  return mass;
}

void ScPhase2TkRecMesonAll::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription mesonDesc;
  mesonDesc.add<std::string>("name");
  mesonDesc.add<double>("minMesonMass");
  mesonDesc.add<double>("maxMesonMass");
  mesonDesc.add<double>("dmass1");
  mesonDesc.add<double>("dmass2");

  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src");
  desc.add<bool>("runStruct", true);
  desc.addVPSet("mesonTypes", mesonDesc);
  desc.add<double>("minDeltaR");
  desc.add<double>("maxDeltaR");
  desc.add<double>("maxDeltaRDaus");
  desc.add<double>("maxDeltaZ");
  desc.add<double>("minPtDau");
  desc.add<double>("maxZIsolation");
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScPhase2TkRecMesonAll);
