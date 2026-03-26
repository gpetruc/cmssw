#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/L1Scouting/interface/OrbitCollection.h"
#include "DataFormats/L1Scouting/interface/OrbitFlatTable.h"
#include "DataFormats/L1TParticleFlow/interface/L1ScoutingPuppi.h"
#include "DataFormats/L1TParticleFlow/interface/L1ScoutingTTrack.h"
#include "DataFormats/L1TParticleFlow/interface/RecMeson.h"
#include "L1TriggerScouting/Utilities/interface/BxOffsetsFiller.h"

#include <ROOT/RVec.hxx>
#include <Math/Vector4D.h>
#include <Math/GenVector/LorentzVector.h>
#include <Math/GenVector/PtEtaPhiM4D.h>

#include <memory>
#include <algorithm>
#include <array>
#include <iostream>

template <typename T>
class ScPhase2RecMesonAll : public edm::stream::EDProducer<> {
public:
  explicit ScPhase2RecMesonAll(const edm::ParameterSet &);
  ~ScPhase2RecMesonAll() override;
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  void produce(edm::Event &, const edm::EventSetup &) override;
  void runObj(const OrbitCollection<T> &src, edm::Event &out);
  int pdgId(const l1Scouting::Puppi &p);
  int pdgId(const l1Scouting::TTrack &p);

  edm::EDGetTokenT<OrbitCollection<T>> structToken_;

  struct mesonTypeStruct {
    std::string name;
    double minMesonMass;
    double maxMesonMass;
    double dauMass1;
    double dauMass2;
    int pdgId;
    bool muonDaughters;
  };

  std::vector<mesonTypeStruct> mesonTypes_;

  double minDeltaR_;
  double maxDeltaR_;
  double maxDeltaRDaus_;
  double maxDeltaZ_;
  double minPtDau_;
  double maxZIsolation_;

  float isolationQ(unsigned int pidex1, unsigned int pidex2, const T *cands, unsigned int size) const;
};

template <typename T>
ScPhase2RecMesonAll<T>::ScPhase2RecMesonAll(const edm::ParameterSet &iConfig)
    : structToken_(consumes<OrbitCollection<T>>(iConfig.getParameter<edm::InputTag>("src"))),
      minDeltaR_(iConfig.getParameter<double>("minDeltaR")),
      maxDeltaR_(iConfig.getParameter<double>("maxDeltaR")),
      maxDeltaRDaus_(iConfig.getParameter<double>("maxDeltaRDaus")),
      maxDeltaZ_(iConfig.getParameter<double>("maxDeltaZ")),
      minPtDau_(iConfig.getParameter<double>("minPtDau")),
      maxZIsolation_(iConfig.getParameter<double>("maxZIsolation")) {
  std::vector<edm::ParameterSet> mesonPsets = iConfig.getParameter<std::vector<edm::ParameterSet>>("mesonTypes");

  for (const auto &pset : mesonPsets) {
    mesonTypeStruct mt;
    mt.name = pset.getParameter<std::string>("name");
    mt.minMesonMass = pset.getParameter<double>("minMesonMass");
    mt.maxMesonMass = pset.getParameter<double>("maxMesonMass");
    mt.dauMass1 = pset.getParameter<double>("dauMass1");
    mt.dauMass2 = pset.getParameter<double>("dauMass2");
    mt.pdgId = pset.getParameter<int>("pdgId");
    if constexpr (std::is_same_v<T, l1Scouting::TTrack>) {
      // for TTracks, the pdgId can never be the muon one so this parameter is not useful
      mt.muonDaughters = false;
    } else {
      mt.muonDaughters = pset.getParameter<bool>("muonDaughters");
    }
    mesonTypes_.push_back(mt);

    // Register output collections with instance labels
    produces<OrbitCollection<l1Scouting::RecMeson<2>>>(mt.name);
  }
}

template <typename T>
ScPhase2RecMesonAll<T>::~ScPhase2RecMesonAll() {}

template <typename T>
void ScPhase2RecMesonAll<T>::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  edm::Handle<OrbitCollection<T>> src;
  iEvent.getByToken(structToken_, src);

  runObj(*src, iEvent);
}

template <typename T>
void ScPhase2RecMesonAll<T>::runObj(const OrbitCollection<T> &src, edm::Event &iEvent) {
  ROOT::RVec<unsigned int> ix;

  std::vector<std::vector<std::vector<l1Scouting::RecMeson<2>>>> all_mesonVec;
  std::vector<unsigned int> all_ntotRecMeson(mesonTypes_.size(), 0);

  for (unsigned int bx = 0; bx <= OrbitCollection<T>::NBX; ++bx) {
    std::vector<std::vector<l1Scouting::RecMeson<2>>> all_mesonVec_thisBx(mesonTypes_.size());

    auto range = src.bxIterator(bx);
    auto size = range.size();
    const T *cands = (size > 0) ? &range.front() : nullptr;

    ix.clear();
    for (unsigned int i = 0; i < size; ++i) {  //make list of all hadrons
      if ((std::abs(pdgId(cands[i])) == 211 or std::abs(pdgId(cands[i])) == 11 or std::abs(pdgId(cands[i])) == 13)) {
        if (cands[i].pt() >= minPtDau_)
          ix.push_back(i);
      }
    }
    unsigned int ndaus = ix.size();

    for (unsigned int i1 = 0; i1 < ndaus; ++i1) {
      for (unsigned int i2 = i1 + 1; i2 < ndaus; ++i2) {
        if (!(cands[ix[i1]].charge() * cands[ix[i2]].charge() < 0))
          continue;

        float drQ = reco::deltaR(cands[ix[i1]], cands[ix[i2]]);
        if (drQ > maxDeltaRDaus_)
          continue;

        float dZ = abs(cands[ix[i1]].z0() - cands[ix[i2]].z0());
        if (dZ > maxDeltaZ_)
          continue;

        if ((std::abs(pdgId(cands[ix[i1]])) == 13) != (std::abs(pdgId(cands[ix[i2]])) == 13)) {
          // if one daughter is a muon, the other must be a muon as well
          continue;
        }

        float isoDR0p25 = isolationQ(ix[i1], ix[i2], cands, size);

        for (unsigned int itype = 0; itype < mesonTypes_.size(); ++itype) {
          if (mesonTypes_[itype].muonDaughters != (std::abs(pdgId(cands[ix[i1]])) == 13)) {
            // if the meson type requires muon daughters, but these are not muons (or vice versa), skip
            continue;
          }
          auto p4_1 = ROOT::Math::PtEtaPhiMVector(
              cands[ix[i1]].pt(), cands[ix[i1]].eta(), cands[ix[i1]].phi(), mesonTypes_[itype].dauMass1);
          auto p4_2 = ROOT::Math::PtEtaPhiMVector(
              cands[ix[i2]].pt(), cands[ix[i2]].eta(), cands[ix[i2]].phi(), mesonTypes_[itype].dauMass2);
          auto recMeson_quad = p4_1 + p4_2;
          auto mass2 = recMeson_quad.mass();
          if (!(mass2 >= mesonTypes_[itype].minMesonMass and mass2 <= mesonTypes_[itype].maxMesonMass))
            continue;

          std::array<double, 2> daughterMasses = {{mesonTypes_[itype].dauMass1, mesonTypes_[itype].dauMass2}};
          std::array<unsigned int, 2> daughterIds = {{ix[i1], ix[i2]}};

          // charge set to 0 because of opposite sign condition
          auto recMeson = l1Scouting::RecMeson<2>(recMeson_quad.pt(),
                                                  recMeson_quad.eta(),
                                                  recMeson_quad.phi(),
                                                  recMeson_quad.mass(),
                                                  0,
                                                  mesonTypes_[itype].pdgId,
                                                  isoDR0p25,
                                                  daughterMasses,
                                                  daughterIds);

          all_mesonVec_thisBx[itype].push_back(recMeson);

          all_ntotRecMeson[itype]++;
        }
      }
    }
    all_mesonVec.push_back(all_mesonVec_thisBx);
  }

  for (unsigned int itype = 0; itype < mesonTypes_.size(); ++itype) {
    std::vector<std::vector<l1Scouting::RecMeson<2>>> mesonVec_perType;
    mesonVec_perType.reserve(all_mesonVec.size());
    for (auto &bxVec : all_mesonVec) {
      mesonVec_perType.push_back(std::move(bxVec[itype]));
    }
    auto outRecMeson =
        std::make_unique<OrbitCollection<l1Scouting::RecMeson<2>>>(mesonVec_perType, all_ntotRecMeson[itype]);
    iEvent.put(std::move(outRecMeson), mesonTypes_[itype].name);
  }
}

template <typename T>
int ScPhase2RecMesonAll<T>::pdgId(const l1Scouting::Puppi &p) {
  return p.pdgId();
}

template <typename T>
int ScPhase2RecMesonAll<T>::pdgId(const l1Scouting::TTrack &p) {
  return p.charge() * 211;
}

template <typename T>
float ScPhase2RecMesonAll<T>::isolationQ(unsigned int pidex1,
                                         unsigned int pidex2,
                                         const T *cands,
                                         unsigned int size) const {
  float psum = 0;
  auto p4_1 = ROOT::Math::PtEtaPhiMVector(cands[pidex1].pt(), cands[pidex1].eta(), cands[pidex1].phi(), 0.0);
  auto p4_2 = ROOT::Math::PtEtaPhiMVector(cands[pidex2].pt(), cands[pidex2].eta(), cands[pidex2].phi(), 0.0);
  float ptQ = (p4_1 + p4_2).pt();
  float etaQ = (p4_1 + p4_2).eta();
  float phiQ = (p4_1 + p4_2).phi();
  for (unsigned int j = 0u; j < size; ++j) {  //loop over other particles
    if (pidex1 == j or pidex2 == j)
      continue;

    // only consider particles with a small distance in z from the candidates for the isolation calculation
    float z_boson = (cands[pidex1].z0() * cands[pidex1].pt() + cands[pidex2].z0() * cands[pidex2].pt()) /
                    (cands[pidex1].pt() + cands[pidex2].pt());
    if (abs(z_boson - cands[j].z0()) > maxZIsolation_)
      continue;

    float deta = etaQ - cands[j].eta(), dphi = ROOT::VecOps::DeltaPhi<float>(phiQ, cands[j].phi());
    float dr2 = deta * deta + dphi * dphi;
    if (dr2 >= minDeltaR_ && dr2 <= maxDeltaR_)
      psum += cands[j].pt();
  }
  return psum / std::max(ptQ, 0.1f);  // protect from 0 division
}

template <typename T>
void ScPhase2RecMesonAll<T>::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription mesonDesc;
  mesonDesc.add<std::string>("name");
  mesonDesc.add<double>("minMesonMass");
  mesonDesc.add<double>("maxMesonMass");
  mesonDesc.add<double>("dauMass1");
  mesonDesc.add<double>("dauMass2");
  mesonDesc.add<int>("pdgId", 0);
  mesonDesc.add<bool>("muonDaughters", false);

  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src");
  desc.addVPSet("mesonTypes", mesonDesc);
  desc.add<double>("minDeltaR");
  desc.add<double>("maxDeltaR");
  desc.add<double>("maxDeltaRDaus");
  desc.add<double>("maxDeltaZ");
  desc.add<double>("minPtDau");
  desc.add<double>("maxZIsolation");
  descriptions.addDefault(desc);
}

typedef ScPhase2RecMesonAll<l1Scouting::Puppi> ScPhase2PuppiRecMesonAll;
typedef ScPhase2RecMesonAll<l1Scouting::TTrack> ScPhase2TTrackRecMesonAll;

DEFINE_FWK_MODULE(ScPhase2PuppiRecMesonAll);
DEFINE_FWK_MODULE(ScPhase2TTrackRecMesonAll);
