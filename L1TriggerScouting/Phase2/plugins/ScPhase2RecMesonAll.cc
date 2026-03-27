#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/L1Scouting/interface/OrbitCollection.h"
#include "DataFormats/L1Scouting/interface/OrbitFlatTable.h"
#include "DataFormats/L1TParticleFlow/interface/L1ScoutingPuppi.h"
#include "DataFormats/L1TParticleFlow/interface/L1ScoutingTTrack.h"
#include "DataFormats/L1TParticleFlow/interface/L1ScoutingTkEm.h"
#include "DataFormats/L1TMuonPhase2/interface/L1ScoutingTrackerMuon.h"
#include "DataFormats/L1TParticleFlow/interface/RecMeson.h"
#include "L1TriggerScouting/Utilities/interface/BxOffsetsFiller.h"

#include "DataFormats/Math/interface/deltaR.h"
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
  int pdgId(const l1Scouting::TrackerMuon &p);
  int pdgId(const l1Scouting::TkEle &p);

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

  double isolationMinDeltaR2_;
  double isolationMaxDeltaR2_;
  double isolationMaxDeltaZ_;
  double isolationMaxRelIso_;
  double maxDeltaR2Daus_;
  double maxDeltaZ_;
  double minPtDau_;

  float isolationQ(unsigned int pidex1, unsigned int pidex2, const T *cands, unsigned int size) const;
};

template <typename T>
ScPhase2RecMesonAll<T>::ScPhase2RecMesonAll(const edm::ParameterSet &iConfig)
    : structToken_(consumes<OrbitCollection<T>>(iConfig.getParameter<edm::InputTag>("src"))),
      isolationMinDeltaR2_(std::pow(iConfig.getParameter<double>("isolationMinDeltaR"), 2)),
      isolationMaxDeltaR2_(std::pow(iConfig.getParameter<double>("isolationMaxDeltaR"), 2)),
      isolationMaxDeltaZ_(iConfig.getParameter<double>("isolationMaxDeltaZ")),
      isolationMaxRelIso_(iConfig.getParameter<double>("isolationMaxRelIso")),
      maxDeltaR2Daus_(std::pow(iConfig.getParameter<double>("maxDeltaRDaus"), 2)),
      maxDeltaZ_(iConfig.getParameter<double>("maxDeltaZ")),
      minPtDau_(iConfig.getParameter<double>("minPtDau")) {
  std::vector<edm::ParameterSet> mesonPsets = iConfig.getParameter<std::vector<edm::ParameterSet>>("mesonTypes");

  for (const auto &pset : mesonPsets) {
    mesonTypeStruct mt;
    mt.name = pset.getParameter<std::string>("name");
    mt.minMesonMass = pset.getParameter<double>("minMesonMass");
    mt.maxMesonMass = pset.getParameter<double>("maxMesonMass");
    mt.dauMass1 = pset.getParameter<double>("dauMass1");
    mt.dauMass2 = pset.getParameter<double>("dauMass2");
    mt.pdgId = pset.getParameter<int>("pdgId");
    if constexpr (std::is_same_v<T, l1Scouting::Puppi>) {
      // only for Puppi we have meaningful PDG IDs
      mt.muonDaughters = pset.getParameter<bool>("muonDaughters");
    } else {
      mt.muonDaughters = true;
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

        float dr2Q = reco::deltaR2(cands[ix[i1]], cands[ix[i2]]);
        if (dr2Q > maxDeltaR2Daus_)
          continue;

        float dZ = std::abs(cands[ix[i1]].z0() - cands[ix[i2]].z0());
        if (dZ > maxDeltaZ_)
          continue;

        if ((std::abs(pdgId(cands[ix[i1]])) == 13) != (std::abs(pdgId(cands[ix[i2]])) == 13)) {
          // if one daughter is a muon, the other must be a muon as well
          continue;
        }

        float isoDR0p25 = isolationQ(ix[i1], ix[i2], cands, size);
        if (isoDR0p25 > isolationMaxRelIso_)
          continue;

        for (unsigned int itype = 0; itype < mesonTypes_.size(); ++itype) {
          if ((std::abs(pdgId(cands[ix[i1]])) == 13) && !mesonTypes_[itype].muonDaughters) {
            // if the meson type requires non-muon daughters, but these are muons, skip
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
int ScPhase2RecMesonAll<T>::pdgId(const l1Scouting::TrackerMuon &p) {
  return p.charge() * (-13);
}

template <typename T>
int ScPhase2RecMesonAll<T>::pdgId(const l1Scouting::TkEle &p) {
  return p.charge() * (-11);
}

template <typename T>
float ScPhase2RecMesonAll<T>::isolationQ(unsigned int pidex1,
                                         unsigned int pidex2,
                                         const T *cands,
                                         unsigned int size) const {
  float psum = 0;
  auto p4_1 = ROOT::Math::PtEtaPhiMVector(cands[pidex1].pt(), cands[pidex1].eta(), cands[pidex1].phi(), 0.0);
  auto p4_2 = ROOT::Math::PtEtaPhiMVector(cands[pidex2].pt(), cands[pidex2].eta(), cands[pidex2].phi(), 0.0);
  float htQ = cands[pidex1].pt() + cands[pidex2].pt();
  float etaQ = (p4_1 + p4_2).eta();
  float phiQ = (p4_1 + p4_2).phi();
  // only consider particles with a small distance in z from the candidates for the isolation calculation
  float z_pair = (cands[pidex1].z0() * cands[pidex1].pt() + cands[pidex2].z0() * cands[pidex2].pt()) / htQ;
  for (unsigned int j = 0u; j < size; ++j) {  //loop over other particles
    if (pidex1 == j or pidex2 == j)
      continue;

    if (cands[j].charge() != 0 && std::abs(z_pair - cands[j].z0()) > isolationMaxDeltaZ_)
      continue;

    float dr2 = reco::deltaR2(etaQ, phiQ, cands[j].eta(), cands[j].phi());
    if (dr2 >= isolationMinDeltaR2_ && dr2 <= isolationMaxDeltaR2_)
      psum += cands[j].pt();
  }
  return psum / htQ;  // htQ can never be zero
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
  desc.add<double>("isolationMinDeltaR");
  desc.add<double>("isolationMaxDeltaR");
  desc.add<double>("isolationMaxDeltaZ");
  desc.add<double>("isolationMaxRelIso");
  desc.add<double>("maxDeltaRDaus");
  desc.add<double>("maxDeltaZ");
  desc.add<double>("minPtDau");
  descriptions.addDefault(desc);
}

typedef ScPhase2RecMesonAll<l1Scouting::Puppi> ScPhase2PuppiRecMesonAll;
typedef ScPhase2RecMesonAll<l1Scouting::TTrack> ScPhase2TTrackRecMesonAll;
typedef ScPhase2RecMesonAll<l1Scouting::TrackerMuon> ScPhase2TrackerMuonRecMesonAll;
typedef ScPhase2RecMesonAll<l1Scouting::TkEle> ScPhase2TkEleRecMesonAll;

DEFINE_FWK_MODULE(ScPhase2PuppiRecMesonAll);
DEFINE_FWK_MODULE(ScPhase2TTrackRecMesonAll);
DEFINE_FWK_MODULE(ScPhase2TrackerMuonRecMesonAll);
DEFINE_FWK_MODULE(ScPhase2TkEleRecMesonAll);
