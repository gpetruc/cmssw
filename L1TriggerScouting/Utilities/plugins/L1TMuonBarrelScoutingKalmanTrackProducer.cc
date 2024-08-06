#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "L1Trigger/L1TMuonBarrel/interface/L1TMuonBarrelKalmanAlgo.h"
#include "L1Trigger/L1TMuonBarrel/interface/L1TMuonBarrelKalmanTrackFinder.h"
#include "L1Trigger/L1TMuonBarrel/interface/L1TMuonBarrelKalmanStubProcessor.h"

#include "DataFormats/L1TMuon/interface/L1MuKBMTCombinedStub.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCandFwd.h"

#include "DataFormats/L1Scouting/interface/OrbitCollection.h"
#include "DataFormats/L1Scouting/interface/L1ScoutingMuon.h"
#include "DataFormats/L1Scouting/interface/L1ScoutingBMTFStub.h"

#include "L1TriggerScouting/Utilities/interface/conversion.h"
#include "L1TriggerScouting/Utilities/interface/convertToL1TFormat.h"

using namespace l1ScoutingRun3;

//
// class declaration
//

class L1TMuonBarrelScoutingKalmanTrackProducer : public edm::stream::EDProducer<> {
public:
  explicit L1TMuonBarrelScoutingKalmanTrackProducer(const edm::ParameterSet&);
  ~L1TMuonBarrelScoutingKalmanTrackProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  unsigned int calcGlobalPhi(const l1t::RegionalMuonCand&);
  double calcDr(const l1t::RegionalMuonCand&, const l1t::Muon&);

private:
  void beginStream(edm::StreamID) override;
  void produce(edm::Event&, const edm::EventSetup&) override;
  void endStream() override;

  edm::EDGetTokenT<L1MuKBMTCombinedStubCollection> src_;
  edm::EDGetTokenT<MuonOrbitCollection> gmtSrc_;
  L1TMuonBarrelKalmanAlgo* algo_;
  L1TMuonBarrelKalmanTrackFinder* trackFinder_;
  bool matchGmt_;
  double drCut_;
  double phiMult_;
  double etaMult_;
  bool debug_;
};
L1TMuonBarrelScoutingKalmanTrackProducer::L1TMuonBarrelScoutingKalmanTrackProducer(const edm::ParameterSet& iConfig)
    : src_(consumes<L1MuKBMTCombinedStubCollection>(iConfig.getParameter<edm::InputTag>("src"))),
      gmtSrc_(consumes<MuonOrbitCollection>(iConfig.getParameter<edm::InputTag>("gmtSrc"))),
      algo_(new L1TMuonBarrelKalmanAlgo(iConfig.getParameter<edm::ParameterSet>("algoSettings"))),
      trackFinder_(new L1TMuonBarrelKalmanTrackFinder(iConfig.getParameter<edm::ParameterSet>("trackFinderSettings"))),
      matchGmt_(iConfig.getParameter<bool>("matchGmt")),
      drCut_(iConfig.getParameter<double>("drCut")),
      phiMult_(iConfig.getParameter<double>("phiMult")),
      etaMult_(iConfig.getParameter<double>("etaMult")),
      debug_(iConfig.getParameter<bool>("debug")) {
  // produces<L1MuKBMTrackBxCollection>("BMTF");
  produces<L1MuKBMTrackOrbitCollection>("L1MuKBMTrack").setBranchAlias("L1MuKBMTrackOrbitCollection");
}

L1TMuonBarrelScoutingKalmanTrackProducer::~L1TMuonBarrelScoutingKalmanTrackProducer() {
  if (algo_ != nullptr)
    delete algo_;

  if (trackFinder_ != nullptr)
    delete trackFinder_;

  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void L1TMuonBarrelScoutingKalmanTrackProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  Handle<L1MuKBMTCombinedStubCollection> stubsCollection;
  Handle<MuonOrbitCollection> muonsCollection;
  iEvent.getByToken(src_, stubsCollection);
  if (matchGmt_) iEvent.getByToken(gmtSrc_, muonsCollection);

  std::vector<int> seenBxs;
  L1MuKBMTCombinedStubRefVector stubs;
  for (size_t i=0; i<stubsCollection->size(); ++i) {
    L1MuKBMTCombinedStubRef r(stubsCollection, i);
    stubs.push_back(r);
    seenBxs.push_back((*stubsCollection)[i].bxNum());
  }

  std::unique_ptr<L1MuKBMTrackOrbitCollection> kbmTrackCollection(new L1MuKBMTrackOrbitCollection);
  std::vector<std::vector<L1MuKBMTrack>> kbmTrackBuffer(3565);
  unsigned nKbmTrack = 0;

  std::sort(seenBxs.begin(), seenBxs.end());
  seenBxs.erase(std::unique(seenBxs.begin(), seenBxs.end()), seenBxs.end());

  double l1dr = 0.0, l1dr_min = 0.0;
  int l1_match_i = -1, l1_i = -1;
  // int64_t n_matches = 0;
  // int64_t n_gmt_m = 0;
  // double l1_reg_m_physPhi, l1_reg_m_physEta;
  // double l1_gmt_m_physPhi, l1_gmt_m_physEta;

  for (const auto& bx : seenBxs) {
    L1MuKBMTrackCollection tmp = trackFinder_->process(algo_, stubs, bx);
    if (tmp.size()==0) continue;

    if (matchGmt_) {
      const auto& muons = muonsCollection->bxIterator(bx);
      for (const auto& muon : muons) {
        const l1t::Muon gmt_m = getL1TMuon(muon);

        l1dr_min = 100.0;
        l1_match_i = -1;
        l1_i = -1;
        // l1_reg_m_physPhi = -100.0;
        // l1_gmt_m_physPhi = gmt_m.hwPhi()/phiMult_;
        // l1_reg_m_physEta = -100.0;
        // l1_gmt_m_physEta = gmt_m.hwEta()/etaMult_;        

        // barrel gmt muons
        if ((gmt_m.tfMuonIndex()>=36) && (gmt_m.tfMuonIndex()<=71)) {
        
          for (const auto& track : tmp) { // loop over KBMTF tracks in same BX
            const l1t::RegionalMuonCand bmtf_m = algo_->convertToBMTF(track);
            ++l1_i;
            l1dr = calcDr(bmtf_m, gmt_m);
            if (l1dr < l1dr_min) {
              l1_match_i = l1_i;
              l1dr_min = l1dr;
              // l1_reg_m_physPhi = calcGlobalPhi(bmtf_m)/phiMult_;
              // l1_reg_m_physEta = bmtf_m.hwEta()/etaMult_;
            }
          }
        }

        if (l1_match_i!=-1) {
          kbmTrackBuffer[bx].push_back(tmp[l1_match_i]);
          nKbmTrack++;
        }
      }
    } else {
      for (const auto& track : tmp) {
        kbmTrackBuffer[bx].push_back(track);
        nKbmTrack++;
      }
    }
    // out->push_back(bx, track);
    // algo_->addBMTFMuon(bx, track, outBMTF);
  }

  kbmTrackCollection->fillAndClear(kbmTrackBuffer, nKbmTrack);
  iEvent.put(std::move(kbmTrackCollection), "L1MuKBMTrack");
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void L1TMuonBarrelScoutingKalmanTrackProducer::beginStream(edm::StreamID) {}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void L1TMuonBarrelScoutingKalmanTrackProducer::endStream() {}

void L1TMuonBarrelScoutingKalmanTrackProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

unsigned L1TMuonBarrelScoutingKalmanTrackProducer::calcGlobalPhi(const l1t::RegionalMuonCand& l1_reg_m) {

  unsigned globalPhi = l1_reg_m.processor()*48 + l1_reg_m.hwPhi();
  globalPhi += 576 - 24;      // first processor starts at -15degrees in cms phi
  globalPhi = globalPhi%576;  // wrap around

  return globalPhi;
}



double L1TMuonBarrelScoutingKalmanTrackProducer::calcDr(const l1t::RegionalMuonCand& l1_reg_m, const l1t::Muon& l1_gmt_m) {

  double l1_reg_m_physPhi = calcGlobalPhi(l1_reg_m)/phiMult_;
  if (l1_reg_m_physPhi > M_PI) { l1_reg_m_physPhi -= 2*M_PI; }

  double l1_gmt_m_physPhi = l1_gmt_m.hwPhi()/phiMult_;
  if (l1_gmt_m_physPhi > M_PI) { l1_gmt_m_physPhi -= 2*M_PI; }

  double dphi = abs(l1_reg_m_physPhi - l1_gmt_m_physPhi);
  if (dphi > M_PI) { dphi = std::abs(dphi - 2*M_PI); }

  double deta = l1_reg_m.hwEta()/etaMult_ - l1_gmt_m.hwEta()/etaMult_;

  double dr = std::sqrt(dphi*dphi + deta*deta);

  if (debug_) {
    std::cout << "****** dr calc debug ******" << std::endl;
    std::cout << "l1_reg_m_phi = " << l1_reg_m_physPhi << std::endl;
    std::cout << "l1_gmt_m_phi = " << l1_gmt_m_physPhi << std::endl;
    std::cout << "l1_reg_m_eta = " << l1_reg_m.hwEta()/etaMult_ << std::endl;
    std::cout << "l1_gmt_m_eta = " << l1_gmt_m.hwEta()/etaMult_ << std::endl;
    std::cout << "dphi         = " << dphi << std::endl;
    std::cout << "deta         = " << deta << std::endl;
    std::cout << "dr           = " << dr   << std::endl;
  }

  return dr;
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1TMuonBarrelScoutingKalmanTrackProducer);
