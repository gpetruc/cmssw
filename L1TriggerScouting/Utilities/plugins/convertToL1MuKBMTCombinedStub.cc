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

#include "DataFormats/L1Scouting/interface/OrbitCollection.h"
#include "DataFormats/L1Scouting/interface/L1ScoutingBMTFStub.h"

using namespace l1ScoutingRun3;

//
// class declaration
//

class convertToL1MuKBMTCombinedStub : public edm::stream::EDProducer<> {
public:
  explicit convertToL1MuKBMTCombinedStub(const edm::ParameterSet&);
  ~convertToL1MuKBMTCombinedStub() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginStream(edm::StreamID) override;
  void produce(edm::Event&, const edm::EventSetup&) override;
  void endStream() override;

  int calculateEta(uint i, int wheel, uint sector, uint station);
  L1MuKBMTCombinedStub buildStub(int wheel, int sector, int station,
                                 int phi, int phiB, bool tag,
                                 int eta, int qeta, int bx,
                                 int quality);
  L1MuKBMTCombinedStub buildStubNoEta(int wheel, int sector, int station,
                                      int phi, int phiB, bool tag,
                                      int bx, int quality);

  edm::EDGetTokenT<BMTFStubOrbitCollection> src_;
  int bxMin_;
  int bxMax_;
  std::vector<int> eta1_;
  std::vector<int> eta2_;
  std::vector<int> eta3_;
  const edm::EDPutTokenT<L1MuKBMTCombinedStubCollection> putToken_;
};
convertToL1MuKBMTCombinedStub::convertToL1MuKBMTCombinedStub(const edm::ParameterSet& iConfig)
    : src_(consumes<BMTFStubOrbitCollection>(iConfig.getParameter<edm::InputTag>("src"))),
      bxMin_(iConfig.getParameter<int>("bxMin")),
      bxMax_(iConfig.getParameter<int>("bxMax")),
      eta1_(iConfig.getParameter<std::vector<int>>("cotTheta_1")),
      eta2_(iConfig.getParameter<std::vector<int>>("cotTheta_2")),
      eta3_(iConfig.getParameter<std::vector<int>>("cotTheta_3")),
      putToken_(produces<L1MuKBMTCombinedStubCollection>()) {}

convertToL1MuKBMTCombinedStub::~convertToL1MuKBMTCombinedStub() {}

//
// member functions
//

// ------------ method called to produce the data  ------------
void convertToL1MuKBMTCombinedStub::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  Handle<BMTFStubOrbitCollection> stubsCollection;
  iEvent.getByToken(src_, stubsCollection);

  L1MuKBMTCombinedStubCollection kbmtStubs;
  for (int bx=bxMin_; bx<=bxMax_; ++bx) {
    const auto& l1dsStubs = stubsCollection->bxIterator(bx);
    for (const auto& l1dsStub : l1dsStubs) {
      if (l1dsStub.hwEta()==0) {
        kbmtStubs.push_back(buildStubNoEta(l1dsStub.wheel(), l1dsStub.sector(), l1dsStub.station(),
                                           l1dsStub.hwPhi(), l1dsStub.hwPhiB(), l1dsStub.tag(),
                                           bx, l1dsStub.hwQual()));
      } else {
        kbmtStubs.push_back(buildStub(l1dsStub.wheel(), l1dsStub.sector(), l1dsStub.station(),
                                      l1dsStub.hwPhi(), l1dsStub.hwPhiB(), l1dsStub.tag(),
                                      l1dsStub.hwEta(), l1dsStub.hwQEta(),
                                      bx, l1dsStub.hwQual()));
      }
    }
  }

  iEvent.emplace(putToken_, std::move(kbmtStubs));
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void convertToL1MuKBMTCombinedStub::beginStream(edm::StreamID) {}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void convertToL1MuKBMTCombinedStub::endStream() {}

void convertToL1MuKBMTCombinedStub::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

// from Kalman stub producer
int convertToL1MuKBMTCombinedStub::calculateEta(uint i, int wheel, uint sector, uint station) {
  int eta = 0;
  if (wheel > 0) {
    eta = 7 * wheel + 3 - i;
  } else if (wheel < 0) {
    eta = 7 * wheel + i - 3;
  } else {
    if (sector == 0 || sector == 3 || sector == 4 || sector == 7 || sector == 8 || sector == 11)
      eta = i - 3;
    else
      eta = 3 - i;
  }

  if (station == 1)
    eta = -eta1_[eta + 17];
  else if (station == 2)
    eta = -eta2_[eta + 17];
  else
    eta = -eta3_[eta + 17];

  return eta;
}

L1MuKBMTCombinedStub convertToL1MuKBMTCombinedStub::buildStub(int wheel, int sector, int station,
                                                              int phi, int phiB, bool tag,
                                                              int eta, int qeta, int bx,
                                                              int quality) {
  // convert eta hw values to global units
  int qeta1 = 0;
  int qeta2 = 0;
  int eta1 = 255;
  int eta2 = 255;
  int mask = 0;

  bool hasEta = false;
  for (uint i = 0; i < 7; ++i) {
    mask = (1<<i);
    if ((eta & mask) == 0)
      continue;
    if (!hasEta) {
      hasEta = true;
      eta1 = calculateEta(i, wheel, sector, station);
      if ((qeta & mask) == 1)
        qeta1 = 2;
      else
        qeta1 = 1;
    } else {
      eta2 = calculateEta(i, wheel, sector, station);
      if ((qeta & mask) == 1)
        qeta2 = 2;
      else
        qeta2 = 1;
    }
  }

  // TODO: tag = true is hardcode for now
  L1MuKBMTCombinedStub stub(wheel, sector, station, phi, phiB, tag, bx, quality, eta1, eta2, qeta1, qeta2);

  return stub;
}

L1MuKBMTCombinedStub convertToL1MuKBMTCombinedStub::buildStubNoEta(int wheel, int sector, int station,
                                                                   int phi, int phiB, bool tag,
                                                                   int bx, int quality) {

  int qeta1 = 0;
  int qeta2 = 0;
  int eta1 = 7;
  int eta2 = 7;
  L1MuKBMTCombinedStub stub(wheel, sector, station, phi, phiB, tag, bx, quality, eta1, eta2, qeta1, qeta2);

  return stub;
}

//define this as a plug-in
DEFINE_FWK_MODULE(convertToL1MuKBMTCombinedStub);
