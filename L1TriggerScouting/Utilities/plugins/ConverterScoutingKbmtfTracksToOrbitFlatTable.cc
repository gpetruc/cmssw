#include "FWCore/Framework/interface/MakerMacros.h"

#include <fstream>
#include <memory>
#include <string>
#include <cmath>
#include <cstdint>

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/MessageLogger/interface/MessageDrop.h"

#include "DataFormats/L1Scouting/interface/OrbitCollection.h"
#include "DataFormats/L1Scouting/interface/L1ScoutingBMTFStub.h"
#include "DataFormats/L1Scouting/interface/OrbitFlatTable.h"

#include "L1Trigger/L1TMuonBarrel/interface/L1TMuonBarrelKalmanAlgo.h"
#include "L1Trigger/L1TMuonBarrel/interface/L1TMuonBarrelKalmanTrackFinder.h"
#include "L1Trigger/L1TMuonBarrel/interface/L1TMuonBarrelKalmanStubProcessor.h"

#include "DataFormats/L1TMuon/interface/L1MuKBMTrack.h"
#include "DataFormats/L1TMuon/interface/L1MuKBMTCombinedStub.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCandFwd.h"

#include "CondFormats/L1TObjects/interface/L1TMuonGlobalParams.h"
#include "CondFormats/DataRecord/interface/L1TMuonGlobalParamsRcd.h"
#include "L1Trigger/L1TMuon/interface/L1TMuonGlobalParamsHelper.h"
#include "L1Trigger/L1TMuon/interface/MicroGMTLUTFactories.h"

#include "L1TriggerScouting/Utilities/interface/convertToL1TFormat.h"

using namespace l1ScoutingRun3;

class ConverterScoutingKbmtfTracksToOrbitFlatTable : public edm::stream::EDProducer<> {
public:
  // constructor and destructor
  explicit ConverterScoutingKbmtfTracksToOrbitFlatTable(const edm::ParameterSet&);
  ~ConverterScoutingKbmtfTracksToOrbitFlatTable() override;

  void produce(edm::Event&, const edm::EventSetup&) override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:

  void beginRun(edm::Run const&, edm::EventSetup const&) override;
  void endRun(edm::Run const&, edm::EventSetup const&) override;

  unsigned calcGlobalPhi(l1t::RegionalMuonCand&);

  // the tokens to access the data
  edm::EDGetTokenT<L1MuKBMTrackOrbitCollection> src_;

  std::string name_, doc_;

  L1TMuonBarrelKalmanAlgo* algo_;

  bool addStubs_;

  std::shared_ptr<l1t::MicroGMTExtrapolationLUT> m_BEtaExtrapolation_;
  std::shared_ptr<l1t::MicroGMTExtrapolationLUT> m_BPhiExtrapolation_;

  std::unique_ptr<L1TMuonGlobalParamsHelper> microGMTParamsHelper_;
  edm::ESGetToken<L1TMuonGlobalParams, L1TMuonGlobalParamsRcd> m_microGMTParamsToken_;

  int fwRev_;
};
// -----------------------------------------------------------------------------

// -------------------------------- constructor  -------------------------------

ConverterScoutingKbmtfTracksToOrbitFlatTable::ConverterScoutingKbmtfTracksToOrbitFlatTable(const edm::ParameterSet& iConfig)
    : src_(consumes<L1MuKBMTrackOrbitCollection>(iConfig.getParameter<edm::InputTag>("src"))),
      name_(iConfig.getParameter<std::string>("name")),
      doc_(iConfig.getParameter<std::string>("doc")),
      algo_(new L1TMuonBarrelKalmanAlgo(iConfig.getParameter<edm::ParameterSet>("algoSettings"))),
      addStubs_(iConfig.getParameter<bool>("addStubs")) {
  produces<OrbitFlatTable>();

  m_microGMTParamsToken_ = esConsumes<L1TMuonGlobalParams, L1TMuonGlobalParamsRcd, edm::Transition::BeginRun>();
  microGMTParamsHelper_ = std::make_unique<L1TMuonGlobalParamsHelper>();
  fwRev_ = 0x8010000;
}
// -----------------------------------------------------------------------------


ConverterScoutingKbmtfTracksToOrbitFlatTable::~ConverterScoutingKbmtfTracksToOrbitFlatTable() {
  if (algo_ != nullptr)
    delete algo_;

  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
}


// ----------------------- method called for each orbit  -----------------------
void ConverterScoutingKbmtfTracksToOrbitFlatTable::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<L1MuKBMTrackOrbitCollection> src;
  iEvent.getByToken(src_, src);
  auto out = std::make_unique<OrbitFlatTable>(src->bxOffsets(), name_);
  out->setDoc(doc_);

  int outputShiftPhi = 3;
  int outputShiftEta = 3;
  if (fwRev_ >= 0x4010000) {
    outputShiftPhi = 2;
    outputShiftEta = 0;
  }

  std::vector<float> pt(out->size());
  std::vector<float> eta(out->size());
  std::vector<float> phi(out->size());
  std::vector<int16_t> charge(out->size());
  std::vector<int16_t> quality(out->size());
  std::vector<int16_t> dxy(out->size());
  std::vector<int16_t> index(out->size());
  std::vector<float> ptUnconstrained(out->size());
  std::vector<float> etaAtVtx(out->size());
  std::vector<float> phiAtVtx(out->size());

  std::vector<int16_t> nStub(out->size());
  std::vector<std::vector<int16_t>> sStation(4, std::vector<int16_t>(out->size(), 0));
  std::vector<std::vector<int16_t>> sSector(4, std::vector<int16_t>(out->size(), 0));
  std::vector<std::vector<int16_t>> sWheel(4, std::vector<int16_t>(out->size(), 0));
  std::vector<std::vector<int16_t>> sHwQual(4, std::vector<int16_t>(out->size(), 0));
  std::vector<std::vector<int16_t>> sHwPhi(4, std::vector<int16_t>(out->size(), 0));
  std::vector<std::vector<int16_t>> sHwPhiB(4, std::vector<int16_t>(out->size(), 0));
  std::vector<std::vector<int16_t>> sHwEta1(4, std::vector<int16_t>(out->size(), 0));
  std::vector<std::vector<int16_t>> sHwQEta1(4, std::vector<int16_t>(out->size(), 0));
  std::vector<std::vector<int16_t>> sHwEta2(4, std::vector<int16_t>(out->size(), 0));
  std::vector<std::vector<int16_t>> sHwQEta2(4, std::vector<int16_t>(out->size(), 0));
  std::vector<std::vector<int16_t>> sTag(4, std::vector<int16_t>(out->size(), 0));

  unsigned int i = 0;
  for (const L1MuKBMTrack& track : *src) {

    l1t::RegionalMuonCand bmtf_m = algo_->convertToBMTF(track);
    pt[i] = ugmt::fPt(bmtf_m.hwPt());
    eta[i] = ugmt::fEta(bmtf_m.hwEta());
    phi[i] = ugmt::fPhi(calcGlobalPhi(bmtf_m));
    charge[i] = bmtf_m.hwSign()==1? -1 : 1;
    quality[i] = bmtf_m.hwQual();
    dxy[i] = track.dxy();
    index[i] = bmtf_m.processor(); // wrong for now
    ptUnconstrained[i] = ugmt::fPtUnconstrained(bmtf_m.hwPtUnconstrained());

    int ptRedInWidth = m_BPhiExtrapolation_->getPtRedInWidth();
    int ptMask = (1 << ptRedInWidth) - 1;
    int etaRedInWidth = m_BPhiExtrapolation_->getEtaRedInWidth();
    int redEtaShift = 8 - etaRedInWidth;

    // only use LSBs of pt:
    int ptRed = bmtf_m.hwPt() & ptMask;
    // here we drop the LSBs and mask the MSB
    int etaAbsRed = (std::abs(bmtf_m.hwEta()) >> redEtaShift) & ((1 << etaRedInWidth) - 1);
    int deltaPhi = 0;
    int deltaEta = 0;

    if (bmtf_m.hwPt() < (1 << ptRedInWidth)) {  // extrapolation only for "low" pT muons
      int sign = 1;
      if (bmtf_m.hwSign() == 1) {
        sign = -1;
      }
      deltaPhi = (m_BPhiExtrapolation_->lookup(etaAbsRed, ptRed) << outputShiftPhi) * sign;
      deltaEta = (m_BEtaExtrapolation_->lookup(etaAbsRed, ptRed) << outputShiftEta);
      if (bmtf_m.hwEta() > 0) {
        deltaEta *= -1;
      }
    }

    etaAtVtx[i] = ugmt::fEta(bmtf_m.hwEta() + deltaEta);
    phiAtVtx[i] = ugmt::fPhi(calcGlobalPhi(bmtf_m) + deltaPhi);

    if (addStubs_) {
      nStub[i] = track.stubs().size();
      unsigned j = 0;
      for (const auto& stub : track.stubs()) {
        sStation[j][i] = (*stub).stNum();
        sSector[j][i] = (*stub).scNum();
        sWheel[j][i] = (*stub).whNum();
        sHwQual[j][i] = (*stub).quality();
        sHwPhi[j][i] = (*stub).phi();
        sHwPhiB[j][i] = (*stub).phiB();
        sHwEta1[j][i] = (*stub).eta1();
        sHwQEta1[j][i] = (*stub).qeta1();
        sHwEta2[j][i] = (*stub).eta2();
        sHwQEta2[j][i] = (*stub).qeta2();
        sTag[j][i] = (*stub).tag();
        ++j;
      }
    }

    ++i;
  }

  out->addColumn<float>("pt", pt, "pt (physical units)");
  out->addColumn<float>("eta", eta, "eta at second muon station (physical units)");
  out->addColumn<float>("phi", phi, "phi at second muon station (physical units)");
  out->addColumn<int16_t>("hwCharge", charge, "hwCharge (hw units)");
  out->addColumn<int16_t>("hwQual", quality, "hwQual (hw units)");
  out->addColumn<int16_t>("hwDXY", dxy, "untruncated transverse impact parameter (hw units)");
  out->addColumn<int16_t>("processor", index, "processor ([0-11])");
  out->addColumn<float>("ptUnconstrained", ptUnconstrained, "pt without vertex constraint (physical units)");
  out->addColumn<float>("etaAtVtx", etaAtVtx, "eta re-extrapolated at vertex (physical units)");
  out->addColumn<float>("phiAtVtx", phiAtVtx, "phi re-extrapolated at vertex (physical units)");

  if (addStubs_) {
    out->addColumn<int16_t>("nStub", nStub, "number of stubs used to reconstruct KBMTF track");
    for (int i=0; i<4; ++i) {
  	  out->addColumn<int16_t>("s"+std::to_string(i+1)+"Station", sStation[i], "stub station");
  	  out->addColumn<int16_t>("s"+std::to_string(i+1)+"Sector", sSector[i], "stub sector");
  	  out->addColumn<int16_t>("s"+std::to_string(i+1)+"Wheel", sWheel[i], "stub wheel");
  	  out->addColumn<int16_t>("s"+std::to_string(i+1)+"HwQual", sHwQual[i], "stub quality (hw units)");
  	  out->addColumn<int16_t>("s"+std::to_string(i+1)+"HwPhi", sHwPhi[i], "stub local phi position (hw units)");
  	  out->addColumn<int16_t>("s"+std::to_string(i+1)+"HwPhiB", sHwPhiB[i], "stub phi bending (hw units)");
  	  out->addColumn<int16_t>("s"+std::to_string(i+1)+"HwEta1", sHwEta1[i], "eta of first stub in chamber (hw units)");
  	  out->addColumn<int16_t>("s"+std::to_string(i+1)+"HwQEta1", sHwQEta1[i], "eta quality of first stub in chamber (hw units)");
  	  out->addColumn<int16_t>("s"+std::to_string(i+1)+"HwEta2", sHwEta2[i], "eta of second stub in chamber (hw units)");
  	  out->addColumn<int16_t>("s"+std::to_string(i+1)+"HwQEta2", sHwQEta2[i], "eta quality of second stub in chamber (hw units)");
  	  out->addColumn<int16_t>("s"+std::to_string(i+1)+"Tag", sTag[i], "tag=0 is for second stub in chamber");
    }
  }

  iEvent.put(std::move(out));
}


// ------------ method called when starting to processes a run  ------------
void ConverterScoutingKbmtfTracksToOrbitFlatTable::beginRun(edm::Run const& run, edm::EventSetup const& iSetup) {

  /*
  edm::ESHandle<L1TMuonGlobalParams> microGMTParamsHandle = iSetup.getHandle(m_microGMTParamsToken);

  std::unique_ptr<L1TMuonGlobalParams_PUBLIC> microGMTParams(
      new L1TMuonGlobalParams_PUBLIC(cast_to_L1TMuonGlobalParams_PUBLIC(*microGMTParamsHandle.product())));
  if (microGMTParams->pnodes_.empty()) {
    edm::ESHandle<L1TMuonGlobalParams> o2oProtoHandle = iSetup.getHandle(m_o2oProtoToken);
    microGMTParamsHelper = std::make_unique<L1TMuonGlobalParamsHelper>(*o2oProtoHandle.product());
  } else
    microGMTParamsHelper =
        std::make_unique<L1TMuonGlobalParamsHelper>(cast_to_L1TMuonGlobalParams(*microGMTParams.get()));
  */

  edm::ESHandle<L1TMuonGlobalParams> microGMTParamsHandle_ = iSetup.getHandle(m_microGMTParamsToken_);
  microGMTParamsHelper_ = std::make_unique<L1TMuonGlobalParamsHelper>(*microGMTParamsHandle_.product());
  if (!microGMTParamsHelper_) {
    edm::LogError("L1TMicroGMTLUTDumper") << "Could not retrieve parameters from Event Setup" << std::endl;
  }
  
  // int fwRev_ = 0x8010000; // microGMTParamsHelper_->fwVersion();

  m_BEtaExtrapolation_ = l1t::MicroGMTExtrapolationLUTFactory::create(microGMTParamsHelper_->bEtaExtrapolationLUT(), l1t::MicroGMTConfiguration::ETA_OUT, fwRev_);
  m_BPhiExtrapolation_ = l1t::MicroGMTExtrapolationLUTFactory::create(microGMTParamsHelper_->bPhiExtrapolationLUT(), l1t::MicroGMTConfiguration::PHI_OUT, fwRev_);
}

// ------------ method called when ending to processes a run  ------------
void ConverterScoutingKbmtfTracksToOrbitFlatTable::endRun(edm::Run const&, edm::EventSetup const&) {}

void ConverterScoutingKbmtfTracksToOrbitFlatTable::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

unsigned ConverterScoutingKbmtfTracksToOrbitFlatTable::calcGlobalPhi(l1t::RegionalMuonCand& l1_reg_m) {

  unsigned globalPhi = l1_reg_m.processor()*48 + l1_reg_m.hwPhi();
  globalPhi += 576 - 24;      // first processor starts at -15degrees in cms phi
  globalPhi = globalPhi%576;  // wrap around

  return globalPhi;
}

DEFINE_FWK_MODULE(ConverterScoutingKbmtfTracksToOrbitFlatTable);
