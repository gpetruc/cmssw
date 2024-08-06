import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
from Configuration.Eras.Modifier_run3_2024_L1T_cff import run3_2024_L1T
import math

options = VarParsing.VarParsing ('analysis')

selbx = "DoubleMuPt0Qual8"

options.parseArguments()

process = cms.Process( "DUMP", run3_2024_L1T )


process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(-1)
)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.load('L1Trigger.L1TMuon.fakeGmtParams_cff')

process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(options.inputFiles)
)




bmtfKalmanTrackingSettings = cms.PSet(
  verbose = cms.bool(False),  #
  lutFile = cms.string("L1Trigger/L1TMuon/data/bmtf_luts/kalmanLUTs_v302.root"),
  initialK = cms.vdouble(-1.196,-1.581,-2.133,-2.263),
  initialK2 = cms.vdouble(-3.26e-4,-7.165e-4,2.305e-3,-5.63e-3),
#  eLoss = cms.vdouble(-2.85e-4,-6.21e-5,-1.26e-4,-1.23e-4),
  eLoss = cms.vdouble(+0.000765,0,0,0),

  aPhi = cms.vdouble(1.942, .01511, .01476, .009799),
  aPhiB = cms.vdouble(-1.508,-0.1237,-0.1496,-0.1333),
  aPhiBNLO = cms.vdouble(0.000331,0,0,0),

  bPhi = cms.vdouble(-1,.18245,.20898,.17286),
  bPhiB = cms.vdouble(-1,1.18245,1.20898,1.17286),
  phiAt2 = cms.double(0.15918),
  etaLUT0 = cms.vdouble(8.946,7.508,6.279,6.399),
  etaLUT1 = cms.vdouble(0.159,0.116,0.088,0.128),
  #generic cuts
  chiSquare = cms.vdouble(0.0,0.109375,0.234375,0.359375),
  chiSquareCutPattern = cms.vint32(7,11,13,14,15),
  chiSquareCutCurvMax = cms.vint32(2500,2500,2500,2500,2500),
  chiSquareCut = cms.vint32(126,126,126,126,126),


  #vertex cuts
  trackComp = cms.vdouble(1.75,1.25,0.625,0.250),
  trackCompErr1 = cms.vdouble(2.0,2.0,2.0,2.0),
  trackCompErr2 = cms.vdouble(0.218750,0.218750,0.218750,0.3125),
  trackCompCutPattern = cms.vint32(3,5,6,9,10,12),
  trackCompCutCurvMax = cms.vint32(34,34,34,34,34,34),   #this is shifted<<4
  trackCompCut        = cms.vint32(15,15,15,15,15,15),
  chiSquareCutTight   = cms.vint32(40,126,60,126,126,126),

  combos4=cms.vint32(9,10,11,12,13,14,15),
  combos3=cms.vint32(5,6,7),
  combos2=cms.vint32(3),
  combos1=cms.vint32(), #for future possible usage

  useOfflineAlgo = cms.bool(False),
  ###Only for the offline algo -not in firmware --------------------
  mScatteringPhi = cms.vdouble(2.49e-3,5.47e-5,3.49e-5,1.37e-5),
  mScatteringPhiB = cms.vdouble(7.22e-3,3.461e-3,4.447e-3,4.12e-3),
  pointResolutionPhi = cms.double(1.),
  pointResolutionPhiB = cms.double(500.),
  pointResolutionPhiBH = cms.vdouble(151., 173., 155., 153.),
  pointResolutionPhiBL = cms.vdouble(17866., 19306., 23984., 23746.),
  pointResolutionVertex = cms.double(1.),

  useNewQualityCalculation = cms.bool(False),
)

bmtfKalmanTrackingOfflineSettings = cms.PSet(
  verbose = cms.bool(False),  #
  lutFile = cms.string("L1Trigger/L1TMuon/data/bmtf_luts/kalmanLUTs_v302.root"),
  initialK = cms.vdouble(-1.196,-1.581,-2.133,-2.263),
  initialK2 = cms.vdouble(-3.26e-4,-7.165e-4,2.305e-3,-5.63e-3),
#  eLoss = cms.vdouble(-2.85e-4,-6.21e-5,-1.26e-4,-1.23e-4),
  eLoss = cms.vdouble(+0.000765,0,0,0),

  aPhi = cms.vdouble(1.942, .01511, .01476, .009799),
  aPhiB = cms.vdouble(-1.508,-0.1237,-0.1496,-0.1333),
  aPhiBNLO = cms.vdouble(0.000331,0,0,0),

  bPhi = cms.vdouble(-1,.18245,.20898,.17286),
  bPhiB = cms.vdouble(-1,1.18245,1.20898,1.17286),
  phiAt2 = cms.double(0.15918),
  etaLUT0 = cms.vdouble(8.946,7.508,6.279,6.399),
  etaLUT1 = cms.vdouble(0.159,0.116,0.088,0.128),
  #generic cuts
  chiSquare = cms.vdouble(0.0,0.109375,0.234375,0.359375),
  chiSquareCutPattern = cms.vint32(7,11,13,14,15),
  chiSquareCutCurvMax = cms.vint32(2500,2500,2500,2500,2500),
  chiSquareCut = cms.vint32(126,126,126,126,126),

  #vertex cuts
  trackComp = cms.vdouble(1.75,1.25,0.625,0.250),
  trackCompErr1 = cms.vdouble(2.0,2.0,2.0,2.0),
  trackCompErr2 = cms.vdouble(0.218750,0.218750,0.218750,0.3125),
  trackCompCutPattern = cms.vint32(3,5,6,9,10,12),
  trackCompCutCurvMax = cms.vint32(34,34,34,34,34,34),   #this is shifted<<4
  trackCompCut        = cms.vint32(15,15,15,15,15,15),
  chiSquareCutTight   = cms.vint32(40,126,60,126,126,126),

  combos4=cms.vint32(9,10,11,12,13,14,15),
  combos3=cms.vint32(5,6,7),
  combos2=cms.vint32(3),
  combos1=cms.vint32(), #for future possible usage

  useOfflineAlgo = cms.bool(True),
  ###Only for the offline algo -not in firmware --------------------
  mScatteringPhi = cms.vdouble(2.49e-3,5.47e-5,3.49e-5,1.37e-5),
  mScatteringPhiB = cms.vdouble(7.22e-3,3.461e-3,4.447e-3,4.12e-3),
  pointResolutionPhi = cms.double(1.),
  pointResolutionPhiB = cms.double(500.),
  pointResolutionPhiBH = cms.vdouble(151., 173., 155., 153.),
  pointResolutionPhiBL = cms.vdouble(17866., 19306., 23984., 23746.),
  pointResolutionVertex = cms.double(1.),

  useNewQualityCalculation = cms.bool(False),
)

run3_2024_L1T.toModify(
  bmtfKalmanTrackingSettings,
  useNewQualityCalculation = cms.bool(True),
)

run3_2024_L1T.toModify(
  bmtfKalmanTrackingOfflineSettings,
  useNewQualityCalculation = cms.bool(True),
)

process.kbmtfConvert = cms.EDProducer("convertToL1MuKBMTCombinedStub",
  src = cms.InputTag("FinalBxSelectorBMTFStub" if selbx else "l1ScBMTFUnpacker", "BMTFStub", "SCHLP"),
  bxMin = cms.int32(1),
  bxMax = cms.int32(3564),
  cotTheta_1 = cms.vint32(105,101,97,93,88,84,79,69,64,58,52,46,40,34,21,14,7,0,-7,-14,-21,-34,-40,-46,-52,-58,-64,-69,-79,-84,-88,-93,-97,-101,-105),
  cotTheta_2 = cms.vint32(93,89,85,81,77,73,68,60,55,50,45,40,34,29,17,12,6,0,-6,-12,-17,-29,-34,-40,-45,-50,-55,-60,-68,-73,-77,-81,-85,-89,-93),
  cotTheta_3 = cms.vint32(81,77,74,70,66,62,58,51,46,42,38,33,29,24,15,10,5,0,-5,-10,-15,-24,-29,-33,-38,-42,-46,-51,-58,-62,-66,-70,-74,-77,-81),
  debug = cms.bool(False)
)

process.kbmtfEmulation = cms.EDProducer("L1TMuonBarrelScoutingKalmanTrackProducer",
  src = cms.InputTag("kbmtfConvert"),
  algoSettings = bmtfKalmanTrackingSettings,
  trackFinderSettings = cms.PSet(
    sectorsToProcess = cms.vint32(0,1,2,3,4,5,6,7,8,9,10,11),
    verbose = cms.int32(0),
    sectorSettings = cms.PSet(
      # verbose = cms.int32(1),
      verbose = cms.int32(0),
      wheelsToProcess = cms.vint32(-2,-1,0,1,2),
      regionSettings = cms.PSet(
        verbose=cms.int32(0)
      )
    )
  ),
  gmtSrc = cms.InputTag("FinalBxSelectorMuon" if selbx else "l1ScGmtUnpacker", "Muon"),
  matchGmt = cms.bool(True),
  drCut = cms.double(0.1),
  phiMult = cms.double(576./(2*math.pi)),
  etaMult = cms.double(1./0.010875),
  debug = cms.bool(False)
)

process.kbmtfOfflineEmulation = cms.EDProducer("L1TMuonBarrelScoutingKalmanTrackProducer",
  src = cms.InputTag("kbmtfConvert"),
  algoSettings = bmtfKalmanTrackingOfflineSettings,
  trackFinderSettings = cms.PSet(
    sectorsToProcess = cms.vint32(0,1,2,3,4,5,6,7,8,9,10,11),
    verbose = cms.int32(0),
    sectorSettings = cms.PSet(
      # verbose = cms.int32(1),
      verbose = cms.int32(0),
      wheelsToProcess = cms.vint32(-2,-1,0,1,2),
      regionSettings = cms.PSet(
        verbose=cms.int32(0)
      )
    )
  ),
  gmtSrc = cms.InputTag("FinalBxSelectorMuon" if selbx else "l1ScGmtUnpacker", "Muon"),
  matchGmt = cms.bool(True),
  drCut = cms.double(0.1),
  phiMult = cms.double(576./(2*math.pi)),
  etaMult = cms.double(1./0.010875),
  debug = cms.bool(False)
)





process.scMuonTable = cms.EDProducer("ConvertScoutingMuonsToOrbitFlatTable",
  src = cms.InputTag("FinalBxSelectorMuon" if selbx else "l1ScGmtUnpacker", "Muon"),
  name = cms.string("L1Mu"),
  doc = cms.string("Muons from GMT"),
)
process.scKbmtfTable = cms.EDProducer("ConverterScoutingKbmtfTracksToOrbitFlatTable",
  src = cms.InputTag("kbmtfEmulation", "L1MuKBMTrack"),
  name = cms.string("L1KBMTF"),
  doc = cms.string("Re-emulated KBMTF muons"),
  algoSettings = bmtfKalmanTrackingSettings,
  addStubs = cms.bool(True)
)
process.scKbmtfOfflineTable = cms.EDProducer("ConverterScoutingKbmtfTracksToOrbitFlatTable",
  src = cms.InputTag("kbmtfOfflineEmulation", "L1MuKBMTrack"),
  name = cms.string("L1KBMTFOff"),
  doc = cms.string("Re-emulated KBMTF muons with offline propagation"),
  algoSettings = bmtfKalmanTrackingOfflineSettings,
  addStubs = cms.bool(True)
)

process.p = cms.Path(
  process.scMuonTable +
  # process.scKbmtfTable
  #process.scJetTable +
  #process.scEgammaTable +
  #process.scTauTable +
  # process.scStubsTable +
  process.kbmtfConvert +
  process.kbmtfEmulation +
  process.kbmtfOfflineEmulation +
  process.scKbmtfTable +
  process.scKbmtfOfflineTable
  #process.scSumTable
)

# process.p.replace(process.scKbmtfTable, process.kbmtfConvert + process.kbmtfEmulation + process.scKbmtfTable)

process.out = cms.OutputModule("OrbitNanoAODOutputModule",
    fileName = cms.untracked.string(options.outputFile),
    skipEmptyBXs = cms.bool(True),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('p')),
    outputCommands = cms.untracked.vstring("drop *", "keep l1ScoutingRun3OrbitFlatTable_*_*_*"),
    compressionLevel = cms.untracked.int32(5),
    compressionAlgorithm = cms.untracked.string("ZSTD"),
)

if selbx:
  process.out.outputCommands += [ "keep uints_*_SelBx_*" ]
  if selbx != "any":
    process.out.selectedBx = cms.InputTag(selbx, "SelBx")
    process.out.skipEmptyBXs = False
  else:
    process.out.selectedBx = cms.InputTag("FinalBxSelector", "SelBx")
    process.out.skipEmptyBXs = False

process.o = cms.EndPath(
  process.out
)
