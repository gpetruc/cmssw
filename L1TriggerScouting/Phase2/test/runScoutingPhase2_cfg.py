from __future__ import print_function
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import os

options = VarParsing.VarParsing ('analysis')
options.register ('runNumber',
                  37,
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.int,          # string, int, or float
                  "Run Number")

options.register ('lumiNumber',
                  1,
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.int,          # string, int, or float
                  "Run Number")

options.register ('daqSourceMode',
                  'ScoutingPhase2', # default value
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "DAQ source data mode")

options.register ('buBaseDir',
                  '/dev/shm/ramdisk', # default value
                  VarParsing.VarParsing.multiplicity.list,
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "BU base directory")

options.register ('buNumStreams',
                  1, # default value
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.int,          # string, int, or float
                  "Number of input streams (i.e. files) used simultaneously for each event")

options.register ('fuBaseDir',
                  '/dev/shm/data', # default value
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "BU base directory")

options.register ('fffBaseDir',
                  '/dev/shm', # default value
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "FFF base directory")

options.register ('numThreads',
                  1, # default value
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.int,          # string, int, or float
                  "Number of CMSSW threads")

options.register ('numFwkStreams',
                  1, # default value
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.int,          # string, int, or float
                  "Number of CMSSW streams")

options.register ('mode',
                  'simple', # default value
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "puppi mode to run (simple, struct, soa)")
                 
options.register ('outMode',
                  'none', # default value
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "puppi mode to run (none, struct, soa)")
                   
options.register ('outFile',
                  "NanoOutput.root",
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.string,
                  "Sub lumisection number to process")
options.register ('doTkTrackerMuons',
                  1, # default value
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.bool,          # string, int, or float
                  "Process tkTrackerMuon paths")


options.parseArguments()
if options.mode not in ("simple", "sparse", "struct", "sparseStruct", "soa", "all", "fast"):
    raise RuntimeError("Unsupported puppiMode %r" %options.puppiMode)

cmsswbase = os.path.expandvars("$CMSSW_BASE/")

process = cms.Process("SCPU")
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

process.options = cms.untracked.PSet(
    numberOfThreads = cms.untracked.uint32(options.numThreads),
    numberOfStreams = cms.untracked.uint32(options.numFwkStreams),
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(1),
    wantSummary = cms.untracked.bool(True)
)
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.EvFDaqDirector = cms.Service("EvFDaqDirector",
    useFileBroker = cms.untracked.bool(False),
    fileBrokerHostFromCfg = cms.untracked.bool(True),
    fileBrokerHost = cms.untracked.string("htcp40.cern.ch"),
    runNumber = cms.untracked.uint32(options.runNumber),
    baseDir = cms.untracked.string(options.fuBaseDir),
    buBaseDir = cms.untracked.string(options.buBaseDir[0]),
    buBaseDirsAll = cms.untracked.vstring(*options.buBaseDir),
    buBaseDirsNumStreams = cms.untracked.vint32([options.buNumStreams for dir in options.buBaseDir]),
    directorIsBU = cms.untracked.bool(False),
)

print(process.EvFDaqDirector.buBaseDirsNumStreams)
fuDir = options.fuBaseDir+("/run%06d" % options.runNumber)
buDirs = [b+("/run%06d" % options.runNumber) for b in options.buBaseDir]
for d in [fuDir, options.fuBaseDir] + buDirs + options.buBaseDir:
  if not os.path.isdir(d):
    os.makedirs(d)

process.source = cms.Source("DAQSource",
    testing = cms.untracked.bool(True),
    dataMode = cms.untracked.string(options.daqSourceMode),
    verifyChecksum = cms.untracked.bool(True),
    useL1EventID = cms.untracked.bool(False),
    eventChunkBlock = cms.untracked.uint32(2 * 1024),
    eventChunkSize = cms.untracked.uint32(2 * 1024),
    maxChunkSize = cms.untracked.uint32(4 * 1024),
    numBuffers = cms.untracked.uint32(4),
    maxBufferedFiles = cms.untracked.uint32(4),
    fileListMode = cms.untracked.bool(True),
    fileNames = cms.untracked.vstring(
        buDirs[0] + "/" + "run%06d_ls%04d_index%06d_stream00.raw" % (options.runNumber, options.lumiNumber, 1)
    )
)

os.system("touch " + buDirs[0] + "/" + "fu.lock")

## test pluging
scPhase2PuppiRawToDigi = cms.EDProducer('ScPhase2PuppiRawToDigi',
  src = cms.InputTag('rawDataCollector'),
  fedIDs = cms.vuint32(0),
  runCandidateUnpacker = cms.bool(False),
  runStructUnpacker = cms.bool(False),
  runSOAUnpacker = cms.bool(False),
)
process.scPhase2PuppiRawToDigiCandidate = scPhase2PuppiRawToDigi.clone(
    runCandidateUnpacker = True
)
process.scPhase2PuppiRawToDigiStruct = scPhase2PuppiRawToDigi.clone(
    runCandidateUnpacker = False,
    runStructUnpacker = True
)
process.scPhase2PuppiRawToDigiSOA = scPhase2PuppiRawToDigi.clone(
    runCandidateUnpacker = False,
    runSOAUnpacker = True
)
process.w3piCandidate = cms.EDProducer("ScPhase2PuppiW3PiDemo",
    src = cms.InputTag("scPhase2PuppiRawToDigiCandidate"),
    runCandidate = cms.bool(True),
    runStruct = cms.bool(False),
    runSOA = cms.bool(False)
)
process.w3piStruct = cms.EDProducer("ScPhase2PuppiW3PiDemo",
    src = cms.InputTag("scPhase2PuppiRawToDigiStruct"),
    runCandidate = cms.bool(False),
    runStruct = cms.bool(True),
    runSOA = cms.bool(False)
)
process.w3piSOA = cms.EDProducer("ScPhase2PuppiW3PiDemo",
    src = cms.InputTag("scPhase2PuppiRawToDigiSOA"),
    runCandidate = cms.bool(False),
    runStruct = cms.bool(False),
    runSOA = cms.bool(True)
)
process.scPhase2PuppiStructToTable = cms.EDProducer("ScPuppiToOrbitFlatTable",
    src = cms.InputTag("scPhase2PuppiRawToDigiStruct"),
    name = cms.string("L1Puppi"),
    doc = cms.string("L1Puppi candidates from Correlator Layer 2"),
)

scPhase2TrackerMuonRawToDigi = cms.EDProducer('ScPhase2TrackerMuonRawToDigi',
  src = cms.InputTag('rawDataCollector'),
  fedIDs = cms.vuint32(0),
  runCandidateUnpacker = cms.bool(False),
  runStructUnpacker = cms.bool(False),
  runSOAUnpacker = cms.bool(False),
)
process.scPhase2TrackerMuonRawToDigiStruct = scPhase2TrackerMuonRawToDigi.clone(
    fedIDs = cms.vuint32(1),
    runStructUnpacker = True
)
process.scPhase2TrackerMuonStructToTable = cms.EDProducer("ScTrackerMuonToOrbitFlatTable",
    src = cms.InputTag("scPhase2TrackerMuonRawToDigiStruct"),
    name = cms.string("L1TrackerMuon"),
    doc = cms.string("L1TrackerMuon candidates from Correlator Layer 2"),
)
process.dimuStruct = cms.EDProducer("ScPhase2TrackerMuonDiMuDemo",
    src = cms.InputTag("scPhase2TrackerMuonRawToDigiStruct"),
    runCandidate = cms.bool(False),
    runStruct = cms.bool(True),
    runSOA = cms.bool(False)
)

process.p_simple = cms.Path(
  process.scPhase2PuppiRawToDigiCandidate
  +process.w3piCandidate
)
process.p_struct = cms.Path(
  process.scPhase2PuppiRawToDigiStruct
  +process.scPhase2TrackerMuonRawToDigiStruct
  +process.w3piStruct
  +process.dimuStruct
  +process.scPhase2PuppiStructToTable
  +process.scPhase2TrackerMuonStructToTable
)
process.p_soa = cms.Path(
  process.scPhase2PuppiRawToDigiSOA
  +process.w3piSOA
)
process.p_all = cms.Path(
  process.scPhase2PuppiRawToDigiCandidate+
  process.scPhase2PuppiRawToDigiStruct+
  process.scPhase2TrackerMuonRawToDigiStruct+
  process.scPhase2PuppiRawToDigiSOA+
  process.scPhase2PuppiStructToTable+
  process.scPhase2TrackerMuonStructToTable+
  process.w3piCandidate+
  process.w3piStruct+
  process.w3piSOA+
  process.dimuStruct
)
process.p_fast = cms.Path(
  process.scPhase2PuppiRawToDigiStruct+
  process.scPhase2TrackerMuonRawToDigiStruct+
  process.scPhase2PuppiRawToDigiSOA+
  process.scPhase2PuppiStructToTable+
  process.scPhase2TrackerMuonStructToTable+
  process.w3piStruct+
  process.w3piSOA+
  process.dimuStruct

)
process.scPhase2PuppiStructNanoAll = cms.OutputModule("OrbitNanoAODOutputModule",
    fileName = cms.untracked.string(options.outFile),
    outputCommands = cms.untracked.vstring("drop *", "keep l1ScoutingRun3OrbitFlatTable_scPhase2PuppiStructToTable_*_*"),
    compressionLevel = cms.untracked.int32(4),
    compressionAlgorithm = cms.untracked.string("LZ4"),
)
process.scPhase2PuppiStructNanoW3pi = cms.OutputModule("OrbitNanoAODOutputModule",
    fileName = cms.untracked.string(options.outFile.replace(".root","")+".w3pi.root"),
    selectedBx = cms.InputTag("w3piStruct","selectedBx"),
    outputCommands = cms.untracked.vstring("drop *", 
        "keep l1ScoutingRun3OrbitFlatTable_scPhase2PuppiStructToTable_*_*",
        "keep l1ScoutingRun3OrbitFlatTable_w3piStruct_*_*"
        ),
    compressionLevel = cms.untracked.int32(4),
    compressionAlgorithm = cms.untracked.string("LZ4"),
)

process.scPhase2TrackerMuonStructNanoAll = cms.OutputModule("OrbitNanoAODOutputModule",
    fileName = cms.untracked.string(options.outFile.replace(".root","")+".tkMu.root"),
    outputCommands = cms.untracked.vstring("drop *", "keep l1ScoutingRun3OrbitFlatTable_scPhase2TrackerMuonStructToTable_*_*"),
    compressionLevel = cms.untracked.int32(4),
    compressionAlgorithm = cms.untracked.string("LZ4"),
)

process.scPhase2TrackerMuonStructNanoDiMu = cms.OutputModule("OrbitNanoAODOutputModule",
    fileName = cms.untracked.string(options.outFile.replace(".root","")+".dimu.root"),
    selectedBx = cms.InputTag("dimuStruct","selectedBx"),
    outputCommands = cms.untracked.vstring("drop *", 
        "keep l1ScoutingRun3OrbitFlatTable_scPhase2TrackerMuonStructToTable_*_*",
        "keep l1ScoutingRun3OrbitFlatTable_dimuStruct_*_*"
        ),
    compressionLevel = cms.untracked.int32(4),
    compressionAlgorithm = cms.untracked.string("LZ4"),
)

process.o_structAll  = cms.EndPath( process.scPhase2PuppiStructNanoAll + process.scPhase2TrackerMuonStructNanoAll )
process.o_structW3pi = cms.EndPath( process.scPhase2PuppiStructNanoW3pi )
process.o_structDiMu = cms.EndPath( process.scPhase2TrackerMuonStructNanoDiMu )

sched = [ getattr(process, "p_"+options.mode) ]
if options.outMode != "none":
    for mode in options.outMode.split(","):
        sched.append(getattr(process, "o_"+mode))

process.schedule = cms.Schedule(*sched)
