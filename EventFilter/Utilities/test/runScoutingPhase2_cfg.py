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
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "BU base directory")

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

options.register ('puppiMode',
                  'simple', # default value
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "puppi mode to run (simple, sparse, struct, sparseStruct, soa)")

options.parseArguments()
if options.puppiMode not in ("simple", "sparse", "struct", "sparseStruct", "soa", "all", "fast"):
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
process.MessageLogger = cms.Service("MessageLogger",
    cout = cms.untracked.PSet(threshold = cms.untracked.string( "INFO" )),
    destinations = cms.untracked.vstring( 'cout' )
)

#process.FastMonitoringService = cms.Service("FastMonitoringService",
#    sleepTime = cms.untracked.int32(1)
#)
#
#process.Timing = cms.Service("Timing",
#  summaryOnly = cms.untracked.bool(False),
#  useJobReport = cms.untracked.bool(True)
#)
#
process.EvFDaqDirector = cms.Service("EvFDaqDirector",
    useFileBroker = cms.untracked.bool(False),
    fileBrokerHostFromCfg = cms.untracked.bool(True),
    fileBrokerHost = cms.untracked.string("htcp40.cern.ch"),
    runNumber = cms.untracked.uint32(options.runNumber),
    baseDir = cms.untracked.string(options.fuBaseDir),
    buBaseDir = cms.untracked.string(options.buBaseDir),
    buBaseDirsAll = cms.untracked.vstring(options.buBaseDir,),
    buBaseDirsNumStreams = cms.untracked.vint32(0),
    directorIsBU = cms.untracked.bool(False),
)

fuDir = options.fuBaseDir+("/run%06d" % options.runNumber)
buDir = options.buBaseDir+("/run%06d" % options.runNumber)
for d in fuDir, buDir, options.fuBaseDir, options.buBaseDir:
  if not os.path.isdir(d):
    os.makedirs(d)

process.source = cms.Source("DAQSource",
    testing = cms.untracked.bool(True),
    dataMode = cms.untracked.string(options.daqSourceMode),
    verifyChecksum = cms.untracked.bool(True),
    useL1EventID = cms.untracked.bool(False),
    eventChunkBlock = cms.untracked.uint32(2 * 1024),
    eventChunkSize = cms.untracked.uint32(3 * 1024),
    maxChunkSize = cms.untracked.uint32(10 * 1024),
    numBuffers = cms.untracked.uint32(2),
    maxBufferedFiles = cms.untracked.uint32(2),
    fileListMode = cms.untracked.bool(True),
    fileNames = cms.untracked.vstring(
        buDir + "/" + "run%06d_ls%04d_index%06d_ts00.raw" % (options.runNumber, options.lumiNumber, 0)
    )
)
os.system("touch " + buDir + "/" + "fu.lock")

## test pluging
process.scPhase2PuppiRawToDigi = cms.EDProducer('ScPhase2PuppiRawToDigi',
  src = cms.InputTag('rawDataCollector'),
  fedIDs = cms.vuint32(0),
  runSimpleUnpacker = cms.bool(True),
  runStructUnpacker = cms.bool(False),
  runSOAUnpacker = cms.bool(False),
  vectorizeSOAUnpacker = cms.bool(False),
  noWrite = cms.bool(False),
  sparseBXVector = cms.bool(False),
)
process.scPhase2PuppiRawToDigiSparse = process.scPhase2PuppiRawToDigi.clone(
    sparseBXVector = True,
)
process.scPhase2PuppiRawToDigiStruct = process.scPhase2PuppiRawToDigi.clone(
    runSimpleUnpacker = False,
    runStructUnpacker = True
)
process.scPhase2PuppiRawToDigiSparseStruct = process.scPhase2PuppiRawToDigiStruct.clone(
    sparseBXVector = True,
)
process.scPhase2PuppiRawToDigiSOA = process.scPhase2PuppiRawToDigi.clone(
    runSimpleUnpacker = False,
    runSOAUnpacker = True
)
process.scPhase2PuppiRawToDigiVSOA = process.scPhase2PuppiRawToDigi.clone(
    runSimpleUnpacker = False,
    runSOAUnpacker = True,
    vectorizeSOAUnpacker = True
)
process.w3piSimple = cms.EDProducer("ScPhase2PuppiW3PiDemo",
    src = cms.InputTag("scPhase2PuppiRawToDigiSparse"),
    runSimple = cms.bool(True),
    runStruct = cms.bool(False),
    runSOA = cms.bool(False),
    noWrite = cms.bool(False)
)
process.w3piStruct = cms.EDProducer("ScPhase2PuppiW3PiDemo",
    src = cms.InputTag("scPhase2PuppiRawToDigiSparseStruct"),
    runSimple = cms.bool(False),
    runStruct = cms.bool(True),
    runSOA = cms.bool(False),
    noWrite = cms.bool(False)
)
process.w3piSOA = cms.EDProducer("ScPhase2PuppiW3PiDemo",
    src = cms.InputTag("scPhase2PuppiRawToDigiSOA"),
    runSimple = cms.bool(False),
    runStruct = cms.bool(False),
    runSOA = cms.bool(True),
    noWrite = cms.bool(False)
)
#irawToDigiTask = cms.Task(
#  process.scPhase2PuppiRawToDigi,
#  process.scPhase2PuppiRawToDigiSparse,
#  process.scPhase2PuppiRawToDigiStruct,
#  process.scPhase2PuppiRawToDigiSparseStruct,
#  process.scPhase2PuppiRawToDigiSOA,
#  process.w3piSimple,
#  process.w3piStruct,
#  process.w3piSOA
#)

process.p_simple = cms.Path(
  process.scPhase2PuppiRawToDigi
)
process.p_sparse = cms.Path(
  process.scPhase2PuppiRawToDigiSparse
  +process.w3piSimple
)
process.p_struct = cms.Path(
  process.scPhase2PuppiRawToDigiStruct
)
process.p_sparseStruct = cms.Path(
  process.scPhase2PuppiRawToDigiSparseStruct
  +process.w3piStruct
)
process.p_soa = cms.Path(
  process.scPhase2PuppiRawToDigiSOA
  +process.w3piSOA
)
process.p_all = cms.Path(
  process.scPhase2PuppiRawToDigi+
  process.scPhase2PuppiRawToDigiSparse+
  process.scPhase2PuppiRawToDigiStruct+
  process.scPhase2PuppiRawToDigiSparseStruct+
  process.scPhase2PuppiRawToDigiSOA+
  process.scPhase2PuppiRawToDigiVSOA+
  process.w3piSimple+
  process.w3piStruct+
  process.w3piSOA
)
process.p_fast = cms.Path(
  process.scPhase2PuppiRawToDigiStruct+
  process.scPhase2PuppiRawToDigiSparseStruct+
  process.scPhase2PuppiRawToDigiSOA+
  process.scPhase2PuppiRawToDigiVSOA+
  process.w3piStruct+
  process.w3piSOA
)
process.schedule = cms.Schedule(getattr(process, "p_"+options.puppiMode))