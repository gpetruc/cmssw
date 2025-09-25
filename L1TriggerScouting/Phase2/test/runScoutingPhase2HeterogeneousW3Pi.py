import os

import FWCore.ParameterSet.Config as cms
from L1TriggerScouting.Phase2.options_heterogeneous_w3pi_cff import args
from L1TriggerScouting.Phase2.modules import (
    l1sc_L1TScPhase2PuppiRawToDigi_alpaka,
    l1sc_L1TScPhase2W3Pi_alpaka,
    l1sc_L1TScPhase2W3PiAnalyzer
)
  

args.parseArguments()
process = cms.Process("L1TScPhase2HeterogeneousW3Pi")

# enable multithreading
process.options.numberOfThreads = args.numberOfThreads if args.numberOfThreads > 1 else 1 
process.options.numberOfStreams = args.numberOfStreams if args.numberOfStreams > 1 else 1 

# enable alpaka and GPU support
process.load("Configuration.StandardSequences.Accelerators_cff")

# logging
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MessageLogger.cerr.threshold = 'ERROR'

# process a limited number of events
process.maxEvents.input = args.numberOfEvents if args.numberOfEvents > 1 else 1 

if len(args.buNumStreams) != len(args.buBaseDir):
    raise RuntimeError("Mismatch between buNumStreams (%d) and buBaseDirs (%d)" % (len(args.buNumStreams), len(args.buBaseDir)))

process.EvFDaqDirector = cms.Service("EvFDaqDirector",
    useFileBroker = cms.untracked.bool(args.broker != "none"),
    fileBrokerHostFromCfg = cms.untracked.bool(False),
    fileBrokerHost = cms.untracked.string(args.broker.split(":")[0] if args.broker != "none" else "htcp40.cern.ch"),
    fileBrokerPort = cms.untracked.string(args.broker.split(":")[1] if args.broker != "none" else "8080"),
    runNumber = cms.untracked.uint32(args.runNumber),
    baseDir = cms.untracked.string(args.fuBaseDir),
    buBaseDir = cms.untracked.string(args.buBaseDir[0]),
    buBaseDirsAll = cms.untracked.vstring(*args.buBaseDir),
    buBaseDirsNumStreams = cms.untracked.vint32(*args.buNumStreams),
    directorIsBU = cms.untracked.bool(False),
)

fuDir = args.fuBaseDir+("/run%06d" % args.runNumber)
buDirs = [b+("/run%06d" % args.runNumber) for b in args.buBaseDir]
for d in [fuDir, args.fuBaseDir] + buDirs + args.buBaseDir:
  if not os.path.isdir(d):
    os.makedirs(d)

process.source = cms.Source("DAQSource",
    testing = cms.untracked.bool(True),
    dataMode = cms.untracked.string(args.daqSourceMode),
    verifyChecksum = cms.untracked.bool(True),
    useL1EventID = cms.untracked.bool(False),
    eventChunkBlock = cms.untracked.uint32(2 * 1024),
    eventChunkSize = cms.untracked.uint32(2 * 1024),
    maxChunkSize = cms.untracked.uint32(4 * 1024),
    numBuffers = cms.untracked.uint32(4),
    maxBufferedFiles = cms.untracked.uint32(4),
    fileListMode = cms.untracked.bool(args.broker == "none"),
    fileNames = cms.untracked.vstring(
        buDirs[0] + "/" + "run%06d_ls%04d_index%06d_stream00.raw" % (args.runNumber, args.lumiNumber, 1),
    )
)
os.system("touch " + buDirs[0] + "/" + "fu.lock")

# setup chain configs
process.L1TScPhase2PuppiRawToDigi = l1sc_L1TScPhase2PuppiRawToDigi_alpaka(
    alpaka = cms.untracked.PSet(
        backend = cms.untracked.string(args.backend)
    ),
    linksIds = cms.vuint32(*list(range(sum(args.buNumStreams))) if args.linksIds == [] else args.linksIds),
    src = cms.InputTag('rawDataCollector'),
    verbose = cms.untracked.bool(args.verbose),
)
process.L1TScPhase2W3Pi = l1sc_L1TScPhase2W3Pi_alpaka(
    alpaka = cms.untracked.PSet(
        backend = cms.untracked.string(args.backend)
    ),
    src = 'L1TScPhase2PuppiRawToDigi',
    verbose = cms.untracked.bool(args.verbose),
)
process.L1TScPhase2W3PiAnalyzer = l1sc_L1TScPhase2W3PiAnalyzer(
    puppi = 'L1TScPhase2W3Pi',
    nbx_map = 'L1TScPhase2W3Pi',
    table = 'L1TScPhase2W3Pi',
    bx_ct = 'L1TScPhase2PuppiRawToDigi',
    verbose = cms.untracked.bool(args.verbose),
    verboseLevel = cms.untracked.int32(args.verboseLevel)
)

# schedule the modules
process.path = cms.Path(
    process.L1TScPhase2PuppiRawToDigi + 
    process.L1TScPhase2W3Pi + 
    process.L1TScPhase2W3PiAnalyzer
)