import os

import FWCore.ParameterSet.Config as cms
from L1TriggerScouting.Phase2.options_heterogeneous_w3pi_cff import args
from L1TriggerScouting.Phase2.modules import (
    l1sc_L1TScPhase2PuppiRawToDigi_alpaka,
    l1sc_L1TScPhase2W3Pi_alpaka,
    l1sc_L1TScPhase2AlpakaAnalyzer
)
  

args.parseArguments()
process = cms.Process("L1TScPhase2HeterogeneousW3Pi")

# summary
process.options.wantSummary = cms.untracked.bool(True)

# enable multithreading
process.options.numberOfThreads = args.numberOfThreads if args.numberOfThreads > 1 else 1 
process.options.numberOfStreams = args.numberOfStreams if args.numberOfStreams > 1 else 1 

# enable alpaka and GPU support
process.load("Configuration.StandardSequences.Accelerators_cff")

# logging
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

# timing service
process.FastMonitoringService = cms.Service("FastMonitoringService")
process.load( "HLTrigger.Timer.FastTimerService_cfi" )
process.FastTimerService.writeJSONSummary = cms.untracked.bool(False)
process.FastTimerService.jsonFileName = cms.untracked.string(f'resources.{os.uname()[1]}.{args.name}.json')

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
    streams = cms.vuint32(*list(range(sum(args.buNumStreams))) if args.streams == [] else args.streams),
    src = cms.InputTag('rawDataCollector'),
    environment = cms.untracked.int32(args.environment),
)
process.L1TScPhase2W3Pi = l1sc_L1TScPhase2W3Pi_alpaka(
    alpaka = cms.untracked.PSet(
        backend = cms.untracked.string(args.backend)
    ),
    src = 'L1TScPhase2PuppiRawToDigi',
    environment = cms.untracked.int32(args.environment),
    # control params
    pT_min = cms.double(args.ptMin),
    pT_int = cms.double(args.ptInt),
    pT_max = cms.double(args.ptMax),
    invariant_mass_lower_bound = cms.double(args.invariantMassLB),
    invariant_mass_upper_bound = cms.double(args.invariantMassUB),
    min_deltar_threshold = cms.double(args.minDeltaR),
    max_deltar_threshold = cms.double(args.maxDeltaR),
    max_isolation_threshold = cms.double(args.maxIso),
    ang_sep_lower_bound = cms.double(args.angSepLB),
    # fast mode
    fast_path = cms.bool(args.fastPath),
)
process.L1TScPhase2W3PiAnalyzer = l1sc_L1TScPhase2AlpakaAnalyzer(
    puppi = 'L1TScPhase2PuppiRawToDigi',
    bx_lookup = 'L1TScPhase2PuppiRawToDigi',
    selected_bxs = 'L1TScPhase2W3Pi',
    w3pi_table = 'L1TScPhase2W3Pi',
    environment = cms.untracked.int32(args.environment),
    fast_path = cms.bool(args.fastPath),
)

# schedule the modules
process.path = cms.Path(
    process.L1TScPhase2PuppiRawToDigi + 
    process.L1TScPhase2W3Pi + 
    process.L1TScPhase2W3PiAnalyzer
)