from __future__ import print_function
import FWCore.ParameterSet.Config as cms
import os

from L1TriggerScouting.Phase2.options_cff import options
options.parseArguments()
if options.buNumStreams == []:
    options.buNumStreams.append(1)
analyses = options.analyses if options.analyses else ["w3pi", "hphijpsi", "h2rho", "h2phi"]
print(f"Analyses set to {analyses}")

if options.run not in ("both", "inclusive", "selected", "candidate", "all", "fast", "alpaka", "unpack", "unpackAlpaka"):
    raise RuntimeError("Unsupported run mode %r" % options.run)

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

if len(options.buNumStreams) != len(options.buBaseDir):
    raise RuntimeError("Mismatch between buNumStreams (%d) and buBaseDirs (%d)" % (len(options.buNumStreams), len(options.buBaseDir)))

if options.puppiStreamIDs == []:
    puppiStreamIDs = list(range(sum(options.buNumStreams))) # take all 
else:
    puppiStreamIDs = options.puppiStreamIDs

process.EvFDaqDirector = cms.Service("EvFDaqDirector",
    useFileBroker = cms.untracked.bool(options.broker != "none"),
    fileBrokerHostFromCfg = cms.untracked.bool(False),
    fileBrokerHost = cms.untracked.string(options.broker.split(":")[0] if options.broker != "none" else "htcp40.cern.ch"),
    fileBrokerPort = cms.untracked.string(options.broker.split(":")[1] if options.broker != "none" else "8080"),
    runNumber = cms.untracked.uint32(options.runNumber),
    baseDir = cms.untracked.string(options.fuBaseDir),
    buBaseDir = cms.untracked.string(options.buBaseDir[0]),
    buBaseDirsAll = cms.untracked.vstring(*options.buBaseDir),
    buBaseDirsNumStreams = cms.untracked.vint32(*options.buNumStreams),
    directorIsBU = cms.untracked.bool(False),
)
process.FastMonitoringService = cms.Service("FastMonitoringService")

process.load( "HLTrigger.Timer.FastTimerService_cfi" )
process.FastTimerService.writeJSONSummary = cms.untracked.bool(True)
process.FastTimerService.jsonFileName = cms.untracked.string(f'resources.{os.uname()[1]}.{options.task}.json')
#process.MessageLogger.cerr.FastReport = cms.untracked.PSet( limit = cms.untracked.int32( 10000000 ) )

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
    fileListMode = cms.untracked.bool(options.broker == "none"),
    fileNames = cms.untracked.vstring(
        buDirs[0] + "/" + "run%06d_ls%04d_index%06d_stream00.raw" % (options.runNumber, options.lumiNumber, 1),
    )
)
os.system("touch " + buDirs[0] + "/" + "fu.lock")

process.load("L1TriggerScouting.Phase2.unpackers_cff")
process.load("L1TriggerScouting.Phase2.rareDecayAnalyses_cff")
process.load("L1TriggerScouting.Phase2.maskedCollections_cff")
process.load("L1TriggerScouting.Phase2.nanoAODOutputs_cff")
if options.run in ("all", "fast", "alpaka", "unpackAlpaka"): 
  process.load("Configuration.StandardSequences.Accelerators_cff")

## Configure unpackers
process.scPhase2PuppiRawToDigiStruct.fedIDs = [*puppiStreamIDs]
process.goodOrbitsByNBX.nbxMin = 3564 * options.timeslices // options.tmuxPeriod
process.goodOrbitsByNBX.unpackers = [ "scPhase2PuppiRawToDigiStruct" ]


## Configure analyses
analysisModules = [getattr(process,f"{a}Struct") for a in analyses]
process.s_analyses = cms.Sequence(sum(analysisModules[1:], analysisModules[0]))


## Configure selected outputs
process.scPhase2SelectedBXs.analysisLabels = [cms.InputTag(f"{a}Struct", "selectedBx") for a in analyses]
process.scPhase2SelectedBXs.nPrint = cms.untracked.uint32(1000)

## Define inclusive processing (ZeroBias)
from FWCore.Modules.preScaler_cfi import preScaler
process.prescaleInclusive = preScaler.clone(prescaleFactor = options.prescaleInclusive)
process.p_inclusive = cms.Path(
  process.scPhase2PuppiRawToDigiStruct +
  process.goodOrbitsByNBX +
  process.prescaleInclusive +
  process.scPhase2PuppiStructToTable
)

## Define selected processing (Physics streams)
process.p_selected = cms.Path(
  process.scPhase2PuppiRawToDigiStruct +
  process.goodOrbitsByNBX +
  process.s_analyses +
  process.scPhase2SelectedBXs +
  process.scPhase2PuppiMasked + 
  process.scPhase2PuppiMaskedStructToTable
)

# Alpaka modules
if options.run in ("all","fast","alpaka", "unpackAlpaka"): 
  from L1TriggerScouting.Phase2.modules import (
      l1sc_L1TScPhase2PuppiRawToDigi_alpaka,
      l1sc_L1TScPhase2W3Pi_alpaka,
      l1sc_L1TScPhase2W3PiAnalyzer
  )
  process.scPhase2PuppiRawToDigiAlpaka = l1sc_L1TScPhase2PuppiRawToDigi_alpaka(
      alpaka = cms.untracked.PSet( backend = cms.untracked.string(options.backend) ),
      linksIds = process.scPhase2PuppiRawToDigiStruct.fedIDs,
      src = process.scPhase2PuppiRawToDigiStruct.src,
      verbose = cms.untracked.bool(options.verbose),
  )
  process.w3piAlpaka = l1sc_L1TScPhase2W3Pi_alpaka(
      alpaka = cms.untracked.PSet( backend = cms.untracked.string(options.backend) ),
      src = 'scPhase2PuppiRawToDigiAlpaka',
      verbose = cms.untracked.bool(options.verbose),
  )
  process.w3piAlpakaAnalyzer = l1sc_L1TScPhase2W3PiAnalyzer(
    puppi = 'w3piAlpaka',
    nbx_map = 'w3piAlpaka',
    table = 'w3piAlpaka',
    bx_ct = 'scPhase2PuppiRawToDigiAlpaka:nbx',
    verbose = cms.untracked.bool(options.verbose),
    verboseLevel = cms.untracked.int32(options.verboseLevel)
  )
  process.goodOrbitsByNBX.unpackersAlpaka = [ "scPhase2PuppiRawToDigiAlpaka" ]
  if options.run in ("alpaka", "unpackAlpaka"):
    process.goodOrbitsByNBX.unpackers = []

# Additional modules and paths for benchmarking different data structures
process.scPhase2PuppiRawToDigiCandidate = process.scPhase2PuppiRawToDigiStruct.clone(
    runStructUnpacker = cms.bool(False),
    runCandidateUnpacker = cms.bool(True),
)

process.w3piCandidate = process.w3piStruct.clone(
    src = 'scPhase2PuppiRawToDigiCandidate',
    runStruct = cms.bool(False),
    runCandidate = cms.bool(True),
)

process.p_candidate = cms.Path(
  process.scPhase2PuppiRawToDigiCandidate +
  process.w3piCandidate
)

if options.run in ("all", "fast", "alpaka", "unpackAlpaka"):
  process.p_all = cms.Path(
    process.scPhase2PuppiRawToDigiCandidate +
    process.scPhase2PuppiRawToDigiStruct +
    process.scPhase2PuppiRawToDigiAlpaka +
    process.goodOrbitsByNBX +
    process.w3piCandidate +
    process.w3piStruct +
    process.w3piAlpaka
  )

  process.p_fast = cms.Path(
    process.scPhase2PuppiRawToDigiStruct +
    process.scPhase2PuppiRawToDigiAlpaka +
    process.goodOrbitsByNBX +
    process.w3piStruct +
    process.w3piAlpaka
  )

  process.p_alpaka = cms.Path(
    process.scPhase2PuppiRawToDigiAlpaka +
    process.goodOrbitsByNBX +
    process.w3piAlpaka
  )

  process.p_unpackAlpaka = cms.Path(
    process.scPhase2PuppiRawToDigiAlpaka +
    process.goodOrbitsByNBX
  )

process.p_unpack = cms.Path(
  process.scPhase2PuppiRawToDigiStruct +
  process.goodOrbitsByNBX
)

process.scPhase2NanoAll.fileName = options.outFile.replace(".root","")+".inclusive.root"
process.scPhase2NanoAll.SelectEvents.SelectEvents = ['p_inclusive']
 
process.scPhase2NanoSelected.fileName = options.outFile.replace(".root","")+".selected.root"
process.scPhase2NanoSelected.SelectEvents.SelectEvents = ['p_selected']
process.scPhase2NanoSelected.outputCommands += [ f"keep *_{a}Struct_*_*" for a in analyses ]

process.o_nanoInclusive = cms.EndPath(process.scPhase2NanoAll)
process.o_nanoSelected = cms.EndPath(process.scPhase2NanoSelected)
process.o_nanoBoth = cms.EndPath(process.scPhase2NanoAll + process.scPhase2NanoSelected)

if options.run not in ("both","inclusive","selected"): 
  sched = [ getattr(process, "p_" + options.run)]
else:
  sched = [ process.p_inclusive, process.p_selected ]
  if options.run in ("inclusive", "selected"):
    sched = [ getattr(process, "p_" + options.run) ]
  if options.outMode != "none":
    sched.append(getattr(process, "o_"+options.outMode))

process.schedule = cms.Schedule(*sched)
