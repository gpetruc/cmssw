from __future__ import print_function
import FWCore.ParameterSet.Config as cms
import os

from L1TriggerScouting.Phase2.options_cff import options
options.parseArguments()
if options.buNumStreams == []:
    options.buNumStreams.append(2)
fullAnalysesList = ["w3pi", "wdsg", "wpig", "zdee",
                    "z2phiRecMeson", "z2rhoRecMeson",
                    "h2phiRecMeson", "h2rhoRecMeson", "hphijpsiRecMeson",
                    "hphigammaRecMeson", "hrhogammaRecMeson", "hjpsigammaRecMeson",
                    "hphijpsiEERecMeson", "hjpsigammaEERecMeson"]
analyses = options.analyses if options.analyses else fullAnalysesList
print(f"Analyses set to {analyses}")

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

if options.puppiStreamIDs == [] and options.tkEmStreamIDs ==  []:
    nStreamsTot = sum(options.buNumStreams)
    puppiStreamIDs = list(range(nStreamsTot//2)) # take first half
    tkEmStreamIDs = list(range(nStreamsTot//2, nStreamsTot)) # take second half
else:
    puppiStreamIDs = options.puppiStreamIDs
    tkEmStreamIDs = options.tkEmStreamIDs

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
process.FastTimerService.enableTimingPaths = cms.untracked.bool(True)
process.FastTimerService.enableTimingModules = cms.untracked.bool(True)
process.FastTimerService.useRealTimeClock = cms.untracked.bool(True)

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
    eventChunkBlock = cms.untracked.uint32(4 * 2 * 1024),
    eventChunkSize = cms.untracked.uint32(4 * 2 * 1024),
    maxChunkSize = cms.untracked.uint32(4 * 4 * 1024),
    numBuffers = cms.untracked.uint32(4),
    maxBufferedFiles = cms.untracked.uint32(4),
    fileListMode = cms.untracked.bool(options.broker == "none"),
    fileNames = cms.untracked.vstring(
        buDirs[0] + "/" + "run%06d_ls%04d_index%06d_stream00.raw" % (options.runNumber, options.lumiNumber, 1),
    )
)
os.system("touch " + buDirs[0] + "/" + "fu.lock")

# Declare analysis to run
process.load("L1TriggerScouting.Phase2.unpackers_cff")
# Declare candidates reconstruction
process.load("L1TriggerScouting.Phase2.candidateReco_cff")
# Declare analyses to run
process.load("L1TriggerScouting.Phase2.rareDecayAnalyses_cff")
process.load("L1TriggerScouting.Phase2.darkPhotonAnalyses_cff")
# Declare the filter of the data that keeps only data
# belonging to selected BXs
process.load("L1TriggerScouting.Phase2.maskedCollections_cff")
# Declare the flat table (ntuples) output
process.load("L1TriggerScouting.Phase2.nanoAODOutputs_cff")

## Configure unpackers
process.scPhase2PuppiRawToDigiStruct.fedIDs = [*puppiStreamIDs]
process.scPhase2TkEmRawToDigiStruct.fedIDs = [*tkEmStreamIDs]
process.goodOrbitsByNBX.nbxMin = 3564 * options.timeslices // options.tmuxPeriod
process.goodOrbitsByNBX.unpackers = [ "scPhase2PuppiRawToDigiStruct", "scPhase2TkEmRawToDigiStruct"]

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
  process.s_unpackers +
  process.prescaleInclusive
)
process.p_inclusive.associate(process.candRecoTasks)
process.p_inclusive.associate(cms.Task(
    process.scPhase2PuppiStructToTable,
    process.scPhase2TkEgTableProducersTask,
    process.scPhase2RecIsoTkEmStructToTable,
    process.scPhase2RecMesonStructToTable,
))

## Define selected processing (Physics streams)
process.p_selected = cms.Path(
  process.s_unpackers +
  process.s_analyses +
  process.scPhase2SelectedBXs +
  process.scPhase2PuppiMasked +
  process.scPhase2TkEmMasked +
  process.scPhase2TkEleMasked + 
  process.scPhase2RecIsoTkEmMasked +
  process.scPhase2RecMesonsMasked
)
process.p_selected.associate(process.candRecoTasks)
process.p_selected.associate(cms.Task(
    process.scPhase2PuppiMaskedStructToTable,
    process.scPhase2TkEgMaskedTableProducersTask,
    process.scPhase2RecIsoTkEmMaskedStructToTable,
    process.scPhase2RecMesonMaskedStructToTable
))

process.scPhase2NanoAll.fileName = options.outFile.replace(".root","")+".inclusive.root"
process.scPhase2NanoAll.SelectEvents.SelectEvents = ['p_inclusive']

process.scPhase2NanoSelected.fileName = options.outFile.replace(".root","")+".selected.root"
process.scPhase2NanoSelected.SelectEvents.SelectEvents = ['p_selected']
process.scPhase2NanoSelected.outputCommands += [ f"keep *_{a}Struct_*_*" for a in analyses ]

process.o_nanoInclusive = cms.EndPath(process.scPhase2NanoAll)
process.o_nanoSelected = cms.EndPath(process.scPhase2NanoSelected)
process.o_nanoBoth = cms.EndPath(process.scPhase2NanoAll + process.scPhase2NanoSelected)

sched = [ process.p_inclusive, process.p_selected ]
if options.run != "both":
    sched = [ getattr(process, "p_" + options.run)]

if options.outMode != "none":
  sched.append(getattr(process, "o_"+options.outMode))
process.schedule = cms.Schedule(*sched)
