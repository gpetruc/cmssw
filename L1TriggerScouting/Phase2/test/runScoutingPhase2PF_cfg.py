from __future__ import print_function
import FWCore.ParameterSet.Config as cms
import os

from L1TriggerScouting.Phase2.options_cff import options
options.parseArguments()
if options.buNumStreams == []:
    options.buNumStreams.append(1)
analyses = options.analyses if options.analyses else ["w3pi", "hphijpsi", "h2rho", "h2phi"]
print(f"Analyses set to {analyses}")

if options.run not in ("both", "inclusive", "selected", "candidate", "soa", "all", "fast", "alpaka", "unpack", "unpackAlpaka"):
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
process.MessageLogger.cerr.FwkReport.reportEvery = 10

if len(options.buNumStreams) != len(options.buBaseDir):
    raise RuntimeError("Mismatch between buNumStreams (%d) and buBaseDirs (%d)" % (len(options.buNumStreams), len(options.buBaseDir)))

if options.pfBarrelStreamIDs == [] and options.pfEndcapStreamIDs == []:
    pfStreamIDs = list(range(sum(options.buNumStreams))) # take all 
else:
    pfStreamIDs = options.pfBarrelStreamIDs + options.pfEndcapStreamIDs

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
if options.run in ("alpaka", "unpackAlpaka"): 
  process.load("Configuration.StandardSequences.Accelerators_cff")

## Configure unpackers
process.scPhase2PFRawToDigiStruct = process.scPhase2PuppiRawToDigiStruct.clone(
  fedIDs = [*pfStreamIDs],
  splitFactor = cms.uint32(len(pfStreamIDs) // options.timeslices)
)
process.goodOrbitsByNBX.nbxMin = 3564 * options.timeslices // options.tmuxPeriod
process.goodOrbitsByNBX.unpackers = [ "scPhase2PFRawToDigiStruct" ]

# Alpaka modules
if options.run in ("alpaka", "unpackAlpaka"): 
  from L1TriggerScouting.Phase2.modules import (
      l1sc_L1TScPhase2PuppiRawToDigi_alpaka
  )
  process.scPhase2PFRawToDigiAlpaka = l1sc_L1TScPhase2PuppiRawToDigi_alpaka(
      alpaka = cms.untracked.PSet( backend = cms.untracked.string(options.backend) ),
      linksIds = process.scPhase2PFRawToDigiStruct.fedIDs,
      src = process.scPhase2PFRawToDigiStruct.src,
      verbose = cms.untracked.bool(options.verbose),
  )
  process.goodOrbitsByNBX.unpackersAlpaka = [ "scPhase2PFRawToDigiAlpaka" ]
  if options.run in ("alpaka", "unpackAlpaka"):
    process.goodOrbitsByNBX.unpackers = []


  process.p_unpackAlpaka = cms.Path(
    process.scPhase2PFRawToDigiAlpaka +
    process.goodOrbitsByNBX
  )

process.p_unpack = cms.Path(
  process.scPhase2PFRawToDigiStruct +
  process.goodOrbitsByNBX
)

if options.run not in ("both","inclusive","selected"): 
  sched = [ getattr(process, "p_" + options.run)]
else:
  sched = [ process.p_inclusive, process.p_selected ]
  if options.run in ("inclusive", "selected"):
    sched = [ getattr(process, "p_" + options.run) ]
  if options.outMode != "none":
    sched.append(getattr(process, "o_"+options.outMode))

process.schedule = cms.Schedule(*sched)
