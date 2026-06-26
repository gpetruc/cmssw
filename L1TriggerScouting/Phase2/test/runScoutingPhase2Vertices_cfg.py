from __future__ import print_function
import FWCore.ParameterSet.Config as cms
import os

from L1TriggerScouting.Phase2.options_cff import options, VarParsing

options.parseArguments()

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

if options.buNumStreams == []:
    options.buNumStreams = [1] * len(options.buBaseDir)

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
    os.makedirs(d, exist_ok=True)

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

open(os.path.join(buDirs[0], "fu.lock"), "a").close()

process.load("L1TriggerScouting.Phase2.unpackers_cff")
process.load("Configuration.StandardSequences.Accelerators_cff")

## Configure unpackers
process.scPhase2VertexRawToDigiStruct = process.scPhase2VertexRawToDigiStruct.clone(
  fedIDs = [1],
  splitFactor = cms.uint32(1)
)
# Alpaka modules
from L1TriggerScouting.Phase2.modules import l1sc_L1TScPhase2VertexRawToDigi_alpaka

process.scPhase2VertexRawToDigiAlpaka = l1sc_L1TScPhase2VertexRawToDigi_alpaka(
    alpaka = cms.untracked.PSet( backend = cms.untracked.string(options.backend) ),
    streams = process.scPhase2VertexRawToDigiStruct.fedIDs,
    splitFactor = cms.uint32(len(pfStreamIDs) // options.timeslices),
    src = process.scPhase2VertexRawToDigiStruct.src,
    environment = cms.untracked.int32(options.environment),
)

process.p_vertexalpaka = cms.Path(
    process.scPhase2VertexRawToDigiAlpaka
)

sched = [ getattr(process, "p_" + options.run)]

process.schedule = cms.Schedule(*sched)
