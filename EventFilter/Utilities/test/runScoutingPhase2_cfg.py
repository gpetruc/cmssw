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


options.parseArguments()

cmsswbase = os.path.expandvars("$CMSSW_BASE/")

process = cms.Process("SCPU")
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(512)
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

fuDir = options.fuBaseDir+"/run"+str(options.runNumber).zfill(6)
buDir = options.buBaseDir+"/run"+str(options.runNumber).zfill(6)
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
        buDir + "/" + "run" + str(options.runNumber) + "_ls0126_index000000.raw"
    )
)
os.system("touch " + buDir + "/" + "fu.lock")
os.system("cp -v Puppi.fedRaw " + buDir + "/" + "run" + str(options.runNumber) + "_ls0126_index000000.raw")

## test pluging
process.scPhase2PuppiRawToDigi = cms.EDProducer('ScPhase2PuppiRawToDigi',
  src = cms.InputTag('rawDataCollector'),
  fedIDs = cms.vuint32(0),
  runSimpleUnpacker = cms.bool(True)
)

rawToDigiTask = cms.Task(
  process.scPhase2PuppiRawToDigi
)

process.p = cms.Path(process.scPhase2PuppiRawToDigi)