from __future__ import print_function
import FWCore.ParameterSet.Config as cms
import os

from L1TriggerScouting.Phase2.options_cff import options, VarParsing
options.register ('njets',
                  16,
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.int,
                  'Number of jets to reconstruct with SCGreedy (iterative seeded cone)'
)
options.register ('minSeedPt',
                  0.0,
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.float,
                  'Minimum pt cut for seed candidates (used by SCNMSWeightedMultiIter / LinkTree)'
)
options.register ('jetR',
                  0.4,
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.float,
                  'Legacy single-radius parameter used by auto / SCGreedy fallback / demo modules'
)
options.register ('dumpClusters',
                  False,
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.bool,
                  'Dump clusters to options.outFile')


# Additional algorithms
# canonical names:
#   SCGreedy                = iterative greedy seeded cone
#   SCNMS                   = old non-iterative seeded cone with split radii:
#                             RSeed for seed finding + old centroid accumulation,
#                             RClu for final nearest-axis assignment
#   SCNMSWeighted           = same old non-iterative seeded cone logic, but final
#                             assignment uses old weighted metric
#   SCNMSWeightedMultiIter  = newer weighted NMS-style seeded cone with explicit
#                             centroid iterations and separate RCen
#   LinkTree                = link-based clustering
#
# backward-compatible aliases kept:
#   iterative               -> SCGreedy
#   seededCone              -> SCNMS
#   seededConeNMSWeighted   -> SCNMSWeighted
#   linkTree                -> LinkTree
options.register('jetAlgo',
                 'SCGreedy',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Jet algorithm: SCGreedy | SCNMS | SCNMSWeighted | SCNMSWeightedMultiIter | LinkTree')
options.register('RSeed',
                 0.3,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 'Seed-finding radius; for SCNMS/SCNMSWeighted also used for old centroid accumulation')
options.register('RCen',
                 0.4,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 'Centroid radius (used only by SCNMSWeightedMultiIter)')
options.register('RClu',
                 0.4,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 'Final particle-assignment radius for SCNMS / SCNMSWeighted / SCNMSWeightedMultiIter')
options.register('RLink',
                 0.3,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 'Link radius for LinkTree')
options.register('alphaSeed',
                 2.0,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 'alpha in score = dr2 / pt_seed^alpha (used by SCNMSWeightedMultiIter)')
options.register('nCentroidIters',
                 1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Number of centroid iterations for SCNMSWeightedMultiIter: 0=use seed axis, >0 refine axis')

options.parseArguments()

if options.buNumStreams == []:
    options.buNumStreams.append(1)

analyses = options.analyses if options.analyses else ["w3pi", "hphijpsi", "h2rho", "h2phi"]
print(f"Analyses set to {analyses}")

if options.run not in ("unpack", "ak4", "sc4", "unpackAlpaka", "clueAlpaka", "sc4Alpaka", "sc4AlpakaTaus"):
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

if "alpaka" in options.run.lower():
  process.load("Configuration.StandardSequences.Accelerators_cff")

## Configure unpackers
process.scPhase2PFRawToDigiStruct = process.scPhase2PuppiRawToDigiStruct.clone(
  fedIDs = [*pfStreamIDs],
  splitFactor = cms.uint32(len(pfStreamIDs) // options.timeslices)
)

process.goodOrbitsByNBX.nbxMin = 3564 * options.timeslices // options.tmuxPeriod
process.goodOrbitsByNBX.unpackers = [ "scPhase2PFRawToDigiStruct" ]

process.scPhase2AK4PFDemo = cms.EDProducer("ScPhase2PuppiAKJetsDemo",
  src = cms.InputTag("scPhase2PFRawToDigiStruct"),
  rParam = cms.double(options.jetR)
)
process.scPhase2SC4PFDemo = cms.EDProducer("ScPhase2PuppiSCJetsDemo",
  src = cms.InputTag("scPhase2PFRawToDigiStruct"),
  rParam = cms.double(options.jetR),
  nJets = cms.uint32(options.njets),
  minSeedPt = cms.double(options.minSeedPt)
)

# Alpaka modules
if "alpaka" in options.run.lower():
  from L1TriggerScouting.Phase2.modules import (
      l1sc_L1TScPhase2PuppiRawToDigi_alpaka,
      l1sc_L1TScPhase2SCJets_alpaka
  )
  from L1TriggerScouting.TauTagging.modules import (
      l1sc_CLUETaus_alpaka,
      l1sc_SoftTauIdML_alpaka,
  )
  process.scPhase2PFRawToDigiAlpaka = l1sc_L1TScPhase2PuppiRawToDigi_alpaka(
      alpaka = cms.untracked.PSet( backend = cms.untracked.string(options.backend) ),
      streams = process.scPhase2PFRawToDigiStruct.fedIDs,
      splitFactor = cms.uint32(len(pfStreamIDs) // options.timeslices),
      src = process.scPhase2PFRawToDigiStruct.src,
      environment = cms.untracked.int32(options.environment),
  )

  process.CLUETaus = l1sc_CLUETaus_alpaka(
      alpaka = cms.untracked.PSet( backend = cms.untracked.string(options.backend) ),
      src = 'scPhase2PFRawToDigiAlpaka',
      dc = cms.double(0.2),
      rhoc = cms.double(5.0),
      dm = cms.double(0.4),
      wrapCoords = cms.bool(False),
      run_scout = cms.bool(True),
  )

  sc4Alpaka_kwargs = {
      "alpaka": cms.untracked.PSet( backend = cms.untracked.string(options.backend) ),
      "src": cms.InputTag("scPhase2PFRawToDigiAlpaka"),
      "algo": cms.string(options.jetAlgo)
  }

  if options.jetAlgo == "SCGreedy":
      sc4Alpaka_kwargs.update({
          "rParam": cms.double(options.jetR),
          "nJets": cms.uint32(options.njets),
      })
  elif options.jetAlgo == "SCNMS" or options.jetAlgo == "SCNMSWeighted":
      sc4Alpaka_kwargs.update({
          "RSeed": cms.double(options.RSeed),
          "RClu": cms.double(options.RClu),
      })
  elif options.jetAlgo == "SCNMSWeightedMultiIter":
      sc4Alpaka_kwargs.update({
          "RSeed": cms.double(options.RSeed),
          "RCen": cms.double(options.RCen),
          "RClu": cms.double(options.RClu),
          "alphaSeed": cms.double(options.alphaSeed),
          "minSeedPt": cms.double(options.minSeedPt),
          "nCentroidIters": cms.uint32(options.nCentroidIters),
      })
  elif options.jetAlgo == "LinkTree":
      sc4Alpaka_kwargs.update({
          "RLink": cms.double(options.RLink),
          "minSeedPt": cms.double(options.minSeedPt),
      })

  process.scPhase2SC4PFAlpaka = l1sc_L1TScPhase2SCJets_alpaka(
      **sc4Alpaka_kwargs
  )

  process.SoftTauIdSC4 = l1sc_SoftTauIdML_alpaka(
      alpaka = cms.untracked.PSet( backend = cms.untracked.string(options.backend) ),
      pf = 'scPhase2PFRawToDigiAlpaka',
      clusters = 'scPhase2SC4PFAlpaka',
      model = cms.FileInPath("L1TriggerScouting/TauTagging/data/softtauid_sigmoid.pt"),
      maxBatchSize = cms.uint32(64),
  )

  process.goodOrbitsByNBX.unpackersAlpaka = [ "scPhase2PFRawToDigiAlpaka" ]
  process.goodOrbitsByNBX.unpackers = []

  process.p_unpackAlpaka = cms.Path(
    process.scPhase2PFRawToDigiAlpaka +
    process.goodOrbitsByNBX
  )
  process.p_clueAlpaka = cms.Path(
    process.scPhase2PFRawToDigiAlpaka +
    process.goodOrbitsByNBX +
    process.CLUETaus
  )
  process.p_sc4Alpaka = cms.Path(
    process.scPhase2PFRawToDigiAlpaka +
    process.goodOrbitsByNBX +
    process.scPhase2SC4PFAlpaka
  )
  process.p_sc4AlpakaTaus = cms.Path(
    process.scPhase2PFRawToDigiAlpaka +
    process.goodOrbitsByNBX +
    process.scPhase2SC4PFAlpaka +
    process.SoftTauIdSC4
  )


process.p_unpack = cms.Path(
  process.scPhase2PFRawToDigiStruct +
  process.goodOrbitsByNBX
)
process.p_ak4 = cms.Path(
  process.scPhase2PFRawToDigiStruct +
  process.goodOrbitsByNBX +
  process.scPhase2AK4PFDemo
)
process.p_sc4 = cms.Path(
  process.scPhase2PFRawToDigiStruct +
  process.goodOrbitsByNBX +
  process.scPhase2SC4PFDemo
)

if options.run not in ("both","inclusive","selected"):
  sched = [ getattr(process, "p_" + options.run)]
  if options.dumpClusters:
    process.scPhase2PFStructToTable = cms.EDProducer("ScPuppiToOrbitFlatTable",
       src = cms.InputTag("scPhase2PFRawToDigiStruct"),
       name = cms.string("L1PF"),
       doc = cms.string("L1PF candidates from Correlator Layer 1"),
    )
    process.p_pfTab = cms.Path(
       process.scPhase2PFRawToDigiStruct +
       process.scPhase2PFStructToTable
    )
    sched.append(process.p_pfTab)
    process.out = cms.OutputModule("OrbitNanoAODOutputModule",
     fileName = cms.untracked.string(options.outFile),
     SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring()),
     outputCommands = cms.untracked.vstring("drop *",
       "keep l1ScoutingRun3OrbitFlatTable_*_*_*")
    )
    process.p_out = cms.EndPath(process.out)

    if options.run in ("sc4Alpaka",):
      process.dumpClusters = cms.EDProducer("ClusterSoAToOrbitFlatTable",
          srcBx = cms.InputTag("scPhase2PFRawToDigiAlpaka"),
          srcClusters = cms.InputTag("scPhase2SC4PFAlpaka"),
          name = cms.string("SC4AlpakaClusters"),
          doc = cms.string("")
      )
      process.dumpJets = cms.EDProducer("ClusterObjSoAToOrbitFlatTable",
          srcBx = cms.InputTag("scPhase2SC4PFAlpaka"),
          srcClusters = cms.InputTag("scPhase2SC4PFAlpaka"),
          name = cms.string("SC4AlpakaJets"),
          doc = cms.string(""),
      )
      process.p_dump = cms.Path(process.dumpClusters + process.dumpJets)
      sched.append(process.p_dump)

    if options.run in ("clueAlpaka",):
      process.dumpClusters = cms.EDProducer("ClusterSoAToOrbitFlatTable",
          srcBx = cms.InputTag("scPhase2PFRawToDigiAlpaka"),
          srcClusters = cms.InputTag("CLUETaus"),
          name = cms.string("ClueClusters"),
          doc = cms.string("")
      )
      process.p_dump = cms.Path(process.dumpClusters)
      sched.append(process.p_dump)

    sched.append(process.p_out)
else:
  sched = [ process.p_inclusive, process.p_selected ]
  if options.run in ("inclusive", "selected"):
    sched = [ getattr(process, "p_" + options.run) ]
  if options.outMode != "none":
    sched.append(getattr(process, "o_"+options.outMode))

process.schedule = cms.Schedule(*sched)
