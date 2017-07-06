import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe2")

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.suppressError = cms.untracked.vstring("patTriggerFull")
process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(
        #'file:debug_Zmm_lostHits.root',
	'/store/express/Run2017B/ExpressPhysics/FEVT/Express-v1/000/297/218/00000/1CDCFB44-DD55-E711-8D01-02163E01450A.root',
	'/store/express/Run2017B/ExpressPhysics/FEVT/Express-v1/000/297/218/00000/242FEC1C-D855-E711-A559-02163E01373C.root',
	'/store/express/Run2017B/ExpressPhysics/FEVT/Express-v1/000/297/218/00000/3410C6CA-DC55-E711-A2F9-02163E012B04.root',
	'/store/express/Run2017B/ExpressPhysics/FEVT/Express-v1/000/297/218/00000/3CEFF073-DB55-E711-9DE5-02163E011CB0.root',
	'/store/express/Run2017B/ExpressPhysics/FEVT/Express-v1/000/297/218/00000/6686B87F-D355-E711-80A0-02163E013806.root',
	'/store/express/Run2017B/ExpressPhysics/FEVT/Express-v1/000/297/218/00000/6EA95BA9-E055-E711-AAC9-02163E012920.root',
	'/store/express/Run2017B/ExpressPhysics/FEVT/Express-v1/000/297/218/00000/96C34217-DB55-E711-A0F9-02163E013937.root',
    )
)
#JSON = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/DCSOnly/json_DCSONLY.txt'
#import FWCore.PythonUtilities.LumiList as LumiList
#process.source.lumisToProcess = LumiList.LumiList(filename = JSON).getVLuminosityBlockRange()

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )    

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '92X_dataRun2_Express_v2', '')

#process.load("HLTrigger.HLTfilters.triggerResultsFilter_cfi")
process.triggerResultsFilter = cms.EDFilter("TriggerResultsFilter",
    daqPartitions = cms.uint32(1),
    hltResults = cms.InputTag("TriggerResults","","HLT"),
    l1tIgnoreMask = cms.bool(False),
    l1tResults = cms.InputTag(""),
    l1techIgnorePrescales = cms.bool(False),
    throw = cms.bool(True),
    triggerConditions = cms.vstring('HLT_IsoMu20_v*','HLT_IsoMu24_v*')
)

process.tagMuons = cms.EDFilter("MuonSelector",
    src = cms.InputTag("muons"),
    cut = cms.string("pt > 20 && numberOfMatchedStations >= 2"+
                     " && pfIsolationR04().sumChargedHadronPt/pt < 0.2"),
)
process.oneTag  = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("tagMuons"), minNumber = cms.uint32(1))

process.probeMuons = cms.EDFilter("MuonSelector",
    src = cms.InputTag("muons"),
    cut = cms.string("pt > 10 && numberOfMatchedStations >= 1 && innerTrack.isNonnull"), 
)
process.tpPairs = cms.EDProducer("CandViewShallowCloneCombiner",
    cut = cms.string('60 < mass < 140'),
    decay = cms.string('tagMuons@+ probeMuons@-')
)
process.onePair = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("tpPairs"), minNumber = cms.uint32(1))

import TrackingTools.KalmanUpdators.Chi2MeasurementEstimator_cfi 
process.HitCollectorForDebug = TrackingTools.KalmanUpdators.Chi2MeasurementEstimator_cfi.Chi2MeasurementEstimator.clone(
    ComponentName = cms.string('HitCollectorForDebug'),
    MaxChi2 = cms.double(30.0), ## was 30 ## TO BE TUNED
    nSigma  = cms.double(3.),   ## was 3  ## TO BE TUNED 
)
process.clusterInfo = cms.EDAnalyzer("DebugPixelHits",
        pairs = cms.InputTag("tpPairs"),
        tracker = cms.InputTag("MeasurementTrackerEvent"),
        vertices = cms.InputTag("offlinePrimaryVertices"),
        lumiScalers = cms.InputTag("scalersRawToDigi"),
        # configuraton for refitter
        DoPredictionsOnly = cms.bool(False),
        Fitter = cms.string('KFFitterForRefitInsideOut'),
        TrackerRecHitBuilder = cms.string('WithAngleAndTemplate'),
        Smoother = cms.string('KFSmootherForRefitInsideOut'),
        MuonRecHitBuilder = cms.string('MuonRecHitBuilder'),
        RefitDirection = cms.string('oppositeToMomentum'),
        RefitRPCHits = cms.bool(True),
        Propagator = cms.string('SmartPropagatorAnyRKOpposite'),
        #Propagators
        PropagatorOpposite = cms.string("RungeKuttaTrackerPropagatorOpposite"),
        Chi2MeasurementEstimator = cms.string("HitCollectorForDebug"),
        #Error rescaling
        rescaleError = cms.double(1),
        #SiPixelQuality = cms.string(''),
        badComponentsFile = cms.string('/afs/cern.ch/work/g/gpetrucc/Tracking/CMSSW_9_2_3_patch2/src/RecoTracker/DebugTools/test/badComponents.txt'),
        ## https://github.com/cms-sw/cmssw/blob/9b7f92a91b55fe1bf3e38435a6afd5b97dea4c9f/RecoLocalTracker/SubCollectionProducers/src/JetCoreClusterSplitter.cc#L139-L153
        debug = cms.untracked.int32(100)
)

process.tagAndProbe = cms.Path( 
    #process.triggerResultsFilter +
    process.tagMuons +
    process.oneTag     +
    process.probeMuons +
    process.tpPairs    +
    process.onePair    +
    process.MeasurementTrackerEvent +
    process.clusterInfo 
)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("debug_Zmm_lostHits.root"),
    outputCommands = cms.untracked.vstring("keep *", "drop *_*_*_TagProbe"),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("tagAndProbe")),
)
#uncomment to skim the selected events
#process.end = cms.EndPath(process.out)

# this below probably not needed
process.TFileService = cms.Service("TFileService", fileName = cms.string("debugHits.root"))

if False:
    del process.clusterInfo.badComponentsFile
    process.TFileService = cms.Service("TFileService", fileName = cms.string("debugHits_MC.root"))
    process.GlobalTag = GlobalTag(process.GlobalTag, '92X_upgrade2017_realistic_v2', '')
    process.source.fileNames = [
        '/store/relval/CMSSW_9_2_3/RelValZMM_13/GEN-SIM-RECO/PUpmx25ns_92X_upgrade2017_realistic_v2_earlyBS2017-v1/10000/1A171A5A-4D51-E711-B321-0025905A6084.root',
        '/store/relval/CMSSW_9_2_3/RelValZMM_13/GEN-SIM-RECO/PUpmx25ns_92X_upgrade2017_realistic_v2_earlyBS2017-v1/10000/5EA8A130-5151-E711-B44C-0025905A6138.root',
        '/store/relval/CMSSW_9_2_3/RelValZMM_13/GEN-SIM-RECO/PUpmx25ns_92X_upgrade2017_realistic_v2_earlyBS2017-v1/10000/7A7BDC00-4D51-E711-B8D2-0CC47A4D7602.root',
        '/store/relval/CMSSW_9_2_3/RelValZMM_13/GEN-SIM-RECO/PUpmx25ns_92X_upgrade2017_realistic_v2_earlyBS2017-v1/10000/8AD53130-4F51-E711-A9FC-0025905A6138.root',
        '/store/relval/CMSSW_9_2_3/RelValZMM_13/GEN-SIM-RECO/PUpmx25ns_92X_upgrade2017_realistic_v2_earlyBS2017-v1/10000/A460431D-4C51-E711-96AA-0CC47A7452DA.root',
        '/store/relval/CMSSW_9_2_3/RelValZMM_13/GEN-SIM-RECO/PUpmx25ns_92X_upgrade2017_realistic_v2_earlyBS2017-v1/10000/AEFA23A6-4851-E711-B314-003048FFD772.root',
        '/store/relval/CMSSW_9_2_3/RelValZMM_13/GEN-SIM-RECO/PUpmx25ns_92X_upgrade2017_realistic_v2_earlyBS2017-v1/10000/C2F6268C-4F51-E711-A9E5-0CC47A7C3572.root',
        '/store/relval/CMSSW_9_2_3/RelValZMM_13/GEN-SIM-RECO/PUpmx25ns_92X_upgrade2017_realistic_v2_earlyBS2017-v1/10000/FC2154C1-5051-E711-ABCA-0CC47A4C8F08.root',
    ]
