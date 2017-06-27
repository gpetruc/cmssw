import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe2")

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
#process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.suppressError = cms.untracked.vstring("patTriggerFull")
process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(
       'file:debug_Zmm_lostHits.root',
       #'/store/express/Run2017B/ExpressPhysics/FEVT/Express-v1/000/297/178/00000/00179E1E-7055-E711-811C-02163E01283D.root',
       #'/store/express/Run2017B/ExpressPhysics/FEVT/Express-v1/000/297/178/00000/0069D1D7-6355-E711-AF34-02163E0143BC.root',
       #'/store/express/Run2017B/ExpressPhysics/FEVT/Express-v1/000/297/178/00000/0085E8F2-6F55-E711-A95E-02163E0133F0.root',
       #'/store/express/Run2017B/ExpressPhysics/FEVT/Express-v1/000/297/178/00000/00A6F68B-6655-E711-8C52-02163E013555.root',
       #'/store/express/Run2017B/ExpressPhysics/FEVT/Express-v1/000/297/178/00000/020B0331-6255-E711-9937-02163E014716.root',
       #'/store/express/Run2017B/ExpressPhysics/FEVT/Express-v1/000/297/178/00000/02501F09-6A55-E711-AFF3-02163E014296.root',
       #'/store/express/Run2017B/ExpressPhysics/FEVT/Express-v1/000/297/178/00000/0279AB4C-7155-E711-81EA-02163E01274E.root',
       #'/store/express/Run2017B/ExpressPhysics/FEVT/Express-v1/000/297/178/00000/0284DD83-5C55-E711-886A-02163E014604.root',
       #'/store/express/Run2017B/ExpressPhysics/FEVT/Express-v1/000/297/178/00000/02B90ABB-5455-E711-83AC-02163E011F46.root',
       #'/store/express/Run2017B/ExpressPhysics/FEVT/Express-v1/000/297/178/00000/02BC23B6-7455-E711-A09D-02163E014200.root',
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
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_Express_v10', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_Prompt_v8', '')

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
    cut = cms.string("pt > 10 && numberOfMatchedStations >= 1 && innerTrack.isNonnull && abs(eta) < 1.0"), 
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
process.clusterInfo = cms.EDProducer("DebugPixelHits",
        pairs = cms.InputTag("tpPairs"),
        tracker = cms.InputTag("MeasurementTrackerEvent"),
        # configuraton for refitter
        DoPredictionsOnly = cms.bool(False),
        Fitter = cms.string('KFFitterForRefitInsideOut'),
        TrackerRecHitBuilder = cms.string('WithAngleAndTemplate'),
        Smoother = cms.string('KFSmootherForRefitInsideOut'),
        MuonRecHitBuilder = cms.string('MuonRecHitBuilder'),
        RefitDirection = cms.string('alongMomentum'),
        RefitRPCHits = cms.bool(True),
        Propagator = cms.string('SmartPropagatorAnyRKOpposite'),
        #Propagators
        PropagatorOpposite = cms.string("RungeKuttaTrackerPropagatorOpposite"),
        Chi2MeasurementEstimator = cms.string("HitCollectorForDebug"),
        #Error rescaling
        rescaleError = cms.double(10),
        #SiPixelQuality = cms.string(''),
        ## https://github.com/cms-sw/cmssw/blob/9b7f92a91b55fe1bf3e38435a6afd5b97dea4c9f/RecoLocalTracker/SubCollectionProducers/src/JetCoreClusterSplitter.cc#L139-L153
)

process.pixelOccupancy = cms.EDProducer("PixelHitOccupancy",
    tracks = cms.InputTag("generalTracks"),
    pixelClusters = cms.InputTag("siPixelClusters"),
)
process.occupancy = cms.Path(process.pixelOccupancy)
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
    outputCommands = cms.untracked.vstring("keep *", "drop *_*_*_TagProbe2"),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("tagAndProbe")),
)
#uncomment to skim the selected events
#process.end = cms.EndPath(process.out)

# this below probably not needed
process.TFileService = cms.Service("TFileService", fileName = cms.string("occupancy.root"))

