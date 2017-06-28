import FWCore.ParameterSet.Config as cms

process = cms.Process("Occupancy")

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(
       #'file:debug_Zmm_lostHits.root',
       '/store/express/Run2017B/ExpressPhysics/FEVT/Express-v1/000/297/178/00000/00179E1E-7055-E711-811C-02163E01283D.root',
       '/store/express/Run2017B/ExpressPhysics/FEVT/Express-v1/000/297/178/00000/0069D1D7-6355-E711-AF34-02163E0143BC.root',
       '/store/express/Run2017B/ExpressPhysics/FEVT/Express-v1/000/297/178/00000/0085E8F2-6F55-E711-A95E-02163E0133F0.root',
       '/store/express/Run2017B/ExpressPhysics/FEVT/Express-v1/000/297/178/00000/00A6F68B-6655-E711-8C52-02163E013555.root',
       '/store/express/Run2017B/ExpressPhysics/FEVT/Express-v1/000/297/178/00000/020B0331-6255-E711-9937-02163E014716.root',
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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )    

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '92X_dataRun2_Express_v2', '')

process.triggerResultsFilter = cms.EDFilter("TriggerResultsFilter",
    daqPartitions = cms.uint32(1),
    hltResults = cms.InputTag("TriggerResults","","HLT"),
    l1tIgnoreMask = cms.bool(False),
    l1tResults = cms.InputTag(""),
    l1techIgnorePrescales = cms.bool(False),
    throw = cms.bool(True),
    triggerConditions = cms.vstring('HLT_IsoMu20_v*','HLT_IsoMu24_v*')
)


process.pixelOccupancy = cms.EDProducer("PixelHitOccupancy",
    tracks = cms.InputTag("generalTracks"),
    pixelClusters = cms.InputTag("siPixelClusters"),
)

process.pixelOccupancyByLumi = cms.EDAnalyzer("PixelHitOccupancyByLumi",
    #tracker = cms.InputTag("MeasurementTrackerEvent"),
    pixelClusters = cms.InputTag("siPixelClusters"),
    maxHitsPerROC = cms.uint32(100), # stop processing LS after I have at least this number of hits on all ROCs that have hits
)
process.pre10 = cms.EDFilter("Prescaler", prescaleFactor = cms.int32(10), prescaleOffset = cms.int32(0) )
process.occupancy = cms.Path(
    process.pre10 +
    #process.MeasurementTrackerEvent +
    #process.triggerResultsFilter +
    process.pixelOccupancyByLumi
)

process.TFileService = cms.Service("TFileService", fileName = cms.string("badComponents.root"))

