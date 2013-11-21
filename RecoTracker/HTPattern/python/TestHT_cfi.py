import FWCore.ParameterSet.Config as cms

import RecoPixelVertexing.PixelVertexFinding.PixelVertexes_cfi
import RecoTracker.IterativeTracking.InitialStep_cff
htCkfTrajectoryBuilder = RecoTracker.IterativeTracking.InitialStep_cff.initialStepTrajectoryBuilder.clone(
    ComponentName = cms.string('htCkfTrajectoryBuilder'),
    minNrOfHitsForRebuild = cms.int32(0),
    requireSeedHitsInRebuild = cms.bool(False),
)
htChi2Est = cms.ESProducer("Chi2MeasurementEstimatorESProducer",
    MaxChi2 = cms.double(30.0),
    nSigma = cms.double(3.0),
    ComponentName = cms.string('htChi2Est')
)

pixelVerticesZero = RecoPixelVertexing.PixelVertexFinding.PixelVertexes_cfi.pixelVertices.clone()
pixelVerticesFake = cms.EDProducer("HTVertexCorrectionStudy",
    tracks = cms.InputTag("generalTracks"),
    trackSelection = cms.string("pt > 0.5 && numberOfValidHits > 6"),
    dz = cms.double(0.3),
    dxy = cms.double(0.0),
    phi = cms.double(0.0),
)
import RecoTracker.TransientTrackingRecHit.TransientTrackingRecHitBuilder_cfi
ttrhbwrBootstrap = RecoTracker.TransientTrackingRecHit.TransientTrackingRecHitBuilder_cfi.ttrhbwr.clone(
    ComputeCoarseLocalPositionFromDisk = True,
    ComponentName = 'WithTrackAngleAndBootstrap'
)    


testHTSeeds = cms.EDProducer("TestHT",
    beamSpot = cms.InputTag("offlineBeamSpot"),
    #vertices = cms.InputTag("pixelVerticesZero"),
    #vertices = cms.InputTag("pixelVerticesFake"),
    vertices = cms.InputTag("pixelVertices"),
    #vertices = cms.InputTag("offlinePrimaryVertices"),
    vertexSelection = cms.string("ndof > 4"),
    pixelHits = cms.InputTag('siPixelRecHits'),
    stripHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
    stripHits2D = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
    measurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    # seeding 
    seed2d = cms.bool(False),
    seed3d = cms.bool(False),
    seedMixed = cms.bool(True),
    # binning
    etabins2d = cms.uint32(32),
    phibins2d = cms.uint32(64),
    etabins3d = cms.uint32(128),
    phibins3d = cms.uint32(64),
    # clustering
    layerSeedCut2d = cms.uint32(4),
    layerSeedCut3d = cms.uint32(3),
    layerCut2d     = cms.uint32(5),
    layerCut3d     = cms.uint32(5),
    layerMoreCut   = cms.uint32(5),
    # pt clustering
    #ptSteps = cms.vdouble( 0.0,       2.5,      1.0,      0.7,      0.5 ), 
    #ptEdges = cms.vdouble( 3.0,-3.0,  1.5,4.0,  0.7,2.0,  0.5,1.0,  0.4,0.8 ), 
    ptSteps = cms.vdouble(0.0, 0.8),
    ptEdges = cms.vdouble(), 
    # seed builder
    seedBuilderConfig = cms.PSet(
        propagator = cms.string('PropagatorWithMaterial'),
        propagatorOpposite = cms.string('PropagatorWithMaterialOpposite'),
        #TTRHBuilder = cms.string('TTRHBuilderWithoutAngle4MixedTriplets'), ## boh?
        #TTRHBuilder = cms.string('WithTrackAngle'), 
        TTRHBuilder = cms.string('WithTrackAngleAndBootstrap'), 
        chi2MeasurementEstimator = cms.string("htChi2Est"),
        NavigationSchool  = cms.string('SimpleNavigationSchool'),
        TrajectoryBuilder = cms.string('htCkfTrajectoryBuilder'),
        TrajectoryCleaner = cms.string('TrajectoryCleanerBySharedHits'),
        cleanTrajectoryAfterInOut = cms.bool(True),
        useHitsSplitting = cms.bool(True),
        splitSeedHits = cms.bool(True),
        # cuts for pair seeding
        detaCutPair = cms.double(0.01),
        dphiCutPair = cms.double(0.02),
        pairsOnSeedCellOnly = cms.bool(True),
        # cuts for additional hits (before and after the refit)
        dphiCutHits = cms.vdouble(0.04, 0.015),
        detaCutHits = cms.vdouble(0.05, 0.02),
        minHits = cms.uint32(4),
        # minimum number of hits in a trajectory for which ckf failed
        minHitsIfNoCkf = cms.uint32(6), 
        # max consec failing microclusters
        maxFailMicroClusters = cms.uint32(5),
        startingCovariance = cms.vdouble(
            10., # 1/pt,
            4., # lambda
            4., # phi
            625., # x_t
            625., # y_t
        ),
        transientInitialStateEstimatorParameters = cms.PSet(
            propagatorAlongTISE = cms.string('PropagatorWithMaterial'),
            numberMeasurementsForFit = cms.int32(4),
            propagatorOppositeTISE = cms.string('PropagatorWithMaterialOpposite')
        ),
    ),
    # for debugging
    debugger = cms.untracked.bool(True),
    tracks = cms.InputTag("generalTracks"),
)

import RecoTracker.TrackProducer.TrackProducer_cfi
testHTTracks = RecoTracker.TrackProducer.TrackProducer_cfi.TrackProducer.clone(
    src = cms.InputTag("testHTSeeds"),
)

testHTSeedsAsTk = cms.EDProducer("FakeTrackProducerFromSeed",
    src = cms.InputTag("testHTSeeds")
)
testHTClustersAsTk = cms.EDProducer("FakeTrackProducerFromSeed",
    src = cms.InputTag("testHTSeeds","clusters")
)
testHTCandiatesAsTk = cms.EDProducer("FakeTrackProducerFromCandidate",
    src = cms.InputTag("testHTSeeds")
)

testHT = cms.Sequence(
    #pixelVerticesZero + 
    #pixelVerticesFake +
    testHTSeeds * testHTTracks +
    testHTSeedsAsTk + testHTClustersAsTk + testHTCandiatesAsTk 
)

binningStudy = cms.EDProducer("HTBinningStudies",
    tracks = cms.InputTag("TrackMCQuality"),
    trackSelection = cms.string("pt > 0.4 && abs(eta) < 2.4 && quality('highPurity') && quality('qualitySize')"),    
)
binStudy = cms.Sequence(
    binningStudy
)
finalStudy = cms.EDProducer("HTFinalOutcomeStudy",
    tracksCkf = cms.InputTag("TrackMCQuality"),
    tracksHT  = cms.InputTag("TrackMCQualityHT"),
    trackSelectionCkf = cms.string("pt > 0.4 && abs(eta) < 2.4 && quality('highPurity')"),
)


from SimTracker.TrackAssociation.TrackAssociatorByHits_cfi import *
from SimTracker.TrackAssociation.TrackMCQuality_cfi import *
TrackMCQualityHT = TrackMCQuality.clone(label_tr = "testHTTracks")
trueTracks = cms.EDFilter("TrackSelector",
    src = cms.InputTag("TrackMCQuality"),
    cut = cms.string("quality('qualitySize')"),
)
testHTMC = cms.Sequence(
    TrackMCQuality + trueTracks +
    testHT +
    TrackMCQualityHT + finalStudy
)
