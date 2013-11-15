import FWCore.ParameterSet.Config as cms

import RecoPixelVertexing.PixelVertexFinding.PixelVertexes_cfi
import RecoTracker.IterativeTracking.InitialStep_cff
htCkfTrajectoryBuilder = RecoTracker.IterativeTracking.InitialStep_cff.initialStepTrajectoryBuilder.clone(
    ComponentName = cms.string('htCkfTrajectoryBuilder'),
    minNrOfHitsForRebuild = cms.int32(0),
    requireSeedHitsInRebuild = cms.bool(False),
)
htChi2Est = cms.ESProducer("Chi2MeasurementEstimatorESProducer",
    MaxChi2 = cms.double(60.0),
    nSigma = cms.double(3.0),
    ComponentName = cms.string('htChi2Est')
)

pixelVerticesZero = RecoPixelVertexing.PixelVertexFinding.PixelVertexes_cfi.pixelVertices.clone()

testHTSeeds = cms.EDProducer("TestHT",
    beamSpot = cms.InputTag("offlineBeamSpot"),
    #vertices = cms.InputTag("pixelVerticesZero"),
    vertices = cms.InputTag("pixelVertices"),
    #vertices = cms.InputTag("offlinePrimaryVertices"),
    vertexSelection = cms.string("ndof > 4"),
    pixelHits = cms.InputTag('siPixelRecHits'),
    stripHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
    stripHits2D = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
    measurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
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
    ptSteps = cms.vdouble(0.0, 0.8), 
    #ptSteps = cms.vdouble(0.0, 1.0, 0.5),
    # seed builder
    seedBuilderConfig = cms.PSet(
        propagator = cms.string('PropagatorWithMaterial'),
        propagatorOpposite = cms.string('PropagatorWithMaterialOpposite'),
        #TTRHBuilder = cms.string('TTRHBuilderWithoutAngle4MixedTriplets'), ## boh?
        TTRHBuilder = cms.string('WithTrackAngle'), 
        chi2MeasurementEstimator = cms.string("htChi2Est"),
        NavigationSchool  = cms.string('SimpleNavigationSchool'),
        TrajectoryBuilder = cms.string('htCkfTrajectoryBuilder'),
        TrajectoryCleaner = cms.string('TrajectoryCleanerBySharedHits'),
        cleanTrajectoryAfterInOut = cms.bool(True),
        useHitsSplitting = cms.bool(True),
        # cuts for pair seeding
        dist2dCut = cms.double(0.10),
        detaCut = cms.double(0.10),
        # cuts for additional hits
        dist2dCorrCut = cms.double(0.10),
        minHits = cms.uint32(4),
        startingCovariance = cms.vdouble(
            100., # 1/pt,
            100., #  lambda
            100., # phi
            100., # x_t
            100., # y_t
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

testHT = cms.Sequence(pixelVerticesZero + 
    testHTSeeds * testHTTracks +
    testHTSeedsAsTk + testHTClustersAsTk + testHTCandiatesAsTk
)
    
