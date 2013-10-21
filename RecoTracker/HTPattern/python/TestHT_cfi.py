import FWCore.ParameterSet.Config as cms

import RecoPixelVertexing.PixelVertexFinding.PixelVertexes_cfi

pixelVerticesZero = RecoPixelVertexing.PixelVertexFinding.PixelVertexes_cfi.pixelVertices.clone()

testHTSeeds = cms.EDProducer("TestHT",
    beamSpot = cms.InputTag("offlineBeamSpot"),
    #vertices = cms.InputTag("pixelVerticesZero"),
    vertices = cms.InputTag("pixelVertices"),
    #vertices = cms.InputTag("offlinePrimaryVertices"),
    pixelHits = cms.InputTag('siPixelRecHits'),
    stripHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
    stripHits2D = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
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
    # seed builder
    seedBuilderConfig = cms.PSet(
        propagator = cms.string('PropagatorWithMaterial'),
        TTRHBuilder = cms.string('TTRHBuilderWithoutAngle4MixedTriplets'), ## boh?
        dist2dCut = cms.double(0.05),
        dist2dCorrCut = cms.double(0.15),
        minHits = cms.uint32(5),
        startingCovariance = cms.vdouble(
            100., # 1/pt,
            100., #  lambda
            100., # phi
            100., # x_t
            100., # y_t
        ),
    ),
    # for debugging
    debugger = cms.untracked.bool(False),
    tracks = cms.InputTag("generalTracks"),
)

import RecoTracker.CkfPattern.CkfTrackCandidates_cfi
testHTCandidates = RecoTracker.CkfPattern.CkfTrackCandidates_cfi.ckfTrackCandidates.clone(
    src = cms.InputTag("testHTSeeds"),
)
import RecoTracker.TrackProducer.TrackProducer_cfi
testHTTracks = RecoTracker.TrackProducer.TrackProducer_cfi.TrackProducer.clone(
    src = cms.InputTag("testHTCandidates"),
)

testHTSeedsAsTk = cms.EDProducer("FakeTrackProducerFromSeed",
    src = cms.InputTag("testHTSeeds")
)
testHTClustersAsTk = cms.EDProducer("FakeTrackProducerFromSeed",
    src = cms.InputTag("testHTSeeds","clusters")
)
testHTCandiatesAsTk = cms.EDProducer("FakeTrackProducerFromCandidate",
    src = cms.InputTag("testHTCandidates")
)

testHT = cms.Sequence(pixelVerticesZero + 
    testHTSeeds * testHTCandidates * testHTTracks +
    testHTSeedsAsTk + testHTClustersAsTk + testHTCandiatesAsTk
)
    
