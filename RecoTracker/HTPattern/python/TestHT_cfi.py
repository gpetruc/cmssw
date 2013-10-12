import FWCore.ParameterSet.Config as cms

import RecoPixelVertexing.PixelVertexFinding.PixelVertexes_cfi

pixelVerticesZero = RecoPixelVertexing.PixelVertexFinding.PixelVertexes_cfi.pixelVertices.clone()

testHT = cms.EDProducer("TestHT",
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
    layerSeedCut2d = cms.uint32(5),
    layerSeedCut3d = cms.uint32(4),
    # for debugging
    tracks = cms.InputTag("generalTracks"),
)


