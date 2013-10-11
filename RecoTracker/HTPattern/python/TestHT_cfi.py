import FWCore.ParameterSet.Config as cms

import RecoPixelVertexing.PixelVertexFinding.PixelVertexes_cfi

pixelVerticesZero = RecoPixelVertexing.PixelVertexFinding.PixelVertexes_cfi.pixelVertices.clone()

testHT = cms.EDProducer("TestHT",
    beamSpot = cms.InputTag("offlineBeamSpot"),
    vertices = cms.InputTag("pixelVerticesZero"),
    #vertices = cms.InputTag("pixelVertices"),
    #vertices = cms.InputTag("offlinePrimaryVertices"),
    pixelHits = cms.InputTag('siPixelRecHits'),
    stripHits2D = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
)


