import FWCore.ParameterSet.Config as cms

recIsoTkEmStruct = cms.EDProducer("ScPhase2RecIsoTkEm",
    src = cms.InputTag("scPhase2PuppiRawToDigiStruct"),
    srcTkEm = cms.InputTag("scPhase2TkEmRawToDigiStruct"),
    runStruct = cms.bool(True),
    minPtGamma = cms.double(20),
    minDeltaR = cms.double(0.05 * 0.05),
    maxDeltaR = cms.double(0.25 * 0.25),
    maxIso = cms.double(0.25)
)

recMesonStruct = cms.EDProducer("ScPhase2RecMesonAll",
    src = cms.InputTag("scPhase2PuppiRawToDigiStruct"),
    mesonTypes = cms.vstring("phi", "rho", "jpsi"),
    runStruct = cms.bool(True)
)

# recMesonStruct = cms.EDProducer("ScPhase2TkRecMesonAll",
#     src = cms.InputTag("scPhase2TrackerTrackRawToDigiStruct"),
#     mesonTypes = cms.vstring("phi", "rho", "jpsi"),
#     runStruct = cms.bool(True)
# )

recMesonPhiStruct = cms.EDProducer("ScPhase2RecMeson",
    src = cms.InputTag("scPhase2PuppiRawToDigiStruct"),
    minPtDau = cms.double(5.0),
    maxDeltaRDaus = cms.double(0.40 * 0.40),
    minDeltaR = cms.double(0.05 * 0.05),
    maxDeltaR = cms.double(0.25 * 0.25),
    mesonType = cms.string("phi"),
    runStruct = cms.bool(True)
)

recMesonRhoStruct = cms.EDProducer("ScPhase2RecMeson",
    src = cms.InputTag("scPhase2PuppiRawToDigiStruct"),
    minPtDau = cms.double(5.0),
    maxDeltaRDaus = cms.double(0.40 * 0.40),
    minDeltaR = cms.double(0.05 * 0.05),
    maxDeltaR = cms.double(0.25 * 0.25),
    mesonType = cms.string("rho"),
    runStruct = cms.bool(True)
)

recMesonJpsiStruct = cms.EDProducer("ScPhase2RecMeson",
    src = cms.InputTag("scPhase2PuppiRawToDigiStruct"),
    minPtDau = cms.double(5.0),
    maxDeltaRDaus = cms.double(0.40 * 0.40),
    minDeltaR = cms.double(0.05 * 0.05),
    maxDeltaR = cms.double(0.25 * 0.25),
    mesonType = cms.string("jpsi"),
    runStruct = cms.bool(True)
)