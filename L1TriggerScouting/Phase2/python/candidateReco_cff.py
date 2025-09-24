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

mesonTypes = cms.VPSet(
    cms.PSet(
        name = cms.string("phi"),
        minMesonMass = cms.double(0.95),
        maxMesonMass = cms.double(1.25),
        dmass1 = cms.double(0.4937),
        dmass2 = cms.double(0.4937)
    ),
    cms.PSet(
        name = cms.string("rho"),
        minMesonMass = cms.double(0.40),
        maxMesonMass = cms.double(1.30),
        dmass1 = cms.double(0.1396),
        dmass2 = cms.double(0.1396)
    ),
    cms.PSet(
        name = cms.string("jpsi"),
        minMesonMass = cms.double(2.50),
        maxMesonMass = cms.double(3.50),
        dmass1 = cms.double(0.1057),
        dmass2 = cms.double(0.1057)
    )
)    

# recMesonStruct = cms.EDProducer("ScPhase2RecMesonAll",
#     src = cms.InputTag("scPhase2PuppiRawToDigiStruct"),
#     mesonTypes = mesonTypes,
#     minDeltaR = cms.double(0.05 * 0.05),
#     maxDeltaR = cms.double(0.25 * 0.25),
#     maxDeltaRDaus = cms.double(0.40 * 0.40),
#     maxDeltaZ = cms.double(1),
#     minPtDau = cms.double(5.0),
#     maxZIsolation = cms.double(1),
#     runStruct = cms.bool(True)
# )

recMesonStruct = cms.EDProducer("ScPhase2TkRecMesonAll",
    src = cms.InputTag("scPhase2TrackerTrackRawToDigiStruct"),
    mesonTypes = mesonTypes,
    minDeltaR = cms.double(0.05 * 0.05),
    maxDeltaR = cms.double(0.25 * 0.25),
    maxDeltaRDaus = cms.double(0.40 * 0.40),
    maxDeltaZ = cms.double(1),
    minPtDau = cms.double(5.0),
    maxZIsolation = cms.double(1),
    runStruct = cms.bool(True)
)

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