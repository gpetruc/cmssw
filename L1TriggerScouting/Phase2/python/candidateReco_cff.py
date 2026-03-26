import FWCore.ParameterSet.Config as cms

recIsoTkEmStruct = cms.EDProducer("ScPhase2RecIsoTkEm",
    src = cms.InputTag("scPhase2PuppiRawToDigiStruct"),
    srcTkEm = cms.InputTag("scPhase2TkEmRawToDigiStruct"),
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
        dauMass1 = cms.double(0.4937),
        dauMass2 = cms.double(0.4937),
        pgdId = cms.int32(333),
    ),
    cms.PSet(
        name = cms.string("rho"),
        minMesonMass = cms.double(0.40),
        maxMesonMass = cms.double(1.30),
        dauMass1 = cms.double(0.1396),
        dauMass2 = cms.double(0.1396),
        pgdId = cms.int32(113),
    ),
    cms.PSet(
        name = cms.string("jpsi"),
        minMesonMass = cms.double(2.50),
        maxMesonMass = cms.double(3.50),
        dauMass1 = cms.double(0.1057),
        dauMass2 = cms.double(0.1057),
        pgdId = cms.int32(443),
        muonDaughters = cms.bool(True),
    )
)

puppiRecMesonStruct = cms.EDProducer("ScPhase2PuppiRecMesonAll",
    src = cms.InputTag("scPhase2PuppiRawToDigiStruct"),
    mesonTypes = mesonTypes.clone(),
    minDeltaR = cms.double(0.05 * 0.05),
    maxDeltaR = cms.double(0.25 * 0.25),
    maxDeltaRDaus = cms.double(0.40 * 0.40),
    maxDeltaZ = cms.double(1),
    minPtDau = cms.double(5.0),
    maxZIsolation = cms.double(1)
)

ttrackRecMesonStruct = cms.EDProducer("ScPhase2TTrackRecMesonAll",
    src = cms.InputTag("scPhase2TrackerTrackRawToDigiStruct"),
    mesonTypes = mesonTypes.clone(),
    minDeltaR = cms.double(0.05 * 0.05),
    maxDeltaR = cms.double(0.25 * 0.25),
    maxDeltaRDaus = cms.double(0.40 * 0.40),
    maxDeltaZ = cms.double(1),
    minPtDau = cms.double(5.0),
    maxZIsolation = cms.double(1)
)