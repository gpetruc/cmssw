import FWCore.ParameterSet.Config as cms

recIsoTkEmStruct = cms.EDProducer("ScPhase2RecIsoTkEm",
    src = cms.InputTag("scPhase2PuppiRawToDigiStruct"),
    srcTkEm = cms.InputTag("scPhase2TkEmRawToDigiStruct"),
    minPtGamma = cms.double(20),
    isolationMinDeltaR = cms.double(0.05),
    isolationMaxDeltaR = cms.double(0.25),
    maxRelIso = cms.double(0.25)
)

mesonTypes = cms.VPSet(
    cms.PSet(
        name = cms.string("phi"),
        minMesonMass = cms.double(0.95),
        maxMesonMass = cms.double(1.25),
        dauMass1 = cms.double(0.4937),
        dauMass2 = cms.double(0.4937),
        pdgId = cms.int32(333),
    ),
    cms.PSet(
        name = cms.string("rho"),
        minMesonMass = cms.double(0.40),
        maxMesonMass = cms.double(1.30),
        dauMass1 = cms.double(0.1396),
        dauMass2 = cms.double(0.1396),
        pdgId = cms.int32(113),
    ),
    cms.PSet(
        name = cms.string("jpsi"),
        minMesonMass = cms.double(2.50),
        maxMesonMass = cms.double(3.50),
        dauMass1 = cms.double(0.1057),
        dauMass2 = cms.double(0.1057),
        pdgId = cms.int32(443),
        muonDaughters = cms.bool(True),
    )
)

puppiRecMesonStruct = cms.EDProducer("ScPhase2PuppiRecMesonAll",
    src = cms.InputTag("scPhase2PuppiRawToDigiStruct"),
    mesonTypes = mesonTypes,
    isolationMinDeltaR = cms.double(0.05),
    isolationMaxDeltaR = cms.double(0.25),
    isolationMaxDeltaZ = cms.double(1),
    isolationMaxRelIso = cms.double(0.5),
    maxDeltaRDaus = cms.double(0.40),
    maxDeltaZ = cms.double(1),
    minPtDau = cms.double(5.0),
)

tkMuonRecMesonStruct = cms.EDProducer("ScPhase2TrackerMuonRecMesonAll",
    src = cms.InputTag("scPhase2TrackerMuonRawToDigiStruct"),
    mesonTypes = cms.VPSet(*[m for m in mesonTypes if m.name.value() == "jpsi"]),
    isolationMinDeltaR = cms.double(0.05),
    isolationMaxDeltaR = cms.double(0.25),
    isolationMaxDeltaZ = cms.double(1),
    isolationMaxRelIso = cms.double(0.5),
    maxDeltaRDaus = cms.double(0.40),
    maxDeltaZ = cms.double(1),
    minPtDau = cms.double(2.0),
)


tkEleRecMesonStruct = cms.EDProducer("ScPhase2TkEleRecMesonAll",
    src = cms.InputTag("scPhase2TkEmRawToDigiStruct"),
    mesonTypes = cms.VPSet(*[m.clone(dauMass1=0.0005, dauMass2=0.0005) for m in mesonTypes if m.name.value() == "jpsi"]),
    isolationMinDeltaR = cms.double(0.05),
    isolationMaxDeltaR = cms.double(0.25),
    isolationMaxDeltaZ = cms.double(1),
    isolationMaxRelIso = cms.double(0.5),
    maxDeltaRDaus = cms.double(0.40),
    maxDeltaZ = cms.double(1),
    minPtDau = cms.double(5.0),
)

ttrackRecMesonStruct = cms.EDProducer("ScPhase2TTrackRecMesonAll",
    src = cms.InputTag("scPhase2TrackerTrackRawToDigiStruct"),
    mesonTypes = mesonTypes,
    isolationMinDeltaR = cms.double(0.05),
    isolationMaxDeltaR = cms.double(0.25),
    isolationMaxDeltaZ = cms.double(1),
    isolationMaxRelIso = cms.double(0.5),
    maxDeltaRDaus = cms.double(0.40),
    maxDeltaZ = cms.double(1),
    minPtDau = cms.double(5.0),
)

candRecoTasks = cms.Task(
    recIsoTkEmStruct,
    puppiRecMesonStruct,
    tkMuonRecMesonStruct,
    tkEleRecMesonStruct,
    ttrackRecMesonStruct,
)