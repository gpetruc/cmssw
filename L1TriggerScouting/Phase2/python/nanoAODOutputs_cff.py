import FWCore.ParameterSet.Config as cms

scPhase2PuppiStructToTable = cms.EDProducer("ScPuppiToOrbitFlatTable",
    src = cms.InputTag("scPhase2PuppiRawToDigiStruct"),
    name = cms.string("L1Puppi"),
    doc = cms.string("L1Puppi candidates from Correlator Layer 2"),
)

scPhase2PuppiMaskedStructToTable = scPhase2PuppiStructToTable.clone(
    src = "scPhase2PuppiMasked"
)

scPhase2PhiRecMesonStructToTable = cms.EDProducer("ScRecMesonToOrbitFlatTable",
    src = cms.InputTag("allRecmesonStruct", "recMeson_phi"),
    name = cms.string("phiRecMeson"),
    doc = cms.string("Reconstructed Phi Meson candidates"),
)

scPhase2RhoRecMesonStructToTable = cms.EDProducer("ScRecMesonToOrbitFlatTable",
    src = cms.InputTag("allRecmesonStruct", "recMeson_rho"),
    name = cms.string("rhoRecMeson"),
    doc = cms.string("Reconstructed Rho Meson candidates"),
)

scPhase2JPsiRecMesonStructToTable = cms.EDProducer("ScRecMesonToOrbitFlatTable",
    src = cms.InputTag("allRecmesonStruct", "recMeson_jpsi"),
    name = cms.string("jpsiRecMeson"),
    doc = cms.string("Reconstructed J/Psi Meson candidates"),
)

scPhase2TkEmStructToTable = cms.EDProducer("ScTkEmToOrbitFlatTable",
    src = cms.InputTag("scPhase2TkEmRawToDigiStruct"),
    name = cms.string("L1TkEm"),
    doc = cms.string("L1TkEm candidates"),
)

scPhase2TkEmMaskedStructToTable = scPhase2TkEmStructToTable.clone(
    src = "scPhase2TkEmMasked"
)

scPhase2TkEleStructToTable = cms.EDProducer("ScTkEleToOrbitFlatTable",
    src = cms.InputTag("scPhase2TkEmRawToDigiStruct"),
    name = cms.string("L1TkEle"),
    doc = cms.string("L1TkEle candidates"),
)

scPhase2TkEleMaskedStructToTable = scPhase2TkEleStructToTable.clone(
    src = "scPhase2TkEleMasked"
)

scPhase2TrackerStructToTable = cms.EDProducer("ScTrackerToOrbitFlatTable",
    src = cms.InputTag("scPhase2TrackerRawToDigiStruct"),
    name = cms.string("L1Tracker"),
    doc = cms.string("L1Tracker candidates from GMT"),
)

scPhase2TrackerMaskedStructToTable = scPhase2TrackerStructToTable.clone(
    src = "scPhase2TrackerMasked"
)

scPhase2TrackerMuonStructToTable = cms.EDProducer("ScTrackerMuonToOrbitFlatTable",
    src = cms.InputTag("scPhase2TrackerMuonRawToDigiStruct"),
    name = cms.string("L1TrackerMuon"),
    doc = cms.string("L1TrackerMuon candidates from GMT"),
)

scPhase2TrackerMuonMaskedStructToTable = scPhase2TrackerMuonStructToTable.clone(
    src = "scPhase2TrackerMuonMasked"
)

tableProducersTkEmTask = cms.Task(
    scPhase2TkEmStructToTable,
    scPhase2TkEleStructToTable,
)

tableProducersTrackerTask = cms.Task(
    scPhase2TrackerStructToTable,
)

tableProducersTask = cms.Task(
    scPhase2PhiRecMesonStructToTable,
    scPhase2RhoRecMesonStructToTable,
    scPhase2JPsiRecMesonStructToTable,
    scPhase2PuppiStructToTable,
    tableProducersTkEmTask,
    scPhase2TrackerMuonStructToTable,
)

maskedTableProducersTkEmTask = cms.Task(
    scPhase2TkEmMaskedStructToTable,
    scPhase2TkEleMaskedStructToTable,
)

maskedTableProducersTrackerTask = cms.Task(
    scPhase2TrackerMaskedStructToTable,
)

maskedTableProducersTask = cms.Task(
    scPhase2PhiRecMesonStructToTable,
    scPhase2RhoRecMesonStructToTable,
    scPhase2JPsiRecMesonStructToTable,
    scPhase2PuppiMaskedStructToTable,
    maskedTableProducersTkEmTask,
    maskedTableProducersTrackerTask,
    scPhase2TrackerMuonMaskedStructToTable,
)

scPhase2NanoAll = cms.OutputModule("OrbitNanoAODOutputModule",
    fileName = cms.untracked.string("all.root"),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring()),
    outputCommands = cms.untracked.vstring("drop *", 
        "keep l1ScoutingRun3OrbitFlatTable_scPhase2PhiRecMesonStructToTable_*_*",
        "keep l1ScoutingRun3OrbitFlatTable_scPhase2RhoRecMesonStructToTable_*_*",
        "keep l1ScoutingRun3OrbitFlatTable_scPhase2JPsiRecMesonStructToTable_*_*",
        # "keep l1ScoutingRun3OrbitFlatTable_scPhase2PuppiStructToTable_*_*", 
        "keep l1ScoutingRun3OrbitFlatTable_scPhase2TkEmStructToTable_*_*", 
        "keep l1ScoutingRun3OrbitFlatTable_scPhase2TkEleStructToTable_*_*",
        "keep l1ScoutingRun3OrbitFlatTable_scPhase2TrackerStructToTable_*_*",
        "keep l1ScoutingRun3OrbitFlatTable_scPhase2TrackerMuonStructToTable_*_*"),
    compressionLevel = cms.untracked.int32(4),
    compressionAlgorithm = cms.untracked.string("LZ4"),
)

scPhase2PuppiNanoSelected = cms.OutputModule("OrbitNanoAODOutputModule",
    fileName = cms.untracked.string("selected.root"),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring()),
    selectedBx = cms.InputTag("scPhase2SelectedBXs","SelBx"),
    outputCommands = cms.untracked.vstring("drop *",
        # "keep l1ScoutingRun3OrbitFlatTable_scPhase2PuppiMaskedStructToTable_*_*",
        "keep l1ScoutingRun3OrbitFlatTable_scPhase2TkEmMaskedStructToTable_*_*",
        "keep l1ScoutingRun3OrbitFlatTable_scPhase2TkEleMaskedStructToTable_*_*",
        "keep l1ScoutingRun3OrbitFlatTable_scPhase2TrackerMuonMaskedStructToTable_*_*",
        "keep *_scPhase2SelectedBXs_*_*"),
    compressionLevel = cms.untracked.int32(4),
    compressionAlgorithm = cms.untracked.string("LZ4"),
)