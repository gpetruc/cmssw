import FWCore.ParameterSet.Config as cms

scPhase2PuppiStructToTable = cms.EDProducer("ScPuppiToOrbitFlatTable",
    src = cms.InputTag("scPhase2PuppiRawToDigiStruct"),
    name = cms.string("L1Puppi"),
    doc = cms.string("L1Puppi candidates from Correlator Layer 2"),
)

scPhase2PuppiMaskedStructToTable = scPhase2PuppiStructToTable.clone(
    src = "scPhase2PuppiMasked"
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

tableProducersTask = cms.Task(
    scPhase2PuppiStructToTable,
    tableProducersTkEmTask,
    scPhase2TrackerMuonStructToTable,
)

maskedTableProducersTkEmTask = cms.Task(
    scPhase2TkEmMaskedStructToTable,
    scPhase2TkEleMaskedStructToTable,
)

maskedTableProducersTask = cms.Task(
    scPhase2PuppiMaskedStructToTable,
    maskedTableProducersTkEmTask,
    scPhase2TrackerMuonMaskedStructToTable,
)

scPhase2NanoAll = cms.OutputModule("OrbitNanoAODOutputModule",
    fileName = cms.untracked.string("all.root"),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring()),
    outputCommands = cms.untracked.vstring("drop *", 
        "keep l1ScoutingRun3OrbitFlatTable_scPhase2PuppiStructToTable_*_*", 
        "keep l1ScoutingRun3OrbitFlatTable_scPhase2TkEmStructToTable_*_*", 
        "keep l1ScoutingRun3OrbitFlatTable_scPhase2TkEleStructToTable_*_*",
        "keep l1ScoutingRun3OrbitFlatTable_scPhase2TrackerMuonStructToTable_*_*"),
    compressionLevel = cms.untracked.int32(4),
    compressionAlgorithm = cms.untracked.string("LZ4"),
)

scPhase2PuppiNanoSelected = cms.OutputModule("OrbitNanoAODOutputModule",
    fileName = cms.untracked.string("selected.root"),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring()),
    selectedBx = cms.InputTag("scPhase2SelectedBXs","SelBx"),
    outputCommands = cms.untracked.vstring("drop *",
        "keep l1ScoutingRun3OrbitFlatTable_scPhase2PuppiMaskedStructToTable_*_*",
        "keep l1ScoutingRun3OrbitFlatTable_scPhase2TkEmMaskedStructToTable_*_*",
        "keep l1ScoutingRun3OrbitFlatTable_scPhase2TkEleMaskedStructToTable_*_*",
        "keep l1ScoutingRun3OrbitFlatTable_scPhase2TrackerMuonMaskedStructToTable_*_*",
        "keep *_scPhase2SelectedBXs_*_*"),
    compressionLevel = cms.untracked.int32(4),
    compressionAlgorithm = cms.untracked.string("LZ4"),
)