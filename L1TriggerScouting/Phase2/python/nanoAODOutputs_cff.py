import FWCore.ParameterSet.Config as cms

# Puppi candidates
scPhase2PuppiStructToTable = cms.EDProducer("ScPuppiToOrbitFlatTable",
    src = cms.InputTag("scPhase2PuppiRawToDigiStruct"),
    name = cms.string("L1Puppi"),
    doc = cms.string("L1Puppi candidates from Correlator Layer 2"),
)

scPhase2PuppiMaskedStructToTable = scPhase2PuppiStructToTable.clone(
    src = "scPhase2PuppiMasked"
)

# Isolated photons
scPhase2RecIsoTkEmStructToTable = cms.EDProducer("ScIsoTkEmToOrbitFlatTable",
    src = cms.InputTag("recIsoTkEmStruct"),
    name = cms.string("recIsoTkEm"),
    doc = cms.string("Selected isolated TkEm candidates"),
)

scPhase2RecIsoTkEmMaskedStructToTable = scPhase2RecIsoTkEmStructToTable.clone(
    src = "scPhase2RecIsoTkEmMasked"
)

# Reconstructed mesons
scPhase2RecMesonPhiStructToTable = cms.EDProducer("ScRecMesonToOrbitFlatTable",
    src = cms.InputTag("puppiRecMesonStruct", "phi"),
    name = cms.string("recMesonPhi"),
    doc = cms.string("Reconstructed Phi Meson candidates"),
)

scPhase2RecMesonPhiMaskedStructToTable = scPhase2RecMesonPhiStructToTable.clone(
    src = "scPhase2RecMesonPhiMasked"
)

scPhase2RecMesonRhoStructToTable = cms.EDProducer("ScRecMesonToOrbitFlatTable",
    src = cms.InputTag("puppiRecMesonStruct", "rho"),
    name = cms.string("recMesonRho"),
    doc = cms.string("Reconstructed Rho Meson candidates"),
)

scPhase2RecMesonRhoMaskedStructToTable = scPhase2RecMesonRhoStructToTable.clone(
    src = "scPhase2RecMesonRhoMasked"
)

scPhase2RecMesonJpsiStructToTable = cms.EDProducer("ScRecMesonToOrbitFlatTable",
    src = cms.InputTag("puppiRecMesonStruct", "jpsi"),
    name = cms.string("recMesonJpsi"),
    doc = cms.string("Reconstructed J/psi Meson candidates"),
)

scPhase2RecMesonJpsiMaskedStructToTable = scPhase2RecMesonJpsiStructToTable.clone(
    src = "scPhase2RecMesonJpsiMasked"
)

# E/Gamma candidates
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

# Tracker tracks
scPhase2TrackerTrackStructToTable = cms.EDProducer("ScTrackerTrackToOrbitFlatTable",
    src = cms.InputTag("scPhase2TrackerTrackRawToDigiStruct"),
    name = cms.string("L1TTrack"),
    doc = cms.string("L1TrackerTrack candidates from GTT"),
)

scPhase2TrackerTrackMaskedStructToTable = scPhase2TrackerTrackStructToTable.clone(
    src = "scPhase2TrackerTrackMasked"
)

# Tracker muons
scPhase2TrackerMuonStructToTable = cms.EDProducer("ScTrackerMuonToOrbitFlatTable",
    src = cms.InputTag("scPhase2TrackerMuonRawToDigiStruct"),
    name = cms.string("L1TrackerMuon"),
    doc = cms.string("L1TrackerMuon candidates from GMT"),
)

scPhase2TrackerMuonMaskedStructToTable = scPhase2TrackerMuonStructToTable.clone(
    src = "scPhase2TrackerMuonMasked"
)

# PF candidates
scPhase2PFStructToTable = cms.EDProducer("ScPuppiToOrbitFlatTable",
    src = cms.InputTag("scPhase2PFRawToDigiStruct"),
    name = cms.string("L1PF"),
    doc = cms.string("L1PF candidates from Correlator Layer 1"),
)

scPhase2PFMaskedStructToTable = scPhase2PFStructToTable.clone(
    src = "scPhase2PFMasked"
)

scPhase2TkEgTableProducersTask = cms.Task(
    scPhase2TkEmStructToTable,
    scPhase2TkEleStructToTable,
)

scPhase2RecMesonStructToTable = cms.Task(
    scPhase2RecMesonPhiStructToTable,
    scPhase2RecMesonRhoStructToTable,
    scPhase2RecMesonJpsiStructToTable
)

tableProducersTask = cms.Task(
    scPhase2PuppiStructToTable,
    scPhase2TkEgTableProducersTask,
    scPhase2TrackerMuonStructToTable,
    scPhase2TrackerTrackStructToTable,
    scPhase2PFStructToTable,
    scPhase2RecIsoTkEmStructToTable,
    scPhase2RecMesonStructToTable,
)

scPhase2TkEgMaskedTableProducersTask = cms.Task(
    scPhase2TkEmMaskedStructToTable,
    scPhase2TkEleMaskedStructToTable,
)

scPhase2RecMesonMaskedStructToTable = cms.Task(
    scPhase2RecMesonPhiMaskedStructToTable,
    scPhase2RecMesonRhoMaskedStructToTable,
    scPhase2RecMesonJpsiMaskedStructToTable
)

maskedTableProducersTask = cms.Task(
    scPhase2PuppiMaskedStructToTable,
    scPhase2TkEgMaskedTableProducersTask,
    scPhase2TrackerMuonMaskedStructToTable,
    scPhase2TrackerTrackMaskedStructToTable,
    scPhase2PFMaskedStructToTable,
    scPhase2RecIsoTkEmMaskedStructToTable,
    scPhase2RecMesonMaskedStructToTable,
)

scPhase2NanoAll = cms.OutputModule("OrbitNanoAODOutputModule",
    fileName = cms.untracked.string("all.root"),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring()),
    outputCommands = cms.untracked.vstring("drop *",
        "keep l1ScoutingRun3OrbitFlatTable_scPhase2PuppiStructToTable_*_*",
        "keep l1ScoutingRun3OrbitFlatTable_scPhase2RecIsoTkEmStructToTable_*_*",
        "keep l1ScoutingRun3OrbitFlatTable_scPhase2RecMesonPhiStructToTable_*_*",
        "keep l1ScoutingRun3OrbitFlatTable_scPhase2RecMesonRhoStructToTable_*_*",
        "keep l1ScoutingRun3OrbitFlatTable_scPhase2RecMesonJpsiStructToTable_*_*",
        "keep l1ScoutingRun3OrbitFlatTable_scPhase2TkEmStructToTable_*_*",
        "keep l1ScoutingRun3OrbitFlatTable_scPhase2TkEleStructToTable_*_*",
        "keep l1ScoutingRun3OrbitFlatTable_scPhase2TrackerMuonStructToTable_*_*",
        "keep l1ScoutingRun3OrbitFlatTable_scPhase2TrackerTrackStructToTable_*_*",
        "keep l1ScoutingRun3OrbitFlatTable_scPhase2PFStructToTable_*_*"
    ),
    compressionLevel = cms.untracked.int32(4),
    compressionAlgorithm = cms.untracked.string("LZ4"),
)

scPhase2NanoSelected = cms.OutputModule("OrbitNanoAODOutputModule",
    fileName = cms.untracked.string("selected.root"),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring()),
    selectedBx = cms.InputTag("scPhase2SelectedBXs","SelBx"),
    outputCommands = cms.untracked.vstring("drop *",
        "keep l1ScoutingRun3OrbitFlatTable_scPhase2PuppiMaskedStructToTable_*_*",
        "keep l1ScoutingRun3OrbitFlatTable_scPhase2RecIsoTkEmMaskedStructToTable_*_*",
        "keep l1ScoutingRun3OrbitFlatTable_scPhase2RecMesonPhiMaskedStructToTable_*_*",
        "keep l1ScoutingRun3OrbitFlatTable_scPhase2RecMesonRhoMaskedStructToTable_*_*",
        "keep l1ScoutingRun3OrbitFlatTable_scPhase2RecMesonJpsiMaskedStructToTable_*_*",
        "keep l1ScoutingRun3OrbitFlatTable_scPhase2TkEmMaskedStructToTable_*_*",
        "keep l1ScoutingRun3OrbitFlatTable_scPhase2TkEleMaskedStructToTable_*_*",
        "keep l1ScoutingRun3OrbitFlatTable_scPhase2TrackerMuonMaskedStructToTable_*_*",
        "keep l1ScoutingRun3OrbitFlatTable_scPhase2TrackerTrackMaskedStructToTable_*_*",
        "keep l1ScoutingRun3OrbitFlatTable_scPhase2PFMaskedStructToTable_*_*",
        "keep *_scPhase2SelectedBXs_*_*"),
    compressionLevel = cms.untracked.int32(4),
    compressionAlgorithm = cms.untracked.string("LZ4"),
)