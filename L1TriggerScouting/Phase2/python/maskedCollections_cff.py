import FWCore.ParameterSet.Config as cms

scPhase2SelectedBXs =  cms.EDFilter("FinalBxSelector",
    analysisLabels = cms.VInputTag(),
)

scPhase2PuppiMasked = cms.EDProducer("MaskOrbitBxScoutingPuppi",
    dataTag = cms.InputTag("scPhase2PuppiRawToDigiStruct"),
    selectBxs = cms.InputTag("scPhase2SelectedBXs","SelBx"),
)

scPhase2TkEmMasked = cms.EDProducer("MaskOrbitBxScoutingTkEm",
    dataTag = cms.InputTag("scPhase2TkEmRawToDigiStruct"),
    selectBxs = cms.InputTag("scPhase2SelectedBXs","SelBx"),
)

scPhase2TkEleMasked = cms.EDProducer("MaskOrbitBxScoutingTkEle",
    dataTag = cms.InputTag("scPhase2TkEmRawToDigiStruct"),
    selectBxs = cms.InputTag("scPhase2SelectedBXs","SelBx"),
)

scPhase2TrackerTrackMasked = cms.EDProducer("MaskOrbitBxScoutingTrackerTrack",
    dataTag = cms.InputTag("scPhase2TrackerTrackRawToDigiStruct"),
    selectBxs = cms.InputTag("scPhase2SelectedBXs","SelBx"),
)

scPhase2TrackerMuonMasked = cms.EDProducer("MaskOrbitBxScoutingTrackerMuon",
    dataTag = cms.InputTag("scPhase2TrackerMuonRawToDigiStruct"),
    selectBxs = cms.InputTag("scPhase2SelectedBXs","SelBx"),
)

scPhase2PFMasked = scPhase2PuppiMasked.clone(
    dataTag = cms.InputTag("scPhase2PFRawToDigiStruct"),
)

scPhase2RecIsoTkEmMasked = cms.EDProducer("MaskOrbitBxScoutingTkEm",
    dataTag = cms.InputTag("recIsoTkEmStruct"),
    selectBxs = cms.InputTag("scPhase2SelectedBXs","SelBx"),
)

scPhase2RecMesonPhiMasked = cms.EDProducer("MaskOrbitBxScoutingRecMeson2",
    dataTag = cms.InputTag("recMesonStruct", "phi"),
    selectBxs = cms.InputTag("scPhase2SelectedBXs","SelBx"),
)

scPhase2RecMesonRhoMasked = cms.EDProducer("MaskOrbitBxScoutingRecMeson2",
    dataTag = cms.InputTag("recMesonStruct", "rho"),
    selectBxs = cms.InputTag("scPhase2SelectedBXs","SelBx"),
)

scPhase2RecMesonJpsiMasked = cms.EDProducer("MaskOrbitBxScoutingRecMeson2",
    dataTag = cms.InputTag("recMesonStruct", "jpsi"),
    selectBxs = cms.InputTag("scPhase2SelectedBXs","SelBx"),
)

scPhase2RecMesonMasked = cms.Sequence(
    scPhase2RecMesonPhiMasked +
    scPhase2RecMesonRhoMasked +
    scPhase2RecMesonJpsiMasked
)

s_maskedCollections = cms.Sequence(
    scPhase2SelectedBXs +
    scPhase2PuppiMasked +
    scPhase2TkEmMasked +
    scPhase2TkEleMasked +
    scPhase2TrackerTrackMasked +
    scPhase2TrackerMuonMasked +
    scPhase2PFMasked +
    scPhase2RecIsoTkEmMasked +
    scPhase2RecMesonMasked
)