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

scPhase2TrackerMuonMasked = cms.EDProducer("MaskOrbitBxScoutingTrackerMuon",
    dataTag = cms.InputTag("scPhase2TrackerMuonRawToDigiStruct"),
    selectBxs = cms.InputTag("scPhase2SelectedBXs","SelBx"),
)

s_maskedCollections = cms.Sequence(
    scPhase2SelectedBXs +
    scPhase2PuppiMasked +
    scPhase2TkEmMasked +
    scPhase2TkEleMasked +
    scPhase2TrackerMuonMasked
)