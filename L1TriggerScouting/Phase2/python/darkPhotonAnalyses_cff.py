import FWCore.ParameterSet.Config as cms

dimuStruct = cms.EDProducer("ScPhase2TrackerMuonDiMuDemo",
    src = cms.InputTag("scPhase2TrackerMuonRawToDigiStruct"),
)