import FWCore.ParameterSet.Config as cms

zdeeStruct = cms.EDProducer("ScPhase2TkEmDarkPhotonDiEle",
    srcTkEm = cms.InputTag("scPhase2TkEmRawToDigiStruct")
)

