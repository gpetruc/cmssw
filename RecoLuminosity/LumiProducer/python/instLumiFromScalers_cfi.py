import FWCore.ParameterSet.Config as cms

instLumiFromScalers = cms.EDProducer("InstLumiAsFloatProducer",
    lumiScalers = cms.InputTag("scalersRawToDigi"),
)
