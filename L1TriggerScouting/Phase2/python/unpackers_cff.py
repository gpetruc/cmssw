import FWCore.ParameterSet.Config as cms

scPhase2PuppiRawToDigiStruct = cms.EDProducer('ScPhase2PuppiRawToDigi',
  src = cms.InputTag('rawDataCollector'),
  fedIDs = cms.vuint32(),
)

scPhase2TkEmRawToDigiStruct = cms.EDProducer('ScPhase2TkEmRawToDigi',
  src = cms.InputTag('rawDataCollector'),
  fedIDs = cms.vuint32(),
)

scPhase2TrackerMuonRawToDigiStruct = cms.EDProducer('ScPhase2TrackerMuonRawToDigi',
  src = cms.InputTag('rawDataCollector'),
  fedIDs = cms.vuint32(),
)

goodOrbitsByNBX = cms.EDFilter("GoodOrbitNBxSelector",
    unpackers = cms.VInputTag(
                    cms.InputTag("scPhase2PuppiRawToDigiStruct"),
                    cms.InputTag("scPhase2TkEmRawToDigiStruct"),
                    cms.InputTag("scPhase2TrackerMuonRawToDigiStruct"),
                ),
    nbxMin = cms.uint32(3564)
)

s_unpackers = cms.Sequence(
   scPhase2PuppiRawToDigiStruct +
   scPhase2TkEmRawToDigiStruct +
   scPhase2TrackerMuonRawToDigiStruct +
   goodOrbitsByNBX 
)