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

scPhase2TrackerTrackRawToDigiStruct = cms.EDProducer('ScPhase2TrackerTrackRawToDigi',
  src = cms.InputTag('rawDataCollector'),
  fedIDs = cms.vuint32(),
  nFitPars = cms.uint32(4)
)

goodOrbitsByNBX = cms.EDFilter("GoodOrbitNBxSelector",
    unpackers = cms.VInputTag(
                    cms.InputTag("scPhase2PuppiRawToDigiStruct"),
                    cms.InputTag("scPhase2TkEmRawToDigiStruct"),
                    cms.InputTag("scPhase2TrackerMuonRawToDigiStruct"),
                    cms.InputTag("scPhase2TrackerTrackRawToDigiStruct"),
                ),
    nbxMin = cms.uint32(3564)
)

s_unpackers = cms.Sequence(
   scPhase2PuppiRawToDigiStruct +
   scPhase2TkEmRawToDigiStruct +
   scPhase2TrackerMuonRawToDigiStruct +
   scPhase2TrackerTrackRawToDigiStruct +
   goodOrbitsByNBX
)