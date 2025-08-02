import FWCore.ParameterSet.Config as cms

#two modules here
#recieves here the parameters to distiguish phi from rho

phiRecmesonStruct = cms.EDProducer("ScPhase2RecMeson",
    src = cms.InputTag("scPhase2PuppiRawToDigiStruct"),
    mesonType = cms.string("phi"),
    runStruct = cms.bool(True)
)

rhoRecmesonStruct = cms.EDProducer("ScPhase2RecMeson",
    src = cms.InputTag("scPhase2PuppiRawToDigiStruct"),
    mesonType = cms.string("rho"),
    runStruct = cms.bool(True)
)

w3piStruct = cms.EDProducer("ScPhase2PuppiW3PiDemo",
    src = cms.InputTag("scPhase2PuppiRawToDigiStruct"),
)

wdsgStruct = cms.EDProducer("ScPhase2PuppiWDsGammaDemo",
    srcPuppi = cms.InputTag("scPhase2PuppiRawToDigiStruct"),
    srcTkEm = cms.InputTag("scPhase2TkEmRawToDigiStruct"),
)

wpigStruct = cms.EDProducer("ScPhase2PuppiWPiGammaDemo",
    srcPuppi = cms.InputTag("scPhase2PuppiRawToDigiStruct"),
    srcTkEm = cms.InputTag("scPhase2TkEmRawToDigiStruct"),
)

hrhogStruct = cms.EDProducer("ScPhase2PuppiHRhoGammaDemo",
    srcPuppi = cms.InputTag("scPhase2PuppiRawToDigiStruct"),
    srcTkEm = cms.InputTag("scPhase2TkEmRawToDigiStruct"),
)

hphigStruct = cms.EDProducer("ScPhase2PuppiHPhiGammaDemo",
    srcPuppi = cms.InputTag("scPhase2PuppiRawToDigiStruct"),
    srcTkEm = cms.InputTag("scPhase2TkEmRawToDigiStruct"),
)

hjpsigStruct = cms.EDProducer("ScPhase2PuppiHJPsiGammaDemo",
    srcPuppi = cms.InputTag("scPhase2PuppiRawToDigiStruct"),
    srcTkEm = cms.InputTag("scPhase2TkEmRawToDigiStruct"),
)

h2rhoStruct = cms.EDProducer("ScPhase2PuppiH2RhoDemo",
    src = cms.InputTag("scPhase2PuppiRawToDigiStruct"),
)

h2phiStruct = cms.EDProducer("ScPhase2PuppiH2PhiDemo",
    src = cms.InputTag("scPhase2PuppiRawToDigiStruct"),
)

hphijpsiStruct = cms.EDProducer("ScPhase2PuppiHPhiJPsiDemo",
    src = cms.InputTag("scPhase2PuppiRawToDigiStruct"),
)