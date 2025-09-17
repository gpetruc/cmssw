import FWCore.ParameterSet.Config as cms


hjpsigammaRecMesonStruct = cms.EDProducer("ScPhase2RecMesonHJPsiGamma",
    # srcMeson = cms.InputTag("jpsiRecmesonStruct"),
    srcMeson = cms.InputTag("recMesonStruct", "jpsi"),
    srcGamma = cms.InputTag("recIsoTkEmStruct"),
)

h2phiRecMesonStruct = cms.EDProducer("ScPhase2RecMesonH2Phi",
    #  src = cms.InputTag("phiRecMesonStruct"),
    src = cms.InputTag("recMesonStruct", "phi"),
)

hphigammaRecMesonStruct = cms.EDProducer("ScPhase2RecMesonHPhiGamma",
    # srcMeson = cms.InputTag("phiRecMesonStruct"),
    srcMeson = cms.InputTag("recMesonStruct", "phi"),
    srcGamma = cms.InputTag("recIsoTkEmStruct"),
)

hphijpsiRecMesonStruct = cms.EDProducer("ScPhase2RecMesonHPhiJPsi",
    # srcPhi = cms.InputTag("phiRecMesonStruct"),
    srcPhi = cms.InputTag("recMesonStruct", "phi"),
    # srcJPsi = cms.InputTag("jpsiRecmesonStruct"),
    srcJPsi = cms.InputTag("recMesonStruct", "jpsi"),
)

h2rhoRecMesonStruct = cms.EDProducer("ScPhase2RecMesonH2Rho",
    #  src = cms.InputTag("rhoRecmesonStruct"),
    src = cms.InputTag("recMesonStruct", "rho"),
)

hrhogammaRecMesonStruct = cms.EDProducer("ScPhase2RecMesonHRhoGamma",
    # srcMeson = cms.InputTag("rhoRecmesonStruct"),
    srcMeson = cms.InputTag("recMesonStruct", "rho"),
    srcGamma = cms.InputTag("recIsoTkEmStruct"),
)

z2phiRecMesonStruct = cms.EDProducer("ScPhase2RecMesonZ2Phi",
    #  src = cms.InputTag("phiRecMesonStruct"),
    src = cms.InputTag("recMesonStruct", "phi"),
)

z2rhoRecMesonStruct = cms.EDProducer("ScPhase2RecMesonZ2Rho",
    #  src = cms.InputTag("rhoRecmesonStruct"),
    src = cms.InputTag("recMesonStruct", "rho"),
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