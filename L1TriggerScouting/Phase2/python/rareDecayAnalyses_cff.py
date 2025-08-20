import FWCore.ParameterSet.Config as cms

#two modules here
#recieves here the parameters to distiguish phi from rho

photonIsolationStruct = cms.EDProducer("ScPhase2PhotonIsolation",
    src = cms.InputTag("scPhase2PuppiRawToDigiStruct"),
    srcTkEm = cms.InputTag("scPhase2TkEmRawToDigiStruct"),
    runStruct = cms.bool(True),
    minPtGamma = cms.double(20),
    minDeltaR = cms.double(0.05 * 0.05),
    maxDeltaR = cms.double(0.25 * 0.25),
    maxIsol = cms.double(0.25)
)

jpsiRecmesonStruct = cms.EDProducer("ScPhase2RecMeson",
    src = cms.InputTag("scPhase2PuppiRawToDigiStruct"),
    minPtDau = cms.double(5.0),
    maxDeltaRDaus = cms.double(0.40 * 0.40),
    minDeltaR = cms.double(0.05 * 0.05),
    maxDeltaR = cms.double(0.25 * 0.25),    
    mesonType = cms.string("jpsi"),
    runStruct = cms.bool(True)
)

hjpsigammaRecmesonStruct = cms.EDProducer("ScPhase2RecMesonHJPsiGamma",
    srcMeson = cms.InputTag("jpsiRecmesonStruct"),
    srcGamma = cms.InputTag("photonIsolationStruct"),
)

phiRecmesonStruct = cms.EDProducer("ScPhase2RecMeson",
    src = cms.InputTag("scPhase2PuppiRawToDigiStruct"),
    minPtDau = cms.double(5.0),
    maxDeltaRDaus = cms.double(0.40 * 0.40),
    minDeltaR = cms.double(0.05 * 0.05),
    maxDeltaR = cms.double(0.25 * 0.25),
    mesonType = cms.string("phi"),
    runStruct = cms.bool(True)
)

h2phiRecmesonStruct = cms.EDProducer("ScPhase2RecMesonH2Phi",
     src = cms.InputTag("phiRecmesonStruct"),
)

hphigammaRecmesonStruct = cms.EDProducer("ScPhase2RecMesonHPhiGamma",
    srcMeson = cms.InputTag("phiRecmesonStruct"),
    srcGamma = cms.InputTag("photonIsolationStruct"),
)

hphijpsiRecmesonStruct = cms.EDProducer("ScPhase2RecMesonHPhiJPsi",
    srcPhi = cms.InputTag("phiRecmesonStruct"),
    srcJPsi = cms.InputTag("jpsiRecmesonStruct"),
)

rhoRecmesonStruct = cms.EDProducer("ScPhase2RecMeson",
    src = cms.InputTag("scPhase2PuppiRawToDigiStruct"),
    minPtDau = cms.double(5.0),
    maxDeltaRDaus = cms.double(0.40 * 0.40),
    minDeltaR = cms.double(0.05 * 0.05),
    maxDeltaR = cms.double(0.25 * 0.25),
    mesonType = cms.string("rho"),
    runStruct = cms.bool(True)
)

h2rhoRecmesonStruct = cms.EDProducer("ScPhase2RecMesonH2Rho",
     src = cms.InputTag("rhoRecmesonStruct"),
)

hrhogammaRecmesonStruct = cms.EDProducer("ScPhase2RecMesonHRhoGamma",
    srcMeson = cms.InputTag("rhoRecmesonStruct"),
    srcGamma = cms.InputTag("photonIsolationStruct"),
)

z2phiRecmesonStruct = cms.EDProducer("ScPhase2RecMesonZ2Phi",
     src = cms.InputTag("phiRecmesonStruct"),
)

z2rhoRecmesonStruct = cms.EDProducer("ScPhase2RecMesonZ2Rho",
     src = cms.InputTag("rhoRecmesonStruct"),
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