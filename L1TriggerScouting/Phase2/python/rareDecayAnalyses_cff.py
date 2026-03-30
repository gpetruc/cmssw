import FWCore.ParameterSet.Config as cms


hjpsigammaRecMesonStruct = cms.EDProducer("ScPhase2BosonToRecMesonGamma",
    srcMeson = cms.InputTag("puppiRecMesonStruct", "jpsi"),
    srcGamma = cms.InputTag("recIsoTkEmStruct"),
    minMassBoson = cms.double(100),
    maxMassBoson = cms.double(150),
    minPtQ = cms.double(30),
    minPtGamma = cms.double(30),
    maxRelIsoQ = cms.double(0.25),
    analysisName = cms.string("HJPsiGamma")
)
hjpsigammaTTrackRecMesonStruct = hjpsigammaRecMesonStruct.clone(
    srcMeson = "ttrackRecMesonStruct:jpsi",
    analysisName = "HJPsiGammaTTrack",
)
hjpsigammaMuMuRecMesonStruct = hjpsigammaRecMesonStruct.clone(
    srcMeson = "puppiRecMesonStruct:jpsi",
    analysisName = "HJPsiGammaMuMu",
)
hjpsigammaEERecMesonStruct = hjpsigammaRecMesonStruct.clone(
    srcMeson = "ttrackRecMesonStruct:jpsi",
    analysisName = "HJPsiGammaEE",
)

h2phiRecMesonStruct = cms.EDProducer("ScPhase2BosonTo2RecMeson",
    srcMeson1 = cms.InputTag("puppiRecMesonStruct", "phi"),
    srcMeson2 = cms.InputTag("puppiRecMesonStruct", "phi"),
    minMassBoson = cms.double(100),
    maxMassBoson = cms.double(150),
    minPtQ = cms.double(1),
    maxIso = cms.double(0.25),
    analysisName = cms.string("H2Phi")
)
h2phiTTrackRecMesonStruct = h2phiRecMesonStruct.clone(
    srcMeson1 = "ttrackRecMesonStruct:phi",
    srcMeson2 = "ttrackRecMesonStruct:phi",
    analysisName = "H2PhiTTrack",
)   

hphigammaRecMesonStruct = cms.EDProducer("ScPhase2BosonToRecMesonGamma",
    srcMeson = cms.InputTag("puppiRecMesonStruct", "phi"),
    srcGamma = cms.InputTag("recIsoTkEmStruct"),
    minMassBoson = cms.double(100),
    maxMassBoson = cms.double(150),
    minPtQ = cms.double(30),
    minPtGamma = cms.double(30),
    maxRelIsoQ = cms.double(0.25),
    analysisName = cms.string("HPhiGamma")
)
hphigammaTTrackRecMesonStruct = hphigammaRecMesonStruct.clone(
    srcMeson = "ttrackRecMesonStruct:phi",
    analysisName = "HPhiGammaTTrack",
)


hphijpsiRecMesonStruct = cms.EDProducer("ScPhase2BosonTo2RecMeson",
    srcMeson1 = cms.InputTag("puppiRecMesonStruct", "phi"),
    srcMeson2 = cms.InputTag("puppiRecMesonStruct", "jpsi"),
    minMassBoson = cms.double(100),
    maxMassBoson = cms.double(150),
    minPtQ = cms.double(30),
    maxIso = cms.double(0.25),
    analysisName = cms.string("HPhiJPsi")
)
hphijpsiTTrackRecMesonStruct = hphijpsiRecMesonStruct.clone(
    srcMeson1 = "ttrackRecMesonStruct:phi",
    srcMeson2 = "ttrackRecMesonStruct:jpsi",
    analysisName = "HPhiJPsiTTrack",
)
hphijpsiMuMuRecMesonStruct = hphijpsiRecMesonStruct.clone(
    srcMeson1 = "puppiRecMesonStruct:phi",
    srcMeson2 = "tkMuonRecMesonStruct:jpsi",
    analysisName = "HPhiJPsiMuMu",
)
hphijpsiEERecMesonStruct = hphijpsiRecMesonStruct.clone(
    srcMeson1 = "puppiRecMesonStruct:phi",
    srcMeson2 = "tkEleRecMesonStruct:jpsi",
    analysisName = "HPhiJPsiEE",
)


h2rhoRecMesonStruct = cms.EDProducer("ScPhase2BosonTo2RecMeson",
    srcMeson1 = cms.InputTag("puppiRecMesonStruct", "rho"),
    srcMeson2 = cms.InputTag("puppiRecMesonStruct", "rho"),
    minMassBoson = cms.double(100),
    maxMassBoson = cms.double(150),
    minPtQ = cms.double(3),
    maxIso = cms.double(0.25),
    analysisName = cms.string("H2Rho")
)
h2rhoTTrackRecMesonStruct = h2rhoRecMesonStruct.clone(
    srcMeson1 = "ttrackRecMesonStruct:rho",
    srcMeson2 = "ttrackRecMesonStruct:rho",
    analysisName = "H2RhoTTrack",
)

hrhogammaRecMesonStruct = cms.EDProducer("ScPhase2BosonToRecMesonGamma",
    srcMeson = cms.InputTag("puppiRecMesonStruct", "rho"),
    srcGamma = cms.InputTag("recIsoTkEmStruct"),
    minMassBoson = cms.double(100),
    maxMassBoson = cms.double(150),
    minPtQ = cms.double(30),
    minPtGamma = cms.double(30),
    maxRelIsoQ = cms.double(0.25),
    analysisName = cms.string("HRhoGamma")
)
hrhogammaTTrackRecMesonStruct = hrhogammaRecMesonStruct.clone(
    srcMeson = "ttrackRecMesonStruct:rho",
    analysisName = "HRhoGammaTTrack",
)


z2phiRecMesonStruct = cms.EDProducer("ScPhase2BosonTo2RecMeson",
    srcMeson1 = cms.InputTag("puppiRecMesonStruct", "phi"),
    srcMeson2 = cms.InputTag("puppiRecMesonStruct", "phi"),
    minMassBoson = cms.double(60),
    maxMassBoson = cms.double(120),
    minPtQ = cms.double(1),
    maxIso = cms.double(0.25),
    analysisName = cms.string("Z2Phi")
)
z2phiTTrackRecMesonStruct = z2phiRecMesonStruct.clone(
    srcMeson1 = "ttrackRecMesonStruct:phi",
    srcMeson2 = "ttrackRecMesonStruct:phi",
    analysisName = "Z2PhiTTrack",
)


z2rhoRecMesonStruct = cms.EDProducer("ScPhase2BosonTo2RecMeson",
    srcMeson1 = cms.InputTag("puppiRecMesonStruct", "rho"),
    srcMeson2 = cms.InputTag("puppiRecMesonStruct", "rho"),
    minMassBoson = cms.double(60),
    maxMassBoson = cms.double(120),
    minPtQ = cms.double(3),
    maxIso = cms.double(0.25),
    analysisName = cms.string("Z2Rho")
)
z2rhoTTrackRecMesonStruct = z2rhoRecMesonStruct.clone(
    srcMeson1 = "ttrackRecMesonStruct:rho",
    srcMeson2 = "ttrackRecMesonStruct:rho",
    analysisName = "Z2RhoTTrack",
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