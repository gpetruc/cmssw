import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
process = cms.Process("TestUnpacking", eras.Phase2C17I13M9)

from PhysicsTools.NanoAOD.common_cff import Var, ExtVar
def LazyVar(expr, valtype, doc=None, precision=-1):
    return Var(expr, valtype, doc, precision, lazyEval=True)


# logging configuration
process.load('Configuration.StandardSequences.Services_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Accelerators_cff")
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.options.wantSummary = cms.untracked.bool(True)
process.maxEvents.input = cms.untracked.int32(100)

process.load('Configuration.Geometry.GeometryExtendedRun4D110Reco_cff')
process.load('Configuration.Geometry.GeometryExtendedRun4D110_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff') # needed to read HCal TPs
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '141X_mcRun4_realistic_v3', '')

process.load('L1Trigger.Phase2L1ParticleFlow.l1ctLayer1_cff')
process.load('L1Trigger.L1TTrackMatch.l1tGTTInputProducer_cfi')
process.load('L1Trigger.L1TTrackMatch.l1tTrackSelectionProducer_cfi')
process.l1tTrackSelectionProducer.processSimulatedTracks = False # these would need stubs, and are not used anyway
process.load('L1Trigger.VertexFinder.l1tVertexProducer_cfi')
from L1Trigger.Configuration.SimL1Emulator_cff import l1tSAMuonsGmt
process.l1tSAMuonsGmt = l1tSAMuonsGmt.clone()
from L1Trigger.L1CaloTrigger.l1tPhase2L1CaloEGammaEmulator_cfi import l1tPhase2L1CaloEGammaEmulator
process.l1tPhase2L1CaloEGammaEmulator = l1tPhase2L1CaloEGammaEmulator.clone()
from L1Trigger.L1CaloTrigger.l1tPhase2CaloPFClusterEmulator_cfi import l1tPhase2CaloPFClusterEmulator
process.l1tPhase2CaloPFClusterEmulator = l1tPhase2CaloPFClusterEmulator.clone()
from L1Trigger.L1CaloTrigger.l1tPhase2GCTBarrelToCorrelatorLayer1Emulator_cfi import l1tPhase2GCTBarrelToCorrelatorLayer1Emulator
process.l1tPhase2GCTBarrelToCorrelatorLayer1Emulator = l1tPhase2GCTBarrelToCorrelatorLayer1Emulator.clone()

process.L1TInputTask = cms.Task(
    process.l1tSAMuonsGmt,
    process.l1tPhase2L1CaloEGammaEmulator,
    process.l1tPhase2CaloPFClusterEmulator,
    process.l1tPhase2GCTBarrelToCorrelatorLayer1Emulator
)

from IOPool.Input.modules import PoolSource
process.source = PoolSource(
    fileNames = [
        #"file:" + cms.FileInPath("L1TriggerScouting/TauTagging/data/pfOnly.root").value()
        #'file:/afs/cern.ch/work/g/gpetrucc/l1p2/l1scout-cmssw/CMSSW_16_0_0_pre1/src/L1TriggerScouting/TauTagging/data/pfOnly.root'
        '/store/cmst3/group/l1tr/vcamagni/L1TauID/DATA/FPinputs/m40/4STEPS/142Xv0/inputs140X_7099346_0.root'
    ]
)

process.runPF = cms.Path( 
        process.l1tGTTInputProducer +
        process.l1tTrackSelectionProducer +
        process.l1tVertexFinderEmulator +
        process.l1tLayer1BarrelExtended + 
        process.l1tLayer1HGCalExtended +
        process.l1tLayer1HGCalNoTK +
        process.l1tLayer1HF +
        process.l1tLayer1Extended 
)
process.runPF.associate(process.L1TLayer1TaskInputsTask)
process.runPF.associate(process.L1TInputTask)

process.pfBarrelTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
                name = cms.string("L1PFBarrel"),
                src = cms.InputTag("l1tLayer1BarrelExtended:PF"),
                cut = cms.string(""),
                doc = cms.string(""),
                singleton = cms.bool(False), # the number of entries is variable
                extension = cms.bool(False), # this is the main table
                variables = cms.PSet(
                    pt  = Var("pt",  float, precision=16),
                    eta  = Var("eta", float, precision=16),
                    phi = Var("phi", float, precision=16),
                    pdgId = Var("pdgId", int,),
                    z0 = Var("vz", float, precision=16),
                    dxy = LazyVar("dxy", float, precision=16),
                    quality = LazyVar("hwQual", int),
                    puppiw = LazyVar("puppiWeight", float, precision=16),
                )
)
process.pfEndcapTable = process.pfBarrelTable.clone(
                name = cms.string("L1PFEndcap"),
                src = cms.InputTag("l1tLayer1HGCalExtended:PF"),
)
process.puppiTable = process.pfBarrelTable.clone(
                name = cms.string("L1Puppi"),
                src = cms.InputTag("l1tLayer1Extended:Puppi"),
)


process.packPFBarrel = cms.EDProducer("ScPhase2PuppiPacker",
    src = cms.InputTag("l1tLayer1BarrelExtended:PF"),
    fedIDs = cms.vuint32(1,2,3),
    splitFactor = cms.uint32(3),
    scoutingHeader = cms.bool(True)
)

process.packPFEndcap = cms.EDProducer("ScPhase2PuppiPacker",
    src = cms.InputTag("l1tLayer1HGCalExtended:PF"),
    fedIDs = cms.vuint32(4,5),
    splitFactor = cms.uint32(2),
    scoutingHeader = cms.bool(True)
)

process.packPuppi = cms.EDProducer("ScPhase2PuppiPacker",
    src = cms.InputTag("l1tLayer1Extended:Puppi"),
    fedIDs = cms.vuint32(0),
    splitFactor = cms.uint32(1),
    scoutingHeader = cms.bool(True)
)

from L1TriggerScouting.Phase2.modules import l1sc_L1TScPhase2PuppiRawToDigi_alpaka
process.unpackPFBarrelAlpaka = l1sc_L1TScPhase2PuppiRawToDigi_alpaka(
    #alpaka = cms.untracked.PSet( backend = cms.untracked.string("serial_sync") ),
    streams = process.packPFBarrel.fedIDs,
    splitFactor = process.packPFBarrel.splitFactor,
    src = cms.InputTag('packPFBarrel'),
)
process.unpackPFEndcapAlpaka = l1sc_L1TScPhase2PuppiRawToDigi_alpaka(
    #alpaka = cms.untracked.PSet( backend = cms.untracked.string("serial_sync") ),
    streams = process.packPFEndcap.fedIDs,
    splitFactor = process.packPFEndcap.splitFactor,
    src = cms.InputTag('packPFEndcap'),
)

process.unpackPFBarrel = cms.EDProducer('ScPhase2PuppiRawToDigi',
    fedIDs = process.packPFBarrel.fedIDs,
    splitFactor = process.packPFBarrel.splitFactor,
    src = cms.InputTag('packPFBarrel'),
)
process.unpackPFEndcap = cms.EDProducer('ScPhase2PuppiRawToDigi',
    fedIDs = process.packPFEndcap.fedIDs,
    splitFactor = process.packPFEndcap.splitFactor,
    src = cms.InputTag('packPFEndcap'),
)
process.unpackPuppi = cms.EDProducer('ScPhase2PuppiRawToDigi',
    fedIDs = process.packPuppi.fedIDs,
    splitFactor = process.packPuppi.splitFactor,
    src = cms.InputTag('packPuppi'),
)
process.scPFBarrelStructToTable = cms.EDProducer("ScPuppiToOrbitFlatTable",
    src = cms.InputTag("unpackPFBarrel"),
    name = cms.string("L1PFUnpackBarrel"),
    doc = cms.string("L1PF candidates from Barrel, unpacked to OrbitCollection"),
)
process.scPFEndcapStructToTable = cms.EDProducer("ScPuppiToOrbitFlatTable",
    src = cms.InputTag("unpackPFEndcap"),
    name = cms.string("L1PFUnpackEndcap"),
    doc = cms.string("L1PF candidates from Endcap, unpacked to OrbitCollection"),
)
process.scPFEndcapStructToTable = cms.EDProducer("ScPuppiToOrbitFlatTable",
    src = cms.InputTag("unpackPFEndcap"),
    name = cms.string("L1PFUnpackEndcap"),
    doc = cms.string("L1PF candidates from Endcap, unpacked to OrbitCollection"),
)
process.scPuppiStructToTable = cms.EDProducer("ScPuppiToOrbitFlatTable",
    src = cms.InputTag("unpackPuppi"),
    name = cms.string("L1PuppiUnpack"),
    doc = cms.string("L1Puppi candidates, unpacked to OrbitCollection"),
)
process.scPFBarrelSoAToTable = cms.EDProducer("PFSoAToOrbitFlatTable",
    srcBx = cms.InputTag("unpackPFBarrelAlpaka"),
    srcPF = cms.InputTag("unpackPFBarrelAlpaka"),
    name = cms.string("L1PFUnpackAlpakaBarrel"),
    doc = cms.string("L1PF candidates from Barrel, unpacked to Alpaka SoA"),
)
process.scPFEndcapSoAToTable = cms.EDProducer("PFSoAToOrbitFlatTable",
    srcBx = cms.InputTag("unpackPFEndcapAlpaka"),
    srcPF = cms.InputTag("unpackPFEndcapAlpaka"),
    name = cms.string("L1PFUnpackAlpakaEndcap"),
    doc = cms.string("L1PF candidates from Endcap, unpacked to Alpaka SoA"),
)

process.p = cms.Path(
    process.packPFBarrel + process.packPFEndcap + process.packPuppi +
    process.unpackPFBarrel + process.unpackPFEndcap + process.unpackPuppi +
    #process.unpackPFBarrelAlpaka + process.unpackPFEndcapAlpaka +
    process.pfBarrelTable + process.pfEndcapTable + process.puppiTable +
    process.scPFBarrelStructToTable + process.scPFEndcapStructToTable + process.scPuppiStructToTable #+
    #process.scPFBarrelSoAToTable + process.scPFEndcapSoAToTable
)
process.outStandard = cms.OutputModule("NanoAODOutputModule",
    fileName = cms.untracked.string("plainNano.root"),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring()),
    outputCommands = cms.untracked.vstring("drop *", "keep nanoaodFlatTable_*Table_*_*"),
    compressionLevel = cms.untracked.int32(4),
    compressionAlgorithm = cms.untracked.string("ZLIB"),
)
process.outScout = cms.OutputModule("OrbitNanoAODOutputModule",
    fileName = cms.untracked.string("scoutingNano.root"),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring()),
    skipEmptyBXs = cms.bool(True),
    outputCommands = cms.untracked.vstring("drop *", "keep l1ScoutingRun3OrbitFlatTable_*_*_*"),
    compressionLevel = cms.untracked.int32(4),
    compressionAlgorithm = cms.untracked.string("ZLIB"),
)
process.e = cms.EndPath(process.outStandard + process.outScout)