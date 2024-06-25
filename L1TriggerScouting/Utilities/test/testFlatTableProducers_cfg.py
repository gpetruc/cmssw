import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing ('analysis')

options.register ('numOrbits',
                  -1,
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.int,
                  "Number of orbits to process")

options.register ('inFile',
                  "file:/dev/shm/PoolOutputTest.root",
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.string,
                  "Sub lumisection number to process")

options.register ('outFile',
                  "NanoOutput.root",
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.string,
                  "Sub lumisection number to process")


options.parseArguments()

process = cms.Process( "DUMP" )


process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(options.numOrbits)
)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(options.inFile)
)

process.scMuonTable = cms.EDProducer("ConvertScoutingMuonsToOrbitFlatTable",
  src = cms.InputTag("l1ScGmtUnpacker", "Muon"),
  name = cms.string("Muon"),
  doc = cms.string("Muons from GMT"),
)
process.scJetTable = cms.EDProducer("ConvertScoutingJetsToOrbitFlatTable",
  src = cms.InputTag("l1ScCaloUnpacker", "Jet"),
  name = cms.string("Jet"),
  doc = cms.string("Jets from Calo Demux"),
)
process.scEgammaTable = cms.EDProducer("ConvertScoutingEGammasToOrbitFlatTable",
  src = cms.InputTag("l1ScCaloUnpacker", "EGamma"),
  name = cms.string("EGamma"),
  doc = cms.string("EGammas from Calo Demux"),
)
process.scTauTable = cms.EDProducer("ConvertScoutingTausToOrbitFlatTable",
  src = cms.InputTag("l1ScCaloUnpacker", "Tau"),
  name = cms.string("Tau"),
  doc = cms.string("Taus from Calo Demux"),
)
process.p = cms.Path(
  process.scMuonTable +
  process.scJetTable +
  process.scEgammaTable +
  process.scTauTable
)
process.out = cms.OutputModule("OrbitNanoAODOutputModule",
    fileName = cms.untracked.string(options.outFile),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('p')),
    outputCommands = cms.untracked.vstring("drop *", "keep l1ScoutingRun3OrbitFlatTable_*_*_*"),
    compressionLevel = cms.untracked.int32(4),
    compressionAlgorithm = cms.untracked.string("LZ4"),
)
process.o = cms.EndPath(
  process.out
)