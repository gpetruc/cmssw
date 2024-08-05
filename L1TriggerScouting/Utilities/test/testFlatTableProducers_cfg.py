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


options.register ('selBx',
                  "none",
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.string,
                 "BX Selection to use")

options.register ('saveStubs',
                  True,
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.bool,
                 "BX Selection to use")

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

selbx = options.selBx if options.selBx != "none" else None

process.scMuonTable = cms.EDProducer("ConvertScoutingMuonsToOrbitFlatTable",
  src = cms.InputTag("FinalBxSelectorMuon" if selbx else "l1ScGmtUnpacker", "Muon"),
  name = cms.string("L1Mu"),
  doc = cms.string("Muons from GMT"),
)
process.scJetTable = cms.EDProducer("ConvertScoutingJetsToOrbitFlatTable",
  src = cms.InputTag("FinalBxSelectorJet" if selbx else "l1ScCaloUnpacker", "Jet"),
  name = cms.string("L1Jet"),
  doc = cms.string("Jets from Calo Demux"),
)
process.scEgammaTable = cms.EDProducer("ConvertScoutingEGammasToOrbitFlatTable",
  src = cms.InputTag("FinalBxSelectorEGamma" if selbx else "l1ScCaloUnpacker", "EGamma"),
  name = cms.string("L1EG"),
  doc = cms.string("EGammas from Calo Demux"),
)
process.scTauTable = cms.EDProducer("ConvertScoutingTausToOrbitFlatTable",
  src = cms.InputTag("l1ScCaloUnpacker", "Tau"),
  name = cms.string("L1Tau"),
  doc = cms.string("Taus from Calo Demux"),
)
process.scStubsTable = cms.EDProducer("ConvertScoutingStubsToOrbitFlatTable",
  src = cms.InputTag("FinalBxSelectorBMTFStub" if selbx else "l1ScBMTFUnpacker", "BMTFStub"),
  name = cms.string("L1BMTFStub"),
  doc = cms.string("Stubs from BMTF"),
)
process.scSumTable = cms.EDProducer("ConvertScoutingSumsToOrbitFlatTable",
  src = cms.InputTag("FinalBxSelectorBxSums" if selbx else "l1ScCaloUnpacker", "EtSum"),
  name = cms.string("L1EtSum"),
  doc = cms.string("Sums from Calo Demux"),
  singleObject = cms.bool(False),
  writeHF = cms.bool(True),
  writeMinBias = cms.bool(False),
  writeCentrality = cms.bool(False),
  writeAsym = cms.bool(False),
)
process.p = cms.Path(
  process.scMuonTable +
  process.scJetTable +
  process.scEgammaTable +
  process.scTauTable +
  process.scStubsTable +
  process.scSumTable
)

process.out = cms.OutputModule("OrbitNanoAODOutputModule",
    fileName = cms.untracked.string(options.outFile),
    skipEmptyBXs = cms.bool(True),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('p')),
    outputCommands = cms.untracked.vstring("drop *", "keep l1ScoutingRun3OrbitFlatTable_*_*_*"),
    compressionLevel = cms.untracked.int32(5),
    compressionAlgorithm = cms.untracked.string("ZSTD"),
)

if not options.saveStubs:
  process.p.remove(process.scStubsTable)
if selbx:
  process.p.remove(process.scTauTable)
  process.out.outputCommands += [ "keep uints_*_SelBx_*" ]
  if selbx != "any":
    process.out.selectedBx = cms.InputTag(selbx, "SelBx")
    process.out.skipEmptyBXs = False
  else:
    process.out.selectedBx = cms.InputTag("FinalBxSelector", "SelBx")
    process.out.skipEmptyBXs = False

process.o = cms.EndPath(
  process.out
)