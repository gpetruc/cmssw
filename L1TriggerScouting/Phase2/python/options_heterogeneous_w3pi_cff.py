import FWCore.ParameterSet.VarParsing as VarParsing


args = VarParsing.VarParsing("analysis")

args.register(
    "numberOfThreads",
    1,
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.int,
    "Number of CMSSW threads"
)

args.register(
    "numberOfStreams",
    1,
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.int,
    "Number of CMSSW streams"
)

args.register(
    "numberOfEvents",
    1,
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.int,
    "Number of events to process"
)

args.register(
    "backend",
    "serial_sync",
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.string,         
    "Hardware accelerator backend: serial_sync, cuda_async, or rocm_async"
)

args.register(
    "environment",
    0, 
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.int,         
    "0 - production, 1 - development, 2 - test"
)

args.register(
    "fastPath",
    True,
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.bool,
    "Use fast path, skip full reconstruction of candidate table and compression of bx selection mask"
)

args.register (
    "runNumber",
    38,
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.int,
    "Run number"
)

args.register (
    "lumiNumber",
    1,
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.int,
    "Lumisection number"
)

args.register (
    "daqSourceMode",
    "ScoutingPhase2",
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.string,
    "DAQ source data mode")

args.register (
    "broker",
    "none",
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.string,
    "Broker: 'none' or 'hostname:port'"
)

args.register (
    "fuBaseDir",
    "/dev/shm/ramdisk",
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.string,
    "BU base directory"
)

args.register (
    "buBaseDir",
    "/dev/shm/ramdisk",
    VarParsing.VarParsing.multiplicity.list,
    VarParsing.VarParsing.varType.string,
    "BU base directory"
)

args.register (
    "buNumStreams",
    [],
    VarParsing.VarParsing.multiplicity.list,
    VarParsing.VarParsing.varType.int,
    "Number of input streams (i.e. files) used simultaneously for each BU directory"
)


args.register (
    "streams",
    [],
    VarParsing.VarParsing.multiplicity.list,
    VarParsing.VarParsing.varType.int,
    "Input links IDs for the inputs"
)

args.register (
    "name",
    "",
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.string,
    "Name for output report file"
)

# W3Pi algorithm parameters
args.register(
    "ptMin",
    7,
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.float,
    "pT min"
)

args.register(
    "ptInt",
    12,
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.float,
    "pT interval"
)

args.register(
    "ptMax",
    15,
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.float,
    "pT max"
)

args.register(
    "invariantMassLB",
    40.0,
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.float,
    "invariant_mass_lower_bound"
)

args.register(
    "invariantMassUB",
    150.0,
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.float,
    "invariant_mass_upper_bound"
)

args.register(
    "minDeltaR",
    0.01 * 0.01,
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.float,
    "min_deltar_threshold"
)

args.register(
    "maxDeltaR",
    0.25 * 0.25,
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.float,
    "max_deltar_threshold"
)

args.register(
    "maxIso",
    2.0,
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.float,
    "max_isolation_threshold"
)

args.register(
    "angSepLB",
    0.5 * 0.5,
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.float,
    "ang_sep_lower_bound"
)