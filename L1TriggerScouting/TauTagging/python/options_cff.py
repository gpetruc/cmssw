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
    VarParsing.VarParsing.varType.bool,         
    "0 - production, 1 - development, 2 - test"
)

args.register(
    "dc",
    0.2,
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.float,         
    "Side of the box inside of which the density of a point is calculated"
)

args.register(
    "rhoc",
    5.0,
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.float,         
    "Minimum rhoc required for a point to be considered a seed candidate "
)

args.register(
    "dm",
    0.4,
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.float,         
    "Side of the box inside of which the followers of a point are searched"
)

args.register(
    "wrapCoords",
    False,
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.bool,         
    "Wrap phi coordinate in CLUEstering"
)

# scouting
args.register(
    "runScouting",
    False,
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.bool,         
    "Run scouting based tagging"
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