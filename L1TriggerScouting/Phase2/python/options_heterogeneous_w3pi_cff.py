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
    "verbose",
    False,
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.bool,         
    "Log debug messages to stdout"
)

args.register(
    "verboseLevel",
    0,
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.int,         
    "0 - debug, 1 - full"
)

args.register (
    "runNumber",
    37,
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
    "linksIds",
    [],
    VarParsing.VarParsing.multiplicity.list,
    VarParsing.VarParsing.varType.int,
    "Input links IDs for the inputs"
)


args.register (
    "fuBaseDir",
    "/dev/shm/data",
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.string,
    "BU base directory"
)