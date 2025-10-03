import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing ('analysis')
options.register ('runNumber',
                  37,
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.int,          # string, int, or float
                  "Run Number")

options.register ('lumiNumber',
                  1,
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.int,          # string, int, or float
                  "Lumisection Number (only for filelist mode)")

options.register ('daqSourceMode',
                  'ScoutingPhase2', # default value
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "DAQ source data mode")

options.register ('broker',
                  'none', # default value
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "Broker: 'none' or 'hostname:port'")

options.register ('buBaseDir',
                  '/dev/shm/ramdisk', # default value
                  VarParsing.VarParsing.multiplicity.list,
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "BU base directory")

options.register ('buNumStreams',
                  [], # default value
                  VarParsing.VarParsing.multiplicity.list,
                  VarParsing.VarParsing.varType.int,          # string, int, or float
                  "Number of input streams (i.e. files) used simultaneously for each BU directory")

options.register ('timeslices',
                  1,
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.int,          # string, int, or float
                  "Number of timeslices")

options.register ('tmuxPeriod',
                  1,
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.int,          # string, int, or float
                  "Time multiplex period")

options.register ('puppiStreamIDs',
                  [], # default value
                  VarParsing.VarParsing.multiplicity.list,
                  VarParsing.VarParsing.varType.int,          # string, int, or float
                  "Stream IDs for the Puppi inputs")

options.register ('tkEmStreamIDs',
                  [], # default value
                  VarParsing.VarParsing.multiplicity.list,
                  VarParsing.VarParsing.varType.int,          # string, int, or float
                  "Stream IDs for the TkEm inputs")

options.register ('tkMuStreamIDs',
                  [], # default value
                  VarParsing.VarParsing.multiplicity.list,
                  VarParsing.VarParsing.varType.int,          # string, int, or float
                  "Stream IDs for the TkMu inputs")

options.register ('pfBarrelStreamIDs',
                  [], # default value
                  VarParsing.VarParsing.multiplicity.list,
                  VarParsing.VarParsing.varType.int,          # string, int, or float
                  "Stream IDs for the PF Barrel inputs")

options.register ('pfEndcapStreamIDs',
                  [], # default value
                  VarParsing.VarParsing.multiplicity.list,
                  VarParsing.VarParsing.varType.int,          # string, int, or float
                  "Stream IDs for the PF Endcap inputs")

options.register ('fuBaseDir',
                  '/dev/shm/data', # default value
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "BU base directory")

options.register ('fffBaseDir',
                  '/dev/shm', # default value
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "FFF base directory")

options.register ('numThreads',
                  1, # default value
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.int,          # string, int, or float
                  "Number of CMSSW threads")

options.register ('numFwkStreams',
                  1, # default value
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.int,          # string, int, or float
                  "Number of CMSSW streams")

options.register ('run',
                  'both', # default value
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "'inclusive', 'selected', 'both' (default).")

options.register ('analyses',
                  [], # default value
                  VarParsing.VarParsing.multiplicity.list,
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "analyses: any list of w3pi, wdsg, wpig, hrhog, hphig, hjpsig, h2rho, h2phi, dimu, zdee")

options.register ('prescaleInclusive',
                  100, # default value
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.int,          # string, int, or float
                  "Prescale factor for the inclusive stream.")

options.register ('outMode',
                  'none', # default value
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "output (none, nanoSelected, nanoInclusive, nanoBoth)")
                   
options.register ('outFile',
                  "NanoOutput.root",
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.string,
                  "Output file name")

options.register ('task',
                  0,
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.int,          # string, int, or float
                  "Task index (used for json outputs)")

options.register ('backend',
                  "serial_sync",
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.string,         
                  "Hardware accelerator backend: serial_sync, cuda_async, or rocm_async")

options.register ('environment',
                  0, 
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.int,         
                  '0 - production, 1 - development, 2 - test'
)

options.register ('fastPath',
                  False,
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.bool,
                  'Use fast path, skip full reconstruction of candidate table and compression of bx selection mask in w3pi portable implementation'
)
