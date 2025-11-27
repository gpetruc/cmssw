import argparse

def parse_args():
    parser = argparse.ArgumentParser()

    # Basic CMSSW settings
    parser.add_argument(
        "-nt", "--numberOfThreads",
        type=int,
        default=1,
        help="Number of CMSSW threads"
    )
    parser.add_argument(
        "-ns", "--numberOfStreams",
        type=int,
        default=1,
        help="Number of CMSSW streams"
    )
    parser.add_argument(
        "-ne", "--numberOfEvents",
        type=int,
        default=1,
        help="Number of events to process"
    )
    
    # pipeline 
    parser.add_argument(
        "-o", "--only",
        nargs="+",
        default=["tagging_inf"],
        choices=["unpacking", "clustering", "tagging_pre", "tagging_inf"],
        help="Run only the specified pipeline stages"
    )

    # dump
    parser.add_argument(
        "-d", "--dump",
        nargs="+",
        default=["none"],
        choices=["none", "candidates", "clusters", "inference", "all_no_cands", "all"],
        help="Select which artifacts to dump on NanoAOD file"
    )

    # Backend and environment
    parser.add_argument(
        "-b", "--backend",
        type=str,
        default="serial_sync",
        choices=["serial_sync", "cuda_async", "rocm_async"],
        help="Hardware accelerator backend"
    )
    parser.add_argument(
        "-e", "--environment",
        type=int,
        default=0,
        choices=[0, 1, 2, 3],
        help="0 - production, 1 - development, 2 - test, 3 - debug"
    )
    parser.add_argument(
        "-ws", "--wantSummary",
        action='store_true',
        help="Show timing report"
    )

    # Clustering parameters
    parser.add_argument(
        "--dc",
        type=float,
        default=0.2,
        help="Side of the box inside which the density of a point is calculated"
    )
    parser.add_argument(
        "--rhoc",
        type=float,
        default=5.0,
        help="Minimum rhoc required for a point to be considered a seed candidate"
    )
    parser.add_argument(
        "--dm",
        type=float,
        default=0.4,
        help="Side of the box inside which the followers of a point are searched"
    )
    parser.add_argument(
        "-wc", "--wrapCoords",
        action="store_true",
        help="Wrap phi coordinate in CLUEstering"
    )

    # Scouting configuration
    parser.add_argument(
        "-scout", "--runScouting",
        action="store_true",
        help="Run scouting-based tagging"
    )
    parser.add_argument(
        "-rn","--runNumber",
        type=int,
        default=38,
        help="Run number"
    )
    parser.add_argument(
        "-ln", "--lumiNumber",
        type=int,
        default=1,
        help="Lumisection number"
    )
    parser.add_argument(
        "-dsm", "--daqSourceMode",
        type=str,
        default="ScoutingPhase2",
        help="DAQ source data mode"
    )
    parser.add_argument(
        "-broker", "--broker",
        type=str,
        default="none",
        help="Broker: 'none' or 'hostname:port'"
    )

    # Tagger
    parser.add_argument(
        "-m","--model",
        type=str,
        default="L1TriggerScouting/TauTagging/data/softtauid_sigmoid.pt",
        help="Path to JIT compiled PyTorch model."
    )

    # Directories and I/O streams
    parser.add_argument(
        "-fbd", "--fuBaseDir",
        type=str,
        default="/dev/shm/ramdisk",
        help="FU base directory"
    )
    parser.add_argument(
        "-bbd", "--buBaseDir",
        nargs="+",
        default=["/dev/shm/ramdisk"],
        help="BU base directory (can specify multiple)"
    )
    parser.add_argument(
        "-bns", "--buNumStreams",
        nargs="+",
        type=int,
        default=[],
        help="Number of input streams (i.e. files) used simultaneously for each BU directory"
    )
    parser.add_argument(
        "-s", "--streams",
        nargs="+",
        type=int,
        default=[],
        help="Input link IDs for the inputs"
    )

    # Output/reporting
    parser.add_argument(
        "-n", "--name",
        type=str,
        default="",
        help="Name for output report file"
    )

    return parser.parse_args()
