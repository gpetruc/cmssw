#!/bin/bash
# to streamline debugging and testing of W3Pi
# usage: bash L1TriggerScouting/Phase2/test/runScoutingPhase2HeterogeneousW3Pi.sh <backend> <num_events> <run_number>
# example usage: bash L1TriggerScouting/Phase2/test/runScoutingPhase2HeterogeneousW3Pi.sh cuda 10 37

function die { echo Failed $1: status $2 ; exit $2 ; }

SCRIPT="L1TriggerScouting/Phase2/test/runScoutingPhase2HeterogeneousW3Pi.py"
DATA=/eos/home-l/lmichals/W3Pi/data/raw/

if [ "$#" != "3" ]; then
    die "Need exactly 3 arguments: 1st ('cpu', 'cuda', or 'rocm'), 2nd ('num_events'), 3rd ('run_number'), got $#" 1
fi
if [[ "$1" =~ ^(cpu|cuda|rocm)$ ]]; then
    TARGET=$1
else
    die "Argument needs to be 'cpu', 'cuda', or 'rocm'; got '$1'" 1
fi

if [ "${TARGET}" == "cpu" ]; then
  echo "Running CPU-only test"
  cmsRun "${SCRIPT}" runNumber=$3 buBaseDir=${DATA} fuBaseDir=${DATA} buNumStreams=6 numberOfEvents=$2 backend=serial_sync
elif [ "${TARGET}" == "cuda" ]; then
  echo "Running CUDA test"
  cmsRun "${SCRIPT}" runNumber=$3 buBaseDir=${DATA} fuBaseDir=${DATA} buNumStreams=6 numberOfEvents=$2 backend=cuda_async
elif [ "${TARGET}" == "rocm" ]; then
  echo "Running ROCm test"
  cmsRun "${SCRIPT}" runNumber=$3 buBaseDir=${DATA} fuBaseDir=${DATA} buNumStreams=1 numberOfEvents=$2 backend=rocm_async
fi