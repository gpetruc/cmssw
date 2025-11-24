#include <cstdio>
#include <random>
#include <iostream>
#include <chrono>

#include <alpaka/alpaka.hpp>

#include "FWCore/Utilities/interface/stringize.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"
#include "L1TriggerScouting/Phase2/interface/alpaka/L1TScPhase2PuppiRawToDigiKernels.h"

using namespace ALPAKA_ACCELERATOR_NAMESPACE;

void testUnpacker(std::vector<uint32_t> const& offsets, std::vector<uint64_t> const& rawdata, unsigned int ntrials = 100) {

  // run the test on each device
  for (auto const& device : cms::alpakatools::devices<Platform>()) {
    std::cout << "Test unpacking on " << alpaka::getName(device) << " over " << offsets.size()-1 << " BXs and " << rawdata.size() << " elements\n";
    auto queue = Queue(device);

    auto puppi = ALPAKA_ACCELERATOR_NAMESPACE::l1sc::PuppiDeviceCollection(rawdata.size(), queue);

    ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels::L1TScPhase2PuppiRawToDigiKernels kernels(queue);
    // move host residing data to device memory space
    alpaka::wait(queue);


    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    for (unsigned int i = 0; i < ntrials; ++i) {
      ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels::decode(queue, const_cast<uint64_t*>(rawdata.data()), puppi);
    }

    // wait for all the operations to complete
    alpaka::wait(queue);
    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Test on " << alpaka::getName(device) << " completed in " << 1000*elapsed.count()/ntrials << " milliseconds.\n";

  }
}

int main(int argc, char** argv) {
  if (argc != 2) {
    std::cout << "Usage: " << argv[0] << " <input file>\n";
    return 2;
  }
  std::cout << "Running alpaka tests for the " EDM_STRINGIZE(ALPAKA_ACCELERATOR_NAMESPACE) " backend on " << argv[1] << "\n";
  std::fstream open(argv[1], std::ios::in | std::ios::binary);
  if (!open.is_open()) {
    std::cout << "Error opening file " << argv[1] << "\n";
    return 3;
  }

  // get the size of the file
  open.seekg(0, std::ios::end);
  std::streamsize file_size = open.tellg();
  open.seekg(0, std::ios::beg);
  std::cout << "File size: " << file_size << " bytes\n";

  std::vector<uint64_t> rawdata(file_size / sizeof(uint64_t), 0u);
  std::vector<uint32_t> offsets(1, 0u);
  for (;;) {
    uint64_t word = 0;
    open.read(reinterpret_cast<char*>(&word), sizeof(uint64_t));
    if (open.eof())
      break;
    uint32_t ncands = word & 0xFFF;
    if (ncands > 0) {
      open.read(reinterpret_cast<char*>(rawdata.data() + offsets.back()), ncands * sizeof(uint64_t));
    }
    offsets.push_back(offsets.back()+ncands);
  }
  rawdata.resize(offsets.back());
  std::cout << "Read " << offsets.size()-1 << " events with a total of " << rawdata.size() << " candidates\n";

  auto const& devices = cms::alpakatools::devices<Platform>();
  if (devices.empty()) {
    std::cout << "No devices available for the " EDM_STRINGIZE(ALPAKA_ACCELERATOR_NAMESPACE) " backend, "
                   "the test will be skipped.\n";
    return 1;
  }

  testUnpacker(offsets, rawdata);

  return 0;
}