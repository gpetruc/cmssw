#include "L1TriggerScouting/Phase2/plugins/alpaka/L1TScPhase2VertexRawToDigiKernels.h"

#include "HeterogeneousCore/AlpakaInterface/interface/prefixScan.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels {

  L1TScPhase2VertexRawToDigiKernels::L1TScPhase2VertexRawToDigiKernels(Queue& queue) { initialize(queue); }

  // Initialize device constant memory for the kernels.
  // Called only once (thread-safe)
  void L1TScPhase2VertexRawToDigiKernels::initialize(Queue& queue) {}

  // Convert raw data to VertexDeviceCollection
  // Takes 64bit words and decodes them into real values for further analysis
  class RawToDigiKernel {
  public:
    ALPAKA_FN_ACC void operator()(Acc1D const& acc, data_t* p_data, VertexDeviceCollection::View Vertex) const {
      for (int32_t idx : cms::alpakatools::uniform_elements(acc, Vertex.metadata().size())) {
        uint64_t data = p_data[idx];

        // hardware values
        // auto valid = decodeBits<uint8_t, 0, l1t::VertexWord::VertexBitWidths::kValidSize>(data);
        // auto z0 = decodeBitsSigned<int16_t, l1t::VertexWord::VertexBitLocations::kZ0LSB, l1t::VertexWord::VertexBitWidths::kZ0Size>(data);
        // auto multIn = decodeBits<uint16_t, l1t::VertexWord::VertexBitLocations::kNTrackInPVLSB, l1t::VertexWord::VertexBitWidths::kNTrackInPVSize>(data);
        // auto sumPt = decodeBits<uint16_t, l1t::VertexWord::VertexBitLocations::kSumPtLSB, l1t::VertexWord::VertexBitWidths::kSumPtSize>(data);
        // auto quality = decodeBits<uint8_t, l1t::VertexWord::VertexBitLocations::kQualityLSB, l1t::VertexWord::VertexBitWidths::kQualitySize>(data);
        // auto multOut = decodeBits<uint16_t, l1t::VertexWord::VertexBitLocations::kNTrackOutPVLSB, l1t::VertexWord::VertexBitWidths::kNTrackOutPVSize>(data);
        // auto unassigned = decodeBits<uint16_t, l1t::VertexWord::VertexBitLocations::kUnassignedLSB, l1t::VertexWord::VertexBitWidths::kUnassignedSize>(data);

        auto valid = decodeBits<uint8_t, 0, 1>(data);
        auto z0 = decodeBitsSigned<int16_t, 1, 15>(data);
        auto multIn = decodeBits<uint16_t, 16, 8>(data);
        auto sumPt = decodeBits<uint16_t, 24, 12>(data);
        auto quality = decodeBits<uint8_t, 36, 3>(data);
        auto multOut = decodeBits<uint16_t, 39, 10>(data);
        auto unassigned = decodeBits<uint16_t, 49, 15>(data);

        // convert to real values
        Vertex.valid()[idx] = valid;
        Vertex.z0()[idx] = static_cast<float>(z0) * (1.0f / 512.0f);
        Vertex.multIn()[idx] = multIn;
        Vertex.sumPt()[idx] = static_cast<float>(sumPt) * 0.25f;
        Vertex.quality()[idx] = quality;
        Vertex.multOut()[idx] = multOut;
        Vertex.unassigned()[idx] = unassigned;
      }
    }
  };

  void decode(Queue& queue, data_t* p_data, VertexDeviceCollection& Vertex) {
    // move host residing data to device memory space
    auto extent = Vec1D{Vertex.const_view().metadata().size()};
    auto p_data_device = alpaka::allocAsyncBuf<data_t, Idx>(queue, extent);
    alpaka::memcpy(queue, p_data_device, createView(cms::alpakatools::host(), p_data, extent));

    // grid dims can be tuned for performance
    uint32_t threads_per_block = 1024;
    uint32_t blocks_per_grid = cms::alpakatools::divide_up_by(Vertex.const_view().metadata().size(), threads_per_block);
    auto grid = cms::alpakatools::make_workdiv<Acc1D>(blocks_per_grid, threads_per_block);

    // decode particles features
    alpaka::exec<Acc1D>(queue, grid, RawToDigiKernel{}, p_data_device.data(), Vertex.view());
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels
