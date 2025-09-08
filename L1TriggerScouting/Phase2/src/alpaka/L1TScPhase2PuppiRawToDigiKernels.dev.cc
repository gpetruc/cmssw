#include "L1TriggerScouting/Phase2/interface/alpaka/L1TScPhase2PuppiRawToDigiKernels.h"

#include "HeterogeneousCore/AlpakaInterface/interface/prefixScan.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels {

  std::once_flag L1TScPhase2PuppiRawToDigiKernels::init_flag_;

  L1TScPhase2PuppiRawToDigiKernels::L1TScPhase2PuppiRawToDigiKernels(Queue& queue) { initialize(queue); }

  /**
   * Initialize device constant memory for the kernels
   * @note called only once (thread-safe)
   */
  void L1TScPhase2PuppiRawToDigiKernels::initialize(Queue& queue) {
    std::call_once(init_flag_, [&]() {
      // pdgid mapping
      constexpr int16_t host_pdgid[8] = {130, 22, -211, 211, 11, -11, 13, -13};
      auto extent = Vec1D{8};
      auto view = createView(cms::alpakatools::host(), host_pdgid, extent);
      alpaka::memcpy(queue, kPdgid<Acc1D>, view);

      // hw to float conversion: 3.14 / 720.0
      constexpr float host_pi_720 = alpaka::math::constants::pi / 720.0f;
      auto extent_var = Vec1D{1};
      auto view_var = createView(cms::alpakatools::host(), &host_pi_720, extent_var);
      alpaka::memcpy(queue, kPi720<Acc1D>, view_var);
    });
  }

  /**
   * Convert raw data to PuppiDeviceCollection
   * Takes 64bit words and decodes them into real values for further analysis
   */
  class RawToDigiKernel {
  public:
    ALPAKA_FN_ACC void operator()(Acc1D const& acc, data_t* p_data, PuppiDeviceCollection::View puppi) const {
      for (int32_t idx : cms::alpakatools::uniform_elements(acc, puppi.metadata().size())) {
        uint64_t data = p_data[idx];

        auto hwPt = decodeBits<uint16_t>(data, 0, 14);
        auto hwEta = decodeBitsSigned<int16_t>(data, 14, 12);
        auto hwPhi = decodeBitsSigned<int16_t>(data, 26, 11);
        auto pid = decodeBits<uint8_t>(data, 37, 3);

        puppi.pt()[idx] = hwPt * 0.25f;
        puppi.eta()[idx] = hwEta * kPi720<Acc1D>.get();
        puppi.phi()[idx] = hwPhi * kPi720<Acc1D>.get();
        puppi.pdgid()[idx] = kPdgid<Acc1D>.get()[pid];

        if (pid > 1) {
          auto hwZ0 = decodeBitsSigned<int16_t>(data, 40, 10);
          auto hwDxy = decodeBitsSigned<int16_t>(data, 50, 8);
          auto hwQual = decodeBits<uint8_t>(data, 58, 3);

          puppi.z0()[idx] = hwZ0 * 0.05f;
          puppi.dxy()[idx] = hwDxy * 0.05f;
          puppi.puppiw()[idx] = 1.0f;
          puppi.quality()[idx] = hwQual;
        } else {
          auto hwPuppiw = decodeBits<uint16_t>(data, 40, 10);
          auto hwQual = decodeBits<uint8_t>(data, 50, 6);

          puppi.z0()[idx] = 0.0f;
          puppi.dxy()[idx] = 0.0f;
          puppi.puppiw()[idx] = hwPuppiw * (1.0f / 256.0f);
          puppi.quality()[idx] = hwQual;
        }

        puppi.selection()[idx] = 0;
      }
    }
  };

  void rawToDigi(Queue& queue, data_t* p_data, PuppiDeviceCollection& puppi) {
    // move host residing data to device memory space
    auto extent = Vec1D(puppi.const_view().metadata().size());
    auto p_data_device = alpaka::allocAsyncBuf<data_t, Idx>(queue, extent);
    alpaka::memcpy(queue, p_data_device, createView(cms::alpakatools::host(), p_data, extent));

    uint32_t threads_per_block = 1024;
    uint32_t blocks_per_grid = cms::alpakatools::divide_up_by(puppi.const_view().metadata().size(), threads_per_block);
    auto grid = cms::alpakatools::make_workdiv<Acc1D>(blocks_per_grid, threads_per_block);

    // decode particles features
    alpaka::exec<Acc1D>(queue, grid, RawToDigiKernel{}, p_data_device.data(), puppi.view());
  }

  void associateOrbitEventIndex(Queue& queue,
                                data_t* h_data,
                                OrbitEventIndexMapDeviceCollection& orbit_association_map) {
    // move host residing data to device memory space
    auto extent = Vec1D(orbit_association_map.const_view().metadata().size());
    auto h_data_device = alpaka::allocAsyncBuf<data_t, Idx>(queue, extent);
    alpaka::memcpy(queue, h_data_device, createView(cms::alpakatools::host(), h_data, extent));

    uint32_t threads_per_block = 1024;
    uint32_t blocks_per_grid =
        cms::alpakatools::divide_up_by(orbit_association_map.const_view().metadata().size(), threads_per_block);
    auto grid = cms::alpakatools::make_workdiv<Acc1D>(blocks_per_grid, threads_per_block);

    // accumulate buffer with events sizes
    alpaka::exec<Acc1D>(
        queue,
        grid,
        [] ALPAKA_FN_ACC(
            Acc1D const& acc, data_t* data, OrbitEventIndexMapDeviceCollection::View orbit_association_map) {
          if (cms::alpakatools::once_per_grid(acc))
            orbit_association_map.offsets()[0] = 0;

          for (int32_t idx : cms::alpakatools::uniform_elements(acc, orbit_association_map.metadata().size() - 1)) {
            auto range = decodeBits<uint32_t>(data[idx], 0, 12);
            orbit_association_map.offsets()[idx + 1] = range;
          }
        },
        h_data_device.data(),
        orbit_association_map.view());

    // prefix sum to build association map used for span extraction and batching
    auto pc = alpaka::allocAsyncBuf<int32_t, Idx>(queue, Vec1D{1});
    alpaka::memset(queue, pc, 0x0);
    alpaka::exec<Acc1D>(queue,
                        grid,
                        cms::alpakatools::multiBlockPrefixScan<uint32_t>{},
                        orbit_association_map.view().offsets() + 1,
                        orbit_association_map.view().offsets() + 1,
                        orbit_association_map.view().metadata().size(),
                        blocks_per_grid,
                        pc.data(),
                        alpaka::getPreferredWarpSize(alpaka::getDev(queue)));
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels