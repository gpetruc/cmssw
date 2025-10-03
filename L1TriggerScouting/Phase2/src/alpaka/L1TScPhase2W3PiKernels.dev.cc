#include "L1TriggerScouting/Phase2/interface/alpaka/L1TScPhase2W3PiKernels.h"

#include <alpaka/alpaka.hpp>

#include "HeterogeneousCore/AlpakaInterface/interface/host.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"
#include "HeterogeneousCore/AlpakaInterface/interface/prefixScan.h"
#include "L1TriggerScouting/Phase2/interface/L1TScPhase2Common.h"
#include "L1TriggerScouting/Phase2/interface/alpaka/PhysicsKernelsUtilities.h"
#include "L1TriggerScouting/Phase2/interface/alpaka/SynchronizingTimer.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels {

  using namespace cms::alpakatools;

  // TODO: - is it possible to vectorize innermost loops?
  class W3PiKernel {
  public:
    template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
    ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                  PuppiDeviceCollection::ConstView puppi,
                                  OffsetsSoA::ConstView offsets,
                                  SelectionBitmaskDeviceCollection::View selection_bitmask,
                                  BufferDevice::View bx_mask,
                                  const W3PiAlgoParams* params,
                                  PortableCounter* w3pi_triplets_ct) const {
      // shared memory
      const uint8_t SHARED_MEM_BLOCK = 128;
      auto& bit_field = alpaka::declareSharedVar<int, __COUNTER__>(acc);
      auto& mask = alpaka::declareSharedVar<int[SHARED_MEM_BLOCK], __COUNTER__>(acc);

      if (once_per_block(acc)) {
        bit_field = 0;
      }

      uint32_t grid_dim = alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0];
      for (uint32_t block_idx : independent_groups(acc, grid_dim)) {
        // get event range
        uint32_t begin = offsets.offsets()[block_idx];
        uint32_t end = offsets.offsets()[block_idx + 1];
        // skip if malformed or empty
        if (end <= begin || end == 0xFFFFFFFF)
          continue;

        uint32_t block_dim = end - begin;
        // match types
        for (uint32_t tid : independent_group_elements(acc, block_dim)) {
          mask[tid] = 0;                  // memset mask
          auto thread_idx = tid + begin;  // global index
          auto cls = alpaka::math::abs(acc, static_cast<int>(puppi.pdgid()[thread_idx]));
          if (cls != 211 && cls != 11)
            continue;
          if (puppi.pt()[thread_idx] < params->pT_min)
            continue;
          // `independent_group_elements` should guarantee the order (if not use atomicExch)
          mask[tid] = 1;
        }
        alpaka::syncBlockThreads(acc);

        // skip combinatorial part if not enough candidates to form a triplet
        if (once_per_block(acc)) {
          int size = 0;
          for (auto idx = 0; idx < SHARED_MEM_BLOCK; ++idx)
            if (mask[idx] == 1)
              ++size;
          if (size < 3)
            continue;
        }

        // combinatorial part
        for (uint32_t tid : independent_group_elements(acc, block_dim)) {
          auto thread_idx = tid + begin;  // global index
          if (mask[tid] == 0)
            continue;
          if (puppi.pt()[thread_idx] < params->pT_max)
            continue;
          if (!coneIsolation(acc,
                             puppi,
                             thread_idx,
                             begin,
                             end,
                             params->min_deltar_threshold,
                             params->max_deltar_threshold,
                             params->max_isolation_threshold))
            continue;
          for (uint32_t i = 0; i < block_dim; i++) {
            auto global_i_idx = i + begin;
            if (mask[i] == 0)
              continue;
            if (global_i_idx == thread_idx || puppi.pt()[global_i_idx] < params->pT_int)
              continue;
            if (puppi.pt()[global_i_idx] > puppi.pt()[thread_idx] ||
                (puppi.pt()[global_i_idx] == puppi.pt()[thread_idx] && global_i_idx < thread_idx))
              continue;
            if (!angularSeparation(acc, puppi, thread_idx, global_i_idx, params->ang_sep_lower_bound))
              continue;
            for (uint32_t j = 0; j < block_dim; j++) {
              auto global_j_idx = j + begin;
              if (mask[j] == 0)
                continue;
              if (global_j_idx == thread_idx || global_j_idx == global_i_idx ||
                  puppi.pt()[global_i_idx] < params->pT_min)
                continue;
              if (puppi.pt()[global_j_idx] > puppi.pt()[thread_idx] ||
                  (puppi.pt()[global_j_idx] == puppi.pt()[thread_idx] && global_j_idx < thread_idx))
                continue;
              if (puppi.pt()[global_j_idx] > puppi.pt()[global_i_idx] ||
                  (puppi.pt()[global_j_idx] == puppi.pt()[global_i_idx] && global_j_idx < global_i_idx))
                continue;
              if (alpaka::math::abs(acc,
                                    static_cast<int>(charge(acc, puppi.pdgid()[thread_idx]) +
                                                     charge(acc, puppi.pdgid()[global_i_idx]) +
                                                     charge(acc, puppi.pdgid()[global_j_idx]))) != 1)
                continue;
              auto mass = massInvariant(acc, puppi, thread_idx, global_i_idx, global_j_idx);
              if (mass < params->invariant_mass_lower_bound || mass > params->invariant_mass_upper_bound)
                continue;
              if (angularSeparation(acc, puppi, thread_idx, global_j_idx, params->ang_sep_lower_bound) &&
                  angularSeparation(acc, puppi, global_i_idx, global_j_idx, params->ang_sep_lower_bound)) {
                if (coneIsolation(acc,
                                  puppi,
                                  global_i_idx,
                                  begin,
                                  end,
                                  params->min_deltar_threshold,
                                  params->max_deltar_threshold,
                                  params->max_isolation_threshold) &&
                    coneIsolation(acc,
                                  puppi,
                                  global_j_idx,
                                  begin,
                                  end,
                                  params->min_deltar_threshold,
                                  params->max_deltar_threshold,
                                  params->max_isolation_threshold)) {
                  // atomic then shift to avoid race conditions
                  int bit_shift = alpaka::atomicAdd(acc, &bit_field, 1);
                  uint64_t bit_mask = (1ull << bit_shift);
                  // encode combinatorial part of triplets (from LSB to MSB)
                  alpaka::atomicOr(acc, &selection_bitmask.bits()[thread_idx], bit_mask);
                  alpaka::atomicOr(acc, &selection_bitmask.bits()[global_i_idx], bit_mask);
                  alpaka::atomicOr(acc, &selection_bitmask.bits()[global_j_idx], bit_mask);
                  // mask to recover selected bxs
                  alpaka::atomicExch(acc, &bx_mask.value()[block_idx], static_cast<uint32_t>(1));
                  // counter to allocate table for decoded triplets
                  alpaka::atomicAdd(acc, &w3pi_triplets_ct->value, 1);

                  // printf("i: %u -> ", thread_idx);
                  // printBits(acc, selection_bitmask.bits()[thread_idx]);
                  // printf("j: %u -> ", global_i_idx);
                  // printBits(acc, selection_bitmask.bits()[global_i_idx]);
                  // printf("k: %u -> ", global_j_idx);
                  // printBits(acc, selection_bitmask.bits()[global_j_idx]);
                  // printf("\n");
                }
              }
            }
          }
        }
      }
    }
  };

  std::tuple<SelectedBxDeviceCollection, W3PiDeviceTable> runW3Pi(Queue& queue,
                                                                  const PuppiDeviceCollection& puppi,
                                                                  const BxLookupDeviceCollection& bx_lookup,
                                                                  const W3PiAlgoParams* params,
                                                                  const bool fast_path) {
    auto bx_mask = BufferDevice(bx_lookup.const_view().metadata().size(), queue);
    bx_mask.zeroInitialise(queue);

    // TODO: selection bitmask should be returned when `fast_path` is enabled, multi collection with masked and standard products?
    auto selection_bitmask = SelectionBitmaskDeviceCollection(puppi.const_view().metadata().size(), queue);
    selection_bitmask.zeroInitialise(queue);
    auto w3pi_triplets_ct = CounterDevice(queue);
    w3pi_triplets_ct.zeroInitialise(queue);

    uint32_t threads_per_block = make_threads_per_block(128);
    uint32_t blocks_per_grid = bx_lookup.const_view<BxIndexSoA>().metadata().size();
    auto grid = make_workdiv<Acc1D>(blocks_per_grid, threads_per_block);

    alpaka::exec<Acc1D>(queue,
                        grid,
                        W3PiKernel{},
                        puppi.const_view(),
                        bx_lookup.const_view<OffsetsSoA>(),
                        selection_bitmask.view(),
                        bx_mask.view(),
                        params,
                        w3pi_triplets_ct.data());

    // TODO: optimize this later
    auto nbx_ct = CounterDevice(queue);
    nbx_ct.zeroInitialise(queue);
    alpaka::exec<Acc1D>(
        queue,
        grid,
        [] ALPAKA_FN_ACC(Acc1D const& acc, BufferDevice::ConstView mask, PortableCounter* ct) {
          // block shared counter (for atomic coarsening)
          auto& ct_shared = alpaka::declareSharedVar<int, __COUNTER__>(acc);
          if (once_per_block(acc))
            ct_shared = 0;

          // accumulate per block memory, to reduce atomic contention
          for (int32_t thread_idx : uniform_elements(acc, mask.metadata().size())) {
            if (mask.value()[thread_idx] == 1)
              alpaka::atomicAdd(acc, &ct_shared, 1);
          }

          // block-level shared memory reduction
          alpaka::syncBlockThreads(acc);
          if (once_per_block(acc))
            alpaka::atomicAdd(acc, &ct->value, ct_shared);
        },
        bx_mask.const_view(),
        nbx_ct.data());

    if (fast_path) {
      auto selected_bxs = SelectedBxDeviceCollection(bx_mask.const_view().metadata().size(), queue);
      auto w3pi_table = W3PiDeviceTable(0, queue);  // empty table
      alpaka::memcpy(queue, selected_bxs.buffer(), bx_mask.buffer());
      return std::make_tuple(std::move(selected_bxs), std::move(w3pi_table));
    }

    // TODO: find a better way to return the requested data, all below kernels should have faster impl on CPU (assuming optimal impl):
    //   - mask might be sufficient instead of compressed indices
    //   - triplets table can be recover later without filling kernel
    auto w3pi_triplets_ct_host = CounterHost(queue);
    auto nbx_ct_host = CounterHost(queue);
    alpaka::memcpy(queue, w3pi_triplets_ct_host.buffer(), w3pi_triplets_ct.buffer());
    alpaka::memcpy(queue, nbx_ct_host.buffer(), nbx_ct.buffer());
    // explicit synchronization in production code should be avoided
    alpaka::wait(queue);

    // allocate buffers for compacted bx and triplets
    auto selected_bxs = SelectedBxDeviceCollection(nbx_ct_host.data()->value, queue);
    auto w3pi_table = W3PiDeviceTable(w3pi_triplets_ct_host.data()->value, queue);

    // return empty buffers
    if (nbx_ct_host.data()->value == 0) {
      return std::make_tuple(std::move(selected_bxs), std::move(w3pi_table));
    }

    // grid dims can be tuned for performance
    uint32_t threads_per_block_inclusive_scan = 512;
    uint32_t blocks_per_grid_inclusive_scan =
        cms::alpakatools::divide_up_by(bx_mask.view().metadata().size(), threads_per_block_inclusive_scan);
    auto grid_inclusive_scan =
        cms::alpakatools::make_workdiv<Acc1D>(blocks_per_grid_inclusive_scan, threads_per_block_inclusive_scan);

    // inclusive scan to recover compressed selected bx indices
    auto prefix_sum = BufferDevice(bx_lookup.const_view().metadata().size(), queue);
    auto pc = alpaka::allocAsyncBuf<int32_t, Idx>(queue, Vec1D{1});
    alpaka::memset(queue, pc, 0x00);
    alpaka::exec<Acc1D>(queue,
                        grid_inclusive_scan,
                        cms::alpakatools::multiBlockPrefixScan<uint32_t>{},
                        bx_mask.view().value(),
                        prefix_sum.view().value(),
                        bx_lookup.view().metadata().size(),
                        blocks_per_grid_inclusive_scan,
                        pc.data(),
                        alpaka::getPreferredWarpSize(alpaka::getDev(queue)));

    // compress mask [0 0 1 0 1 0 ... BX] into selected bxs only e.g. [124, ... 2456].
    alpaka::exec<Acc1D>(
        queue,
        grid,
        [] ALPAKA_FN_ACC(Acc1D const& acc,
                         BufferDevice::ConstView mask,
                         BufferDevice::ConstView prefixsum,
                         SelectedBxDeviceCollection::View compressed) {
          for (uint32_t thread_idx : uniform_elements(acc, mask.metadata().size())) {
            if (mask.value()[thread_idx] == 1) {
              uint32_t dest_idx = prefixsum.value()[thread_idx] - 1;  // inclusive scan!
              compressed.bx()[dest_idx] = thread_idx;
            }
          }
        },
        bx_mask.const_view(),
        prefix_sum.const_view(),
        selected_bxs.view());

    // grid dims can be tuned for performance
    uint32_t threads_per_block_fill_table = 32;
    uint32_t blocks_per_grid_fill_table =
        cms::alpakatools::divide_up_by(selected_bxs.view().metadata().size(), threads_per_block_fill_table);
    auto grid_fill_table =
        cms::alpakatools::make_workdiv<Acc1D>(blocks_per_grid_fill_table, threads_per_block_fill_table);

    auto table_pos = CounterDevice(queue);
    table_pos.zeroInitialise(queue);
    // kernel to fill w3pi table
    alpaka::exec<Acc1D>(
        queue,
        grid_fill_table,
        [] ALPAKA_FN_ACC(Acc1D const& acc,
                         PuppiDeviceCollection::ConstView puppi,
                         SelectedBxDeviceCollection::ConstView selected_bxs,
                         SelectionBitmaskDeviceCollection::ConstView bitmask,
                         OffsetsSoA::ConstView offsets,
                         W3PiDeviceTable::View w3pi_table,
                         PortableCounter* table_pos) {
          for (uint32_t thread_idx : uniform_elements(acc, selected_bxs.metadata().size())) {
            uint32_t bx_idx = selected_bxs.bx()[thread_idx];
            uint32_t begin = offsets.offsets()[bx_idx];
            uint32_t end = offsets.offsets()[bx_idx + 1];
            if (end <= begin || end == 0xFFFFFFFF)
              continue;

            // assumed 64bit words will not be exceeded at any time
            for (int bit = 0; bit < std::numeric_limits<uint64_t>::digits; ++bit) {
              int found = 0;
              uint32_t i = 0u, j = 0u, k = 0u;
              for (uint32_t p = begin; p < end; ++p) {
                if ((bitmask.bits()[p] & (1ull << bit)) != 0ull) {
                  ++found;
                  if (found == 1)
                    i = p;
                  else if (found == 2)
                    j = p;
                  else if (found == 3) {
                    k = p;
                    break;
                  }
                }
              }

              // reject combinations with less than 3 particles
              // continue to search for rest of combination
              if (found < 3)
                continue;

              // assert to be sure encoding is correct and bits are set correctly
              ALPAKA_ASSERT_ACC(found == 3);

              // fill table
              int32_t pos = alpaka::atomicAdd(acc, &table_pos->value, 1);
              w3pi_table.i()[pos] = i;
              w3pi_table.j()[pos] = j;
              w3pi_table.k()[pos] = k;
            }
          }
        },
        puppi.const_view(),
        selected_bxs.const_view(),
        selection_bitmask.const_view(),
        bx_lookup.const_view<OffsetsSoA>(),
        w3pi_table.view(),
        table_pos.data());

    return std::make_tuple(std::move(selected_bxs), std::move(w3pi_table));
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels