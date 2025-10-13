#include "L1TriggerScouting/Phase2/interface/alpaka/L1TScPhase2SCJetsKernels.h"

#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"
#include "HeterogeneousCore/AlpakaInterface/interface/host.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"
#include "HeterogeneousCore/AlpakaInterface/interface/prefixScan.h"
#include "HeterogeneousCore/AlpakaMath/interface/deltaPhi.h"

//#define L1TSC_VERBOSE_DEBUG

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels {

  using namespace cms::alpakatools;

  class JetKernel {
  public:
    template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
    ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                  PuppiDeviceCollection::ConstView puppi,
                                  OffsetsSoA::ConstView bx_lookup,
                                  float R2,
                                  ClustersDeviceCollection::View clusters,
                                  ClusterObjDeviceCollection::View jets) const {
      uint32_t grid_dim = alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0];
      for (uint32_t block_idx : independent_groups(acc, grid_dim)) {
        // get event range
        uint32_t begin = bx_lookup.offsets()[block_idx];
        uint32_t end = bx_lookup.offsets()[block_idx + 1];
        // skip if malformed or empty
        if (end <= begin)
          continue;

        uint32_t block_dim = end - begin;
        // pre-cluster
        for (uint32_t tid : independent_group_elements(acc, block_dim)) {
          // try this as a seed
          uint32_t iseed = tid + begin;  // global index
          float seed_pt = puppi.pt()[iseed], seed_eta = puppi.eta()[iseed], seed_phi = puppi.phi()[iseed];
          float sum_pt = seed_pt, sum_eta = 0, sum_phi = 0;
          clusters.is_seed()[iseed] = 1;
          for (uint32_t j = 0; j < block_dim; ++j) {
            if (j == tid)
              continue;
            float deta = puppi.eta()[j + begin] - seed_eta;
            float dphi = cms::alpakatools::deltaPhi(acc, puppi.phi()[j + begin], seed_phi);
            if (deta * deta + dphi * dphi < R2) {
              if (puppi.pt()[j + begin] > seed_pt || (puppi.pt()[j + begin] == seed_pt && j < tid)) {
                clusters.is_seed()[iseed] = 0;
                break;
              } else {
                sum_pt += puppi.pt()[j + begin];
                sum_eta += puppi.pt()[j + begin] * deta;
                sum_phi += puppi.pt()[j + begin] * dphi;
              }
            }
          }
          sum_eta = seed_eta + sum_eta / sum_pt;
          sum_phi = cms::alpakatools::reducePhiRange(acc, seed_phi + sum_phi / sum_pt);
          jets.pt()[iseed] = clusters.is_seed()[iseed] ? sum_pt : 0;
          jets.eta()[iseed] = clusters.is_seed()[iseed] ? sum_eta : 0;
          jets.phi()[iseed] = clusters.is_seed()[iseed] ? sum_phi : 0;
          jets.cluster()[iseed] = clusters.is_seed()[iseed] ? iseed : 0;
#ifdef L1TSC_VERBOSE_DEBUG
          if (block_idx <= 2)
            if (clusters.is_seed()[iseed])
              printf("Jet pt %7.2f eta %+6.3f phi %+6.3f, seed %d\n\n",
                    jets.pt()[iseed],
                    jets.eta()[iseed],
                    jets.phi()[iseed],
                    iseed);
#endif
        }
        alpaka::syncBlockThreads(acc);

        // reassociate
        for (uint32_t tid : independent_group_elements(acc, block_dim)) {
          auto ipart = tid + begin;  // global index
          float nearest = R2;
          clusters.cluster()[ipart] = -1;
          for (uint32_t j = 0; j < block_dim; ++j) {
            auto jseed = j + begin;  // global index
            if (!clusters.is_seed()[jseed])
              continue;
            float deta = puppi.eta()[ipart] - jets.eta()[jseed];
            float dphi = cms::alpakatools::deltaPhi(acc, puppi.phi()[ipart], jets.phi()[jseed]);
            float dr2 = deta * deta + dphi * dphi;
            if (dr2 < nearest) {
              clusters.cluster()[ipart] = jseed;
              nearest = dr2;
            }
          }
        }
      }  // block
    }  // operator()
  };  // class

  class JetIterKernel {
  public:
    template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
    ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                  ClusterObjDeviceCollection::View puppi,
                                  OffsetsSoA::ConstView bx_lookup,
                                  float R2,
                                  unsigned int nIters,
                                  ClustersDeviceCollection::View clusters,
                                  uint32_t* tag,
                                  ClusterObjDeviceCollection::View work2,
                                  ClusterObjDeviceCollection::View jets) const {
      // for prefix scan (only on GPU)
      uint32_t* ws = nullptr;
      [[maybe_unused]] constexpr bool single_thread = requires_single_thread_per_block<TAcc>::value;
      if constexpr (!requires_single_thread_per_block_v<TAcc>) {
        ws = alpaka::getDynSharedMem<uint32_t>(acc);
      }
      uint32_t grid_dim = alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0];
      for (uint32_t block_idx : independent_groups(acc, grid_dim)) {
        // get event range
        uint32_t begin = bx_lookup.offsets()[block_idx];
        uint32_t end = bx_lookup.offsets()[block_idx + 1];

        // skip if empty
        if (end <= begin)
          continue;

        auto& size = alpaka::declareSharedVar<uint32_t, __COUNTER__>(acc);
        size = end - begin;

#ifdef L1TSC_VERBOSE_DEBUG
        if (once_per_block(acc) && (block_idx <= 2))
            printf("In BX %u begin with %u PF candidates: \n", block_idx + 1, end - begin);
#endif
        // running sums (accumulating on multiple threads)
        auto& seed_pt = alpaka::declareSharedVar<float, __COUNTER__>(acc);
        auto& seed_eta = alpaka::declareSharedVar<float, __COUNTER__>(acc);
        auto& seed_phi = alpaka::declareSharedVar<float, __COUNTER__>(acc);
        auto& seed_i = alpaka::declareSharedVar<unsigned int, __COUNTER__>(acc);
        auto& sum_pt = alpaka::declareSharedVar<float, __COUNTER__>(acc);
        auto& sum_eta = alpaka::declareSharedVar<float, __COUNTER__>(acc);
        auto& sum_phi = alpaka::declareSharedVar<float, __COUNTER__>(acc);

        for (unsigned int iter = 0; iter < nIters; ++iter) {
          bool even = (iter % 2 == 0);
          auto pt = even ? puppi.pt() : work2.pt();
          auto eta = even ? puppi.eta() : work2.eta();
          auto phi = even ? puppi.phi() : work2.phi();
          auto cluster = even ? puppi.cluster() : work2.cluster();
          auto pt2 = !even ? puppi.pt() : work2.pt();
          auto eta2 = !even ? puppi.eta() : work2.eta();
          auto phi2 = !even ? puppi.phi() : work2.phi();
          auto cluster2 = !even ? puppi.cluster() : work2.cluster();

          // seeding (identical on all threads)
          if (once_per_block(acc)) {
            float spt = 0, seta = 0, sphi = 0;
            unsigned int iseed = end;
            for (unsigned int j = begin, myend = begin + size; j < myend; ++j) {
  #ifdef L1TSC_VERBOSE_DEBUG
              if ((block_idx <= 2) && (iter == 0) && single_thread) {
                printf("  %4u: pt %7.2f eta %+6.3f phi %+6.3f  index %7d\n",
                      j,
                      pt[j],
                      eta[j],
                      puppi.phi()[j],
                      puppi.cluster()[j] - begin);
              }
#endif
              if (pt[j] > spt) {
                spt = pt[j];
                seta = eta[j];
                sphi = phi[j];
                iseed = j;
              }
            }
            seed_pt = spt;
            seed_eta = seta;
            seed_phi = sphi;
            seed_i = iseed;
            sum_pt = 0;
            sum_eta = 0;
            sum_phi = 0;
#ifdef L1TSC_VERBOSE_DEBUG
            if (block_idx <= 2)
              printf(
                  "In BX %u selected %u (pt %7.2f eta %+6.3f phi %+6.3f) at %u as seed for iteration %u (%u/%u particles left)\n",
                  block_idx + 1,
                  iseed - begin,
                  bestpt,
                  seed_eta,
                  seed_phi,
                  jseed,
                  iter,
                  size,
                  end - begin);            
#endif
          }

          alpaka::syncBlockThreads(acc);

          if (seed_pt == 0)
            break;

          for (uint32_t tid : independent_group_elements(acc, size)) {
            auto ipart = tid + begin;  // global index
            float deta = eta[ipart] - seed_eta;
            float dphi = cms::alpakatools::deltaPhi(acc, phi[ipart], seed_phi);
            float dr2 = deta * deta + dphi * dphi;
            if (dr2 < R2) {
              clusters.is_seed()[cluster[ipart]] = (ipart == seed_i ? 1 : 0);
              clusters.cluster()[cluster[ipart]] = iter;
              tag[ipart] = 0;
              alpaka::atomicAdd(acc, &sum_pt, pt[ipart], alpaka::hierarchy::Blocks{});
              alpaka::atomicAdd(acc, &sum_eta, deta * pt[ipart], alpaka::hierarchy::Blocks{});
              alpaka::atomicAdd(acc, &sum_phi, dphi * pt[ipart], alpaka::hierarchy::Blocks{});
#ifdef L1TSC_VERBOSE_DEBUG
              if (block_idx <= 2 && single_thread)
                printf("  %4u: pt %7.2f eta %+6.3f phi %+6.3f  cluster %7d <<= selected (dr %7.4f)\n",
                       ipart,
                       pt[ipart],
                       eta[ipart],
                       phi[ipart],
                       cluster[ipart] - begin,
                       alpaka::math::sqrt(acc, dr2));
#endif
            } else {
              tag[ipart] = 1;
            }
          }  // elements

          alpaka::syncBlockThreads(acc);

          if (once_per_block(acc)) {
            jets.pt()[begin + iter] = sum_pt;
            jets.eta()[begin + iter] = seed_eta + sum_eta / sum_pt;
            jets.phi()[begin + iter] = cms::alpakatools::reducePhiRange(acc, seed_phi + sum_phi / sum_pt);
#ifdef L1TSC_VERBOSE_DEBUG
            if (block_idx <= 2)
              printf("In BX %u Jet pt %7.2f eta %+6.3f phi %+6.3f, seed %d\n\n",
                    block_idx + 1,
                    jets.pt()[begin + iter],
                    jets.eta()[begin + iter],
                    jets.phi()[begin + iter],
                    iseed);
#endif
          }

          blockPrefixScan(acc, tag + begin, size, ws);

#ifdef L1TSC_VERBOSE_DEBUG
          if (once_per_block(acc) && (block_idx <= 2) && single_thread)
            printf("Reordering candidates\n");
#endif
          for (uint32_t tid : independent_group_elements(acc, size)) {
            auto ipart = tid + begin;  // global index
            if (tag[ipart] > (tid == 0 ? 0 : tag[ipart - 1])) {
              int dest = begin + tag[ipart] - 1;
              pt2[dest] = pt[ipart];
              eta2[dest] = eta[ipart];
              phi2[dest] = phi[ipart];
              cluster2[dest] = cluster[ipart];
#ifdef L1TSC_VERBOSE_DEBUG
              if (block_idx <= 2 && single_thread)
                printf("  %4u: pt %7.2f eta %+6.3f phi %+6.3f  cluster %7d tag %u --> %8u\n",
                       ipart,
                       pt[ipart],
                       eta[ipart],
                       phi[ipart],
                       cluster[ipart] - begin,
                       tag[ipart],
                       dest);
#endif
            }
          }
          if (once_per_block(acc)) {
            size = tag[begin + size - 1];
#ifdef L1TSC_VERBOSE_DEBUG
            if (block_idx <= 2)
              printf("In BX %u size updated to %u\n\n", block_idx + 1, size);
#endif
          }
          alpaka::syncBlockThreads(acc);

        }  // iter

      }  // block
    }  // operator()
  };  // class

  L1TScPhase2SCJetsKernels::L1TScPhase2SCJetsKernels() {}

  void L1TScPhase2SCJetsKernels::run(Queue& queue,
                                     const PuppiDeviceCollection& src,
                                     const BxLookupDeviceCollection& bx_lookup,
                                     float R2,
                                     ClustersDeviceCollection& clusters,
                                     ClusterObjDeviceCollection& jets) const {
    uint32_t threads_per_block = 256;
    uint32_t blocks_per_grid = bx_lookup.const_view<OffsetsSoA>().metadata().size() - 1;
    auto grid = make_workdiv<Acc1D>(blocks_per_grid, threads_per_block);

    clusters.zeroInitialise(queue);
    jets.zeroInitialise(queue);

    alpaka::exec<Acc1D>(queue,
                        grid,
                        JetKernel{},
                        src.const_view(),
                        bx_lookup.const_view<OffsetsSoA>(),
                        R2,
                        clusters.view(),
                        jets.view());
  }

  void L1TScPhase2SCJetsKernels::run(Queue& queue,
                                     const PuppiDeviceCollection& src,
                                     const BxLookupDeviceCollection& bx_lookup,
                                     float R2,
                                     unsigned int nJets,
                                     ClustersDeviceCollection& clusters,
                                     ClusterObjDeviceCollection& jets) const {
    // one grid per particles in blocks per BX
    uint32_t threads_per_block = 256;
    uint32_t blocks_per_grid = bx_lookup.const_view<OffsetsSoA>().metadata().size() - 1;
    auto grid = make_workdiv<Acc1D>(blocks_per_grid, threads_per_block);

    // space for output and for reordering inputs
    auto work = ClusterObjDeviceCollection(src.const_view().metadata().size(), queue);
    auto work2 = ClusterObjDeviceCollection(src.const_view().metadata().size(), queue);

    // a buffer space for counting items
    auto partExtent = Vec1D(src.const_view().metadata().size());
    auto h_tag_device = alpaka::allocAsyncBuf<uint32_t, Idx>(queue, partExtent);
    alpaka::memset(queue, h_tag_device, 0x00);

    // one flat grid per particle, with arbitrary block size
    uint32_t threads_per_flatblock = 1024;
    uint32_t blocks_per_flatgrid =
        cms::alpakatools::divide_up_by(src.const_view().metadata().size(), threads_per_flatblock);
    auto flatgrid = cms::alpakatools::make_workdiv<Acc1D>(blocks_per_flatgrid, threads_per_flatblock);

    jets.zeroInitialise(queue);

    alpaka::exec<Acc1D>(
        queue,
        flatgrid,
        [] ALPAKA_FN_ACC(Acc1D const& acc,
                         PuppiDeviceCollection::ConstView puppi,
                         ClusterObjDeviceCollection::View work,
                         ClustersDeviceCollection::View clusters) {
          for (int32_t idx : cms::alpakatools::uniform_elements(acc, clusters.metadata().size())) {
            work.pt()[idx] = puppi.pt()[idx];
            work.eta()[idx] = puppi.eta()[idx];
            work.phi()[idx] = puppi.phi()[idx];
            work.cluster()[idx] = idx;
            clusters.cluster()[idx] = -1;
          }
        },
        src.const_view(),
        work.view(),
        clusters.view());

    alpaka::exec<Acc1D>(queue,
                        grid,
                        JetIterKernel{},
                        work.view(),
                        bx_lookup.const_view<OffsetsSoA>(),
                        R2,
                        nJets,
                        clusters.view(),
                        h_tag_device.data(),
                        work2.view(),
                        jets.view());
  }
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels