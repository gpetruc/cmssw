#include "L1TriggerScouting/TauTagging/plugins/alpaka/CLUEsteringAlgo.h"

#include <any>

#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"
#include "HeterogeneousCore/AlpakaInterface/interface/radixSort.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels {

  template<typename TAcc, typename T>
  inline ALPAKA_FN_ACC void swap(TAcc const& acc, T &a, T &b) {
    T temp = a;
    a = b;
    b = temp;
  }

  using namespace cms::alpakatools;

  class SortClustersKernel {
  public:
    ALPAKA_FN_ACC void operator()(
        Acc1D const& acc,
        const float* weights,
        IndexSoA::View indexes, 
        OffsetsSoA::View offsets) const {
      const uint8_t kSharedMemSize = 128;
      auto& indices_shared = alpaka::declareSharedVar<int[kSharedMemSize], __COUNTER__>(acc);
      auto& weights_shared = alpaka::declareSharedVar<float[kSharedMemSize], __COUNTER__>(acc); 

      // loop over clusters in parallel
      for (uint32_t block_idx: independent_groups(acc, offsets.metadata().size() - 1)) {
        // bind range to hw block
        uint32_t begin = offsets.offsets()[block_idx];
        uint32_t end = offsets.offsets()[block_idx + 1];
        // define block dimensions
        uint32_t block_dim = end - begin;
        if (block_dim == 0)
          continue;

        // load global to shared memory with EOF sentinels
        for (uint32_t tid : independent_group_elements(acc, kSharedMemSize)) {
          if (tid < block_dim) {
            uint32_t thread_idx = begin + tid;
            auto p_index = indexes.indexes()[thread_idx];
            indices_shared[tid] = p_index;
            weights_shared[tid] = weights[p_index];
          } else {
            // sentinel so unused slots never win
            indices_shared[tid] = -1;
            weights_shared[tid] = -1.0f;
          }
        }
        alpaka::syncBlockThreads(acc);

        // odd-even sort algorithm
        // this should be replaced by radixSortMulti?
        // but in the near future CLUE will provide sorted clusters directly
        for (uint32_t i = 0; i < block_dim; i++) {
          for (uint32_t tid : independent_group_elements(acc, block_dim - 1)) {
            if (tid + 1 < block_dim) {
              if ((i % 2 == 0 && tid % 2 == 0) || (i % 2 == 1 && tid % 2 == 1)) {
                if (weights_shared[tid] < weights_shared[tid + 1]) {
                  swap(acc, weights_shared[tid], weights_shared[tid + 1]);
                  swap(acc, indices_shared[tid], indices_shared[tid + 1]);
                }
              }
            }
            // sync tree
            alpaka::syncBlockThreads(acc);
          }
        }

        // write back shared to global memory
        for (uint32_t tid : independent_group_elements(acc, block_dim)) {
          uint32_t thread_idx = tid + begin;
          indexes.indexes()[thread_idx] = indices_shared[tid];
        }
      }
    }
  };

  class CopyAssociatorKernel {
  public:
    ALPAKA_FN_ACC void operator()(
        Acc1D const& acc, 
        clue::AssociationMapView associator, 
        IndexSoA::View indexes, 
        OffsetsSoA::View offsets, 
        const float* weights) const {
      if (once_per_grid(acc)) {
        offsets.offsets()[0] = 0;
        for (int c_id = 0; c_id < offsets.metadata().size() - 1; c_id++) {
          auto span = associator[c_id];
          offsets.offsets()[c_id+1] = span.size() + offsets.offsets()[c_id];
          auto begin = offsets.offsets()[c_id];
          for (int i = 0; i < span.size(); i++) {
            indexes.indexes()[i+begin] = span[i];
          }
        }
      }
    }
  };

  CLUEsteringAlgo::CLUEsteringAlgo(float dc, float rhoc, float dm, bool wrap_coords)
      : dc_(dc), rhoc_(rhoc), dm_(dm), wrap_coords_(wrap_coords) {}

  AssociationMapDevice CLUEsteringAlgo::run(Queue& queue, const PFCandidateDeviceCollection& pf, ClustersDeviceCollection& clusters) const {
    const uint32_t n_points = pf.const_view().metadata().size();

    // buffers
    // CLUEstering call internally reinterpret_cast<T*> to non-const ptr
    auto* eta_coord_ptr = const_cast<float*>(pf.const_view().eta().data());
    auto* phi_coord_ptr = const_cast<float*>(pf.const_view().phi().data());
    auto* weights_ptr = const_cast<float*>(pf.const_view().pt().data());
    auto* clusters_ptr = clusters.view().cluster().data();

    // wrap device buffers
    auto points_device =
        clue::PointsDevice<kDims, Device>(queue, n_points, eta_coord_ptr, phi_coord_ptr, weights_ptr, clusters_ptr);
    // run (wrap coords if enabled)
    auto clue_algo = clue::Clusterer<kDims>(queue, dc_, rhoc_, dm_);
    if (wrap_coords_)
      clue_algo.setWrappedCoordinates({{0, 1}});
    clue_algo.make_clusters(queue, points_device);
    auto associator = clue_algo.getClusters(queue, points_device);

    auto association_soa = AssociationMapDevice({{static_cast<int>(n_points), static_cast<int>(associator.size()+1)}}, queue);
    association_soa.zeroInitialise(queue);

    // copy clue::AssociationMapView (should be parallelized inside kernel)
    alpaka::exec<Acc1D>(queue, 
      make_workdiv<Acc1D>(1, 1), 
      CopyAssociatorKernel{}, 
      associator.view(), 
      association_soa.view<IndexSoA>(), 
      association_soa.view<OffsetsSoA>(),
      pf.const_view().pt().data());

    // sort clusters by pt
    // TODO: this will be done inside CLUE in the future, efficient copy required
    // right now all underlying buffers are private
    alpaka::exec<Acc1D>(queue, 
      make_workdiv<Acc1D>(associator.size(), 128), 
      SortClustersKernel{}, 
      pf.const_view().pt().data(),
      association_soa.view<IndexSoA>(), 
      association_soa.view<OffsetsSoA>());

    return association_soa;
  }

  AssociationMapDevice CLUEsteringAlgo::run(Queue& queue,
                            const PFCandidateDeviceCollection& pf,
                            const BxLookupDeviceCollection& bx_lookup,
                            ClustersDeviceCollection& clusters) const {
    // TODO: CLUE is not yet adapted to run on multiple BXs and batch efficiently, for loop required.
    const auto nbx = static_cast<int32_t>(bx_lookup.const_view<BxIndexSoA>().metadata().size());
    auto bx_lookup_host = BxLookupHostCollection({{nbx, nbx + 1}}, queue);
    alpaka::memcpy(queue, bx_lookup_host.buffer(), bx_lookup.buffer());
    alpaka::wait(queue);

    // allocate a big map to store all BXs associations, will be compacted at the end
    // by design the number of clusters is <= number of PF candidates
    const auto n_pf = pf.const_view().metadata().size();
    auto association_nonzs = AssociationMapDevice({{n_pf, n_pf+1}}, queue);
    association_nonzs.zeroInitialise(queue);

    unsigned int clustersTotal = 0, clusteredTotal = 0;
    for (int32_t idx = 0; idx < bx_lookup_host.const_view<BxIndexSoA>().metadata().size(); idx++) {
      const auto begin = bx_lookup_host.const_view<OffsetsSoA>().offsets()[idx];
      const auto end = bx_lookup_host.const_view<OffsetsSoA>().offsets()[idx + 1];
      const uint32_t n_points = end - begin;

      if (n_points == 0) {
        continue;
      }

      // buffers
      // CLUEstering call internally reinterpret_cast<T*> to non-const ptr
      auto* eta_coord_ptr = const_cast<float*>(pf.const_view().eta().data() + begin);
      auto* phi_coord_ptr = const_cast<float*>(pf.const_view().phi().data() + begin);
      auto* weights_ptr = const_cast<float*>(pf.const_view().pt().data() + begin);
      auto* clusters_ptr = clusters.view().cluster().data() + begin;

      // wrap device buffers
      auto points_device =
          clue::PointsDevice<kDims, Device>(queue, n_points, eta_coord_ptr, phi_coord_ptr, weights_ptr, clusters_ptr);
      auto clue_algo = clue::Clusterer<kDims>(queue, dc_, rhoc_, dm_);
      if (wrap_coords_)
        clue_algo.setWrappedCoordinates({{0, 1}});
      clue_algo.make_clusters(queue, points_device);
#if 0
      // temporary workaround: check if any cluster found
      auto p_found = alpaka::allocAsyncBuf<uint32_t, Idx>(queue, Vec1D{1});
      alpaka::exec<Acc1D>(queue, 
        make_workdiv<Acc1D>(1, 1), 
        [] ALPAKA_FN_ACC(Acc1D const& acc,
                        const int * clusters,
                        uint32_t npoints,
                        uint32_t* nfound) {
          if (once_per_grid(acc)) {
            for (int i = 0; i < npoints; i++) {
              if (clusters[i] != -1) {
                (*nfound)++;
                break;
              }
            }
          }
        }, 
        clusters.view().cluster().data() + begin,
        n_points,
        p_found);
      auto p_found_host = alpaka::allocBuf<uint32_t, Idx>(alpaka::getDevHost(), Vec1D{1});
      alpaka::memcpy(queue, p_found_host, p_found);
      alpaka::wait(queue);
      if (p_found_host.data()[0] == 0)
        continue;
#endif      
      // get clusters
      auto associator = clue_algo.getClusters(queue, points_device);

      // concatenate the association maps (could be done better!)
      auto nclusters = associator.size();
      auto nclustered = associator.extents().values;
      //std::cout << "BX " << idx << ": found " << nclusters << " clusters from " << n_points << " PF candidates of which " << nclustered << " are clustered." << std::endl;
      alpaka::exec<Acc1D>(queue, 
        make_workdiv<Acc1D>(1, 1), 
        [] ALPAKA_FN_ACC(Acc1D const& acc,
                         clue::AssociationMapView associator, 
                         IndexSoA::View indexes, 
                         OffsetsSoA::View offsets,
                         int nclusters,
                         unsigned int begin_clusters,
                         unsigned int begin_clustered) {
          if (once_per_grid(acc)) {
            offsets.offsets()[begin_clusters] = begin_clustered;
            for (int c_id = 0; c_id < nclusters; c_id++) {
              auto span = associator[c_id];
              offsets.offsets()[begin_clusters+c_id+1] = span.size() + offsets.offsets()[begin_clusters+c_id];
              auto begin = offsets.offsets()[begin_clusters+c_id];
              for (int i = 0; i < span.size(); i++) {
                indexes.indexes()[i+begin] = span[i];
              }
            }
          }
        },        
        associator.view(), 
        association_nonzs.view<IndexSoA>(), 
        association_nonzs.view<OffsetsSoA>(),
        nclusters,
        clustersTotal,
        clusteredTotal);
        clustersTotal += nclusters;
        clusteredTotal += nclustered;
    }

    // make final association map copying only used entries
    // could be done better
    auto association_soa = AssociationMapDevice({{int(clusteredTotal), int(clustersTotal+1)}}, queue);
    association_soa.zeroInitialise(queue);
    auto srcIdx = alpaka::createView(alpaka::getDev(queue), 
                             association_nonzs.view<IndexSoA>().indexes().data(),
                             Vec1D{clusteredTotal});
    auto dstIdx = alpaka::createView(alpaka::getDev(queue), 
                             association_soa.view<IndexSoA>().indexes().data(),
                             Vec1D{clusteredTotal});                             
    auto srcOffs = alpaka::createView(alpaka::getDev(queue), 
                             association_nonzs.view<OffsetsSoA>().offsets().data(),
                             Vec1D{clustersTotal});
    auto dstOffs = alpaka::createView(alpaka::getDev(queue), 
                             association_soa.view<OffsetsSoA>().offsets().data(),
                             Vec1D{clustersTotal});                             
    alpaka::memcpy(queue, dstIdx, srcIdx);
    alpaka::memcpy(queue, dstOffs, srcOffs);
    return association_soa;
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels
