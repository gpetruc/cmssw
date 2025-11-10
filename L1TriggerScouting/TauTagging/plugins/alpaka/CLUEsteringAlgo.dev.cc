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

  // class GenerateAssociatorKernel {
  // public:
  //   ALPAKA_FN_ACC void operator()(
  //       Acc1D const& acc, 
  //       clue::AssociationMapView associator, 
  //       IndexSoA::View indexes, 
  //       OffsetsSoA::View offsets, 
  //       unsigned int begin_clustered) const {
  //     if (once_per_grid(acc)) {
  //       // update first offset
  //       auto span = associator[0];
  //       offsets.offsets()[0] = begin_clustered + span.size();
        
  //       // this variable keeps track of the local indexes entry that is filled
  //       unsigned int begin = 0;

  //       // update indexes corresponding to the first cluster
  //       for (unsigned int i = 0; i < span.size(); i++) {
  //           indexes.indexes()[i] = span[i];
  //       }
        
  //       for (unsigned int c_id = 1; c_id < offsets.metadata().size(); c_id++) { // ATTENTION HERE: AT LEAST TWO CLUSTERS ARE REQUIRED
  //         // update offsets
  //         span = associator[c_id];
  //         offsets.offsets()[c_id] = span.size() + offsets.offsets()[c_id - 1];

  //         // update indexes
  //         begin += associator[c_id - 1].size();
  //         for (unsigned int i = 0; i < span.size(); i++) {
  //           indexes.indexes()[i + begin] = span[i];
  //         }
  //       }
  //     }
  //   }
  // };

  class GetAssociatorKernel {
  public:
    ALPAKA_FN_ACC void operator()(
        Acc1D const& acc, 
        clue::AssociationMapView associator, 
        IndexSoA::View indexes, 
        OffsetsSoA::View offsets, 
        unsigned int begin_cands,
        unsigned int begin_clusters,
        unsigned int begin_clustered) const {
        if (once_per_grid(acc)) {
          offsets.offsets()[begin_clusters] = begin_clustered;
          for (int c_id = 0; c_id < offsets.metadata().size() - 1; c_id++) {
            auto span = associator[c_id];
            offsets.offsets()[begin_clusters + c_id + 1] = span.size() + offsets.offsets()[begin_clusters + c_id];
            auto begin = offsets.offsets()[begin_clusters + c_id];
            for (int i = 0; i < span.size(); i++) {
              indexes.indexes()[i + begin] = span[i] + begin_cands;
            }
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

    // create vector to store all association maps for each event
    std::vector<clue::AssociationMap<>> association_collection;
    association_collection.reserve(nbx);

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
      auto associator = clue_algo.getClusters(queue, points_device);

      // print stuff for control and debug
      auto nclusters = associator.size();
      auto nclustered = associator.extents().values;
      std::cout << "BX " << idx << ": found " << nclusters << " clusters from " << n_points << " PF candidates of which " << nclustered << " are clustered." << std::endl;

      // set entry of the vector
      association_collection.push_back(std::move(associator));

      // update total number of clusters and total number of clustered candidates
      clustersTotal += nclusters;
      clusteredTotal += nclustered;

      // gpetruc kernel was here
    }

    // *** ALLOCATE MULTICOLLECTION ***
    CandsClusterBxDeviceCollection ccbMap({{clusteredTotal, clustersTotal + 1, nbx, nbx + 1}}, queue);
    ccbMap.zeroInitialise();

    // * FILL BX INDEX SOA *
    alpaka::memcpy(queue, ccbMap.view<BxIndexSoA>(), bx_lookup.const_view<BxIndexSoA>());

    // * FILL BX OFFSET SOA *
    std::vector<std::size_t> bx_offsets;
    bx_offsets.reserve(association_collection.size() + 1);
    bx_offsets.push_back(0); 

    for (auto const& m : association_collection) {
        bx_offsets.push_back(bx_offsets.back() + static_cast<std::size_t>(m.size()));
    }
    
    auto bx_offsets_device = alpaka::allocAsyncBuf<data_t, Idx>(queue, Vec1D{bx_offsets.size()});
    auto bx_offsets_view = alpaka::createView(cms::alpakatools::host(), 
                                              bx_offsets.data(), 
                                              Vec1D{bx_offsets.size()});

    alpaka::memcpy(queue, bx_offsets_device, bx_offsets_view);

    // * FILL CANDIDATE INDEX SOA AND CLUSTER OFFSET SOA *
    // get initial pointers to where to insert data inside the SoAs
    auto idxsInsertPtr = ccbMap.view<IndexSoA>().indexes().data();
    auto cluOffsInsertPtr = ccbMap.view<ClusterOffsetSoA>().indexes().data() + 1; // the +1 here is fundamental

    // define begin of clusters and clustered candidates
    unsigned int begin_clusters = 0;
    unsigned int begin_clustered = 0;

    // loop over the association maps (i.e. loop over the events)
    for (unsigned int idx = 0; idx < association_collection.size(); ++idx) {
      // get current map
      auto const &asmap = association_collection[idx];

      // get offset of the candidates for the current event (! it is NOT begin_clustered)
      // this is needed for two reasons:
      // 1) to get global (i.e. orbit-wise) PF candidate indexes, since the association map
      //    returned by getClusters includes local indexes (i.e. each event starts from zero)
      // 2) as a consequence, to know which weights (pt values) to use to sort the 
      //    current association map
      const auto begin_cands = bx_lookup_host.const_view<OffsetsSoA>().offsets()[idx];

      // copy current association map into association SoA 
      // pay attention to the asmap.size() + 1 because it is very important
      auto asmap_soa = AssociationMapDevice({{asmap.extents().values, asmap.size() + 1}}, queue);
      asmap_soa.zeroInitialise(queue);

      alpaka::exec<Acc1D>(queue, 
        make_workdiv<Acc1D>(1, 1), 
        GetAssociatorKernel{}, 
        asmap.view(), 
        asmap_soa.view<IndexSoA>(), 
        asmap_soa.view<OffsetsSoA>(),
        begin_cands, 
        begin_clusters,
        begin_clustered);

      // sort asmap_soa by pt
      alpaka::exec<Acc1D>(queue, 
        make_workdiv<Acc1D>(associator.size(), 128), 
        SortClustersKernel{}, 
        pf.const_view().pt().data() + begin_cands,
        asmap_soa.view<IndexSoA>(), 
        asmap_soa.view<OffsetsSoA>());

      // create view of the buffers to be updated with the content of the current map
      auto dstIdx = alpaka::createView(alpaka::getDev(queue),
                            idxsInsertPtr, 
                            Vec1D{asmap.extents().values});                  
      auto dstOffs = alpaka::createView(alpaka::getDev(queue), 
                            offsInsertPtr,
                            Vec1D{asmap.size()}); // also here pay attention that threre is not the +1

      // create "custom" view of the offset SoA since the first element has to be discarded
      // indeed the starting point is ...data() +1
      auto srcOffs = alpaka::createView(alpaka::getDev(queue), 
                                        asmap_soa.view<OffsetsSoA>().offsets().data() + 1, 
                                        Vec1D{asmap.size()});

      // copy data from the current event asmap_soa to the global SoAs
      alpaka::memcpy(queue, dstIdx, asmap_soa.view<IndexSoA>());
      alpaka::memcpy(queue, dstOffs, srcOffs);
      
      // update pointers to where new data has to be written inside the global SoAs
      idxsInsertPtr += asmap.extents().values;
      offsInsertPtr += asmap.size(); // here there is not the +1

      // update beginners
      begin_clusters += asmap.size();
      begin_clustered += asmap.extents().values; 
    }

    return ccbMap;
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels
