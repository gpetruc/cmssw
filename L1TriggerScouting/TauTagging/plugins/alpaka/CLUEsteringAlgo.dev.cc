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
        LongIndexSoA::View indexes, 
        LongOffsetsSoA::View offsets) const {
      const uint8_t kSharedMemSize = 128;
      auto& indices_shared = alpaka::declareSharedVar<int[kSharedMemSize], __COUNTER__>(acc);
      auto& weights_shared = alpaka::declareSharedVar<float[kSharedMemSize], __COUNTER__>(acc); 

      // loop over clusters in parallel
      for (auto block_idx: independent_groups(acc, offsets.metadata().size() - 1)) {
        // bind range to hw block
        auto begin = offsets.offsets()[block_idx];
        auto end = offsets.offsets()[block_idx + 1];
        // define block dimensions
        auto block_dim = end - begin;
        if (block_dim == 0)
          continue;

        // load global to shared memory with EOF sentinels
        for (auto tid : independent_group_elements(acc, kSharedMemSize)) {
          if (tid < block_dim) {
            auto thread_idx = begin + tid;
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
        for (auto i = 0; i < block_dim; i++) {
          for (auto tid : independent_group_elements(acc, block_dim - 1)) {
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
        for (auto tid : independent_group_elements(acc, block_dim)) {
          auto thread_idx = tid + begin;
          indexes.indexes()[thread_idx] = indices_shared[tid];
        }
      }
    }
  };

  // class UpdateAssociatorKernel {
  // public:
  //   ALPAKA_FN_ACC void operator()(
  //     Acc1D const& acc, 
  //     LongIndexSoA::View indexes, 
  //     LongOffsetsSoA::View offsets, 
  //     int32_t begin_indexes,
  //     int32_t begin_offsets, 
  //     int32_t num_indexes, 
  //     int32_t num_offsets) const {
  //       if (once_per_grid(acc)) {
  //         for (auto ii = 0; ii < num_indexes; ++ii) { 
  //           indexes.indexes()[ii] += begin_indexes;
  //         }
          
  //         for (auto ii = 0; ii < num_offsets; ++ii) {
  //           offsets.offsets()[ii] += begin_offsets;
  //         }
  //       }
  //     } 
  // };

  class UpdateAssociatorKernel {
  public:
    ALPAKA_FN_ACC void operator()(
      Acc1D const& acc, 
      LongIndexSoA::View indexes, 
      LongOffsetsSoA::View offsets, 
      int32_t begin_indexes,
      int32_t begin_offsets, 
      int32_t num_indexes, 
      int32_t num_offsets) const {
        for (auto thread_idx : alpaka::uniformElements(acc, num_indexes)) {
          indexes.indexes()[thread_idx] += begin_indexes;
        }

        for (auto thread_idx : alpaka::uniformElements(acc, num_offsets)) {
          offsets.offsets()[thread_idx] += begin_indexes;
        }
      } 
  };

  struct GetLastEntryKernel {
    uint32_t dst_idx; 
    uint32_t src_idx; 

    ALPAKA_FN_ACC void operator()(
      Acc1D const& acc, 
      int32_t* dst_data,
      const int32_t* src_data) const {
        if (once_per_grid(acc)) {
          auto val = src_data[src_idx];
          dst_data[dst_idx] = val;
        }
    }
  };

  CLUEsteringAlgo::CLUEsteringAlgo(float dc, float rhoc, float dm, bool wrap_coords)
      : dc_(dc), rhoc_(rhoc), dm_(dm), wrap_coords_(wrap_coords) {}

  CandsClusterBxDeviceCollection CLUEsteringAlgo::run(Queue& queue,
                            const PFCandidateDeviceCollection& pf,
                            const BxLookupDeviceCollection& bx_lookup,
                            ClustersDeviceCollection& clusters) const {
    const auto nbx = static_cast<int32_t>(bx_lookup.const_view<BxIndexSoA>().metadata().size());
    auto bx_lookup_host = BxLookupHostCollection({{nbx, nbx + 1}}, queue);
    alpaka::memcpy(queue, bx_lookup_host.buffer(), bx_lookup.buffer());
    alpaka::wait(queue);

    // create vector to store all association maps for each event
    std::vector<clue::AssociationMap<>> association_collection;

    // create host vector where to store the TRUE number of clustered candidates per event
    // and also create an index to fill it
    auto num_clustered_d = alpaka::allocAsyncBuf<int32_t, int32_t>(queue, nbx);
    // keep track of the total number of clusters in the current orbit
    int32_t clustersTotal = 0;

    for (uint32_t idx = 0; idx < bx_lookup_host.const_view<BxIndexSoA>().metadata().size(); idx++) {

#if defined(__DEBUG__) || defined(__DEBUGLITE__)
      std::cout << "\n\n------------------------- ITERATION INDEX " << idx << "------------------------------" << std::endl;
#endif

      const auto begin = bx_lookup_host.const_view<OffsetsSoA>().offsets()[idx];
      const auto end = bx_lookup_host.const_view<OffsetsSoA>().offsets()[idx + 1];
      const uint32_t n_points = end - begin;

      if (n_points == 0) {
#if defined(__DEBUG__) || defined(__DEBUGLITE__)
        std::cout << "WARNING! NO POINTS FOUND IN THE EVENT" << std::endl;
#endif
        throw std::runtime_error("Encountered an event with no candidates inside");
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

      if (associator.size() == 0) {
#if defined(__DEBUG__) || defined(__DEBUGLITE__)
        std::cout << "WARNING! associatior.size()IN THE EVENT" << std::endl;
        // maybe here it could be necessary to handcraft a "mock" AssociationMap
#endif
      }
      
      // fetch true number of clustered points from the buffers of the association map
      auto keys_size = alpaka::getExtents(associator.extract().keys);
      alpaka::exec<Acc1D>(queue,
                          make_workdiv<Acc1D>(1, 1), 
                          GetLastEntryKernel{idx, keys_size[0] - 1}, 
                          num_clustered_d.data(),
                          associator.extract().keys.data());
      
      // set entry of the association map vector
      association_collection.push_back(std::move(associator));
      // update total number of clusters in the current orbit
      clustersTotal += associator.size();
    }

    // bring number of clustered candidates per event back to the host
    std::vector<int32_t> num_clustered(nbx);
    alpaka::memcpy(queue, num_clustered, num_clustered_d);

    // accumulate over the entries of the num_clustered vector in order to get the 
    // total number of clustered candidates in the current orbit
    int32_t clusteredTotal = std::accumulate(num_clustered.begin(), num_clustered.end(), 0);

#if defined(__DEBUG__) || defined(__DEBUGLITE__)
    std::cout << "\n\n\nTOTAL CLUSTERED: " << clusteredTotal << " TOTAL CLUSTERS: " << clustersTotal << std::endl;
    std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
#endif
    
    // *** ALLOCATE MULTICOLLECTION ***
    CandsClusterBxDeviceCollection ccbMap({{clusteredTotal, clustersTotal + 1, clustersTotal, nbx + 1}}, queue);
    //                                      LongIndexSoA    ClusterOffsetSoA   ClusterIndexSoA LongOffsetSoA
    ccbMap.zeroInitialise(queue);

    // * FILL CLUSTER INDEX SOA *
    std::vector<int32_t> cluster_indexes(clustersTotal);
    std::iota(cluster_indexes.begin(), cluster_indexes.end(), 0);

    auto srcClusterIndex = alpaka::createView(cms::alpakatools::host(), 
                                              cluster_indexes.data(), 
                                              Vec1D{cluster_indexes.size()});
    
    auto dstClusterIndex = alpaka::createView(alpaka::getDev(queue), 
                                              ccbMap.view<ClusterIndexSoA>().indexes().data(), 
                                              Vec1D{clustersTotal});
                                                    
    alpaka::memcpy(queue, dstClusterIndex, srcClusterIndex);

    // * FILL LAST LONG OFFSET SOA *
    std::vector<int32_t> bx_offsets;
    bx_offsets.reserve(association_collection.size() + 1);
    bx_offsets.push_back(0); 
    
    for (auto const& m : association_collection) {
      bx_offsets.push_back(bx_offsets.back() + m.size());
    }
    
    auto dstClusterOffsets = alpaka::createView(alpaka::getDev(queue), 
                                                ccbMap.view<LongOffsetsSoA>().offsets().data(), 
                                                Vec1D{nbx + 1});

    auto srcClusterOffsets = alpaka::createView(cms::alpakatools::host(), 
                                                bx_offsets.data(), 
                                                Vec1D{bx_offsets.size()});

    alpaka::memcpy(queue, dstClusterOffsets, srcClusterOffsets);
    
    // * FILL FIRST LONG INDEX SOA AND CLUSTER OFFSET SOA *
    // get initial pointers to where to insert data inside the SoAs
    auto idxsInsertPtr = ccbMap.view<LongIndexSoA>().indexes().data();
    auto offsInsertPtr = ccbMap.view<ClusterOffsetsSoA>().offsets().data() + 1; // the +1 here is fundamental
    
    // define begin of clusters and clustered candidates
    int32_t begin_offsets = 0;
    
    // loop over the association maps (i.e. loop over the events)
    for (uint32_t idx = 0; idx < association_collection.size(); ++idx) {

#if defined(__DEBUG__) || defined(__DEBUGLITE__)
      std::cout << "\n\n------------------------- ITERATION INDEX " << idx << "------------------------------" << std::endl;
#endif

      // get current map
      auto const &asmap = association_collection[idx];
      
      // get offset of the candidates for the current event (! it is NOT begin_clustered)
      // this is needed for two reasons:
      // 1) to get global (i.e. orbit-wise) PF candidate indexes, since the association map
      //    returned by getClusters includes local indexes (i.e. each event starts from zero)
      // 2) as a consequence, to know which weights (pt values) to use to sort the 
      //    current association map
      const auto begin_indexes = static_cast<int32_t>(bx_lookup_host.const_view<OffsetsSoA>().offsets()[idx]);

#if defined(__DEBUG__) || defined(__DEBUGLITE__)
      std::cout << "NOW CREATING LOCAL AssociationMapDevice WITH INDEXES LENGTH = " << num_clustered[idx]
                << " OFFSETS LENGTH = " << asmap.size() + 1 << std::endl;
#endif

      // copy current association map into association SoA 
      // pay attention to the asmap.size() + 1 because it is very important
      auto asmap_soa = LongAssociationMapDevice({{num_clustered[idx], 
                                                  static_cast<int32_t>(asmap.size() + 1)}}, 
                                                  queue);
      asmap_soa.zeroInitialise(queue);

      // copy local indexes inside local association map
      auto dstLocalIdx = alpaka::createView(alpaka::getDev(queue), 
                                      asmap_soa.view<LongIndexSoA>().indexes().data(), 
                                      Vec1D{asmap_soa.view<LongIndexSoA>().metadata().size()});

      auto srcLocalIdx = alpaka::createView(alpaka::getDev(queue), // here creating srcLocalIdx is needed since we need to copy only the TRUE number of indexes
                                      asmap.extract().values.data(), 
                                      Vec1D{num_clustered[idx]});
        
      alpaka::memcpy(queue, dstLocalIdx, srcLocalIdx);  
      
      // copy local offsets inside local association map
      auto dstLocalOffs = alpaka::createView(alpaka::getDev(queue), 
                                      asmap_soa.view<LongOffsetsSoA>().offsets().data(), 
                                      Vec1D{asmap_soa.view<LongOffsetsSoA>().metadata().size()});

      alpaka::memcpy(queue, dstLocalOffs, asmap.extract().keys);

      // sort local association map
      alpaka::exec<Acc1D>(queue, 
        make_workdiv<Acc1D>(asmap.size(), 128), 
        SortClustersKernel{}, 
        pf.const_view().pt().data() + begin_indexes, // importanto to use begin indexes
        asmap_soa.view<LongIndexSoA>(), 
        asmap_soa.view<LongOffsetsSoA>());
      
      // update the indexes of the local association map in order to
      // prepare it to be copied to the global association map
      alpaka::exec<Acc1D>(queue, 
        make_workdiv<Acc1D>(1, 1), 
        UpdateAssociatorKernel{}, 
        asmap_soa.view<LongIndexSoA>(), 
        asmap_soa.view<LongOffsetsSoA>(),
        begin_indexes, 
        begin_offsets,
        asmap_soa.view<LongIndexSoA>().metadata().size(), 
        asmap_soa.view<LongOffsetsSoA>().metadata().size());

#if defined(__DEBUG__)
      std::cout << "begin_indexes: " << begin_indexes << std::endl;
      
      std::cout << "Original Indexes" << std::endl;                         
      std::vector<int32_t> debug_clue_indexes(asmap.extents().values);
      alpaka::memcpy(queue, debug_clue_indexes, asmap.extract().values);
      
      for (auto el : debug_clue_indexes) 
      std::cout << el << " : ";
      std::cout << "\n";
      
      std::cout << "Indexes (sorted) copied inside the local association map" << std::endl;
      std::vector<int32_t> debug_soa_indexes(asmap_soa.view<LongIndexSoA>().metadata().size());
      
      auto debugSrc = alpaka::createView(alpaka::getDev(queue),
                                    asmap_soa.view<LongIndexSoA>().indexes().data(), 
                                    Vec1D{asmap_soa.view<LongIndexSoA>().metadata().size()});
      
      alpaka::memcpy(queue, debug_soa_indexes, debugSrc);

      for (auto el : debug_soa_indexes) 
        std::cout << el << " : ";
      std::cout << "\n";
      
      std::cout << "begin_offsets: " << begin_offsets << std::endl;

      std::cout << "Original Offsets" << std::endl;  
      std::vector<int32_t> debug_clue_offsets(asmap.extents().keys);
      alpaka::memcpy(queue, debug_clue_offsets, asmap.extract().keys);

      for (auto el : debug_clue_offsets) 
        std::cout << el << " : ";
      std::cout << "\n";
      
      std::cout << "Offsets copied inside the local association map" << std::endl;
      std::vector<int32_t> debug_soa_offsets(asmap_soa.view<LongOffsetsSoA>().metadata().size());
      debugSrc = alpaka::createView(alpaka::getDev(queue),
                                    asmap_soa.view<LongOffsetsSoA>().offsets().data(), 
                                    Vec1D{asmap_soa.view<LongOffsetsSoA>().metadata().size()}); 
      
      alpaka::memcpy(queue, debug_soa_offsets, debugSrc);

      for (auto el : debug_soa_offsets) 
        std::cout << el << " : ";
      std::cout << "\n";
#endif

      // create view of the buffers to be updated with the content of the current map
      auto dstIdx = alpaka::createView(alpaka::getDev(queue),
                            idxsInsertPtr, 
                            Vec1D{asmap_soa.view<LongIndexSoA>().metadata().size()});      
      
      auto dstOffs = alpaka::createView(alpaka::getDev(queue), 
                            offsInsertPtr,
                            Vec1D{asmap_soa.view<LongOffsetsSoA>().metadata().size() - 1});
      
      // create view of the buffers to be copied
      auto srcIdx = alpaka::createView(alpaka::getDev(queue),
                                      asmap_soa.view<LongIndexSoA>().indexes().data(), 
                                      Vec1D{asmap_soa.view<LongIndexSoA>().metadata().size()});                          
      
      // create "custom" view of the offset SoA since the first element has to be discarded
      // indeed the starting point is ...data() +1
      auto srcOffs = alpaka::createView(alpaka::getDev(queue), 
                                        asmap_soa.view<LongOffsetsSoA>().offsets().data() + 1, 
                                        Vec1D{asmap_soa.view<LongOffsetsSoA>().metadata().size() - 1});
      
      alpaka::memcpy(queue, dstIdx, srcIdx);
      alpaka::memcpy(queue, dstOffs, srcOffs);
      

      // update pointers to where new data has to be written inside the global SoAs
      idxsInsertPtr += num_clustered[idx];
      offsInsertPtr += asmap.size(); // here there is not the + 1

      // update beginners
      begin_offsets += num_clustered[idx];
    }

    return ccbMap;
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels
