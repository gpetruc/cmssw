#include "L1TriggerScouting/TauTagging/plugins/alpaka/CLUEsteringAlgo.h"

#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels {

  using namespace cms::alpakatools;

  CLUEsteringAlgo::CLUEsteringAlgo(float dc, float rhoc, float dm, bool wrap_coords)
      : dc_(dc), rhoc_(rhoc), dm_(dm), wrap_coords_(wrap_coords) {}

  void CLUEsteringAlgo::run(Queue& queue,
                            const PFCandidateDeviceCollection& pf,
                            ClustersDeviceCollection& clusters) const {
    // CLUEstering call internally reinterpret_cast<T*> to non-const ptr
    auto& pf_cast = const_cast<PFCandidateDeviceCollection&>(pf);
    const uint32_t n_points = pf_cast.view().metadata().size();

    // buffers
    auto* coords_ptr = pf_cast.view().eta();
    auto* weights_ptr = pf_cast.view().pt();
    auto* cluster_indices_ptr = clusters.view().cluster();
    auto* is_seed_ptr = clusters.view().is_seed();

    // wrap device buffers
    auto points_device =
        clue::PointsDevice<kDims, Device>(queue, n_points, coords_ptr, weights_ptr, cluster_indices_ptr, is_seed_ptr);

    // run (wrap coords if enabled)
    auto clue_algo = clue::Clusterer<kDims>(queue, dc_, rhoc_, dm_);
    if (wrap_coords_)
      clue_algo.setWrappedCoordinates({{0, 1}});
    clue_algo.make_clusters(queue, points_device);
  }

  void CLUEsteringAlgo::run(Queue& queue,
                            const PFCandidateDeviceCollection& pf,
                            const BxLookupDeviceCollection& bx_lookup,
                            ClustersDeviceCollection& clusters) const {
    // CLUEstering call internally reinterpret_cast<T*> to non-const ptr, alpaka::memcpy also
    auto& pf_cast = const_cast<PFCandidateDeviceCollection&>(pf);
    auto& bx_lookup_cast = const_cast<BxLookupDeviceCollection&>(bx_lookup);

    const auto nbx = static_cast<int32_t>(bx_lookup_cast.const_view<BxIndexSoA>().metadata().size());
    auto bx_lookup_host = BxLookupHostCollection({{nbx, nbx + 1}}, queue);
    alpaka::memcpy(queue, bx_lookup_host.buffer(), bx_lookup_cast.buffer());
    alpaka::wait(queue);

    for (int32_t idx = 0; idx < bx_lookup_host.const_view<BxIndexSoA>().metadata().size(); idx++) {
      const auto begin = bx_lookup_host.const_view<OffsetsSoA>().offsets()[idx];
      const auto end = bx_lookup_host.const_view<OffsetsSoA>().offsets()[idx + 1];
      const uint32_t n_points = end - begin;

      if (n_points == 0)
        continue;

      // types
      using coords_dtype_t = std::remove_pointer_t<decltype(pf_cast.view().eta())>;
      using weights_dtype_t = std::remove_pointer_t<decltype(pf_cast.view().pt())>;
      using clusters_dtype_t = std::remove_pointer_t<decltype(clusters.view().cluster())>;
      using seeds_dtype_t = std::remove_pointer_t<decltype(clusters.view().is_seed())>;

      // base extent
      auto extent = Vec1D{n_points};

      // CLUEstering expects contiguous memory buffers which cannot be guaranteed
      // by the currect scouting SOA design, so slow D2D copy is done instead.
      // This should be optimized along with the talks with the CLUEstering developers (batching support)
      auto coords_buffer = make_device_buffer<coords_dtype_t[]>(queue, n_points * kDims);

      // spans to memory layouts
      std::span<coords_dtype_t> eta_span(pf_cast.view().eta() + begin, pf_cast.view().eta() + end);
      std::span<coords_dtype_t> phi_span(pf_cast.view().phi() + begin, pf_cast.view().phi() + end);
      std::span<weights_dtype_t> pt_span(pf_cast.view().pt() + begin, pf_cast.view().pt() + end);
      std::span<clusters_dtype_t> cluster_span(clusters.view().cluster() + begin, clusters.view().cluster() + end);
      std::span<seeds_dtype_t> seed_span(clusters.view().is_seed() + begin, clusters.view().is_seed() + end);

      // memcpy
      auto eta_coord = createView(alpaka::getDev(queue), alpaka::getPtrNative(coords_buffer), extent);
      auto phi_coord = createView(alpaka::getDev(queue), alpaka::getPtrNative(coords_buffer) + n_points, extent);
      alpaka::memcpy(queue, eta_coord, createView(alpaka::getDev(queue), eta_span.data(), extent));
      alpaka::memcpy(queue, phi_coord, createView(alpaka::getDev(queue), phi_span.data(), extent));

      // // buffers
      auto* coords_ptr = alpaka::getPtrNative(coords_buffer);
      auto* weights_ptr = pt_span.data();
      auto* cluster_indices_ptr = cluster_span.data();
      auto* is_seed_ptr = seed_span.data();

      // wrap device buffers
      auto points_device =
          clue::PointsDevice<kDims, Device>(queue, n_points, coords_ptr, weights_ptr, cluster_indices_ptr, is_seed_ptr);

      // run (wrap coords if enabled)
      auto clue_algo = clue::Clusterer<kDims>(queue, dc_, rhoc_, dm_);
      if (wrap_coords_)
        clue_algo.setWrappedCoordinates({{0, 1}});
      clue_algo.make_clusters(queue, points_device);
    }
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels