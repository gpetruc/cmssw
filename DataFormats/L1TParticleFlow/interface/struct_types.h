#ifndef DataFormats_L1TParticleFlow_struct_types_h
#define DataFormats_L1TParticleFlow_struct_types_h

#include <vector>
#include <utility>
#include <cstdint>
#include <boost/align/aligned_allocator.hpp>

namespace l1ct {
    namespace structs {
        struct Puppi {
            float pt, eta, phi, z0, dxy, puppiw;
            int16_t pdgId;
            uint8_t quality;
        };
        struct PuppiSOA {
            std::vector<uint16_t> bx;
            std::vector<uint32_t> offsets;
            std::vector<float> pt, eta, phi, z0, dxy, puppiw;
            std::vector<int16_t> pdgId;
            std::vector<uint8_t> quality;
            PuppiSOA() : bx(), offsets(), pt(), eta(), phi(), z0(), dxy(), puppiw(), pdgId(), quality() {}
            PuppiSOA(const PuppiSOA& other) = default;
            PuppiSOA(PuppiSOA&& other) = default;
            PuppiSOA & operator=(const PuppiSOA& other) = default;
            PuppiSOA & operator=(PuppiSOA&& other) = default;
            void swap(PuppiSOA& other) {
              using std::swap;
              swap(bx, other.bx);
              swap(offsets, other.offsets);
              swap(pt, other.pt);
              swap(eta, other.eta);
              swap(phi, other.phi);
              swap(z0, other.z0);
              swap(dxy, other.dxy);
              swap(puppiw, other.puppiw);
              swap(pdgId, other.pdgId);
              swap(quality, other.quality);
            }
        };
        inline void swap(PuppiSOA & a, PuppiSOA & b) {
            a.swap(b);
        }
        struct PuppiASOA {
          template <class T>
          using avector = std::vector<T, boost::alignment::aligned_allocator<T, 64>>;
          avector<uint16_t> bx;
          avector<uint32_t> offsets;
          avector<float> pt, eta, phi, z0, dxy, puppiw;
          avector<int16_t> pdgId;
          avector<uint8_t> quality;
          PuppiASOA() : bx(), offsets(), pt(), eta(), phi(), z0(), dxy(), puppiw(), pdgId(), quality() {}
          PuppiASOA(const PuppiASOA& other) = default;
          PuppiASOA(PuppiASOA&& other) = default;
          PuppiASOA& operator=(const PuppiASOA& other) = default;
          PuppiASOA& operator=(PuppiASOA&& other) = default;
          void swap(PuppiASOA& other) {
            using std::swap;
            swap(bx, other.bx);
            swap(offsets, other.offsets);
            swap(pt, other.pt);
            swap(eta, other.eta);
            swap(phi, other.phi);
            swap(z0, other.z0);
            swap(dxy, other.dxy);
            swap(puppiw, other.puppiw);
            swap(pdgId, other.pdgId);
            swap(quality, other.quality);
            }
        };
        inline void swap(PuppiASOA & a, PuppiASOA & b) {
            a.swap(b);
        }
    }
}
#endif
