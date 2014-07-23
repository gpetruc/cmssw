#ifndef _RecoTracker_DebugTools_SlidingPeakFinder_
#define _RecoTracker_DebugTools_SlidingPeakFinder_
#include <algorithm>
#include <vector>
#include <cstdint>
#include <cstdio>

namespace utils { 

class SlidingPeakFinder {
    public:
        SlidingPeakFinder(unsigned int size) :
            size_(size), half_((size+1)/2) {}

        template<typename Test>
        bool apply(const uint8_t *x, const uint8_t *begin, const uint8_t *end, const Test & test, bool verbose=false, int firststrip=0) {
            const uint8_t * ileft  = (x != begin)      ? std::min_element(x-1,x+half_)                     : begin - 1;
            const uint8_t * iright = ((x+size_) < end) ? std::min_element(x+half_,std::min(x+size_+1,end)) : end   + 1;
            uint8_t left   = (ileft  <  begin ? 0 : *ileft );
            uint8_t right =  (iright >= end   ? 0 : *iright);
            uint8_t center = *std::max_element(x, std::min(x+size_,end));
            uint8_t maxmin = std::max(left,right);
            if (maxmin < center) {
                bool ret = test(center, maxmin);
                if (ret) {
                    //static int _ndebug2 = 0;
                    //if (++_ndebug2 < 1000) {
                    //    printf(" test x = %4d, begin = %4d, end = %4d; peak size %d\n", int(x-begin)+firststrip, firststrip, int(end-begin)+firststrip, size_);
                    //    printf("   left   min[%4d,%4d] = %3d (%d) @ strip %4d\n", int((x-1)-begin)+firststrip, int(x+half_-begin)-1+firststrip, int(left), (x != begin), int(ileft-begin)+firststrip);
                    //    printf("   right  min[%4d,%4d] = %3d (%d) @ strip %4d\n", int(x+half_-begin)+firststrip, int(std::min(x+size_+1,end)-begin)-1+firststrip, int(right), ((x+half_) < end), int(iright-begin)+firststrip);
                    //    printf("   center max[%4d,%4d] = %3d\n", int(x-begin)+firststrip, int(std::min(x+size_,end)-begin)-1+firststrip, int(center));
                    //}
                    ret = test(ileft,iright,begin,end);
                }
                return ret;
            } else {
                return false;
            }
        }

        template<typename Test>
        bool apply(const std::vector<uint8_t> & ampls, const Test & test, bool verbose=false, int firststrip=0) {
            const uint8_t *begin = &*ampls.begin();
            const uint8_t *end = &*ampls.end();
            for (const uint8_t *x = begin; x < end - (half_-1); ++x) {
                if (apply(x,begin,end,test,verbose,firststrip)) {
                    return true;
                }
            }
            return false;
        }

    private:
        unsigned int size_, half_;
            
};

}

#endif
