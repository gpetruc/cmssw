#ifndef BadComponents_h
#define BadComponents_h

#include <cstdint>
#include <vector>
#include <array>
#include <unordered_map>
#include <string>
#include <bitset>

namespace impl {
    struct BadSingleComponent {
        std::vector<bool> lumi;
        // return true if lumi[i] is true for i any ls, ls+1
        bool test(unsigned int ls) const {
            if (ls-1 < lumi.size()) return lumi[ls] || lumi[ls+1];
            else if (ls < lumi.size()) return lumi[ls];
            else return false;
        }
   };

    struct BadModule {
        uint32_t detid;
        std::array<BadSingleComponent,4>   tbm;
        std::array<BadSingleComponent,16> roc;
        BadModule(uint32_t id = 0) : detid(id) { 
            for (auto & c : tbm) c = BadSingleComponent();
            for (auto & c : roc) c = BadSingleComponent();
        }
        void fetch(unsigned int lumi, std::vector<int> & rocs) const {
            rocs.clear();
            std::bitset<16> byroc(false);
            for (int itbm = 0; itbm < 4; ++itbm) {
                if (tbm[itbm].test(lumi)) { 
                    for (int i = 0; i < 4; ++i) byroc[4*itbm+i] = true;
                }
            }
            for (int iroc = 0; iroc < 16; ++iroc) {
                if (!byroc[iroc] &&  roc[iroc].test(lumi)) byroc[iroc] = true;
            }
            for (int iroc = 0; iroc < 16; ++iroc) {
                if (byroc[iroc]) rocs.push_back(iroc);
            }
        }
    };

    struct BadRun {
        std::unordered_map<unsigned,BadModule> modules;

        void fill(uint32_t detid, int tbm, int roc, unsigned int lumi) {
            auto &mod = modules[detid];
            if (mod.detid != detid) mod.detid = detid;
            auto & ls = (tbm == -1) ? mod.roc[roc].lumi : mod.tbm[tbm].lumi;
            if (ls.size() <= lumi) ls.resize(std::max<unsigned>(lumi+1, ls.size()*2), false);
            ls[lumi] = true;
        }

        void fetch(uint32_t detid, unsigned int lumi, std::vector<int> & rocs) const {
            auto it = modules.find(detid);
            if (it != modules.end()) it->second.fetch(lumi,rocs);
        }
    };
}

class BadComponents {
    public:
        BadComponents() {}
        void init(const std::string &filename) ;
        void fetch(int run, int lumi, uint32_t detid, std::vector<int> & badrocs) const {
            auto it = byrun.find(run);
            if (it != byrun.end()) it->second.fetch(detid, lumi, badrocs);
        }
    private:
        std::unordered_map<unsigned,impl::BadRun> byrun;
        
};

#endif
