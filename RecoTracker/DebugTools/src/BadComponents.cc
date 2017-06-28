#include "../interface/BadComponents.h"
#include "FWCore/Utilities/interface/Exception.h"

#include <cstdlib>
#include <cstdio>
void BadComponents::init(const std::string &filename) {
    FILE *fin = fopen(filename.c_str(), "r");
    if (!fin) throw cms::Exception("Cannot read input filename", filename);
    unsigned int run, lumi, detid;
    int tbm, roc;
    unsigned int rows = 0;
    while (fscanf(fin, "%u %u %u %d %d", &run, &lumi, &detid, &tbm, &roc) == 5) {
        byrun[run].fill(detid,tbm,roc,lumi);
        rows++;
    }
    printf("Read %d records\n",rows);
}
