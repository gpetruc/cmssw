#include "DataFormats/Portable/interface/PortableHostCollectionReadRules.h"
#include "DataFormats/L1ScoutingSoA/interface/OrbitEventIndexMapHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/PuppiHostCollection.h"

SET_PORTABLEHOSTCOLLECTION_READ_RULES(l1sc::OrbitEventIndexMapHostCollection);
SET_PORTABLEHOSTCOLLECTION_READ_RULES(l1sc::PuppiHostCollection);
