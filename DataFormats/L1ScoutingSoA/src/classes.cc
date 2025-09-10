#include "DataFormats/Portable/interface/PortableHostCollectionReadRules.h"
#include "DataFormats/Portable/interface/PortableHostObjectReadRules.h"
#include "DataFormats/L1ScoutingSoA/interface/OrbitEventIndexMapHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/PuppiHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/CounterHost.h"

SET_PORTABLEHOSTCOLLECTION_READ_RULES(l1sc::OrbitEventIndexMapHostCollection);
SET_PORTABLEHOSTCOLLECTION_READ_RULES(l1sc::PuppiHostCollection);
SET_PORTABLEHOSTOBJECT_READ_RULES(l1sc::CounterHost);
