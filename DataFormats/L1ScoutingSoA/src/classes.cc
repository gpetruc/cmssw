#include "DataFormats/Portable/interface/PortableHostCollectionReadRules.h"
#include "DataFormats/Portable/interface/PortableHostObjectReadRules.h"
#include "DataFormats/L1ScoutingSoA/interface/BxLookupHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/ClustersHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/PuppiHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/PFCandidateHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/SelectedBxHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/W3PiHostTable.h"
#include "DataFormats/L1ScoutingSoA/interface/CounterHost.h"

SET_PORTABLEHOSTMULTICOLLECTION_READ_RULES(l1sc::BxLookupHostCollection);
SET_PORTABLEHOSTCOLLECTION_READ_RULES(l1sc::ClustersHostCollection);
SET_PORTABLEHOSTCOLLECTION_READ_RULES(l1sc::PuppiHostCollection);
SET_PORTABLEHOSTCOLLECTION_READ_RULES(l1sc::SelectedBxHostCollection);
SET_PORTABLEHOSTCOLLECTION_READ_RULES(l1sc::W3PiHostTable);
SET_PORTABLEHOSTCOLLECTION_READ_RULES(l1sc::PFCandidateHostCollection);
SET_PORTABLEHOSTOBJECT_READ_RULES(l1sc::CounterHost);
