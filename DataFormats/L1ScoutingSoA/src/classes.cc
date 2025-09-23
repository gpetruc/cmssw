#include "DataFormats/Portable/interface/PortableHostCollectionReadRules.h"
#include "DataFormats/Portable/interface/PortableHostObjectReadRules.h"
#include "DataFormats/L1ScoutingSoA/interface/HostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/HostObject.h"
#include "DataFormats/L1ScoutingSoA/interface/CounterHost.h"

SET_PORTABLEHOSTCOLLECTION_READ_RULES(l1sc::PuppiHostCollection);
SET_PORTABLEHOSTMULTICOLLECTION_READ_RULES(l1sc::NbxMapHostCollection);
SET_PORTABLEHOSTOBJECT_READ_RULES(l1sc::W3PiTripletHostObject);
SET_PORTABLEHOSTCOLLECTION_READ_RULES(l1sc::W3PiPuppiTableHostCollection);
SET_PORTABLEHOSTOBJECT_READ_RULES(l1sc::CounterHost);
