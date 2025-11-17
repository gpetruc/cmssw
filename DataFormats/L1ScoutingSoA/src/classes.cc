#include "DataFormats/L1ScoutingSoA/src/classes.h"
#include "DataFormats/Portable/interface/PortableHostCollectionReadRules.h"
#include "DataFormats/Portable/interface/PortableHostObjectReadRules.h"
#include "DataFormats/L1ScoutingSoA/interface/AssociationMapHost.h"
#include "DataFormats/L1ScoutingSoA/interface/BxLookupHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/ClustersHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/PuppiHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/TkEmHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/TkEleHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/SelectedBxHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/W3PiHostTable.h"
#include "DataFormats/L1ScoutingSoA/interface/CounterHost.h"
#include "DataFormats/L1ScoutingSoA/interface/SoftTauHostTensor.h"

SET_PORTABLEHOSTMULTICOLLECTION_READ_RULES(l1sc::AssociationMapHost);
SET_PORTABLEHOSTMULTICOLLECTION_READ_RULES(l1sc::BxLookupHostCollection);
SET_PORTABLEHOSTCOLLECTION_READ_RULES(l1sc::ClustersHostCollection);
SET_PORTABLEHOSTCOLLECTION_READ_RULES(l1sc::ClusterObjHostCollection);
SET_PORTABLEHOSTCOLLECTION_READ_RULES(l1sc::PuppiHostCollection);
SET_PORTABLEHOSTCOLLECTION_READ_RULES(l1sc::TkEmHostCollection);
SET_PORTABLEHOSTCOLLECTION_READ_RULES(l1sc::TkEleHostCollection);
SET_PORTABLEHOSTCOLLECTION_READ_RULES(l1sc::SelectedBxHostCollection);
SET_PORTABLEHOSTCOLLECTION_READ_RULES(l1sc::W3PiHostTable);
SET_PORTABLEHOSTCOLLECTION_READ_RULES(l1sc::SoftTauInputHostTensor);
SET_PORTABLEHOSTCOLLECTION_READ_RULES(l1sc::SoftTauOutputHostTensor);
SET_PORTABLEHOSTOBJECT_READ_RULES(l1sc::CounterHost);
