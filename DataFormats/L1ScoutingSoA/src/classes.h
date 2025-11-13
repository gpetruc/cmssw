#ifndef DataFormats_L1ScoutingSoA_src_classes_h
#define DataFormats_L1ScoutingSoA_src_classes_h

// these first to make sure they get included before any SoA header
#include <Eigen/Core>
#include <Eigen/Dense>

#include "DataFormats/Common/interface/Wrapper.h"

#include "DataFormats/L1ScoutingSoA/interface/IndexSoA.h"
#include "DataFormats/L1ScoutingSoA/interface/LongIndexSoA.h"
#include "DataFormats/L1ScoutingSoA/interface/BxIndexSoA.h"
#include "DataFormats/L1ScoutingSoA/interface/ClustersSoA.h"
#include "DataFormats/L1ScoutingSoA/interface/OffsetsSoA.h"
#include "DataFormats/L1ScoutingSoA/interface/ClusterIndexSoA.h"
#include "DataFormats/L1ScoutingSoA/interface/ClusterOffsetsSoA.h"
#include "DataFormats/L1ScoutingSoA/interface/PuppiSoA.h"
#include "DataFormats/L1ScoutingSoA/interface/TkEmSoA.h"
#include "DataFormats/L1ScoutingSoA/interface/TkEleSoA.h"
#include "DataFormats/L1ScoutingSoA/interface/SelectedBxSoA.h"
#include "DataFormats/L1ScoutingSoA/interface/W3PiTable.h"
#include "DataFormats/L1ScoutingSoA/interface/SoftTauTensorSoA.h"

#include "DataFormats/L1ScoutingSoA/interface/AssociationMapHost.h"
#include "DataFormats/L1ScoutingSoA/interface/LongAssociationMapHost.h"
#include "DataFormats/L1ScoutingSoA/interface/BxLookupHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/CandsClusterBxHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/ClustersHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/PuppiHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/TkEmHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/TkEleHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/SelectedBxHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/W3PiHostTable.h"
#include "DataFormats/L1ScoutingSoA/interface/CounterHost.h"
#include "DataFormats/L1ScoutingSoA/interface/SoftTauHostTensor.h"

#endif  // DataFormats_L1ScoutingSoA_src_classes_h