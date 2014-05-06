// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "RecoTracker/DebugTools/interface/ClusterShapeDebugTrajectoryFilter.h"

//
// class decleration
//

class ClusterShapeDebugTrajectoryFilterESProducer : public edm::ESProducer
{
 public:
  ClusterShapeDebugTrajectoryFilterESProducer(const edm::ParameterSet&);
  ~ClusterShapeDebugTrajectoryFilterESProducer();

  typedef std::auto_ptr<TrajectoryFilter> ReturnType;

  ReturnType produce(const TrajectoryFilter::Record &);

 private:
  std::string componentName;
  edm::ParameterSet filterPSet;
};

ClusterShapeDebugTrajectoryFilterESProducer::ClusterShapeDebugTrajectoryFilterESProducer
  (const edm::ParameterSet& iConfig) :
  filterPSet(iConfig)
{
  componentName = iConfig.getParameter<std::string>("ComponentName");
  
  setWhatProduced(this, componentName);
}


/*****************************************************************************/
ClusterShapeDebugTrajectoryFilterESProducer::~ClusterShapeDebugTrajectoryFilterESProducer
  ()
{
}

/*****************************************************************************/
ClusterShapeDebugTrajectoryFilterESProducer::ReturnType
ClusterShapeDebugTrajectoryFilterESProducer::produce
(const TrajectoryFilter::Record &iRecord)
{
  ReturnType aFilter(new ClusterShapeDebugTrajectoryFilter(filterPSet));
  return aFilter;
}

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/typelookup.h"
DEFINE_FWK_EVENTSETUP_MODULE(ClusterShapeDebugTrajectoryFilterESProducer);
