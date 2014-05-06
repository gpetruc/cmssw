#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "RecoTracker/DebugTools/interface/ClusterShapeDebugTrajectoryFilter.h"

class ClusterShapeDebugTrajectoryFilterUpdator : public edm::EDAnalyzer {
   public:
      ClusterShapeDebugTrajectoryFilterUpdator(const edm::ParameterSet &iConfig) :
         label_(iConfig.getParameter<std::string>("label")),
         isInit_(iConfig.getParameter<std::string>("op") == "init")
      {
         if (!isInit_ && iConfig.getParameter<std::string>("op") != "done") {
            throw cms::Exception("Configuration", "op must be 'init' or 'done'");
         }
      }

      ~ClusterShapeDebugTrajectoryFilterUpdator() 
      {
      }

   private:
      virtual void analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup) override
      {
         iSetup.get<TrajectoryFilter::Record>().get(label_, handle);
         const ClusterShapeDebugTrajectoryFilter & filter = dynamic_cast<const ClusterShapeDebugTrajectoryFilter &>(*handle);
         if (isInit_) filter.init(iEvent, iSetup);
         else         filter.done(); 
      }

      std::string label_;
      bool        isInit_;
      edm::ESHandle<TrajectoryFilter> handle;
};

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(ClusterShapeDebugTrajectoryFilterUpdator);
