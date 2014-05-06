#ifndef RecoPixelVertexing_PixelLowPtUtilities_ClusterShapeDebugTrajectoryFilterESProducer_H
#define RecoPixelVertexing_PixelLowPtUtilities_ClusterShapeDebugTrajectoryFilterESProducer_H


// -*- C++ -*-
//
// Package:    ClusterShapeDebugTrajectoryFilterESProducer
// Class:      ClusterShapeDebugTrajectoryFilterESProducer
// 
/**\class ClusterShapeDebugTrajectoryFilterESProducer ClusterShapeDebugTrajectoryFilterESProducer.h TrackingTools/ClusterShapeDebugTrajectoryFilterESProducer/src/ClusterShapeDebugTrajectoryFilterESProducer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Jean-Roch Vlimant
//         Created:  Fri Sep 28 18:07:52 CEST 2007
//
//


// system include files
#include <memory>
#include "boost/shared_ptr.hpp"

// user include files
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Framework/interface/ESProducer.h"

#include "FWCore/Framework/interface/ESHandle.h"

//#include "RecoTracker/DebugTools/interface/ClusterShapeDebugTrajectoryFilter.h"

//class TrajectoryFilter;
#include "TrackingTools/TrajectoryFiltering/interface/TrajectoryFilter.h"

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
  std::string componentType;
  edm::ParameterSet filterPset;
};

#endif
