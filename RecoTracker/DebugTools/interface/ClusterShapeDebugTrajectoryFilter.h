#ifndef _ClusterShapeDebugTrajectoryFilter_h_
#define _ClusterShapeDebugTrajectoryFilter_h_

#include <vector>
#include <unordered_map>
#include "TrackingTools/TrajectoryFiltering/interface/TrajectoryFilter.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class ClusterShapeHitFilter;
class TrackerTopology;
class TrackerGeometry;
class TrackerHitAssociator;
class TrajectoryMeasurement;
class TrajectoryStateOnSurface;
class MeasurementTrackerEvent;
class TTree;
namespace edm { class Event; class EventSetup; }


class ClusterShapeDebugTrajectoryFilter : public TrajectoryFilter {
 public:
  //  ClusterShapeDebugTrajectoryFilter(const edm::EventSetup& es);

  ClusterShapeDebugTrajectoryFilter(const edm::ParameterSet &iConfig);

  virtual ~ClusterShapeDebugTrajectoryFilter();

  virtual bool qualityFilter(const TempTrajectory&) const;
  virtual bool qualityFilter(const Trajectory&) const;
 
  virtual bool toBeContinued(TempTrajectory&) const;
  virtual bool toBeContinued(Trajectory&) const;

  virtual std::string name() const { return "ClusterShapeDebugTrajectoryFilter"; }

  void init(const edm::Event &, const edm::EventSetup &) const ;
  void done() const;

 private:
  edm::ParameterSet              theConfig;
  mutable TrackerHitAssociator * theAssociator;
  mutable const ClusterShapeHitFilter   * theFilter;
  mutable const TrackerTopology         * theTopology;
  mutable const TrackerGeometry         * theTracker;
  mutable const MeasurementTrackerEvent * theMTEvent;

  mutable TTree * theTree;
  mutable float trackPt_, trackEta_; 
  mutable int   trackAssoc_;
  mutable int   hitDet_, hitLayer_, hitAssoc_, hitSingleSim_, hitUsable_, hitStrips_, hitCompatPos_, hitCompatNoPos_, hitBadBoundaries_;
  mutable float hitCharge_,  hitChargeNorm_, hitChargeRms_;
  mutable float hitPredPos_, hitPredNoPos_;

  typedef std::pair<unsigned int, unsigned int> Id2;
  void fillAssociations(const TrackingRecHit *hit, std::vector<Id2> &out) const  ;
  void fillCluster(const TrackingRecHit *hit, const TrajectoryStateOnSurface &tsos, Id2 simtk, bool mustProject=false) const ;
  bool hasBadStrip(unsigned int detid, int strip) const ;
  mutable std::unordered_map<unsigned int, unsigned int> theStripDetLookup;
  void    indexStripDets() const ;
};

#endif
