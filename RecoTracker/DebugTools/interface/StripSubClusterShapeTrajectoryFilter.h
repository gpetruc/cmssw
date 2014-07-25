#ifndef _StripSubClusterShapeTrajectoryFilter_h_
#define _StripSubClusterShapeTrajectoryFilter_h_

#include <vector>
#include <unordered_map>
#include "TrackingTools/TrajectoryFiltering/interface/TrajectoryFilter.h"
#include "RecoTracker/MeasurementDet/interface/MeasurementTrackerEvent.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class ClusterShapeHitFilter;
class TrackerTopology;
class TrackerGeometry;
class TrajectoryMeasurement;
class TrajectoryStateOnSurface;
class MeasurementTrackerEvent;
class SiStripNoises;
class TTree;
namespace edm { class Event; class EventSetup; class ConsumesCollector; }


class StripSubClusterShapeTrajectoryFilter : public TrajectoryFilter {
 public:
  //  StripSubClusterShapeTrajectoryFilter(const edm::EventSetup& es);

  StripSubClusterShapeTrajectoryFilter(const edm::ParameterSet &iConfig, edm::ConsumesCollector& iC);

  virtual ~StripSubClusterShapeTrajectoryFilter();

  virtual bool qualityFilter(const TempTrajectory&) const override;
  virtual bool qualityFilter(const Trajectory&) const override;
 
  virtual bool toBeContinued(TempTrajectory&) const override;
  virtual bool toBeContinued(Trajectory&) const override;

  virtual std::string name() const { return "StripSubClusterShapeTrajectoryFilter"; }

  virtual void setEvent(const edm::Event &, const edm::EventSetup &) override;

 protected:

  virtual bool testLastHit(const TrackingRecHit *hit, const TrajectoryStateOnSurface &tsos, bool mustProject=false) const ;

  // pass-through of clusters with too many consecutive saturated strips
  uint32_t maxNSat_;

  // trimming parameters
  uint8_t trimMaxADC_;
  float   trimMaxFracTotal_, trimMaxFracNeigh_;

  // maximum difference after peak finding
  float   maxTrimmedSizeDiffPos_, maxTrimmedSizeDiffNeg_;

  // peak finding parameters
  float subclusterWindow_;
  float seedCutMIPs_, seedCutSN_;
  float subclusterCutMIPs_, subclusterCutSN_;
 
  mutable uint64_t called_, saturated_, test_, passTrim_, failTooLarge_, passSC_, failTooNarrow_;

  edm::EDGetTokenT<MeasurementTrackerEvent> mteToken_;
  const ClusterShapeHitFilter   * theFilter;
  const TrackerTopology         * theTopology;
  const TrackerGeometry         * theTracker;
  const MeasurementTrackerEvent * theMTEvent;
  const SiStripNoises           * theNoise;

};

#endif
