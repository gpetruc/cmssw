#ifndef _ClusterShapeDebugTrajectoryFilter_h_
#define _ClusterShapeDebugTrajectoryFilter_h_

#include <vector>
#include <unordered_map>
#include "TrackingTools/TrajectoryFiltering/interface/TrajectoryFilter.h"
#include "RecoTracker/MeasurementDet/interface/MeasurementTrackerEvent.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "RecoTracker/DebugTools/interface/TrackMixingAssociator.h"

class ClusterShapeHitFilter;
class TrackerTopology;
class TrackerGeometry;
class TrackerHitAssociator;
class TrajectoryMeasurement;
class TrajectoryStateOnSurface;
class MeasurementTrackerEvent;
class SiStripNoises;
class TTree;
namespace edm { class Event; class EventSetup; class ConsumesCollector; }


class ClusterShapeDebugTrajectoryFilter : public TrajectoryFilter {
 public:
  //  ClusterShapeDebugTrajectoryFilter(const edm::EventSetup& es);

  ClusterShapeDebugTrajectoryFilter(const edm::ParameterSet &iConfig, edm::ConsumesCollector& iC);

  virtual ~ClusterShapeDebugTrajectoryFilter();

  virtual bool qualityFilter(const TempTrajectory&) const override;
  virtual bool qualityFilter(const Trajectory&) const override;
 
  virtual bool toBeContinued(TempTrajectory&) const override;
  virtual bool toBeContinued(Trajectory&) const override;

  virtual std::string name() const { return "ClusterShapeDebugTrajectoryFilter"; }

  virtual void setEvent(const edm::Event &, const edm::EventSetup &) override;
  void initTree() const;
 protected:
  edm::ParameterSet              theConfig;
  edm::EDGetTokenT<MeasurementTrackerEvent> mteToken_;
  std::vector<edm::EDGetTokenT<edmNew::DetSetVector<SiStripCluster>>> clusterTokens_;
  std::unique_ptr<TrackerHitAssociator>  theAssociator;
  std::unique_ptr<TrackMixingAssociator> theRecoAssociator;
  mutable const ClusterShapeHitFilter   * theFilter;
  mutable const TrackerTopology         * theTopology;
  mutable const TrackerGeometry         * theTracker;
  mutable const MeasurementTrackerEvent * theMTEvent;
  mutable const SiStripNoises           * theNoise;

  mutable TTree * theTree;
  mutable float trackPt_, trackEta_; 
  mutable int   trackAssoc_;
  mutable int   hitDet_, hitLayer_, hitAssoc_, hitSingleSim_, hitUsable_, hitStrips_, hitFirstStrip_, hitCompatPos_, hitCompatNoPos_, hitBadBoundaries_, hitSaturatedStrips_;
  mutable float hitCharge_,  hitChargeNorm_, hitChargeRms_;
  mutable float hitPredPos_, hitPredNoPos_, hitLocalDyDz_, hitLocalDxDz_, hitDetPitch_, hitDetThickness_;
  mutable int   hitPFNearestStrips_, hitPFNearestNSat_; mutable float hitPFNearestCharge_;
  mutable int   hitPFSmallestStrips_, hitPFSmallestNSat_; mutable float hitPFSmallestCharge_;
  mutable int   hitPFLargestStrips_, hitPFLargestNSat_; mutable float hitPFLargestCharge_;
  mutable int   hitTRNearestStrips_, hitTRNearestNSat_; mutable float hitTRNearestCharge_;
  mutable int   hitTRSmallestStrips_, hitTRSmallestNSat_; mutable float hitTRSmallestCharge_;
  mutable int   hitTRLargestStrips_, hitTRLargestNSat_; mutable float hitTRLargestCharge_;
  mutable float hitPFTBestMip_, hitPFTBestSN_;

  typedef std::pair<unsigned int, unsigned int> Id2;
  virtual void fillAssociations(const TrackingRecHit *hit, std::vector<Id2> &out) const  ;
  void fillCluster(const TrackingRecHit *hit, const TrajectoryStateOnSurface &tsos, Id2 simtk, bool mustProject=false) const ;
  bool hasBadStrip(unsigned int detid, int strip) const ;
  mutable std::unordered_map<unsigned int, unsigned int> theStripDetLookup;
  void    indexStripDets() const ;
  virtual void initAssociator(const edm::Event &, const edm::EventSetup &); 

  struct Peak { 
        uint16_t charge, first; uint8_t maxval, nstrips, nsat;  
        Peak() : charge(0), first(0), maxval(0), nstrips(0), nsat(0) {}
        Peak(unsigned int seed, uint8_t val) : charge(val), first(seed), maxval(val), nstrips(1), nsat(val >= 254) {}
        bool operator<(const Peak &p2) const { return maxval > p2.maxval; } 
  };
  void peakFinder(int detid, int ifirst, const std::vector<uint8_t> & ampls, std::vector<Peak> &peaks, float alpha, float seedThr, float clustThr, float stepThr, bool debug=false) const ;
  void subtractPeaks(int detid, int ifirst, const std::vector<uint8_t> & ampls, const std::vector<Peak> &peaksIn, std::vector<uint8_t> & amplsOut, std::vector<Peak> &peaksOut, int maxSize, float stripThr, float seedThr, float clustThr, bool debug=false) const ;
  void thresholdRaiser(int detid, int ifirst, const std::vector<uint8_t> & ampls, std::vector<Peak> &peaks, float stripThr, float seedThr, float clustThr, bool debug=false) const ;
  void dumpCluster(int detid, int ifirst, const std::vector<uint8_t> & ampls, double mip) const ;
  void dumpAssociatedClusters(int detid, const SiStripCluster &cluster, double mip) const ;
  enum FailureCode { Success=0, Trimmable=1, Splittable=2, Overlapping=3, Saturated=4, Helpless=5 };
};

#endif
