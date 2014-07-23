#include "RecoTracker/DebugTools/interface/StripSubClusterShapeTrajectoryFilter.h"

#include <map>
#include <set>
#include <algorithm>

#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/TrackerRecHit2D/interface/ProjectedSiStripRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "RecoPixelVertexing/PixelLowPtUtilities/interface/ClusterShapeHitFilter.h"
#include "RecoTracker/Record/interface/CkfComponentsRecord.h"
#include "RecoTracker/MeasurementDet/src/TkMeasurementDetSet.h"
#include "RecoTracker/MeasurementDet/interface/MeasurementTrackerEvent.h"
#include "TrackingTools/PatternTools/interface/TempTrajectory.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "RecoTracker/DebugTools/interface/SlidingPeakFinder.h"

#include "CondFormats/SiStripObjects/interface/SiStripNoises.h"
#include "CondFormats/DataRecord/interface/SiStripNoisesRcd.h"

namespace {
    struct PeakFinderTest {
        PeakFinderTest(uint8_t cut, uint32_t detid, uint32_t firstStrip, float mip, const StripSubClusterShapeTrajectoryFilter &main) :
            cut_(cut), detid_(detid), firstStrip_(firstStrip), mip_(mip), main_(main) {}
   
        bool operator()(uint8_t max, uint8_t min) const {
            return max-min > cut_;
            //printf("cluster candidate with max = %3d, maxmin = %3d.  (max-min)/mip = %.2f, (max-min)/noise = %.2f\n", int(max), int(min), int(max-min)/mip_, int(max-min)/noise_);
        } 
        bool operator()(const uint8_t *left, const uint8_t *right, const uint8_t *begin, const uint8_t *end) const {
            return main_.testSubCluster(detid_, firstStrip_, left, right, begin, end, mip_);
        }
        private:
            uint8_t cut_;
            unsigned int detid_; int firstStrip_;
            float mip_;
            const StripSubClusterShapeTrajectoryFilter & main_;
    };
}

/*****************************************************************************/
StripSubClusterShapeTrajectoryFilter::StripSubClusterShapeTrajectoryFilter
   (const edm::ParameterSet &iCfg, edm::ConsumesCollector& iC) :
   maxNSat_(iCfg.getParameter<uint32_t>("maxNSat")),
   maxPredSize_(iCfg.getParameter<double>("maxPredSize")),
   sizeCut_(iCfg.getParameter<double>("sizeCut")),
   called_(0),test_(0),pass_(0),fullpass_(0),
   mteToken_(iC.consumes<MeasurementTrackerEvent>(iCfg.getParameter<edm::InputTag>("MeasurementTrackerEvent"))),
   theFilter(0),
   theTopology(0),
   theTracker(0)
{
}

/*****************************************************************************/
StripSubClusterShapeTrajectoryFilter::~StripSubClusterShapeTrajectoryFilter()
{
    std::cout << "StripSubClusterShapeTrajectoryFilter:  called " << called_ << ", tested " << test_ << ", passing " << pass_ << " + " << fullpass_ << " = " << (pass_+fullpass_) << std::endl;
}


/*****************************************************************************/
bool StripSubClusterShapeTrajectoryFilter::toBeContinued
(Trajectory& trajectory) const 
{
   throw cms::Exception("toBeContinued(Traj) instead of toBeContinued(TempTraj)");
}

/*****************************************************************************/
bool StripSubClusterShapeTrajectoryFilter::toBeContinued
(TempTrajectory& trajectory) const 
{
   const TempTrajectory::DataContainer & tms = trajectory.measurements();

   const TrajectoryMeasurement &last = *tms.rbegin();
   const TrackingRecHit* hit = last.recHit()->hit();
   if (!last.updatedState().isValid()) return true;
   if (hit == 0 || !hit->isValid()) return true;
   if (hit->geographicalId().subdetId() <= 2) return true; // we look only at strips for now
   return testLastHit(hit, last.updatedState(), false);
}

 
bool StripSubClusterShapeTrajectoryFilter::testLastHit
   (const TrackingRecHit *hit, const TrajectoryStateOnSurface &tsos, bool mustProject) const
{
   const TrackerSingleRecHit *stripHit = 0;
   if (typeid(*hit) == typeid(SiStripMatchedRecHit2D)) {
      const SiStripMatchedRecHit2D & mhit = static_cast<const SiStripMatchedRecHit2D &>(*hit);
      SiStripRecHit2D mono = mhit.monoHit();
      SiStripRecHit2D stereo = mhit.stereoHit();
      return testLastHit(&mono, tsos, true) && testLastHit(&stereo, tsos, true);
   } else if (typeid(*hit) == typeid(ProjectedSiStripRecHit2D)) {
      const ProjectedSiStripRecHit2D & mhit = static_cast<const ProjectedSiStripRecHit2D &>(*hit);
      const SiStripRecHit2D &orig = mhit.originalHit();
      return testLastHit(&orig,   tsos, true);
   } else if ((stripHit = dynamic_cast<const TrackerSingleRecHit *>(hit)) != 0) {
      DetId detId = hit->geographicalId();

      GlobalVector gdir = tsos.globalDirection();
      GlobalPoint  gpos = tsos.globalPosition();
      const GeomDet *det = theTracker->idToDet(detId);
      LocalVector ldir = det->toLocal(gdir);
      LocalPoint  lpos;
      if (mustProject) {
         lpos = det->toLocal(gpos);
         lpos -= ldir * lpos.z()/ldir.z();
      } else {
         lpos = tsos.localPosition();
      }
      int hitStrips; float hitPredPos;
      const SiStripCluster &cluster = stripHit->stripCluster();
      bool usable = theFilter->getSizes(detId, cluster, lpos, ldir, hitStrips, hitPredPos);
      if (!usable) return true;

      called_++;
      if (fabs(hitPredPos) > maxPredSize_) return true;
      uint32_t nsat = std::count_if(cluster.amplitudes().begin(), cluster.amplitudes().end(), [](const uint8_t v) { return v >= 254; });
      if (nsat > maxNSat_) return true;

      test_++; 
      if (hitStrips - fabs(hitPredPos) < sizeCut_) {
          pass_++;
          return true;
      } else {
          const StripGeomDetUnit* stripDetUnit = dynamic_cast<const StripGeomDetUnit *>(det);
          if (det == 0) throw cms::Exception("Strip not a StripGeomDetUnit?") << " on " << detId.rawId() << " aka " << theTopology->print(detId) << "\n";

          float MeVperADCStrip =  3.61e-06*265; 
          float mip = 3.9 / ( MeVperADCStrip/stripDetUnit->surface().bounds().thickness() ); // 3.9 MeV/cm = ionization in silicon 
          uint8_t cut = std::min(0.4*mip, 20.);
          PeakFinderTest test(cut, detId(), cluster.firstStrip(), mip, *this);
          utils::SlidingPeakFinder pf(ceil(fabs(hitPredPos)+0.7));
          if (pf.apply(cluster.amplitudes(), test)) {
              fullpass_++;
              return true;
          }
      }
        
      return false;
    }
    return true; 
}

bool StripSubClusterShapeTrajectoryFilter::testSubCluster(uint32_t detid, int firstStrip, const uint8_t *left, const uint8_t *right, const uint8_t *begin, const uint8_t *end, float mip) const {
    SiStripNoises::Range noises = theNoise->getRange(detid);
    int yleft  = (left  <  begin ? 0 : *left);
    int yright = (right >= end   ? 0 : *right);
    unsigned int sum = 0, strips = 0; int maxval = 0; float noise = 0;
    for (const uint8_t *x = left+1; x < right; ++x) {
        int baseline = (yleft * int(right-x) + yright * int(x-left)) / int(right-left);
        sum += int(*x) - baseline;
        noise += std::pow(theNoise->getNoise(firstStrip + int(x-begin), noises), 2);
        maxval = std::max(maxval, int(*x) - baseline);
        strips++;
    }
    printf("refined candidate with sum %d, strips = %d: sum/mip = %.2f, sum/noise = %.2f\n", sum, strips, sum/mip, sum/std::sqrt(noise));
    return sum/mip > 0.6 && sum/std::sqrt(noise) > 9;
}
/*****************************************************************************/
bool StripSubClusterShapeTrajectoryFilter::qualityFilter
  (const Trajectory& trajectory) const
{
   const Trajectory::DataContainer & tms = trajectory.measurements();

   const TrajectoryMeasurement &last = *tms.rbegin();
   const TrackingRecHit* hit = last.recHit()->hit();
   if (!last.updatedState().isValid()) return true;
   if (hit == 0 || !hit->isValid()) return true;
   if (hit->geographicalId().subdetId() <= 2) return true; // we look only at strips for now
   return testLastHit(hit, last.updatedState(), false);

}

/*****************************************************************************/
bool StripSubClusterShapeTrajectoryFilter::qualityFilter
  (const TempTrajectory& trajectory) const
{
   const TempTrajectory::DataContainer & tms = trajectory.measurements();
   const TrajectoryMeasurement &last = *tms.rbegin();
   const TrackingRecHit* hit = last.recHit()->hit();
   if (!last.updatedState().isValid()) return true;
   if (hit == 0 || !hit->isValid()) return true;
   if (hit->geographicalId().subdetId() <= 2) return true; // we look only at strips for now
   return testLastHit(hit, last.updatedState(), false);

}

void StripSubClusterShapeTrajectoryFilter::setEvent
  (const edm::Event &event, const edm::EventSetup &es) 
{
  // Get tracker geometry
  edm::ESHandle<TrackerGeometry> tracker;
  es.get<TrackerDigiGeometryRecord>().get(tracker);
  theTracker = tracker.product();

  edm::ESHandle<ClusterShapeHitFilter> shape;
  es.get<CkfComponentsRecord>().get("ClusterShapeHitFilter",shape);
  theFilter = shape.product();

  //Retrieve tracker topology from geometry
  edm::ESHandle<TrackerTopology> tTopo;
  es.get<IdealGeometryRecord>().get(tTopo);
  theTopology = tTopo.product();

  edm::ESHandle<SiStripNoises>  noise;
  es.get<SiStripNoisesRcd>().get(noise);
  theNoise = noise.product();
  
  edm::Handle<MeasurementTrackerEvent> mte;
  event.getByToken(mteToken_, mte);
  theMTEvent = mte.product();
}

