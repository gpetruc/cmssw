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
        PeakFinderTest(float mip, uint32_t detid, uint32_t firstStrip, const SiStripNoises *theNoise,
                       float seedCutMIPs, float seedCutSN, 
                       float subclusterCutMIPs, float subclusterCutSN, bool verbose=false) :
            mip_(mip),  detid_(detid), firstStrip_(firstStrip), noiseObj_(theNoise), noises_(theNoise->getRange(detid)),
            subclusterCutMIPs_(subclusterCutMIPs), subclusterCutSN_(subclusterCutSN),
            verbose_(verbose)
        {
            cut_ = std::min<float>(seedCutMIPs*mip, seedCutSN*noiseObj_->getNoise(firstStrip+1, noises_));
        }
   
        bool operator()(uint8_t max, uint8_t min) const {
            return max-min > cut_;
        } 
        bool operator()(const uint8_t *left, const uint8_t *right, const uint8_t *begin, const uint8_t *end) const {
            int yleft  = (left  <  begin ? 0 : *left);
            int yright = (right >= end   ? 0 : *right);
            unsigned int sum = 0, sumall = 0, strips = 0; int maxval = 0; float noise = 0;
            for (const uint8_t *x = left+1; x < right; ++x) {
                int baseline = (yleft * int(right-x) + yright * int(x-left)) / int(right-left);
                sumall += int(*x);
                sum += int(*x) - baseline;
                noise += std::pow(noiseObj_->getNoise(firstStrip_ + int(x-begin), noises_), 2);
                maxval = std::max(maxval, int(*x) - baseline);
                strips++;
            }
            if (verbose_) {
                printf("refined candidate [%4d,%4d] with sum %4d, strips = %d: sum/mip = %5.2f (%+5.2f), sumall/mip = %5.2f (%+5.2f), sum/noise = %6.2f, pede/mip = %5.2f, pede/noise = %6.2f\n", 
                            int((left+1)-begin)+firstStrip_, int((right-1)-begin)+firstStrip_, sum, strips, sum/mip_, log(sum/mip_), sumall/mip_, log(sumall/mip_), sum/std::sqrt(noise), (sumall-sum)/mip_, (sumall-sum)/std::sqrt(noise));
            }
            if (sum/mip_ > subclusterCutMIPs_ && sum/std::sqrt(noise) > subclusterCutSN_) return true;
            //if ((sumall-sum) < 0.3*sum || ( (sumall-sum)/mip_ < 0.2 && (sumall-sum) < 5*std::sqrt(noise) )) {
            //    if (sumall/mip_ > subclusterCutMIPs_ && sumall/std::sqrt(noise) > subclusterCutSN_) return true;
            //}
            return false; 
        }

        private:
            float mip_;
            unsigned int detid_; int firstStrip_;
            const SiStripNoises *noiseObj_;
            SiStripNoises::Range noises_;
            uint8_t cut_;
            float subclusterCutMIPs_, subclusterCutSN_;
            bool verbose_;
    };

}

/*****************************************************************************/
StripSubClusterShapeFilterBase::StripSubClusterShapeFilterBase
   (const edm::ParameterSet &iCfg, edm::ConsumesCollector& iC) :
   label_(iCfg.getUntrackedParameter<std::string>("label","")),
   maxNSat_(iCfg.getParameter<uint32_t>("maxNSat")),
   trimMaxADC_(iCfg.getParameter<double>("trimMaxADC")),
   trimMaxFracTotal_(iCfg.getParameter<double>("trimMaxFracTotal")),
   trimMaxFracNeigh_(iCfg.getParameter<double>("trimMaxFracNeigh")),
   maxTrimmedSizeDiffPos_(iCfg.getParameter<double>("maxTrimmedSizeDiffPos")),
   maxTrimmedSizeDiffNeg_(iCfg.getParameter<double>("maxTrimmedSizeDiffNeg")),
   subclusterWindow_(iCfg.getParameter<double>("subclusterWindow")),
   seedCutMIPs_(iCfg.getParameter<double>("seedCutMIPs")),
   seedCutSN_(iCfg.getParameter<double>("seedCutSN")),
   subclusterCutMIPs_(iCfg.getParameter<double>("subclusterCutMIPs")),
   subclusterCutSN_(iCfg.getParameter<double>("subclusterCutSN")),
   called_(0),saturated_(0),test_(0),passTrim_(0),failTooLarge_(0),passSC_(0),failTooNarrow_(0),
   theFilter(0),
   theTracker(0)
{
}

/*****************************************************************************/
StripSubClusterShapeFilterBase::~StripSubClusterShapeFilterBase()
{
    std::cout << "StripSubClusterShapeFilterBase " << label_ <<": called        " << called_ << std::endl;
    std::cout << "StripSubClusterShapeFilterBase " << label_ <<": saturated     " << saturated_ << std::endl;
    std::cout << "StripSubClusterShapeFilterBase " << label_ <<": test          " << test_ << std::endl;
    std::cout << "StripSubClusterShapeFilterBase " << label_ <<": failTooNarrow " << failTooNarrow_ << std::endl;
    std::cout << "StripSubClusterShapeFilterBase " << label_ <<": passTrim      " << passTrim_ << std::endl;
    std::cout << "StripSubClusterShapeFilterBase " << label_ <<": passSC        " << passSC_ << std::endl;
    std::cout << "StripSubClusterShapeFilterBase " << label_ <<": failTooLarge  " << failTooLarge_ << std::endl;
}


bool StripSubClusterShapeFilterBase::testLastHit
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
      const std::vector<uint8_t> &ampls = cluster.amplitudes();
      
      // pass-through of trivial case 
      if (std::abs(hitPredPos) < 1.5f && hitStrips <= 2) {
        return true;
      }


      // compute number of consecutive saturated strips
      unsigned int thisSat = (ampls[0] >= 254), maxSat = thisSat;
      for (unsigned int i = 1, n = ampls.size(); i < n; ++i) {
          if (ampls[i] >= 254) {
              thisSat++;
          } else if (thisSat > 0) {
              maxSat = std::max<int>(maxSat, thisSat);
              thisSat = 0;
          }
      }
      if (thisSat > 0) {
          maxSat = std::max<int>(maxSat, thisSat);
      }
      if (maxSat >= maxNSat_) {
          saturated_++;
          return true;
      }
     
      // trimming
      test_++;
      unsigned int hitStripsTrim = ampls.size();
      int sum = std::accumulate(ampls.begin(), ampls.end(), 0);
      uint8_t trimCut = std::min<uint8_t>(trimMaxADC_, std::floor(trimMaxFracTotal_ * sum));
      std::vector<uint8_t>::const_iterator begin = ampls.begin();
      std::vector<uint8_t>::const_iterator last = ampls.end()-1;
      while (hitStripsTrim > 1 && (*begin < std::max<uint8_t>(trimCut, trimMaxFracNeigh_*(*(begin+1)))) ) { hitStripsTrim--; ++begin; }
      while (hitStripsTrim > 1 && (*last  < std::max<uint8_t>(trimCut, trimMaxFracNeigh_*(*(last -1)))) ) { hitStripsTrim--; --last; }

      if (hitStripsTrim < std::floor(std::abs(hitPredPos)-maxTrimmedSizeDiffNeg_)) {
          failTooNarrow_++;
          return false;
      } else if (hitStripsTrim <= std::ceil(std::abs(hitPredPos)+maxTrimmedSizeDiffPos_)) {
          passTrim_++;
          return true;
      } 

      const StripGeomDetUnit* stripDetUnit = dynamic_cast<const StripGeomDetUnit *>(det);
      if (det == 0) throw cms::Exception("Strip not a StripGeomDetUnit?") << " on " << detId.rawId() << "\n";

      float MeVperADCStrip =  3.61e-06*265; 
      float mip = 3.9 / ( MeVperADCStrip/stripDetUnit->surface().bounds().thickness() ); // 3.9 MeV/cm = ionization in silicon 
      float mipnorm = mip/std::abs(ldir.z());
      utils::SlidingPeakFinder pf(std::max<int>(2,std::ceil(std::abs(hitPredPos)+subclusterWindow_)));
      PeakFinderTest test(mipnorm, detId(), cluster.firstStrip(), theNoise, seedCutMIPs_, seedCutSN_, subclusterCutMIPs_, subclusterCutSN_);
      if (pf.apply(cluster.amplitudes(), test)) {
          passSC_++;
          return true;
      } else {
          failTooLarge_++;
          return false;
      }

    }
    return true; 
}

void StripSubClusterShapeFilterBase::setEventBase
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

  edm::ESHandle<SiStripNoises>  noise;
  es.get<SiStripNoisesRcd>().get(noise);
  theNoise = noise.product();
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

/*****************************************************************************/
bool StripSubClusterShapeSeedFilter::compatible
  (const  TrajectoryStateOnSurface & tsos, SeedingHitSet::ConstRecHitPointer thit) const
{
   const TrackingRecHit* hit = thit->hit();
   if (hit == 0 || !hit->isValid()) return true;
   if (hit->geographicalId().subdetId() <= 2) return true; // we look only at strips for now
   return testLastHit(hit, tsos, false);
}

