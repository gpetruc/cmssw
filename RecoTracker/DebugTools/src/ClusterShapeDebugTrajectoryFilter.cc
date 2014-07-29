#include "RecoTracker/DebugTools/interface/ClusterShapeDebugTrajectoryFilter.h"

#include <map>
#include <set>
#include <algorithm>

#include "CommonTools/UtilAlgos/interface/TFileService.h"
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
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "RecoTracker/DebugTools/interface/SlidingPeakFinder.h"

#include "CondFormats/SiStripObjects/interface/SiStripNoises.h"
#include "CondFormats/DataRecord/interface/SiStripNoisesRcd.h"

#include <TTree.h>

namespace {
    struct PeakFinderTest {
        PeakFinderTest(float mip, uint32_t detid, uint32_t firstStrip, const SiStripNoises *theNoise,
                       float seedCutMIPs, float seedCutSN, 
                       float subclusterCutMIPs, float subclusterCutSN, bool verbose=false, bool docut=true) :
            mip_(mip),  detid_(detid), firstStrip_(firstStrip), noiseObj_(theNoise), noises_(theNoise->getRange(detid)),
            subclusterCutMIPs_(subclusterCutMIPs), subclusterCutSN_(subclusterCutSN),
            verbose_(verbose), docut_(docut)
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
            if (!docut_) return false;
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
            bool verbose_, docut_;
    };


}

/*****************************************************************************/
ClusterShapeDebugTrajectoryFilter::ClusterShapeDebugTrajectoryFilter
   (const edm::ParameterSet &iCfg, edm::ConsumesCollector& iC) :
   theConfig(iCfg.getParameter<edm::ParameterSet>("HitAssociatorPSet")),
   mteToken_(iC.consumes<MeasurementTrackerEvent>(iCfg.getParameter<edm::InputTag>("MeasurementTrackerEvent"))),
   theAssociator(),
   theRecoAssociator(),
   theFilter(0),
   theTopology(0),
   theTracker(0),
   maxNSat_(iCfg.getParameter<uint32_t>("maxNSat")),
   trimMaxADC_(iCfg.getParameter<double>("trimMaxADC")),
   trimMaxADCTight_(iCfg.getParameter<double>("trimMaxADCTight")),
   trimMaxFracTotal_(iCfg.getParameter<double>("trimMaxFracTotal")),
   trimMaxFracNeigh_(iCfg.getParameter<double>("trimMaxFracNeigh")),
   trimMaxFracTotalTight_(iCfg.getParameter<double>("trimMaxFracTotalTight")),
   trimMaxFracNeighTight_(iCfg.getParameter<double>("trimMaxFracNeighTight")),
   maxTrimmedSizeDiffPos_(iCfg.getParameter<double>("maxTrimmedSizeDiffPos")),
   maxTrimmedSizeDiffNeg_(iCfg.getParameter<double>("maxTrimmedSizeDiffNeg")),
   maxTrimmedSizeDiffPosTight_(iCfg.getParameter<double>("maxTrimmedSizeDiffPosTight")),
   maxTrimmedSizeDiffNegTight_(iCfg.getParameter<double>("maxTrimmedSizeDiffNegTight")),
   subclusterWindow_(iCfg.getParameter<double>("subclusterWindow")),
   seedCutMIPs_(iCfg.getParameter<double>("seedCutMIPs")),
   seedCutSN_(iCfg.getParameter<double>("seedCutSN")),
   subclusterCutMIPs_(iCfg.getParameter<double>("subclusterCutMIPs")),
   subclusterCutSN_(iCfg.getParameter<double>("subclusterCutSN")),
   subclusterCutMIPsTight_(iCfg.getParameter<double>("subclusterCutMIPsTight")),
   subclusterCutSNTight_(iCfg.getParameter<double>("subclusterCutSNTight")),
   maxBadHits_(iCfg.getParameter<uint32_t>("maxBadHits")),
   theTree(0)
{
    if (iCfg.existsAs<std::vector<edm::InputTag>>("referenceStripClusters")) {
        std::vector<edm::InputTag> tags = iCfg.getParameter<std::vector<edm::InputTag>>("referenceStripClusters");
        for (const auto &t : tags) {
            clusterTokens_.push_back(iC.consumes<edmNew::DetSetVector<SiStripCluster>>(t));
        }
    }
    if (iCfg.exists("layerMask")) {
        const edm::ParameterSet &iLM = iCfg.getParameter<edm::ParameterSet>("layerMask");
        const char *ndets[4] =  { "TIB", "TID", "TOB", "TEC" };
        const int   idets[4] =  {   3,     4,     5,     6   };
        for (unsigned int i = 0; i < 4; ++i) {
            if (iLM.existsAs<bool>(ndets[i])) {
                std::fill(layerMask_[idets[i]].begin(), layerMask_[idets[i]].end(), iLM.getParameter<bool>(ndets[i]));
            } else {
                layerMask_[idets[i]][0] = 2;
                std::fill(layerMask_[idets[i]].begin()+1, layerMask_[idets[i]].end(), 0);
                for (uint32_t lay : iLM.getParameter<std::vector<uint32_t>>(ndets[i])) {
                    layerMask_[idets[i]][lay] = 1;
                }
            }
        }
    } else {
        for (auto & arr : layerMask_) {
            std::fill(arr.begin(), arr.end(), 1);
        }
    }
    printf("Layer mask: \n");
    for (unsigned int i = 3; i <= 6; ++i) {
        printf("Subdetector %d mode = %d: layers: ", i, layerMask_[i][0]);
        for (unsigned int j = 1; j < layerMask_[i].size(); ++j) {
            printf("[%d] = %d   ", j, layerMask_[i][j]);
        }
        printf("\n");
    }
    printf("\n");

}

/*****************************************************************************/
ClusterShapeDebugTrajectoryFilter::~ClusterShapeDebugTrajectoryFilter()
{
}

void
ClusterShapeDebugTrajectoryFilter::initTree() const 
{
   edm::Service<TFileService> fs;
   theTree = fs->make<TTree>("t","t");
   theTree->Branch("trackPt", &trackPt_, "trackPt/F");
   theTree->Branch("trackEta", &trackEta_, "trackEta/F");
   theTree->Branch("trackAssoc", &trackAssoc_, "trackAssoc/I");
   theTree->Branch("hitDet", &hitDet_, "hitDet/I");
   theTree->Branch("hitLayer", &hitLayer_, "hitLayer/I");
   theTree->Branch("hitAssoc", &hitAssoc_, "hitAssoc/I");
   theTree->Branch("hitSingleSim", &hitSingleSim_, "hitSingleSim/I");
   theTree->Branch("hitUsable", &hitUsable_, "hitUsable/I");
   theTree->Branch("hitCompatPos", &hitCompatPos_, "hitCompatPos/I");
   theTree->Branch("hitCompatNoPos", &hitCompatNoPos_, "hitCompatNoPos/I");
   theTree->Branch("hitBadBoundaries", &hitBadBoundaries_, "hitBadBoundaries/I");
   theTree->Branch("hitStrips", &hitStrips_, "hitStrips/I");
   theTree->Branch("hitStripsTrim", &hitStripsTrim_, "hitStripsTrim/I");
   theTree->Branch("hitFirstStrip", &hitFirstStrip_, "hitFirstStrip/I");
   theTree->Branch("hitCharge", &hitCharge_, "hitCharge/F");
   theTree->Branch("hitChargeNorm", &hitChargeNorm_, "hitChargeNorm/F");
   theTree->Branch("hitChargeRms", &hitChargeRms_, "hitChargeRms/F");
   theTree->Branch("hitPredPos", &hitPredPos_, "hitPredPos/F");
   theTree->Branch("hitPredNoPos", &hitPredNoPos_, "hitPredNoPos/F");
   theTree->Branch("hitLocalDxDz", &hitLocalDxDz_, "hitLocalDxDz/F");
   theTree->Branch("hitLocalDyDz", &hitLocalDyDz_, "hitLocalDyDz/F");
   theTree->Branch("hitDetPitch", &hitDetPitch_, "hitDetPitch/F");
   theTree->Branch("hitDetThickness", &hitDetThickness_, "hitDetThickness/F");
   theTree->Branch("hitSaturatedStrips", &hitSaturatedStrips_, "hitSaturatedStrips/I");
   theTree->Branch("hitPFNearestStrips", &hitPFNearestStrips_, "hitPFNearestStrips/I");
   theTree->Branch("hitPFNearestNSat",   &hitPFNearestNSat_,   "hitPFNearestNSat/I");
   theTree->Branch("hitPFNearestCharge", &hitPFNearestCharge_, "hitPFNearestCharge/F");
   theTree->Branch("hitPFSmallestStrips", &hitPFSmallestStrips_, "hitPFSmallestStrips/I");
   theTree->Branch("hitPFSmallestNSat",   &hitPFSmallestNSat_,   "hitPFSmallestNSat/I");
   theTree->Branch("hitPFSmallestCharge", &hitPFSmallestCharge_, "hitPFSmallestCharge/F");
   theTree->Branch("hitPFLargestStrips", &hitPFLargestStrips_, "hitPFLargestStrips/I");
   theTree->Branch("hitPFLargestNSat",   &hitPFLargestNSat_,   "hitPFLargestNSat/I");
   theTree->Branch("hitPFLargestCharge", &hitPFLargestCharge_, "hitPFLargestCharge/F");
   theTree->Branch("hitTRNearestStrips", &hitTRNearestStrips_, "hitTRNearestStrips/I");
   theTree->Branch("hitTRNearestNSat",   &hitTRNearestNSat_,   "hitTRNearestNSat/I");
   theTree->Branch("hitTRNearestCharge", &hitTRNearestCharge_, "hitTRNearestCharge/F");
   theTree->Branch("hitTRSmallestStrips", &hitTRSmallestStrips_, "hitTRSmallestStrips/I");
   theTree->Branch("hitTRSmallestNSat",   &hitTRSmallestNSat_,   "hitTRSmallestNSat/I");
   theTree->Branch("hitTRSmallestCharge", &hitTRSmallestCharge_, "hitTRSmallestCharge/F");
   theTree->Branch("hitTRLargestStrips", &hitTRLargestStrips_, "hitTRLargestStrips/I");
   theTree->Branch("hitTRLargestNSat",   &hitTRLargestNSat_,   "hitTRLargestNSat/I");
   theTree->Branch("hitTRLargestCharge", &hitTRLargestCharge_, "hitTRLargestCharge/F");
   theTree->Branch("hitPFTBestMip", &hitPFTBestMip_, "hitPFTBestMip/F");
   theTree->Branch("hitPFTBestSN", &hitPFTBestSN_, "hitPFTBestSN/F");
}

/*****************************************************************************/
bool ClusterShapeDebugTrajectoryFilter::toBeContinued
(Trajectory& trajectory) const 
{
   throw cms::Exception("toBeContinued(Traj) instead of toBeContinued(TempTraj)");
}

/*****************************************************************************/
bool ClusterShapeDebugTrajectoryFilter::toBeContinued
(TempTrajectory& trajectory) const 
{
   if (theTree == 0) initTree();
   TempTrajectory::DataContainer tms = trajectory.measurements();

   std::vector<Id2> work;
   std::map<Id2,std::set<Id2> > scoreboard;
   std::set<Id2> countlayers;
   for(TempTrajectory::DataContainer::const_iterator
         tm = tms.rbegin(); tm!= tms.rend(); --tm)
   {
      if (tm == tms.rbegin()) continue; // skip last measurement
      const TrackingRecHit* hit = tm->recHit()->hit();

      if(hit == 0 || !hit->isValid() || hit->geographicalId().rawId() == 0) continue;

      Id2 mylayer = Id2(hit->geographicalId().subdetId(), theTopology->layer(hit->geographicalId()));
      countlayers.insert(mylayer);

      fillAssociations(hit, work);
      //std::cout << "  on layer " << mylayer.first << "/" << mylayer.second << " found match with " << work.size() << " simtracks." << std::endl;
      for (Id2 simtk: work) {
         //std::cout << "       simtk  " << simtk.first << "/" << simtk.second << std::endl;
         scoreboard[simtk].insert(mylayer);
      }
   }

   Id2 bestSim; unsigned int layers = 0;
   for (const auto &p : scoreboard) {
      if (p.second.size() > layers) {
         layers = p.second.size();
         bestSim = p.first;
      }
   }
   EncodedEventId evid(bestSim.first);
   //if (layers > 0) {
   //}

   if (layers >= 3 && layers >= 0.67*countlayers.size()) {
      trackAssoc_ = (evid.bunchCrossing() == 0 ? 1 : -1); 
      //std::cout << "best match: " << layers << " layers (out of " << countlayers.size() << "), id " << bestSim.first << ", tk " << bestSim.second << ", event/bx = " << evid.event() << "/" << evid.bunchCrossing() << std::endl;
   } else {
      trackAssoc_ = 0;
   }

   const TrajectoryMeasurement &last = *tms.rbegin();
   const TrackingRecHit* hit = last.recHit()->hit();
   if (!last.updatedState().isValid()) return true;
   if (hit == 0 || !hit->isValid()) return true;
   DetId detId = hit->geographicalId();
   if (detId.subdetId() <= 2) return true; // we look only at strips for now

   if (layerMask_[detId.subdetId()][0] == 0) {
       std::cout << "Not filtering on " << detId.subdetId() << std::endl;
       return true; // no filtering here
   } else if (layerMask_[detId.subdetId()][0] == 2) {
       std::cout << "Selective filtering on " << detId.subdetId() << std::endl;
       unsigned int ilayer = theTopology->layer(detId);          
       if (layerMask_[detId.subdetId()][ilayer] == 0) {
           std::cout << "Not filtering on " << detId.subdetId() << "/" << ilayer << std::endl;
           return true; // no filtering here
       }
   }

   bool theOk = fillCluster(hit, last.updatedState(), bestSim);
   return theOk ? true : true; // stupid
}

void ClusterShapeDebugTrajectoryFilter::fillAssociations
   (const TrackingRecHit *hit, std::vector<Id2> &out) const
{
   assert(theAssociator);
   out.clear();
   std::vector<PSimHit> simHits = theAssociator->associateHit(*hit);
   hitSingleSim_ = (simHits.size() == 1);
   for (const PSimHit &hit : simHits) {
      Id2 id2(hit.eventId().rawId(), hit.trackId());
      out.push_back(id2);
   }
}
 
bool ClusterShapeDebugTrajectoryFilter::fillCluster
   (const TrackingRecHit *hit, const TrajectoryStateOnSurface &tsos, Id2 simtk, bool mustProject) const
{
   const TrackerSingleRecHit *stripHit = 0;
   if (typeid(*hit) == typeid(SiStripMatchedRecHit2D)) {
      const SiStripMatchedRecHit2D & mhit = static_cast<const SiStripMatchedRecHit2D &>(*hit);
      SiStripRecHit2D mono = mhit.monoHit();
      SiStripRecHit2D stereo = mhit.stereoHit();
      return  fillCluster(&mono,   tsos, simtk, true) && 
              fillCluster(&stereo, tsos, simtk, true);
   } else if (typeid(*hit) == typeid(ProjectedSiStripRecHit2D)) {
      const ProjectedSiStripRecHit2D & mhit = static_cast<const ProjectedSiStripRecHit2D &>(*hit);
      const SiStripRecHit2D &orig = mhit.originalHit();
      return fillCluster(&orig,   tsos, simtk, true);
   } else if ((stripHit = dynamic_cast<const TrackerSingleRecHit *>(hit)) != 0) {
      DetId detId = hit->geographicalId();
      Id2 mylayer = Id2(detId.subdetId(), theTopology->layer(detId));
      hitDet_   = mylayer.first;
      hitLayer_ = mylayer.second;
      std::vector<Id2> work;
      fillAssociations(hit, work);
      if (work.empty()) {
         hitAssoc_ = 0;
      } else if (simtk.second != 0 && std::find(work.begin(), work.end(), simtk) != work.end()) {
         hitAssoc_ = 1;
      } else {
         hitAssoc_ = -1;
      }
     
      trackPt_  = tsos.globalMomentum().perp();
      trackEta_ = tsos.globalMomentum().eta();
 
      GlobalVector gdir = tsos.globalDirection();
      GlobalPoint  gpos = tsos.globalPosition();
      const GeomDet *det = theTracker->idToDet(detId);
      LocalVector ldir = det->toLocal(gdir);
      LocalPoint  lpos, lp0(0.,0.,0.);
      if (mustProject) {
         lpos = det->toLocal(gpos);
         // now here we do the transformation
         lpos -= ldir * lpos.z()/ldir.z();
      } else {
         lpos = tsos.localPosition();
      }
      const SiStripCluster &cluster = stripHit->stripCluster();
      hitUsable_ = theFilter->getSizes(detId, cluster, lpos, ldir, hitStrips_, hitPredPos_);
      hitUsable_ = theFilter->getSizes(detId, cluster, lp0,  ldir, hitStrips_, hitPredNoPos_);
      hitCompatPos_   = theFilter->isCompatible(detId, cluster, lpos, ldir);
      hitCompatNoPos_ = theFilter->isCompatible(detId, cluster, lp0,  ldir);
      hitLocalDxDz_ = ldir.y()/ldir.z();
      hitLocalDyDz_ = ldir.x()/ldir.z();
      hitFirstStrip_ = cluster.firstStrip();
      const StripGeomDetUnit * stripDet = dynamic_cast<const StripGeomDetUnit *>(det);
      if (stripDet) {
          hitDetPitch_ = stripDet->specificTopology().localPitch(lpos);
          hitDetThickness_ = stripDet->surface().bounds().thickness();
      } else {
        hitDetPitch_ = -1;
        hitDetThickness_ = -1;
      }
      int sfirst = cluster.firstStrip(), slast = sfirst + cluster.amplitudes().size()-1;
      hitBadBoundaries_ = hasBadStrip(detId(), sfirst-1) || hasBadStrip(detId(), slast+1);

      const std::vector<uint8_t> &ampls = cluster.amplitudes();
      double sum = 0, sumx = 0, sumx2 = 0; 
      for (unsigned int i = 0, n = ampls.size(); i < n; ++i) {
         sum   +=         ampls[i];
         sumx  +=     i * ampls[i];
         sumx2 += i * i * ampls[i];
      }
      sumx /= sum; sumx2 /= sum;
      hitChargeRms_ = sumx2 - sumx*sumx;
    
      uint8_t trimCut = std::min<uint8_t>(trimMaxADC_, std::floor(trimMaxFracTotal_ * sum));
      hitStripsTrim_ = hitStrips_;
      std::vector<uint8_t>::const_iterator begin = cluster.amplitudes().begin();
      std::vector<uint8_t>::const_iterator last = cluster.amplitudes().end()-1;
      while (hitStripsTrim_ > 1 && (*begin < std::max<uint8_t>(trimCut, trimMaxFracNeigh_*(*(begin+1)))) ) { hitStripsTrim_--; ++begin; }
      while (hitStripsTrim_ > 1 && (*last  < std::max<uint8_t>(trimCut, trimMaxFracNeigh_*(*(last -1)))) ) { hitStripsTrim_--; --last; }
    
      uint8_t trimCutTight = std::min<uint8_t>(trimMaxADCTight_, std::floor(trimMaxFracTotalTight_ * sum));
      unsigned int hitStripsTrimTight = hitStrips_;
      begin = cluster.amplitudes().begin();
      last = cluster.amplitudes().end()-1;
      while (hitStripsTrimTight > 1 && (*begin < std::max<uint8_t>(trimCutTight, trimMaxFracNeighTight_*(*(begin+1)))) ) { hitStripsTrimTight--; ++begin; }
      while (hitStripsTrimTight > 1 && (*last  < std::max<uint8_t>(trimCutTight, trimMaxFracNeighTight_*(*(last -1)))) ) { hitStripsTrimTight--; --last; }



      const StripGeomDetUnit* stripDetUnit = dynamic_cast<const StripGeomDetUnit *>(det);
      if (det == 0) throw cms::Exception("Strip not a StripGeomDetUnit?") << " on " << detId.rawId() << " aka " << theTopology->print(detId) << "\n";

      float MeVperADCStrip =  3.61e-06*265; 
      float norm = MeVperADCStrip/stripDetUnit->surface().bounds().thickness();
      float mip = 3.9 / norm; // 3.9 MeV/cm = ionization in silicon 
      hitCharge_ = sum * norm;
      hitChargeNorm_ = hitCharge_ / std::abs(ldir.z());
      std::vector<Peak> peaksPF, peaksTR;
      peakFinder(detId(),sfirst,ampls,peaksPF,0.4, 5., 9., 4., false);
      thresholdRaiser(detId(),sfirst,ampls,peaksTR,4., 5., 9., false);
      if (!peaksPF.empty()) {
        std::sort(peaksPF.begin(),peaksPF.end(), [](const Peak &a, const Peak &b) { return a.nstrips < b.nstrips; });
        hitPFSmallestStrips_ = peaksPF.front().nstrips;
        hitPFSmallestNSat_   = peaksPF.front().nsat;
        hitPFSmallestCharge_ = peaksPF.front().charge * norm; 
        hitPFLargestStrips_ = peaksPF.back().nstrips;
        hitPFLargestNSat_   = peaksPF.back().nsat;
        hitPFLargestCharge_ = peaksPF.back().charge * norm; 
        std::sort(peaksPF.begin(),peaksPF.end(), [=](const Peak &a, const Peak &b) { return std::abs(a.nstrips-std::abs(hitPredPos_)) <  std::abs(b.nstrips-std::abs(hitPredPos_)); });  
        hitPFNearestStrips_ = peaksPF.front().nstrips;
        hitPFNearestNSat_   = peaksPF.front().nsat;
        hitPFNearestCharge_ = peaksPF.front().charge * norm; 
      } else {
        hitPFNearestStrips_ = hitPFSmallestStrips_ = hitPFLargestStrips_ = -1;
      }
      if (!peaksTR.empty()) {
        std::sort(peaksTR.begin(),peaksTR.end(), [](const Peak &a, const Peak &b) { return a.nstrips < b.nstrips; });
        hitTRSmallestStrips_ = peaksTR.front().nstrips;
        hitTRSmallestNSat_   = peaksTR.front().nsat;
        hitTRSmallestCharge_ = peaksTR.front().charge * norm; 
        hitTRLargestStrips_ = peaksTR.back().nstrips;
        hitTRLargestNSat_   = peaksTR.back().nsat;
        hitTRLargestCharge_ = peaksTR.back().charge * norm; 
        std::sort(peaksTR.begin(),peaksTR.end(), [=](const Peak &a, const Peak &b) { return std::abs(a.nstrips-std::abs(hitPredPos_)) <  std::abs(b.nstrips-std::abs(hitPredPos_)); });  
        hitTRNearestStrips_ = peaksTR.front().nstrips;
        hitTRNearestNSat_   = peaksTR.front().nsat;
        hitTRNearestCharge_ = peaksTR.front().charge * norm; 
      } else {
        hitTRNearestStrips_ = hitTRSmallestStrips_ = hitTRLargestStrips_ = -1;
      }

      unsigned int thisSat = (ampls[0] >= 254);
      hitSaturatedStrips_ = thisSat;
      for (unsigned int i = 1, n = ampls.size(); i < n; ++i) {
        if (ampls[i] >= 254) {
            thisSat++;
        } else if (thisSat > 0) {
            hitSaturatedStrips_ = std::max<int>(hitSaturatedStrips_, thisSat);
            thisSat = 0;
        }
      }
      if (thisSat > 0) {
          hitSaturatedStrips_ = std::max<int>(hitSaturatedStrips_, thisSat);
      }

      if (hitUsable_) {
          std::vector<Peak> peaks;
         
          if (hitSaturatedStrips_ < int(maxNSat_) && !(std::abs(hitPredPos_)<1.5f && hitStrips_ <= 2)) {
              int tooNarrowTight = std::floor(std::abs(hitPredPos_)-maxTrimmedSizeDiffNegTight_);
              int tooNarrowLoose = std::floor(std::abs(hitPredPos_)-maxTrimmedSizeDiffNeg_);
              int narrowTight = std::max<int>(2,std::ceil(std::abs(hitPredPos_)+maxTrimmedSizeDiffPosTight_));
              int narrowLoose = std::max<int>(2,std::ceil(std::abs(hitPredPos_)+maxTrimmedSizeDiffPos_));
              printf("\nCluster with predicted strips %.2f, observed strips %d (consec. saturated strips: %d), trimmed size %d (tight %d); tooNarrow if < %d(%d), narrow if <= %d(%d) \n", hitPredPos_, hitStrips_, hitSaturatedStrips_, hitStripsTrim_, hitStripsTrimTight, tooNarrowLoose, tooNarrowTight, narrowLoose, narrowTight);
              std::string verdict;
              if (trackAssoc_ == 1 && hitAssoc_ == 1) {
                verdict = "GOOD";
              } else {
                verdict = "BAD";
              }
              printf("  cluster corresponds to a %s match\n", verdict.c_str());
              if (hitStripsTrim_ != int(hitStripsTrimTight)) {
                  dumpCluster(detId(),sfirst,ampls,mip);
              }

              //dumpAssociatedClusters(detId(),cluster,mip);
              if (hitStripsTrim_ < std::floor(std::abs(hitPredPos_)-maxTrimmedSizeDiffNegTight_)) {
                  printf(" --> %s REJECTED (too narrow, tight cut)\n", verdict.c_str());
                  theTree->Fill(); return false;
              } else if (hitStripsTrim_ < std::floor(std::abs(hitPredPos_)-maxTrimmedSizeDiffNeg_)) {
                  printf(" --> %s REJECTED (too narrow)\n", verdict.c_str());
                  dumpCluster(detId(),sfirst,ampls,mip);
                  theTree->Fill(); return false;
              } else if (hitStripsTrim_ <= std::max<int>(2,std::ceil(std::abs(hitPredPos_)+maxTrimmedSizeDiffPosTight_))) {
                  printf(" --> %s ACCEPTED (after trimming, tight cut)\n", verdict.c_str());
                  theTree->Fill(); return true;
              } else if (hitStripsTrim_ <= std::max<int>(2,std::ceil(std::abs(hitPredPos_)+maxTrimmedSizeDiffPos_))) {
                  printf(" --> %s ACCEPTED (after trimming)\n", verdict.c_str());
                  dumpCluster(detId(),sfirst,ampls,mip);
                  theTree->Fill(); return true;
              } else {
                  float mipnorm = mip/std::abs(ldir.z());
                  utils::SlidingPeakFinder pf(std::max<int>(2,std::ceil(std::abs(hitPredPos_)+subclusterWindow_)));
                  PeakFinderTest pfCutTight(mipnorm, detId(), cluster.firstStrip(), theNoise, seedCutMIPs_, seedCutSN_, subclusterCutMIPsTight_, subclusterCutSNTight_);
                  PeakFinderTest pfCut(mipnorm, detId(), cluster.firstStrip(), theNoise, seedCutMIPs_, seedCutSN_, subclusterCutMIPs_, subclusterCutSN_);

                  if (!pf.apply(cluster.amplitudes(), pfCutTight)) {
                      dumpCluster(detId(),sfirst,ampls,mip);
                      printf("running sliding peak finder with window %d strips\n", std::max<int>(2,std::ceil(std::abs(hitPredPos_)+subclusterWindow_)));
                      printf("  expected charge due to track angle: %.2f normal MIPs = %.1f\n", 1.0/std::abs(ldir.z()), mipnorm);
                      PeakFinderTest pfTest(mipnorm, detId(), cluster.firstStrip(), theNoise, seedCutMIPs_, seedCutSN_, subclusterCutMIPs_, subclusterCutSN_, true, false);
                      pf.apply(cluster.amplitudes(), pfTest, false, sfirst); 

                      if (!pf.apply(cluster.amplitudes(), pfCut)) {
                        printf(" --> %s REJECTED (too large)\n", verdict.c_str());
                        theTree->Fill(); return false;
                      } else {
                        printf(" --> %s ACCEPTED (subcluster)\n", verdict.c_str());
                        theTree->Fill(); return true;
                      }
                  } else {
                      printf(" --> %s ACCEPTED (subcluster, tight cut)\n", verdict.c_str());
                      theTree->Fill(); return true;
                  }
              }

              //peakFinder(detId(),sfirst,ampls,peaks,0.4, 5., 9., 4., true);
              //thresholdRaiser(detId(),sfirst,ampls,peaks,4., 5., 9., true);
          }
      }

      theTree->Fill();
   }
   return true;
}

/*****************************************************************************/
bool ClusterShapeDebugTrajectoryFilter::qualityFilter
  (const Trajectory& trajectory) const
{
   Trajectory::DataContainer tms = trajectory.measurements();

   std::vector<Id2> work;
   std::map<Id2,std::set<Id2> > scoreboard;
   std::set<Id2> countlayers;
   for(Trajectory::DataContainer::const_iterator tm = tms.begin(); tm!= tms.end(); ++tm)
   {
       const TrackingRecHit* hit = tm->recHit()->hit();
       if(hit == 0 || !hit->isValid() || hit->geographicalId().rawId() == 0) continue;

       Id2 mylayer = Id2(hit->geographicalId().subdetId(), theTopology->layer(hit->geographicalId()));
       countlayers.insert(mylayer);

       fillAssociations(hit, work);
       for (Id2 simtk: work) {
           scoreboard[simtk].insert(mylayer);
       }
   }

   Id2 bestSim; unsigned int layers = 0;
   for (const auto &p : scoreboard) {
      if (p.second.size() > layers) {
         layers = p.second.size();
         bestSim = p.first;
      }
   }
   EncodedEventId evid(bestSim.first);

   if (layers >= 3 && layers >= 0.67*countlayers.size()) {
      trackAssoc_ = (evid.bunchCrossing() == 0 ? 1 : -1); 
   } else {
      trackAssoc_ = 0;
   }

   if (maxBadHits_ > 1) {
   printf("======== QUALITYFILTER called for Trajectory with %d/%d associated layers ======\n", layers, int(countlayers.size()));
   unsigned int badHits = 0;
   for (auto i = tms.begin(), e = tms.end(); i != e; ++i) {
       const TrackingRecHit* hit = i->recHit()->hit();
       if (!i->updatedState().isValid()) continue;
       if (hit == 0 || !hit->isValid()) continue;
       if (hit->geographicalId().subdetId() <= 2) continue; // we look only at strips for now
       bool theOk = fillCluster(hit, i->updatedState(), bestSim);
       if (!theOk) badHits++;
   }
   printf("======== QUALITYFILTER called for Trajectory with %d/%d associated layers --> %d bad hits ======\n", layers, int(countlayers.size()), badHits);
   }
  
   return true;
}

/*****************************************************************************/
bool ClusterShapeDebugTrajectoryFilter::qualityFilter
  (const TempTrajectory& trajectory) const
{
   TempTrajectory::DataContainer tms = trajectory.measurements();

   std::vector<Id2> work;
   std::map<Id2,std::set<Id2> > scoreboard;
   std::set<Id2> countlayers;
   for(TempTrajectory::DataContainer::const_iterator tm = tms.rbegin(); tm!= tms.rend(); --tm)
   {
       const TrackingRecHit* hit = tm->recHit()->hit();
       if(hit == 0 || !hit->isValid() || hit->geographicalId().rawId() == 0) continue;

       Id2 mylayer = Id2(hit->geographicalId().subdetId(), theTopology->layer(hit->geographicalId()));
       countlayers.insert(mylayer);

       fillAssociations(hit, work);
       for (Id2 simtk: work) {
           scoreboard[simtk].insert(mylayer);
       }
   }

   Id2 bestSim; unsigned int layers = 0;
   for (const auto &p : scoreboard) {
      if (p.second.size() > layers) {
         layers = p.second.size();
         bestSim = p.first;
      }
   }
   EncodedEventId evid(bestSim.first);

   if (layers >= 3 && layers >= 0.67*countlayers.size()) {
      trackAssoc_ = (evid.bunchCrossing() == 0 ? 1 : -1); 
   } else {
      trackAssoc_ = 0;
   }

   if (maxBadHits_ > 1) {
   printf("======== QUALITYFILTER called for TempTrajectory with %d/%d associated layers ======\n", layers, int(countlayers.size()));
   unsigned int badHits = 0;
   for (auto i = tms.rbegin(), e = tms.rend(); i != e; --i) {
       const TrackingRecHit* hit = i->recHit()->hit();
       if (!i->updatedState().isValid()) continue;
       if (hit == 0 || !hit->isValid()) continue;
       if (hit->geographicalId().subdetId() <= 2) continue; // we look only at strips for now
       bool theOk = fillCluster(hit, i->updatedState(), bestSim);
       if (!theOk) badHits++;
   }
   printf("======== QUALITYFILTER called for TempTrajectory with %d/%d associated layers --> %d bad hits ======\n", layers, int(countlayers.size()), badHits);
   }
  
   return true;
}

void ClusterShapeDebugTrajectoryFilter::setEvent
  (const edm::Event &event, const edm::EventSetup &es) 
{
  // Get tracker geometry
  edm::ESHandle<TrackerGeometry> tracker;
  es.get<TrackerDigiGeometryRecord>().get(tracker);
  theTracker = tracker.product();

  //
  //  theClusterShape = new ClusterShapeHitFilter(es);
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
  bool reindex = (theMTEvent == 0);
  theMTEvent = mte.product();
  if (reindex) indexStripDets();

  theRecoAssociator.reset(new TrackMixingAssociator());
  edmNew::DetSetVector<SiPixelCluster> dummy;
  edm::Handle<edmNew::DetSetVector<SiStripCluster>> handle;
  int icluster = 0;
  for (const auto & t : clusterTokens_) {
    event.getByToken(t, handle);
    theRecoAssociator->registerClusterEvent(++icluster, dummy, *handle);
  }
  initAssociator(event, es);
}

void ClusterShapeDebugTrajectoryFilter::initAssociator(const edm::Event &event, const edm::EventSetup &)
{
  theAssociator.reset(new TrackerHitAssociator(event, theConfig));
}


void ClusterShapeDebugTrajectoryFilter::indexStripDets
   () const 
{
   theStripDetLookup.clear();
   const StMeasurementConditionSet & stripData = theMTEvent->stripData().conditions();
   for (unsigned int i = 0, n = stripData.nDet(); i < n; ++i) {
      theStripDetLookup[stripData.id(i)] = i;
   }
}

bool ClusterShapeDebugTrajectoryFilter::hasBadStrip
   (unsigned int detid, int strip) const 
{
   const StMeasurementConditionSet & stripData = theMTEvent->stripData().conditions();
   int index = theStripDetLookup[detid];
   if (stripData.id(index) != detid) { indexStripDets(); index = theStripDetLookup[detid]; }
   if (strip < 0 || strip >= stripData.totalStrips(index)) return false;
   if (stripData.bad128Strip(index, strip)) return true;
   for (const StMeasurementConditionSet::BadStripBlock & bsb : stripData.badStripBlocks(index)) {
      if (bsb.first <= strip && bsb.last >= strip) return true;
   }
   return false;
}

void ClusterShapeDebugTrajectoryFilter::peakFinder(int detid, int ifirst, const std::vector<uint8_t> & ampls, std::vector<Peak> &peaks, float alpha, float seedThr, float clustThr, float stepThr, bool debug) const {
    int total = 0;
    for (uint8_t a : ampls) total += a;
    SiStripNoises::Range noises = theNoise->getRange(detid);
    peaks.clear();
    if (debug) printf("running peak finder with alpha = %.2f, thresholds = %.1f (seed), %.1f (cluster), %.1f (step)\n", alpha, seedThr, clustThr, stepThr);
    if (debug) printf("  1) make list of local maxima\n");
    for (unsigned int i = 0, n = ampls.size(); i < n; ++i) {
        if ((i == 0 || ampls[i-1] < ampls[i]) && (i == n-1 || ampls[i+1] <= ampls[i])) {
            float ston = unsigned(ampls[i])/float(theNoise->getNoise(i+ifirst,noises));
            if (ston > seedThr) {
                peaks.emplace_back(Peak(i,ampls[i]));
                if (debug) printf("       - strip at %d (S/N = %.1f)\n",ifirst+i, ston);
            }
        }
    }
    if (debug) printf("  2) sort and process list\n");
    std::sort(peaks.begin(), peaks.end());
    for (unsigned int ip = 0, np = peaks.size(); ip  < np; ++ip) {
        Peak &me = peaks[ip];
        if (me.maxval == 0) continue; // ignore peaks that have been already killed
        uint8_t floor = alpha*std::min<uint16_t>(me.maxval,170);
        //printf("       - starting from strip at %d (adc %d, floor %d)\n",ifirst+int(me.first),int(me.maxval),int(floor));
        unsigned int lbound = 0, rbound = ampls.size();
        for (unsigned int ip2 = 0; ip2 < ip; ++ip2) {
            if (peaks[ip2].first < me.first) lbound = std::max<unsigned int>(lbound, peaks[ip2].first + peaks[ip2].nstrips);
            else if (peaks[ip2].first > me.first) rbound = std::min<unsigned int>(rbound, peaks[ip2].first);
        }
        //printf("       - bounds set to [%d,%d[\n",ifirst+int(lbound),ifirst+int(rbound));
        unsigned int ileft = me.first, iright = me.first;
        while (ileft > lbound && ampls[ileft-1] > floor) {
            // check if I have moved upwards significantly
            if (ampls[ileft-1] > ampls[ileft] && ampls[ileft-1] > ampls[ileft] + stepThr * std::max(theNoise->getNoise(ileft-1+ifirst, noises),theNoise->getNoise(ileft+ifirst, noises))) {
                //printf("       - stop left scan at %d: gain is %d, gain/noise = %.1f\n",ifirst+int(ileft),int(ampls[ileft-1]-ampls[ileft]), float(ampls[ileft-1]-ampls[ileft])/std::max(theNoise->getNoise(ileft-1+ifirst, noises),theNoise->getNoise(ileft+ifirst, noises)));
                break; 
            }
            ileft--;
            me.charge += ampls[ileft];
            me.nsat += (ampls[ileft] >= 254);
        }
        //printf("       - end of left scan at %d\n",ifirst+int(ileft));
        while (iright+1 < rbound && ampls[iright+1] > floor) {
            // check if I have moved upwards significantly
            if (ampls[iright+1] > ampls[iright] && ampls[iright+1] > ampls[iright] + stepThr * std::max(theNoise->getNoise(iright+ifirst, noises), theNoise->getNoise(iright+1+ifirst,noises))) {
                //printf("       - stop right scan at %d: gain is %d, gain/noise = %.1f\n",ifirst+int(iright),int(ampls[iright+1]-ampls[iright]), float(ampls[iright+1]-ampls[iright])/std::max(theNoise->getNoise(iright+1+ifirst, noises),theNoise->getNoise(iright+ifirst, noises)));
                break; 
            }
            iright++;
            me.charge += ampls[iright];
            me.nsat += (ampls[iright] >= 254);
        }
        //printf("       - end of right scan at %d\n",ifirst+int(iright));
        me.nstrips = iright - ileft + 1;
        me.first = ileft;
        for (unsigned int ip2 = ip+1; ip2 < np; ++ip2) {
            if (ileft <= peaks[ip2].first  && peaks[ip2].first <= iright) peaks[ip2].maxval = 0; // kill
        }
    }
    // compute S/N for the clusters
    for (Peak &p : peaks) {
        if (p.maxval == 0) continue;
        double c = 0, n2 = 0;
        for (unsigned int i = p.first, n = i+p.nstrips; i < n; ++i) {
            c += ampls[i];
            n2 += std::pow(theNoise->getNoise(i+ifirst,noises),2);
        }
        if (debug) printf("   - peak [%d,%d]: strips %d, maxval %d, charge %d, nsat %d, S/N = %.1f, peak/cluster = %.2f\n", ifirst+int(p.first), ifirst+int(p.first+p.nstrips-1), int(p.nstrips), int(p.maxval), int(p.charge), int(p.nsat), c/std::sqrt(n2), p.charge/float(total));
        if (c/std::sqrt(n2) < clustThr) p.maxval = 0;
    }
    // resort to bring killed clusters to the bottom
    std::sort(peaks.begin(), peaks.end());
    while(!peaks.empty() && peaks.back().maxval == 0) peaks.pop_back();
    int peak_total = 0;
    for (const Peak &p : peaks) peak_total += p.charge;
    if (debug) printf("   - total charge in peaks / cluster charge = %.2f\n", peak_total/float(total));
}

void ClusterShapeDebugTrajectoryFilter::thresholdRaiser(int detid, int ifirst, const std::vector<uint8_t> & ampls, std::vector<Peak> &peaks, float stripThr, float seedThr, float clustThr, bool debug) const {
    SiStripNoises::Range noises = theNoise->getRange(detid);
    peaks.clear();
    if (debug) printf("running thredholdRaiser with thresholds = %.1f (strip), %.1f (seed), %.1f (cluster)\n", stripThr, seedThr, clustThr);
    Peak current; double noise2 = 0; float seedStoN = 0;
    for (unsigned int i = 0, n = ampls.size(); i < n; ++i) {
        float mynoise = float(theNoise->getNoise(i+ifirst,noises));
        float ston = unsigned(ampls[i])/mynoise;
        if (ston < stripThr) {
            if (current.nstrips != 0) { 
                if (seedStoN > seedThr && current.charge > clustThr*std::sqrt(noise2)) {
                    peaks.push_back(current); 
                    if (debug) printf("   - peak [%d,%d]: strips %d, maxval %d, charge %d, nsat %d, S/N = %.1f\n", ifirst+int(current.first), ifirst+int(current.first+current.nstrips-1), int(current.nstrips), int(current.maxval), int(current.charge), int(current.nsat), current.charge/std::sqrt(noise2));
                }
                current = Peak(); noise2 = 0; seedStoN = 0;
            }
        } else if (current.nstrips == 0){
            current  = Peak(i, ampls[i]);
            noise2   = mynoise*mynoise;
            seedStoN = ston;
        } else {
            current.nstrips++;
            current.charge += ampls[i];
            current.nsat   += (ampls[i] >= 254);
            current.maxval = std::max(current.maxval, ampls[i]);
            seedStoN = std::max(seedStoN, ston);
            noise2 += mynoise*mynoise;
        }
    }
    if (current.nstrips != 0) {
        if (seedStoN > seedThr && current.charge > clustThr*std::sqrt(noise2)) {
            peaks.push_back(current);
            if (debug) printf("   - peak [%d,%d]: strips %d, maxval %d, charge %d, nsat %d, S/N = %.1f\n", ifirst+int(current.first), ifirst+int(current.first+current.nstrips-1), int(current.nstrips), int(current.maxval), int(current.charge), int(current.nsat), current.charge/std::sqrt(noise2));
        }
    }
}

void ClusterShapeDebugTrajectoryFilter::subtractPeaks(int detid, int ifirst, const std::vector<uint8_t> & ampls, const std::vector<Peak> &peaksIn, std::vector<uint8_t> & amplsOut, std::vector<Peak> &peaksOut, int maxSize, float stripThr, float seedThr, float clustThr, bool debug) const 
{
    if (debug) printf("subtracting peaks from the PeakFinder run\n");
    int total = 0;
    for (uint8_t a : ampls) total += a;
    peaksOut = peaksIn;
    amplsOut = ampls;
    for (unsigned int i = 0, n = peaksIn.size(); i < n; ++i) {
        const Peak &p = peaksIn[i];
        if (p.nstrips > maxSize) continue;
        Peak &pOut = peaksOut[i];
        uint8_t baseline = ((p.first == 0                      ? 0 : std::min<uint8_t>(ampls[p.first-1],ampls[p.first])) + 
                            (p.first+p.nstrips == ampls.size() ? 0 : std::min<uint8_t>(ampls[p.first+p.nstrips], ampls[p.first+p.nstrips-1])))/2;
        pOut.charge = 0; 
        pOut.maxval = 0;
        for (int j = 0; j < p.nstrips; ++j) {
            if (amplsOut[p.first+j] >= 254) {
                pOut.charge += amplsOut[p.first+j];
                pOut.maxval = std::max<uint8_t>(pOut.maxval,amplsOut[p.first+j]);
            } else {
                pOut.charge += amplsOut[p.first+j] - baseline;
                pOut.maxval = std::max<uint8_t>(pOut.maxval,amplsOut[p.first+j] - baseline);
            } 
            amplsOut[p.first+j] = baseline;
        }
        if (debug) printf("   - peak [%d,%d]: strips %d, baseline %d, charge %d -> %d, peak/cluster = %.2f -> %.2f\n", ifirst+int(p.first), ifirst+int(p.first+p.nstrips-1), int(p.nstrips), int(baseline), int(p.charge), int(pOut.charge), p.charge/float(total), pOut.charge/float(total));
    }
    std::sort(peaksOut.begin(), peaksOut.end());
    int peak_total = 0;
    for (const Peak &p : peaksOut) peak_total += p.charge;
    if (debug) printf("   - total charge in peaksOut / cluster charge = %.2f\n", peak_total/float(total));
}

void ClusterShapeDebugTrajectoryFilter::dumpCluster(int detid, int ifirst, const std::vector<uint8_t> & ampls, double mip) const
{
    SiStripNoises::Range noises = theNoise->getRange(detid);
    double ctotal = 0, n2 = 0, bestSN=0;
    for (unsigned int i = 0, n = ampls.size(); i < n; ++i) {
        double myn = theNoise->getNoise(i+ifirst, noises);
        printf("\t%4u: %4d (S/N = %6.1f): ", i+ifirst, int(ampls[i]), float(ampls[i])/myn);
        for (uint8_t j = 0; j < ampls[i]/2; ++j) putchar('=');
        if (ampls[i] % 2 == 1) putchar('-');
        if (ampls[i] >= 254) putchar('X');
        printf("\n");
        n2 += myn*myn;
        ctotal += ampls[i];
        bestSN = std::max(bestSN, float(ampls[i])/myn);
    }
    printf("    Total charge: %.0f (%.2f normal MIPs), S/N = %.1f (highest strip: %.1f)\n", ctotal, ctotal/mip, ctotal/std::sqrt(n2), bestSN);
}

void ClusterShapeDebugTrajectoryFilter::dumpAssociatedClusters(int detid, const SiStripCluster &cluster, double mip) const
{
    std::vector<std::pair<int, const SiStripCluster *>> assocs;
    theRecoAssociator->getAssociatedClusters(detid, cluster, assocs);
    printf("    Associated clusters (%ld)\n", assocs.size());
    for (const auto &p : assocs) {
        printf("    - from %d:\n", p.first);
        dumpCluster(detid, p.second->firstStrip(), p.second->amplitudes(), mip);
    }
}

