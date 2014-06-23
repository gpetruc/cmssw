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

#include "CondFormats/SiStripObjects/interface/SiStripNoises.h"
#include "CondFormats/DataRecord/interface/SiStripNoisesRcd.h"

#include <TTree.h>

/*****************************************************************************/
ClusterShapeDebugTrajectoryFilter::ClusterShapeDebugTrajectoryFilter
   (const edm::ParameterSet &iCfg, edm::ConsumesCollector& iC) :
   theConfig(iCfg.getParameter<edm::ParameterSet>("HitAssociatorPSet")),
   mteToken_(iC.consumes<MeasurementTrackerEvent>(iCfg.getParameter<edm::InputTag>("MeasurementTrackerEvent"))),
   theAssociator(0),
   theFilter(0),
   theTopology(0),
   theTracker(0),
   theTree(0)
{
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
   assert(theAssociator != 0);
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
   if (hit->geographicalId().subdetId() <= 2) return true; // we look only at strips for now

   fillCluster(hit, last.updatedState(), bestSim);
   return true;
}

void ClusterShapeDebugTrajectoryFilter::fillAssociations
   (const TrackingRecHit *hit, std::vector<Id2> &out) const
{
   out.clear();
   std::vector<PSimHit> simHits = theAssociator->associateHit(*hit);
   hitSingleSim_ = (simHits.size() == 1);
   for (const PSimHit &hit : simHits) {
      Id2 id2(hit.eventId().rawId(), hit.trackId());
      out.push_back(id2);
#if 0
      // Outgoing?
      DetId id = DetId(simHit.detUnitId());
      const GeomDetUnit *gdu = theTracker->idToDetUnit(id);
      if (gdu == 0) throw cms::Exception("MissingData") << "Missing DetUnit for detid " << id.rawId() << "\n" << std::endl;
      GlobalVector gvec = theTracker->idToDetUnit(id)->position() - GlobalPoint(0,0,0);
      LocalVector  lvec = theTracker->idToDetUnit(id)->toLocal(gvec);
      LocalVector  ldir = simHit.exitPoint() - simHit.entryPoint();
      bool isOutgoing = (lvec.z()*ldir.z() > 0);

      // From a relevant process? primary or decay
      bool isRelevant = (simHit.processType() == 2 || simHit.processType() == 4);

      // Fast enough? pt > 50 MeV/c
      bool isFast = (simHit.momentumAtEntry().perp() > 0.050);

      if (isOutgoing && isRelevant && isFast)
#endif
   }
}
 
void ClusterShapeDebugTrajectoryFilter::fillCluster
   (const TrackingRecHit *hit, const TrajectoryStateOnSurface &tsos, Id2 simtk, bool mustProject) const
{
   const TrackerSingleRecHit *stripHit = 0;
   if (typeid(*hit) == typeid(SiStripMatchedRecHit2D)) {
      const SiStripMatchedRecHit2D & mhit = static_cast<const SiStripMatchedRecHit2D &>(*hit);
      SiStripRecHit2D mono = mhit.monoHit();
      SiStripRecHit2D stereo = mhit.stereoHit();
      fillCluster(&mono,   tsos, simtk, true);
      fillCluster(&stereo, tsos, simtk, true);
   } else if (typeid(*hit) == typeid(ProjectedSiStripRecHit2D)) {
      const ProjectedSiStripRecHit2D & mhit = static_cast<const ProjectedSiStripRecHit2D &>(*hit);
      const SiStripRecHit2D &orig = mhit.originalHit();
      fillCluster(&orig,   tsos, simtk, true);
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
      

      const StripGeomDetUnit* stripDetUnit = dynamic_cast<const StripGeomDetUnit *>(det);
      if (det == 0) throw cms::Exception("Strip not a StripGeomDetUnit?") << " on " << detId.rawId() << " aka " << theTopology->print(detId) << "\n";

      float MeVperADCStrip =  3.61e-06*265; 
      float norm = MeVperADCStrip/stripDetUnit->surface().bounds().thickness();
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


      theTree->Fill();
      static int i_good = 0, i_bad = 0;
      if (hitUsable_ && trackAssoc_ == 1 && hitAssoc_ == 1) {
          SiStripNoises::Range noises = theNoise->getRange(detId());
          std::vector<Peak> peaks;
          if (hitCompatPos_) {
              if (i_good < 200) {
                  i_good++;
                  printf("\nCluster with predicted strips %.2f, observed strips %d: GOOD\n", hitPredPos_, hitStrips_);
                  for (unsigned int i = 0, n = ampls.size(); i < n; ++i) {
                      printf("\t%3u: %4d (S/N = %6.1f): ", i, int(ampls[i]), float(ampls[i])/theNoise->getNoise(i+sfirst, noises));
                      for (uint8_t j = 0; j < ampls[i]/2; ++j) putchar('=');
                      if (ampls[i] % 2 == 1) putchar('-');
                      if (ampls[i] >= 254) putchar('X');
                      printf("\n");
                  }
                  peakFinder(detId(),sfirst,ampls,peaks,0.4, 5., 9., 4., true);
                  thresholdRaiser(detId(),sfirst,ampls,peaks,4., 5., 9., true);
              }
          } else {    
              if (i_bad < 200) {
                  i_bad++;
                  printf("\nCluster with predicted strips %.2f, observed strips %d: FAIL\n", hitPredPos_, hitStrips_);
                  for (unsigned int i = 0, n = ampls.size(); i < n; ++i) {
                      printf("\t%3u: %4d (S/N = %6.1f): ", i, int(ampls[i]), float(ampls[i])/theNoise->getNoise(i+sfirst, noises));
                      for (uint8_t j = 0; j < ampls[i]/2; ++j) putchar('=');
                      if (ampls[i] % 2 == 1) putchar('-');
                      if (ampls[i] >= 254) putchar('X');
                      printf("\n");
                  }
                  peakFinder(detId(),sfirst,ampls,peaks,0.4,  5., 9.,  4., true);
                  thresholdRaiser(detId(),sfirst,ampls,peaks, 5., 7., 12., true);
              }
          } 
      }
   }
}

/*****************************************************************************/
bool ClusterShapeDebugTrajectoryFilter::qualityFilter
  (const Trajectory& trajectory) const
{
  return true;
}

/*****************************************************************************/
bool ClusterShapeDebugTrajectoryFilter::qualityFilter
  (const TempTrajectory& trajectory) const
{
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
  
  // create the hit associator
  delete theAssociator;
  theAssociator = new TrackerHitAssociator(event, theConfig);

  edm::Handle<MeasurementTrackerEvent> mte;
  event.getByToken(mteToken_, mte);
  bool reindex = (theMTEvent == 0);
  theMTEvent = mte.product();
  if (reindex) indexStripDets();
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
    SiStripNoises::Range noises = theNoise->getRange(detid);
    peaks.clear();
    if (debug) printf("running peak finder with alpha = %.2f, thresholds = %.1f (seed), %.1f (cluster), %.1f (step)\n", alpha, seedThr, clustThr, stepThr);
    if (debug) printf("  1) make list of local maxima\n");
    for (unsigned int i = 0, n = ampls.size(); i < n; ++i) {
        if ((i == 0 || ampls[i-1] < ampls[i]) && (i == n-1 || ampls[i+1] <= ampls[i])) {
            float ston = unsigned(ampls[i])/float(theNoise->getNoise(i+ifirst,noises));
            if (ston > seedThr) {
                peaks.emplace_back(Peak(i,ampls[i]));
                if (debug) printf("       - strip at %d (S/N = %.1f)\n",i, ston);
            }
        }
    }
    if (debug) printf("  2) sort and process list\n");
    std::sort(peaks.begin(), peaks.end());
    for (unsigned int ip = 0, np = peaks.size(); ip  < np; ++ip) {
        Peak &me = peaks[ip];
        if (me.maxval == 0) continue; // ignore peaks that have been already killed
        uint8_t floor = alpha*std::min<uint16_t>(me.maxval,170);
        //printf("       - starting from strip at %d (adc %d, floor %d)\n",int(me.first),int(me.maxval),int(floor));
        unsigned int lbound = 0, rbound = ampls.size();
        for (unsigned int ip2 = 0; ip2 < ip; ++ip2) {
            if (peaks[ip2].first < me.first) lbound = std::max<unsigned int>(lbound, peaks[ip2].first + peaks[ip2].nstrips);
            else if (peaks[ip2].first > me.first) rbound = std::min<unsigned int>(rbound, peaks[ip2].first);
        }
        //printf("       - bounds set to [%d,%d[\n",int(lbound),int(rbound));
        unsigned int ileft = me.first, iright = me.first;
        while (ileft > lbound && ampls[ileft-1] > floor) {
            // check if I have moved upwards significantly
            if (ampls[ileft-1] > ampls[ileft] && ampls[ileft-1] > ampls[ileft] + stepThr * std::max(theNoise->getNoise(ileft-1+ifirst, noises),theNoise->getNoise(ileft+ifirst, noises))) {
                //printf("       - stop left scan at %d: gain is %d, gain/noise = %.1f\n",int(ileft),int(ampls[ileft-1]-ampls[ileft]), float(ampls[ileft-1]-ampls[ileft])/std::max(theNoise->getNoise(ileft-1+ifirst, noises),theNoise->getNoise(ileft+ifirst, noises)));
                break; 
            }
            ileft--;
            me.charge += ampls[ileft];
            me.nsat += (ampls[ileft] >= 254);
        }
        //printf("       - end of left scan at %d\n",int(ileft));
        while (iright+1 < rbound && ampls[iright+1] > floor) {
            // check if I have moved upwards significantly
            if (ampls[iright+1] > ampls[iright] && ampls[iright+1] > ampls[iright] + stepThr * std::max(theNoise->getNoise(iright+ifirst, noises), theNoise->getNoise(iright+1+ifirst,noises))) {
                //printf("       - stop right scan at %d: gain is %d, gain/noise = %.1f\n",int(iright),int(ampls[iright+1]-ampls[iright]), float(ampls[iright+1]-ampls[iright])/std::max(theNoise->getNoise(iright+1+ifirst, noises),theNoise->getNoise(iright+ifirst, noises)));
                break; 
            }
            iright++;
            me.charge += ampls[iright];
            me.nsat += (ampls[iright] >= 254);
        }
        //printf("       - end of right scan at %d\n",int(iright));
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
        if (debug) printf("   - peak [%d,%d]: strips %d, maxval %d, charge %d, nsat %d, S/N = %.1f\n", int(p.first), int(p.first+p.nstrips-1), int(p.nstrips), int(p.maxval), int(p.charge), int(p.nsat), c/std::sqrt(n2));
        if (c/std::sqrt(n2) < clustThr) p.maxval = 0;
    }
    // resort to bring killed clusters to the bottom
    std::sort(peaks.begin(), peaks.end());
    while(!peaks.empty() && peaks.back().maxval == 0) peaks.pop_back();
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
                    if (debug) printf("   - peak [%d,%d]: strips %d, maxval %d, charge %d, nsat %d, S/N = %.1f\n", int(current.first), int(current.first+current.nstrips-1), int(current.nstrips), int(current.maxval), int(current.charge), int(current.nsat), current.charge/std::sqrt(noise2));
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
            if (debug) printf("   - peak [%d,%d]: strips %d, maxval %d, charge %d, nsat %d, S/N = %.1f\n", int(current.first), int(current.first+current.nstrips-1), int(current.nstrips), int(current.maxval), int(current.charge), int(current.nsat), current.charge/std::sqrt(noise2));
        }
    }
}
