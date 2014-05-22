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
   theTree->Branch("hitCharge", &hitCharge_, "hitCharge/F");
   theTree->Branch("hitChargeNorm", &hitChargeNorm_, "hitChargeNorm/F");
   theTree->Branch("hitChargeRms", &hitChargeRms_, "hitChargeRms/F");
   theTree->Branch("hitPredPos", &hitPredPos_, "hitPredPos/F");
   theTree->Branch("hitPredNoPos", &hitPredNoPos_, "hitPredNoPos/F");
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
      hitChargeNorm_ = hitCharge_ / ldir.z();
      theTree->Fill();
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
