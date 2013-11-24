//
//

/**
  \class    TestHT TestHT.h "RecoTracker/HTPattern/interface/TestHT.h"
  \brief    Matcher of reconstructed objects to other reconstructed objects using the tracks inside them 
            
  \author   Giovanni Petrucciani
*/


#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "RecoTracker/HTPattern/interface/HTHits.h"
#include "RecoTracker/HTPattern/interface/HTHitMap.h"
#include "RecoTracker/HTPattern/interface/HTDebugger.h"
#include "RecoTracker/HTPattern/interface/TrackCandidateBuilderFromCluster.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#ifdef  HTDEBUG
#define DEBUG HTDEBUG
#else
#define DEBUG 0.5
#endif
#define DEBUG1_printf  if (DEBUG>=1) printf
#define DEBUG2_printf  if (DEBUG>=2) printf
#define DEBUG3_printf  if (DEBUG>=3) printf


class TestHT : public edm::EDProducer {
    public:
      explicit TestHT(const edm::ParameterSet & iConfig);
      virtual ~TestHT() { }

      virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup) override;

      void getPixelHits3D(const edm::Event & iEvent, const reco::BeamSpot &bsp, HTHits3D &hits3d) const ;
      void getStripHits3D(const edm::Event & iEvent, const reco::BeamSpot &bsp, HTHits3D &hits3d) const ;
      void getStripHits2D(const edm::Event & iEvent, const reco::BeamSpot &bsp, HTHits3D &hits3d) const ;
      //std::auto_ptr<HTHitMap> getHitMap(const edm::Event & iEvent, const std::vector<bool> & toSkip) const ;
      //void fillHitMap(HTHitMap &map, const HTHitsSpher &hitsSpher, const std::vector<bool> & toSkip) const ;
      
    private:
      typedef edm::ContainerMask<edmNew::DetSetVector<SiPixelCluster> > PixelMask;
      typedef edm::ContainerMask<edmNew::DetSetVector<SiStripCluster> > StripMask;

      /// Labels for input collections
      edm::EDGetTokenT<edmNew::DetSetVector<SiPixelRecHit> > pixelHits_;
      edm::EDGetTokenT<edmNew::DetSetVector<SiStripRecHit2D> > stripHits_;
      edm::EDGetTokenT<edmNew::DetSetVector<SiStripMatchedRecHit2D> > stripHits2D_;
      bool hasMasks_;
      edm::EDGetTokenT<std::vector<reco::Vertex> > vertices_;
      edm::EDGetTokenT<reco::BeamSpot> beamSpot_;
      edm::EDGetTokenT<std::vector<reco::Track> > tracks_;
      edm::EDGetTokenT<MeasurementTrackerEvent> mtEvent_;
      edm::EDGetTokenT<PixelMask> pixelHitMask_;  
      edm::EDGetTokenT<StripMask> stripHitMask_;  


      // Configurables
      bool     seed2d_, seed3d_, seedMixed_;
      uint32_t etabins2d_, etabins3d_, phibins2d_, phibins3d_;
      uint32_t layerSeedCut2d_, layerSeedCut3d_;
      uint32_t layerCut2d_, layerCut3d_;
      uint32_t layerMoreCut_;
      std::vector<double> ptSteps_;
      std::vector<double> ptEdges_;

      // Helpers
      StringCutObjectSelector<reco::Vertex> vertexSelection_;
      TrackCandidateBuilderFromCluster tcBuilder_;

      // EventSetup stuff
      edm::ESHandle<TrackerGeometry> geometry_;
      edm::ESHandle<TrackerTopology> tTopo_;
      edm::ESHandle<MagneticField> bfield_;

      // Debugging
      bool debugger_;
};

TestHT::TestHT(const edm::ParameterSet & iConfig) :
    pixelHits_(consumes<edmNew::DetSetVector<SiPixelRecHit> >(iConfig.getParameter<edm::InputTag>("pixelHits"))),
    stripHits_(consumes<edmNew::DetSetVector<SiStripRecHit2D> >(iConfig.getParameter<edm::InputTag>("stripHits"))),
    stripHits2D_(consumes<edmNew::DetSetVector<SiStripMatchedRecHit2D> >(iConfig.getParameter<edm::InputTag>("stripHits2D"))),
    hasMasks_(iConfig.existsAs<edm::InputTag>("clustersToSkip")),
    vertices_(consumes<std::vector<reco::Vertex> >(iConfig.getParameter<edm::InputTag>("vertices"))),
    beamSpot_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
    tracks_(consumes<std::vector<reco::Track> >(iConfig.getParameter<edm::InputTag>("tracks"))),
    mtEvent_(consumes<MeasurementTrackerEvent>(iConfig.getParameter<edm::InputTag>("measurementTrackerEvent"))),
    seed2d_(iConfig.getParameter<bool>("seed2d")),
    seed3d_(iConfig.getParameter<bool>("seed3d")),
    seedMixed_(iConfig.getParameter<bool>("seedMixed")),
    etabins2d_(iConfig.getParameter<uint32_t>("etabins2d")),
    etabins3d_(iConfig.getParameter<uint32_t>("etabins3d")),
    phibins2d_(iConfig.getParameter<uint32_t>("phibins2d")),
    phibins3d_(iConfig.getParameter<uint32_t>("phibins3d")),
    layerSeedCut2d_(iConfig.getParameter<uint32_t>("layerSeedCut2d")),
    layerSeedCut3d_(iConfig.getParameter<uint32_t>("layerSeedCut3d")),
    layerCut2d_(iConfig.getParameter<uint32_t>("layerCut2d")),
    layerCut3d_(iConfig.getParameter<uint32_t>("layerCut3d")),
    layerMoreCut_(iConfig.getParameter<uint32_t>("layerMoreCut")),
    ptSteps_(iConfig.getParameter<std::vector<double> >("ptSteps")),
    ptEdges_(iConfig.getParameter<std::vector<double> >("ptEdges")),
    vertexSelection_(iConfig.getParameter<std::string>("vertexSelection")),
    tcBuilder_(iConfig.getParameter<edm::ParameterSet>("seedBuilderConfig")),
    debugger_(iConfig.getUntrackedParameter<bool>("debugger",false))
{
    if (hasMasks_) {
        pixelHitMask_ = consumes<PixelMask>(iConfig.getParameter<edm::InputTag>("clustersToSkip"));
        stripHitMask_ = consumes<StripMask>(iConfig.getParameter<edm::InputTag>("clustersToSkip"));
    }
    produces<TrajectorySeedCollection>();  
    produces<TrajectorySeedCollection>("clusters");  
    produces<TrackCandidateCollection>();  
}

void 
TestHT::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {
    using namespace edm;
    using namespace std;

    char evid[50];
    sprintf(evid, "r%d_l%d_e%d_", iEvent.id().run(), iEvent.id().luminosityBlock(), iEvent.id().event());
    std::string prefix(evid);

    //Get the geometry
    iSetup.get<TrackerDigiGeometryRecord>().get(geometry_);

    //Retrieve tracker topology from geometry
    iSetup.get<IdealGeometryRecord>().get(tTopo_);

    //Retrieve magnetic field
    iSetup.get<IdealMagneticFieldRecord>().get(bfield_);
    double bfield = bfield_->inTesla(GlobalPoint(0.,0.,0.)).z();

    edm::Handle<reco::BeamSpot> beamSpot;
    iEvent.getByToken(beamSpot_, beamSpot);

    Handle<vector<reco::Vertex> > vertices;
    iEvent.getByToken(vertices_, vertices);

    Handle<vector<reco::Track> > tracks;
    iEvent.getByToken(tracks_, tracks);
    if (debugger_) { if (DEBUG>=2) HTDebugger::dumpTracks(prefix+"tk", *tracks); } 
    if (DEBUG>0.1) HTDebugger::registerTracks(*tracks); 

    edm::Handle<MeasurementTrackerEvent> mtEvent;
    iEvent.getByToken(mtEvent_, mtEvent);
    std::auto_ptr<MeasurementTrackerEvent> myMtEvent;
    if (hasMasks_) {
        // create a masked measurement tracker
        edm::Handle<StripMask> stripHitMask;
        iEvent.getByToken(stripHitMask_, stripHitMask);
        edm::Handle<PixelMask> pixelHitMask;
        iEvent.getByToken(pixelHitMask_, pixelHitMask);
        myMtEvent.reset(new MeasurementTrackerEvent(*mtEvent, *stripHitMask, *pixelHitMask));
        // initialize the builder
        tcBuilder_.init(iSetup, *myMtEvent, nullptr);
    } else {
        // initialize the builder
        tcBuilder_.init(iSetup, *mtEvent, nullptr);
    }

    std::auto_ptr<TrajectorySeedCollection> out(new TrajectorySeedCollection());
    std::auto_ptr<TrajectorySeedCollection> outCl(new TrajectorySeedCollection());
    std::auto_ptr<TrackCandidateCollection> outTC(new TrackCandidateCollection());

    bool do3d = seed3d_ || seedMixed_;
    bool do2d = seed2d_ || seedMixed_;
    HTHits3D hits3d(beamSpot->position().x(), beamSpot->position().y()); // make pixels
    if (do3d) getPixelHits3D(iEvent, *beamSpot, hits3d);

    HTHits3D hits2d(beamSpot->position().x(), beamSpot->position().y());
    if (do2d) hits2d = hits3d; // start as copy
    if (do2d) getStripHits2D(iEvent, *beamSpot, hits2d);

    if (do3d) getStripHits3D(iEvent, *beamSpot, hits3d);
    if (debugger_) if (DEBUG>=2) HTDebugger::dumpHTHits3D(prefix+"3d", hits3d, vertices->empty() ? 0.0 : vertices->front().z());
    if (debugger_) if (DEBUG>=2) HTDebugger::dumpHTHits3D(prefix+"2d", hits2d, vertices->empty() ? 0.0 : vertices->front().z());

    DEBUG1_printf("Number of hits: %d 3d, %d 2d\n", hits3d.size(), hits2d.size());
    // Now we just do one step
    HTHitMap map3d(etabins3d_,phibins3d_);
    HTHitMap map2d(etabins2d_,phibins2d_);
    std::vector<bool> mask3d(hits3d.size(), false);
    std::vector<bool> mask2d(hits2d.size(), false);
    int etashift = std::round(log(double(etabins3d_)/etabins2d_)/std::log(2));
    int phishift = std::round(log(double(phibins3d_)/phibins2d_)/std::log(2));
    for (unsigned int iptStep = 0, nptSteps = ptSteps_.size(); iptStep < nptSteps; ++iptStep) {
        double ptStep = ptSteps_[iptStep];
        if (ptStep != 0 && bfield == 0) continue;

        for (int ptsign = +1; ptsign > -2; ptsign -= 2) {
            if (ptStep == 0 && ptsign < 0) continue;

            bool firstVertex = true;
            for(const reco::Vertex &vtx : *vertices) {
                // always process first vertex, but apply selection to others
                if (!firstVertex && !vertexSelection_(vtx)) continue; else firstVertex = false;

                sprintf(evid, "r%d_l%d_e%d_pT%c%03d_", iEvent.id().run(), iEvent.id().luminosityBlock(), iEvent.id().event(), ptsign > 0 ? 'p' : 'm', int(ptStep*100));
                std::string prefix(evid);

                float alpha = ptStep ? 0.5 * 0.003 * bfield / ptStep * ptsign : 0;
                if (!ptEdges_.empty()) {
                    float alpha1 =  0.5 * 0.003 * bfield / ptEdges_[2*iptStep] * ptsign;
                    float alpha2 = 0.5 * 0.003 * bfield / ptEdges_[2*iptStep+1] * ptsign;
                    tcBuilder_.setAlphaRange( alpha1, alpha2 );
                    printf("\n========== HT iteration with pT = %+.2f, alpha = %+.5f [ %+.5f : %+.5f ] ==========\n", ptStep*ptsign, alpha, alpha1, alpha2);
                } else {
                    printf("\n========== HT iteration with pT = %+.2f, alpha = %+.5f ==========\n", ptStep*ptsign, alpha);
                }

                HTHitsSpher hits3ds(hits3d, vtx.z(), etabins3d_);
                hits3ds.filliphi(alpha, phibins3d_);
                HTHitsSpher hits2ds(hits2d, vtx.z(), etabins2d_);
                hits2ds.filliphi(alpha, phibins2d_);

                if (debugger_) if (DEBUG>=2) HTDebugger::dumpHTHitsSpher(prefix+"spher3d", hits3ds);
                if (debugger_) if (DEBUG>=2) HTDebugger::dumpHTHitsSpher(prefix+"spher2d", hits2ds);
                // fill the maps
                for (unsigned int i = 0; i < hits3ds.size(); ++i) {
                    if (mask3d[i]) continue;
                    map3d.addHit(hits3ds.ieta(i), hits3ds.iphi(i), i, hits3ds.layermask(i), layerSeedCut3d_);
                }
                // fill the maps
                for (unsigned int i = 0; i < hits2ds.size(); ++i) {
                    if (mask2d[i]) continue;
                    map2d.addHit(hits2ds.ieta(i), hits2ds.iphi(i), i, hits2ds.layermask(i), layerSeedCut2d_);
                }
                // clusterize
                if (do3d) map3d.clusterize();
                if (debugger_) if (DEBUG>=2) HTDebugger::dumpHTHitMap(prefix+"map3d", map3d);

                if (seed3d_) {
                    if (debugger_) if (DEBUG>=2) HTDebugger::dumpHTClusters(prefix+"c3d", map3d, hits3ds, layerCut3d_);
                    tcBuilder_.setHits(hits3ds, hits2ds, map3d, map2d, mask3d, mask2d, etashift, phishift, false);
                    for (auto & cluster : map3d.clusters()) {
                        if (cluster.nmorelayers() >= layerCut3d_) tcBuilder_.run(cluster, layerSeedCut3d_, layerCut3d_, *out, &*outCl);
                    }
                }

                if (seed2d_) {
                    if (debugger_) if (DEBUG>=2) HTDebugger::dumpHTHitMap(prefix+"map2d", map2d);
                    map2d.clusterize();
                    if (debugger_) if (DEBUG>=2) HTDebugger::dumpHTClusters(prefix+"c2d", map2d, hits2ds, layerCut2d_);
                    tcBuilder_.setHits(hits2ds, hits2ds, map2d, map2d, mask2d, mask2d, etashift, phishift, false);
                    for (auto & cluster : map3d.clusters()) {
                        if (cluster.nmorelayers() >= layerCut2d_) tcBuilder_.run(cluster, layerSeedCut2d_, layerCut2d_, *out, &*outCl);
                    }
                }

                if (seedMixed_) {
                    // Recover more layers from 2D map
                    tcBuilder_.setHits(hits3ds, hits2ds, map3d, map2d, mask3d, mask2d, etashift, phishift, true);
                    for (auto & cluster : map3d.clusters()) {
                        cluster.addMoreCell(map2d.get(cluster.ieta() >> etashift, cluster.iphi() >> phishift));
                        const HTCell &cell = map3d.get(cluster.ieta(),cluster.iphi());
                        DEBUG2_printf("cluster #%03d: eta %+5.3f, phi %+5.3f, %d seed layers, %d hi-res layers, %d layers, front detid %d\n",
                                cell.icluster(),
                                hits3ds.eta(cell.hits().front()), hits3d.phi(cell.hits().front()),
                                cluster.nseedlayers(), cluster.nlayers(), cluster.nmorelayers(),
                                hits3ds.hit(cell.hits().front())->geographicalId().rawId());
                    }
                    std::vector<const HTCluster *> clustersSorted;
                    for (auto & cluster : map3d.clusters()) {
                        if (cluster.nmorelayers() >= layerMoreCut_) clustersSorted.push_back(&cluster);
                    } 
                    std::sort(clustersSorted.begin(), clustersSorted.end(),
                              [&map3d](const HTCluster *c1, const HTCluster *c2) -> bool 
                              { return (c1->nseedlayers() == c2->nseedlayers() ? c1->nlayers() > c2->nlayers() : c1->nseedlayers() > c2->nseedlayers()); });
                    for (const HTCluster *cluster : clustersSorted) {
                        tcBuilder_.run(*cluster, layerSeedCut3d_, layerMoreCut_, *out, &*outCl);
                    }
                    if (debugger_) if (DEBUG>=1) HTDebugger::dumpHTClusters(prefix+"c3dp", map3d, hits3ds, layerSeedCut3d_, layerMoreCut_);
                }

                // clear the map
                for (unsigned int i = 0; i < hits3ds.size(); ++i) {
                    map3d.clear(hits3ds.ieta(i), hits3ds.iphi(i));
                }
                map3d.clearClusters();
                for (unsigned int i = 0; i < hits2ds.size(); ++i) {
                    map2d.clear(hits2ds.ieta(i), hits2ds.iphi(i));
                }
                map2d.clearClusters();
            }
        }
    }

    tcBuilder_.done(*outTC);

    if (DEBUG>=0.1) {
        HTDebugger::registerClustersAsSeeds(2, *out);
        HTDebugger::registerClustersAsSeeds(1, *outCl);
        HTDebugger::registerTrackCandidates(*outTC);
        float ptMin = ptSteps_.back() == 0 ? 0.7 :  0.7*std::abs(ptSteps_.back());
        HTDebugger::printBackAssociation(*tracks, *outTC, *geometry_, *bfield_, ptMin);
    }
    
    iEvent.put(outTC);
    iEvent.put(out);
    iEvent.put(outCl, "clusters");
}

void 
TestHT::getPixelHits3D(const edm::Event & iEvent,  const reco::BeamSpot &bspot, HTHits3D & hits3d)  const 
{
   edm::Handle<edmNew::DetSetVector<SiPixelRecHit> > pixelHits;
   iEvent.getByToken(pixelHits_, pixelHits);
   GlobalPoint bspotPosition(bspot.position().x(), bspot.position().y(), 0);

   edm::Handle<PixelMask> pixelHitMask;
   if (hasMasks_) iEvent.getByToken(pixelHitMask_, pixelHitMask);

   for (auto ds : *pixelHits) {
        if (ds.empty()) continue;
        uint32_t id = ds.detId();
        //if (DetId(id).subdetId() != PixelSubdetector::PixelBarrel) continue;
        const GeomDet* geomDet = geometry_->idToDet(DetId(id));
	unsigned int layer = tTopo_->layer(DetId(id));
        unsigned int layermask = 1 << (layer-1);
        for (auto const &hit : ds) {
            if (hasMasks_ && pixelHitMask->mask(hit.cluster().key())) continue;
            GlobalVector gp = geomDet->surface().toGlobal( hit.localPosition() ) - bspotPosition;
            hits3d.push_back(&hit, gp.perp(), gp.phi(), gp.z(), layermask);
        }
   } 
}
void 
TestHT::getStripHits3D(const edm::Event & iEvent,  const reco::BeamSpot &bspot, HTHits3D & hits3d) const 
{
   edm::Handle<edmNew::DetSetVector<SiStripMatchedRecHit2D> > stripHits2D;
   iEvent.getByToken(stripHits2D_, stripHits2D);
   GlobalPoint bspotPosition(bspot.position().x(), bspot.position().y(), 0);

   edm::Handle<StripMask> stripHitMask;
   if (hasMasks_) iEvent.getByToken(stripHitMask_, stripHitMask);

   for (auto ds : *stripHits2D) {
        if (ds.empty()) continue;
        DetId id(ds.detId());
        const GeomDet* geomDet = geometry_->idToDet(id);
	unsigned int layer = tTopo_->layer(id) + 3;
        if (id.subdetId() == StripSubdetector::TOB || id.subdetId() == StripSubdetector::TEC) layer += 4;
        unsigned int layermask = 1 << (layer-1);
        for (auto const &hit : ds) {
            if (hasMasks_ && (stripHitMask->mask(hit.monoClusterRef().index()) || stripHitMask->mask(hit.stereoClusterRef().index()))) continue;
            GlobalVector gp = geomDet->surface().toGlobal( hit.localPosition() ) - bspotPosition;
            hits3d.push_back(&hit, gp.perp(), gp.phi(), gp.z(), layermask);
        }
   }
}
void 
TestHT::getStripHits2D(const edm::Event & iEvent,  const reco::BeamSpot &bspot, HTHits3D & hits3d) const 
{
   edm::Handle<edmNew::DetSetVector<SiStripRecHit2D> > stripHits;
   iEvent.getByToken(stripHits_, stripHits);
   GlobalPoint bspotPosition(bspot.position().x(), bspot.position().y(), 0);

   edm::Handle<StripMask> stripHitMask;
   if (hasMasks_) iEvent.getByToken(stripHitMask_, stripHitMask);

   for (auto ds : *stripHits) {
        if (ds.empty()) continue;
        DetId id(ds.detId());
        const GeomDet* geomDet = geometry_->idToDet(id);
	unsigned int layer = tTopo_->layer(id) + 3;
        if (id.subdetId() == StripSubdetector::TOB || id.subdetId() == StripSubdetector::TEC) layer += 4;
        unsigned int layermask = 1 << (layer-1);
        for (auto const &hit : ds) {
            if (hasMasks_ && stripHitMask->mask(hit.cluster().key())) continue;
            GlobalVector gp = geomDet->surface().toGlobal( hit.localPosition() ) - bspotPosition;
            hits3d.push_back(&hit, gp.perp(), gp.phi(), gp.z(), layermask);
        }
   }
}



#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TestHT);
