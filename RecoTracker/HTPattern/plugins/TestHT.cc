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
#include "RecoTracker/HTPattern/plugins/HTDebugger.h"




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
      /// Labels for input collections
      edm::EDGetTokenT<edmNew::DetSetVector<SiPixelRecHit> > pixelHits_;
      edm::EDGetTokenT<edmNew::DetSetVector<SiStripRecHit2D> > stripHits_;
      edm::EDGetTokenT<edmNew::DetSetVector<SiStripMatchedRecHit2D> > stripHits2D_;
      edm::EDGetTokenT<std::vector<reco::Vertex> > vertices_;
      edm::EDGetTokenT<reco::BeamSpot> beamSpot_;
      edm::EDGetTokenT<std::vector<reco::Track> > tracks_;

      uint32_t etabins2d_, etabins3d_, phibins2d_, phibins3d_;
      uint32_t layerSeedCut2d_, layerSeedCut3d_;
      edm::ESHandle<TrackerGeometry> geometry_;
      edm::ESHandle<TrackerTopology> tTopo_;
};

TestHT::TestHT(const edm::ParameterSet & iConfig) :
    pixelHits_(consumes<edmNew::DetSetVector<SiPixelRecHit> >(iConfig.getParameter<edm::InputTag>("pixelHits"))),
    stripHits_(consumes<edmNew::DetSetVector<SiStripRecHit2D> >(iConfig.getParameter<edm::InputTag>("stripHits"))),
    stripHits2D_(consumes<edmNew::DetSetVector<SiStripMatchedRecHit2D> >(iConfig.getParameter<edm::InputTag>("stripHits2D"))),
    vertices_(consumes<std::vector<reco::Vertex> >(iConfig.getParameter<edm::InputTag>("vertices"))),
    beamSpot_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
    tracks_(consumes<std::vector<reco::Track> >(iConfig.getParameter<edm::InputTag>("tracks"))),
    etabins2d_(iConfig.getParameter<uint32_t>("etabins2d")),
    etabins3d_(iConfig.getParameter<uint32_t>("etabins3d")),
    phibins2d_(iConfig.getParameter<uint32_t>("phibins2d")),
    phibins3d_(iConfig.getParameter<uint32_t>("phibins3d")),
    layerSeedCut2d_(iConfig.getParameter<uint32_t>("layerSeedCut2d")),
    layerSeedCut3d_(iConfig.getParameter<uint32_t>("layerSeedCut3d"))
{
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

    edm::Handle<reco::BeamSpot> beamSpot;
    iEvent.getByToken(beamSpot_, beamSpot);

    //Retrieve tracker topology from geometry
    iSetup.get<IdealGeometryRecord>().get(tTopo_);

    Handle<vector<reco::Vertex> > vertices;
    iEvent.getByToken(vertices_, vertices);

    Handle<vector<reco::Track> > tracks;
    iEvent.getByToken(tracks_, tracks);
    HTDebugger::dumpTracks(prefix+"tk", *tracks);

    HTHits3D hits3d; // make pixels
    getPixelHits3D(iEvent, *beamSpot, hits3d);
    HTHits3D hits2d(hits3d); // then copy, and make strips
    getStripHits2D(iEvent, *beamSpot, hits2d);
    getStripHits3D(iEvent, *beamSpot, hits3d);
    HTDebugger::dumpHTHits3D(prefix+"3d", hits3d, vertices->empty() ? 0.0 : vertices->front().z());
    HTDebugger::dumpHTHits3D(prefix+"2d", hits2d, vertices->empty() ? 0.0 : vertices->front().z());


    // Now we just do one step
    HTHitMap map3d(etabins3d_,phibins3d_);
    HTHitMap map2d(etabins2d_,phibins2d_);
    std::vector<bool> mask(hits3d.size(), false);
    for(const reco::Vertex &vtx : *vertices) {
        HTHitsSpher hits3ds(hits3d, vtx.z(), etabins3d_);
        hits3ds.filliphi(0.f, phibins3d_);
        HTHitsSpher hits2ds(hits2d, vtx.z(), etabins2d_);
        hits2ds.filliphi(0.f, phibins2d_);
        HTDebugger::dumpHTHitsSpher(prefix+"spher3d", hits3ds);
        HTDebugger::dumpHTHitsSpher(prefix+"spher2d", hits2ds);
        // fill the maps
        for (unsigned int i = 0; i < hits3ds.size(); ++i) {
            if (mask[i]) continue;
            map3d.addHit(hits3ds.ieta(i), hits3ds.iphi(i), i, hits3ds.layermask(i), layerSeedCut3d_);
        }
        // fill the maps
        for (unsigned int i = 0; i < hits2ds.size(); ++i) {
            if (mask[i]) continue;
            map2d.addHit(hits2ds.ieta(i), hits2ds.iphi(i), i, hits2ds.layermask(i), layerSeedCut2d_);
        }
        // analzye
        HTDebugger::dumpHTHitMap(prefix+"map3d", map3d);
        HTDebugger::dumpHTHitMap(prefix+"map2d", map2d);
        HTDebugger::dumpHTClusters(prefix+"c3d", map3d, hits3ds);
        HTDebugger::dumpHTClusters(prefix+"c2d", map2d, hits2ds);
        // clear the map
        for (unsigned int i = 0; i < hits3ds.size(); ++i) {
             map3d.clear(hits3ds.ieta(i), hits3ds.iphi(i));
        }
        map3d.clearClusters();
        for (unsigned int i = 0; i < hits2ds.size(); ++i) {
             map2d.clear(hits2ds.ieta(i), hits2ds.iphi(i));
        }
        map2d.clearClusters();
        break;
    }

    std::auto_ptr<TrackCandidateCollection> out(new TrackCandidateCollection());
    iEvent.put(out);
}

void 
TestHT::getPixelHits3D(const edm::Event & iEvent,  const reco::BeamSpot &bspot, HTHits3D & hits3d)  const 
{
   edm::Handle<edmNew::DetSetVector<SiPixelRecHit> > pixelHits;
   iEvent.getByToken(pixelHits_, pixelHits);
   GlobalPoint bspotPosition(bspot.position().x(), bspot.position().y(), 0);

   for (auto ds : *pixelHits) {
        if (ds.empty()) continue;
        uint32_t id = ds.detId();
        //if (DetId(id).subdetId() != PixelSubdetector::PixelBarrel) continue;
        const GeomDet* geomDet = geometry_->idToDet(DetId(id));
	unsigned int layer = tTopo_->layer(DetId(id));
        unsigned int layermask = 1 << (layer-1);
        for (auto hit : ds) {
            GlobalVector gp = geomDet->surface().toGlobal( hit.localPosition() ) - bspotPosition;
            hits3d.push_back(hit, gp.perp(), gp.phi(), gp.z(), layermask);
        }
   } 
}
void 
TestHT::getStripHits3D(const edm::Event & iEvent,  const reco::BeamSpot &bspot, HTHits3D & hits3d) const 
{
   edm::Handle<edmNew::DetSetVector<SiStripMatchedRecHit2D> > stripHits2D;
   iEvent.getByToken(stripHits2D_, stripHits2D);
   GlobalPoint bspotPosition(bspot.position().x(), bspot.position().y(), 0);

   for (auto ds : *stripHits2D) {
        if (ds.empty()) continue;
        uint32_t id = ds.detId();
        const GeomDet* geomDet = geometry_->idToDet(DetId(id));
	unsigned int layer = tTopo_->layer(DetId(id)) + 3;
        unsigned int layermask = 1 << (layer-1);
        for (auto hit : ds) {
            GlobalVector gp = geomDet->surface().toGlobal( hit.localPosition() ) - bspotPosition;
            hits3d.push_back(hit, gp.perp(), gp.phi(), gp.z(), layermask);
        }
   }
}
void 
TestHT::getStripHits2D(const edm::Event & iEvent,  const reco::BeamSpot &bspot, HTHits3D & hits3d) const 
{
   edm::Handle<edmNew::DetSetVector<SiStripRecHit2D> > stripHits;
   iEvent.getByToken(stripHits_, stripHits);
   GlobalPoint bspotPosition(bspot.position().x(), bspot.position().y(), 0);

   for (auto ds : *stripHits) {
        if (ds.empty()) continue;
        uint32_t id = ds.detId();
        const GeomDet* geomDet = geometry_->idToDet(DetId(id));
	unsigned int layer = tTopo_->layer(DetId(id)) + 3;
        unsigned int layermask = 1 << (layer-1);
        for (auto hit : ds) {
            GlobalVector gp = geomDet->surface().toGlobal( hit.localPosition() ) - bspotPosition;
            hits3d.push_back(hit, gp.perp(), gp.phi(), gp.z(), layermask);
        }
   }
}



#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TestHT);
