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
#include "RecoTracker/HTPattern/interface/TrackCandidateBuilderFromCluster.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"




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
      edm::EDGetTokenT<MeasurementTrackerEvent> mtEvent_;


      // Configurables
      uint32_t etabins2d_, etabins3d_, phibins2d_, phibins3d_;
      uint32_t layerSeedCut2d_, layerSeedCut3d_;
      uint32_t layerCut2d_, layerCut3d_;
      uint32_t layerMoreCut_;
      std::vector<double> ptSteps_;

      // Helpers
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
    vertices_(consumes<std::vector<reco::Vertex> >(iConfig.getParameter<edm::InputTag>("vertices"))),
    beamSpot_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
    tracks_(consumes<std::vector<reco::Track> >(iConfig.getParameter<edm::InputTag>("tracks"))),
    mtEvent_(consumes<MeasurementTrackerEvent>(iConfig.getParameter<edm::InputTag>("measurementTrackerEvent"))),
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
    tcBuilder_(iConfig.getParameter<edm::ParameterSet>("seedBuilderConfig")),
    debugger_(iConfig.getUntrackedParameter<bool>("debugger",false))
{
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
    std::cout << "Magnetic field: " << bfield << std::endl;

    edm::Handle<reco::BeamSpot> beamSpot;
    iEvent.getByToken(beamSpot_, beamSpot);

    Handle<vector<reco::Vertex> > vertices;
    iEvent.getByToken(vertices_, vertices);

    Handle<vector<reco::Track> > tracks;
    iEvent.getByToken(tracks_, tracks);
    if (debugger_) HTDebugger::dumpTracks(prefix+"tk", *tracks);

    edm::Handle<MeasurementTrackerEvent> mtEvent;
    iEvent.getByToken(mtEvent_, mtEvent);

    // initialize the builder
    tcBuilder_.init(iSetup, *mtEvent, nullptr);

    std::auto_ptr<TrajectorySeedCollection> out(new TrajectorySeedCollection());
    std::auto_ptr<TrajectorySeedCollection> outCl(new TrajectorySeedCollection());
    std::auto_ptr<TrackCandidateCollection> outTC(new TrackCandidateCollection());

    HTHits3D hits3d(beamSpot->position().x(), beamSpot->position().y()); // make pixels
    getPixelHits3D(iEvent, *beamSpot, hits3d);
    HTHits3D hits2d(hits3d); // then copy, and make strips
    getStripHits2D(iEvent, *beamSpot, hits2d);
    getStripHits3D(iEvent, *beamSpot, hits3d);
    if (debugger_) HTDebugger::dumpHTHits3D(prefix+"3d", hits3d, vertices->empty() ? 0.0 : vertices->front().z());
    if (debugger_) HTDebugger::dumpHTHits3D(prefix+"2d", hits2d, vertices->empty() ? 0.0 : vertices->front().z());


    // Now we just do one step
    HTHitMap map3d(etabins3d_,phibins3d_);
    HTHitMap map2d(etabins2d_,phibins2d_);
    std::vector<bool> mask3d(hits3d.size(), false);
    std::vector<bool> mask2d(hits2d.size(), false);
    int etashift = std::round(log(double(etabins3d_)/etabins2d_)/std::log(2));
    int phishift = std::round(log(double(phibins3d_)/phibins2d_)/std::log(2));
    for(const reco::Vertex &vtx : *vertices) {

        for (double ptStep : ptSteps_) {
            if (ptStep != 0 && bfield == 0) continue;

            for (int ptsign = +1; ptsign > -2; ptsign -= 2) {
                if (ptStep == 0 && ptsign < 0) continue;

                sprintf(evid, "r%d_l%d_e%d_pT%c%03d_", iEvent.id().run(), iEvent.id().luminosityBlock(), iEvent.id().event(), ptsign > 0 ? 'p' : 'm', int(ptStep*100));
                std::string prefix(evid);

                float alpha = ptStep ? 0.5 * 0.003 * bfield / ptStep * ptsign : 0;

                printf("\n========== HT iteration with pT = %+.2f, alpha = %+.5f ==========\n", ptStep*ptsign, alpha);

                HTHitsSpher hits3ds(hits3d, vtx.z(), etabins3d_);
                hits3ds.filliphi(alpha, phibins3d_);
                HTHitsSpher hits2ds(hits2d, vtx.z(), etabins2d_);
                hits2ds.filliphi(alpha, phibins2d_);

                tcBuilder_.setHits(hits3ds, hits2ds, map3d, map2d, mask3d, mask2d, etashift, phishift);

                if (debugger_) HTDebugger::dumpHTHitsSpher(prefix+"spher3d", hits3ds);
                if (debugger_) HTDebugger::dumpHTHitsSpher(prefix+"spher2d", hits2ds);
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
                map3d.clusterize();
                map2d.clusterize();

                // analzye
                if (debugger_) HTDebugger::dumpHTHitMap(prefix+"map3d", map3d);
                if (debugger_) HTDebugger::dumpHTHitMap(prefix+"map2d", map2d);
                if (debugger_) HTDebugger::dumpHTClusters(prefix+"c3d", map3d, hits3ds, layerCut3d_);
                if (debugger_) HTDebugger::dumpHTClusters(prefix+"c2d", map2d, hits2ds, layerCut2d_);

                // Recover more layers from 2D map
                for (auto & cluster : map3d.clusters()) {
                    cluster.addMoreCell(map2d.get(cluster.ieta() >> etashift, cluster.iphi() >> phishift));
                    const HTCell &cell = map3d.get(cluster.ieta(),cluster.iphi());
                    printf("cluster #%03d: eta %+5.3f, phi %+5.3f, %d seed layers, %d hi-res layers, %d layers\n",
                        cell.icluster(),
                        hits3ds.eta(cell.hits().front()), hits3d.phi(cell.hits().front()),
                        cluster.nseedlayers(), cluster.nlayers(), cluster.nmorelayers());
                }
                for (auto & cluster : map3d.clusters()) {
                    if (cluster.nmorelayers() >= layerMoreCut_) tcBuilder_.run(cluster, *outTC, *out, &*outCl);
                }
                if (debugger_) HTDebugger::dumpHTClusters(prefix+"c3dp", map3d, hits3ds, layerSeedCut3d_, layerMoreCut_);

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
        break;
    }

    tcBuilder_.done();
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

   for (auto ds : *pixelHits) {
        if (ds.empty()) continue;
        uint32_t id = ds.detId();
        //if (DetId(id).subdetId() != PixelSubdetector::PixelBarrel) continue;
        const GeomDet* geomDet = geometry_->idToDet(DetId(id));
	unsigned int layer = tTopo_->layer(DetId(id));
        unsigned int layermask = 1 << (layer-1);
        for (auto const &hit : ds) {
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

   for (auto ds : *stripHits2D) {
        if (ds.empty()) continue;
        uint32_t id = ds.detId();
        const GeomDet* geomDet = geometry_->idToDet(DetId(id));
	unsigned int layer = tTopo_->layer(DetId(id)) + 3;
        unsigned int layermask = 1 << (layer-1);
        for (auto const &hit : ds) {
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

   for (auto ds : *stripHits) {
        if (ds.empty()) continue;
        uint32_t id = ds.detId();
        const GeomDet* geomDet = geometry_->idToDet(DetId(id));
	unsigned int layer = tTopo_->layer(DetId(id)) + 3;
        unsigned int layermask = 1 << (layer-1);
        for (auto const &hit : ds) {
            GlobalVector gp = geomDet->surface().toGlobal( hit.localPosition() ) - bspotPosition;
            hits3d.push_back(&hit, gp.perp(), gp.phi(), gp.z(), layermask);
        }
   }
}



#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TestHT);
