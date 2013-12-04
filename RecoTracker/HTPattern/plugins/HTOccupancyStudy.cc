//
//

/**
  \class    HTOccupancyStudy HTOccupancyStudy.h "RecoTracker/HTPattern/interface/HTOccupancyStudy.h"
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
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include <TTree.h>


class HTOccupancyStudy : public edm::EDProducer {
    public:
      explicit HTOccupancyStudy(const edm::ParameterSet & iConfig);
      virtual ~HTOccupancyStudy() { }

      virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup) override;

      void getPixelHits3D(const edm::Event & iEvent, const reco::BeamSpot &bsp, HTHits3D &hits3d) const ;
      void getStripHits3D(const edm::Event & iEvent, const reco::BeamSpot &bsp, HTHits3D &hits3d) const ;
      void getStripHits2D(const edm::Event & iEvent, const reco::BeamSpot &bsp, HTHits3D &hits3d) const ;
      
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

      // EventSetup stuff
      edm::ESHandle<TrackerGeometry> geometry_;
      edm::ESHandle<TrackerTopology> tTopo_;
      edm::ESHandle<MagneticField> bfield_;

      // Debugging
      TTree *tree_;
      Float_t ptStep_;
      Int_t firstvertex_;
      Int_t seedLayers_, hiresLayers_, moreLayers_;
};

HTOccupancyStudy::HTOccupancyStudy(const edm::ParameterSet & iConfig) :
    pixelHits_(consumes<edmNew::DetSetVector<SiPixelRecHit> >(iConfig.getParameter<edm::InputTag>("pixelHits"))),
    stripHits_(consumes<edmNew::DetSetVector<SiStripRecHit2D> >(iConfig.getParameter<edm::InputTag>("stripHits"))),
    stripHits2D_(consumes<edmNew::DetSetVector<SiStripMatchedRecHit2D> >(iConfig.getParameter<edm::InputTag>("stripHits2D"))),
    hasMasks_(iConfig.existsAs<edm::InputTag>("clustersToSkip")),
    vertices_(consumes<std::vector<reco::Vertex> >(iConfig.getParameter<edm::InputTag>("vertices"))),
    beamSpot_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
    tracks_(consumes<std::vector<reco::Track> >(iConfig.getParameter<edm::InputTag>("tracks"))),
    seed2d_(iConfig.getParameter<bool>("seed2d")),
    seed3d_(iConfig.getParameter<bool>("seed3d")),
    seedMixed_(iConfig.getParameter<bool>("seedMixed")),
    etabins2d_(iConfig.getParameter<uint32_t>("etabins2d")),
    etabins3d_(iConfig.getParameter<uint32_t>("etabins3d")),
    phibins2d_(iConfig.getParameter<uint32_t>("phibins2d")),
    phibins3d_(iConfig.getParameter<uint32_t>("phibins3d")),
    ptSteps_(iConfig.getParameter<std::vector<double> >("ptSteps")),
    vertexSelection_(iConfig.getParameter<std::string>("vertexSelection"))
{
    if (hasMasks_) {
        pixelHitMask_ = consumes<PixelMask>(iConfig.getParameter<edm::InputTag>("clustersToSkip"));
        stripHitMask_ = consumes<StripMask>(iConfig.getParameter<edm::InputTag>("clustersToSkip"));
    }
    edm::Service<TFileService> fs;
    tree_ = fs->make<TTree>("t","t");
    tree_->Branch("ptStep",&ptStep_,"ptStep/F");
    tree_->Branch("seedLayers",&seedLayers_,"seedLayers/I");
    tree_->Branch("hiresLayers",&hiresLayers_,"hiresLayers/I");
    tree_->Branch("moreLayers",&moreLayers_,"moreLayers/I");
}

void 
HTOccupancyStudy::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {
    using namespace edm;
    using namespace std;

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

    bool do3d = seed3d_ || seedMixed_;
    bool do2d = seed2d_ || seedMixed_;
    HTHits3D hits3d(beamSpot->position().x(), beamSpot->position().y()); // make pixels
    if (do3d) getPixelHits3D(iEvent, *beamSpot, hits3d);

    HTHits3D hits2d(beamSpot->position().x(), beamSpot->position().y());
    if (do2d) hits2d = hits3d; // start as copy
    if (do2d) getStripHits2D(iEvent, *beamSpot, hits2d);
    if (do3d) getStripHits3D(iEvent, *beamSpot, hits3d);

    printf("Number of hits: %d 3d, %d 2d\n", hits3d.size(), hits2d.size());
    // Now we just do one step
    HTHitMap map3d(etabins3d_,phibins3d_);
    HTHitMap map2d(etabins2d_,phibins2d_);
    int etashift = std::round(log(double(etabins3d_)/etabins2d_)/std::log(2));
    int phishift = std::round(log(double(phibins3d_)/phibins2d_)/std::log(2));
    for (unsigned int iptStep = 0, nptSteps = ptSteps_.size(); iptStep < nptSteps; ++iptStep) {
        double ptStep = ptSteps_[iptStep]; ptStep_ = ptStep;
        if (ptStep != 0 && bfield == 0) continue;

        for (int ptsign = +1; ptsign > -2; ptsign -= 2) {
            if (ptStep == 0 && ptsign < 0) continue;

            bool firstvertex_ = true;
            for(const reco::Vertex &vtx : *vertices) {
                // always process first vertex, but apply selection to others
                if (!firstvertex_ && !vertexSelection_(vtx)) continue; else firstvertex_ = false;

                float alpha = ptStep ? 0.5 * 0.003 * bfield / ptStep * ptsign : 0;
                HTHitsSpher hits3ds(hits3d, vtx.z(), etabins3d_);
                hits3ds.filliphi(alpha, phibins3d_);
                HTHitsSpher hits2ds(hits2d, vtx.z(), etabins2d_);
                hits2ds.filliphi(alpha, phibins2d_);

                for (unsigned int i = 0; i < hits3ds.size(); ++i) {
                    map3d.addHit(hits3ds.ieta(i), hits3ds.iphi(i), i, hits3ds.layermask(i), 3);
                }
                for (unsigned int i = 0; i < hits2ds.size(); ++i) {
                    map2d.addHit(hits2ds.ieta(i), hits2ds.iphi(i), i, hits2ds.layermask(i), 3);
                }
                // clusterize
                if (do3d) map3d.clusterize();

                if (seed3d_) {
                    for (auto & cluster : map3d.clusters()) {
                        seedLayers_  = cluster.nseedlayers();
                        hiresLayers_ = cluster.nlayers();
                        tree_->Fill();
                    }
                }

                if (seed2d_) {
                    map2d.clusterize();
                    for (auto & cluster : map3d.clusters()) {
                        seedLayers_  = cluster.nseedlayers();
                        hiresLayers_ = cluster.nlayers();
                        tree_->Fill();
                    }
                }

                if (seedMixed_) {
                    for (auto & cluster : map3d.clusters()) {
                        cluster.addMoreCell(map2d.get(cluster.ieta() >> etashift, cluster.iphi() >> phishift));
                        seedLayers_  = cluster.nseedlayers();
                        hiresLayers_ = cluster.nlayers();
                        moreLayers_ = cluster.nmorelayers();
                        tree_->Fill();
                    }
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


}

void 
HTOccupancyStudy::getPixelHits3D(const edm::Event & iEvent,  const reco::BeamSpot &bspot, HTHits3D & hits3d)  const 
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
HTOccupancyStudy::getStripHits3D(const edm::Event & iEvent,  const reco::BeamSpot &bspot, HTHits3D & hits3d) const 
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
HTOccupancyStudy::getStripHits2D(const edm::Event & iEvent,  const reco::BeamSpot &bspot, HTHits3D & hits3d) const 
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
DEFINE_FWK_MODULE(HTOccupancyStudy);
