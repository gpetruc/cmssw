
//

/**
  \class    HTBinningStudies HTBinningStudies.h "RecoTracker/HTPattern/interface/HTBinningStudies.h"
            
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
#include "DataFormats/Math/interface/deltaPhi.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include <TTree.h>

class HTBinningStudies : public edm::EDProducer {
    public:
      explicit HTBinningStudies(const edm::ParameterSet & iConfig);
      virtual ~HTBinningStudies() { }

      virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup) override;

     
    private:
      /// Labels for input collections
      edm::EDGetTokenT<std::vector<reco::Track> > tracks_;

      // Configurables
      StringCutObjectSelector<reco::Track> trackSelection_;

      // EventSetup stuff
      edm::ESHandle<TrackerGeometry> geometry_;
      edm::ESHandle<TrackerTopology> tTopo_;
      edm::ESHandle<MagneticField> bfield_;

      TTree *tree_;
      Float_t  trackPt_, trackEta_;
      Float_t  deta_, dphi_; Int_t layer_;
};

HTBinningStudies::HTBinningStudies(const edm::ParameterSet & iConfig) :
    tracks_(consumes<std::vector<reco::Track> >(iConfig.getParameter<edm::InputTag>("tracks"))),
    trackSelection_(iConfig.getParameter<std::string>("trackSelection"))
{
    edm::Service<TFileService> fs;
    tree_ = fs->make<TTree>("t","t");
    tree_->Branch("pt",&trackPt_,"pt/F");
    tree_->Branch("eta",&trackEta_,"eta/F");
    tree_->Branch("deta",&deta_,"deta/F");
    tree_->Branch("dphi",&dphi_,"dphi/F");
    tree_->Branch("layer",&layer_,"layer/I");
}

void 
HTBinningStudies::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {
    //Get the geometry
    iSetup.get<TrackerDigiGeometryRecord>().get(geometry_);

    //Retrieve tracker topology from geometry
    iSetup.get<IdealGeometryRecord>().get(tTopo_);

    //Retrieve magnetic field
    iSetup.get<IdealMagneticFieldRecord>().get(bfield_);
    double bfield = bfield_->inTesla(GlobalPoint(0.,0.,0.)).z();

    edm::Handle<std::vector<reco::Track> > tracks;
    iEvent.getByToken(tracks_, tracks);
    for (const reco::Track & tk : *tracks) {
        if (!trackSelection_(tk)) continue;
        trackPt_ = tk.pt(); trackEta_ = tk.eta();
        GlobalPoint p0(tk.vx(), tk.vy(), tk.vz());
        float alpha = 0.5 * 0.003 * bfield / trackPt_;

        for (trackingRecHit_iterator ihit = tk.recHitsBegin(), ehit = tk.recHitsEnd(); ihit != ehit; ++ihit) {
            const TrackingRecHit *hit = &**ihit;
            if (!hit->isValid()) continue;
            DetId id = hit->geographicalId(); 
            const GeomDet* geomDet = geometry_->idToDet(id);
            if (geomDet == 0) continue;
            layer_ = tTopo_->layer(id);
            if (id.subdetId() > 2) layer_ += 3; // strip
            if (id.subdetId() == StripSubdetector::TOB || id.subdetId() == StripSubdetector::TEC) layer_ += 4;
            GlobalVector vec = geomDet->surface().toGlobal( hit->localPosition() ) - p0;
            deta_ = fabs(vec.eta() - tk.eta());
            dphi_ = fabs(deltaPhi(vec.phi()-alpha*tk.phi(), tk.phi()));
            printf("hit on layer %d: dphi = %7.5f\n", layer_, dphi_);
            tree_->Fill();
        }
    }
}

class HTBinningEfficiencyStudy : public edm::EDProducer {
    public:
      explicit HTBinningEfficiencyStudy(const edm::ParameterSet & iConfig);
      virtual ~HTBinningEfficiencyStudy() { }

      virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup) override;

     
    private:
      /// Labels for input collections
      edm::EDGetTokenT<std::vector<reco::Track> > tracks_;

      // Configurables
      StringCutObjectSelector<reco::Track> trackSelection_;

      // Vertices and beam spot
      bool useVertices_; 
      edm::EDGetTokenT<std::vector<reco::Vertex> > vertices_;
      edm::EDGetTokenT<reco::BeamSpot> beamSpot_;

      // EventSetup stuff
      edm::ESHandle<TrackerGeometry> geometry_;
      edm::ESHandle<TrackerTopology> tTopo_;
      edm::ESHandle<MagneticField>   bfield_;

      // Configurables
      Int_t  etabins2d_, etabins3d_, phibins2d_, phibins3d_;
      std::vector<double> ptSteps_,ptEdges_;

      TTree *tree_;
      Float_t  trackPt_, trackEta_, trackPhi_;
      Float_t  trackDz_, trackDxy_;
      Int_t    trackHits_, trackPixelHits_, trackLayers_, track3dLayers_;
      Int_t    trackAlgo_;
      Int_t    trackHP_;
      Int_t    seedLayers_, hiresLayers_, lowresSeedLayers_, lowresLayers_, moreLayers_;
      Float_t  ptStep_;
      Float_t  minDEta_, minDPhi_;
};

HTBinningEfficiencyStudy::HTBinningEfficiencyStudy(const edm::ParameterSet & iConfig) :
    tracks_(consumes<std::vector<reco::Track> >(iConfig.getParameter<edm::InputTag>("tracks"))),
    trackSelection_(iConfig.getParameter<std::string>("trackSelection")),
    useVertices_(iConfig.getParameter<bool>("useVertices")),
    vertices_(consumes<std::vector<reco::Vertex> >(iConfig.getParameter<edm::InputTag>("vertices"))),
    beamSpot_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
    etabins2d_(iConfig.getParameter<uint32_t>("etabins2d")),
    etabins3d_(iConfig.getParameter<uint32_t>("etabins3d")),
    phibins2d_(iConfig.getParameter<uint32_t>("phibins2d")),
    phibins3d_(iConfig.getParameter<uint32_t>("phibins3d")),
    ptSteps_(iConfig.getParameter<std::vector<double> >("ptSteps"))
{
    edm::Service<TFileService> fs;
    tree_ = fs->make<TTree>("t","t");
    tree_->Branch("pt",&trackPt_,"pt/F");
    tree_->Branch("eta",&trackEta_,"eta/F");
    tree_->Branch("phi",&trackPhi_,"phi/F");
    tree_->Branch("dxy",&trackDxy_,"dxy/F");
    tree_->Branch("dz",&trackDz_,"dz/F");
    tree_->Branch("hits",&trackHits_,"hits/F");
    tree_->Branch("pixelHits",&trackPixelHits_,"pixelHits/F");
    tree_->Branch("layers",&trackLayers_,"layers/F");
    tree_->Branch("layers3d",&track3dLayers_,"layers3dHits/F");
    tree_->Branch("algo",&trackAlgo_,"algo/F");
    tree_->Branch("hp",&trackHP_,"hp/F");
    tree_->Branch("ptStep",&ptStep_,"ptStep/F");
    tree_->Branch("minDEta",&minDEta_,"minDEta/F");
    tree_->Branch("minDPhi",&minDPhi_,"minDPhi/F");
    tree_->Branch("seedLayers",&seedLayers_,"seedLayers/I");
    tree_->Branch("hiresLayers",&hiresLayers_,"hiresLayers/I");
    tree_->Branch("lowresLayers",&lowresLayers_,"lowresLayers/I");
    tree_->Branch("lowresSeedLayers",&lowresSeedLayers_,"lowresSeedLayers/I");
    tree_->Branch("moreLayers",&moreLayers_,"moreLayers/I");
}

void 
HTBinningEfficiencyStudy::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {
    using std::max; using std::abs;

    //Get the geometry
    iSetup.get<TrackerDigiGeometryRecord>().get(geometry_);

    //Retrieve tracker topology from geometry
    iSetup.get<IdealGeometryRecord>().get(tTopo_);

    //Retrieve magnetic field
    iSetup.get<IdealMagneticFieldRecord>().get(bfield_);
    double bfield = bfield_->inTesla(GlobalPoint(0.,0.,0.)).z();

    edm::Handle<std::vector<reco::Track> > tracks;
    iEvent.getByToken(tracks_, tracks);

    edm::Handle<reco::BeamSpot> beamSpot;
    iEvent.getByToken(beamSpot_, beamSpot);

    edm::Handle<std::vector<reco::Vertex> > vertices;
    iEvent.getByToken(vertices_, vertices);

    HTHitMap map3d(etabins3d_,phibins3d_);
    HTHitMap map2d(etabins2d_,phibins2d_);
    int etashift = std::round(log(double(etabins3d_)/etabins2d_)/std::log(2));
    int phishift = std::round(log(double(phibins3d_)/phibins2d_)/std::log(2));

    for (const reco::Track & tk : *tracks) {
        if (!trackSelection_(tk)) continue;
        trackPt_ = tk.pt(); trackEta_ = tk.eta(); trackPhi_ = tk.phi();
        trackDz_  = abs(tk.dz(vertices->front().position()));
        trackDxy_ = abs(tk.dz(vertices->front().position()));
        trackHits_ = tk.hitPattern().numberOfValidHits();
        trackPixelHits_ = tk.hitPattern().numberOfValidPixelHits();
        trackLayers_    = tk.hitPattern().trackerLayersWithMeasurement();
        track3dLayers_  = tk.hitPattern().pixelLayersWithMeasurement() +  tk.hitPattern().numberOfValidStripLayersWithMonoAndStereo();
        trackAlgo_      = tk.algo();
        trackHP_        = tk.quality(reco::Track::highPurity);
        GlobalPoint p0;
        if (useVertices_) {
            p0 = GlobalPoint(beamSpot->position().x(), beamSpot->position().y(), vertices->front().z());
        } else {
            p0 = GlobalPoint(tk.vx(), tk.vy(), tk.vz());
        }

        HTHits3D hits3d(p0.x(), p0.y()); 
        HTHits3D hits2d(p0.x(), p0.y());
        for (trackingRecHit_iterator ihit = tk.recHitsBegin(), ehit = tk.recHitsEnd(); ihit != ehit; ++ihit) {
            const TrackingRecHit *hit = &**ihit;
            if (!hit->isValid()) continue;
            DetId id = hit->geographicalId(); 
            const GeomDet* geomDet = geometry_->idToDet(id);
            if (geomDet == 0) continue;
            int layer = tTopo_->layer(id);
            if (id.subdetId() > 2) layer += 3; // strip
            if (id.subdetId() == StripSubdetector::TOB || id.subdetId() == StripSubdetector::TEC) layer += 4;
            unsigned int layermask = 1 << (layer-1);
            GlobalVector vec = geomDet->surface().toGlobal( hit->localPosition() ) - p0;
            bool is3d = false;
            if (id.subdetId() <= 2) is3d = true;
            else if (id.subdetId() == StripSubdetector::TIB) is3d = tTopo_->tibIsDoubleSide(id);
            else if (id.subdetId() == StripSubdetector::TID) is3d = tTopo_->tidIsDoubleSide(id);
            else if (id.subdetId() == StripSubdetector::TOB) is3d = tTopo_->tobIsDoubleSide(id);
            else if (id.subdetId() == StripSubdetector::TEC) is3d = tTopo_->tecIsDoubleSide(id);
            if (is3d) {
                hits3d.push_back(hit, vec.perp(), vec.phi(), vec.z(), layermask);
            } 
            hits2d.push_back(hit, vec.perp(), vec.phi(), vec.z(), layermask);
        }

        unsigned int seedLayersMax = 0, hiresLayersMax = 0, lowresSeedLayersMax = 0, lowresLayersMax = 0, moreLayersMax = 0;
        for (unsigned int i = 0, n = ptSteps_.size(); i < n; ++i) {
            ptStep_ = ptSteps_[i];
            for (int ptsign = +1; ptsign > -2; ptsign -= 2) {
                if (ptStep_ == 0 && ptsign < 0) continue;

                float alpha = ptStep_ ? 0.5 * 0.003 * bfield / ptStep_ * ptsign : 0;

                HTHitsSpher hits3ds(hits3d, p0.z(), etabins3d_);
                hits3ds.filliphi(alpha, phibins3d_);
                HTHitsSpher hits2ds(hits2d, p0.z(), etabins2d_);
                hits2ds.filliphi(alpha, phibins2d_);

                // fill the maps
                for (unsigned int i = 0; i < hits3ds.size(); ++i) {
                    map3d.addHit(hits3ds.ieta(i), hits3ds.iphi(i), i, hits3ds.layermask(i), 2);
                }
                for (unsigned int i = 0; i < hits2ds.size(); ++i) {
                    map2d.addHit(hits2ds.ieta(i), hits2ds.iphi(i), i, hits2ds.layermask(i), 2);
                }

                // clusterize
                map3d.clusterize();
                map2d.clusterize();

                // count what found
                seedLayers_ = 0; hiresLayers_ = 0; lowresSeedLayers_ = 0; lowresLayers_ = 0; moreLayers_ = 0;

                // High res
                for (auto & cluster : map3d.clusters()) {
                    seedLayers_  = max<int>(seedLayers_,  cluster.nseedlayers());
                    hiresLayers_ = max<int>(hiresLayers_, cluster.nlayers());
                }
                // Low res
                for (auto & cluster : map2d.clusters()) {
                    lowresSeedLayers_ = max<int>(lowresSeedLayers_,  cluster.nseedlayers());
                    lowresLayers_     = max<int>(lowresLayers_,      cluster.nlayers());
                }
                // Mixed
                for (auto & cluster : map3d.clusters()) {
                    cluster.addMoreCell(map2d.get(cluster.ieta() >> etashift, cluster.iphi() >> phishift));
                    moreLayers_ = max<int>(moreLayers_, cluster.nmorelayers());
                }
                tree_->Fill();

                // store the best among the pt steps
                seedLayersMax       = max<int>(seedLayers_,seedLayersMax);
                hiresLayersMax      = max<int>(hiresLayers_,hiresLayersMax);
                lowresSeedLayersMax = max<int>(lowresSeedLayers_,lowresSeedLayersMax);
                lowresLayersMax     = max<int>(lowresLayers_,lowresLayersMax);
                moreLayersMax       = max<int>(moreLayers_,moreLayersMax);

                // clear the map
                for (unsigned int i = 0; i < hits3ds.size(); ++i) {
                    map3d.clear(hits3ds.ieta(i), hits3ds.iphi(i));
                }
                map3d.clearClusters();
                for (unsigned int i = 0; i < hits2ds.size(); ++i) {
                    map2d.clear(hits2ds.ieta(i), hits2ds.iphi(i));
                }
                map2d.clearClusters();
            } // ptSign
        } // ptStep
        ptStep_ = -1;
        seedLayers_ = seedLayersMax; hiresLayers_ = hiresLayersMax; lowresSeedLayers_ = lowresSeedLayersMax; lowresLayers_ = lowresLayersMax; moreLayers_ = moreLayersMax;
        tree_->Fill();
    }
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HTBinningStudies);
DEFINE_FWK_MODULE(HTBinningEfficiencyStudy);
