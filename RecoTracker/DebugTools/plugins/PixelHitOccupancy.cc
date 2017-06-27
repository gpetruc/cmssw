
// system include files
#include <memory>
#include <cmath>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"

#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include <TTree.h>
//#include "CalibFormats/SiPixelObjects/interface/SiPixelQuality.h"
//#include "CalibTracker/Records/interface/SiPixelQualityRcd.h"


class PixelHitOccupancy : public edm::one::EDProducer<edm::one::SharedResources> {
    public:
        explicit PixelHitOccupancy(const edm::ParameterSet&);
        ~PixelHitOccupancy();

    private:
        virtual void produce(edm::Event&, const edm::EventSetup&);

        // ----------member data ---------------------------
        const edm::EDGetTokenT<std::vector<reco::Track>> tracks_;    
        const edm::EDGetTokenT<edmNew::DetSetVector<SiPixelCluster>> pixelClusterLabel_;

        edm::ESHandle<TrackerTopology> theTrkTopo;
        //edm::ESHandle<SiPixelQuality> thePixelQuality;
        edm::ESHandle<TrackerGeometry> theGeometry;

        // output
        TTree *tree_;
        int run_, lumi_;    
        int detid_, ontrack_, charge_;
        float module_z_, module_phi_, hit_z_, hit_phi_, hit_lx_, hit_ly_, hit_px_, hit_py_;

        const PixelGeomDetUnit*  doModule_(DetId detid) {
            const PixelGeomDetUnit* det = dynamic_cast<const PixelGeomDetUnit*>(theGeometry->idToDet(detid));
            if (det) {
                const auto & detgp = det->position();
                detid_ = detid.rawId();
                module_z_ = detgp.z();
                module_phi_ = detgp.phi();
            }
            return det;
        }
        void fill_(const SiPixelCluster & cluster, const PixelGeomDetUnit * det) {
            const PixelTopology& topol = det->specificTopology();
            LocalPoint  clustlp = topol.localPosition(MeasurementPoint(cluster.x(), cluster.y()));
            GlobalPoint clustgp = det->surface().toGlobal(clustlp);
            hit_px_ = cluster.x();
            hit_py_ = cluster.y();
            hit_lx_ = clustlp.x();
            hit_ly_ = clustlp.y();
            hit_z_ = clustgp.z();
            hit_phi_ = clustgp.phi();
            tree_->Fill();
        }
        
};

//
// constructors and destructor
//
PixelHitOccupancy::PixelHitOccupancy(const edm::ParameterSet& iConfig):
    tracks_(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("tracks"))),
    pixelClusterLabel_(consumes<edmNew::DetSetVector<SiPixelCluster>>(iConfig.getParameter<edm::InputTag>("pixelClusters")))
{
    usesResource("TFileService");
    edm::Service<TFileService> fs;
    tree_ = fs->make<TTree>("tree","tree");
    tree_->Branch("run",  &run_,  "run/i");
    tree_->Branch("lumi", &lumi_, "lumi/i");
    tree_->Branch("detid", &detid_, "detid/I");
    tree_->Branch("ontrack", &ontrack_, "ontrack/I");
    tree_->Branch("charge", &charge_, "charge/I");
    tree_->Branch("module_z", &module_z_, "module_z/F");
    tree_->Branch("module_phi", &module_phi_, "module_phi/F");
    tree_->Branch("hit_z", &hit_z_, "hit_z/F");
    tree_->Branch("hit_phi", &hit_phi_, "hit_phi/F");
    tree_->Branch("hit_lx", &hit_lx_, "hit_lx/F");
    tree_->Branch("hit_ly", &hit_ly_, "hit_ly/F");
    tree_->Branch("hit_px", &hit_px_, "hit_px/F");
    tree_->Branch("hit_py", &hit_py_, "hit_py/F");
}


PixelHitOccupancy::~PixelHitOccupancy()
{
}

void
PixelHitOccupancy::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;

    run_  = iEvent.id().run();
    lumi_ = iEvent.id().luminosityBlock();

    // read input
    Handle<std::vector<reco::Track>> tracks;
    iEvent.getByToken(tracks_, tracks);

    iSetup.get<TrackerTopologyRcd>().get(theTrkTopo);
    iSetup.get<TrackerDigiGeometryRecord>().get(theGeometry);


    edm::Handle<edmNew::DetSetVector<SiPixelCluster> > pixelC; 
    iEvent.getByToken(pixelClusterLabel_, pixelC); 

    ontrack_ = 1;
    for (const reco::Track & track : *tracks) {
        if (!track.hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 1)) continue;
        for (unsigned int i = 0, nhits = track.recHitsSize(); i < nhits; ++i) {
            const TrackingRecHit *hit = &* track.recHit(i);
            DetId id = hit->geographicalId();
            if (hit->isValid() && id.subdetId() == PixelSubdetector::SubDetector::PixelBarrel && theTrkTopo->layer(id) == 1) { 
                const auto * pixhit = dynamic_cast<const SiPixelRecHit*>(hit->hit()); if (!pixhit) throw cms::Exception("CorruptData", "Valid PXB1 hit that is not a SiPixelRecHit");
                const PixelGeomDetUnit*  det = doModule_(id); if (!det) throw cms::Exception("CorruptData", "Missing DetID in geometry??");
                auto  clustref = pixhit->cluster(); if (clustref.isNull()) throw cms::Exception("CorruptData", "Valid PXB1 SiPixelRecHit with null cluster ref");
                fill_(*clustref, det);
            }
        }
    }

    ontrack_ = 0;
    for (const auto & detset : *pixelC) {
        DetId id(detset.detId());
        if (id.subdetId() == PixelSubdetector::SubDetector::PixelBarrel && theTrkTopo->layer(id) == 1) {
            const PixelGeomDetUnit * det = doModule_(id); if (!det) throw cms::Exception("CorruptData", "Missing DetID in geometry??");
            for (const SiPixelCluster & cluster : detset) { 
                fill_(cluster, det);
            }
        }
    }
}


//define this as a plug-in
DEFINE_FWK_MODULE(PixelHitOccupancy);
