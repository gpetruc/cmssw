
// system include files
#include <memory>
#include <cmath>
#include <unordered_map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"

#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"

#include "RecoTracker/MeasurementDet/interface/MeasurementTrackerEvent.h"
#include "RecoTracker/MeasurementDet/src/TkMeasurementDetSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include <TTree.h>
//#include "CalibFormats/SiPixelObjects/interface/SiPixelQuality.h"
//#include "CalibTracker/Records/interface/SiPixelQualityRcd.h"


class PixelHitOccupancyByLumi : public edm::one::EDAnalyzer<edm::one::WatchLuminosityBlocks,edm::one::SharedResources> {
    public:
        explicit PixelHitOccupancyByLumi(const edm::ParameterSet&);
        ~PixelHitOccupancyByLumi();

    private:
        virtual void analyze(const edm::Event&, const edm::EventSetup&) ;
        void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) { firstEventInLS_ = true; } // waitNextLS_ = false; }
        void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) ;

        // ----------member data ---------------------------
        //const edm::EDGetTokenT<MeasurementTrackerEvent> tracker_;
        const edm::EDGetTokenT<edmNew::DetSetVector<SiPixelCluster>> pixelClusterLabel_;
        //unsigned int maxHitsPerROC_;

        edm::ESHandle<TrackerTopology> theTrkTopo;
        edm::ESHandle<TrackerGeometry> theGeometry;
        uint32_t theGeometryCacheID;

        std::unordered_map<unsigned int, unsigned int> detIdToIndexMap_;
        std::vector<unsigned int> detIds_;
        std::vector<std::array<unsigned int, 4>> tpmHits_;
        std::vector<std::array<unsigned int,16>> rocHits_;
        bool firstEventInLS_; //, waitNextLS_;

        TTree *tree_;
        int run_, lumi_;    
        int detid_, tpm_, roc_;        

};

//
// constructors and destructor
//
PixelHitOccupancyByLumi::PixelHitOccupancyByLumi(const edm::ParameterSet& iConfig):
    //tracker_(consumes<MeasurementTrackerEvent>(iConfig.getParameter<edm::InputTag>("tracker"))),
    pixelClusterLabel_(consumes<edmNew::DetSetVector<SiPixelCluster>>(iConfig.getParameter<edm::InputTag>("pixelClusters"))),
    //maxHitsPerROC_(iConfig.getParameter<uint32_t>("maxHitsPerROC")),
    theGeometryCacheID(0),
    firstEventInLS_(true)//, waitNextLS_(false)
{
    usesResource("TFileService");
    edm::Service<TFileService> fs;
    tree_ = fs->make<TTree>("tree","tree");
    tree_->Branch("run",  &run_,  "run/i");
    tree_->Branch("lumi", &lumi_, "lumi/i");
    tree_->Branch("detid", &detid_, "detid/I");
    tree_->Branch("tpm", &tpm_, "tpm/I");
    tree_->Branch("roc", &roc_, "roc/I");
}


PixelHitOccupancyByLumi::~PixelHitOccupancyByLumi()
{
}

void
PixelHitOccupancyByLumi::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace std;
    if (firstEventInLS_) {
        unsigned int geocache = iSetup.get<TrackerDigiGeometryRecord>().cacheIdentifier();
        if (geocache != theGeometryCacheID) {
            iSetup.get<TrackerTopologyRcd>().get(theTrkTopo);
            iSetup.get<TrackerDigiGeometryRecord>().get(theGeometry);
            detIdToIndexMap_.clear();
            detIds_.clear();
            for (const auto * det : theGeometry->detsPXB()) {
                DetId id = det->geographicalId();
                if (theTrkTopo->layer(id) == 1) {
                    detIds_.push_back(id());
                    detIdToIndexMap_[id()] = detIds_.size()-1;
                }
            }
            tpmHits_.resize(detIds_.size());
            rocHits_.resize(detIds_.size());
            theGeometryCacheID = geocache;
        }
        // zero out existing, if any
        for (auto &p : rocHits_) fill(p.begin(), p.end(), 0);
        for (auto &p : tpmHits_) fill(p.begin(), p.end(), 0);
        firstEventInLS_ = false;
        //waitNextLS_ = false;
    }
    //if (waitNextLS_) return;

    //unsigned int minNonZero = 0;
    iSetup.get<TrackerTopologyRcd>().get(theTrkTopo);
    edm::Handle<edmNew::DetSetVector<SiPixelCluster> > pixelClusters; 
    iEvent.getByToken(pixelClusterLabel_, pixelClusters); 
    for (const auto & detset : *pixelClusters) {
        DetId id(detset.detId());
        if (id.subdetId() != PixelSubdetector::SubDetector::PixelBarrel || theTrkTopo->layer(id) != 1) continue;
        int index = detIdToIndexMap_[id()];
        for (const SiPixelCluster & cluster : detset) { 
            int ix1 = (cluster.maxPixelRow()/80);
            int iy1 = (cluster.maxPixelCol()/52);
            int ix2 = (cluster.minPixelRow()/80);
            int iy2 = (cluster.minPixelCol()/52);
            if (ix1 != ix2 || iy1 != iy2) continue; // skip ones on the border for the moment
            int itpm = iy1/4 + ix1*2;
            int iroc = iy1 + ix1*8;
            tpmHits_[index][itpm]++;
            rocHits_[index][iroc]++;
            //if (rocHits_[index][iroc]> 0 && (minNonZero == 0 || minNonZero > rocHits_[index][iroc])) minNonZero = rocHits_[index][iroc];
        }
    }
    //printf("Min number of hits on non-empty ROC: %u\n", minNonZero);
    //if (minNonZero > maxHitsPerROC_) waitNextLS_ = true;
}

void 
PixelHitOccupancyByLumi::endLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const&) 
{
    run_  = iLumi.run();
    lumi_ = iLumi.luminosityBlock();
    unsigned int entries = 0;
    for (unsigned int i = 0, n = detIds_.size(); i < n; ++i) {
        detid_ = detIds_[i];
        for (tpm_ = 0; tpm_ < 4; ++tpm_) { 
            if (tpmHits_[i][tpm_] == 0) {
                roc_ = -1;
                tree_->Fill(); entries++;
            } else {
                for (roc_ = tpm_*4; roc_ < (tpm_+1)*4; ++roc_) { 
                    if (rocHits_[i][roc_] == 0) {
                        tree_->Fill();
                        entries++;
                    }
                }
            }
        }
    }
    printf("Wrote %u entries for this lumi\n", entries);
}




//define this as a plug-in
DEFINE_FWK_MODULE(PixelHitOccupancyByLumi);
