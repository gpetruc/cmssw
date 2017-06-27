
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
        void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) { firstEventInLS_ = true; waitNextLS_ = false; }
        void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) ;

        // ----------member data ---------------------------
        const edm::EDGetTokenT<MeasurementTrackerEvent> tracker_;
        unsigned int maxHitsPerROC_;

        edm::ESHandle<TrackerTopology> theTrkTopo;

        std::vector<unsigned int> pxb1ModuleIndices_;
        std::map<std::pair<int,int>,unsigned int> rocHits_;
        std::map<std::pair<int,int>,unsigned int> tpmHits_;
        bool firstEventInLS_, waitNextLS_;

        TTree *tree_;
        int run_, lumi_;    
        int detid_, tpm_, roc_;        

};

//
// constructors and destructor
//
PixelHitOccupancyByLumi::PixelHitOccupancyByLumi(const edm::ParameterSet& iConfig):
    tracker_(consumes<MeasurementTrackerEvent>(iConfig.getParameter<edm::InputTag>("tracker"))),
    maxHitsPerROC_(iConfig.getParameter<uint32_t>("maxHitsPerROC")),
    firstEventInLS_(true), waitNextLS_(false)
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
    using namespace edm;
    using namespace std;

    iSetup.get<TrackerTopologyRcd>().get(theTrkTopo);
    Handle<MeasurementTrackerEvent> tracker;

    iEvent.getByToken(tracker_, tracker);
    const auto & pixelData = tracker->pixelData();
    unsigned int ndet = pixelData.nDet();
    if (firstEventInLS_) {
        pxb1ModuleIndices_.clear();
        // zero out existing, if any
        for (auto &p : rocHits_) p.second = 0;
        for (auto &p : tpmHits_) p.second = 0;
        // make sure all ids are here
        for (unsigned int idet = 0; idet < ndet; ++idet) {
            DetId id(pixelData.id(idet));
            if (id.subdetId() == PixelSubdetector::SubDetector::PixelBarrel && theTrkTopo->layer(id) == 1) {
                for (unsigned int i = 0; i < 16; ++i) {
                    rocHits_[make_pair(id(),i)] = 0;
                }
                for (unsigned int i = 0; i < 4; ++i) {
                    tpmHits_[make_pair(id(),i)] = 0;
                }
            }
            pxb1ModuleIndices_.push_back(idet);
        }
        firstEventInLS_ = false;
        waitNextLS_ = false;
    }
    if (waitNextLS_) return;
    unsigned int minNonZero = 0;
    for (unsigned int i = 0, n = pxb1ModuleIndices_.size(); i < n; ++i) {
        int idet = pxb1ModuleIndices_[i];
        if (pixelData.empty(idet) || !pixelData.isActive(idet)) continue;
        int detid = pixelData.id(idet);
        for (const SiPixelCluster & cluster : pixelData.detSet(idet)) { 
            int ix1 = (cluster.maxPixelRow()/80);
            int iy1 = (cluster.maxPixelCol()/52);
            int ix2 = (cluster.minPixelRow()/80);
            int iy2 = (cluster.minPixelCol()/52);
            if (ix1 != ix2 || iy1 != iy2) continue; // skip ones on the border for the moment
            int itpm = iy1/4 + ix1*2;
            int iroc = iy1 + ix1*8;
            tpmHits_[make_pair(detid,itpm)]++;
            unsigned int & hits = rocHits_[make_pair(detid,iroc)]; hits++;
            if (hits > 0 && (minNonZero == 0 || minNonZero > hits)) minNonZero = hits;
        }
    }
    printf("Min number of hits on non-empty ROC: %u\n", minNonZero);
    if (minNonZero > maxHitsPerROC_) waitNextLS_ = true;
}

void 
PixelHitOccupancyByLumi::endLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const&) 
{
    run_  = iLumi.run();
    lumi_ = iLumi.luminosityBlock();

    for (const auto &p : tpmHits_) {
        detid_ =  p.first.first;
        tpm_ = p.first.second;
        if (p.second == 0) {
            roc_ = -1;
            tree_->Fill();
        } else {
            tpm_ = -1;
            for (unsigned int i = 0; i < 4; ++i) {
                unsigned int hits = rocHits_[std::make_pair(p.first.first,p.first.second*4+i)];
                if (hits == 0) { roc_ = i; tree_->Fill(); }
            }
        }
    }
}




//define this as a plug-in
DEFINE_FWK_MODULE(PixelHitOccupancyByLumi);
