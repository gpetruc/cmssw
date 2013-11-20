#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "RecoTracker/HTPattern/interface/HTDebugger.h"
#include "RecoTracker/HTPattern/interface/TrackCandidateBuilderFromCluster.h"
#include "RecoTracker/DebugTools/interface/TrackMixingAssociator.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateAccessor.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2D.h"

#include <TGraph.h>
#include <TH2.h>
#include <TTree.h>
#include <map>

namespace {
    std::map<std::string, TTree *> trees3D_;
    Float_t rho_, z_, z0_, phi_; 
    Int_t layermask_;

    std::map<std::string, TTree *> treesTk_;
    Float_t ptTk_,etaTk_,phiTk_,zTk_;
    Int_t hitsTk_,layers3dTk_,hpTk_; 

    std::auto_ptr<TrackMixingAssociator> associator_;

    TTree *microClusterTree_ = 0;
    Int_t ntracks_; Float_t trackPt_; Int_t trackHits_, sharedHits_; Int_t seedhits_, hits_, morehits_; Float_t alpha_;
    Int_t npair_; Int_t hit1layer_, hit1class_, hit2layer_, hit2class_; Float_t hit12dr_, hit12deta_, hit12dphi_;
    Int_t hit12match_; Float_t hit12trackPt_, hit12trackEta_; const reco::Track * hit12track_;
    Int_t hit3layer_, hit3class_; Float_t hit3dr_, hit3deta_, hit3dphi_, thisalpha_; Int_t hit3pass_, hit3final_, hit3good_; 
    
}

void
HTDebugger::dumpHTHits3D(const std::string &name, const HTHits3D &hits, double zref)
{
    TTree * & tree = trees3D_[name];
    if (tree == 0) {
        edm::Service<TFileService> fs;
        tree = fs->make<TTree>( name.c_str()  , name.c_str() );
        tree->Branch("rho",&rho_,"rho/F");
        tree->Branch("z",&z_,"z/F");
        tree->Branch("z0",&z0_,"z0/F");
        tree->Branch("phi",&phi_,"phi/F");
        tree->Branch("layermask",&layermask_,"layermask/I");
    }
    z0_ = zref;
    for (unsigned int i = 0, n = hits.size(); i < n; ++i) {
        rho_ = hits.rho(i);
        z_ = hits.z(i);
        phi_ = hits.phi(i);
        layermask_ = hits.layermask(i);
        tree->Fill();
    }
}
void
HTDebugger::dumpTracks(const std::string &name, const std::vector<reco::Track> & tracks, double pT)
{
    
    TTree * & tree = treesTk_[name];
    if (tree == 0) {
        edm::Service<TFileService> fs;
        tree = fs->make<TTree>( name.c_str()  , name.c_str() );
        tree->Branch("z",&zTk_,"z/F");
        tree->Branch("phi",&phiTk_,"phi/F");
        tree->Branch("pt",&ptTk_,"pt/F");
        tree->Branch("eta",&etaTk_,"eta/F");
        tree->Branch("hits",&hitsTk_,"hits/I");
        tree->Branch("layers3d",&layers3dTk_,"layers3d/I");
        tree->Branch("hp",&hpTk_,"hp/I");
    }
    std::vector<std::pair<float, const reco::Track *> > tksorted;
    for (const reco::Track & tk : tracks) tksorted.push_back(std::make_pair(tk.eta(),&tk));
    std::sort(tksorted.begin(),tksorted.end());
    for (auto tkp : tksorted) {
        const reco::Track &tk = *tkp.second;
        ptTk_ = tk.pt();
        etaTk_ = tk.eta();
        phiTk_ = tk.phi(); if (phiTk_ < 0) phiTk_ += 2*M_PI;
        zTk_ = tk.vertex().Z();
        hitsTk_ = tk.hitPattern().numberOfValidHits();
        layers3dTk_ = tk.hitPattern().pixelLayersWithMeasurement() + tk.hitPattern().numberOfValidStripLayersWithMonoAndStereo(); 
        hpTk_ = tk.quality(reco::Track::highPurity); 
        if (pT == 0) {
            tree->Fill();
            if (ptTk_ > 0.250 && hpTk_) {
                printf("Track pt %8.4f +/- %7.2f, eta %+5.3f, phi %+5.3f, q %+2d, hits %2d, 3d layers %2d, vz %+5.3f, algo %2d (rho = %5.1f, z = %+5.1f --> rho = %5.1f, z = %+5.1f)\n", 
                        ptTk_, tk.ptError(), etaTk_, phiTk_, tk.charge(), hitsTk_, layers3dTk_, zTk_, tk.algo(), tk.innerPosition().Rho(), tk.innerPosition().Z(), tk.outerPosition().Rho(), tk.outerPosition().Z()); 
                for (trackingRecHit_iterator it = tk.recHitsBegin(), ed = tk.recHitsEnd(); it != ed; ++it) {
                    printf("\thit on detid %10d/%7d (valid? %1d)\n", (*it)->geographicalId().rawId(), TrackCandidateBuilderFromCluster::hitid(&**it), (*it)->isValid());
                }
            }
        }
    }
}



void
HTDebugger::dumpHTHitsSpher(const std::string &name, const HTHitsSpher &hits) {
    edm::Service<TFileService> fs;
    char ptitle[1024];
    char pname[1024];
    sprintf(ptitle, "%s (etabins %d, phibins %d, z0 %.3f, alpha %.3f);#eta;#phi", name.c_str(), hits.etabins(), hits.phibins(), hits.z0(), hits.alpha());
    sprintf(pname, "%s_etaphi", name.c_str());
    TGraph *g = fs->make<TGraph>(hits.size());
    g->SetNameTitle(pname,ptitle);
    for (unsigned int i = 0; i < hits.size(); ++i) {
        g->SetPoint(i, hits.eta(i), hits.phi(i));
    } 
    sprintf(pname, "%s_ietaiphi", name.c_str());
    g = fs->make<TGraph>(hits.size());
    g->SetNameTitle(pname,ptitle);
    for (unsigned int i = 0; i < hits.size(); ++i) {
        g->SetPoint(i, hits.ieta(i), hits.iphi(i));
    } 
 
    //fs->getBareDirectory()->WriteTObject(g, pname);
}

void
HTDebugger::dumpHTHitMap(const std::string &name, const HTHitMap &map) 
{
    edm::Service<TFileService> fs;
    char ptitle[1024];
    char pname[1024];
    sprintf(ptitle, "%s (etabins %d, phibins %d);i#eta;i#phi", name.c_str(), map.etabins(), map.phibins());
    sprintf(pname, "%s_hits", name.c_str());
    TH2D *hits = fs->make<TH2D>(pname, ptitle, map.etabins(), -3., 3., map.phibins(), 0, 2*M_PI );
    sprintf(pname, "%s_layers", name.c_str());
    TH2D *layers = fs->make<TH2D>(pname, ptitle, map.etabins(), -3., 3., map.phibins(), 0, 2*M_PI );
    for (unsigned int i = 0; i < map.etabins(); ++i) {
        for (unsigned int j = 0; j < map.phibins(); ++j) {
            hits->SetBinContent(i+1, j+1, map.get(i,j).nhits());
            layers->SetBinContent(i+1, j+1, map.get(i,j).nlayers());
        }
    } 

}

namespace {
    void addHits(const HTCell &cell, const HTHitsSpher &hits, TGraph *gep, TGraph *grz, unsigned int &nhits)  {
        gep->Set(nhits+cell.nhits());
        grz->Set(nhits+cell.nhits());
        for (unsigned int h : cell.hits()) {
            gep->SetPoint(nhits, hits.eta(h), hits.phi(h));
            grz->SetPoint(nhits, hits.rho(h), hits.z(h));
            nhits++;
        }
    }
}
void
HTDebugger::dumpHTClusters(const std::string &name, const HTHitMap &map, const HTHitsSpher &hits, unsigned int minlayers, unsigned int minmorelayers) 
{
    edm::Service<TFileService> fs;
    char ptitle[1024];
    char pname[1024];
    sprintf(ptitle, "%s (etabins %d, phibins %d);i#eta;i#phi", name.c_str(), map.etabins(), map.phibins());
    sprintf(pname, "%s_layers", name.c_str());
    TH2D *layers = fs->make<TH2D>(pname, ptitle, map.etabins(), -3., 3., map.phibins(), 0, 2*M_PI );
    sprintf(pname, "%s_hits_etaphi", name.c_str());
    TObjArray *geps = new TObjArray(); geps->SetName(pname);
    sprintf(pname, "%s_hits_rhoz", name.c_str());
    TObjArray *grzs = new TObjArray(); grzs->SetName(pname);
    geps->SetOwner(true);
    grzs->SetOwner(true);
    unsigned int color = 1;
    std::array<const HTCell *, 8> cells; unsigned int ncells;
    for (auto & cluster : map.clusters()) {
        if (cluster.nlayers() < minlayers) continue;
        if (cluster.nmorelayers() < minmorelayers) continue;
        unsigned int ieta = cluster.ieta(), iphi = cluster.iphi();
        layers->SetBinContent(ieta+1, iphi+1, minlayers ? cluster.nlayers() : cluster.nmorelayers());
        unsigned int nhits = 0;
        TGraph *gep = new TGraph(); 
        TGraph *grz = new TGraph(); 
        addHits(map.get(ieta, iphi), hits, gep, grz, nhits);
        map.getNeighbours(ieta, iphi, cells, ncells);
        for (unsigned int i = 0; i < ncells; ++i) {
            addHits(*cells[i], hits, gep, grz, nhits);
        }
        grz->SetMarkerStyle(7);
        gep->SetMarkerStyle(7);
        grz->SetMarkerColor(color);
        gep->SetMarkerColor(color);
        grzs->AddLast(grz);
        geps->AddLast(gep);
        color++; if (color == 10) color = 1;
    }
    fs->getBareDirectory()->WriteTObject(geps, 0, "SingleKey");
    fs->getBareDirectory()->WriteTObject(grzs, 0, "SingleKey");
}

void
HTDebugger::registerTracks(const std::vector<reco::Track> & tracks) 
{
    associator_.reset(new TrackMixingAssociator());
    associator_->registerTrackEvent(1, tracks);
}

void
HTDebugger::registerClustersAsSeeds(int step, const TrajectorySeedCollection &seeds)
{
    if (associator_.get() == 0) return;
    associator_->registerSeedEvent(step, seeds);
}

void 
HTDebugger::registerTrackCandidates(const std::vector<TrackCandidate> &tcs)
{
    if (associator_.get() == 0) return;
    associator_->registerTrackCandEvent(3, tcs);
}

void 
HTDebugger::printAssociatedTracks(const TrajectorySeed &cluster)
{
    if (associator_.get() == 0) return;
    std::vector<TrackMixingAssociator::TrackAssociation> assoc;
    associator_->associateToTracks(cluster, assoc);
    for (const auto &a : assoc) {
        if (a.sharedHits <= 2 && assoc.front().sharedHits > 2) break;
        float phi = a.track->phi(); if (phi < 0) phi += 2*M_PI;
        printf("\t\t-> track pt %8.4f, eta %+5.3f, phi %+5.3f, q %+2d, hits %2d, 3d layers %2d: associated layers %2d, hits %2d (%5.1f%%)\n", 
            a.track->pt(), a.track->eta(), phi, a.track->charge(), a.track->hitPattern().numberOfValidHits(), a.track->hitPattern().pixelLayersWithMeasurement() + a.track->hitPattern().numberOfValidStripLayersWithMonoAndStereo(), a.sharedHitPattern.trackerLayersWithMeasurement(), a.sharedHits, (100.0*a.sharedHits)/a.track->found());
    }
}

void 
HTDebugger::printAssociatedTracks(const Trajectory &traj, int minHits)
{
    if (associator_.get() == 0) return;
    std::vector<TrackMixingAssociator::TrackAssociation> assoc;
    associator_->associateToTracks(traj, assoc);
    for (const auto &a : assoc) { if (int(a.sharedHits) < minHits) continue;
        if (a.sharedHits <= 2 && assoc.front().sharedHits > 2) break;
        float phi = a.track->phi(); if (phi < 0) phi += 2*M_PI;
        printf("\t\t-> track pt %8.4f, eta %+5.3f, phi %+5.3f, q %+2d, hits %2d, 3d layers %2d: associated layers %2d, hits %2d (%5.1f%%)\n", 
            a.track->pt(), a.track->eta(), phi, a.track->charge(), a.track->hitPattern().numberOfValidHits(), a.track->hitPattern().pixelLayersWithMeasurement() + a.track->hitPattern().numberOfValidStripLayersWithMonoAndStereo(), a.sharedHitPattern.trackerLayersWithMeasurement(), a.sharedHits, (100.0*a.sharedHits)/a.track->found());
    }
}

void
HTDebugger::printBackAssociation(const std::vector<reco::Track> & tracks, const std::vector<TrackCandidate> &tcs, const TrackingGeometry &geometry, const MagneticField &bfield, float ptMin) 
{
    if (associator_.get() == 0) return;
    std::vector<TrackMixingAssociator::SeedAssociation> assoc;
    std::vector<TrackMixingAssociator::TrackCandAssociation> tcassoc;
    printf("\n\n========================================== CKF -> HT SUMMARY ==================================");
    for (const reco::Track &t : tracks) {
        if (t.pt() < ptMin) continue;
        float phi = t.phi(); if (phi < 0) phi += 2*M_PI;
        printf("\ntrack pt %8.4f, eta %+5.3f, phi %+5.3f, q %+2d, hits %2d, 3d layers %2d. Associations:\n", 
                    t.pt(), t.eta(), phi, t.charge(), t.hitPattern().numberOfValidHits(), t.hitPattern().pixelLayersWithMeasurement() + t.hitPattern().numberOfValidStripLayersWithMonoAndStereo());
        associator_->associateToSeed(t, assoc);
        for (const auto &a : assoc) {
            if (a.eventId != 1) continue; if (a.sharedHits <= 2) continue;
            printf("\t-> seed from step %d with front detid %d (%2d associated hits)\n", a.eventId, a.track->startingState().detId(), a.sharedHits);
        }  
        for (const auto &a : assoc) {
            if (a.eventId == 1) continue; if (a.sharedHits <= 2) continue;
            printf("\t-> seed from step %d with front detid %d (%2d associated hits)\n", a.eventId, a.track->startingState().detId(), a.sharedHits);
        }  
        associator_->associateToTrackCands(t, tcassoc);
        for (const auto &a : tcassoc) {
            if (a.sharedHits <= 2) continue;
            printf("\t-> track candidate with %d hits, starting detid %d (%2d associated hits, %3.0f%%)\n", int(std::distance(a.track->recHits().first, a.track->recHits().second)), a.track->trajectoryStateOnDet().detId(), a.sharedHits, a.sharedHits*100.0/t.numberOfValidHits());
        }  
    }

    printf("\n========================================== HT -> CKF SUMMARY ==================================");
    std::vector<TrackMixingAssociator::TrackAssociation> tassoc;
    for (const TrackCandidate &tc : tcs) {
        const GeomDet *det = geometry.idToDet(DetId(tc.trajectoryStateOnDet().detId())); if (det == 0) continue;
        TrajectoryStateOnSurface tsos = trajectoryStateTransform::transientState(tc.trajectoryStateOnDet(), &det->surface(), &bfield);
        float phipos = tsos.globalMomentum().phi(); if (phipos < 0) phipos += 2*M_PI;
        printf("\ntrack candidate pt %8.4f +/- %7.2f, eta %+5.3f, phi %+5.3f, q %+2d, hits %2d, starting detid %d\n",  
                tsos.globalMomentum().perp(), TrajectoryStateAccessor(*tsos.freeState()).inversePtError()*tsos.globalMomentum().perp2(), 
                tsos.globalMomentum().eta(), phipos, 
                tsos.charge(), int(std::distance(tc.recHits().first, tc.recHits().second)), tc.trajectoryStateOnDet().detId()),
        associator_->associateToTracks(tc, tassoc);
        for (const auto &a : tassoc) { if (a.sharedHits <= 2) continue; 
            float phi = a.track->phi(); if (phi < 0) phi += 2*M_PI;
            printf("\t-> track pt %8.4f, eta %+5.3f, phi %+5.3f, q %+2d, hits %2d, 3d layers %2d: associated hits: %2d (%3.0f%%)\n", 
                    a.track->pt(), a.track->eta(), phi, a.track->charge(), a.track->hitPattern().numberOfValidHits(), a.track->hitPattern().pixelLayersWithMeasurement() + a.track->hitPattern().numberOfValidStripLayersWithMonoAndStereo(), a.sharedHits, a.sharedHits*100.0/a.track->numberOfValidHits());
        }
    }
}

void HTDebugger::printAssociatedTracks(const TrackingRecHit *hit1, const TrackingRecHit *hit2)
{
    if (associator_.get() == 0) return;
    std::vector<const TrackingRecHit *> hits(2, hit1);
    hits[1] = hit2;
    std::vector<TrackMixingAssociator::TrackAssociation> assoc;
    associator_->associateToTracks(hits, assoc);
    for (const auto &a : assoc) {
        float phi = a.track->phi(); if (phi < 0) phi += 2*M_PI;
        if (a.sharedHits < 2) break;
        if (a.sharedHitPattern.trackerLayersWithMeasurement() < 2) continue;
        printf("\t\t-> track pt %8.4f, eta %+5.3f, phi %+5.3f, q %+2d, hits %2d, 3d layers %2d: associated layers %2d, hits %2d\n", 
            a.track->pt(), a.track->eta(), phi, a.track->charge(), a.track->hitPattern().numberOfValidHits(), a.track->hitPattern().pixelLayersWithMeasurement() + a.track->hitPattern().numberOfValidStripLayersWithMonoAndStereo(), a.sharedHitPattern.trackerLayersWithMeasurement(), a.sharedHits);
    }
}

bool 
HTDebugger::thirdHitOnSameTrack(const TrackingRecHit *hit1, const TrackingRecHit *hit2, const TrackingRecHit *hit3)
{
    if (associator_.get() == 0) return false;
    std::vector<const TrackingRecHit *> hits(2, hit1);
    hits[1] = hit2;
    std::vector<TrackMixingAssociator::TrackAssociation> assoc, assoc2;
    associator_->associateToTracks(hits, assoc);
    std::vector<const reco::Track *> assoc2l;
    for (const auto &a : assoc) {
        if (a.sharedHitPattern.trackerLayersWithMeasurement() >= 2) assoc2l.push_back(a.track);
    }
    if (assoc2l.empty()) return false;
    hits.clear(); hits.push_back(hit3);
    associator_->associateToTracks(hits, assoc2);
    for (const auto &a : assoc2) {
        if (std::find(assoc2l.begin(), assoc2l.end(), a.track) != assoc2l.end()) {
            return true;
        }
    }
    return false;
}

void HTDebugger::beginLoggingCluster(const TrajectorySeed &cluster, unsigned int seedhits, unsigned int hits, unsigned int morehits, double alpha)
{
    if (associator_.get() == 0) return;
    if (microClusterTree_ == 0) {
        edm::Service<TFileService> fs;
        microClusterTree_ = fs->make<TTree>( "microClusters", "" );
        microClusterTree_->Branch("ntracks", &ntracks_, "ntracks/I");
        microClusterTree_->Branch("trackPt", &trackPt_, "trackPt/F");
        microClusterTree_->Branch("trackHits", &trackHits_, "trackHits/I");
        microClusterTree_->Branch("sharedHits", &sharedHits_, "sharedHits/I");
        microClusterTree_->Branch("seedhits", &seedhits_, "seedhits/I");
        microClusterTree_->Branch("hits", &hits_, "hits/I");
        microClusterTree_->Branch("morehits", &morehits_, "morehits/I");
        microClusterTree_->Branch("alpha", &alpha_, "alpha/F");
        microClusterTree_->Branch("npair", &npair_, "npair/I");
        microClusterTree_->Branch("hit1layer", &hit1layer_, "hit1layer/I");
        microClusterTree_->Branch("hit1class", &hit1class_, "hit1class/I");
        microClusterTree_->Branch("hit2layer", &hit2layer_, "hit2layer/I");
        microClusterTree_->Branch("hit2class", &hit2class_, "hit2class/I");
        microClusterTree_->Branch("hit12dr", &hit12dr_, "hit12dr/F");
        microClusterTree_->Branch("hit12deta", &hit12deta_, "hit12deta/F");
        microClusterTree_->Branch("hit12dphi", &hit12dphi_, "hit12dphi/F");
        microClusterTree_->Branch("hit12match", &hit12match_, "hit12match/I");
        microClusterTree_->Branch("hit12trackPt", &hit12trackPt_, "hit12trackPt/F");
        microClusterTree_->Branch("hit12trackEta", &hit12trackEta_, "hit12trackEta/F");
        microClusterTree_->Branch("hit3layer", &hit3layer_, "hit3layer/I");
        microClusterTree_->Branch("hit3class", &hit3class_, "hit3class/I");
        microClusterTree_->Branch("hit3dr", &hit3dr_, "hit3dr/F");
        microClusterTree_->Branch("hit3dphi", &hit3dphi_, "hit3dphi/F");
        microClusterTree_->Branch("hit3deta", &hit3deta_, "hit3deta/F");
        microClusterTree_->Branch("thisalpha", &thisalpha_, "thisalpha/F");
        microClusterTree_->Branch("hit3pass", &hit3pass_, "hit3pass/I");
        microClusterTree_->Branch("hit3final", &hit3final_, "hit3final/I");
        microClusterTree_->Branch("hit3good", &hit3good_, "hit3good/I");
        
    }
    ntracks_ = 0; npair_ = 0; seedhits_ = seedhits; hits_ = hits; morehits_ = morehits_; alpha_ = alpha;
    std::vector<TrackMixingAssociator::TrackAssociation> assoc;
    associator_->associateToTracks(cluster, assoc);
    for (const auto &a : assoc) {
        if (a.sharedHits <= 2) continue;
        ntracks_++;
        if (ntracks_ == 1) {
            trackPt_ = a.track->pt();
        }
    }
    hit1layer_ = -1; microClusterTree_->Fill();
}
void HTDebugger::logSeedingPair(const TrackingRecHit *hit1, int layer1, int class1, const TrackingRecHit *hit2, int layer2, int class2, float dr, float dphi, float deta)
{
    if (associator_.get() == 0 || microClusterTree_ == 0) return;
    std::vector<const TrackingRecHit *> hits(2, hit1);
    hits[1] = hit2;
    std::vector<TrackMixingAssociator::TrackAssociation> assoc;
    associator_->associateToTracks(hits, assoc);
    npair_++; hit1layer_ = layer1; hit1class_ = class1; hit2layer_ = layer2; hit2class_ = class2; hit12dr_ = dr; hit12deta_ = deta; hit12dphi_ = dphi;
    hit12track_ = 0; hit12match_ = 0; 
    for (const auto &a : assoc) {
        if (a.sharedHitPattern.trackerLayersWithMeasurement() >= 2) {
            hit12track_ = a.track; hit12match_ = 1; hit12trackPt_ = a.track->pt(); hit12trackEta_ = a.track->eta();
        }
    }
    hit3layer_ = -1; microClusterTree_->Fill();
}
void HTDebugger::logThirdHit(const TrackingRecHit *hit3, int layer3, int class3, float d2c, float dphic, float detac, int pass, float thisalpha, bool finalsel)
{
    if (associator_.get() == 0 || microClusterTree_ == 0) return;
    hit3layer_ = layer3; hit3class_ = class3; hit3dr_ = d2c; hit3dphi_ = dphic; hit3deta_ = detac; thisalpha_ = thisalpha; hit3final_ = finalsel; hit3pass_ = pass;
    std::vector<const TrackingRecHit *> hits(1, hit3);
    std::vector<TrackMixingAssociator::TrackAssociation> assoc;
    associator_->associateToTracks(hits, assoc);
    hit3good_ = 0;
    for (const auto &a : assoc) {
        if (a.track == hit12track_ && hit12track_ != 0) {
            hit3good_ = 1;
            break;
        }
    }
    microClusterTree_->Fill();
    
}


void HTDebugger::debugHelixParameters(const std::vector<const TrackingRecHit *> &hits, float alpha, float beta, float phi0, float eta0, float bfieldAtZero, GlobalPoint vertex) 
{
    if (associator_.get() == 0) return;
    std::vector<TrackMixingAssociator::TrackAssociation> assoc;
    associator_->associateToTracks(hits, assoc);
    for (const auto &a : assoc) {
        if (a.sharedHits >= 3) {
            float pt = alpha ? std::abs(0.5f*0.003f*bfieldAtZero/alpha) : 100.0f;
            float dz = beta * cosh(eta0); 
            float tphi = a.track->phi(); if (tphi < 0) tphi += 2*M_PI;
            float dphi = std::abs(phi0 - tphi); if (dphi > M_PI) dphi = 2*M_PI - dphi;
            float tdz  = a.track->dz(reco::Track::Point(vertex));
            printf("helix  pt = %8.4f, eta = %+5.3f, phi = %+5.3f, dz = %+5.3f; track pt %8.4f, eta = %+5.3f, phi = %+5.3f, dz = %+5.3f;  d(1/pt) = %6.4f, deta = %6.4f, dphi = %6.4f, dz = %6.4f\n",
                    pt, eta0, phi0, dz, a.track->pt(), a.track->eta(), tphi, tdz, std::abs(1/pt - 1/a.track->pt()), std::abs(eta0-a.track->eta()), dphi, std::abs(dz-tdz));
            break;
        }
    }
}

const reco::Track * HTDebugger::matchCandidate(const std::vector<const TrackingRecHit *> &hits) 
{
    if (associator_.get() == 0) return 0;
    std::vector<TrackMixingAssociator::TrackAssociation> assoc;
    associator_->associateToTracks(hits, assoc);
    for (const auto &a : assoc) {
        if (a.sharedHits >= 3 && a.sharedHitPattern.trackerLayersWithMeasurement() >= 3) {
            float phi = a.track->phi(); if (phi < 0) phi += 2*M_PI;
            printf("\t\t-> track pt %8.4f, eta %+5.3f, phi %+5.3f, q %+2d, hits %2d, 3d layers %2d: associated layers %2d, hits %2d (%5.1f%%)\n", 
                a.track->pt(), a.track->eta(), phi, a.track->charge(), a.track->hitPattern().numberOfValidHits(), a.track->hitPattern().pixelLayersWithMeasurement() + a.track->hitPattern().numberOfValidStripLayersWithMonoAndStereo(), a.sharedHitPattern.trackerLayersWithMeasurement(), a.sharedHits, (100.0*a.sharedHits)/a.track->found());
            return a.track;
        }
    }
    return 0;
}
int HTDebugger::isMatchedToTrack(const TrackingRecHit *hit, const reco::Track *tk) {
    if (associator_.get() == 0) return 0;
    std::vector<const TrackingRecHit *> hits(1, hit); std::vector<SiStripRecHit2D> split;
    if (typeid(*hit) == typeid(SiStripMatchedRecHit2D)) {
        const SiStripMatchedRecHit2D &m = static_cast<const SiStripMatchedRecHit2D &>(*hit);
        split.push_back(m.monoHit());
        split.push_back(m.stereoHit());
        hits.resize(2);
        hits[0] = & split[0];
        hits[1] = & split[1];
    }
    std::vector<TrackMixingAssociator::TrackAssociation> assoc;
    associator_->associateToTracks(hits, assoc);
    int ret = 0;
    for (const auto &a : assoc) {
        if (a.track == tk) { 
            if (a.sharedHits == hits.size()) return 2;
            else ret = 1;
        }
    }
    return ret;
}
