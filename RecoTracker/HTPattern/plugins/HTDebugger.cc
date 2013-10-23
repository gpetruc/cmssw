#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "RecoTracker/HTPattern/plugins/HTDebugger.h"
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
    for (auto tk : tracks) {
        ptTk_ = tk.pt();
        etaTk_ = tk.eta();
        phiTk_ = tk.phi(); if (phiTk_ < 0) phiTk_ += 2*M_PI;
        zTk_ = tk.vertex().Z();
        hitsTk_ = tk.hitPattern().numberOfValidHits();
        layers3dTk_ = tk.hitPattern().pixelLayersWithMeasurement() + tk.hitPattern().numberOfValidStripLayersWithMonoAndStereo(); 
        hpTk_ = tk.quality(reco::Track::highPurity); 
        if (pT == 0) {
            tree->Fill();
            if (ptTk_ > 0.5 && hpTk_) {
                printf("Track pt %8.4f, eta %+5.3f, phi %+5.3f, q %+2d, hits %2d, 3d layers %2d, vz %+5.3f\n", ptTk_, etaTk_, phiTk_, tk.charge(), hitsTk_, layers3dTk_, zTk_); 
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

