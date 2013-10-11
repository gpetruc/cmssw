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
