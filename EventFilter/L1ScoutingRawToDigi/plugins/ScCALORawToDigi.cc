#include "EventFilter/L1ScoutingRawToDigi/plugins/ScCALORawToDigi.h"

ScCaloRawToDigi::ScCaloRawToDigi(const edm::ParameterSet& iConfig) {
  using namespace edm;
  using namespace scoutingRun3;
  srcInputTag  = iConfig.getParameter<InputTag>( "srcInputTag" );
  debug = iConfig.getUntrackedParameter<bool>("debug", false);

  // produces<l1t::JetBxCollection>().setBranchAlias( "JetBxCollection" );
  // produces<l1t::TauBxCollection>().setBranchAlias( "TauBxCollection" );
  // produces<l1t::EGammaBxCollection>().setBranchAlias( "EGammaBxCollection" );
  // produces<l1t::EtSumBxCollection>().setBranchAlias( "EtSumBxCollection" );

  produces<JetOrbitCollection>().setBranchAlias( "JetOrbitCollection" );
  produces<TauOrbitCollection>().setBranchAlias( "TauOrbitCollection" );
  produces<EGammaOrbitCollection>().setBranchAlias( "EGammaOrbitCollection" );
  produces<EtSumOrbitCollection>().setBranchAlias( "EtSumOrbitCollection" );
  
  rawToken = consumes<SRDCollection>(srcInputTag);
  
  // bx_jets.reserve(12);
  // bx_taus.reserve(12);
  // bx_eGammas.reserve(12);
  // bx_etSums.reserve(12);
  dummyLVec_.reset( new ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>() );
}

ScCaloRawToDigi::~ScCaloRawToDigi() {};

void ScCaloRawToDigi::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace scoutingRun3;

  Handle<SRDCollection> ScoutingRawDataCollection;
  iEvent.getByToken( rawToken, ScoutingRawDataCollection );

  const FEDRawData& sourceRawData = ScoutingRawDataCollection->FEDData(SDSNumbering::CaloSDSID);
  size_t orbitSize = sourceRawData.size();

  // std::unique_ptr<l1t::JetBxCollection> unpackedJets(new l1t::JetBxCollection);
  // std::unique_ptr<l1t::TauBxCollection> unpackedTaus(new l1t::TauBxCollection);
  // std::unique_ptr<l1t::EGammaBxCollection> unpackedEGammas(new l1t::EGammaBxCollection);
  // std::unique_ptr<l1t::EtSumBxCollection> unpackedEtSums(new l1t::EtSumBxCollection);

  std::unique_ptr<JetOrbitCollection> unpackedJets(new JetOrbitCollection);
  std::unique_ptr<TauOrbitCollection> unpackedTaus(new TauOrbitCollection);
  std::unique_ptr<EGammaOrbitCollection> unpackedEGammas(new EGammaOrbitCollection);
  std::unique_ptr<EtSumOrbitCollection> unpackedEtSums(new EtSumOrbitCollection);
  
  if((sourceRawData.size()==0) && debug ){
    std::cout << "No raw data for CALO FED\n";  
  }

  unpackOrbit(
    unpackedJets.get(), unpackedTaus.get(),
    unpackedEGammas.get(), unpackedEtSums.get(),
    sourceRawData.data(), orbitSize
  ); 

  // store collections in the event
  iEvent.put( std::move(unpackedJets) );
  iEvent.put( std::move(unpackedTaus) );
  iEvent.put( std::move(unpackedEGammas) );
  iEvent.put( std::move(unpackedEtSums) );
}

void ScCaloRawToDigi::unpackOrbit(
  // l1t::JetBxCollection* jets, l1t::TauBxCollection* taus,
  // l1t::EGammaBxCollection* eGammas, l1t::EtSumBxCollection* etSums,
  scoutingRun3::JetOrbitCollection* jets, scoutingRun3::TauOrbitCollection* taus,
  scoutingRun3::EGammaOrbitCollection* eGammas, scoutingRun3::EtSumOrbitCollection* etSums,
  const unsigned char* buf, size_t len
  ){
  
  using namespace scoutingRun3;
  
  size_t pos = 0;

  // jets->setBXRange(0,3565);
  // taus->setBXRange(0,3565);
  // eGammas->setBXRange(0,3565);
  // etSums->setBXRange(0,3565);

  while (pos < len) {
    
    assert(pos+ (4+4+4+56*4) <= len); //sizeof(demux::block)

    demux::block *bl = (demux::block *)(buf + pos);
    pos += 4+4+4+56*4;

    assert(pos <= len);
    uint32_t orbit = bl->orbit & 0x7FFFFFFF;
    uint32_t bx = bl->bx;

    if(debug) {
      std::cout << " CALO Orbit " << orbit << ", BX -> "<< bx << std::endl;
    }
    
    // reset vectors
    // bx_jets.clear();
    // bx_taus.clear();
    // bx_eGammas.clear();
    // bx_etSums.clear(); 

    int32_t ET(0), Eta(0), Phi(0), Iso(0);

    // unpack jets from first link
    for (uint32_t i=0; i<6; i++) {
      ET = ((bl->jet1[i] >> demux::shiftsJet::ET)  & demux::masksJet::ET);
      
      if (ET != 0) {
        Eta = ((bl->jet1[i] >> demux::shiftsJet::eta) & demux::masksJet::eta);
        Phi = ((bl->jet1[i] >> demux::shiftsJet::phi) & demux::masksJet::phi);
        Iso = 0;

        if (Eta > 127) Eta = Eta - 256;
        
        l1t::Jet jet(*dummyLVec_, ET, Eta, Phi, Iso);
        jets->push_back(bx,jet);
        // bx_jets.emplace_back(l1t::Jet(*dummyLVec_, ET, Eta, Phi, Iso));
      } 
    } // end link1 jet unpacking loop
    
    // unpack jets from second link
    for (uint32_t i=0; i<6; i++) {
      ET = ((bl->jet2[i] >> demux::shiftsJet::ET)  & demux::masksJet::ET);
      
      if (ET != 0) {
        Eta = ((bl->jet2[i] >> demux::shiftsJet::eta) & demux::masksJet::eta);
        Phi = ((bl->jet2[i] >> demux::shiftsJet::phi) & demux::masksJet::phi);
        Iso = 0;

        if (Eta > 127) Eta = Eta - 256;
        
        l1t::Jet jet(*dummyLVec_, ET, Eta, Phi, Iso);
        jets->push_back(bx,jet);
        // bx_jets.emplace_back(l1t::Jet(*dummyLVec_, ET, Eta, Phi, Iso));
      } 
    } // end link1 jet unpacking loop

    // add jets to the event collection 
    //jets->push_back(bx, bx_jets);


    // unpack eg from first link
    for (uint32_t i=0; i<6; i++) {
      ET = ((bl->egamma1[i] >> demux::shiftsEGamma::ET)  & demux::masksEGamma::ET);
      if (ET != 0) {
        Eta   = ((bl->egamma1[i] >> demux::shiftsEGamma::eta) & demux::masksEGamma::eta);
        Phi   = ((bl->egamma1[i] >> demux::shiftsEGamma::phi) & demux::masksEGamma::phi);
        Iso   = ((bl->egamma1[i] >> demux::shiftsEGamma::iso) & demux::masksEGamma::iso);

        if (Eta > 127) Eta = Eta - 256;
        
        l1t::EGamma eGamma(*dummyLVec_, ET, Eta, Phi, 0, Iso);
        eGammas->push_back(bx, eGamma);
        // bx_eGammas.emplace_back(l1t::EGamma(*dummyLVec_, ET, Eta, Phi, 0, Iso));
      }
    } // end eg link 1

    // unpack eg from second link link
    for (uint32_t i=0; i<6; i++) {
      ET = ((bl->egamma2[i] >> demux::shiftsEGamma::ET)  & demux::masksEGamma::ET);
      if (ET != 0) {
        Eta   = ((bl->egamma2[i] >> demux::shiftsEGamma::eta) & demux::masksEGamma::eta);
        Phi   = ((bl->egamma2[i] >> demux::shiftsEGamma::phi) & demux::masksEGamma::phi);
        Iso   = ((bl->egamma2[i] >> demux::shiftsEGamma::iso) & demux::masksEGamma::iso);

        if (Eta > 127) Eta = Eta - 256;
        
        l1t::EGamma eGamma(*dummyLVec_, ET, Eta, Phi, 0, Iso);
        eGammas->push_back(bx, eGamma);
        // bx_eGammas.emplace_back(l1t::EGamma(*dummyLVec_, ET, Eta, Phi, 0, Iso));
      }

    } // end of eg unpacker

    // add eg to the event
    // eGammas->push_back(bx, bx_eGammas);


    // unpack taus from first link
    for (uint32_t i=0; i<6; i++) { 
      ET = ((bl->tau1[i] >> demux::shiftsTau::ET)  & demux::masksTau::ET);
      if (ET != 0) {
          Eta   = ((bl->tau1[i] >> demux::shiftsTau::eta) & demux::masksTau::eta);
          Phi   = ((bl->tau1[i] >> demux::shiftsTau::phi) & demux::masksTau::phi);
          Iso   = ((bl->tau1[i] >> demux::shiftsTau::iso) & demux::masksTau::iso);

          if (Eta > 127) Eta = Eta - 256; 

          l1t::Tau tau(*dummyLVec_, ET, Eta, Phi, 0, Iso);
          taus->push_back(bx, tau);
          //bx_taus.emplace_back(l1t::Tau(*dummyLVec_, ET, Eta, Phi, 0, Iso));
      }
    } // end tau link 1

    // unpack taus from second link
    for (uint32_t i=0; i<6; i++) { 
      ET = ((bl->tau2[i] >> demux::shiftsTau::ET)  & demux::masksTau::ET);
      if (ET != 0) {
          Eta   = ((bl->tau2[i] >> demux::shiftsTau::eta) & demux::masksTau::eta);
          Phi   = ((bl->tau2[i] >> demux::shiftsTau::phi) & demux::masksTau::phi);
          Iso   = ((bl->tau2[i] >> demux::shiftsTau::iso) & demux::masksTau::iso);

          if (Eta > 127) Eta = Eta - 256;   

          l1t::Tau tau(*dummyLVec_, ET, Eta, Phi, 0, Iso);
          taus->push_back(bx, tau);
          //bx_taus.emplace_back(l1t::Tau(*dummyLVec_, ET, Eta, Phi, 0, Iso));
      }
    } // end tau unpacker

    // add taus to event
    // taus->push_back(bx, bx_taus);   

    // unpack et sums
    int32_t ETEt(0), ETEttem(0), //ETMinBiasHF(0),
            HTEt(0), HTtowerCount(0), //HTMinBiasHF(0),
            ETmissEt(0), ETmissPhi(0), ETmissASYMET(0), //ETmissMinBiasHF(0),
            HTmissEt(0), HTmissPhi(0), HTmissASYMHT(0), //HTmissMinBiasHF(0),
            ETHFmissEt(0), ETHFmissPhi(0), ETHFmissASYMETHF(0), //ETHFmissCENT(0),
            HTHFmissEt(0), HTHFmissPhi(0), HTHFmissASYMHTHF(0); //HTHFmissCENT(0);

    // ET
    ETEt        = ((bl->sum[0] >> demux::shiftsESums::ETEt)        & demux::masksESums::ETEt);
    ETEttem     = ((bl->sum[0] >> demux::shiftsESums::ETEttem)     & demux::masksESums::ETEttem);
    //ETMinBiasHF = ((bl->sum[0] >> demux::shiftsESums::ETMinBiasHF) & demux::masksESums::ETMinBiasHF);
    
    l1t::EtSum sum = l1t::EtSum(*dummyLVec_, l1t::EtSum::EtSumType::kTotalEt, ETEt);
    etSums->push_back(bx, sum);
    sum  = l1t::EtSum(*dummyLVec_, l1t::EtSum::EtSumType::kTotalEtEm, ETEttem);
    etSums->push_back(bx, sum);
    
    // HT
    HTEt         = ((bl->sum[1] >> demux::shiftsESums::HTEt)         & demux::masksESums::HTEt);
    HTtowerCount = ((bl->sum[1] >> demux::shiftsESums::HTtowerCount) & demux::masksESums::HTtowerCount);
    //HTMinBiasHF  = ((bl->sum[1] >> demux::shiftsESums::HTMinBiasHF)  & demux::masksESums::HTMinBiasHF);

    sum = l1t::EtSum(*dummyLVec_, l1t::EtSum::EtSumType::kTotalHt, HTEt);
    etSums->push_back(bx, sum);
    sum = l1t::EtSum(*dummyLVec_, l1t::EtSum::EtSumType::kTowerCount, HTtowerCount);
    etSums->push_back(bx, sum);

    // ETMiss
    ETmissEt        = ((bl->sum[2] >> demux::shiftsESums::ETmissEt)        & demux::masksESums::ETmissEt);
    ETmissPhi       = ((bl->sum[2] >> demux::shiftsESums::ETmissPhi)       & demux::masksESums::ETmissPhi);
    ETmissASYMET    = ((bl->sum[2] >> demux::shiftsESums::ETmissASYMET)    & demux::masksESums::ETmissASYMET);
    //ETmissMinBiasHF = ((bl->sum[2] >> demux::shiftsESums::ETmissMinBiasHF) & demux::masksESums::ETmissMinBiasHF);
    
    sum = l1t::EtSum(*dummyLVec_, l1t::EtSum::EtSumType::kMissingEt, ETmissEt, 0, ETmissPhi);
    etSums->push_back(bx, sum);
    sum = l1t::EtSum(*dummyLVec_, l1t::EtSum::EtSumType::kAsymEt, ETmissASYMET);
    etSums->push_back(bx, sum);
    

    // HTMiss
    HTmissEt        = ((bl->sum[3] >> demux::shiftsESums::HTmissEt)        & demux::masksESums::HTmissEt);
    HTmissPhi       = ((bl->sum[3] >> demux::shiftsESums::HTmissPhi)       & demux::masksESums::HTmissPhi);
    HTmissASYMHT    = ((bl->sum[3] >> demux::shiftsESums::HTmissASYMHT)    & demux::masksESums::HTmissASYMHT);
    //HTmissMinBiasHF = ((bl->sum[3] >> demux::shiftsESums::HTmissMinBiasHF) & demux::masksESums::HTmissMinBiasHF);

    sum = l1t::EtSum(*dummyLVec_, l1t::EtSum::EtSumType::kMissingHt, HTmissEt, 0, HTmissPhi);
    etSums->push_back(bx, sum);
    sum = l1t::EtSum(*dummyLVec_, l1t::EtSum::EtSumType::kAsymHt, HTmissASYMHT);
    etSums->push_back(bx, sum);

    // ETHFMiss
    ETHFmissEt       = ((bl->sum[4] >> demux::shiftsESums::ETHFmissEt)       & demux::masksESums::ETHFmissEt);
    ETHFmissPhi      = ((bl->sum[4] >> demux::shiftsESums::ETHFmissPhi)      & demux::masksESums::ETHFmissPhi);
    ETHFmissASYMETHF = ((bl->sum[4] >> demux::shiftsESums::ETHFmissASYMETHF) & demux::masksESums::ETHFmissASYMETHF);
    //ETHFmissCENT     = ((bl.sum[4] >> demux::shiftsESums::ETHFmissCENT)     & demux::masksESums::ETHFmissCENT);
    
    sum = l1t::EtSum(*dummyLVec_, l1t::EtSum::EtSumType::kMissingEtHF, ETHFmissEt, 0, ETHFmissPhi);
    etSums->push_back(bx, sum);
    sum =  l1t::EtSum(*dummyLVec_, l1t::EtSum::EtSumType::kAsymEtHF, ETHFmissASYMETHF);
    etSums->push_back(bx, sum);

   
    // HTHFMiss
    HTHFmissEt       = ((bl->sum[5] >> demux::shiftsESums::HTHFmissEt)       & demux::masksESums::HTHFmissEt);
    HTHFmissPhi      = ((bl->sum[5] >> demux::shiftsESums::HTHFmissPhi)      & demux::masksESums::HTHFmissPhi);
    HTHFmissASYMHTHF = ((bl->sum[5] >> demux::shiftsESums::HTHFmissASYMHTHF) & demux::masksESums::HTHFmissASYMHTHF);
    //HTHFmissCENT     = ((bl->sum[5] >> demux::shiftsESums::HTHFmissCENT)     & demux::masksESums::HTHFmissCENT);

    sum = l1t::EtSum(*dummyLVec_, l1t::EtSum::EtSumType::kMissingHtHF, HTHFmissEt, 0, HTHFmissPhi);
    etSums->push_back(bx, sum);
    sum =  l1t::EtSum(*dummyLVec_, l1t::EtSum::EtSumType::kAsymHtHF, HTHFmissASYMHTHF);
    etSums->push_back(bx, sum);

    // add sums to event
    // etSums->push_back(bx, bx_etSums);

  } // end of orbit loop

  jets->flatten();
  eGammas->flatten();
  taus->flatten();
  etSums->flatten();
  
}

void ScCaloRawToDigi::unpackRawJet(std::vector<l1t::Jet>& jets, uint32_t *rawData){
  
 return; 
}

void ScCaloRawToDigi::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScCaloRawToDigi);
