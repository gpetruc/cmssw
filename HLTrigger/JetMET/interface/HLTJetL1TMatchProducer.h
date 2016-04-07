#ifndef HLTJetL1TMatchProducer_h
#define HLTJetL1TMatchProducer_h

#include <string>
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/L1Trigger/interface/L1JetParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/Jet.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/TrackJetCollection.h"
#include "DataFormats/JetReco/interface/BasicJetCollection.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"

template<typename T>
class HLTJetL1TMatchProducer : public edm::stream::EDProducer<> {
 public:
  explicit HLTJetL1TMatchProducer(const edm::ParameterSet&);
  ~HLTJetL1TMatchProducer();
  static  void fillDescriptions(edm::ConfigurationDescriptions & descriptions);
  virtual void beginJob() ; 
  virtual void produce(edm::Event &, const edm::EventSetup&);
 private:
  edm::EDGetTokenT<std::vector<T>> m_theJetToken;
  edm::EDGetTokenT<l1t::JetBxCollection> m_theL1JetToken;
  edm::EDGetTokenT<trigger::TriggerFilterObjectWithRefs> m_theL1FilterToken;
  edm::InputTag jetsInput_;
  edm::InputTag L1Jets_;
  edm::InputTag L1Filter_;
  //  std::string jetType_;
  double DeltaR_;         // DeltaR(HLT,L1)
  bool useL1Filter_;
};

#endif
