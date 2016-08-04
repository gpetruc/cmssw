/**
  \class    InstLumiAsFloatProducer 
  \brief    Copy instantaneous luminosity from the LumiScalers to a float
  \Author   Giovanni Petrucciani (from code written by Carlo Battilana for the T&P)
*/


#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Scalers/interface/LumiScalers.h"
#include "DataFormats/Luminosity/interface/LumiDetails.h"

class InstLumiAsFloatProducer : public edm::stream::EDProducer<>
{
  
public:
    explicit InstLumiAsFloatProducer(const edm::ParameterSet & iConfig);
    virtual ~InstLumiAsFloatProducer() { }

    virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup) final override;
private:
    edm::EDGetTokenT<LumiScalersCollection> m_lumiScalerTag;

};

InstLumiAsFloatProducer::InstLumiAsFloatProducer(const edm::ParameterSet & iConfig) :
  m_lumiScalerTag(consumes<LumiScalersCollection>(iConfig.getParameter<edm::InputTag>("lumiScalers")))
{
    produces<float>("instantLumi");
}

void InstLumiAsFloatProducer::produce(edm::Event & ev, const edm::EventSetup & iSetup)
{
    std::auto_ptr<float> instLumi(new float(-999.));

    if (ev.isRealData()) {
        edm::Handle<LumiScalersCollection> lumiScaler;
        ev.getByToken(m_lumiScalerTag, lumiScaler);
        if (!lumiScaler->empty()) {
            *instLumi= lumiScaler->front().instantLumi();
        }
    }

    ev.put(instLumi, "instantLumi");
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(InstLumiAsFloatProducer);
