#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDFilter.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/EDPutToken.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "DataFormats/L1ScoutingSoA/interface/CounterHost.h"

#include <vector>
#include <set>

/*
 * Filter orbits that don't contain at least one selected BX
 * from a BxSelector module and produce a vector of selected BXs
 */
class GoodOrbitNBxSelector : public edm::global::EDFilter<> {
public:
  explicit GoodOrbitNBxSelector(const edm::ParameterSet&);
  ~GoodOrbitNBxSelector() override {}
  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  bool filter(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;
  void endJob() override;

  // tokens for BX selected by each analysis
  std::vector<edm::EDGetTokenT<unsigned int>> nbxTokens_;
  std::vector<edm::EDGetTokenT<l1sc::CounterHost>> nbxPortableTokens_;
  unsigned int threshold_;
  unsigned long long nPrint_;
  mutable std::atomic<unsigned long long> goodOrbits_, badOrbits_, events_;
};

GoodOrbitNBxSelector::GoodOrbitNBxSelector(const edm::ParameterSet& iPSet)
    : threshold_(iPSet.getParameter<unsigned int>("nbxMin")),
      nPrint_(iPSet.getUntrackedParameter<unsigned int>("nPrint")),
      goodOrbits_(0),
      badOrbits_(0),
      events_(0) {
  // get the list of selected BXs
  std::vector<edm::InputTag> bxLabels = iPSet.getParameter<std::vector<edm::InputTag>>("unpackers");
  for (const auto& bxLabel : bxLabels) {
    edm::InputTag tag{bxLabel.label(), "nbx", bxLabel.process()};
    nbxTokens_.push_back(consumes<unsigned int>(tag));
  }
  bxLabels = iPSet.getParameter<std::vector<edm::InputTag>>("unpackersAlpaka");
  for (const auto& bxLabel : bxLabels) {
    edm::InputTag tag{bxLabel.label(), "nbx", bxLabel.process()};
    nbxPortableTokens_.push_back(consumes<l1sc::CounterHost>(tag));
  }
}

// ------------ method called for each ORBIT  ------------
bool GoodOrbitNBxSelector::filter(edm::StreamID, edm::Event& iEvent, const edm::EventSetup&) const {
  unsigned nbxMin = threshold_;
  for (const auto& token : nbxTokens_) {
    edm::Handle<unsigned> nbx;
    iEvent.getByToken(token, nbx);
    if (*nbx < threshold_) {
      ++badOrbits_;
      return false;
    }
    nbxMin = std::min(*nbx, nbxMin);
  }
  for (const auto& token : nbxPortableTokens_) {
    edm::Handle<l1sc::CounterHost> nbx;
    iEvent.getByToken(token, nbx);
    if (**nbx < threshold_) {
      ++badOrbits_;
      return false;
    }
    nbxMin = std::min(**nbx, nbxMin);
  }
  events_ += nbxMin;
  if ((goodOrbits_++) % nPrint_ == 0) {
    edm::LogImportant("GoodOrbitNBxSelector")
        << "Processed " << (goodOrbits_.load() + badOrbits_.load()) << " orbits, of which " << badOrbits_.load()
        << " truncated, and " << events_.load() << " events.\n";
  }
  return true;
}

void GoodOrbitNBxSelector::endJob() {
  edm::LogImportant("GoodOrbitNBxSelector")
      << "Processed " << (goodOrbits_.load() + badOrbits_.load()) << " orbits, of which " << badOrbits_.load()
      << " truncated, and " << events_.load() << " events.\n";
}

void GoodOrbitNBxSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::vector<edm::InputTag>>("unpackers", std::vector<edm::InputTag>{});
  desc.add<std::vector<edm::InputTag>>("unpackersAlpaka", std::vector<edm::InputTag>{});
  desc.add<unsigned int>("nbxMin", 3564);           // BXs in one orbit
  desc.addUntracked<unsigned int>("nPrint", 1000);  // Number of orbits between printouts
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(GoodOrbitNBxSelector);
