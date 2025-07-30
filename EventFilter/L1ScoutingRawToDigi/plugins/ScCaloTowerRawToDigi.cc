#include "EventFilter/L1ScoutingRawToDigi/plugins/ScCaloTowerRawToDigi.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

ScCaloTowerRawToDigi::ScCaloTowerRawToDigi(const edm::ParameterSet& iConfig)
    : orbitBuffer_(NBX + 1),
      nCaloTowersOrbit_(0),
      sourceIdList_(iConfig.getParameter<std::vector<int>>("sourceIdList")),
      srcInputTag_(iConfig.getParameter<edm::InputTag>("srcInputTag")),
      rawToken_(consumes<SDSRawDataCollection>(srcInputTag_)),
      debug_(iConfig.getUntrackedParameter<bool>("debug", false)) {
  for (auto& bxVec : orbitBuffer_) {
    bxVec.reserve(4096);  // reasonable upper estimate
  }
  for (const auto& sdsId : sourceIdList_) {
    if ((sdsId < SDSNumbering::CaloTowerMinSDSID) || (sdsId > SDSNumbering::CaloTowerMaxSDSID))
      edm::LogError("ScCaloTowerRawToDigi")
          << "Provided a source ID outside the expected range: " << sdsId << ", expected range ["
          << SDSNumbering::CaloTowerMinSDSID << ", " << SDSNumbering::CaloTowerMaxSDSID;
  }
  produces<l1ScoutingRun3::CaloTowerOrbitCollection>("CaloTower").setBranchAlias("CaloTowerOrbitCollection");
}

ScCaloTowerRawToDigi::~ScCaloTowerRawToDigi() {}

void ScCaloTowerRawToDigi::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<SDSRawDataCollection> ScoutingRawDataCollection;
  iEvent.getByToken(rawToken_, ScoutingRawDataCollection);

  std::unique_ptr<l1ScoutingRun3::CaloTowerOrbitCollection> unpackedCaloTowers(
      new l1ScoutingRun3::CaloTowerOrbitCollection);

  nCaloTowersOrbit_ = 0;

  for (const auto& sdsId : sourceIdList_) {
    const FEDRawData& sourceRawData = ScoutingRawDataCollection->FEDData(sdsId);
    size_t orbitSize = sourceRawData.size();

    if (debug_ && (orbitSize == 0)) {
      edm::LogWarning("ScCaloTowerRawToDigi::produce") << "No raw data for CaloTower FED " << sdsId << std::endl;
    }

    // unpack current orbit and store data into the orbitBufferr
    unpackOrbit(sourceRawData.data(), orbitSize, sdsId);
  }

  // fill orbit collection and clear the Bx buffer vector
  unpackedCaloTowers->fillAndClear(orbitBuffer_, nCaloTowersOrbit_);

  // store collection in the event
  iEvent.put(std::move(unpackedCaloTowers), "CaloTower");
}

void ScCaloTowerRawToDigi::unpackOrbit(const unsigned char* buf, size_t len, int sdsId) {
  using namespace l1ScoutingRun3;

  size_t pos = 0;
  const size_t blockHeaderSize = 3 * sizeof(uint32_t);

  while (pos < len) {
    if (pos + blockHeaderSize > len) {
      edm::LogError("ScCaloTowerRawToDigi") << "Corrupt data in sourceId " << sdsId << ", incomplete header";
      break;  // no sense trying to unpack further
    }
    const calol2::block* bl = reinterpret_cast<const calol2::block*>(buf + pos);

    unsigned bx = bl->bx;
    unsigned orbit = (bl->orbit) & 0x7FFFFFFF;
    unsigned ctCount = bl->header;

    pos += blockHeaderSize;

    if (pos + 4 * ctCount > len) {
      edm::LogError("ScCaloTowerRawToDigi")
          << "Corrupt data in sourceId " << sdsId << ", orbit " << orbit << ", BX " << bx << ": expecting " << ctCount
          << " towers but only " << (len - pos) << " bytes left in the block.";
      break;  // no sense trying to unpack further
    }
    if (bx > NBX) {  // need this check as otherwise the code will crash later accessing orbitBuffer_
      edm::LogError("ScCaloTowerRawToDigi")
          << "Corrupt data in sourceId " << sdsId << ", orbit " << orbit << ", invalid BX " << bx;
      break;  // we could go to the next block, but if the data is corrupted it probably doesn't help
    }
    if (debug_) {
      edm::LogPrint("ScCaloTowerRawToDigi") << " CaloTower #" << sdsId << " Orbit " << orbit << ", BX -> " << bx
                                            << ", nCaloTowers -> " << ctCount << std::endl;
    }

    // Unpack calo towers
    auto& bufferThisBX = orbitBuffer_[bx];
    int32_t ET, erBits, miscBits, eta, phi;
    const uint32_t* towerPtr = reinterpret_cast<const uint32_t*>(buf + pos);
    for (unsigned int i = 0; i < ctCount; ++i, ++towerPtr) {
      uint32_t ct_raw = *towerPtr;

      ET = ((ct_raw >> calol2::shiftsCaloTowers::ET) & calol2::masksCaloTowers::ET);
      erBits = ((ct_raw >> calol2::shiftsCaloTowers::erBits) & calol2::masksCaloTowers::erBits);
      miscBits = ((ct_raw >> calol2::shiftsCaloTowers::miscBits) & calol2::masksCaloTowers::miscBits);
      eta = ((ct_raw >> calol2::shiftsCaloTowers::eta) & calol2::masksCaloTowers::eta);
      phi = ((ct_raw >> calol2::shiftsCaloTowers::phi) & calol2::masksCaloTowers::phi);

      eta = eta >= 128 ? eta - 256 : eta;

      bufferThisBX.emplace_back(ET, erBits, miscBits, eta, phi);

      if (debug_) {
        edm::LogPrint("LogPrint") << "Calo Tower " << i << ", raw: 0x" << std::hex << ct_raw << std::dec
                                  << "\n\tET: " << ET << "\n\tER bits: " << erBits << "\n\tMisc bits: " << miscBits
                                  << "\n\tEta: " << eta << "\n\tPhi: " << phi;
      }
    }
    pos += 4 * ctCount;
    nCaloTowersOrbit_ += ctCount;

  }  // end orbit while loop
}

void ScCaloTowerRawToDigi::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("srcInputTag", edm::InputTag("rawDataCollector"));

  std::vector<int> sourceIds;
  for (int id = SDSNumbering::CaloTowerMinSDSID; id <= SDSNumbering::CaloTowerMaxSDSID; ++id)
    sourceIds.emplace_back(id);
  desc.add<std::vector<int>>("sourceIdList", sourceIds);

  desc.addUntracked<bool>("debug", false);

  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScCaloTowerRawToDigi);