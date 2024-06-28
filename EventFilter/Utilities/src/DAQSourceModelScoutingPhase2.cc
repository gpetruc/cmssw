#include "EventFilter//Utilities/interface/DAQSourceModelsScoutingPhase2.h"

using namespace edm::streamer;

void DataModeScoutingPhase2::makeDirectoryEntries(std::vector<std::string> const& baseDirs,
                                                  std::vector<int> const& numSources,
                                                  std::string const& runDir) {
  std::filesystem::path runDirP(runDir);
  for (auto& baseDir : baseDirs) {
    std::filesystem::path baseDirP(baseDir);
    buPaths_.emplace_back(baseDirP / runDirP);
  }

  // store the number of sources in each BU
  buNumSources_ = numSources;
}

std::pair<bool, std::vector<std::string>> DataModeScoutingPhase2::defineAdditionalFiles(std::string const& primaryName,
                                                                                        bool fileListMode) const {
  std::vector<std::string> additionalFiles;
  if (fileListMode) {
    assert(!buNumSources_.empty());
    //for the unit test
    for (int j = 1; j < buNumSources_[0]; j++) {
      additionalFiles.push_back(primaryName + "_" + std::to_string(j));
    }
    return std::make_pair(true, additionalFiles);
  }

  auto fullpath = std::filesystem::path(primaryName);
  auto fullname = fullpath.filename();

  for (size_t i = 0; i < buPaths_.size(); i++) {
    std::filesystem::path newPath = buPaths_[i] / fullname;

    if (i != 0) {
      // secondary files from other ramdisks
      additionalFiles.push_back(newPath.generic_string());
    }

    // add extra sources from the same ramdisk
    for (int j = 1; j < buNumSources_[i]; j++) {
      additionalFiles.push_back(newPath.generic_string() + "_" + std::to_string(j));
    }
  }
  return std::make_pair(true, additionalFiles);
}

void DataModeScoutingPhase2::readEvent(edm::EventPrincipal& eventPrincipal) {
  assert(!events_.empty());

  edm::TimeValue_t time;
  timeval stv;
  gettimeofday(&stv, nullptr);
  time = stv.tv_sec;
  time = (time << 32) + stv.tv_usec;
  edm::Timestamp tstamp(time);

  // set provenance helpers
  uint32_t hdrEventID = events_.front()->event();
  edm::EventID eventID = edm::EventID(daqSource_->eventRunNumber(), daqSource_->currentLumiSection(), hdrEventID);
  edm::EventAuxiliary aux(
      eventID, daqSource_->processGUID(), tstamp, events_[0]->isRealData(), edm::EventAuxiliary::PhysicsTrigger);

  aux.setProcessHistoryID(daqSource_->processHistoryID());
  daqSource_->makeEventWrapper(eventPrincipal, aux);

  // create scouting raw data collection
  //std::cout << "Called here with " << events_.size() << " events" << std::endl;
  std::unique_ptr<SDSRawDataCollection> rawData(new SDSRawDataCollection);
  for (unsigned int id = 0, n = events_.size(); id < n; ++id) {
    const auto& e = *events_[id];
    auto size = e.eventSize();
    //std::cout << "Event " << id << " of size " << size << " payload at " << (void*)(e.payload()) << std::endl;
    auto& fedData = rawData->FEDData(id);
    fedData.resize(size);
    memcpy(fedData.data(), e.payload(), size);
  }
  std::unique_ptr<edm::WrapperBase> edp(new edm::Wrapper<SDSRawDataCollection>(std::move(rawData)));
  eventPrincipal.put(
      daqProvenanceHelpers_[0]->branchDescription(), std::move(edp), daqProvenanceHelpers_[0]->dummyProvenance());

  eventCached_ = false;
}

std::vector<std::shared_ptr<const edm::DaqProvenanceHelper>>& DataModeScoutingPhase2::makeDaqProvenanceHelpers() {
  //set FRD data collection
  daqProvenanceHelpers_.clear();
  daqProvenanceHelpers_.emplace_back(std::make_shared<const edm::DaqProvenanceHelper>(
      edm::TypeID(typeid(SDSRawDataCollection)), "SDSRawDataCollection", "SDSRawDataCollection", "DAQSource"));
  return daqProvenanceHelpers_;
}

bool DataModeScoutingPhase2::nextEventView() {
  blockCompleted_ = false;
  if (eventCached_)
    return true;
  for (unsigned int i = 0; i < events_.size(); i++) {
    //add last event length to each stripe
    dataBlockAddrs_[i] += events_[i]->size();
  }
  return makeEvents();
}

bool DataModeScoutingPhase2::makeEvents() {
  events_.clear();
  assert(!blockCompleted_);
  //std::cout << "In makeEvents(), numFiles_ = " << numFiles_ << std::endl;
  for (int i = 0; i < numFiles_; i++) {
    /*
    std::cout << " i = " << i << ", dataBlockAddrs_[i] = " << (void*)(dataBlockAddrs_[i])
              << ", dataBlockMaxAddrs_[i] = " << (void*)(dataBlockMaxAddrs_[i])
              << ", blockCompleted_ = " << blockCompleted_ << std::endl;
    */
    if (dataBlockAddrs_[i] >= dataBlockMaxAddrs_[i]) {
      //must be exact
      assert(dataBlockAddrs_[i] == dataBlockMaxAddrs_[i]);
      blockCompleted_ = true;
      return false;
    } else {
      if (blockCompleted_)
        throw cms::Exception("DataModeFRDStriped::makeEvents")
            << "not all striped blocks were completed at the same time";
    }
    if (blockCompleted_)
      continue;
    events_.emplace_back(std::make_unique<FRDEventMsgView>(dataBlockAddrs_[i]));
    /*
    std::cout << "Emplaced event " << i << " (check: " << (events_.size() - 1) << ") of size " << events_[i]->size()
              << " eventSize " << events_[i]->eventSize() << ", payload at " << (void*)(events_[i]->payload())
              << std::endl;
    */
    if (dataBlockAddrs_[i] + events_[i]->size() > dataBlockMaxAddrs_[i])
      throw cms::Exception("DAQSource::getNextEvent")
          << " event id:" << events_[i]->event() << " lumi:" << events_[i]->lumi() << " run:" << events_[i]->run()
          << " of size:" << events_[i]->size() << " bytes does not fit into the buffer or has corrupted header";
  }
  return !blockCompleted_;
}

bool DataModeScoutingPhase2::checksumValid() { return true; }

std::string DataModeScoutingPhase2::getChecksumError() const { return std::string(); }
