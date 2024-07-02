#ifndef L1TriggerScouting_Utilities_OrbitTableOutputBranches_h
#define L1TriggerScouting_Utilities_OrbitTableOutputBranches_h

#include <string>
#include <vector>
#include <TTree.h>
#include "FWCore/Framework/interface/OccurrenceForOutput.h"
#include "DataFormats/L1Scouting/interface/OrbitFlatTable.h"
#include "DataFormats/Provenance/interface/BranchDescription.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

class OrbitTableOutputBranches {
public:
  OrbitTableOutputBranches(const edm::BranchDescription *desc, const edm::EDGetToken &token)
      : m_token(token), m_extension(DontKnowYetIfMainOrExtension), m_branchesBooked(false) {
    if (desc->className() != "l1ScoutingRun3::OrbitFlatTable")
      throw cms::Exception("Configuration",
                           "NanoAODOutputModule can only write out l1ScoutingRun3::OrbitFlatTable objects");
  }

  void defineBranchesFromFirstEvent(const l1ScoutingRun3::OrbitFlatTable &tab);
  void branch(TTree &tree);

  /// Fill the current table, if extensions == table.extension().
  /// This parameter is used so that the fill is called first for non-extensions and then for extensions
  void beginFill(const edm::OccurrenceForOutput &iWhatever, TTree &tree, bool extensions);
  bool hasBx(uint32_t bx);
  void fillBx(uint32_t bx, bool skipReadingSize=false);
  void endFill();

private:
  edm::EDGetToken m_token;
  std::string m_baseName;
  bool m_singleton = false;
  enum { IsMain = 0, IsExtension = 1, DontKnowYetIfMainOrExtension = 2 } m_extension;
  std::string m_doc;
  typedef Int_t CounterType;
  CounterType m_counter;
  struct NamedBranchPtr {
    std::string name, title, rootTypeCode;
    int columIndex;
    TBranch *branch;
    NamedBranchPtr(const std::string &aname,
                   const std::string &atitle,
                   const std::string &rootType,
                   TBranch *branchptr = nullptr)
        : name(aname), title(atitle), rootTypeCode(rootType), columIndex(-1), branch(branchptr) {}
    void getIndex(const l1ScoutingRun3::OrbitFlatTable &tab, const std::string &baseName) {
      columIndex = tab.columnIndex(name);
      if (columIndex == -1)
        throw cms::Exception("LogicError", "Missing column in input for " + baseName + "_" + name);
    }
  };
  TBranch *m_counterBranch = nullptr;
  std::vector<NamedBranchPtr> m_uint8Branches;
  std::vector<NamedBranchPtr> m_int16Branches;
  std::vector<NamedBranchPtr> m_uint16Branches;
  std::vector<NamedBranchPtr> m_int32Branches;
  std::vector<NamedBranchPtr> m_uint32Branches;
  std::vector<NamedBranchPtr> m_floatBranches;
  std::vector<NamedBranchPtr> m_doubleBranches;
  bool m_branchesBooked;

  edm::Handle<l1ScoutingRun3::OrbitFlatTable> m_handle;
  const l1ScoutingRun3::OrbitFlatTable *m_table;

  template <typename T>
  void fillColumn(NamedBranchPtr &pair, uint32_t bx, bool empty) {
    pair.branch->SetAddress(
        empty ? static_cast<T *>(nullptr)
              : const_cast<T *>(
                    &m_table->columnData<T>(pair.columIndex, bx).front()));  // SetAddress should take a const * !
  }
};

#endif
