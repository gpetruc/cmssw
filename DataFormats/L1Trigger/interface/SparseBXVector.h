#ifndef DataFormats_L1Trigger_SparseBXVector_h
#define DataFormats_L1Trigger_SparseBXVector_h

// this class is an extension of std::vector
// designed to store objects corresponding to several time-samples (BX)
// the lists of BXs and objects in each BX are flattened into linear vectors
// BXs should be added in order

#include "DataFormats/Common/interface/FillView.h"
#include "DataFormats/Common/interface/fillPtrVector.h"
#include "DataFormats/Common/interface/setPtr.h"
#include "DataFormats/Common/interface/traits.h"
#include <vector>
#include <algorithm>

template <class T>
class SparseBXVector {
public:
  typedef typename std::vector<T>::iterator iterator;
  typedef typename std::vector<T>::const_iterator const_iterator;
  typedef typename std::vector<int>::const_iterator bx_const_iterator;
  typedef T value_type;
  typedef typename std::vector<T>::size_type size_type;

  class BxAndRangeIterator {
  public:
    BxAndRangeIterator(const SparseBXVector<T>* src = nullptr, unsigned int i = 0)
        : src_(src), i_(i), n_(src != nullptr ? src->nbx() : 0) {}
    bool operator==(const BxAndRangeIterator& other) const { return src_ == other.src_ && i_ == other.i_; }
    bool operator!=(const BxAndRangeIterator& other) const { return src_ != other.src_ || i_ != other.i_; }
    BxAndRangeIterator& operator++() {
      ++i_;
      return *this;
    }
    BxAndRangeIterator& operator--() {
      --i_;
      return *this;
    }
    unsigned int ibx() const { return i_; }
    int bx() const { return src_->bxs_[i_]; }
    const T * data() const { return &src_->data_[src_->itrs_[i_]]; }
    const_iterator begin() const { return src_->data_.begin() + src_->itrs_[i_]; }
    const_iterator end() const { return i_ < n_ - 1 ? src_->data_.begin() + src_->itrs_[i_ + 1] : src_->data_.end(); }
    unsigned int size() const { return (i_ < n_ - 1 ? src_->itrs_[i_ + 1] : src_->data_.size()) - src_->itrs_[i_]; }
    const T& operator[](unsigned int i) const { return *(begin() + i); }
    bool good() const { return i_ < n_; }

  private:
    const SparseBXVector<T>* src_;
    unsigned int i_, n_;
  };

  // default ctor
  SparseBXVector() {}

  // default dtor
  ~SparseBXVector() {}

  // copy ctor
  SparseBXVector(const SparseBXVector<T> & other) = default;
  // copy assign
  SparseBXVector<T> & operator=(const SparseBXVector<T> & other) = default;
  // move ctor
  SparseBXVector(SparseBXVector<T> && other) = default;
  // move assign
  SparseBXVector<T> & operator=(SparseBXVector<T> && other) = default;

  // swap
  void swap(SparseBXVector<T> & other) {
    using std::swap;
    swap(bxs_, other.bxs_);
    swap(data_, other.data_);
    swap(itrs_, other.itrs_);
  }

  // add one BX to end of a SparseBXVector
  template <typename VI>
  void addBX(int bx, VI objectsBegin, VI objectsEnd) {
    assert(bxs_.empty() || bx > bxs_.back());
    bxs_.push_back(bx);
    itrs_.push_back(data_.size());
    data_.insert(data_.end(), objectsBegin, objectsEnd);
  }

  // iterate on BX and data ranges (preferred option)
  BxAndRangeIterator iter() const { return BxAndRangeIterator(this); }

  // iterator access for BX list
  bx_const_iterator beginBXs() const { return bxs_.begin(); }

  // iterator access for BX list
  bx_const_iterator endBXs() const { return bxs_.end(); }

  // iterator access by BX
  const_iterator begin(int bx) const { return range(bx).first; }

  // iterator access by BX
  const_iterator end(int bx) const { return range(bx).second; }

  std::pair<const_iterator, const_iterator> range(int bx) const {
    auto bxmatch = std::lower_bound(bxs_.begin(), bxs_.end(), bx);
    if ((bxmatch != bxs_.end()) && (*bxmatch == bx)) {
      auto index = std::distance(bxs_.begin(), bxmatch);
      auto tail = (index + 1 < itrs_.size() ? itrs_[index + 1] : data_.size());
      return std::make_pair(begin() + itrs_[index], begin() + next);
    } else {
      return std::make_pair(end(), end());
    }
  }

  // check whether the specified BX is available
  bool hasBX(int bx) const { return std::binary_search(bxs_.begin(), bxs_.end(), bx); }

  // get N objects for a given BX
  unsigned size(int bx) const {
    auto r = range(bx);
    return std::distance(r.first, r.second);
  }

  // get N BXs 
  unsigned nbx() const { return bxs_.size(); }

  // get N objects for all BXs together
  unsigned size() const { return data_.size(); }

  // add element with given BX index
  // slow unless this is the last BX in the vector
  void push_back(int bx, const T& object) {
    assert(bxs_.empty() || bxs_.back() <= bx);
    if (bxs_.empty() || bxs_.back() != bx) {
      bxs_.push_back(bx);
      itrs_.push_back(data_.size());
      data_.push_back(object);
    } else if (bxs_.back() == bx) {
      data_.push_back(object);
    }
  }

  // clear entire BXVector
  void clear() {
    bxs_.clear();
    data_.clear();
    itrs_.clear();
  }

  // check if it's empty
  bool empty() const { return bxs_.empty(); }

  // check if data has empty location
  bool isEmpty(int bx) const {
    auto r = range(bx);
    return r.first == r.second;
  }

  // support looping over entire collection (note also that begin() is needed by edm::Ref)
  const_iterator begin() const { return data_.begin(); }
  const_iterator end() const { return data_.end(); }
  //int bx(const_iterator & iter) const; (potentially useful)
  unsigned int key(const_iterator& iter) const { return iter - begin(); }

  // array subscript operator (incited by TriggerSummaryProducerAOD::fillTriggerObject...)
  T& operator[](std::size_t i) { return data_[i]; }
  const T& operator[](std::size_t i) const { return data_[i]; }

  // edm::View support
  void fillView(edm::ProductID const& id,
                std::vector<void const*>& pointers,
                edm::FillViewHelperVector& helpers) const {
    edm::detail::reallyFillView(*this, id, pointers, helpers);
  }
  // edm::Ptr support
  void setPtr(std::type_info const& toType, unsigned long index, void const*& ptr) const {
    edm::detail::reallySetPtr<SparseBXVector<T> >(*this, toType, index, ptr);
  }
  void fillPtrVector(std::type_info const& toType,
                     std::vector<unsigned long> const& indices,
                     std::vector<void const*>& ptrs) const {
    edm::detail::reallyfillPtrVector(*this, toType, indices, ptrs);
  }

private:
  std::vector<int> bxs_;
  std::vector<T> data_;
  std::vector<unsigned> itrs_;
};

// edm::View support
namespace edm {
  template <class T>
  inline void fillView(SparseBXVector<T> const& obj,
                       edm::ProductID const& id,
                       std::vector<void const*>& pointers,
                       edm::FillViewHelperVector& helpers) {
    obj.fillView(id, pointers, helpers);
  }
  template <class T>
  struct has_fillView<SparseBXVector<T> > {
    static bool const value = true;
  };
}  // namespace edm
// edm::Ptr support
template <class T>
inline void setPtr(SparseBXVector<T> const& obj, std::type_info const& toType, unsigned long index, void const*& ptr) {
  obj.setPtr(toType, index, ptr);
}
template <class T>
inline void fillPtrVector(SparseBXVector<T> const& obj,
                          std::type_info const& toType,
                          std::vector<unsigned long> const& indices,
                          std::vector<void const*>& ptrs) {
  obj.fillPtrVector(toType, indices, ptrs);
}
namespace edm {
  template <class T>
  struct has_setPtr<SparseBXVector<T> > {
    static bool const value = true;
  };
}  // namespace edm
#endif
