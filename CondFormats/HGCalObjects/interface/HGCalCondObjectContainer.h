#ifndef _condformats_hgcalobjects_hgcalcondobjectcontainer_h_
#define _condformats_hgcalobjects_hgcalcondobjectcontainer_h_

#include "CondFormats/Serialization/interface/Serializable.h"
#include <vector>

/**
   @class HGCalCondObjectContainer
   @short it defines a templated container indexed by a uint32_t
   (the uint32_t is the raw value of HGCal DetId or ElectronicsId)
   Templated items can be added, modified and retrieved from the collection.
   The container is serializable as cond format
 */
template <typename T>
class HGCalCondObjectContainer {
public:

  HGCalCondObjectContainer(){};
  ~HGCalCondObjectContainer(){};

  inline void clear() { objColl_.clear(); }
  inline const std::vector<T> &items() const { return objColl_.items(); }
  inline const T &item(size_t hashedIndex) const { return objColl_.item(hashedIndex); }
  inline void insert(std::pair<uint32_t, T> const &a) { objColl_.insert(a); }
  inline typename std::vector<T>::const_iterator find(uint32_t rawId) const { return objColl_.find(rawId); }
  inline typename std::vector<T>::const_iterator end() const { return objColl_.end(); }
  inline void setValue(const uint32_t id, const T &item) { (*this)[id] = item; }
  inline const HGCalCondObjectContainer<T> &getMap() const { return *this; }
  inline size_t size() const { return objColl_.size(); }
  inline T &operator[](uint32_t rawId) { return objColl_[rawId]; }
  inline T operator[](uint32_t rawId) const { return objColl_[rawId]; }

private:

  std::map<uint32_t,T> objColl_;

  COND_SERIALIZABLE;
};

typedef HGCalCondObjectContainer<uint32_t> HGCalCondUInt32Container;

#endif
