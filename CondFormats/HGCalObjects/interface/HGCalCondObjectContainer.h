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
template <typename K,typename V>
class HGCalCondObjectContainer {
public:

  HGCalCondObjectContainer(){};
  ~HGCalCondObjectContainer(){};

  inline void clear() { objColl_.clear(); }
  inline const std::vector<V> &items() const { return objColl_.items(); }
  inline const V &item(size_t hashedIndex) const { return objColl_.item(hashedIndex); }
  inline void insert(std::pair<K, V> const &a) { objColl_.insert(a); }
  inline typename std::vector<V>::const_iterator find(K rawId) const { return objColl_.find(rawId); }
  inline typename std::vector<V>::const_iterator end() const { return objColl_.end(); }
  inline void setValue(const K id, const V &item) { (*this)[id] = item; }
  inline const HGCalCondObjectContainer<K,V> &getMap() const { return *this; }
  inline size_t size() const { return objColl_.size(); }
  inline V &operator[](K rawId) { return objColl_[rawId]; }
  inline V operator[](K rawId) const { return objColl_[rawId]; }

private:

  std::map<K,V> objColl_;

  COND_SERIALIZABLE;
};

#endif
