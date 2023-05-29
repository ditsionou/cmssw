#ifndef CondFormats_HGCalObjects_HGCalCondSerializableModuleInfo_h
#define CondFormats_HGCalObjects_HGCalCondSerializableModuleInfo_h

#include <vector>
#include <cstdint>

#include "CondFormats/Serialization/interface/Serializable.h"
/**
   @short representation of a si-cell channel information read from a txt file or db
*/
struct HGCalModuleInfo {
  bool zside,isSiPM,isHD;
  int plane, u, v;
  std::string modtype;
  int econdidx, captureblock, slink, captureblockidx, fedid;
  std::string DAQ;  
  COND_SERIALIZABLE;
};

/**
   @holder for the si cell channel mappeing
 */
class HGCalCondSerializableModuleInfo {
  
public:
  
  HGCalCondSerializableModuleInfo() {}
  virtual ~HGCalCondSerializableModuleInfo() {}
  inline HGCalCondSerializableModuleInfo& addParameter(HGCalModuleInfo &i) {
    params_.push_back(i);
    return *this;
  }
  std::vector<HGCalModuleInfo> params_;
  
  COND_SERIALIZABLE;
};

#endif
