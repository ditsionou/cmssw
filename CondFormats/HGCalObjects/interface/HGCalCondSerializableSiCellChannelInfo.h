#ifndef CondFormats_HGCalObjects_HGCalCondSerializableSiCellChannelInfo_h
#define CondFormats_HGCalObjects_HGCalCondSerializableSiCellChannelInfo_h

#include <vector>
#include <map>

#include "CondFormats/Serialization/interface/Serializable.h"

/**
   @short cond. serializable class holding information on a Si sensor cell
 */
class HGCalCondSerializableSiCellChannelInfo {
public:

  struct HGCalSiCellChannelInfo {
    bool isHD,iscalib;
    uint8_t wafType, chip, half;
    uint16_t seq,rocpin;
    int sicell,iu,iv,t;
    float trace;
    COND_SERIALIZABLE;
  };
  
  HGCalCondSerializableSiCellChannelInfo() {}
  virtual ~HGCalCondSerializableSiCellChannelInfo() {}
  HGCalCondSerializableSiCellChannelInfo& addParameter(HGCalSiCellChannelInfo&);
  std::vector<HGCalSiCellChannelInfo> params_;

  COND_SERIALIZABLE;
};

#endif
