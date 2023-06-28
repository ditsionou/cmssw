#ifndef _geometry_hgcalmapping_hgcalsicelllocator_h_
#define _geometry_hgcalmapping_hgcalsicelllocator_h_

#include <string>
#include <vector>
#include "CondFormats/HGCalObjects/interface/HGCalCondSerializableSiCellChannelInfo.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetId.h"
#include "DataFormats/HGCalDigi/interface/HGCalElectronicsId.h"

class HGCalSiCellLocator {
public:
  HGCalSiCellLocator();
  void buildLocatorFrom(std::string url, bool append = false, bool usefip = false);
  HGCalSiCellChannelInfo locateCellByGeom(int iu, int iv, uint8_t wafType, bool isHD);
  HGCalSiCellChannelInfo locateCellByChannel(uint8_t roc, uint8_t rocHalf, uint16_t rocPin, uint8_t wafType, bool isHD);
  ~HGCalSiCellLocator();

private:
  HGCalCondSerializableSiCellChannelInfo cellColl_;
};

#endif
