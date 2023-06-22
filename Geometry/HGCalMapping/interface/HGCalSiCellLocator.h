#ifndef _geometry_hgcalmapping_hgcalsicelllocator_h_
#define _geometry_hgcalmapping_hgcalsicelllocator_h_

#include <string>
#include <vector>
#include "CondFormats/HGCalObjects/interface/HGCalCondSerializableSiCellChannelInfo.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetId.h"
#include "DataFormats/HGCalDigi/interface/HGCalElectronicsId.h"

class HGCalSiCellLocator
{

public:

  HGCalSiCellLocator();
  void buildLocatorFrom(std::string url,bool append=false,bool usefip=false);
  HGCalSiCellChannelInfo locateCellByGeom(int iu,int iv,uint8_t wafType, bool isHD);
  DetId getDetId(HGCalElectronicsId& id, int z, int layer, int modU, int modZ) const;  
  ~HGCalSiCellLocator();

private:

  HGCalCondSerializableSiCellChannelInfo cellColl_;
};



#endif
