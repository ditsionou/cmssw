#ifndef _geometry_hgcalmapping_hgcalsicelllocator_h_
#define _geometry_hgcalmapping_hgcalsicelllocator_h_

#include <string>
#include <vector>

class HGCalSiCellLocator
{
  struct HGCalSiCellChannel{
    bool isHD,iscalib;
    uint8_t wafType, chip, half;
    uint16_t seq,rocpin;
    int sicell,iu,iv,t;
  };

public:

  HGCalSiCellLocator();
  void buildLocatorFrom(std::string url,bool append=false);
  HGCalSiCellChannel locateCellByGeom(int iu,int iv,uint8_t wafType, bool isHD);
  ~HGCalSiCellLocator();

private:

  std::vector<HGCalSiCellChannel> cellColl_;

};



#endif
