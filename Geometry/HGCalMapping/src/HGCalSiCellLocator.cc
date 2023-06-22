#include "Geometry/HGCalMapping/interface/HGCalSiCellLocator.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

//
HGCalSiCellLocator::HGCalSiCellLocator() {
}

//
void HGCalSiCellLocator::buildLocatorFrom(std::string url,bool append,bool usefip) {

  if(!append) cellColl_.params_.clear();
  //open file and parse each line
  if(usefip){
    edm::FileInPath fip(url);
    url=fip.fullPath();
  }
  std::ifstream inF(url);
  std::string line;
  //skip the first line
  std::getline(inF, line);
  while (std::getline(inF, line)) {
    std::istringstream strm(line);
    //density and rocpin  need to be decoded from the string, other fields are read directly
    HGCalSiCellChannelInfo c;
    std::string denscol,rocpincol;
    strm >> denscol;
    c.isHD = denscol=="LD" ? false : true;
    std::string wafType;
    strm >> wafType >> c.chip >> c.half >> c.seq;
    c.wafType = std::stoi(wafType);
    strm >> rocpincol;
    if(rocpincol.find("CALIB")!=std::string::npos) {
      c.iscalib=true;
      c.rocpin=uint16_t(rocpincol[rocpincol.size()-1]);
    }
    else {
      c.iscalib=false;
      c.rocpin=std::stoi(rocpincol);
    }
    strm >> c.sicell >> c.triglink >> c.trigcell >> c.iu >> c.iv >> c.trace >> c.t;
    cellColl_.addParameter(c);
  }

}

//
HGCalSiCellChannelInfo HGCalSiCellLocator::locateCellByGeom(int iu,int iv,uint8_t wafType, bool isHD) {
  auto _matchesByGeom = [iu,iv,wafType,isHD](HGCalSiCellChannelInfo c){ return c.iu == iu && c.iv == iv; };
  auto it = std::find_if(begin(cellColl_.params_), end(cellColl_.params_), _matchesByGeom);
  if(it==cellColl_.params_.end()) {
    std::cout << "iu: " << iu << ", iv: " << iv << ", wafType: " << int(wafType) << ", isHD: " << isHD << std::endl; 
    edm::Exception e(edm::errors::NotFound,"Failed to match Si cell to channel by geometry");
    throw e;
  }
  return *it;
}

DetId HGCalSiCellLocator::getDetId(HGCalElectronicsId& id, int z, int layer, int modU, int modV) const
{
  int chip = (int)(id.econdeRx()/2.);
  int half = id.econdeRx()%2;
  int rocpin = id.halfrocChannel();
  int seq = id.sequentialHalfrocChannel();
  std::cout << "chip: " << chip << ", half: " << half << ", rocpin: " << rocpin << ", seq: " << seq << std::endl; 
  auto _matchesByChannel = [rocpin, seq, chip, half](HGCalSiCellChannelInfo c){
    return c.rocpin == rocpin && c.seq == seq && c.chip == chip && c.half == half;
  };
  auto it = std::find_if(begin(cellColl_.params_), end(cellColl_.params_), _matchesByChannel);
  if(it==cellColl_.params_.end()) {
    edm::Exception e(edm::errors::NotFound,"Failed to match Si cell to geometry by channel");
    throw e;
  }

  HGCSiliconDetId detId(DetId::HGCalHSc, z, it->wafType, layer, modU, modV, it->iu, it->iv);
  return detId;
}

//
HGCalSiCellLocator::~HGCalSiCellLocator() {
}
