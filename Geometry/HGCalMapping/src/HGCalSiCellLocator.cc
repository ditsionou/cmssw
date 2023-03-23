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
void HGCalSiCellLocator::buildLocatorFrom(std::string url,bool append) {

  if(!append) cellColl_.clear();

  //open file and parse each line
  edm::FileInPath fip(url);
  std::ifstream inF(fip.fullPath());
  std::string line;
  while (std::getline(inF, line)) {
    
    HGCalSiCellChannel c;
          
    std::istringstream strm(line);

    //density and rocpin  need to be decoded from the string
    //other fields are read directly
    std::string denscol,rocpincol;
    strm >> denscol;
    c.isHD = denscol=="LD" ? false : true;
    strm >> c.wafType >> c.chip >> c.half >> c.seq;
    strm >> rocpincol;
    if(rocpincol.find("CALIB")!=std::string::npos) {
      c.iscalib=true;
      c.rocpin=uint16_t(rocpincol[rocpincol.size()-1]);
    }
    else {
      c.iscalib=false;
      c.rocpin=std::stoi(rocpincol);
    }
    strm >> c.sicell >> c.iu >> c.iv >> c.t;
    cellColl_.push_back(c);
  }

}

//
HGCalSiCellLocator::HGCalSiCellChannel HGCalSiCellLocator::locateCellByGeom(int iu,int iv,uint8_t wafType, bool isHD) {

  auto _matchesByGeom = [iu,iv,wafType,isHD](HGCalSiCellChannel c){ return c.iu==iu && c.iv==iv && c.wafType==wafType && c.isHD==isHD; };
  auto it = std::find_if(begin(cellColl_), end(cellColl_), _matchesByGeom);
  if(it==cellColl_.end()) {
    edm::Exception e(edm::errors::NotFound,"Failed to match Si cell to channel by geometry");
    throw e;
  }
  return *it;
  
}


//
HGCalSiCellLocator::~HGCalSiCellLocator() {
}
