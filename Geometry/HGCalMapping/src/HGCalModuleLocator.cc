#include "Geometry/HGCalMapping/interface/HGCalModuleLocator.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/EDMException.h"

//
void HGCalModuleLocator::buildLocatorFrom(std::string path,bool usefip)
{
  if(usefip) {
    edm::FileInPath fip(path);
    path=fip.fullPath();
  }
  std::ifstream file(path);

  std::string line;
  //Skip the first line of the file, which contains the column labels
  std::getline(file,line);

  while(std::getline(file, line))
  {
    std::istringstream stream(line);
    int zside;
    HGCalModuleInfo m;
    std::string itype;
    stream >> m.plane >> m.u >> m.v >> m.isSiPM >> m.isHD;
    //Some modules have an itype column, those that don't have two sequential spaces. 
    //See which we have by eating the first space, then looking for the second.
    char c;
    stream.get(c);
    std::string modtype;
    if(!isspace(stream.peek())){stream >> modtype;}
    stream >> m.econdidx >> m.captureblock >> m.slink >> m.captureblockidx >> m.fedid >> m.DAQ >> zside;
    //zside true for -1, false for +1 (matching convention from sc det ID)
    m.zside = (zside<0) ? true : false;
    mod2loc_.addParameter(m);
  }
}

HGCalModuleInfo HGCalModuleLocator::getModule(int econdidx, int captureblock, int fedid) const
{
  auto _electronicsMatch = [econdidx, captureblock, fedid](HGCalModuleInfo m){ 
    return m.econdidx == econdidx && m.captureblock == captureblock && m.fedid == fedid;
  };
  auto it = std::find_if(begin(mod2loc_.params_), end(mod2loc_.params_), _electronicsMatch);
  if(it==mod2loc_.params_.end()){
     edm::LogError ("HGCalMapping") << "econidx: " << int(econdidx) << ", captblock: " << int(captureblock) << ", fed ID: " << fedid <<std::endl;
     edm::Exception e(edm::errors::NotFound,"Failed to find module from channel info");
     throw e;
  }
  return *it;
}

HGCalModuleInfo HGCalModuleLocator::getModuleFromGeom(int plane, int u, int v, int zside, bool isSiPM) const
{
  bool binZside = (zside<0) ? true : false;
  auto _geometryMatch = [plane, u, v, binZside, isSiPM](HGCalModuleInfo m){ 
    return m.plane == plane && m.u == u && m.v == v && m.zside == binZside && m.isSiPM == isSiPM;
  };
  auto it = std::find_if(begin(mod2loc_.params_), end(mod2loc_.params_), _geometryMatch);
  if(it==mod2loc_.params_.end()){
     edm::Exception e(edm::errors::NotFound,"Failed to find module from Geom info");
     throw e;
  }

  return *it;
}

std::tuple<int,int,int> HGCalModuleLocator::getModuleLocation(int econdidx, int captureblockidx, int fedid) const
{
  auto _electronicsMatch = [econdidx, captureblockidx, fedid](HGCalModuleInfo m){ 
    return m.econdidx == econdidx && m.captureblockidx == captureblockidx && m.fedid == fedid;
  };
  auto it = std::find_if(begin(mod2loc_.params_), end(mod2loc_.params_), _electronicsMatch);
  return std::make_tuple(it->plane, it->u, it->v);
}

std::tuple<int,int,int> HGCalModuleLocator::getModuleLocation(HGCalElectronicsId& id) const
{
  uint8_t econdidx = id.econdIdx();
  uint8_t captureblockidx = id.captureBlock();
  uint16_t fedid = id.fedId();
  
  auto _electronicsMatch = [econdidx, captureblockidx, fedid](HGCalModuleInfo m){ 
    return m.econdidx == econdidx && m.captureblockidx == captureblockidx && m.fedid == fedid;
  };
  auto it = std::find_if(begin(mod2loc_.params_), end(mod2loc_.params_), _electronicsMatch);

  return std::make_tuple(it->plane, it->u, it->v);
}

int HGCalModuleLocator::getEcondIdx(int plane, int u, int v, int isSiPM) const
{
  auto _geometryMatch = [plane, u, v, isSiPM](HGCalModuleInfo m){ 
    return m.plane == plane && m.u == u && m.v == v && m.isSiPM == isSiPM;
  };
  auto it = std::find_if(begin(mod2loc_.params_), end(mod2loc_.params_), _geometryMatch);

  return it->econdidx;
}

int HGCalModuleLocator::getCaptureBlockIdx(int plane, int u, int v, int isSiPM) const
{
  auto _geometryMatch = [plane, u, v, isSiPM](HGCalModuleInfo m){ 
    return m.plane == plane && m.u == u && m.v == v && m.isSiPM == isSiPM;
  };
  auto it = std::find_if(begin(mod2loc_.params_), end(mod2loc_.params_), _geometryMatch);

  return it->captureblockidx;
}

int HGCalModuleLocator::getFedId(int plane, int u, int v, int isSiPM) const
{
  auto _geometryMatch = [plane, u, v, isSiPM](HGCalModuleInfo m){ 
    return m.plane == plane && m.u == u && m.v == v && m.isSiPM == isSiPM;
  };
  auto it = std::find_if(begin(mod2loc_.params_), end(mod2loc_.params_), _geometryMatch);

  return it->fedid;
}

std::string HGCalModuleLocator::getDAQ(int plane, int u, int v, int isSiPM) const
{
  auto _geometryMatch = [plane, u, v, isSiPM](HGCalModuleInfo m){ 
    return m.plane == plane && m.u == u && m.v == v && m.isSiPM == isSiPM;
  };
  auto it = std::find_if(begin(mod2loc_.params_), end(mod2loc_.params_), _geometryMatch);

  return it->DAQ;
}
