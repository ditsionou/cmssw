#include "Geometry/HGCalMapping/interface/HGCalModuleLocator.h"

#include "FWCore/Utilities/interface/Exception.h"

void HGCalModuleLocator::buildLocatorFrom(std::string path)
{
  edm::FileInPath fip(input);
  std::ifstream file(fip.fullPath());

  std::string line;
  while(std::getline(file, line))
  {
    HGCalModule m;

    int econdidx, captureblock, slink, captureblockidx;

    stream >> m.plane >> m.u >> m.v >> m.isSiPM >> econdidx >> captureblock >> slink >> captureblockidx >> m.fedid >> m.DAQ;

    m.econdidx = econdidx;
    m.captureblock = captureblock;
    m.slink = slink;
    m.captureblockidx = captureblockidx;
    ModuleToLocation_.push_back(m);
  }
}

HGCalModuleLocator::HGCalModule HGCalModuleLocator::getModule(int econdidx, int captureblockidx, int fedid) const
{
  auto _electronicsMatch = [econdidx, captureblockidx, fedid](HGCalModule m){ 
    return m.econdidx == econdidx && m.captureblockidx == captureblockidx && m.fedid == fedid;
  };
  auto it = std::find_if(begin(ModuleToLocation_), end(ModuleToLocation_), _electronicsMatch);
  return *it;
}

std::tuple<int,int,int> HGCalModuleLocator::getModuleLocation(int econdidx, int captureblockidx, int fedid) const
{
  auto _electronicsMatch = [econdidx, captureblockidx, fedid](HGCalModule m){ 
    return m.econdidx == econdidx && m.captureblockidx == captureblockidx && m.fedid == fedid;
  };
  auto it = std::find_if(begin(ModuleToLocation_), end(ModuleToLocation_), _electronicsMatch);
  return std::make_tuple(it->plane, it->u, it->v);
}

std::tuple<int,int,int> HGCalModuleLocator::getModuleLocation(const HGCalElectronicsId& id) const
{
  int econdidx = (int)id.econdIdx();
  int captureblockidx = (int)id.captureBlock();
  int fedid = (int)id.fedId();
  
  auto _electronicsMatch = [econdidx, captureblockidx, fedid](HGCalModule m){ 
    return m.econdidx == econdidx && m.captureblockidx == captureblockidx && m.fedid == fedid;
  };
  auto it = std::find_if(begin(ModuleToLocation_), end(ModuleToLocation_), _electronicsMatch);

  return std::make_tuple(it->plane, it->u, it->v);
}

int HGCalModuleLocator::getEcondIdx(int plane, int u, int v, int isSiPM) const
{
  auto _geometryMatch = [plane, u, v, isSiPM](HGCalModule m){ 
    return m.plane == plane && m.u == u && m.v == v && m.isSiPM == isSiPM;
  };
  auto it = std::find_if(begin(ModuleToLocation_), end(ModuleToLocation_), _geometryMatch);

  return it->econdidx;
}

int HGCalModuleLocator::getCaptureBlockIdx(int plane, int u, int v, int isSiPM) const
{
  auto _geometryMatch = [plane, u, v, isSiPM](HGCalModule m){ 
    return m.plane == plane && m.u == u && m.v == v && m.isSiPM == isSiPM;
  };
  auto it = std::find_if(begin(ModuleToLocation_), end(ModuleToLocation_), _geometryMatch);

  return it->captureblockidx;
}

int HGCalModuleLocator::getFedId(int plane, int u, int v, int isSiPM) const
{
  auto _geometryMatch = [plane, u, v, isSiPM](HGCalModule m){ 
    return m.plane == plane && m.u == u && m.v == v && m.isSiPM == isSiPM;
  };
  auto it = std::find_if(begin(ModuleToLocation_), end(ModuleToLocation_), _geometryMatch);

  return it->fedid;
}

std::string HGCalModuleLocator::getDAQ(int plane, int u, int v, int isSiPM) const
{
  auto _geometryMatch = [plane, u, v, isSiPM](HGCalModule m){ 
    return m.plane == plane && m.u == u && m.v == v && m.isSiPM == isSiPM;
  };
  auto it = std::find_if(begin(ModuleToLocation_), end(ModuleToLocation_), _geometryMatch);

  return it->DAQ;
}