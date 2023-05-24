#include "Geometry/HGCalMapping/interface/HGCalCellLocator.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/Utilities/interface/EDMException.h"

void HGCalCellLocator::buildLocatorFrom(std::string channelpath)
{
  edm::FileInPath fip(channelpath);
  std::ifstream file(fip.fullPath());
  std::string line;
  if (file.is_open())
  {
    while (std::getline(file, line))
    {
      std::istringstream stream(line);

      HGCalSiPMCellChannel c; 

      stream >> c.sipmcell >> c.plane >> c.iring >> c.iphi >> c.t >> c.trigch >> c.trigsum >> c.modiring;
      cellColl_.push_back(c);
    }
  }
  else
  {
    edm::Exception e(edm::errors::FileOpenError, "HGCalSiPMCellLocator::buildLocatorFrom : SiPM channel mapping file can not be found.");
    throw e;
  }
  file.close();

}

int HGCalCellLocator::getSiPMchannel(int seq, int econderx, int halfrocch) const
{
  return seq + 72*(int)(econderx/2) + 36*halfrocch;
}

std::tuple<int,int,int> HGCalCellLocator::getCellLocation(int seq, int econderx, int halfrocch, int layer, int modiring, int modiphi) const
{
  int sipmcell = getSiPMchannel(seq, econderx, halfrocch);

  auto _matchesByChannel = [sipmcell, layer, modiring](HGCalSiPMCellChannel c){
    return c.sipmcell == sipmcell && c.plane == layer && c.modiring == modiring;
  };
  auto it = std::find_if(begin(cellColl_), end(cellColl_), _matchesByChannel);
  if(it==cellColl_.end()){
    edm::Exception e(edm::errors::NotFound,"HGCalSiPMCellLocator::getCellLocation: Failed to match cell by channel number");
    throw e;
  }
  return std::make_tuple(it->plane, it->iring, it->iphi);
}


std::tuple<int,int> HGCalCellLocator::getCellLocation(const HGCalElectronicsId& id, int seq, int layer, int modiring, int modiphi) const
{
  int sipmcell = getSiPMchannel(seq, (int)id.econdeRx(), (int)id.halfrocChannel());

  auto _matchesByChannel = [sipmcell, layer, modiring](HGCalSiPMCellChannel c){
    return c.sipmcell == sipmcell && c.plane == layer && c.modiring == modiring;
  };
  auto it = std::find_if(begin(cellColl_), end(cellColl_), _matchesByChannel);
  if(it==cellColl_.end()){
    edm::Exception e(edm::errors::NotFound,"HGCalSiPMCellLocator::getCellLocation: Failed to match cell by channel number");
    throw e;
  }
  int celliring = it->iring;
  int celliphi = it->iphi;

  return std::make_tuple(celliring, celliphi);
}

DetId HGCalCellLocator::getDetId(const HGCalElectronicsId& id, int seq, int z, int layer, int modiring, int modiphi) const
{
  int sipmcell = getSiPMchannel(seq, (int)id.econdeRx(), (int)id.halfrocChannel());

  auto _matchesByChannel = [sipmcell, layer, modiring](HGCalSiPMCellChannel c){
    return c.sipmcell == sipmcell && c.plane == layer && c.modiring == modiring;
  };
  auto it = std::find_if(begin(cellColl_), end(cellColl_), _matchesByChannel);
  if(it==cellColl_.end()){
    edm::Exception e(edm::errors::NotFound,"HGCalSiPMCellLocator::getDetId: Failed to match cell by channel number");
    throw e;
  }

  int celliring = it->iring;
  int celliphi = it->iphi;

  // Layer desription has an offset of 25
  int idlayer = layer - 25;
  int idtype = ((idlayer <= 8) ? 0 : ((idlayer <= 17) ? 1 : 2));
  int ring = ((z == 0) ? celliring : (-1)*celliring);
  int iphi = modiphi*10 + celliphi;

  HGCScintillatorDetId detid(idtype, idlayer, ring, iphi, false, true);
  return detid;
}

std::tuple<int,int,int> HGCalCellLocator::getModuleLocation(DetId& id) const
{
  DetId::Detector subdet = id.det();

  if (subdet == DetId::Detector::HGCalHSc)
  {
    HGCScintillatorDetId detid(id); 
    
    // Layer desription has an offset of 25
    int layer = detid.layer() + 25;
    int cellring = detid.ring();
    int modiring;
    if (layer <= 37)
    {
      modiring = ((cellring <= 25) ? 0 : 1);
    }
    else if (layer <= 43)
    {
      modiring = ((cellring <= 17) ? 0 : ((cellring <= 25 ) ? 1 : ((cellring <= 33) ? 2 : 3)));
    }
    else
    {
      modiring = ((cellring <= 5) ? 0 : ((cellring <= 17 ) ? 1 : ((cellring <= 25) ? 2 : ((cellring <= 33) ? 3 : 4))));
    }
    int modiphi = detid.iphi()/10;

    return std::make_tuple(layer,modiring,modiphi);
  }
  else
  {
    throw cms::Exception("InvalidDetId") << "Wrong HGCal DetId::Detector in HGCalCellLocator::getModuleLocation.";
  
  }
}