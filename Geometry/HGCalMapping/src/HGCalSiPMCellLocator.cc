#include "Geometry/HGCalMapping/interface/HGCalCellLocator.h"

#include "FWCore/Utilities/interface/EDMException.h"

void HGCalCellLocator::buildLocatorFrom(std::string channelpath, std::string geometrypath)
{
  std::ifstream file;
  std::string line;
  file.open(channelpath.c_str());
  if (file.is_open())
  {
    while (std::getline(file, line))
    {
      std::istringstream stream(line);

      HGCalSiPMCellChannel c; 

      stream >> c.sipmcell >> c.plane >> c.iring >> c.iphi >> c.t >> c.trigch >> c.trigsum;
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
    edm::Exception e(edm::errors::NotFound,"Failed to match cell by channel number");
    throw e;
  }
  return std::make_tuple(it->plane, it->iring, it->iphi);
}

DetId HGCalCellLocator::getDetId(const HGCalElectronicsId& id, int seq, int z, int layer, int modiring, int modiphi) const
{
  int sipmcell = getSiPMchannel(seq, (int)id.econdeRx(), (int)id.halfrocChannel());

  auto _matchesByChannel = [sipmcell, layer, modiring](HGCalSiPMCellChannel c){
    return c.sipmcell == sipmcell && c.plane == layer && c.modiring == modiring;
  };
  auto it = std::find_if(begin(cellColl_), end(cellColl_), _matchesByChannel);
  if(it==cellColl_.end()){
    edm::Exception e(edm::errors::NotFound,"Failed to match cell by channel number");
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
