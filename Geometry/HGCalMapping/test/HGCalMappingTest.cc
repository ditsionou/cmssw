#include "Geometry/HGCalMapping/interface/HGCalCellLocator.h"

#include "DataFormats/ForwardDetId/interface/HGCScintillatorDetId.h"

#include <iostream>
#include <cassert>
#include <string>
#include <vector>
#include <cmath>

#include "FWCore/Utilities/interface/EDMException.h"

void testSiPMCellLocator(std::string channelmap, std::string geometrymap, std::string modulemap)
{
  std::cout << "Testing of HGCalSiPMCellLocator class" << std::endl;

  HGCalCellLocator celllocator;
  celllocator.buildLocatorFrom(channelmap, geometrymap);

  int plane,modiu,modiv,isSiPM,z(0);
  int econdidx,captureblock,slink,captureblockidx,fedid,econderx(0),halfrocch(0),seq;
  std::string DAQ;

  std::ifstream file;
  std::string line;
  file.open(modulemap.c_str());
  if (file.is_open())
  {
    std::getline(file, line);
    while(std::getline(file, line))
    {
      std::istringstream stream(line);
      stream >> plane >> modiu >> modiv >> isSiPM >> econdidx >> captureblock >> slink >> captureblockidx >> fedid >> DAQ;

      if(isSiPM)
      {
        for(seq = 0; seq < 20; seq++) {
          // Calibration and common mode channels
          if (seq == 8 || seq == 17 || seq == 18) continue;

          HGCalElectronicsId eid(fedid, captureblockidx, econdidx, econderx, halfrocch);                 
          
          HGCScintillatorDetId detid = celllocator.getDetId(eid, seq, z, plane, modiu, modiv);
          assert(detid.sipm());
          assert(!std::isnan(detid.iphi()));
          assert(!std::isnan(detid.ring()));
          assert(detid.ring()>0 && detid.layer()>0);
          assert(detid.iphi()<=360);
        }
      }
    }
  }
  std::cout << "HGCalSiPMCellLocator done" << std::endl;
}

int main(int argc, char** argv) {

  if (argc<3) {
      std::cout << "Usage: HGCalMappingTest path_to_channels_map path_to_goemetry_map path_to_module_map" << std::endl;
      return -1;
  }
  std::string channelsmap(argv[1]);
  std::string geometrymap(argv[2]);
  std::string modulemap(argv[3]);
  testSiPMCellLocator(channlesmap, geometrymap, modulemap);
  return 0;
}