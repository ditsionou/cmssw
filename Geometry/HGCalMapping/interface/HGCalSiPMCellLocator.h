#ifndef Geometry_HGCalMapping_HGCalCellLocator_H
#define Geometry_HGCalMapping_HGCalCellLocator_H 1

#include "DataFormats/HGCalDigi/interface/HGCalElectronicsId.h"
#include "DataFormats/ForwardDetId/interface/HGCScintillatorDetId.h"

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <tuple>
#include <iterator>
#include <algorithm>

class HGCalCellLocator {

    struct HGCalSiPMCellChannel{
        int sipmcell,plane,iring,iphi,trigch,trigsum;
        std::string t;
    };

    struct HGCalSiPMCellType{
        int modplane,modiring,modiphi;
        std::string t;
    };

    public:

        HGCalCellLocator(){};

        void buildLocatorFrom(std::string channelpath, std::string geometrypath);

        // // Cell location from ROC fields and Module location
        std::tuple<int,int,int> getCellLocation(int seq, int econderx, int halfrocch, int layer, int modiring, int modiphi) const;

        // // DetId from ElectronicsId and Module location, including z-side
        DetId getDetId(const HGCalElectronicsId& id, int seq, int z, int layer, int modiring, int modiphi) const;

    private:
        std::vector<HGCalSiPMCellChannel> cellColl_;
        std::vector<HGCalSiPMCellType> cellType_;

        // Module type from module location
        std::string getModuleType(int layer, int modiring, int modiphi) const;

        // SiPM channel number from sequence number and ROC fields
        int getSiPMchannel(int seq, int econderx, int halfrocch) const;
};

#endif