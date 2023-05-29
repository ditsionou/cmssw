#include "CondFormats/Serialization/interface/Test.h"
#include "CondFormats/HGCalObjects/src/headers.h"

int main()
{
  //generic configurables
  testSerialization<HGCalCondSerializableGenericConfig>();

  //si cell
  testSerialization<HGCalCondSerializableSiCellChannelInfo::HGCalSiCellChannelInfo>();
  testSerialization<std::vector<HGCalCondSerializableSiCellChannelInfo::HGCalSiCellChannelInfo>>();
  testSerialization<HGCalCondSerializableSiCellChannelInfo>();
  testSerialization<std::vector<HGCalCondSerializableSiCellChannelInfo>>();

  //sipm-on-tile cell
  testSerialization<HGCalCondSerializableSiPMTileInfo::HGCalSiPMTileInfo>();
  testSerialization<std::vector<HGCalCondSerializableSiPMTileInfo::HGCalSiPMTileInfo>>();
  testSerialization<HGCalCondSerializableSiPMTileInfo>();
  testSerialization<std::vector<HGCalCondSerializableSiPMTileInfo>>();

  return 0;
}
