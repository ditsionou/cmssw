#include "CondFormats/HGCalObjects/src/headers.h"

namespace CondFormats_HGCalObjects {

  std::vector<int> v_i;
  std::pair<std::string, std::vector<int> > p_s_v_i;
  std::map<std::string, std::vector<int> > m_s_v_i;
  HGCalCondSerializableGenericConfig h_csgc();
  
  HGCalCondSerializableSiCellChannelInfo::HGCalSiCellChannelInfo hscci;
  std::vector<HGCalCondSerializableSiCellChannelInfo::HGCalSiCellChannelInfo> v_hscci;
  HGCalCondSerializableSiCellChannelInfo h_cscci();

  HGCalCondSerializableSiPMTileInfo::HGCalSiPMTileInfo hsti;
  std::vector<HGCalCondSerializableSiPMTileInfo::HGCalSiPMTileInfo> v_hsti;
  HGCalCondSerializableSiPMTileInfo h_csti();
  
}  // namespace CondFormats_HGCalObjects
