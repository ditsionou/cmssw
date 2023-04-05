#ifndef simcalorimetry_hgcalsimalgos_hgcalsiconfigurationbyalgo_h
#define simcalorimetry_hgcalsimalgos_hgcalsiconfigurationbyalgo_h

#include "SimCalorimetry/HGCalSimAlgos/interface/HGCalConfigurationByAlgoWrapper.h"

/**
   @short implements a si section configurator based on groups of DetIds
 */
class HGCalGeomSiSectionConfigurationByAlgo :
public HGCalConfigurationByAlgoWrapper<HGCalSiConditionsByAlgo::SiCellOpCharacteristics,HGCSiliconDetId>
{
 public:
  HGCalGeomSiSectionConfigurationByAlgo() { }
  ~HGCalGeomSiSectionConfigurationByAlgo() { }
  uint32_t toConfigurableKey(HGCSiliconDetId &d);
};

#endif
