#ifndef simcalorimetry_hgcalsimalgos_hgcalconfigurationbyalgowrapper_h
#define simcalorimetry_hgcalsimalgos_hgcalconfigurationbyalgowrapper_h

#include "DataFormats/ForwardDetId/interface/HGCScintillatorDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetId.h"
#include "DataFormats/ForwardDetId/interface/HFNoseDetId.h"
#include "DataFormats/HGCalDigi/interface/HGCalElectronicsId.h"

#include "SimCalorimetry/HGCalSimAlgos/interface/HGCROCEmulator.h"
#include "SimCalorimetry/HGCalSimAlgos/interface/HGCalSiConditionsByAlgo.h"
#include "SimCalorimetry/HGCalSimAlgos/interface/HGCalSiPMonTileConditionsByAlgo.h"


/**
   @class HGCalConfigurationByAlgoWrapper
   @short this class finds the electronics configuration algorithmically given a set of DetIds
   it allows the user to retrieve the front-end configuration for a given DetId
   given the specifities of the Si and SiPM-on-tile sections and of the HFNose this class is templated
   <Conditions, DetId> or <Conditions, LogicalId>
   The baseline configuration will be encoded in a struct (HGCROCConfiguration) containing the gain, a ZS threshold, the leakage to the next bunch and the noise (total, common mode and series+parallel contributions to total)
*/
template <class C, class D>
class HGCalConfigurationByAlgoWrapper {
public:

  enum HGCalConfigurationAlgos { BYMAJORITY, BYMAXGAIN };
  
  HGCalConfigurationByAlgoWrapper() { }
  
  virtual ~HGCalConfigurationByAlgoWrapper() = default;

  /**
  @short virtual method which maps a DetId/LogicalId to a configurable key 
  and which is used to group several to sensors (e.g. to a module or a ROC)
  different classes should override this method
  */
  virtual uint32_t toConfigurableKey(D &d) = 0;

  /**
   @short adds to an internal map of separate entities which should have a common configuraiton
  */
  inline void addConfigurableToEntity(D &d, C &conf) { 
    uint32_t key = toConfigurableKey(d);
    if(configurableEntities_.count(key)==0 ) {
      std::vector<C> newentry;
      configurableEntities_[key] = newentry;
    }
    configurableEntities_[key].push_back(conf);
  }

  /**
   @short clears current maps
  */
  inline void resetConfigurableEntities() { configurableEntities_.clear(); }

  /**
     @short finds the best electronics configuration fom a map of conditions
     @param algo - the algorithm version can take the values 
                   BYMAJORITY - the configuration with largest counts is used
                   BYMAXGAIN - the channel which requires maximal gain (min dynamical range) is used
     @param resetAfter - can be used to avoid resetting the map of configurables
   */
  void findFEConfigurationByAlgo(int algo,bool resetAfter=true);

  /**
     @short returns the configuration for a given identifier
     if not found returns a default config where the gain is HGCROCDynamicRange_t::NULL
   */
  inline HGCROCConfiguration getConfigurationFor(const D &d) {
    uint32_t key = toConfigurableKey(d);
    auto it = confCache_.begin();
    if(it!=confCache_.end()) return it->second;
    return  HGCROCConfiguration();
  }

  /**
    @short fills configurations from file
  */
  void fillConfigurationsFrom(std::string);

  /**
    @short saves configurations to file
  */
  void saveConfigurationsTo(std::string);
  
private:

  //configurable entities map
  std::map<uint32_t, std::vector<C> > configurableEntities_;

  //cache of configurations associated by identifier
  std::map<uint32_t, HGCROCConfiguration> confCache_;
};

//define typedef for the different flavours of templated classes
typedef HGCalConfigurationByAlgoWrapper<HGCalSiConditionsByAlgo::SiCellOpCharacteristics,HGCSiliconDetId> HGCalSiGeoConfigurationByAlgo;
//typedef HGCalConfigurationByAlgoWrapper<HGCalSiConditionsByAlgo::SiCellOpCharacteristics,HGCalElectronicsId> HGCalSiLogConfigurationByAlgo;
//typedef HGCalConfigurationByAlgoWrapper<HGCalSiConditionsByAlgo::SiCellOpCharacteristics,HFNoseDetId> HGCalHFNoseGeoConfigurationByAlgo;
//typedef HGCalConfigurationByAlgoWrapper<HGCalSiPMonTileConditionsByAlgo::SiPMonTileCharacteristicsCore,HGCScintillatorDetId> HGCalSiPMonTileGeoConfigurationByAlgo;
//typedef HGCalConfigurationByAlgoWrapper<HGCalSiPMonTileConditionsByAlgo::SiPMonTileCharacteristicsCore,HGCalElectronicsId> HGCalSiPMonTileLogConfigurationByAlgo;

#endif
