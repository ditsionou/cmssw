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
enum HGCalConfigurationAlgos { BYMAJORITY, BYMAXGAIN };

template <class C, class D>
class HGCalConfigurationByAlgoWrapper {
public:
  
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
  inline void findFEConfigurationByAlgo(int algo=HGCalConfigurationAlgos::BYMAXGAIN,
                                        bool resetAfter=false) {
    confCache_.clear();
    HGCROCEmulator<HGCROCChannelDataFrameSpec> roc;
    
    //loop over the current map of configurables
    for(auto it : configurableEntities_) {
      
      auto key = it.first;
      auto condColl = it.second;
      
      //for each sensor get the best configuration from the ROC based on the expected mip equivalent position
      std::map<HGCROCDynamicRange,int> counts;
      for(auto cond : condColl) {
        float S(cond.mipEqfC);
        HGCROCDynamicRange gain = roc.proposeConfig(S).gain;
        if(counts.count(gain)==0) counts[gain]=0;
        counts[gain]+=1;
      }
      
      HGCROCConfiguration cfg;
      cfg.opMode=HGCROCOperationMode::DEFAULT;
      cfg.gain=HGCROCDynamicRange::NULLGAIN;
      
      //the first entry corresponds to the max. gain possible
      if(algo==HGCalConfigurationAlgos::BYMAXGAIN) {
        cfg.gain=counts.begin()->first;
      }
      
      //use the configuration preferred by the majority instead
      else if(algo==HGCalConfigurationAlgos::BYMAJORITY){
        using pair_type = decltype(counts)::value_type;
        auto itr = std::max_element(counts.begin(),
                                    counts.end(),
                                    [](const pair_type a, const pair_type b) { return a.second < b.second; }
                                    );      
        cfg.gain=itr->first;      
      }
      
      confCache_[key]=cfg;
    }

    if(resetAfter) configurableEntities_.clear();
  }

  /**
     @short returns the configuration for a given identifier
     if not found returns a default config where the gain is HGCROCDynamicRange_t::NULL
   */
  inline HGCROCConfiguration getConfigurationFor(D &d) {
    uint32_t key = toConfigurableKey(d);
    auto it = confCache_.find(key);
    if(it!=confCache_.end()) return it->second;
    std::cout << "Can't find config for " << d << std::endl;
    return  HGCROCConfiguration();
  }

  /**
    @short fills configurations from file
  */
  void fillConfigurationsFrom(std::string) {
    //FIXME
  }

  /**
    @short saves configurations to file
  */
  void saveConfigurationsTo(std::string) {
    //FIXME
  }

  
private:

  //configurable entities map
  std::map<uint32_t, std::vector<C> > configurableEntities_;

  //cache of configurations associated by identifier
  std::map<uint32_t, HGCROCConfiguration> confCache_;
};

#endif
