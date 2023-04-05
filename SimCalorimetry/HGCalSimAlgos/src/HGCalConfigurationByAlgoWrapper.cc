#include "SimCalorimetry/HGCalSimAlgos/interface/HGCalConfigurationByAlgoWrapper.h"
#include "SimCalorimetry/HGCalSimAlgos/interface/HGCROCEmulator.h"
#include <algorithm>

//
template <class C, class D>
void HGCalConfigurationByAlgoWrapper<C,D>::findFEConfigurationByAlgo(int algo,bool resetAfter) {

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
}

//
template <class C, class D>
void HGCalConfigurationByAlgoWrapper<C,D>::fillConfigurationsFrom(std::string url){
  //FIXME
}


//
template <class C, class D>
void HGCalConfigurationByAlgoWrapper<C,D>::saveConfigurationsTo(std::string url){
  //FIXME
}
