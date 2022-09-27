#include "SimCalorimetry/HGCalSimAlgos/interface/HGCalRadiationMap.h"
#include "DataFormats/DetId/interface/DetId.h"
#include <iostream>
#include <cassert>
#include <string>
#include <cmath>
#include <chrono>
#include <random>

// run for instance with:
//
//                         time HGCalRadiationMapTest n_trials path_to_radmap
//                         time HGCalRadiationMapTest 1000 SimCalorimetry/HGCalSimProducers/data/doseParams_3000fb_fluka-6.2.0.1.txt
// for a measureble amount of time taken

int main(int argc, char** argv) {
  
  std::cout << "Basic test of radiation map parser " << std::endl;

  if (argc<2) {
    std::cout << "Usage: HGCalRadiationMapTest n_trials path_to_params" << std::endl;
    return -1;
  }
  
  //parse the arguments from command line
  unsigned long int ntries = std::stoul(argv[1], nullptr, 0);
  std::cout << "\t will test " << ntries << " random tries" << std::endl;
  
  std::string filepath(argv[2]);
  std::cout << "\t + input file: " << filepath << std::endl;

  //start a linear congruential engine
  unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
  std::minstd_rand0 myrand(seed1);
  
  //start the radiation map object
  HGCalRadiationMap mapper;
  mapper.setDoseMap(filepath);
  const HGCalRadiationMap::DoseParametersMap_t &paramMap = mapper.getDoseMap();

  //loop over available sub-detectors and layers
  for(auto entry : paramMap) {
    
    int subdet=entry.first.first;
    int layer=entry.first.second;
    HGCalRadiationMap::DoseParameters params=entry.second;
    float min_radius=std::max(params.doff_,params.foff_);
    std::cout << "Testing subdet=" << subdet << " @ layer=" << layer;
    
    //iterate over ntries
    for(unsigned long i=0; i<ntries; ++i) {

      //test radius in the [min_radius-10,300+min_radius] cm range
      double radius = myrand() % 310 + (min_radius-10); 
    
      //test nominal predictions
      mapper.setFluenceScaleFactor(1.0);
      float dose_eol=mapper.getDoseValue(subdet,layer,radius);
      float fluence_eol=mapper.getFluenceValue(subdet,layer,radius);
      if(subdet==DetId::Detector::HGCalHSc) {
        assert(!std::isnan(dose_eol));
        assert(dose_eol>=0);
      }
      assert(!std::isnan(fluence_eol));
      assert(fluence_eol>=0);
         
      //for meaning values (we expect fluences >1e14 neq/cm2)
      //assert fluence scaling for a random scenario in the [0.5,2.0] range
      if(fluence_eol>1e12) {
        float sf = 1.5*(myrand() % 100)/100.+0.5;    
        mapper.setFluenceScaleFactor(sf);
        float fluence=mapper.getFluenceValue(subdet,layer,radius);
        assert(!std::isnan(fluence));
        assert(fluence>=0);
        assert(fabs(fluence-sf*fluence_eol)<=(1e-3)*fluence_eol);
      }
    }
    
    std::cout << "...OK" << std::endl;
  }
  
  return 0;
}
