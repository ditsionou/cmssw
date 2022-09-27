#ifndef simcalorimetry_hgcalsimalgos_hgcalconditionsbyalgowrapper_h
#define simcalorimetry_hgcalsimalgos_hgcalconditionsbyalgowrapper_h

#include "SimCalorimetry/HGCalSimAlgos/interface/HGCalSiPMonTileConditionsByAlgo.h"
#include "SimCalorimetry/HGCalSimAlgos/interface/HGCalSiConditionsByAlgo.h"
#include "DataFormats/ForwardDetId/interface/HGCScintillatorDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetId.h"
#include "DataFormats/ForwardDetId/interface/HFNoseDetId.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include <string>

/**
   @class HGCalConditionsByAlgoWrapper
   @short this class interfaces the geometry with the base HGCal*ConditionsByAlgo classes
   it contains the geometry/topology pointers as well as the pointer to the corresponding HGCal*Conditions algo class
   it allows the user to retrieve the conditions for a given DetId
   given the specifities of the Si and SiPM-on-tile sections and of the HFNose this class is templated
   <Algo, Conditions, DetId>
*/
template <class A, class C, class D = DetId>
class HGCalConditionsByAlgoWrapper {
public:

  HGCalConditionsByAlgoWrapper() { }
  
  ~HGCalConditionsByAlgoWrapper() { }

  /**
     @short sets the geometry 
     a cache mapping the DetIds to conditions is created the first time this method is called
     the filling of the cache is ignored if it has already been filled
   */
  void setGeometry(const CaloSubdetectorGeometry *);

  /**
     @short fills the conditions cache
     by default it will skip if there are already entries in the cache
   */
  void fillCache(bool forceFill=false);
  
  /**
     @short returns the conditions for a given det id
   */
  C getConditionsForDetId(const D &detID);

  /**
     @short wraps the procedure of getting the radius for a given DetId
   */
  double computeRadius(const D &detId);
  
  /**
     @short retrieve conditions algorithm being used
   */
  A &getConditionsAlgo();

  /**
     @short retrieve the cache (DetId to conditions struct)
   */
  std::map<uint32_t, C> &getConditionsCache();

  /**
     @short gets the geometry being used by the wrapper
  */
  const HGCalGeometry *geom();

  /**
     @short gets the topology associated to the geometry being used by the wrapper
   */
  const HGCalTopology *topo();

  /**
     @short returns the DDD constants associated to the geometry being used by the wrapper
   */
  const HGCalDDDConstants *ddd();
  
private:

  /**
     @short substitutes the original detId word by another equivalent
     (profits from the detector geometry and it is used to reduce the number of keys stored in the cache)
   */
  unsigned int getEquivalentDetId(DetId rawId);
  
  //geometry pointers
  const HGCalGeometry *hgcalGeom_;
  const HGCalTopology *hgcalTopology_;
  const HGCalDDDConstants *hgcalDDD_;

  //conditions algorithm class
  A condAlgo_;
  
  //cache of conditions associated to a detID
  std::map<uint32_t, C> condCache_;
};


#include "HGCalConditionsByAlgoWrapper.icc"

template class HGCalConditionsByAlgoWrapper<HGCalSiConditionsByAlgo,
                                            HGCalSiConditionsByAlgo::SiCellOpCharacteristicsCore,
                                            HGCSiliconDetId>;
template class HGCalConditionsByAlgoWrapper<HGCalSiConditionsByAlgo,
                                            HGCalSiConditionsByAlgo::SiCellOpCharacteristicsCore,
                                            HFNoseDetId>;
template class HGCalConditionsByAlgoWrapper<HGCalSiPMonTileConditionsByAlgo,
                                            HGCalSiPMonTileConditionsByAlgo::SiPMonTileCharacteristicsCore,
                                            HGCScintillatorDetId>;

typedef HGCalConditionsByAlgoWrapper<HGCalSiConditionsByAlgo,
                                     HGCalSiConditionsByAlgo::SiCellOpCharacteristicsCore,
                                     HGCSiliconDetId> HGCSiConditionsByAlgoWrapper;
typedef HGCalConditionsByAlgoWrapper<HGCalSiConditionsByAlgo,
                                     HGCalSiConditionsByAlgo::SiCellOpCharacteristicsCore,
                                     HFNoseDetId> HFConditionsByAlgoWrapper;
typedef HGCalConditionsByAlgoWrapper<HGCalSiPMonTileConditionsByAlgo,
                                     HGCalSiPMonTileConditionsByAlgo::SiPMonTileCharacteristicsCore,
                                     HGCScintillatorDetId> HGCSiPMonTileConditionsByAlgoWrapper;

#endif
