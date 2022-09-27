#ifndef simcalorimetry_hgcalsimalgos_hgcalconfigurationbyalgowrapper_h
#define simcalorimetry_hgcalsimalgos_hgcalconfigurationbyalgowrapper_h

#include "DataFormats/ForwardDetId/interface/HGCScintillatorDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetId.h"
#include "DataFormats/ForwardDetId/interface/HFNoseDetId.h"

#include "SimCalorimetry/HGCalSimAlgos/interface/HGCalSiConditionsByAlgo.h"
#include "SimCalorimetry/HGCalSimAlgos/interface/HGCalSiPMonTileConditionsByAlgo.h"

//this should move to electronics class (temporary for now)
namespace hgcroc {
  struct HGCROCConfiguration {
    HGCROCConfiguration() : gain(0), thrADC(0), bxp1leak(0.), enc(0.), enc_s(0.), enc_p(0.), enc_cm(0.) {}
    unsigned short gain, thrADC;
    float bxp1leak,enc,enc_s,enc_p,enc_cm;
  };
  
  enum HGCROCDynamicRange_t { q80fC, q160fC, q320fC, AUTO };
  
  typedef std::array<float, 6> HGCROCPreampPulseShape_t;
}


/**
   @class HGCalConfigurationByAlgoWrapper
   @short this class finds the electronics configuration algorithmically given a set of DetIds
   it allows the user to retrieve the front-end configuration for a given DetId
   given the specifities of the Si and SiPM-on-tile sections and of the HFNose this class is templated
   <Conditions, DetId>
   The baseline configuration is encoded in a struct (HGCROCConfiguration) containing the gain, a ZS threshold, the leakage to the next bunch and the noise (total, common mode and series+parallel contributions to total)
*/
template <class C, class D = DetId>
class HGCalConfigurationByAlgoWrapper {
public:
  
  HGCalConfigurationByAlgoWrapper();
  
  ~HGCalConfigurationByAlgoWrapper() { }

  /**
     @short finds the best electronics configuration fom a map of conditions
   */
  void findFEConfigurationByAlgo(std::map<uint32_t, C> &);

  /**
     @short returns the configuration for a given det id
   */
  hgcroc::HGCROCConfiguration getConfigurationForDetId(const D &detID);

  /**
     @short getters for private members
   */
  const hgcroc::HGCROCPreampPulseShape_t &preampPulseForGain(hgcroc::HGCROCDynamicRange_t gain) { return preampPulses_[gain];  };
  const std::vector<hgcroc::HGCROCPreampPulseShape_t> &getPreampPulses() { return preampPulses_; }
  const double &getFSCForGain(hgcroc::HGCROCDynamicRange_t gain) { return chargeAtFullScaleADCPerGain_[gain]; }
  const std::vector<double> &getFSC() { return chargeAtFullScaleADCPerGain_; }
  const double &getLSBForGain(hgcroc::HGCROCDynamicRange_t gain) { return lsbPerGain_[gain]; }
  const std::vector<double> &getLSB() { return lsbPerGain_; }
  const std::vector<double> &getENCsParamForGain(hgcroc::HGCROCDynamicRange_t gain) { return encsParam_[gain]; }
  const std::vector<std::vector<double> > &getENCsParam() { return encsParam_; }

  /**
     @short returns the TOT onset threshold (a small tolerance is subtracted)
   */
  double getTOTOnset(hgcroc::HGCROCDynamicRange_t gainIdx,double vtol=1e-6);

  /**
     @short given the gain and the cell capacitance returns an estimate of the series noise
   */
  double estimateENCs(hgcroc::HGCROCDynamicRange_t gainIdx,double cap);
  
  /**
     @short returns the parallel noise corresponding the most favorable DAQ configuraiton
   */
  double estimateOptimalENCp(double ileak);
  
private:
  
  //cache of configurations associated by identifier
  std::map<uint32_t, hgcroc::HGCROCConfiguration> confCache_;

  //electronics series noise polynomial coeffs (per gain)
  std::vector<std::vector<double> > encsParam_;

  //ADC pulses (per gain)
  std::vector<hgcroc::HGCROCPreampPulseShape_t> preampPulses_;

  //lsb and dynamic range (per gain)
  std::vector<double> lsbPerGain_, chargeAtFullScaleADCPerGain_;

};

#include "HGCalConfigurationByAlgoWrapper.icc"

template class HGCalConfigurationByAlgoWrapper<HGCalSiConditionsByAlgo::SiCellOpCharacteristics,HGCSiliconDetId>;
template class HGCalConfigurationByAlgoWrapper<HGCalSiPMonTileConditionsByAlgo::SiPMonTileCharacteristicsCore,HGCScintillatorDetId>;

#endif
