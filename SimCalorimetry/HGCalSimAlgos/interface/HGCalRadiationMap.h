#ifndef simcalorimetry_hgcalsimalgos_hgcalradiationmap_h
#define simcalorimetry_hgcalsimalgos_hgcalradiationmap_h

#include "FWCore/Utilities/interface/Exception.h"
#include "vdt/vdtMath.h"
#include <string>
#include <map>

/**
   @class HGCalRadiationMap
   @short parses the txt file with dose/fluence parameters and implements the functions needed to 
   compute these values at a given point in space specified by the HGCal subdetector index, the layer,
   and the radius of interest
   The txt file is created outside CMSSW with the tool available in the hgcal-integration gitlab area
   https://gitlab.cern.ch/hgcal-integration/hgcal_modmap/-/blob/main/notebooks/RadiationMapFits.ipynb
 */
class HGCalRadiationMap {

public:

  //a struct to hold the parameters to be used for fluence and dose
  struct DoseParameters {
    DoseParameters()
        : a_(0.), b_(0.), c_(0.), d_(0.), e_(0.), doff_(0.), f_(0.), g_(0.), h_(0.), i_(0.), j_(0.), foff_(0.) {}
    double a_, b_, c_, d_, e_, doff_, f_, g_, h_, i_, j_, foff_;
  };

  //the key of the parameters map is a pair <sub det index, layer>
  typedef std::map<std::pair<int, int>, DoseParameters> DoseParametersMap_t;
  
  HGCalRadiationMap();
  ~HGCalRadiationMap(){};

  /**
     @shor triggers the parsing of the dose parameters from a file
   */
  void setDoseMap(const std::string &);

  /**
     @short retrieves the dose (D) value for a given subdet index, layer and radius
     if logVal enabled returns log10(D)
   */
  double getDoseValue(const int, const int, const double, bool logVal = false);


  /**
     @short retrieves the fluence (f) value for a given subdet index, layer and radius
     if logVal enabled returns log10(f)
  */
  double getFluenceValue(const int, const int, const double, bool logVal = false);

  /**
     @short retrives the full list of parameters
   */
  inline const DoseParametersMap_t &getDoseMap() { return doseMap_; }

  /**
     @short sets a common scaling factor for fluence, useful to generate scenarios with arbitrary integrated luminosity
     by default the parameters correspond to an end-of-life scenario (3/ab)
   */
  inline void setFluenceScaleFactor(double val) { fluenceSFlog10_ = log10(val); }

private:

  DoseParametersMap_t readDosePars(const std::string &);

  DoseParametersMap_t doseMap_;

  //conversion from gray to krad (1Gy=100rad=0.1krad)
  const double grayToKrad_ = 0.1;

  double fluenceSFlog10_;
};

#endif
