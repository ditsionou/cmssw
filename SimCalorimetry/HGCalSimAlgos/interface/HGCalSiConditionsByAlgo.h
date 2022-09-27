#ifndef simcalorimetry_hgcalsimalgos_hgcalsiconditionsbyalgo_h
#define simcalorimetry_hgcalsimalgos_hgcalsiconditionsbyalgo_h

#include "SimCalorimetry/HGCalSimAlgos/interface/HGCalRadiationMap.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetId.h"
#include "DataFormats/ForwardDetId/interface/HFNoseDetId.h"
#include <array>
#include <vector>

/**
   @class HGCalSiConditionsByAlgo
   @short this class derives from HGCalRadiation map to parse fluence parameters and provides Si-specific calculations to be used 
   in emulating realistic conditions namely charge collection efficiency and leakage current (see DN-19-045 for details)
   By configuration the class can ignore the fluence effects on both leakage current and charge collection efficiency (HGCalSiConditions_t::FLUENCE)
   and/or the charge collection efficiency (HGCalSiConditions_t::CCE)
*/
class HGCalSiConditionsByAlgo : public HGCalRadiationMap {
public:

  enum HGCalSiConditionsAlgoBits_t { FLUENCE, CCE };

  struct SiCellOpCharacteristicsCore {
    SiCellOpCharacteristicsCore() : cce(0.), ileak(0.) {};
    float cce, ileak;
  };
  
  struct SiCellOpCharacteristics {
    SiCellOpCharacteristics() : fluence(0.), lnfluence(0.) {}
    SiCellOpCharacteristicsCore core;
    float fluence, lnfluence;
  };

  HGCalSiConditionsByAlgo();
  ~HGCalSiConditionsByAlgo(){};

  /**
     @short set the ileak parameters to use
  */
  void setIleakParam(const std::vector<double> &pars) { ileakParam_ = pars; }

  /**
     @short set the cce parameters to use
  */
  void setCceParam(const std::vector<double> &parsFine,
                   const std::vector<double> &parsThin,
                   const std::vector<double> &parsThick) {
    cceParam_.push_back(parsFine);   //120
    cceParam_.push_back(parsThin);   //200
    cceParam_.push_back(parsThick);  //300
  }

  /**
     @short overrides base class method with specifics for the configuration of the algo
  */
  void setDoseMap(const std::string &, unsigned int);

  /**
     @short returns the conditions for a pad placed a given layer,radius for a given sub-detector
  */
  SiCellOpCharacteristics getConditionsByAlgo(DetId::Detector &subdet,
                                              int &layer,
                                              double &radius,
                                              unsigned int &cellThick);
  SiCellOpCharacteristicsCore getCoreConditionsByAlgo(DetId::Detector &subdet,
                                                      int &layer,
                                                      double &radius,
                                                      unsigned int &cellThick) {
    return getConditionsByAlgo(subdet,layer,radius,cellThick).core;
  }
  SiCellOpCharacteristicsCore getCoreConditionsByAlgo(unsigned int &rawId, double &radius) {

    DetId did(rawId);
    DetId::Detector det = did.det();
    int subdet = did.subdetId();

    //start by assuming this is HGCAL if not correct for HFNose (although the outcome should be the same)
    HGCSiliconDetId hgcsi(rawId);
    int layer=hgcsi.layer();
    unsigned int cellThick=hgcsi.type();
    if(subdet==ForwardSubdetector::HFNose) {
      HFNoseDetId hfnose(rawId);
      layer=hfnose.layer();
      cellThick=hfnose.type();
    }

    return getCoreConditionsByAlgo(det,layer,radius,cellThick);
  }
  
  /**
     @short determines the conditions using all the relevant parameterizations
   */
  SiCellOpCharacteristics getSiCellOpCharacteristics(double &cellVol,
                                                     std::vector<double> &cceParam,
                                                     DetId::Detector &subdet,
                                                     int &layer,
                                                     double &radius);

  std::array<double, 3> &getMipEqfC() { return mipEqfC_; }
  std::array<double, 3> &getCellCapacitance() { return cellCapacitance_; }
  std::array<double, 3> &getCellVolume() { return cellVolume_; }
  std::vector<std::vector<double> > &getCCEParam() { return cceParam_; }
  std::vector<double> &getIleakParam() { return ileakParam_; }

private:

  //flags used to disable specific components of the Si operation parameter
  bool ignoreFluence_, ignoreCCE_;
  
  //vector of three params, per sensor type: 0:120 [mum], 1:200, 2:300
  std::array<double, 3> mipEqfC_, cellCapacitance_, cellVolume_;
  std::vector<std::vector<double> > cceParam_;

  //leakage current/volume vs fluence
  std::vector<double> ileakParam_;

  //electron charge in fC
  const double qe2fc_ = 1.60217646E-4;

  //conversions
  const double unitToMicro_ = 1.e6;
  const double unitToMicroLog_ = log(unitToMicro_);
};

#endif
