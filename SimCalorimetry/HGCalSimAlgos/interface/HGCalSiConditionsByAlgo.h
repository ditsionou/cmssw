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
  enum HGCalSiSensorTypes_t { HD120, LD200, LD300, HD200 };

  struct SiCellOpCharacteristicsCore {
    SiCellOpCharacteristicsCore() : cce(0.), ileak(0.) {};
    float cce, ileak;
  };
  
  struct SiCellOpCharacteristics {
    SiCellOpCharacteristics() : fluence(0.), lnfluence(0.) {}
    SiCellOpCharacteristicsCore core;
    float fluence, lnfluence, mipEqfC, capacitance;
  };

  HGCalSiConditionsByAlgo();
  ~HGCalSiConditionsByAlgo(){};

  /**
     @short set the ileak parameters to use
  */
  void setIleakParam(const std::vector<double> &pars) { ileakParam_ = pars; }

  /**
     @short set the charge collection parameters to use
  */
  void setCceParam(HGCalSiSensorTypes_t sensType,const std::vector<double> &pars) { 
    cceParam_[sensType]=pars; 
    if(sensType==HGCalSiSensorTypes_t::LD200) cceParam_[HGCalSiSensorTypes_t::HD200]=cceParam_[sensType];
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
                                              HGCalSiSensorTypes_t &sensType);
  SiCellOpCharacteristicsCore getCoreConditionsByAlgo(DetId::Detector &subdet,
                                                      int &layer,
                                                      double &radius,
                                                      HGCalSiSensorTypes_t &sensType) {
    return getConditionsByAlgo(subdet,layer,radius,sensType).core;
  }
  SiCellOpCharacteristicsCore getCoreConditionsByAlgo(unsigned int &rawId, double &radius) {

    DetId did(rawId);
    DetId::Detector det = did.det();
    int subdet = did.subdetId();

    //start by assuming this is HGCAL if not correct for HFNose (although the outcome should be the same)
    HGCSiliconDetId hgcsi(rawId);
    int layer=hgcsi.layer();
    HGCalSiSensorTypes_t sensType=(HGCalSiSensorTypes_t)hgcsi.type();
    if(subdet==ForwardSubdetector::HFNose) {
      HFNoseDetId hfnose(rawId);
      layer=hfnose.layer();
      sensType=(HGCalSiSensorTypes_t)hfnose.type();
    }

    return getCoreConditionsByAlgo(det,layer,radius,sensType);
  }
  
  /**
     @short determines the conditions using all the relevant parameterizations
   */
  SiCellOpCharacteristics getSiCellOpCharacteristics(double &cellVol,
                                                     double &capacitance,
                                                     double &mipEqfC,
                                                     std::vector<double> &cceParam,
                                                     DetId::Detector &subdet,
                                                     int &layer,
                                                     double &radius);

  std::map<HGCalSiSensorTypes_t,double> &getMipEqfC() { return mipEqfC_; }
  std::map<HGCalSiSensorTypes_t,double> &getCellCapacitance() { return cellCapacitance_; }
  std::map<HGCalSiSensorTypes_t,double> &getCellVolume() { return cellVolume_; }
  std::map<HGCalSiSensorTypes_t,std::vector<double> > &getCCEParam() { return cceParam_; }
  std::vector<double> &getIleakParam() { return ileakParam_; }

private:

  //flags used to disable specific components of the Si operation parameter
  bool ignoreFluence_, ignoreCCE_;
  
  //maps parameters for sensor types
  std::map<HGCalSiSensorTypes_t,double> mipEqfC_, cellCapacitance_, cellVolume_;
  std::map<HGCalSiSensorTypes_t,std::vector<double> > cceParam_;

  //leakage current/volume vs fluence
  std::vector<double> ileakParam_;

  //electron charge in fC
  const double qe2fc_ = 1.60217646E-4;

  //conversions
  const double unitToMicro_ = 1.e6;
  const double unitToMicroLog_ = log(unitToMicro_);
};

#endif
