#include "SimCalorimetry/HGCalSimAlgos/interface/HGCalSiConditionsByAlgo.h"
#include "TMath.h"
#include<iostream>

//
HGCalSiConditionsByAlgo::HGCalSiConditionsByAlgo()
  : ignoreFluence_(false),
    ignoreCCE_(false) {
  
  //120 mum -  67: MPV of charge[number of e-]/mum for a mip in silicon; source PDG
  const double mipEqfC_120 = 120. * 67. * qe2fc_;
  mipEqfC_[HGCalSiSensorTypes_t::HD120] = mipEqfC_120;
  const double cellCapacitance_120 = 50;
  cellCapacitance_[HGCalSiSensorTypes_t::HD120] = cellCapacitance_120;
  const double cellVolume_120 = 0.56 * (120.e-4);
  cellVolume_[HGCalSiSensorTypes_t::HD120] = cellVolume_120;

  //200 mum
  const double mipEqfC_200 = 200. * 70. * qe2fc_;
  mipEqfC_[HGCalSiSensorTypes_t::LD200] = mipEqfC_200;
  mipEqfC_[HGCalSiSensorTypes_t::HD200] = mipEqfC_200;
  const double cellCapacitance_LD200 = 65;
  const double kfactLD2HD = 0.56/1.26;
  cellCapacitance_[HGCalSiSensorTypes_t::LD200] = cellCapacitance_LD200;
  cellCapacitance_[HGCalSiSensorTypes_t::HD200] = cellCapacitance_LD200 * kfactLD2HD;
  const double cellVolume_200 = 1.26 * (200.e-4);
  cellVolume_[HGCalSiSensorTypes_t::LD200] = cellVolume_200;
  cellVolume_[HGCalSiSensorTypes_t::HD200] = cellVolume_200 * kfactLD2HD;

  //300 mum
  const double mipEqfC_300 = 300. * 73. * qe2fc_;
  mipEqfC_[HGCalSiSensorTypes_t::LD300] = mipEqfC_300;
  const double cellCapacitance_300 = 45;
  cellCapacitance_[HGCalSiSensorTypes_t::LD300] = cellCapacitance_300;
  const double cellVolume_300 = 1.26 * (300.e-4);
  cellVolume_[HGCalSiSensorTypes_t::LD300] = cellVolume_300;
}

//
void HGCalSiConditionsByAlgo::setDoseMap(const std::string &fullpath, unsigned int algo) {

  //decode bits in the algo word
  ignoreFluence_ = ((algo >> FLUENCE) & 0x1);
  ignoreCCE_ = ((algo >> CCE) & 0x1);

  //call base class method
  HGCalRadiationMap::setDoseMap(fullpath);
}

//
HGCalSiConditionsByAlgo::SiCellOpCharacteristics HGCalSiConditionsByAlgo::getConditionsByAlgo(DetId::Detector &subdet,
                                                                                              int &layer,
                                                                                              double &radius,
                                                                                              HGCalSiSensorTypes_t &sensType) {
  //get the appropriate parameters
  double cellVol(cellVolume_[sensType]);
  double capacitance(cellCapacitance_[sensType]);
  double mipEqfC(mipEqfC_[sensType]);
  std::vector<double> &cceParam = cceParam_[sensType];
  
  //call baseline method and add to cache
  return getSiCellOpCharacteristics(cellVol, capacitance, mipEqfC, cceParam, subdet, layer, radius);
}

//
HGCalSiConditionsByAlgo::SiCellOpCharacteristics HGCalSiConditionsByAlgo::getSiCellOpCharacteristics(double &cellVol,
                                                                                                     double &capacitance,
                                                                                                     double &mipEqfC,
                                                                                                     std::vector<double> &cceParam,
                                                                                                     DetId::Detector &subdet,
                                                                                                     int &layer,
                                                                                                     double &radius) {
  SiCellOpCharacteristics siop;  
  siop.capacitance = capacitance;

  //leakage current and CCE [muA]
  if (ignoreFluence_) {
    siop.fluence = 0;
    siop.lnfluence = -1;
    siop.core.ileak = 0;
    siop.core.cce = 1;
  } else {
    
    if (getDoseMap().empty()) {
      throw cms::Exception("BadConfiguration")
          << " Fluence is required but no DoseMap has been passed to HGCalSiNoiseMap";
      return siop;
    }

    siop.lnfluence = getFluenceValue(subdet, layer, radius, true);
    siop.fluence = exp(siop.lnfluence);
    siop.core.ileak = evaluateLeakageCurrent(cellVol,siop.lnfluence);

    //charge collection efficiency
    siop.core.cce = ignoreCCE_ ? 1.0 : evaluateCCE(siop.fluence,siop.lnfluence,cceParam);

    //apply CCE to the mip equivalent
    siop.mipEqfC = mipEqfC * siop.core.cce;
  }

  return siop;
}

//
HGCalSiConditionsByAlgo::SiCellOpCharacteristics HGCalSiConditionsByAlgo::getSiCellOpCharacteristics(double &fluence,
                                                                                                     double &cellVol,
                                                                                                     double &capacitance,
                                                                                                     double &mipEqfC,
                                                                                                     std::vector<double> &cceParam) {
  SiCellOpCharacteristics siop;  
  siop.capacitance = capacitance;
  siop.lnfluence = log(fluence);
  siop.fluence = fluence;
  siop.core.ileak = evaluateLeakageCurrent(cellVol,siop.lnfluence);
  siop.core.cce = evaluateCCE(siop.fluence,siop.lnfluence,cceParam);
  siop.mipEqfC = mipEqfC * siop.core.cce;
  return siop;
}

//
HGCalSiConditionsByAlgo::SiCellOpCharacteristics HGCalSiConditionsByAlgo::getSiCellOpCharacteristics(double &fluence,HGCalSiSensorTypes_t &sens) {
  return getSiCellOpCharacteristics(fluence, cellVolume_[sens], cellCapacitance_[sens], mipEqfC_[sens], cceParam_[sens]);
}

//
float HGCalSiConditionsByAlgo::evaluateLeakageCurrent(float cellvol, float lnfluence) {  
  return cellvol*vdt::fast_exp(lnAlpha_ + lnfluence + lnUnitToMicro_);
}

//
float HGCalSiConditionsByAlgo::evaluateCCE(float fluence,float lnfluence,std::vector<double> &cceParam) {
  
  double a=cceParam[0];
  double b=cceParam[1];
  double cce = a*lnfluence+b;

  //regularize near the extremities (if needed)
  double x=fluence;
  double xm=TMath::Exp((95-b)/a);
  double xp=TMath::Exp((5-b)/a);
  if(x<xm) {
    double m=50;
    double tm=TMath::ErfInverse(1-(a/m)*TMath::Log(xm)-b/m);
    double sigmam=-(2*m*xm)/(a*TMath::Sqrt(TMath::Pi()))*TMath::Exp(-pow(tm,2));
    double x0m=xm-sigmam*tm;
    cce = m*(1-TMath::Erf((x-x0m)/sigmam));
  } else if(x>xp) {
    double p=50;
    double tp=TMath::ErfInverse(1-(a/p)*TMath::Log(xp)-b/p);
    double sigmap=-(2*p*xp)/(a*TMath::Sqrt(TMath::Pi()))*TMath::Exp(-pow(tp,2));
    double x0p=xp-sigmap*tp;
    cce = p*(1-TMath::Erf((x-x0p)/sigmap));
  }
  
  return cce/100.;
}
