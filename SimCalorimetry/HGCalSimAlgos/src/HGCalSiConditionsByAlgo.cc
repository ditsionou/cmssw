#include "SimCalorimetry/HGCalSimAlgos/interface/HGCalSiConditionsByAlgo.h"
#include "TMath.h"
#include<iostream>

//
HGCalSiConditionsByAlgo::HGCalSiConditionsByAlgo()
  : ignoreFluence_(false),
    ignoreCCE_(false) {
  
  //fine sensors: 120 mum -  67: MPV of charge[number of e-]/mum for a mip in silicon; srouce PDG
  const double mipEqfC_120 = 120. * 67. * qe2fc_;
  mipEqfC_[0] = mipEqfC_120;
  const double cellCapacitance_120 = 50;
  cellCapacitance_[0] = cellCapacitance_120;
  const double cellVolume_120 = 0.56 * (120.e-4);
  cellVolume_[0] = cellVolume_120;

  //thin sensors: 200 mum
  const double mipEqfC_200 = 200. * 70. * qe2fc_;
  mipEqfC_[1] = mipEqfC_200;
  const double cellCapacitance_200 = 65;
  cellCapacitance_[1] = cellCapacitance_200;
  const double cellVolume_200 = 1.26 * (200.e-4);
  cellVolume_[1] = cellVolume_200;

  //thick sensors: 300 mum
  const double mipEqfC_300 = 300. * 73. * qe2fc_;
  mipEqfC_[2] = mipEqfC_300;
  const double cellCapacitance_300 = 45;
  cellCapacitance_[2] = cellCapacitance_300;
  const double cellVolume_300 = 1.26 * (300.e-4);
  cellVolume_[2] = cellVolume_300;
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
                                                                                              unsigned int &cellThick) {
  //get the appropriate parameters
  double cellVol(cellVolume_[cellThick]);
  std::vector<double> &cceParam = cceParam_[cellThick];
  
  //call baseline method and add to cache
  return getSiCellOpCharacteristics(cellVol, cceParam, subdet, layer, radius);
}

//
HGCalSiConditionsByAlgo::SiCellOpCharacteristics HGCalSiConditionsByAlgo::getSiCellOpCharacteristics(double &cellVol,
                                                                                                     std::vector<double> &cceParam,
                                                                                                     DetId::Detector &subdet,
                                                                                                     int &layer,
                                                                                                     double &radius) {
  SiCellOpCharacteristics siop;

  //leakage current and CCE [muA]
  if (ignoreFluence_) {
    siop.fluence = 0;
    siop.lnfluence = -1;
    siop.core.ileak = exp(ileakParam_[1]) * cellVol * unitToMicro_;
    siop.core.cce = 1;
  } else {
    
    if (getDoseMap().empty()) {
      throw cms::Exception("BadConfiguration")
          << " Fluence is required but no DoseMap has been passed to HGCalSiNoiseMap";
      return siop;
    }

    siop.lnfluence = getFluenceValue(subdet, layer, radius, true);
    siop.fluence = exp(siop.lnfluence);

    double conv(log(cellVol) + unitToMicroLog_);
    siop.core.ileak = exp(ileakParam_[0] * siop.lnfluence + ileakParam_[1] + conv);

    //charge collection efficiency
    if (ignoreCCE_) {
      siop.core.cce = 1.0;
    } else {
      
      float a=cceParam[0];
      float b=cceParam[1];
      siop.core.cce = a*siop.lnfluence+b;

      //regularize near the extremities (if needed)
      float x=siop.fluence;
      float xm=TMath::Exp((95-b)/a);
      float xp=TMath::Exp((5-b)/a);
      if(x<xm) {
        float m=50;
        float tm=TMath::ErfInverse(1-(a/m)*TMath::Log(xm)-b/m);
        float sigmam=-(2*m*xm)/(a*TMath::Sqrt(TMath::Pi()))*TMath::Exp(-pow(tm,2));
        float x0m=xm-sigmam*tm;
        siop.core.cce = m*(1-TMath::Erf((x-x0m)/sigmam));
      } else if(x>xp) {
        float p=50;
        float tp=TMath::ErfInverse(1-(a/p)*TMath::Log(xp)-b/p);
        float sigmap=-(2*p*xp)/(a*TMath::Sqrt(TMath::Pi()))*TMath::Exp(-pow(tp,2));
        float x0p=xp-sigmap*tp;
        siop.core.cce = p*(1-TMath::Erf((x-x0p)/sigmap));
      }
      
      siop.core.cce = siop.core.cce/100.;
    }
  }

  return siop;
}
