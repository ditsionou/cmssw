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
      //in the region where measurements are available we have used a linear approximation in log (f)
      //the extrapolation to the lower fluence regions is done with an error function matched at the boundary
      if (siop.fluence > cceParam[0]) {
        siop.core.cce = cceParam[2] * siop.lnfluence + cceParam[1];
      } else {
        double at = cceParam[2] * log(cceParam[0]) + cceParam[1];
        double bt = -cceParam[0] / TMath::ErfcInverse(1. / at);
        siop.core.cce = at * TMath::Erfc((siop.fluence - cceParam[0]) / bt);
      }

      siop.core.cce = std::max((float)0., siop.core.cce);
    }
  }

  return siop;
}
