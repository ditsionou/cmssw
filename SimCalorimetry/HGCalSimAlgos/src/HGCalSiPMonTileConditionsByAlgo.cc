#include "SimCalorimetry/HGCalSimAlgos/interface/HGCalSiPMonTileConditionsByAlgo.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "CLHEP/Units/GlobalPhysicalConstants.h"

//
HGCalSiPMonTileConditionsByAlgo::HGCalSiPMonTileConditionsByAlgo()
    : refEdge_(3.),
      refDarkCurrent_(0.5),
      ignoreSiPMarea_(false),
      ignoreTileArea_(false),
      ignoreDoseScale_(false),
      ignoreFluenceScale_(false),
      ignoreNoise_(false)
{
  //number of photo electrons per MIP per scintillator type (irradiated, based on testbeam results)
  //reference is a 30*30 mm^2 tile and 2 mm^2 SiPM (with 15um pixels), at the 2 V over-voltage
  //based on https://indico.cern.ch/event/927798/contributions/3900921/attachments/2054679/3444966/2020Jun10_sn_scenes.pdf
  nPEperMIP_[HGCalTypes::TileType::TileCoarseCast] = 80.31 / 2;   //1.25deg cast
  nPEperMIP_[HGCalTypes::TileType::TileCoarseMould] = 57.35 / 2;  //1.25deg moulded
  nPEperMIP_[HGCalTypes::TileType::TileFine] = nPEperMIP_[HGCalTypes::TileType::TileCoarseCast]; //fine 1deg 
  
  maxsipmPE_[HGCalTypes::TileSiPMType::SiPMSmall] = 8888;
  maxsipmPE_[HGCalTypes::TileSiPMType::SiPMLarge] = 2*maxsipmPE_[HGCalTypes::TileSiPMType::SiPMSmall];
  maxsipmPE_[HGCalTypes::TileSiPMType::SiPMUnknown] = maxsipmPE_[HGCalTypes::TileSiPMType::SiPMSmall];
}

//
void HGCalSiPMonTileConditionsByAlgo::setDoseMap(const std::string& fullpath, const unsigned int algo) {

  //decode bits of the algo word
  ignoreSiPMarea_ = ((algo >> IGNORE_SIPMAREA) & 0x1);
  ignoreTileArea_ = ((algo >> IGNORE_TILEAREA) & 0x1);
  ignoreDoseScale_ = ((algo >> IGNORE_DOSESCALE) & 0x1);
  ignoreFluenceScale_ = ((algo >> IGNORE_FLUENCESCALE) & 0x1);
  ignoreNoise_ = ((algo >> IGNORE_NOISE) & 0x1);
  ignoreTileType_ = ((algo >> IGNORE_TILETYPE) & 0x1);

  //call base class method
  HGCalRadiationMap::setDoseMap(fullpath);
}


//
void HGCalSiPMonTileConditionsByAlgo::setNpePerMIP(HGCalTypes::TileType tile,float npe) {
  nPEperMIP_[tile] = npe;
}

//
void HGCalSiPMonTileConditionsByAlgo::setReferenceDarkCurrent(double idark) {
  refDarkCurrent_ = idark;
}

//
HGCalSiPMonTileConditionsByAlgo::SiPMonTileCharacteristics HGCalSiPMonTileConditionsByAlgo::getConditionsByAlgo(int &layer,
                                                                                                                double &radius,
                                                                                                                HGCalTypes::TileType &tileType,
                                                                                                                HGCalTypes::TileSiPMType &sipmArea) {

  bool hasDoseMap(!(getDoseMap().empty()));

  //LIGHT YIELD SCALING
  //formula is: A = A0 * exp( -D^0.65 / 199.6)
  //where A0 is the response of the undamaged detector, D is the dose
  //it will be applied if dosemap exists and this is not to be ignored
  double lyScaleFactor(1.f);
  if (!ignoreDoseScale_ && hasDoseMap) {
    double cellDose = getDoseValue(DetId::HGCalHSc, layer, radius);  //in kRad
    constexpr double expofactor = 1. / 199.6;
    const double dosespower = 0.65;
    lyScaleFactor = std::exp(-std::pow(cellDose, dosespower) * expofactor);
  }

  //SiPM NOISE
  //starting from the reference dark current, a scaling for fluence is applied
  //formula is  3.16 *  sqrt( (Idark * 1e-12) / (qe * gain) * (F / F0) )
  //where F is the fluence (neq/cm2), gain is the SiPM gain, qe is the electron charge (C), Idark is dark current (mA)
  double noise(0.f),Rdark(0.f);
  if (!ignoreNoise_ && hasDoseMap) {
 
    constexpr double refFluence(2.0E+13);
    constexpr double refGain(235000.);
    Rdark = (refDarkCurrent_ * 1E-12) / (CLHEP::e_SI * refGain);
    if (!ignoreFluenceScale_ && hasDoseMap) {
      double cellFluence = getFluenceValue(DetId::HGCalHSc, layer, radius);  //in 1-Mev-equivalent neutrons per cm2
      Rdark *= (cellFluence / refFluence);
    }
    noise = 3.16 * sqrt(Rdark);
  }
  
  //TILE and SIPM AREA scalings
  double tileAreaSF = scaleByTileArea(tileType, radius);
  double sipmAreaSF = scaleBySipmArea(sipmArea);

  lyScaleFactor *= tileAreaSF * sipmAreaSF;
  noise *= sqrt(sipmAreaSF);

  //final signal depending on scintillator type
  double S(nPEperMIP_[HGCalTypes::TileType::TileCoarseCast]);
  if (!ignoreTileType_)
    S = nPEperMIP_[tileType];
  
  HGCalSiPMonTileConditionsByAlgo::SiPMonTileCharacteristics sipmChar;
  sipmChar.core.npe = S;
  sipmChar.core.darkPx = noise;
  sipmChar.lysf = lyScaleFactor;
  sipmChar.Rdark = Rdark;
  return sipmChar;
}

//
double HGCalSiPMonTileConditionsByAlgo::scaleByTileArea(const HGCalTypes::TileType &tileType, const double radius) {

  double scaleFactor(1.f);
  
  if (ignoreTileArea_)
    return scaleFactor;

  double edge(refEdge_);  //start with reference 3cm of edge
  if (tileType==HGCalTypes::TileType::TileFine) {
    constexpr double factor = 2 * M_PI * 1. / 360.;
    edge = radius * factor;  //1 degree
  } else {
    constexpr double factor = 2 * M_PI * 1. / 288.;
    edge = radius * factor;  //1.25 degrees
  }
  scaleFactor = refEdge_ / edge;
  return scaleFactor;
}

//
double HGCalSiPMonTileConditionsByAlgo::scaleBySipmArea(const HGCalTypes::TileSiPMType &sipmArea) {
  double scaleFactor(1.f);

  if (ignoreSiPMarea_)
    return scaleFactor;

  if(sipmArea==HGCalTypes::TileSiPMType::SiPMLarge)
    scaleFactor = 2.f;
  
  return scaleFactor;
}
