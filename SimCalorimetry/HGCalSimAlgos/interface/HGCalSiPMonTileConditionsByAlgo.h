#ifndef simcalorimetry_hgcalsimalgos_hgcalscinoisemap
#define simcalorimetry_hgcalsimalgos_hgcalscinoisemap

#include "SimCalorimetry/HGCalSimAlgos/interface/HGCalRadiationMap.h"
#include "Geometry/HGCalCommonData/interface/HGCalTypes.h"
#include "DataFormats/ForwardDetId/interface/HGCScintillatorDetId.h"
#include <string>
#include <array>

/**
   @class HGCalSiPMonTileConditionsByAlgo
   @short derives from HGCalRadiation map to parse fluence/dose parameters and provides Sci-specific functions
   to compute the conditions for each SiPM-on-tile
   the algo word (set at configuration level) is used to control which aspects are simulated
   bit 1 - ignores the scaling of signal and noise with SIPMAREA
   bit 2 - instead of the geometry-based SiPM area (from detId, if available) use the boundaries read from a txt file
   bit 3 - ignores the scaling of the signal light yield with the tile area
   bit 4 - ignores the scaling of the light yield with the dose
   bit 5 - ignores the scaling of the noise with the fluence (=constant noise scenario)
   bit 6 - ignores noise
   bit 7 - ignore tile type (fallback on CAST)
   bit 8 - ignore pedestal subtraction
*/

class HGCalSiPMonTileConditionsByAlgo : public HGCalRadiationMap {
public:
  
  enum HGCalSiPMonTileConditionsAlgoBits_t {
    IGNORE_SIPMAREA,
    IGNORE_TILEAREA,
    IGNORE_DOSESCALE,
    IGNORE_FLUENCESCALE,
    IGNORE_NOISE,
    IGNORE_TILETYPE,
    IGNORE_AUTO_PEDESTALSUB
  };

  struct SiPMonTileCharacteristicsCore {
    SiPMonTileCharacteristicsCore() : npe(0.), darkPx(0.) {}
    float npe,darkPx;
  };
  
  struct SiPMonTileCharacteristics {
    SiPMonTileCharacteristics() : lysf(0.), Rdark(0.) {}
    SiPMonTileCharacteristicsCore core;
    float lysf, Rdark;
  };

  HGCalSiPMonTileConditionsByAlgo();
  ~HGCalSiPMonTileConditionsByAlgo(){};

  /**
     @short overrides base class method with specifics for the configuration of the algo 
  */
  void setDoseMap(const std::string &, const unsigned int);

  /**
     @short returns the conditions for a pad placed a given layer,radius for a given sub-detector
  */
  SiPMonTileCharacteristics getConditionsByAlgo(int &layer,
                                                double &radius,
                                                HGCalTypes::TileType &tileType,
                                                HGCalTypes::TileSiPMType &sipmArea);
  SiPMonTileCharacteristicsCore getCoreConditionsByAlgo(int &layer,
                                                        double &radius,
                                                        HGCalTypes::TileType &tileType,
                                                        HGCalTypes::TileSiPMType &sipmArea) {
    return getConditionsByAlgo(layer,radius,tileType,sipmArea).core;
  }
  SiPMonTileCharacteristicsCore getCoreConditionsByAlgo(unsigned int &rawId,double &radius) {
    HGCScintillatorDetId did(rawId);
    int layer=did.layer();
    HGCalTypes::TileType tileType=(HGCalTypes::TileType)did.type();
    HGCalTypes::TileSiPMType sipmArea=(HGCalTypes::TileSiPMType)did.sipm();
    return getConditionsByAlgo(layer,radius,tileType,sipmArea).core;
  }
  

  /**
     @short reconfigures the nPE expected for a given tile type
   */
  void setNpePerMIP(HGCalTypes::TileType tile,float npe);

  /**
     @short reconfigures the reference dark current
   */
  void setReferenceDarkCurrent(double idark);

private:
  
  /**
     @short returns the scale factor to apply for a tile of a given type placed at a given radius
  */
  double scaleByTileArea(const HGCalTypes::TileType &, const double);

  /**
     @short returns the scale factor to apply for a given SiPM area
   */
  double scaleBySipmArea(const HGCalTypes::TileSiPMType &);
    
  //size of the reference scintillator tile
  const double refEdge_;

  //reference dark current for the noise (mA)
  double refDarkCurrent_;
  
  //flags used to disable/override specific SiPM-on-tile operation parameters
  bool ignoreSiPMarea_, ignoreTileArea_, ignoreDoseScale_, ignoreFluenceScale_, ignoreNoise_,
    ignoreTileType_, ignoreAutoPedestalSub_;
    
  //reference values
  std::map<HGCalTypes::TileType, float> nPEperMIP_;
  std::map<HGCalTypes::TileSiPMType, float> maxsipmPE_;  
};

#endif
