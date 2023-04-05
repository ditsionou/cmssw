#include "SimCalorimetry/HGCalSimAlgos/interface/HGCalConfigurationByAlgoWrapper.h"
#include "DataFormats/HGCalDigi/interface/HGCalDigiCollections.h"

//
template<>
uint32_t HGCalSiGeoConfigurationByAlgo::toConfigurableKey(HGCSiliconDetId &d) {

  //reset cell (u,v) in detId and return new value
  HGCSiliconDetId modId(d.det(),d.zside(),d.type(),d.layer(),d.waferU(),d.waferV(),0,0);
  return modId.rawId();

}





//LEFT OVERS START HERE
//#include "SimCalorimetry/HGCalSimProducers/interface/HGCFEElectronics.h"
//#include "DataFormats/ForwardDetId/interface/HGCSiliconDetId.h"
//#include "DataFormats/ForwardDetId/interface/HFNoseDetId.h"
//#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
//#include <string>
//#include <array>
//#include <unordered_map>
//
//
//          enum SiPMGainRange_t { GAIN_2, GAIN_4, AUTO, GAINRANGE_N };  //roc gain for 2mm2 and 4mm2
//
//    SiPMonTileCharacteristics() : s(0.), lySF(0.), n(0.), xtalk(0), gain(0), thrADC(0), ntotalPE(0) {}
//    float s, lySF, n, xtalk;
//    unsigned short gain, thrADC, ntotalPE;
//  };
//
//  /**
//     @short returns the signal scaling and the noise
//  */
//  SiPMonTileCharacteristics scaleByDose(const HGCScintillatorDetId &,
//                                        const double,
//                                        const int aimMIPtoADC = 15,
//                                        const GainRange_t gainPreChoice = GainRange_t::AUTO);
//
//
//  void setSipmMap(const std::string &);
//  void setReferenceDarkCurrent(double idark);
//  void setReferenceCrossTalk(double xtalk) { refXtalk_ = xtalk; }
//  void setNpePerMIP(float npePerMIP);
//  std::array<double, GAINRANGE_N> getLSBPerGain() { return lsbPerGain_; }
//  std::array<double, GAINRANGE_N> getMaxADCPerGain() { return fscADCPerGain_; }
//  std::array<double, TILETYPE_N> getNpePerMIP() { return nPEperMIP_; }
//  float getNPeInSiPM() { return maxSiPMPE_; }
//  bool ignoreAutoPedestalSubtraction() { return ignoreAutoPedestalSub_; }
//
//  
//  /**
//     @short parses the radius boundaries for the SiPM area assignment from a custom file
//   */
//  std::unordered_map<int, float> readSipmPars(const std::string &);
//
//  
//  //lsb and fsc per gain
//  std::array<double, GAINRANGE_N> lsbPerGain_, fscADCPerGain_;
//
//  
//  //reference cross talk parameter (0,1) - if -1 it will be used to ignore effect in the digitization step
//  double refXtalk_;
//
//  //reference ADC counts for the MIP peak
//  int aimMIPtoADC_;
//
//  //sipm size boundaries
//  std::unordered_map<int, float> sipmMap_;
//      aimMIPtoADC_(15),
//
// //full scale charge per gain in nPE
//  //this is chosen for now such that the ref. MIP peak is at N ADC counts
//  //to be changed once the specs are fully defined such that the algorithm chooses the best gain
//  //for the MIP peak to be close to N ADC counts (similar to Si) or applying other specific criteria
//  fscADCPerGain_[GAIN_2] = nPEperMIP_[CAST] * 1024. / aimMIPtoADC_;      //2 mm^2  SiPM
//  fscADCPerGain_[GAIN_4] = 2 * nPEperMIP_[CAST] * 1024. / aimMIPtoADC_;  //4mm^2   SiPM
//
//  //lsb: adc has 10 bits -> 1024 counts at max
//  for (size_t i = 0; i < GAINRANGE_N; i++)
//    lsbPerGain_[i] = fscADCPerGain_[i] / 1024.f;
//}
//
////
//void HGCalSiPMonTileConditionsByAlgo::setSipmMap(const std::string& fullpath) { sipmMap_ = readSipmPars(fullpath); }
//
////
//std::unordered_map<int, float> HGCalSiPMonTileConditionsByAlgo::readSipmPars(const std::string& fullpath) {
//  std::unordered_map<int, float> result;
//  //no file means default sipm size
//  if (fullpath.empty())
//    return result;
//
//  edm::FileInPath fp(fullpath);
//  std::ifstream infile(fp.fullPath());
//  if (!infile.is_open()) {
//    throw cms::Exception("FileNotFound") << "Unable to open '" << fullpath << "'" << std::endl;
//  }
//  std::string line;
//  while (getline(infile, line)) {
//    int layer;
//    float boundary;
//
//    //space-separated
//    std::stringstream linestream(line);
//    linestream >> layer >> boundary;
//
//    result[layer] = boundary;
//  }
//  return result;
//}
//
//
//                                                                          int aimMIPtoADC,
//                                                                          GainRange_t gainPreChoice)
//  HGCalSiPMonTileConditionsByAlgo::GainRange_t gain = sipm.second;
//
//  sipmChar.thrADC = std::floor(0.5 * S / lsbPerGain_[gain]);
//    sipmChar.xtalk = refXtalk_;
//
////
//std::pair<double, HGCalSiPMonTileConditionsByAlgo::GainRange_t> HGCalSiPMonTileConditionsByAlgo::scaleBySipmArea(
//    const HGCScintillatorDetId& cellId, const double radius, const HGCalSiPMonTileConditionsByAlgo::GainRange_t& gainPreChoice) {
//  //start with the prechosen gain
//  //if auto then override it according to the SiPM area
//  HGCalSiPMonTileConditionsByAlgo::GainRange_t gain(gainPreChoice);
//  if (gainPreChoice == HGCalSiPMonTileConditionsByAlgo::GainRange_t::AUTO)
//    gain = GainRange_t::GAIN_2;
//
//  double scaleFactor(1.f);
//
//  if (ignoreSiPMarea_)
//    return std::pair<double, HGCalSiPMonTileConditionsByAlgo::GainRange_t>(scaleFactor, gain);
//
//  //use sipm area boundary map
//  if (overrideSiPMarea_) {
//    int layer = cellId.layer();
//    if (sipmMap_.count(layer) > 0 && radius < sipmMap_[layer]) {
//      scaleFactor = 2.f;
//      if (gainPreChoice == HGCalSiPMonTileConditionsByAlgo::GainRange_t::AUTO)
//        gain = GainRange_t::GAIN_4;
//    }
//  }
//  //read from DetId
//  else {
//    int sipm = cellId.sipm();
//    if (sipm == 0) {
//      scaleFactor = 2.f;
//      if (gainPreChoice == HGCalSiPMonTileConditionsByAlgo::GainRange_t::AUTO)
//        gain = GainRange_t::GAIN_4;
//    }
//  }
//
//  return std::pair<double, HGCalSiPMonTileConditionsByAlgo::GainRange_t>(scaleFactor, gain);
//}
