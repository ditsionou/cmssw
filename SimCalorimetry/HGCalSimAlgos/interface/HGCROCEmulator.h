#ifndef SimCalorimetry_HGCalSimAlgos_HGCROCEmulator_h
#define SimCalorimetry_HGCalSimAlgos_HGCROCEmulator_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/HGCalDigi/interface/HGCalDigiCollections.h"
#include "SimCalorimetry/HGCalSimProducers/interface/HGCDigitizerTypes.h"

#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Random/RandFlat.h"

typedef std::array<float, 15> HGCROCSimHitData;
typedef std::array<bool, 15> HGCROCSimHitFlags;
typedef std::array<float, 6> HGCROCPreampPulseShape;
typedef std::array<float, 2> HGCROCTOAJitterParam;
typedef std::array<float, 3> HGCROCTDCChargeDrainParam;
typedef std::array<float, 2> HGCROCTDCChargeDrainJitterParam;


/**
   @class HGCROCEmulator
   @short models the behavior of the front-end electronics
 */
template <class DFr>
class HGCROCEmulator {
public:

  enum HGCROCOperationMode { DEFAULT, CHARACTERIZATION, TRIVIAL };

  /**
     @short CTOR
   */
  HGCROCEmulator(const edm::ParameterSet& ps);

  /**
     @short 
     runs the emulation according to the operation mode
     dataFrame - structure to be filled
     chargeColl - simulated charge per bunch
     toaColl - time of arrival per bunch
     engine - a random number generator engine
     addNoise - include noise generation
     itbx - the in-time index to use
   */
  inline void run(DFr& dataFrame,
                  HGCROCSimHitData& chargeColl,
                  HGCROCSimHitData& toaColl,
                  CLHEP::HepRandomEngine* engine,
                  bool addNoise,
                  short itbx = 9) {
    if (opMode_ == TRIVIAL)
      digitizeTrivial(dataFrame, chargeColl, toaColl, engine, addNoise, itbx);
    else
      digitize(dataFrame, chargeColl, toaColl, engine, addNoise, itbx);
  }

  /**
     @short runs a trivial digitization routine ignoring out-of-time pileup, pulse shapes etc.
     the floating numbers are simply digitized according to the LSB and FSC
     dataFrame - structure to be filled
     chargeColl - simulated charge per bunch
     toa - time of arrival per bunch
     engine - used to draw random numbers
     addNoise - by default no noise is added
     itbx - the in-time index to use (9 by default)
   */
  void digitizeTrivial(DFr& dataFrame, HGCROCSimHitData& chargeColl, HGCROCSimHitData& toaColl, CLHEP::HepRandomEngine* engine, bool addNoise=false,short itbx=9);

  /**
     @short runs the digitization routine to fill the dataframe
     switch to time over threshold including deadtime if needed
     the values used for LSBs, dynamical range, pulse shape are the ones configured in the class
     dataFrame - structure to be filled
     chargeColl - simulated charge per bunch
     toaColl - time of arrival per bunch
     engine - a random number generator engine
     addNoise - by default noise is added
     itbx - the in-time index to use (9 by default)
  */
  void digitize(DFr& dataFrame,
                HGCROCSimHitData& chargeColl,
                HGCROCSimHitData& toaColl,
                CLHEP::HepRandomEngine* engine,
                bool addNoise=true,
                short itbx=9);

  /**
     @short returns the current charge estimations including leakage and busy state effects
   */
  const HGCROCSimHitData &currentCharges() { return newCharge_; }

  /**
     @short returns the current TOT integration time estimations
   */
  const HGCROCSimHitData &currentIntegrationTimes() { return integTime_; }

  /**
     @short returns the generated noise vector
  */
  const HGCROCSimHitData &currentNoise() { return noiseCharge_; }
  

  /**
     @short operation mode
   */
  HGCROCOperationMode opMode() const { return opMode_; }
  
  /**
     @short returns least-significant bit of the ADC
   */
  float adcLSB() const { return adcLSB_; }
 
  /**
     @short returns full scale charge (dyn. range) of the ADC
  */
  float adcFSC() const { return adcFSC_; }

  /**
     @short return number of bits in the ADC
   */
  uint16_t adcNbits() const { return adcNbits_; }

  /**
     @short returns the max = 2^adcNbits_
   */
  uint16_t adcMax() const { return adcMax_; }

  /**
     @short returns pulse shape leakage to neighboring bunches
   */
  HGCROCPreampPulseShape adcPulse() const { return adcPulse_; }

  /**
     @short returns least-significant bit of the TDC used for time-of-arrival measurement
   */
  float toaLSB() const { return toaLSB_; }

  /**
     @short returns full scale charge of the TDC used for time-of-arrival
   */
  float toaFSC() const { return toaFSC_; }

  /**
     @short charge onset to trigger the measurement of time-of-arrival
   */
  float toaOnset() const { return toaOnset_; }

  /**
     @short number of bits of the TDC used to measure time-of-arrival
   */
  uint16_t toaNbits() const { return toaNbits_; }

  /**
     @short returns max value of the TDC used to measure time-of-arrival = 2^toaNbits_
   */
  uint16_t toaMax() const { return toaMax_; }
  
  /**
     @short returns least-significant bit of the TDC used for time-over-threshold
   */
  float totLSB() const { return totLSB_; }

  /**
     @short returns full scale charge of the TDC for time-over-threshold
   */
  float totFSC() const { return totFSC_; }

  /**
     @short returns number of bits in the TDC for time-over-threshold
   */
  uint16_t totNbits() const { return totNbits_; }

  /**
     @short returns max in TDC for time-over-threshold = 2^totNbits_
   */
  uint16_t totMax() const { return totMax_; }

  /**
     @short returns charge onset to trigger time-over-threshold
   */
  float totOnset() const { return totOnset_; }

  /**
     @short returns the parameters used to model the charge drain after the time-over-threshold period
   */
  HGCROCTDCChargeDrainParam totChargeDrainParam() const { return totChargeDrainParam_; }

  /**
     @short returns the parameters use to introduce a jitter in the charge drain after time-over-threshold
   */
  HGCROCTDCChargeDrainJitterParam totChargeDrainJitterParam() const {  return totChargeDrainJitterParam_; }
  
  /**
     @short configure operation mode
  */
  void configureMode(HGCROCOperationMode mode);

  /**
     @short configure ADC
     fsc = full scale charge (triggers recomputation of the LSB)
     shape = pulse shape to use
  */
  void configureADC(float fsc, HGCROCPreampPulseShape shape);

  /**
     @short configure time of arrival
     fsc = dynamical range (ns) (triggers recomputation of the LSB)
     onset = charge after which toa is computed
     toa{jitter,clkoff} = {stochastic, constant} for the smearing of the toa
  */
  void configureTOA(float fsc, float onset, float toajitter, float toaclkoff);

  /**
     @short configures the noise parameters
     ped - pedestal
     jitter - stochastic component
     cm - common mode
   */
  void configureNoise(float ped, float jitter,float cm);
  
  /**
     @short configure time over threshold
     fsc = dynamical range (ns) (triggers recomputation of the LSB)
     chargeDrainParam = parameters to use when evaluating how many bunches TOT is busy
     bxUndershoot = additional bunches where the system undershoots the charge
  */
  void configureTOT(float fsc,
                    float onset,
                    HGCROCTDCChargeDrainParam chargeDrainParam,
                    HGCROCTDCChargeDrainJitterParam chargeDrainJitterParam,
                    short bxUndershoot);

  /**
     @short returns charge integration time, in units of bunch crossings, when in TOT mode
  */
  int16_t getChargeIntegrationTime(float charge, CLHEP::HepRandomEngine* engine = nullptr);
  
  /**
     @short getter for configuration
   */
  edm::ParameterSet cfg() const { return myCfg_; }

  /**
     @short evaluates the expected parallel component to the noise for a given gain and sensor capacitance
   */
  float getENCs(float gain,float cap);

  /**
     @short evalutes the expected series component to the noise for a given gain and leakage current
   */
  float getENCp(float gain,float ileak);
  
  /**
     @short DTOR
   */
  ~HGCROCEmulator() {}

private:

  /**
     @short generates the noise bunch-by-bunch
   */
  void generateNoise(CLHEP::HepRandomEngine* engine);
  
  /**
     @short generates a common mode word
     At the moment this is a rather simple call to the noise smearing including a pedestal and 
     common mode value configured.
   */
  float generateCMWord(CLHEP::HepRandomEngine* engine);
  
  /**
     @short a simple implementation of the time of arrival based on a smearing of the MC truth
     @params the array of collected charge, time of arrival simulated, the random number generator and the index for the in-time-bunch
  */
  float measureToA(HGCROCSimHitData& chargeColl,
                   HGCROCSimHitData& toaColl,
                   CLHEP::HepRandomEngine* engine,
                   short itbx);

  /**
     @short determines which bunches will trigger ToT
     The algorithm is iterative and integration time needs to be updated depending on the charge which is accumulated while in ToT mode or which leaks from previous bunches
     The cache for newCharge_ will be updated with the final measurements
   */
  void measureChargeWithTOT(HGCROCSimHitData& chargeColl);
  
  /**
     @short convolve the pulse shape with the charge accumulated in each bunch
     The cache for newCharge_ will be updated with the final measurements
   */
  void measureChargeWithADCPreamp(HGCROCSimHitData& chargeColl);

  /**
     @short estimates leakage to the bunch crossing it using the configured pre-amp charge
     All busy/tot flags which have been set will be used to skip bunch crossings in those states
   */
  float estimateLeakage(HGCROCSimHitData& chargeColl, size_t it);
  
  /**
     @short resets the caches used in the digitization
     it's also used to trigger the generation of noise
  */
  void resetCaches(CLHEP::HepRandomEngine* engine);

  //configuration parameters
  HGCROCOperationMode opMode_;
  uint16_t adcNbits_, adcMax_, toaNbits_, toaMax_, totNbits_, totMax_, totBxUndershoot_;
  float adcLSB_, adcFSC_, toaLSB_, toaFSC_, toaOnset_, totLSB_, totFSC_, totOnset_;
  float pedestal_, noiseJitter_, commonNoise_;
  float toaJitter_, toaClockOffset_;
  HGCROCPreampPulseShape adcPulse_;
  HGCROCTDCChargeDrainParam totChargeDrainParam_;
  HGCROCTDCChargeDrainJitterParam totChargeDrainJitterParam_;
  edm::ParameterSet myCfg_;

  //caches
  HGCROCSimHitFlags busyFlags_, totFlags_, toaFlags_;
  HGCROCSimHitData noiseCharge_, newCharge_,integTime_;
};

#endif
