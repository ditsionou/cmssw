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
typedef std::array<float, 5> HGCROCTDCChargeDrainParam;


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
     itbx - the in-time index to use
   */
  inline void run(DFr& dataFrame,
                  HGCROCSimHitData& chargeColl,
                  HGCROCSimHitData& toaColl,
                  CLHEP::HepRandomEngine* engine,
                  short itbx = 9) {
    if (opMode_ == TRIVIAL)
      digitizeTrivial(dataFrame, chargeColl, toaColl, itbx);
    else
      digitize(dataFrame, chargeColl, toaColl, engine, itbx);
  }

  /**
     @short runs a trivial digitization routine, the dataframe is filled as in characterization mode
     dataFrame - structure to be filled
     chargeColl - simulated charge per bunch
     toa - time of arrival per bunch
     itbx - the in-time index to use
   */
  void digitizeTrivial(DFr& dataFrame, HGCROCSimHitData& chargeColl, HGCROCSimHitData& toaColl, short itbx);

  /**
     @short runs the digitization routine to fill the dataframe
     switch to time over threshold including deadtime if needed
     the values used for LSBs, dynamical range, pulse shape are the ones configured in the class
     dataFrame - structure to be filled
     chargeColl - simulated charge per bunch
     toaColl - time of arrival per bunch
     engine - a random number generator engine
     itbx - the in-time index to use
  */
  void digitize(DFr& dataFrame,
                HGCROCSimHitData& chargeColl,
                HGCROCSimHitData& toaColl,
                CLHEP::HepRandomEngine* engine,
                short itbx);

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
     @short configure noise components
     note: this should evolve once the common mode noise model is understood
   */
  void configureNoise(float noise) { noise_=noise; }

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
     @short configure time over threshold
     fsc = dynamical range (ns) (triggers recomputation of the LSB)
     chargeDrainParam = parameters to use when evaluating how many bunches TOT is busy
     bxUndershoot = additional bunches where the system undershoots the charge
  */
  void configureTOT(float fsc, float onset, HGCROCTDCChargeDrainParam chargeDrainParam, short bxUndershoot);

  /**
     @short returns charge integration time in TOT mode
   */
  float getChargeIntegrationTime(float charge, CLHEP::HepRandomEngine* engine = nullptr);

  /**
     @short returns the post-TOT charge leakage to a neighboring bunch
     deltaT2nextBx = the deltaT (ns) to the bunch crossing of interest
  */
  float getToTOnsetLeakage(float deltaT2nextBx);
  
  /**
     @short getter for configuration
   */
  edm::ParameterSet cfg() const { return myCfg_; }

  /**
     @short DTOR
   */
  ~HGCROCEmulator() {}

private:

  /**
     @short a simple implementation of the time of arrival based on a smearing of the MC truth
     @params the array of collected charge, time of arrival simulated, the random number generator and the index for the in-time-bunch
  */
  float measureToA(HGCROCSimHitData& chargeColl,
                   HGCROCSimHitData& toaColl,
                   CLHEP::HepRandomEngine* engine,
                   short itbx);
  
  HGCROCOperationMode opMode_;
  uint16_t adcNbits_, adcMax_, toaNbits_, toaMax_, totNbits_, totMax_, totBxUndershoot_;
  float adcLSB_, adcFSC_, toaLSB_, toaFSC_, toaOnset_, totLSB_, totFSC_, totOnset_;
  float noise_;
  float toaJitter_, toaClockOffset_;
  HGCROCPreampPulseShape adcPulse_;
  HGCROCTDCChargeDrainParam totChargeDrainParam_;
  edm::ParameterSet myCfg_;

  //caches
  HGCROCSimHitFlags busyFlags_, undershootFlags_, totFlags_, toaFlags_;
  HGCROCSimHitData newCharge_;
};

#endif
