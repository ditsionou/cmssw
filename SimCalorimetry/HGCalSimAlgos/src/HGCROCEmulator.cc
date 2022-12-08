#include "SimCalorimetry/HGCalSimAlgos/interface/HGCROCEmulator.h"
#include "FWCore/Utilities/interface/transform.h"
#include "vdt/vdtMath.h"

//
template <class DFr>
HGCROCEmulator<DFr>::HGCROCEmulator(const edm::ParameterSet& ps) {

  myCfg_=ps;

  //operation mode
  configureMode((HGCROCOperationMode)ps.getParameter<uint32_t>("opMode"));

  //adc configuration
  adcNbits_ = 10;
  adcMax_ = pow(2, adcNbits_);
  adcFSC_ = ps.getParameter<double>("adcFSC");
  auto adcPulseCfg = ps.getParameter<std::vector<double> >("adcPulse");
  for (size_t i = 0; i < adcPulseCfg.size(); ++i) {
    adcPulse_[i] = (float)adcPulseCfg[i];
  }
  configureADC(adcFSC_, adcPulse_);

  //tot configuration
  totNbits_ = 12;
  totMax_ = pow(2, totNbits_);
  totFSC_ = ps.getParameter<double>("totFSC");
  totOnset_ = ps.getParameter<double>("totOnset");
  totBxUndershoot_ = ps.getParameter<uint32_t>("totBxUndershoot");
  auto totChargeDrainParamCfg = ps.getParameter<std::vector<double> >("totChargeDrainParam");
  for (size_t i = 0; i < totChargeDrainParamCfg.size(); i++) {
    totChargeDrainParam_[i] = totChargeDrainParamCfg[i];
  }
  auto totChargeDrainJitterParamCfg = ps.getParameter<std::vector<double> >("totChargeDrainJitterParam");
  for (size_t i = 0; i < totChargeDrainJitterParamCfg.size(); i++) {
    totChargeDrainJitterParam_[i] = totChargeDrainJitterParamCfg[i];
  }
  configureTOT(totFSC_, totOnset_, totChargeDrainParam_, totChargeDrainJitterParam_, totBxUndershoot_);

  //toa configuration
  toaNbits_ = 10;
  toaMax_ = pow(2, toaNbits_);
  toaFSC_ = ps.getParameter<double>("toaFSC");
  toaOnset_ = ps.getParameter<double>("toaOnset");
  toaJitter_ = ps.getParameter<double>("toaJitter");
  toaClockOffset_ = ps.getParameter<double>("toaClockOffset");
  configureTOA(toaFSC_, toaOnset_, toaJitter_, toaClockOffset_);
}

//
template <class DFr>
void HGCROCEmulator<DFr>::configureMode(HGCROCOperationMode mode) {
  opMode_ = mode;
}

//
template <class DFr>
void HGCROCEmulator<DFr>::configureADC(float fsc, HGCROCPreampPulseShape shape) {
  adcFSC_ = fsc;
  adcLSB_ = (adcFSC_ / adcMax_);
  adcPulse_ = shape;
}

//
template <class DFr>
void HGCROCEmulator<DFr>::configureTOA(float fsc, float onset, float toajitter, float toaclkoff) {
  toaFSC_ = fsc;
  toaLSB_ = (toaFSC_ / toaMax_);
  toaOnset_ = onset;
  toaJitter_ = toajitter;
  toaClockOffset_ = toaclkoff;
}

//
template <class DFr>
void HGCROCEmulator<DFr>::configureTOT(float fsc,
                                       float onset,
                                       HGCROCTDCChargeDrainParam chargeDrainParam,
                                       HGCROCTDCChargeDrainJitterParam chargeDrainJitterParam,
                                       short bxUndershoot) {
  totFSC_ = fsc;
  totLSB_ = (totFSC_ / totMax_);
  totOnset_ = onset;
  totBxUndershoot_ = bxUndershoot;
  totChargeDrainParam_ = chargeDrainParam;
  totChargeDrainJitterParam_ = chargeDrainJitterParam;
}

//
template <class DFr>
int16_t HGCROCEmulator<DFr>::getChargeIntegrationTime(float charge, CLHEP::HepRandomEngine* engine) {

  //average busy time
  float isq_charge=vdt::fast_isqrtf(charge);  
  float busyBx = totChargeDrainParam_[0]*charge + \
    totChargeDrainParam_[1] / isq_charge + \
    totChargeDrainParam_[2];
  
  //smear time if an engine has been passed
  if (engine) {
    float busySigma=totChargeDrainJitterParam_[0]*charge + totChargeDrainJitterParam_[1];
    busyBx = CLHEP::RandGaussQ::shoot(engine, busyBx, busySigma);
  }

  busyBx += totBxUndershoot_;
  
  return busyBx;
};

//
template <class DFr>
void HGCROCEmulator<DFr>::digitizeTrivial(DFr& dataFrame,
                                          HGCROCSimHitData& chargeColl,
                                          HGCROCSimHitData& toaColl,
                                          short itbx) {
#ifdef EDM_ML_DEBUG
  edm::LogVerbatim("HGCROCEmulator::digitizeTrivial") << "[digitizeTrivial]" << std::endl;
#endif

  resetCaches();
  
  //check if tot is to be triggered based on the charge
  //digitize charge in ADC and TOT modes and time of arrival
  float charge(chargeColl[itbx]);
  bool tp(false), tc(false);
  if (charge > totOnset_) {
    tc=true;
    tp=true;
  }
  
  uint32_t adcm1(0);
  uint32_t adc = std::floor(std::min(charge, adcFSC_ - adcLSB_) / adcLSB_);
  uint32_t tot = tc ? std::floor(std::min(charge, totFSC_ - totLSB_) / totLSB_) : 0;
  uint32_t toa = charge > toaOnset_ ? std::floor(std::min(toaColl[itbx], toaFSC_ - toaLSB_) / toaLSB_) : 0;

  //fill the dataframe
  dataFrame.fill(opMode_ == HGCROCOperationMode::CHARACTERIZATION, tc, tp, adcm1, adc, tot, toa);

#ifdef EDM_ML_DEBUG
    std::ostringstream msg;
    dataFrame.print(msg);
    edm::LogVerbatim("HGCROCEmulator::digitizeTrivial") << msg.str() << std::endl;
#endif
    
}

//
template <class DFr>
float HGCROCEmulator<DFr>::measureToA(HGCROCSimHitData& chargeColl,
                                      HGCROCSimHitData& toaColl,
                                      CLHEP::HepRandomEngine* engine,
                                      short itbx
                                      ) {
  
  //the simulation is simplified for the moment
  //in all cases look only at in-time signals: ToA will be computed for the central BX only
  //to be done properly with realistic ToA shaper and jitter: for the moment accounted in the smearing
  toaFlags_.fill(false);
  float timeToA = 0.f;
  if (toaColl[itbx] != 0.f) {

    timeToA = toaColl[itbx];    

    //emulate stochastic term ~ A / (S/N)
    if(noise_>0) {
      float sovern = chargeColl[itbx] / noise_;
      float jitter = toaJitter_ / sovern;
      timeToA = CLHEP::RandGaussQ::shoot(engine, timeToA, jitter);
    }

    //add a clock jump (constant term)
    timeToA += toaClockOffset_;
    
    if (timeToA >= 0.f && timeToA <= toaFSC_)
      toaFlags_[itbx] = true;
  }

  return timeToA;
}

//
template <class DFr>
void HGCROCEmulator<DFr>::digitize(DFr& dataFrame,
                                   HGCROCSimHitData& chargeColl,
                                   HGCROCSimHitData& toaColl,
                                   CLHEP::HepRandomEngine* engine,
                                   short itbx) {

  resetCaches();
  
  //get time-of-arrival
  float timeToA = this->measureToA(chargeColl,toaColl,engine,itbx);

  //TOT charge measurements must be estimated first
  this->measureChargeWithTOT(chargeColl);

  //ADC charge measurement for non-busy and non-ToT bunches 
  this->measureChargeWithADCPreamp(chargeColl);
  
  //fill the ROC dataframe for this channel

  //Tc/Tp flags
  bool tc(totFlags_[itbx] && !busyFlags_[itbx]);
  bool tp(totFlags_[itbx] || busyFlags_[itbx]);

  //ADC BX-1 is only non zero if it was not busy
  //TOT mode is not considered for previous bunch, so ADC will be saturated if TOT was triggered
  uint16_t adcm1(0);
  int itbxm1(itbx - 1);
  if (itbxm1 >= 0 && !busyFlags_[itbxm1]) {
    adcm1 = std::floor(std::min(newCharge_[itbxm1], adcFSC_ - adcLSB_) / adcLSB_);
  }

  //ADC, TOT, TOA in-time BX
  //note: in principle the TOT has an intrinsic offset but this is ignored for the moment
  uint16_t adc(0), toa(0), tot(0);
  if (!busyFlags_[itbx]) {
    adc = std::floor(std::min(newCharge_[itbx], adcFSC_ - adcLSB_) / adcLSB_);
    if (totFlags_[itbx] || opMode_ == HGCROCOperationMode::CHARACTERIZATION)
      tot = std::floor(std::min(newCharge_[itbx], totFSC_ - totLSB_) / totLSB_);
    toa = newCharge_[itbx] > toaOnset_ ? std::floor(std::min(timeToA, toaFSC_ - toaLSB_) / toaLSB_) : 0;
  }

  dataFrame.fill(opMode_ == HGCROCOperationMode::CHARACTERIZATION, tc, tp, adcm1, adc, tot, toa);
  
#ifdef EDM_ML_DEBUG
  std::ostringstream msg;
  dataFrame.print(msg);
  edm::LogVerbatim("HGCROCEmulator") << msg.str() << std::endl;
#endif  

}

//
template <class DFr>
void HGCROCEmulator<DFr>::measureChargeWithTOT(HGCROCSimHitData& chargeColl) {
  
  //charge measurement (requires identifying bunches which will trigger ToT and then run charge sharing)
  //reset the caches as the algorithm requires some iterations over the initial values
  for (size_t it = 0; it < chargeColl.size(); ++it) {

    //if already flagged as busy it can't be re-used to trigger the ToT
    if (busyFlags_[it]) continue;
    
    float charge = chargeColl[it];
    
    //if below TDC onset will be handled by the ADC later
    if (charge < totOnset_) continue;
    
    //raise TDC mode for charge computation
    totFlags_[it] = true;
    
#ifdef EDM_ML_DEBUG
    edm::LogVerbatim("HGCROCEmulator") << "HGCROCEmulator::digitize"
      << "\t q=" << charge << " fC with triggers ToT @ " << it << std::endl;
#endif
    
    //compute total charge to be integrated and integration time
    //needs a loop as ToT will last as long as there is charge to dissipate
    int16_t busyBxs(0);
    float totalCharge(charge);
    bool addLeak(true);
    while (true) {
      
      //compute integration time in #bunches
      //if no update is needed regarding the number of bunches, then the ToT integration time has converged
      //and the infinite while is made collapse
      const int16_t newBusyBxs = getChargeIntegrationTime(totalCharge);
      if (newBusyBxs == busyBxs) break;
#ifdef EDM_ML_DEBUG
      if (busyBxs == 0)
        edm::LogVerbatim("HGCROCEmulator") << "\t Intial busy estimate=" << newBusyBxs << " bxs" << std::endl;
      else
        edm::LogVerbatim("HGCROCEmulator")
          << "\t ...updating initial busy estimate to " << newBusyBxs << " bxs, iterating again" << std::endl;
#endif
    
      //update number of busy bunches
      busyBxs = newBusyBxs;
      
      //recompute charge to be integrated adding leakage from previous bunches in SARS ADC mode
      //this only needs to be done once
      if(addLeak) {
        for (size_t jt = 0; jt < it; ++jt) {
          const size_t deltaT = (it - jt);
          if ((deltaT + 2) >= adcPulse_.size() || chargeColl[jt] == 0.f || totFlags_[jt] || busyFlags_[jt])
            continue;
          totalCharge += chargeColl[jt] * adcPulse_[deltaT + 2];
        }
#ifdef EDM_ML_DEBUG
        edm::LogVerbatim("HGCROCEmulator") << "\t\t adding leakage charge new estimate is " << totalCharge << "fC" << std::endl;
#endif
        addLeak=false;
      }

      //add contamination from posterior bunches if they haven't yet been flagged as busy in a previous iteration
      for (size_t jt = it + 1; jt < it + busyBxs && jt < chargeColl.size(); ++jt) {
        if(busyFlags_[jt]) continue;
        busyFlags_[jt] = true;
        totalCharge += chargeColl[jt];
#ifdef EDM_ML_DEBUG
        edm::LogVerbatim("HGCROCEmulator") << "\t\t adding charge @ bx +" << (jt-it) << " new estimate is" << totalCharge << "fC" << std::endl;
#endif
      } 
    }

    //final integrated charge in ToT
    newCharge_[it] = totalCharge;
#ifdef EDM_ML_DEBUG
    std::cout << "HGCROCEmulator::digitize"
      << "\t Final busy estimate=" << busyBxs << " bxs" << std::endl
      << "\t Total integrated=" << totalCharge << " fC" << std::endl;
#endif
  }
  
}


//
template <class DFr>
void HGCROCEmulator<DFr>::measureChargeWithADCPreamp(HGCROCSimHitData& chargeColl) {

  //convolve the pulse shape with the charge accumulated in each bunch
  int ipulse = 0;
  for (int it = 0; it < (int)(chargeColl.size()); ++it) {

    //if busy, charge has been already integrated
    if (!totFlags_[it] & !busyFlags_[it]) {
      const int start = std::max(0, 2 - it);
      const int stop = std::min((int)adcPulse_.size(), (int)newCharge_.size() - it + 2);
      for (ipulse = start; ipulse < stop; ++ipulse) {
        const int itoffset = it + ipulse - 2;
        //notice that if the channel is already busy,
        //it has already been affected by the leakage of the SARS ADC
        if (!totFlags_[itoffset] & !busyFlags_[itoffset]) {
          newCharge_[itoffset] += chargeColl[it] * adcPulse_[ipulse];
        }
      }
    }
  }
};



//
template <class DFr>
void HGCROCEmulator<DFr>::resetCaches() {
  toaFlags_.fill(false);
  busyFlags_.fill(false);
  totFlags_.fill(false);
  newCharge_.fill(0.f);  
}


// trigger the compiler to generate the appropriate code
#include "DataFormats/HGCalDigi/interface/HGCalDigiCollections.h"
template class HGCROCEmulator<HGCROCChannelDataFrameSpec>;
