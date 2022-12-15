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

  //noise parameters
  pedestal_ = ps.getParameter<double>("pedestal");
  noiseJitter_ = ps.getParameter<double>("noiseJitter");
  commonNoise_ = ps.getParameter<double>("commonNoise");
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
void HGCROCEmulator<DFr>::configureNoise(float ped, float jitter, float cm) {
  pedestal_=ped;
  noiseJitter_=jitter;
  commonNoise_=cm;
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

  //block additional bunches where the charge is expected to undershoot the baseline
  busyBx += totBxUndershoot_;
  
  return busyBx;
};

//
template <class DFr>
void HGCROCEmulator<DFr>::digitizeTrivial(DFr& dataFrame,
                                          HGCROCSimHitData& chargeColl,
                                          HGCROCSimHitData& toaColl,
                                          CLHEP::HepRandomEngine* engine,
                                          bool addNoise,
                                          short itbx) {
#ifdef EDM_ML_DEBUG
  edm::LogVerbatim("HGCROCEmulator::digitizeTrivial") << "[digitizeTrivial]" << std::endl;
#endif

  resetCaches(engine);
  
  //check if tot is to be triggered based on the charge
  //digitize charge in ADC and TOT modes and time of arrival
  float charge(chargeColl[itbx] + addNoise*noiseCharge_[itbx]);
  bool tp(false), tc(false);
  if (charge > totOnset_) {
    tc=true;
    tp=true;
  }
  
  uint32_t adcm1(0);
  uint32_t adc = std::floor(std::min(charge, adcFSC_ - adcLSB_) / adcLSB_);
  uint32_t tot = tc ? std::floor(std::min(charge, totFSC_ - totLSB_) / totLSB_) : 0;
  uint32_t digiToA(std::floor(toaColl[itbx]/toaLSB_));
  uint32_t toa = charge > toaOnset_ ? digiToA%1024 : 0;

  //update caches (only in-time bunch as this is trivial)
  newCharge_[itbx] = charge;
  integTime_[itbx] = float(tc ? getChargeIntegrationTime(charge) : 0.);
  totFlags_[itbx] = true;
  
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
    if(noiseJitter_>0) {
      float sovern = chargeColl[itbx] / noiseJitter_;
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
void HGCROCEmulator<DFr>::generateNoise(CLHEP::HepRandomEngine* engine) {
  for(size_t i=0; i<noiseCharge_.size(); i++) {
    noiseCharge_[i] = CLHEP::RandGaussQ::shoot(engine, pedestal_, noiseJitter_) + commonNoise_;
  }
}

//
template <class DFr>
float HGCROCEmulator<DFr>::estimateLeakage(HGCROCSimHitData& chargeColl, size_t it) {

  //recompute charge to be integrated adding leakage from previous bunches in SARS ADC mode
  float leakCharge(0.f);
  for (size_t jt = 0; jt < it; ++jt) {
    if(chargeColl[jt] == 0.f || totFlags_[jt] || busyFlags_[jt]) continue;
    const size_t deltaT = (it - jt);
    if ((deltaT + 2) >= adcPulse_.size()) continue;
    leakCharge += chargeColl[jt] * adcPulse_[deltaT + 2];
  }
#ifdef EDM_ML_DEBUG
  edm::LogVerbatim("HGCROCEmulator::estimateLeakage") << "Leakage estimated is " << leakCharge << "fC" << std::endl;
#endif
  
  return leakCharge;
}



//
template <class DFr>
void HGCROCEmulator<DFr>::digitize(DFr& dataFrame,
                                   HGCROCSimHitData& simChargeColl,
                                   HGCROCSimHitData& toaColl,
                                   CLHEP::HepRandomEngine* engine,                                   
                                   bool addNoise,
                                   short itbx) {

  resetCaches(engine);

  //add noise
  HGCROCSimHitData chargeColl;
  for(size_t i=0; i<chargeColl.size(); i++)
    chargeColl[i] = simChargeColl[i] + addNoise*noiseCharge_[i];
  
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

    //adc and tot saturate
    adc = std::floor(std::min(newCharge_[itbx], adcFSC_ - adcLSB_) / adcLSB_);
    if (totFlags_[itbx] || opMode_ == HGCROCOperationMode::CHARACTERIZATION)
      tot = std::floor(std::min(newCharge_[itbx], totFSC_ - totLSB_) / totLSB_);

    //toa cycles...
    uint32_t digiToA(std::floor(timeToA/toaLSB_));
    toa = newCharge_[itbx] > toaOnset_ ? digiToA%1024 : 0;
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

    //first estimate of the charge only takes into account leakage charge and the in-time charge
    float leakCharge=this->estimateLeakage(chargeColl,it);
    float charge=chargeColl[it] + leakCharge;
    
    //if below TDC onset, this will be handled by the ADC later
    if (charge < totOnset_) continue;
        
#ifdef EDM_ML_DEBUG
    edm::LogVerbatim("HGCROCEmulator::digitize")
      << "\t q=" << charge << " fC (leak="<< leakCharge << " fC) "
      << " triggers ToT @ " << it << "-th bunch" << std::endl;
#endif
    
    //compute total charge to be integrated and integration time
    //raise TOT mode for charge computation
    totFlags_[it] = true;
    int16_t busyBxs = getChargeIntegrationTime(charge);    
#ifdef EDM_ML_DEBUG
    edm::LogVerbatim("HGCROCEmulator") << "\t Intial busy estimate=" << newBusyBxs << " bxs" << std::endl;    
#endif

    //add contamination from posterior bunches
    //if they haven't yet been flagged as busy in a previous iteration
    //the loop ends when the charge integration time has stabilized
    while(true) {
      
      for (size_t jt = it + 1; jt < it + busyBxs && jt < chargeColl.size(); ++jt) {
        if(busyFlags_[jt]) continue;
        busyFlags_[jt] = true;
        charge += chargeColl[jt];
#ifdef EDM_ML_DEBUG
        edm::LogVerbatim("HGCROCEmulator")
          << "\t\t adding charge @ bx +" << (jt-it)
          << " new estimate is " << charge << " fC" << std::endl;
#endif
      }

      //check if this has converged
      uint16_t updatedBusyBxs = getChargeIntegrationTime(charge);
      if(updatedBusyBxs==busyBxs) break;
      busyBxs=updatedBusyBxs;
    }
    
    //final integrated charge in ToT
    newCharge_[it] = charge;
    integTime_[it] = busyBxs;

#ifdef EDM_ML_DEBUG
    edm::LogVerbatim("HGCROCEmulator")
      << "\t Final busy estimate=" << busyBxs << " bxs" << std::endl
      << "\t Total integrated=" << charge << " fC" << std::endl;
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
void HGCROCEmulator<DFr>::resetCaches(CLHEP::HepRandomEngine* engine) {
  toaFlags_.fill(false);
  busyFlags_.fill(false);
  totFlags_.fill(false);
  newCharge_.fill(0.f);
  integTime_.fill(0.f);
  generateNoise(engine);
}


//
template <class DFr>
float HGCROCEmulator<DFr>::getENCs(float gain,float cap) {
  if(gain==80)  return 0.000017*pow(cap,2)+0.0033*cap+0.1269;
  if(gain==160) return 0.000017*pow(cap,2)+0.0021*cap+0.1903;
  if(gain==320) return 0.000023*pow(cap,2)+0.0011*cap+0.3451;
  return 0.;
}

//
template <class DFr>
float HGCROCEmulator<DFr>::getENCp(float gain,float ileak) {

  if(gain==80){
    if(ileak>34.22) return 0.00228*pow(ileak,2)+-0.15429*ileak+2.92155;
    else if(ileak>29.57) return 0.00180*pow(ileak,2)+-0.10225*ileak+1.74070;
    else if(ileak>24.68) return 0.00122*pow(ileak,2)+-0.05270*ileak+0.82095;
    else if(ileak>19.79) return 0.00097*pow(ileak,2)+-0.03050*ileak+0.45769;
    else if(ileak>14.58) return 0.00063*pow(ileak,2)+-0.00800*ileak+0.18280;
    else if(ileak>9.83) return 0.00020*pow(ileak,2)+0.01021*ileak+0.04443;
    else if(ileak>4.80) return -0.00023*pow(ileak,2)+0.02150*ileak+0.01669;
    return -0.00149*pow(ileak,2)+0.03648*ileak+0.02429;
  }
  if(gain==160){
    if(ileak>34.22) return 0.00081*pow(ileak,2)+-0.04680*ileak+0.95977;
    else if(ileak>29.47) return 0.00064*pow(ileak,2)+-0.02779*ileak+0.55068;
    else if(ileak>24.63) return 0.00060*pow(ileak,2)+-0.02028*ileak+0.39420;
    else if(ileak>19.75) return 0.00054*pow(ileak,2)+-0.01162*ileak+0.25168;
    else if(ileak>14.63) return 0.00034*pow(ileak,2)+0.00200*ileak+0.10092;
    else if(ileak>9.74) return 0.00016*pow(ileak,2)+0.01155*ileak+0.03837;
    else if(ileak>4.85) return -0.00015*pow(ileak,2)+0.02084*ileak+0.01906;
    return -0.00172*pow(ileak,2)+0.03924*ileak+0.02344;
  }
  if(gain==320){
    if(ileak>34.17) return 0.00031*pow(ileak,2)+-0.01268*ileak+0.36704;
    else if(ileak>29.33) return 0.00026*pow(ileak,2)+-0.00615*ileak+0.23206;
    else if(ileak>24.58) return 0.00021*pow(ileak,2)+-0.00019*ileak+0.12965;
    else if(ileak>19.70) return 0.00015*pow(ileak,2)+0.00501*ileak+0.06792;
    else if(ileak>14.58) return 0.00014*pow(ileak,2)+0.00797*ileak+0.04938;
    else if(ileak>9.74) return -0.00002*pow(ileak,2)+0.01479*ileak+0.01662;
    else if(ileak>4.94) return -0.00021*pow(ileak,2)+0.02039*ileak+0.01977;
    return -0.00166*pow(ileak,2)+0.03746*ileak+0.02337;
  }
  return 0.f;
}

// trigger the compiler to generate the appropriate code
#include "DataFormats/HGCalDigi/interface/HGCalDigiCollections.h"
template class HGCROCEmulator<HGCROCChannelDataFrameSpec>;
