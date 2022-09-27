#include "SimCalorimetry/HGCalSimAlgos/interface/HGCROCAlgos.h"
#include "FWCore/Utilities/interface/transform.h"
#include "vdt/vdtMath.h"

using namespace hgc_digi;

//
template <class DFr>
HGCROCAlgos<DFr>::HGCROCAlgos(const edm::ParameterSet& ps) : verbose_(false) {

  myCfg_=ps;

#ifdef EDM_ML_DEBUG
  verbose_ = true;
#else
  verbose_ = ps.getParameter<bool>("verbose");
#endif

  //operation mode
  configureMode((HGCROCAlgosOperationMode)ps.getParameter<uint32_t>("opMode"));

  //adc configuration
  adcNbits_ = 10;
  adcMax_ = pow(2, adcNbits_);
  adcFSC_ = ps.getParameter<double>("adcFSC");
  auto adcPulseCfg = ps.getParameter<std::vector<double> >("adcPulse");
  for (size_t i = 0; i < adcPulseCfg.size(); ++i) {
    adcPulse_[i] = (float)adcPulseCfg[i];
  }
  configureADC(adcFSC_, adcPulse_);

  //toa configuration
  toaNbits_ = 10;
  toaMax_ = pow(2, toaNbits_);
  toaFSC_ = ps.getParameter<double>("toaFSC");
  toaOnset_ = ps.getParameter<double>("toaOnset");
  auto toaJitterCfg = ps.getParameter<std::vector<double> >("toaJitterParam");
  for (size_t i = 0; i < toaJitterCfg.size(); ++i) {
    toaJitterParam_[i] = toaJitterCfg[i];
  }
  configureTOA(toaFSC_, toaOnset_, toaJitterParam_);

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
  configureTOT(totFSC_, totOnset_, totChargeDrainParam_, totBxUndershoot_);
}

//
template <class DFr>
void HGCROCAlgos<DFr>::configureMode(HGCROCAlgosOperationMode mode) {
  opMode_ = mode;
}

//
template <class DFr>
void HGCROCAlgos<DFr>::configureADC(float fsc, hgc_digi::FEADCPulseShape shape) {
  adcFSC_ = fsc;
  adcLSB_ = (adcFSC_ / adcMax_);
  adcPulse_ = shape;
}

//
template <class DFr>
void HGCROCAlgos<DFr>::configureTOA(float fsc, float onset, hgc_digi::FETOAJitterParam toaJitterParam) {
  toaFSC_ = fsc;
  toaLSB_ = (toaFSC_ / toaMax_);
  toaOnset_ = onset;
  toaJitterParam_ = toaJitterParam;
}

//
template <class DFr>
void HGCROCAlgos<DFr>::configureTOT(float fsc,
                                         float onset,
                                         hgc_digi::FETDCChargeDrainParam chargeDrainParam,
                                         short bxUndershoot) {
  totFSC_ = fsc;
  totLSB_ = (totFSC_ / totMax_);
  totOnset_ = onset;
  totBxUndershoot_ = bxUndershoot;
  totChargeDrainParam_ = chargeDrainParam;
}

//
template <class DFr>
float HGCROCAlgos<DFr>::getChargeIntegrationTime(float charge, CLHEP::HepRandomEngine* engine) {
  //compute the average
  short poffsetMean = 0;
  const float timeMean = totChargeDrainParam_[poffsetMean] * charge + totChargeDrainParam_[poffsetMean + 1];
  float integTime(timeMean);

  //smear time if an engine has been passed
  if (engine) {
    short poffsetSigma = 2;
    const float timeSigma = totChargeDrainParam_[poffsetSigma] * charge + totChargeDrainParam_[poffsetSigma + 1];
    integTime = CLHEP::RandGaussQ::shoot(engine, timeMean, timeSigma);
  }

  return integTime;
};

//
template <class DFr>
float HGCROCAlgos<DFr>::getToTOnsetLeakage(float deltaT2nextBx) {
  float tau(totChargeDrainParam_[4]);
  float chargeLeakage(totOnset_ * vdt::fast_expf(-deltaT2nextBx / tau));
  return chargeLeakage;
}

//
template <class DFr>
void HGCROCAlgos<DFr>::digitizeTrivial(DFr& dataFrame,
                                            HGCSimHitData& chargeColl,
                                            hgc::HGCSimHitData& toaColl,
                                            short itbx) {
  if (verbose_) {
    edm::LogVerbatim("HGCROCAlgos::digitizeTrivial") << "[digitizeTrivial]" << std::endl;
  }

  //check if tot is to be triggered based on the charge
  //digitize charge in ADC and TOT modes and time of arrival
  float charge(chargeColl[itbx]);
  bool tp(false), tc(false);
  if (charge > totOnset_)
    tp = true;

  uint32_t adcm1(0);
  uint32_t adc = std::floor(std::min(charge, adcFSC_ - adcLSB_) / adcLSB_);
  uint32_t tot = tc ? std::floor(std::min(charge, totFSC_ - totLSB_) / totLSB_) : 0;
  uint32_t toa = charge > toaOnset_ ? std::floor(std::min(toaColl[itbx], toaFSC_ - toaLSB_) / toaLSB_) : 0;

  //fill the dataframe
  dataFrame.fill(true, tc, tp, adcm1, adc, tot, toa);

  if (verbose_) {
    std::ostringstream msg;
    dataFrame.print(msg);
    edm::LogVerbatim("HGCROCAlgos::digitizeTrivial") << msg.str() << std::endl;
  }
}

//
template <class DFr>
void HGCROCAlgos<DFr>::digitize(
    DFr& dataFrame, HGCSimHitData& chargeColl, HGCSimHitData& toaColl, CLHEP::HepRandomEngine* engine, short itbx) {
  //first look at time as the simulation is simplified for the moment
  //in all cases look only at in-time signals: ToA will be computed for the central BX only
  //to be done properly with realistic ToA shaper and jitter: for the moment accounted in the smearing
  toaFlags_.fill(false);
  float timeToA = 0.f;
  if (toaColl[itbx] != 0.f) {
    timeToA = toaColl[itbx];
    float jitter = getTimeJitter(chargeColl[itbx]);
    if (jitter != 0)
      timeToA = CLHEP::RandGaussQ::shoot(engine, timeToA, jitter);
    if (timeToA >= 0.f && timeToA <= toaFSC_)
      toaFlags_[itbx] = true;
  }

  //second charge (requires identifying bunches which will trigger ToT and then run charge sharing)
  //reset the caches as the algorithm requires some iterations over the initial values
  busyFlags_.fill(false);
  undershootFlags_.fill(false);
  totFlags_.fill(false);
  newCharge_.fill(0.f);
  toaFromToT_.fill(0.f);
  for (size_t it = 0; it < chargeColl.size(); ++it) {
    //if already flagged as busy it can't be re-used to trigger the ToT
    if (busyFlags_[it])
      continue;

    float charge = chargeColl[it];

    //if below TDC onset will be handled by SARS ADC later
    if (charge < totOnset_)
      continue;

    //raise TDC mode for charge computation
    totFlags_[it] = true;
    if (verbose_)
      edm::LogVerbatim("HGCROCAlgos::digitize")
          << "\t q=" << charge << " fC with triggers ToT @ " << it << std::endl;

    //compute total charge to be integrated and integration time
    //needs a loop as ToT will last as long as there is charge to dissipate
    uint16_t busyBxs(0);
    float totalCharge(charge), integTime(0);
    while (true) {
      //update integration time in ns and # bunches
      integTime = getChargeIntegrationTime(charge);
      const int newBusyBxs = std::floor(integTime / 25.f) + 1;

      //if no update is needed regarding the number of bunches, then the ToT integration time has converged
      if (newBusyBxs == busyBxs)
        break;

      //update charge integrated during ToT
      if (verbose_) {
        if (busyBxs == 0)
          edm::LogVerbatim("HGCROCAlgos::digitize")
              << "\t Intial busy estimate=" << integTime << " ns = " << newBusyBxs << " bxs" << std::endl;
        else
          edm::LogVerbatim("HGCROCAlgos::digitize")
              << "\t ...integrated charge overflows initial busy estimate, interating again" << std::endl;
      }

      //update number of busy bunches
      busyBxs = newBusyBxs;

      //recompute charge to be integrated adding leakage from previous bunches in SARS ADC mode
      totalCharge = charge;
      for (size_t jt = 0; jt < it; ++jt) {
        const size_t deltaT = (it - jt);
        if ((deltaT + 2) >= adcPulse_.size() || chargeColl[jt] == 0.f || totFlags_[jt] || busyFlags_[jt])
          continue;

        const float leakCharge = chargeColl[jt] * adcPulse_[deltaT + 2];
        totalCharge += leakCharge;

        if (verbose_)
          edm::LogVerbatim("HGCROCAlgos::digitize")
              << "\t\t leaking " << chargeColl[jt] << " fC @ deltaT=-" << deltaT << " -> +" << leakCharge << std::endl;
      }

      //add contamination from posterior bunches
      for (size_t jt = it + 1; jt < it + busyBxs && jt < chargeColl.size(); ++jt) {
        //this charge will be integrated in TDC mode disable for SARS ADC
        busyFlags_[jt] = true;

        const float extraCharge = chargeColl[jt];
        totalCharge += extraCharge;

        if (verbose_)
          edm::LogVerbatim("HGCROCAlgos::digitize")
              << "\t\t adding " << extraCharge << " fC @ deltaT=+" << (jt - it) << std::endl;
      }
    }

    newCharge_[it] = totalCharge;
    if (verbose_)
      edm::LogVerbatim("HGCROCAlgos::digitize")
          << "\t Final busy estimate=" << integTime << " ns = " << busyBxs << " bxs" << std::endl
          << "\t Total integrated=" << totalCharge << " fC" << std::endl;

    //last fC (tdcOnset) are dissipated through pulse with a decay constant
    //note in principle the leakage should be multiplied by (1-e^(-25/tau)) but this is very close to 1 so it's neglected
    if (it + busyBxs < newCharge_.size()) {
      const float deltaT2nextBx(busyBxs * 25 - integTime);
      const float totOnsetLeakage = getToTOnsetLeakage(deltaT2nextBx);
      newCharge_[it + busyBxs] += totOnsetLeakage;
      if (verbose_)
        edm::LogVerbatim("HGCROCAlgos::digitize")
            << "\t Leaking remainder of TOT onset " << totOnset_ << " fC, to be dissipated in " << deltaT2nextBx
            << " adds " << totOnsetLeakage << " fC @ " << it + busyBxs << " bx (first free bx)" << std::endl;
    }

    //flag bunches which will be in undershoot mode after TOT
    for (size_t jt = it + busyBxs; jt < it + busyBxs + totBxUndershoot_ && jt < chargeColl.size(); ++jt) {
      //this charge will be integrated in TDC mode disable for SARS ADC
      undershootFlags_[jt] = true;

      if (verbose_)
        edm::LogVerbatim("HGCROCAlgos::digitize")
            << "\t\t flagging undershoot for deltaT=+" << (jt - it) << std::endl;
    }
  }

  //including the leakage from bunches in SARS ADC when not declared busy or in ToT
  auto runChargeSharing = [&]() {
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
  runChargeSharing();

  //fill the dataframe

  //Tc/Tp flags
  bool tc(totFlags_[itbx] && !busyFlags_[itbx]);
  bool tp(totFlags_[itbx] || busyFlags_[itbx] || undershootFlags_[itbx]);

  //ADC BX-1 is only non zero if it was not busy or in undershoot mode
  //TOT mode is not considered for previous bunch, so ADC will be saturated if TOT was triggered
  uint16_t adcm1(0);
  int itbxm1(itbx - 1);
  if (itbxm1 >= 0 && !busyFlags_[itbxm1] && !undershootFlags_[itbxm1]) {
    adcm1 = std::floor(std::min(newCharge_[itbxm1], adcFSC_ - adcLSB_) / adcLSB_);
  }

  //ADC, TOT, TOA in-time BX
  //note: in principle the TOT has an intrinsic offset but this is ignored for the moment
  uint16_t adc(0), toa(0), tot(0);
  if (!busyFlags_[itbx] && !undershootFlags_[itbx]) {
    adc = std::floor(std::min(newCharge_[itbx], adcFSC_ - adcLSB_) / adcLSB_);
    if (totFlags_[itbx] || opMode_ == HGCROCAlgosOperationMode::CHARACTERIZATION)
      tot = std::floor(std::min(newCharge_[itbx], totFSC_ - totLSB_) / totLSB_);
    toa = newCharge_[itbx] > toaOnset_ ? std::floor(std::min(timeToA, toaFSC_ - toaLSB_) / toaLSB_) : 0;
  }

  dataFrame.fill(opMode_ == HGCROCAlgosOperationMode::CHARACTERIZATION, tc, tp, adcm1, adc, tot, toa);

  if (verbose_) {
    std::ostringstream msg;
    dataFrame.print(msg);
    edm::LogVerbatim("HGCROCAlgos") << msg.str() << std::endl;
  }
}

// trigger the compiler to generate the appropriate code
#include "DataFormats/HGCalDigi/interface/HGCalDigiCollections.h"
template class HGCROCAlgos<HGCROCChannelDataFrameSpec>;
