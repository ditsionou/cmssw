#include "DataFormats/HGCalDigi/interface/HGCalDigiCollections.h"
#include "SimCalorimetry/HGCalSimProducers/interface/HGCFEElectronics.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSetReader/interface/ParameterSetReader.h"
#include "FWCore/Utilities/interface/FileInPath.h"

#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/RandFlat.h"

#include <iostream>
#include <cassert>
#include <string>
#include <chrono>
#include <random>

#define NTRIALS 10000000

//helper method to print arrays
template <class T, std::size_t N>
std::ostream &operator<<(std::ostream &o, const std::array<T, N> &arr) {
  std::copy(arr.cbegin(), arr.cend(), std::ostream_iterator<T>(o, ","));
  return o;
}

//print a custom debug message before assert bombs
void massert(bool cond,
             HGCROCChannelDataFrameSpec &dfr,
             hgc_digi::HGCSimHitData &chargeColl,
             hgc_digi::HGCSimHitData &toa,
             std::string extra = "") {
  if (!cond) {
    std::cout << "Tests fail with" << std::endl
              << "Simulated | charge: " << chargeColl << std::endl
              << "          | toa: " << toa << std::endl;
    dfr.print(std::cout);
    std::cout << extra << std::endl;
  }
  assert(cond);
}

//
int main(int argc, char **argv) {
  //if passed at command line use new cfi to configure the electronics emulator
  std::string url("SimCalorimetry/HGCalSimProducers/python/hgcROCParameters_cfi.py");
  if (argc > 1)
    url = argv[1];
  url = edm::FileInPath(url).fullPath();

  //get configuration and instantiate the ROC emulator
  const std::shared_ptr<edm::ParameterSet> &pset = edm::readPSetsFrom(url);
  const edm::ParameterSet &cfg = pset->getParameter<edm::ParameterSet>("hgcROCSettings");
  float adcFSC = cfg.getParameter<double>("adcFSC");
  float toaOnset = cfg.getParameter<double>("toaOnset");
  float toaFSC = cfg.getParameter<double>("toaFSC");
  float totOnset = cfg.getParameter<double>("totOnset");
  float totFSC = cfg.getParameter<double>("totFSC");

  std::cout << "Testing HGCFEElectronics class with the following configuration:" << std::endl;
  std::cout << cfg << std::endl << std::endl;

  //simulated data arrays
  hgc_digi::HGCSimHitData chargeColl, toa;
  size_t itbx(9);

  //output dataframe
  HGCROCChannelDataFrameSpec dfr;

  //random number generator
  CLHEP::HepJamesRandom *hre = new CLHEP::HepJamesRandom();
  hre->setSeed(0);

  //digitizer
  HGCFEElectronics<HGCROCChannelDataFrameSpec> roc(cfg);
  roc.configureNoise(0.f);
  float adcLSB(roc.adcLSB()), totLSB(roc.totLSB()), toaLSB(roc.toaLSB());
  hgc_digi::FEADCPulseShape pulseShape = roc.adcPulse();

  //TEST: in-time charge deposits only
  //do once in default mode and another in characterization mode
  chargeColl.fill(0.f);
  toa.fill(0.f);
  for (size_t m = 0; m < 2; m++) {
    roc.configureMode(m == 0 ? roc.HGCFEElectronicsOperationMode::DEFAULT
                             : roc.HGCFEElectronicsOperationMode::CHARACTERIZATION);

    for (int i = 0; i < NTRIALS; i++) {
      //prepare charge/time injection
      float qin = CLHEP::RandFlat::shoot(hre, 0, totFSC);
      chargeColl[itbx] = qin;
      toa[itbx] = CLHEP::RandFlat::shoot(hre, 0, toaFSC);

      //run and check digitization output
      roc.run(dfr, chargeColl, toa, hre, itbx);
      uint16_t digitoa = dfr.toa();
      uint16_t adc = dfr.adc(roc.opMode() == roc.HGCFEElectronicsOperationMode::CHARACTERIZATION);
      uint16_t adcm1 = dfr.adcm1(roc.opMode() == roc.HGCFEElectronicsOperationMode::CHARACTERIZATION);
      uint16_t tot = dfr.tot(roc.opMode() == roc.HGCFEElectronicsOperationMode::CHARACTERIZATION);

      //no ADC for BX-1 in all cases
      massert(adcm1 == 0, dfr, chargeColl, toa);

      //TOA: check resolution when triggered, otherwise require 0
      //with noise 0 ToA has perfect resolution
      if (qin > toaOnset) {
        int delta = floor(fabs(toa[itbx] / toaLSB - (digitoa + 0.5)));
        massert(delta < 1, dfr, chargeColl, toa, "Reconstructed ToA has poor resolution");
      } else {
        massert(digitoa == 0, dfr, chargeColl, toa, "Expect ToA=0");
      }

      //ADC mode: no Tc/Tp, TOT=0, check ADC resolution
      if (qin < totOnset) {
        massert(!dfr.tc(), dfr, chargeColl, toa, "Expect Tc=0");
        massert(!dfr.tp(), dfr, chargeColl, toa, "Expect Tp=0");
        if (roc.opMode() == roc.HGCFEElectronicsOperationMode::DEFAULT) {
          massert(dfr.tot() == 0, dfr, chargeColl, toa, "Expect ToT=0");
        }

        int delta = floor(fabs(qin / adcLSB - (adc + 0.5)));
        massert(delta < 1, dfr, chargeColl, toa, "Reconstructed charge (ADC) has poor resolution");
      }
      //TOT mode: Tc=Tp=True, ADC=0, check TOT resolution
      //TOT resolution is spoiled by the 12b->10b compression
      else {
        massert(dfr.tc(), dfr, chargeColl, toa, "Expect Tc=1");
        massert(dfr.tp(), dfr, chargeColl, toa, "Expect Tp=1");
        if (roc.opMode() == roc.HGCFEElectronicsOperationMode::DEFAULT) {
          massert(adc == 0, dfr, chargeColl, toa, "Expect ADC=0");
        }
        int delta = floor(fabs(qin / totLSB - (tot + 0.5)));
        massert(delta < 5, dfr, chargeColl, toa, "Reconstructed charge (TOT) has poor resolution");
      }
    }

    std::cout << "[In-time only charge deposits] OK for opMode=" << roc.opMode() << std::endl;
  }

  //test ADC <-> TOT transition in default mode
  //scan the neighborhood of the totOnset threshold
  //firs in steps of adcLSB, second in steps of tdcLSB
  roc.configureMode(roc.HGCFEElectronicsOperationMode::DEFAULT);
  toa.fill(0.f);
  chargeColl.fill(0.f);
  for (size_t k = 0; k < 2; k++) {
    float deltaqin(k == 0 ? adcLSB : totLSB);

    for (int i = -2; i <= 2; i++) {
      float qin = totOnset + i * deltaqin;
      chargeColl[itbx] = qin;
      roc.run(dfr, chargeColl, toa, hre, itbx);
      uint16_t adc = dfr.adc(roc.opMode() == roc.HGCFEElectronicsOperationMode::CHARACTERIZATION);
      uint16_t tot = dfr.tot(roc.opMode() == roc.HGCFEElectronicsOperationMode::CHARACTERIZATION);

      if (i < 0) {
        massert(tot == 0, dfr, chargeColl, toa, "Unexpected ToT measurement");
        int delta = floor(fabs(qin / adcLSB - (adc + 0.5)));
        massert(delta < 1, dfr, chargeColl, toa, "Reconstructed charge (ADC) has poor resolution");
      } else {
        massert(adc == 0, dfr, chargeColl, toa, "Unexpected ADC measurement");
        int delta = floor(fabs(qin / totLSB - (tot + 0.5)));
        massert(delta < 5, dfr, chargeColl, toa, "Reconstructed charge (TOT) has poor resolution");
      }
    }
    std::cout << "[ADC <-> TOT transition] in steps of " << deltaqin << " OK" << std::endl;
  }

  //test charge injection @ BX-1 (within ADC range)
  //inject random charge in BX-1 and check the leakage
  for (size_t m = 0; m < 2; m++) {
    roc.configureMode(m == 0 ? roc.HGCFEElectronicsOperationMode::DEFAULT
                             : roc.HGCFEElectronicsOperationMode::CHARACTERIZATION);

    toa.fill(0.f);
    chargeColl.fill(0.f);

    for (int i = 0; i <= NTRIALS; i++) {
      float qin = CLHEP::RandFlat::shoot(hre, 0, adcFSC);
      chargeColl[itbx - 1] = qin;
      roc.run(dfr, chargeColl, toa, hre, itbx);

      uint16_t adc = dfr.adc(roc.opMode() == roc.HGCFEElectronicsOperationMode::CHARACTERIZATION);
      uint16_t adcm1 = dfr.adcm1(roc.opMode() == roc.HGCFEElectronicsOperationMode::CHARACTERIZATION);
      uint16_t tot = dfr.tot(roc.opMode() == roc.HGCFEElectronicsOperationMode::CHARACTERIZATION);

      //both operation modes check resolution of leaked pulse in BX
      float qleak = qin * pulseShape[3];
      int delta = floor(fabs(qleak / adcLSB - (adc + 0.5)));
      massert(delta < 1, dfr, chargeColl, toa, "Leaked charge has poor resolution");

      //mode==CHAR : check ADC BX-1 = 0 and ADC-TOT resolution
      if (roc.opMode() == roc.HGCFEElectronicsOperationMode::CHARACTERIZATION) {
        massert(adcm1 == 0, dfr, chargeColl, toa, "Unexpected ADC BX-1 reading");

        float delta = fabs((adc + 0.5) * adcLSB - (tot + 0.5) * totLSB);
        massert(delta < totLSB, dfr, chargeColl, toa, "ADC - TOT difference has poor resolution");
      }
      //mode==DEFAULT : check resolution of ADC BX-1
      else {
        int delta = floor(fabs(qin / adcLSB - (adcm1 + 0.5)));
        massert(delta < 1, dfr, chargeColl, toa, "ADC-1 has poor resolution");
        massert(tot == 0, dfr, chargeColl, toa, "Unexpected TOT measurement");
      }
    }
    std::cout << "[BX-1 charge injection] OK with mode=" << roc.opMode() << std::endl;
  }

  //all tests done
  return 0;
}
