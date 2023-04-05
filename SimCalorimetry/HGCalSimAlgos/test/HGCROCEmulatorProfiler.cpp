#include "DataFormats/HGCalDigi/interface/HGCalDigiCollections.h"
#include "SimCalorimetry/HGCalSimAlgos/interface/HGCROCEmulator.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSetReader/interface/ParameterSetReader.h"
#include "FWCore/Utilities/interface/FileInPath.h"

#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/RandFlat.h"

#include "TFile.h"
#include "TTree.h"

#include <iostream>
#include <cassert>
#include <string>
#include <chrono>
#include <random>

//
int main(int argc, char **argv) {

  //if passed at command line use new cfi to configure the electronics emulator
  std::string url("SimCalorimetry/HGCalSimAlgos/python/hgcrocEmulator_cfi.py");
  if (argc > 1) url = argv[1];
  url = edm::FileInPath(url).fullPath();
  
  //get configuration and instantiate the ROC emulator
  const std::shared_ptr<edm::ParameterSet> &pset = edm::readPSetsFrom(url);
  const edm::ParameterSet &cfg = pset->getParameter<edm::ParameterSet>("hgcrocEmulator");
  std::cout << "Profiling HGCROCEmulator behavior with the following configuration:" << std::endl;
  std::cout << cfg << std::endl << std::endl;

  //digitizer
  HGCROCEmulator<HGCROCChannelDataFrameSpec> roc(cfg);
  float adcFSC(roc.adcFSC()), totFSC(roc.totFSC()), toaFSC(roc.toaFSC());
  //float toaOnset(roc.toaOnset()), totOnset(roc.totOnset());
  float toaLSB(roc.toaLSB());
  HGCROCOperationMode opMode(roc.opMode());
    
  //simulated data arrays
  HGCROCSimHitData_t chargeColl, chargeToA;
  size_t itbx(9);

  //output dataframe
  HGCROCChannelDataFrameSpec dfr;

  //random number generator
  CLHEP::HepJamesRandom *hre = new CLHEP::HepJamesRandom();
  hre->setSeed(0);

  //book some profiles and start output root file
  TFile *fOut=TFile::Open("hgrocemulatorprofile.root","RECREATE");
  TTree *t=new TTree("ch","HGCROC channel data");
  t->SetDirectory(fOut);
  float qinj,toainj;
  t->Branch("qinj",&qinj, "qinj/F");
  t->Branch("toainj",&toainj, "toainj/F");
  uint16_t toa,adc,adcm1,tot,tctp,busy;
  t->Branch("adc",&adc, "adc/s");
  t->Branch("adcm1",&adcm1, "adcm1/s");
  t->Branch("tot",&tot, "tot/s");
  t->Branch("toa",&toa, "toa/s");
  t->Branch("tctp",&tctp, "tctp/s");
  t->Branch("busy",&busy, "busy/s");

  //profile responses
  for(unsigned long i=0; i<10000; i++) {
    
    //random values (full range every 5)
    float fsc(i%5==0 ? totFSC*1.1 : adcFSC*1.1);
    chargeColl[itbx] = CLHEP::RandFlat::shoot(hre, 0, fsc);
    chargeToA[itbx] = CLHEP::RandFlat::shoot(hre, 0, toaFSC*1.1);
    
    //run and check digitization output
    roc.run(dfr, chargeColl, chargeToA, hre, true, itbx);
    toa = dfr.toa();
    adc = dfr.adc(opMode);
    adcm1 = dfr.adcm1(opMode);
    tot = dfr.tot(opMode);
    tctp = (dfr.tc() << 1) | dfr.tp();
    busy = roc.currentIntegrationTimes()[itbx];

    //truth
    qinj=chargeColl[itbx];
    toainj=chargeToA[itbx]/toaLSB;

    if(i<10)
      std::cout << i << " | q=" << qinj << " tcptp=" << tctp
                << " tot=" << tot << " busy=" << busy
                << " adc=" << adc << " adcm1=" << adcm1
                << " | t=" << toainj << " toa=" << toa << std::endl;
        

    t->Fill();    
  }
  
  //close root file and write profiles to it
  t->Write();
  fOut->Close();
  
  //all done
  return 0;
}
