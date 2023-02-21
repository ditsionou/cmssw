#include <iostream>
#include <algorithm>

#include "Rivet/Tools/ParticleIdUtils.hh"

#include "FWCore/Utilities/interface/FileInPath.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSetReader/interface/ParameterSetReader.h"

#include "fastjet/ClusterSequence.hh"

#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/RandPoisson.h"

#include "TRandom.h"
#include "TMath.h"
#include "TString.h"
#include "TChain.h"
#include "TFile.h"
#include "Math/Vector4D.h"
#include "TH2F.h"

using namespace std;
using namespace fastjet;

const UInt_t nToAThresholds=10;
const Float_t toaThresholds[nToAThresholds]={12,20,24,30,40,50,60,70,80,100};
std::map<TString, TH1F *> histos;

/**
   @short reads the required entries from the chain and returns a vector of particles for clustering
 */
std::vector<PseudoJet> getParticlesFrom(TChain *t,int ievt,int nevts,bool ispu, TString tag="GenPart") {

  //atach variables to read from
  UInt_t maxPseudoParticles(10000);
  UInt_t nPseudoParticle;
  Int_t PseudoParticle_pdgId[maxPseudoParticles], PseudoParticle_status[maxPseudoParticles];
  Float_t PseudoParticle_eta[maxPseudoParticles], PseudoParticle_mass[maxPseudoParticles], PseudoParticle_phi[maxPseudoParticles], PseudoParticle_pt[maxPseudoParticles];
  Bool_t PseudoParticle_crossedBoundary[maxPseudoParticles];
  t->SetBranchAddress("n"+tag,&nPseudoParticle);
  t->SetBranchAddress(tag+"_pdgId",PseudoParticle_pdgId);
  bool requireStatus(tag.Contains("GenPart"));
  if(requireStatus)
    t->SetBranchAddress(tag+"_status",PseudoParticle_status);
  t->SetBranchAddress(tag+"_pt",PseudoParticle_pt);
  t->SetBranchAddress(tag+"_eta",PseudoParticle_eta);
  t->SetBranchAddress(tag+"_phi",PseudoParticle_phi);
  t->SetBranchAddress(tag+"_mass",PseudoParticle_mass);
  bool requireCrossing(tag.Contains("SimClus"));
  if(requireCrossing)
    t->SetBranchAddress(tag+"_crossedBoundary",PseudoParticle_crossedBoundary);

  std::vector<PseudoJet> pseudoParticles;

  //loop over the required events
  UInt_t ientries(t->GetEntries());
  for(int i=ievt; i<ievt+nevts; i++){

    int ientry(i);
    if(ispu) ientry=gRandom->Integer(ientries-1);
    t->GetEntry(ientry);
    
    //convert kinematics to pseudojets
    for(UInt_t n=0; n<nPseudoParticle; n++) {

      if(requireStatus && PseudoParticle_status[n]!=1) continue;
      if(Rivet::PID::isNeutrino(PseudoParticle_pdgId[n])) continue;
      if(requireCrossing && PseudoParticle_crossedBoundary[n]==false) continue;

      ROOT::Math::PtEtaPhiMVector p4(PseudoParticle_pt[n],PseudoParticle_eta[n],PseudoParticle_phi[n],PseudoParticle_mass[n]);
      auto ip = PseudoJet(p4.px(),p4.py(),p4.pz(),p4.energy());

      //check nature of this particle
      bool isHadron( Rivet::PID::isHadron(PseudoParticle_pdgId[n]) );
      bool isCharged( Rivet::PID::isCharged(PseudoParticle_pdgId[n]) );
      bool isMuon( Rivet::PID::isMuon(PseudoParticle_pdgId[n]) );
      int pfid(0);
      if(isMuon) pfid=13;
      else {
        if(isHadron) {
          if(isCharged) pfid=211;
          else pfid=130;
        } else {
          pfid=22;
        }
      }
      
      //identify as pileup or signal and type of particle (e/g, charged hadron, neutral hadron or muon)
      ip.set_user_index(ispu ? -pfid : pfid);
      pseudoParticles.push_back( ip );
    }
  }
  

  return pseudoParticles;
}

/**
   @short returns a weight which is ~1 if shower can be tagged by HGCAL and ~0 if not
 */
float timeTaggedByHGCAL(float en, float thr=12,bool isem=false){

  if(thr<0) return 1.;

  float enthr(0.),eff(1.0);
  if(isem) {
    eff=1.0;
    enthr=TMath::CosH(2)*(0.000060309567449866875*pow(thr,2)+0.022640337461370418*thr+0.2636313805211345);
    enthr=TMath::Max((Double_t)(enthr),(Double_t)1.0*TMath::CosH(2));
  }else {
    eff=1.0;
    enthr=TMath::CosH(2)*(0.00035331364772303635*pow(thr,2)+0.1273691858065309*thr-0.47506821606965366);
    enthr=TMath::Max((Double_t)(enthr),(Double_t)(1.0*TMath::CosH(2)));
  }

  float x=(en-enthr)/0.5;
  return 0.5*eff*(1+TMath::Erf(x));
}


/**
   @short returns the expected time resolution by HGCAL for a shower of a given energy and using a given ToA threshold
 */
float timeResolutionByHGCAL(float en, float thr=12,bool isem=false){

  if(thr<0) return 0.;
  if(en<0.5) return 12500.;

  float A(0.), B(0.), C(20.);
  if(isem) {
    B=15.01;
  }else {
    B=-0.014548962655619186*pow(thr,2)+3.500515522985467*thr+167.6628226353696;
  }

  float x(en/TMath::CosH(2));  //the original sample was at eta=2
  float resolsq=pow(A,2)/x + pow(B/x,2) + pow(C,2);
  return resolsq;
}


float timeResolutionWeight(float resol, float thr=90.,float cte=20.){
  float r=(thr-resol)/20;
  return 0.5*(1+TMath::Erf(r));
}




/**
   @fill jet constituents based on the energy fraction
*/
struct JetConstituents_t {
  JetConstituents_t() : 
    nhf(0), chf(0), emf(0), muf(0), 
    punhf(0), puchf(0), puemf(0), pumuf(0)
  {
  };
  float nhf, chf, emf, muf;
  float punhf, puchf, puemf, pumuf;
};
JetConstituents_t fillJetConstituents(const PseudoJet &j,TString tag="",bool fillHistos=false) {

  float en(j.e());
  JetConstituents_t jc;
  for(auto c : j.constituents()) {

    float cen(c.e());
    float cenf(cen/en);
    int pfid=c.user_index();

    //compute the time tagging weights

    if(pfid==22)   {
      jc.emf += cenf;
      if(fillHistos) histos[tag+"emf_en"]->Fill(cen);
    }
    else if(pfid==-22)  {
      jc.puemf += cenf;
      if(fillHistos) histos[tag+"puemf_en"]->Fill(cen);
    }
    else if(pfid==130)  {
      jc.nhf+= cenf;
      if(fillHistos) histos[tag+"nhf_en"]->Fill(cen);
    }
    else if(pfid==-130) {
      jc.punhf += cenf;
      if(fillHistos) histos[tag+"punhf_en"]->Fill(cen);
    }
    else if(pfid==13)   jc.muf+= cenf;
    else if(pfid==-13)  jc.pumuf += cenf;
    else if(pfid==211)  {
      jc.chf+= cenf;
      if(fillHistos) histos[tag+"chf_en"]->Fill(cen);
    }
    else if(pfid==-211) {
      jc.puchf += cenf;
      if(fillHistos) histos[tag+"puchf_en"]->Fill(cen);
    }
    else 
      std::cout << "PFid = " << pfid << " ?? " << std::endl;
  }
  
  return jc;
};


int main(int argc, char** argv) {

  //if passed at command line use new cfi
  std::string url("Configuration/GenProduction/test/mixandcluster_cfi.py");
  if (argc > 1) url = argv[1];
  url = edm::FileInPath(url).fullPath();

  //get configuration
  const std::shared_ptr<edm::ParameterSet> &pset = edm::readPSetsFrom(url);
  const edm::ParameterSet &cfg = pset->getParameter<edm::ParameterSet>("mixandcluster");
  std::vector<std::string> pufiles = cfg.getParameter<std::vector<std::string> >("pu");
  std::vector<std::string> sigfiles = cfg.getParameter<std::vector<std::string> >("sig");
  int avgpu = cfg.getParameter<int>("avgpu");
  int maxevts = cfg.getParameter<int>("maxevts");

  //fast jet definition
  int jetAlgo = cfg.getParameter<int>("jetAlgo");
  double jetR = cfg.getParameter<double>("jetR");

  //
  TFile *fout=TFile::Open(Form("jets_ak%d_vbfhgg.root",int(jetR*10)),"RECREATE");
  TTree *tree = new TTree("data","data");
  UInt_t nPU,puMode;
  Float_t toaThr;
  Float_t GenJet_pt,GenJet_en, GenJet_eta;
  Float_t Jet_pt,Jet_eta,Jet_en;
  Float_t PuJet_en,PuJet_eta;
  Float_t PuJet_chf, PuJet_puchf, PuJet_nhf, PuJet_punhf, PuJet_emf, PuJet_puemf;
         
  tree->Branch("nPU",&nPU);
  tree->Branch("GenJet_pt",&GenJet_pt);
  tree->Branch("GenJet_en",&GenJet_en);
  tree->Branch("GenJet_eta",&GenJet_eta);
  tree->Branch("Jet_pt",&Jet_pt);
  tree->Branch("Jet_eta",&Jet_eta);
  tree->Branch("Jet_en",&Jet_en);
  tree->Branch("PuJet_en",&PuJet_en);
  tree->Branch("PuJet_eta",&PuJet_eta);
  tree->Branch("toaThr",&toaThr);
  tree->Branch("puMode",&puMode);
  tree->Branch("PuJet_chf",&PuJet_chf);
  tree->Branch("PuJet_puchf",&PuJet_puchf);
  tree->Branch("PuJet_nhf",&PuJet_nhf);
  tree->Branch("PuJet_punhf",&PuJet_punhf);
  tree->Branch("PuJet_emf",&PuJet_emf);
  tree->Branch("PuJet_puemf",&PuJet_puemf);
  
  TString parts[]={"chf","puchf","nhf","punhf","emf","puemf"};
  for(size_t i=0; i<sizeof(parts)/sizeof(TString); i++) {
    for(UInt_t ithr=0; ithr<nToAThresholds+1; ithr++) {
      TString tag(Form("pu%d_",ithr));
      histos[tag+parts[i]+"_en"] = new TH1F(tag+parts[i]+"_en",";Energy [GeV]; Particles",100,0,25);
    }
  }

  //fastjet definition
  JetDefinition jet_def( (JetAlgorithm)(jetAlgo), jetR);

  //start chains
  TChain *pu=new TChain("Events");
  for(size_t i=0; i<pufiles.size(); i++) pu->AddFile(pufiles[i].c_str());
  int npuEvts(pu->GetEntries());

  TChain *sig=new TChain("Events");
  for(size_t i=0; i<sigfiles.size(); i++) sig->AddFile(sigfiles[i].c_str());
  int nsigEvts(sig->GetEntries());

  std::cout << "Generating an average pileup of " << avgpu 
            << " from " << pufiles.size() << " pileup files with " << npuEvts << " events,"
            << " and signal injected from " << sigfiles.size() << " files with " << nsigEvts << " events." << std::endl;

  //random number generator
  CLHEP::HepJamesRandom *hre = new CLHEP::HepJamesRandom();
  hre->setSeed(0);
  
  maxevts = maxevts<0 || maxevts>nsigEvts ? nsigEvts : maxevts;
  for(int i=0; i<maxevts; i++) {

    nPU = UInt_t( CLHEP::RandPoisson::shoot(hre,avgpu) );

    //pure signal from gen particles
    auto sigParticles = getParticlesFrom(sig,i,1,false);
    ClusterSequence sigCS(sigParticles, jet_def);
    const auto sigGenJets = sorted_by_pt(sigCS.inclusive_jets());
    
    //pure signal from sim clusters crossing HCAL boundary
    auto sigSimClusters = getParticlesFrom(sig,i,1,false,"SimCluster");
    ClusterSequence sigsimCS(sigSimClusters, jet_def);
    const auto sigJets = sorted_by_pt(sigsimCS.inclusive_jets());
    
    //pileup
    auto puParticles = getParticlesFrom(pu,0,nPU,true,"SimCluster");

    //pileup+signal under different time-tagging hypothesis
    //puMode=0 all particles
    //puMode=1 only neutrals
    //puMode=2 only neutral hadrons
    auto sigpuParticles(puParticles);
    sigpuParticles.insert(sigpuParticles.end(), sigParticles.begin(), sigParticles.end());
    for(puMode=0; puMode<3; puMode++) {
      
      for(UInt_t ithr=0; ithr<nToAThresholds+1; ithr++) {

        toaThr=(ithr==0 ? -1 : toaThresholds[ithr-1]);
        
        std::vector<PseudoJet> sel_sigpuParticles;
        for(size_t ipart=0; ipart<sigpuParticles.size(); ipart++) {
          
          PseudoJet ijet(sigpuParticles[ipart]);
          int pfid=ijet.user_index();
          float resol=timeResolutionByHGCAL(ijet.e(),toaThr,abs(pfid)==22);
          float wgt=timeTaggedByHGCAL(ijet.e(),toaThr,abs(pfid)==22)*timeResolutionWeight(resol,90.,20.);

          if(pfid<0) wgt=1.0-wgt;
          if(puMode==0) ijet *= wgt;
          if(puMode==1 && (pfid==22 || pfid==130)) ijet*=wgt;
          if(puMode==2 && pfid==130) ijet*=wgt;
          
          sel_sigpuParticles.push_back(ijet);
        }
        
        //run clustering
        ClusterSequence sigpuCS(sel_sigpuParticles, jet_def);
        const auto sigpuJets = sorted_by_pt(sigpuCS.inclusive_jets());
        
        //select up to two signal jets and pair with the closest pileup jet
        int nSelJets(0);
        for(size_t ij=0; ij<sigGenJets.size(); ij++) {
          
          //select
          if(sigGenJets[ij].pt()<30) continue;
          if(fabs(sigGenJets[ij].eta())>3) continue;
          if(fabs(sigGenJets[ij].eta())<1.5) continue;
          if(sigGenJets[ij].constituents().size()<3) continue;
          nSelJets+=1;
          if(nSelJets>1) break;
          
          //find sim cluster jet closest in eta-phi space
          auto ijsc = min_element(begin(sigJets), end(sigJets), [=] (PseudoJet x, PseudoJet y)
          {
            return sigGenJets[ij].plain_distance(x) < sigGenJets[ij].plain_distance(y);
          });
          auto scidx = std::distance(begin(sigJets),ijsc);

          //find pileup closest in eta-phi space
          auto ijpu = min_element(begin(sigpuJets), end(sigpuJets), [=] (PseudoJet x, PseudoJet y)
          {
            return sigJets[scidx].plain_distance(x) < sigJets[scidx].plain_distance(y);
          });
          auto puidx = std::distance(begin(sigpuJets),ijpu);
          JetConstituents_t sigpujc = fillJetConstituents(sigpuJets[puidx],Form("pu%d_",ithr),puMode==0);

          //fill tree variables
          GenJet_pt=sigGenJets[ij].pt();
          GenJet_en=sigGenJets[ij].e();
          GenJet_eta=sigGenJets[ij].eta();
          Jet_pt=sigJets[scidx].pt();
          Jet_eta=sigJets[scidx].eta();
          Jet_en=sigJets[scidx].e();          
          PuJet_eta=sigpuJets[puidx].eta();
          PuJet_en=sigpuJets[puidx].e();
          PuJet_chf=sigpujc.chf;
          PuJet_puchf=sigpujc.puchf;
          PuJet_nhf=sigpujc.nhf;
          PuJet_punhf=sigpujc.punhf;
          PuJet_emf=sigpujc.emf;
          PuJet_puemf=sigpujc.puemf;
          tree->Fill();          
        }

        //debug
        if(i%10==0 && puMode==0 && ithr==0) 
          {
            std::cout << "Event: " << i << " pu=" << nPU 
                      << " signal/pu particles: " << sigParticles.size() << "/" << puParticles.size()
                      << " signal/signal+pu jets: " << sigJets.size() << "/" << sigpuJets.size() << std::endl;
          }
        
      }
    }

  }

  //save TTree
  tree->Write();
  for(auto it : histos) it.second->Write();
  fout->Close();

  return 0;
}
