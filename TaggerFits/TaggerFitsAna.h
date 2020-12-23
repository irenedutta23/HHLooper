/////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Sep 10 04:03:58 2018 by ROOT version 6.10/09
// from TTree tree/tree
// found on file: ../condor/condor_output/condor_logs/ElData16_job0.root
//////////////////////////////////////////////////////////

#ifndef TaggerFitsAna_h
#define TaggerFitsAna_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <cmath>
// Header file for the classes stored in the TTree if any.
#include "vector"
#include "NtupleVariables.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TLorentzVector.h"
#include<string>
#include "TF1.h"
#include "TAxis.h"
class TaggerFitsAna:public NtupleVariables  {
public :

  TaggerFitsAna(const TString &inputFileList="foo.txt", const char *outFileName="histo.root",const char *dataset="data",bool isData = 0, const TString &year="2018");
  virtual ~TaggerFitsAna();
  Bool_t   FillChain(TChain *chain, const TString &inputFileList);
  Long64_t LoadTree(Long64_t entry);
  void     EventLoop(TString, bool, double, const TString, double);
  void     BookHistogram(const char *);
  double   kFactor_Lookup(float genpt, TString procname);

  TFile *oFile;
  TH1D *hQCD_Z;
  TH1D *hQCD_W;
  TH1D *hEWK_Z;
  TH1D *hEWK_W;
  TH1D *h_znlo;
  TH1D *h_wnlo;

  //define histograms
  TH1D *histMET;
  TH1D *histJet1Mass;
  TH1D *histJet1Pt;
  TH1D *histJet1Eta;
  TH1D *histJet1Phi;
  TH1D *histJet1DDB;
  TH1D *histJet1PNetXbb;
  TH1D *histJet1Tau3OverTau2;
  TH1D *histJet1n2b1;
  TH1D *histJet1HasMuon;
  TH1D *histJet1HasElectron;
  TH1D *histJet1HasBJetCSVLoose;
  TH1D *histJet1HasBJetCSVMedium;
  TH1D *histJet1HasBJetCSVTight;
  TH1D *histJet1OppositeHemisphereHasBJet;
  TH1D *histJet1Rho;

  TH1D *histMET_Pass;
  TH1D *histJet1Mass_Pass;
  TH1D *histJet1Pt_Pass;
  TH1D *histJet1Eta_Pass;
  TH1D *histJet1Phi_Pass;
  TH1D *histJet1DDB_Pass;
  TH1D *histJet1PNetXbb_Pass;
  TH1D *histJet1Tau3OverTau2_Pass;
  TH1D *histJet1n2b1_Pass;
  TH1D *histJet1HasMuon_Pass;
  TH1D *histJet1HasElectron_Pass;
  TH1D *histJet1HasBJetCSVLoose_Pass;
  TH1D *histJet1HasBJetCSVMedium_Pass;
  TH1D *histJet1HasBJetCSVTight_Pass;
  TH1D *histJet1OppositeHemisphereHasBJet_Pass;
  TH1D *histJet1Rho_Pass;

  TH1D *histMET_Fail;
  TH1D *histJet1Mass_Fail;
  TH1D *histJet1Pt_Fail;
  TH1D *histJet1Eta_Fail;
  TH1D *histJet1Phi_Fail;
  TH1D *histJet1DDB_Fail;
  TH1D *histJet1PNetXbb_Fail;
  TH1D *histJet1Tau3OverTau2_Fail;
  TH1D *histJet1n2b1_Fail;
  TH1D *histJet1HasMuon_Fail;
  TH1D *histJet1HasElectron_Fail;
  TH1D *histJet1HasBJetCSVLoose_Fail;
  TH1D *histJet1HasBJetCSVMedium_Fail;
  TH1D *histJet1HasBJetCSVTight_Fail;
  TH1D *histJet1OppositeHemisphereHasBJet_Fail;
  TH1D *histJet1Rho_Fail;
};

#endif

#ifdef TaggerFitsAna_cxx

void TaggerFitsAna::BookHistogram(const char *outFileName) {

//  char hname[200], htit[200];
//  float xlow = 0.0,  xhigh = 2000.0;
//  int nbins = 2000;0

  oFile = new TFile(outFileName, "recreate");
  oFile->mkdir("Base");
  oFile->cd("Base");
  histMET = new TH1D("h_MET", "; MET [GeV] ; Number of Events", 100, 0, 200);
  histJet1Mass = new TH1D("h_Jet1Mass", "; Leading Jet Mass [GeV] ; Number of Events", 75, 50, 200);
  histJet1Pt = new TH1D("h_Jet1Pt", "; Leading Jet p_{T} [GeV] ; Number of Events", 50, 0, 2000);
  histJet1Eta = new TH1D("h_Jet1Eta", "; Leading Jet #eta ; Number of Events", 20, -5, 5);
  histJet1Phi = new TH1D("h_Jet1Phi", "; Leading Jet #phi ; Number of Events", 20, -3.2, 3.2);
  histJet1DDB = new TH1D("h_Jet1DDB", "; Leading Jet DDB ; Number of Events", 25, 0, 1.0);
  histJet1PNetXbb = new TH1D("h_Jet1PNetXbb", "; Leading Jet PNetXbb ; Number of Events", 25, 0, 1.0);
  histJet1Tau3OverTau2 = new TH1D("h_Jet1Tau3OverTau2", "; Leading Jet #tau_{32} ; Number of Events", 25, 0, 1.0);
  histJet1n2b1 = new TH1D("h_Jet1DDB", "; Leading Jet N2b1 ; Number of Events", 25, 0, 0.5);
  histJet1HasMuon = new TH1D("h_Jet1HasMuon", "; Leading Jet Has Muon; Number of Events", 2, 0, 2.0);
  histJet1HasElectron = new TH1D("h_Jet1HasElectron", "; Leading Jet Has Electron; Number of Events", 2, 0, 2.0);
  histJet1HasBJetCSVLoose = new TH1D("h_Jet1HasBJetCSVLoose", "; Leading Jet is loose CSV; Number of Events", 2, 0, 2.0);
  histJet1HasBJetCSVMedium = new TH1D("h_Jet1HasBJetCSVMedium", "; Leading Jet is Medium CSV; Number of Events", 2, 0, 2.0);
  histJet1HasBJetCSVTight = new TH1D("h_Jet1HasBJetCSVTight", "; Leading Jet is Tight CSV ; Number of Events", 2, 0, 2.0);
  histJet1OppositeHemisphereHasBJet = new TH1D("h_Jet1OppositeHemisphereHasBJet", "; Leading Jet has a b-jet in #Delta#phi>2.5; Number of Events", 2, 0, 2.0);
  histJet1Rho = new TH1D("h_Jet1Rho", "; Leading Jet #rho ; Number of Events", 60, -7, -1);

  oFile->mkdir("Pass");
  oFile->cd("Pass");
  histMET_Pass = new TH1D("h_MET", "; MET [GeV] ; Number of Events", 100, 0, 200);
  histJet1Mass_Pass = new TH1D("h_Jet1Mass", "; Leading Jet Mass [GeV] ; Number of Events", 75, 50, 200);
  histJet1Pt_Pass = new TH1D("h_Jet1Pt", "; Leading Jet p_{T} [GeV] ; Number of Events", 50, 0, 2000);
  histJet1Eta_Pass = new TH1D("h_Jet1Eta", "; Leading Jet #eta ; Number of Events", 20, -5, 5);
  histJet1Phi_Pass = new TH1D("h_Jet1Phi", "; Leading Jet #phi ; Number of Events", 20, -3.2, 3.2);
  histJet1DDB_Pass = new TH1D("h_Jet1DDB", "; Leading Jet DDB ; Number of Events", 25, 0, 1.0);
  histJet1PNetXbb_Pass = new TH1D("h_Jet1PNetXbb", "; Leading Jet PNetXbb ; Number of Events", 25, 0, 1.0);
  histJet1Tau3OverTau2_Pass = new TH1D("h_Jet1Tau3OverTau2", "; Leading Jet #tau_{32} ; Number of Events", 25, 0, 1.0);
  histJet1n2b1_Pass = new TH1D("h_Jet1DDB", "; Leading Jet N2b1 ; Number of Events", 25, 0, 0.5);
  histJet1HasMuon_Pass = new TH1D("h_Jet1HasMuon", "; Leading Jet Has Muon; Number of Events", 2, 0, 2.0);
  histJet1HasElectron_Pass = new TH1D("h_Jet1HasElectron", "; Leading Jet Has Electron; Number of Events", 2, 0, 2.0);
  histJet1HasBJetCSVLoose_Pass = new TH1D("h_Jet1HasBJetCSVLoose", "; Leading Jet is loose CSV; Number of Events", 2, 0, 2.0);
  histJet1HasBJetCSVMedium_Pass = new TH1D("h_Jet1HasBJetCSVMedium", "; Leading Jet is Medium CSV; Number of Events", 2, 0, 2.0);
  histJet1HasBJetCSVTight_Pass = new TH1D("h_Jet1HasBJetCSVTight", "; Leading Jet is Tight CSV ; Number of Events", 2, 0, 2.0);
  histJet1OppositeHemisphereHasBJet_Pass = new TH1D("h_Jet1OppositeHemisphereHasBJet", "; Leading Jet has a b-jet in #Delta#phi>2.5; Number of Events", 2, 0, 2.0);
  histJet1Rho_Pass = new TH1D("h_Jet1Rho", "; Leading Jet #rho ; Number of Events", 60, -7, -1);

  oFile->mkdir("Fail");
  oFile->cd("Fail");
  histMET_Fail = new TH1D("h_MET", "; MET [GeV] ; Number of Events", 100, 0, 200);
  histJet1Mass_Fail = new TH1D("h_Jet1Mass", "; Leading Jet Mass [GeV] ; Number of Events", 75, 50, 200);
  histJet1Pt_Fail = new TH1D("h_Jet1Pt", "; Leading Jet p_{T} [GeV] ; Number of Events", 50, 0, 2000);
  histJet1Eta_Fail = new TH1D("h_Jet1Eta", "; Leading Jet #eta ; Number of Events", 20, -5, 5);
  histJet1Phi_Fail = new TH1D("h_Jet1Phi", "; Leading Jet #phi ; Number of Events", 20, -3.2, 3.2);
  histJet1DDB_Fail = new TH1D("h_Jet1DDB", "; Leading Jet DDB ; Number of Events", 25, 0, 1.0);
  histJet1PNetXbb_Fail = new TH1D("h_Jet1PNetXbb", "; Leading Jet PNetXbb ; Number of Events", 25, 0, 1.0);
  histJet1Tau3OverTau2_Fail = new TH1D("h_Jet1Tau3OverTau2", "; Leading Jet #tau_{32} ; Number of Events", 25, 0, 1.0);
  histJet1n2b1_Fail = new TH1D("h_Jet1DDB", "; Leading Jet N2b1 ; Number of Events", 25, 0, 0.5);
  histJet1HasMuon_Fail = new TH1D("h_Jet1HasMuon", "; Leading Jet Has Muon; Number of Events", 2, 0, 2.0);
  histJet1HasElectron_Fail = new TH1D("h_Jet1HasElectron", "; Leading Jet Has Electron; Number of Events", 2, 0, 2.0);
  histJet1HasBJetCSVLoose_Fail = new TH1D("h_Jet1HasBJetCSVLoose", "; Leading Jet is loose CSV; Number of Events", 2, 0, 2.0);
  histJet1HasBJetCSVMedium_Fail = new TH1D("h_Jet1HasBJetCSVMedium", "; Leading Jet is Medium CSV; Number of Events", 2, 0, 2.0);
  histJet1HasBJetCSVTight_Fail = new TH1D("h_Jet1HasBJetCSVTight", "; Leading Jet is Tight CSV ; Number of Events", 2, 0, 2.0);
  histJet1OppositeHemisphereHasBJet_Fail = new TH1D("h_Jet1OppositeHemisphereHasBJet", "; Leading Jet has a b-jet in #Delta#phi>2.5; Number of Events", 2, 0, 2.0);
  histJet1Rho_Fail = new TH1D("h_Jet1Rho", "; Leading Jet #rho ; Number of Events", 60, -7, -1);

}

TaggerFitsAna::TaggerFitsAna(const TString &inputFileList, const char *outFileName, const char* dataset, bool isData, const TString &year)
{

  const char* kfactorsfile_path = "kfactors.root";
  const char* zjetsnlo_path = "ZJets_QCD_NLO.root";
  const char* wjetsnlo_path = "WJets_QCD_NLO.root";
  TFile* kfactorsFile = TFile::Open(kfactorsfile_path);
  TFile* zjetsNLOFile = TFile::Open(zjetsnlo_path);
  TFile* wjetsNLOFile = TFile::Open(wjetsnlo_path);

  hQCD_Z = (TH1D*)kfactorsFile->Get("ZJets_012j_NLO/nominal");
  hQCD_W = (TH1D*)kfactorsFile->Get("WJets_012j_NLO/nominal");
  hEWK_Z = (TH1D*)kfactorsFile->Get("EWKcorr/Z");
  hEWK_W = (TH1D*)kfactorsFile->Get("EWKcorr/W");
  hEWK_Z->Divide(hQCD_Z);
  hEWK_W->Divide(hQCD_W);


  if(year=="2017" || year=="2018"){
    h_znlo = (TH1D*)zjetsNLOFile->Get("Z_NLO_QCD_2017");
    h_wnlo = (TH1D*)wjetsNLOFile->Get("W_NLO_QCD_2017");
  }
  else{
    if(year=="2016"){
      h_znlo = (TH1D*)zjetsNLOFile->Get("Z_NLO_QCD_2016");
      h_wnlo = (TH1D*)wjetsNLOFile->Get("W_NLO_QCD_2016");
    }
  }

  TChain *tree = new TChain("tree");

  if( ! FillChain(tree, inputFileList) ) {
    std::cerr << "Cannot get the tree " << std::endl;
  } else {
    std::cout << "Initiating analysis of dataset " << dataset << std::endl;
    
    if(isData)std::cout<<"Initiating analysis on Data"<<endl;
    else std::cout<<"Initiating analysis on MC"<<endl;
  }
   
  NtupleVariables::Init(tree);
  BookHistogram(outFileName);
  }
Bool_t TaggerFitsAna::FillChain(TChain *chain, const TString &inputFileList) {

  ifstream infile(inputFileList, ifstream::in);
  std::string buffer;

  if(!infile.is_open()) {
    std::cerr << "** ERROR: Can't open '" << inputFileList << "' for input" << std::endl;
    return kFALSE;
  }

  std::cout << "TreeUtilities : FillChain " << std::endl;
  while(1) {
    infile >> buffer;
    if(!infile.good()) break;
    std::cout << "Adding tree from " << buffer.c_str() << std::endl;                                                              
    chain->Add(buffer.c_str());
  }
  std::cout << "No. of Entries in this tree : " << chain->GetEntries() << std::endl;
  
  return kTRUE;
}

TaggerFitsAna::~TaggerFitsAna()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
   oFile->cd();
   oFile->Write();
   oFile->Close();
   
}

Long64_t TaggerFitsAna::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
    if (!fChain->InheritsFrom(TChain::Class()))  return centry;
  TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      // Notify();
   }
   return centry;
}

double TaggerFitsAna::kFactor_Lookup(float genpt, TString procname){
  double vjetsKF = 1.;
  double iEWKKF = 1.;
  double iQCDKF = 1.;

  if(procname.Contains("zqq")){
    iEWKKF = hEWK_Z->GetBinContent(hEWK_Z->FindBin(genpt));
    iQCDKF = h_znlo->GetBinContent(h_znlo->FindBin(genpt));
  }
  else{
    if( procname.Contains("wqq")){
      iEWKKF = hEWK_W->GetBinContent(hEWK_W->FindBin(genpt));
      iQCDKF = h_wnlo->GetBinContent(h_wnlo->FindBin(genpt));
    }
  }
  vjetsKF = iQCDKF*iEWKKF;
  return vjetsKF;
}
#endif // #ifdef TaggerFitsAna_cxx
