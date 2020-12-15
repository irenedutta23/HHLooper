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

  TaggerFitsAna(const TString &inputFileList="foo.txt", const char *outFileName="histo.root",const char *dataset="data",const char *isData="F");
  virtual ~TaggerFitsAna();
  Bool_t   FillChain(TChain *chain, const TString &inputFileList);
  Long64_t LoadTree(Long64_t entry);
  void     EventLoop(const char *, const char *, double);
  void     Categorization(const char *, const char *, float , float,double, TString procname); 
 //void     Categorization(const char *, const char *, float , float, double);
  void     GenInfo(const char *, const char *, double);
  void     BookHistogram(const char *);
 
  TFile *oFile;
  //define histograms
  TH1D *catyield;
};

#endif

#ifdef TaggerFitsAna_cxx

void TaggerFitsAna::BookHistogram(const char *outFileName) {

//  char hname[200], htit[200];
//  float xlow = 0.0,  xhigh = 2000.0;
//  int nbins = 2000;0

  oFile = new TFile(outFileName, "recreate");
  //oFile->mkdir("Cutflow");
  //oFile->cd("Cutflow");
}

TaggerFitsAna::TaggerFitsAna(const TString &inputFileList, const char *outFileName, const char* dataset, const char *isData)
{
TChain *tree = new TChain("tree");

  if( ! FillChain(tree, inputFileList) ) {
    std::cerr << "Cannot get the tree " << std::endl;
  } else {
    std::cout << "Initiating analysis of dataset " << dataset << std::endl;
    char temp[]="T";

    if(strcmp(temp,isData)==0)std::cout<<"Initiating analysis on Data"<<endl;
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
#endif // #ifdef TaggerFitsAna_cxx
