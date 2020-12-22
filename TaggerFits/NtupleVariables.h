//////////////////////////////////////////////////////////
// 
// Original Author : Irene Dutta
//                   Caltech
// Date Created    : Mon 27 Aug, 2018
//////////////////////////////////////////////////////////
#ifndef NtupleVariables_h
#define NtupleVariables_h

#include <TROOT.h>
#include <TChain.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TSelector.h>
#include "TH1F.h"
#include "TH2F.h"

// Header file for the classes stored in the TTree if any.
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.
#ifdef __CINT__

#pragma link C++ class vector<float>+;
#pragma link C++ class vector<int>+;
#pragma link C++ class vector<bool>+;
#endif

using namespace std;
class NtupleVariables : public TSelector {
public :

   NtupleVariables(TTree * /*tree*/ =0) : fChain(0) { }
   ~NtupleVariables() { }
   void    Init(TTree *tree);
   Bool_t  Notify();
   Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   double  DeltaPhi(double, double);
   double  DeltaR(double eta1, double phi1, double eta2, double phi2);
   double  HT (std::vector<TLorentzVector>);
   std::pair<double,double> helicityAngles( TLorentzVector v1, TLorentzVector v2 );
   std::pair<double,double> CSAngles( TLorentzVector v1, TLorentzVector v2, int charge);
   TLorentzVector MHT(std::vector<TLorentzVector>);
   
   // void LeptHist(TH1D *h[20],const char *leptype, int Zcut);
   //ClassDef(NtupleVariables,0);

   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   
   // Declaration of leaf types
      Float_t         weight;
   Float_t         triggerEffWeight;
   Float_t         pileupWeight;
   Float_t         totalWeight;
   UInt_t          run;
   UInt_t          lumi;
   ULong64_t       event;
   Float_t         npu;
   Float_t         rho;
   Int_t           NJets;
   Float_t         MET;
   Float_t         fatJet1Pt;
   Float_t         fatJet1Eta;
   Float_t         fatJet1Phi;
   Float_t         fatJet1Mass;
   Float_t         fatJet1MassSD;
   Float_t         fatJet1DDBTagger;
   Float_t         fatJet1PNetXbb;
   Float_t         fatJet1PNetQCDb;
   Float_t         fatJet1PNetQCDbb;
   Float_t         fatJet1PNetQCDc;
   Float_t         fatJet1PNetQCDcc;
   Float_t         fatJet1PNetQCDothers;
   Int_t           fatJet1GenMatchIndex;
   Float_t         fatJet1Tau3OverTau2;
   Float_t         fatJet1_n2b1;
   Bool_t          fatJet1HasMuon;
   Bool_t          fatJet1HasElectron;
   Bool_t          fatJet1HasBJetCSVLoose;
   Bool_t          fatJet1HasBJetCSVMedium;
   Bool_t          fatJet1HasBJetCSVTight;
   Bool_t          fatJet1OppositeHemisphereHasBJet;
   Float_t         genWPt;
   Float_t         genWEta;
   Float_t         genWPhi;
   Float_t         genZPt;
   Float_t         genZEta;
   Float_t         genZPhi;
   Float_t         pho1Pt;
   Float_t         pho1Eta;
   Float_t         pho1Phi;
   Bool_t          HLT_Ele27_WPTight_Gsf;
   Bool_t          HLT_Ele28_WPTight_Gsf;
   Bool_t          HLT_Ele30_WPTight_Gsf;
   Bool_t          HLT_Ele32_WPTight_Gsf;
   Bool_t          HLT_Ele35_WPTight_Gsf;
   Bool_t          HLT_Ele38_WPTight_Gsf;
   Bool_t          HLT_Ele40_WPTight_Gsf;
   Bool_t          HLT_IsoMu20;
   Bool_t          HLT_IsoMu24;
   Bool_t          HLT_IsoMu24_eta2p1;
   Bool_t          HLT_IsoMu27;
   Bool_t          HLT_IsoMu30;
   Bool_t          HLT_Mu50;
   Bool_t          HLT_Mu55;
   Bool_t          HLT_Photon175;
   Bool_t          HLT_PFHT780;
   Bool_t          HLT_PFHT890;
   Bool_t          HLT_PFHT1050;
   Bool_t          HLT_AK8PFJet360_TrimMass30;
   Bool_t          HLT_AK8PFJet380_TrimMass30;
   Bool_t          HLT_AK8PFJet400_TrimMass30;
   Bool_t          HLT_AK8PFJet420_TrimMass30;
   Bool_t          HLT_AK8PFHT750_TrimMass50;
   Bool_t          HLT_AK8PFHT800_TrimMass50;
   Bool_t          HLT_AK8PFHT850_TrimMass50;
   Bool_t          HLT_AK8PFHT900_TrimMass50;
   Bool_t          HLT_PFJet450;
   Bool_t          HLT_PFJet500;
   Bool_t          HLT_PFJet550;
   Bool_t          HLT_AK8PFJet450;
   Bool_t          HLT_AK8PFJet500;
   Bool_t          HLT_AK8PFJet550;
   Bool_t          HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17;
   Bool_t          HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1;
   Bool_t          HLT_AK8PFJet330_PFAK8BTagCSV_p17;
   Bool_t          HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02;
   Bool_t          HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2;
   Bool_t          HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4;
   Bool_t          HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20;
   Bool_t          HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p087;
   Bool_t          HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p087;
   Bool_t          HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20;
   Bool_t          HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20;
   Bool_t          HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20;

      TBranch        *b_weight;   //!
   TBranch        *b_triggerEffWeight;   //!
   TBranch        *b_pileupWeight;   //!
   TBranch        *b_totalWeight;   //!
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_event;   //!
   TBranch        *b_npu;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_NJets;   //!
   TBranch        *b_MET;   //!
   TBranch        *b_fatJet1Pt;   //!
   TBranch        *b_fatJet1Eta;   //!
   TBranch        *b_fatJet1Phi;   //!
   TBranch        *b_fatJet1Mass;   //!
   TBranch        *b_fatJet1MassSD;   //!
   TBranch        *b_fatJet1DDBTagger;   //!
   TBranch        *b_fatJet1PNetXbb;   //!
   TBranch        *b_fatJet1PNetQCDb;   //!
   TBranch        *b_fatJet1PNetQCDbb;   //!
   TBranch        *b_fatJet1PNetQCDc;   //!
   TBranch        *b_fatJet1PNetQCDcc;   //!
   TBranch        *b_fatJet1PNetQCDothers;   //!
   TBranch        *b_fatJet1GenMatchIndex;   //!
   TBranch        *b_fatJet1Tau3OverTau2;   //!
   TBranch        *b_fatJet1_n2b1;   //!
   TBranch        *b_fatJet1HasMuon;   //!
   TBranch        *b_fatJet1HasElectron;   //!
   TBranch        *b_fatJet1HasBJetCSVLoose;   //!
   TBranch        *b_fatJet1HasBJetCSVMedium;   //!
   TBranch        *b_fatJet1HasBJetCSVTight;   //!
   TBranch        *b_fatJet1OppositeHemisphereHasBJet;   //!
   TBranch        *b_genWPt;   //!
   TBranch        *b_genWEta;   //!
   TBranch        *b_genWPhi;   //!
   TBranch        *b_genZPt;   //!
   TBranch        *b_genZEta;   //!
   TBranch        *b_genZPhi;   //!
   TBranch        *b_pho1Pt;   //!
   TBranch        *b_pho1Eta;   //!
   TBranch        *b_pho1Phi;   //!
   TBranch        *b_HLT_Ele27_WPTight_Gsf;   //!
   TBranch        *b_HLT_Ele28_WPTight_Gsf;   //!
   TBranch        *b_HLT_Ele30_WPTight_Gsf;   //!
   TBranch        *b_HLT_Ele32_WPTight_Gsf;   //!
   TBranch        *b_HLT_Ele35_WPTight_Gsf;   //!
   TBranch        *b_HLT_Ele38_WPTight_Gsf;   //!
   TBranch        *b_HLT_Ele40_WPTight_Gsf;   //!
   TBranch        *b_HLT_IsoMu20;   //!
   TBranch        *b_HLT_IsoMu24;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1;   //!
   TBranch        *b_HLT_IsoMu27;   //!
   TBranch        *b_HLT_IsoMu30;   //!
   TBranch        *b_HLT_Mu50;   //!
   TBranch        *b_HLT_Mu55;   //!
   TBranch        *b_HLT_Photon175;   //!
   TBranch        *b_HLT_PFHT780;   //!
   TBranch        *b_HLT_PFHT890;   //!
   TBranch        *b_HLT_PFHT1050;   //!
   TBranch        *b_HLT_AK8PFJet360_TrimMass30;   //!
   TBranch        *b_HLT_AK8PFJet380_TrimMass30;   //!
   TBranch        *b_HLT_AK8PFJet400_TrimMass30;   //!
   TBranch        *b_HLT_AK8PFJet420_TrimMass30;   //!
   TBranch        *b_HLT_AK8PFHT750_TrimMass50;   //!
   TBranch        *b_HLT_AK8PFHT800_TrimMass50;   //!
   TBranch        *b_HLT_AK8PFHT850_TrimMass50;   //!
   TBranch        *b_HLT_AK8PFHT900_TrimMass50;   //!
   TBranch        *b_HLT_PFJet450;   //!
   TBranch        *b_HLT_PFJet500;   //!
   TBranch        *b_HLT_PFJet550;   //!
   TBranch        *b_HLT_AK8PFJet450;   //!
   TBranch        *b_HLT_AK8PFJet500;   //!
   TBranch        *b_HLT_AK8PFJet550;   //!
   TBranch        *b_HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17;   //!
   TBranch        *b_HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1;   //!
   TBranch        *b_HLT_AK8PFJet330_PFAK8BTagCSV_p17;   //!
   TBranch        *b_HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02;   //!
   TBranch        *b_HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2;   //!
   TBranch        *b_HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4;   //!
   TBranch        *b_HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20;   //!
   TBranch        *b_HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p087;   //!
   TBranch        *b_HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p087;   //!
   TBranch        *b_HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20;   //!
   TBranch        *b_HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20;   //!
   TBranch        *b_HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20;   //!

    

};
#endif

#ifdef NtupleVariables_cxx
void NtupleVariables::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

  // Set object pointer
  

 // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   // Set object pointer
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("triggerEffWeight", &triggerEffWeight, &b_triggerEffWeight);
   fChain->SetBranchAddress("pileupWeight", &pileupWeight, &b_pileupWeight);
   fChain->SetBranchAddress("totalWeight", &totalWeight, &b_totalWeight);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("npu", &npu, &b_npu);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("NJets", &NJets, &b_NJets);
   fChain->SetBranchAddress("MET", &MET, &b_MET);
   fChain->SetBranchAddress("fatJet1Pt", &fatJet1Pt, &b_fatJet1Pt);
   fChain->SetBranchAddress("fatJet1Eta", &fatJet1Eta, &b_fatJet1Eta);
   fChain->SetBranchAddress("fatJet1Phi", &fatJet1Phi, &b_fatJet1Phi);
   fChain->SetBranchAddress("fatJet1Mass", &fatJet1Mass, &b_fatJet1Mass);
   fChain->SetBranchAddress("fatJet1MassSD", &fatJet1MassSD, &b_fatJet1MassSD);
   fChain->SetBranchAddress("fatJet1DDBTagger", &fatJet1DDBTagger, &b_fatJet1DDBTagger);
   fChain->SetBranchAddress("fatJet1PNetXbb", &fatJet1PNetXbb, &b_fatJet1PNetXbb);
   fChain->SetBranchAddress("fatJet1PNetQCDb", &fatJet1PNetQCDb, &b_fatJet1PNetQCDb);
   fChain->SetBranchAddress("fatJet1PNetQCDbb", &fatJet1PNetQCDbb, &b_fatJet1PNetQCDbb);
   fChain->SetBranchAddress("fatJet1PNetQCDc", &fatJet1PNetQCDc, &b_fatJet1PNetQCDc);
   fChain->SetBranchAddress("fatJet1PNetQCDcc", &fatJet1PNetQCDcc, &b_fatJet1PNetQCDcc);
   fChain->SetBranchAddress("fatJet1PNetQCDothers", &fatJet1PNetQCDothers, &b_fatJet1PNetQCDothers);
   fChain->SetBranchAddress("fatJet1GenMatchIndex", &fatJet1GenMatchIndex, &b_fatJet1GenMatchIndex);
   fChain->SetBranchAddress("fatJet1Tau3OverTau2", &fatJet1Tau3OverTau2, &b_fatJet1Tau3OverTau2);
   fChain->SetBranchAddress("fatJet1_n2b1", &fatJet1_n2b1, &b_fatJet1_n2b1);
   fChain->SetBranchAddress("fatJet1HasMuon", &fatJet1HasMuon, &b_fatJet1HasMuon);
   fChain->SetBranchAddress("fatJet1HasElectron", &fatJet1HasElectron, &b_fatJet1HasElectron);
   fChain->SetBranchAddress("fatJet1HasBJetCSVLoose", &fatJet1HasBJetCSVLoose, &b_fatJet1HasBJetCSVLoose);
   fChain->SetBranchAddress("fatJet1HasBJetCSVMedium", &fatJet1HasBJetCSVMedium, &b_fatJet1HasBJetCSVMedium);
   fChain->SetBranchAddress("fatJet1HasBJetCSVTight", &fatJet1HasBJetCSVTight, &b_fatJet1HasBJetCSVTight);
   fChain->SetBranchAddress("fatJet1OppositeHemisphereHasBJet", &fatJet1OppositeHemisphereHasBJet, &b_fatJet1OppositeHemisphereHasBJet);
   fChain->SetBranchAddress("genWPt", &genWPt, &b_genWPt);
   fChain->SetBranchAddress("genWEta", &genWEta, &b_genWEta);
   fChain->SetBranchAddress("genWPhi", &genWPhi, &b_genWPhi);
   fChain->SetBranchAddress("genZPt", &genZPt, &b_genZPt);
   fChain->SetBranchAddress("genZEta", &genZEta, &b_genZEta);
   fChain->SetBranchAddress("genZPhi", &genZPhi, &b_genZPhi);
   fChain->SetBranchAddress("pho1Pt", &pho1Pt, &b_pho1Pt);
   fChain->SetBranchAddress("pho1Eta", &pho1Eta, &b_pho1Eta);
   fChain->SetBranchAddress("pho1Phi", &pho1Phi, &b_pho1Phi);
   fChain->SetBranchAddress("HLT_Ele27_WPTight_Gsf", &HLT_Ele27_WPTight_Gsf, &b_HLT_Ele27_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_Ele28_WPTight_Gsf", &HLT_Ele28_WPTight_Gsf, &b_HLT_Ele28_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_Ele30_WPTight_Gsf", &HLT_Ele30_WPTight_Gsf, &b_HLT_Ele30_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_Ele32_WPTight_Gsf", &HLT_Ele32_WPTight_Gsf, &b_HLT_Ele32_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_Ele35_WPTight_Gsf", &HLT_Ele35_WPTight_Gsf, &b_HLT_Ele35_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_Ele38_WPTight_Gsf", &HLT_Ele38_WPTight_Gsf, &b_HLT_Ele38_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_Ele40_WPTight_Gsf", &HLT_Ele40_WPTight_Gsf, &b_HLT_Ele40_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_IsoMu20", &HLT_IsoMu20, &b_HLT_IsoMu20);
   fChain->SetBranchAddress("HLT_IsoMu24", &HLT_IsoMu24, &b_HLT_IsoMu24);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1", &HLT_IsoMu24_eta2p1, &b_HLT_IsoMu24_eta2p1);
   fChain->SetBranchAddress("HLT_IsoMu27", &HLT_IsoMu27, &b_HLT_IsoMu27);
   fChain->SetBranchAddress("HLT_IsoMu30", &HLT_IsoMu30, &b_HLT_IsoMu30);
   fChain->SetBranchAddress("HLT_Mu50", &HLT_Mu50, &b_HLT_Mu50);
   fChain->SetBranchAddress("HLT_Mu55", &HLT_Mu55, &b_HLT_Mu55);
   fChain->SetBranchAddress("HLT_Photon175", &HLT_Photon175, &b_HLT_Photon175);
   fChain->SetBranchAddress("HLT_PFHT780", &HLT_PFHT780, &b_HLT_PFHT780);
   fChain->SetBranchAddress("HLT_PFHT890", &HLT_PFHT890, &b_HLT_PFHT890);
   fChain->SetBranchAddress("HLT_PFHT1050", &HLT_PFHT1050, &b_HLT_PFHT1050);
   fChain->SetBranchAddress("HLT_AK8PFJet360_TrimMass30", &HLT_AK8PFJet360_TrimMass30, &b_HLT_AK8PFJet360_TrimMass30);
   fChain->SetBranchAddress("HLT_AK8PFJet380_TrimMass30", &HLT_AK8PFJet380_TrimMass30, &b_HLT_AK8PFJet380_TrimMass30);
   fChain->SetBranchAddress("HLT_AK8PFJet400_TrimMass30", &HLT_AK8PFJet400_TrimMass30, &b_HLT_AK8PFJet400_TrimMass30);
   fChain->SetBranchAddress("HLT_AK8PFJet420_TrimMass30", &HLT_AK8PFJet420_TrimMass30, &b_HLT_AK8PFJet420_TrimMass30);
   fChain->SetBranchAddress("HLT_AK8PFHT750_TrimMass50", &HLT_AK8PFHT750_TrimMass50, &b_HLT_AK8PFHT750_TrimMass50);
   fChain->SetBranchAddress("HLT_AK8PFHT800_TrimMass50", &HLT_AK8PFHT800_TrimMass50, &b_HLT_AK8PFHT800_TrimMass50);
   fChain->SetBranchAddress("HLT_AK8PFHT850_TrimMass50", &HLT_AK8PFHT850_TrimMass50, &b_HLT_AK8PFHT850_TrimMass50);
   fChain->SetBranchAddress("HLT_AK8PFHT900_TrimMass50", &HLT_AK8PFHT900_TrimMass50, &b_HLT_AK8PFHT900_TrimMass50);
   fChain->SetBranchAddress("HLT_PFJet450", &HLT_PFJet450, &b_HLT_PFJet450);
   fChain->SetBranchAddress("HLT_PFJet500", &HLT_PFJet500, &b_HLT_PFJet500);
   fChain->SetBranchAddress("HLT_PFJet550", &HLT_PFJet550, &b_HLT_PFJet550);
   fChain->SetBranchAddress("HLT_AK8PFJet450", &HLT_AK8PFJet450, &b_HLT_AK8PFJet450);
   fChain->SetBranchAddress("HLT_AK8PFJet500", &HLT_AK8PFJet500, &b_HLT_AK8PFJet500);
   fChain->SetBranchAddress("HLT_AK8PFJet550", &HLT_AK8PFJet550, &b_HLT_AK8PFJet550);
   fChain->SetBranchAddress("HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17", &HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17, &b_HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17);
   fChain->SetBranchAddress("HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1", &HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1, &b_HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1);
   fChain->SetBranchAddress("HLT_AK8PFJet330_PFAK8BTagCSV_p17", &HLT_AK8PFJet330_PFAK8BTagCSV_p17, &b_HLT_AK8PFJet330_PFAK8BTagCSV_p17);
   fChain->SetBranchAddress("HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02", &HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02, &b_HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02);
   fChain->SetBranchAddress("HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2", &HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2, &b_HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2);
   fChain->SetBranchAddress("HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4", &HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4, &b_HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4);
   fChain->SetBranchAddress("HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20", &HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20, &b_HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20);
   fChain->SetBranchAddress("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p087", &HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p087, &b_HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p087);
   fChain->SetBranchAddress("HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p087", &HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p087, &b_HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p087);
   fChain->SetBranchAddress("HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20", &HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20, &b_HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20);
   fChain->SetBranchAddress("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20", &HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20, &b_HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20);
   fChain->SetBranchAddress("HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20", &HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20, &b_HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20);
   Notify();
}

Bool_t NtupleVariables::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef NtupleVariables_cxx
