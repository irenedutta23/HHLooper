#define TaggerFitsAna_cxx
#include "TaggerFitsAna.h"
#include <TH2.h>
#include <TMath.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <map>
#include <vector>
#include <cstring>
#include<string>
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
int main(int argc, char* argv[])
{

  if(argc < 3) {
    cerr << "Please give 4 arguments " << "runList " << " " << "outputFileName" << " " << "dataset" << "data type"<< " year"<<endl;
    return -1;
  }
  /*
  const char *inputFileList = argv[1];
  const char *outFileName   = argv[2];
  const char *data          = argv[3];
  const char *isData        = argv[4];
    TString scales        = argv[5];
  TaggerFitsAna hmm(inputFileList, outFileName, data, isData);
  cout << "dataset " << data << " " << scales<<endl;
  float scalef = scales.Atof();
  */
    
    const char *inputFileList = argv[1];
    const char *outFileName   = argv[2];
    const char *data          = argv[3];
    const char *isdata        = argv[4];
    const TString year          = argv[5];
    bool isData = false;
    if (*isdata=='1') isData = true;
    
    TaggerFitsAna tfit(inputFileList, outFileName, data, isData, year);
    cout << "dataset " << data << " " << endl;
    double lumi;
    if(year=="2016"){lumi  = 35.9*1000.;}
    else{
      if(year=="2017"){lumi  = 41.529*1000.;} 
      else{
	if(year=="2018"){lumi  = 59.741*1000.;}
      }
    }
    TString procname   = argv[3];
    if(!isData) cout <<"process name: "<<procname<<endl;
    double scale =1.0;
    tfit.EventLoop(procname, isData, scale,  year, lumi);
  return 0;
}

void TaggerFitsAna::EventLoop(TString procname, bool isData, double scale, TString year, double lumi)
{  if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if(jentry%5000 ==0) cout<<jentry<<endl;
      double puWeight = 1;      
      double myWeight = 1;
      if (!isData) {	 
	myWeight = lumi * weight * triggerEffWeight * pileupWeight * scale;
      }

      //******************************
      //Trigger Selection
      //******************************
      bool passTrigger = false;
      
      if (year == "2016") {
	    passTrigger = 
	      (0 == 1)
	      || HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20
	      || HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20
	      || HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20	 	    
	      ;       

	    // apply trigger efficiency correction for some triggers that were not enabled for full run
	    if (!isData) {
	      // double triggerSF = 1.0;
	      // if (HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20)                         triggerSF = 1.0;
	      // else if (HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20 
	      // 	     || HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20)                 triggerSF = 19.9 / 35.9;	      
	      // else                                                                      triggerSF = 0;
	      // myWeight = myWeight * triggerSF;	  

	      passTrigger = true;
	      // double triggerEff = 1.0 - 
	      //   (1 - getTriggerEff( triggerEff2016Hist , fatJet1Pt, fatJet1MassSD )) * 
	      //   (1 - getTriggerEff( triggerEff2016Hist , fatJet2Pt, fatJet2MassSD ))
	      //   ;
	      // myWeight = myWeight * triggerEff;
	    }

	  }
	  if (year == "2017") {
	    passTrigger = 
	      (0 == 1) 
	      || HLT_PFJet500    
	      || HLT_AK8PFJet500 
	      || HLT_AK8PFJet360_TrimMass30
	      || HLT_AK8PFJet380_TrimMass30
	      || HLT_AK8PFJet400_TrimMass30   
	      || HLT_AK8PFHT800_TrimMass50 
	      || HLT_AK8PFJet330_PFAK8BTagCSV_p17	  
	      ;       

	    // apply trigger efficiency correction for some triggers that were not enabled for full run
	    if (!isData) {
	      // double triggerSF = 1.0;
	      // if (HLT_PFJet500 || HLT_AK8PFJet500)                                    triggerSF = 1.0;
	      // else {
	      //   //cout << "fail\n";
	      //   if (HLT_AK8PFJet400_TrimMass30 || HLT_AK8PFHT800_TrimMass50)          triggerSF = 36.42 / 41.48;
	      //   else if (HLT_AK8PFJet380_TrimMass30)                                    triggerSF = 31.15 / 41.48;            
	      //   else if (HLT_AK8PFJet360_TrimMass30)                                    triggerSF = 28.23 / 41.48;
	      //   else if (HLT_AK8PFJet330_PFAK8BTagCSV_p17)                              triggerSF = 7.73 / 41.48;	    
	      //   else                                                                    triggerSF = 0;
	      //   //cout << "triggerSF = " << triggerSF << "\n";
	      // }
	      // myWeight = myWeight * triggerSF;
	    
	      passTrigger = true;
	      // double triggerEff = 1.0 - 
	      //   (1 - getTriggerEff( triggerEff2017Hist , fatJet1Pt, fatJet1MassSD )) * 
	      //   (1 - getTriggerEff( triggerEff2017Hist , fatJet2Pt, fatJet2MassSD ))
	      //   ;
	      // myWeight = myWeight * triggerEff;
	    }


	  }


	  if (year == "2018") {
	    passTrigger = 
	      (0 == 1) 
	      || HLT_AK8PFJet400_TrimMass30 
	      || HLT_AK8PFHT800_TrimMass50     
	      || HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4 
	      || HLT_AK8PFJet420_TrimMass30
	      || HLT_PFHT1050
	      || HLT_PFJet500
	      || HLT_AK8PFJet500
	      || HLT_AK8PFJet330_PFAK8BTagCSV_p17
	      
	      ;    


	    // apply trigger efficiency correction for some triggers that were not enabled for full run
	    if (!isData) {
	      // double triggerSF = 1.0;
	      // if (HLT_AK8PFJet400_TrimMass30 || HLT_AK8PFHT800_TrimMass50)              triggerSF = 1.0;
	      // else if (HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4)              triggerSF = 54.5 / 59.7;	        
	      // else                                                                      triggerSF = 0;

	      // myWeight = myWeight * triggerSF;	  

	      passTrigger = true;
	      // double triggerEff = 1.0 - 
	      //   (1 - getTriggerEff( triggerEff2018Hist , fatJet1Pt, fatJet1MassSD )) * 
	      //   (1 - getTriggerEff( triggerEff2018Hist , fatJet2Pt, fatJet2MassSD ))
	      //   ;
	      // myWeight = myWeight * triggerEff;
	    }   
	  }


	  // if (isData) {
	  //   cout << "year = " << year << " : " << passTrigger << "\n";
	  // }

	  if (!passTrigger) continue;


	  //if (isData) continue;

	  //******************************
	  //Selection Cuts 
	  //******************************
	  if ( !(fatJet1Pt > 450 )) continue;
	  if ( !(fatJet1MassSD > 30)) continue;

	  double vjetsKF = 1.;
	  float dR = 1.;
	  if(procname.Contains("zqq")){
	    dR = DeltaR(fatJet1Eta,fatJet1Phi, genZEta, genZPhi);
	  }
	  else{
	    if( procname.Contains("wqq")){
	      dR = DeltaR(fatJet1Eta,fatJet1Phi, genWEta, genWPhi);  
	    }
	  }
	  if(dR<0.8){
	   if(procname.Contains("zqq")){
	     vjetsKF = kFactor_Lookup(genZPt, procname);
	   }
	   else{
	     if( procname.Contains("wqq")){
	       vjetsKF = kFactor_Lookup(genWPt, procname);  
	     }
	   } 
	  }
	  if(!isData){myWeight*= vjetsKF;}

	  //******************************
	  //Fill histograms
	  //******************************
      
	  if (fatJet1HasMuon || fatJet1HasElectron) continue;
	  double fatJet1Rho = 2*log(fatJet1MassSD/fatJet1Pt);
	  histMET->Fill(MET, myWeight);    
	  histJet1Mass->Fill(fatJet1MassSD, myWeight);      
	  histJet1Pt->Fill(fatJet1Pt, myWeight);
	  histJet1Eta->Fill(fatJet1Eta, myWeight);
	  histJet1Phi->Fill(fatJet1Phi, myWeight);
	  histJet1DDB->Fill(fatJet1DDBTagger, myWeight);
	  histJet1PNetXbb->Fill(fatJet1PNetXbb, myWeight);
	  histJet1Tau3OverTau2->Fill(fatJet1Tau3OverTau2, myWeight);	    
	  histJet1n2b1->Fill(fatJet1_n2b1,myWeight);
	  histJet1HasMuon->Fill(fatJet1HasMuon, myWeight);
	  histJet1HasElectron->Fill(fatJet1HasElectron, myWeight);
	  histJet1HasBJetCSVLoose->Fill(fatJet1HasBJetCSVLoose, myWeight);
	  histJet1HasBJetCSVMedium->Fill(fatJet1HasBJetCSVMedium, myWeight);
	  histJet1HasBJetCSVTight->Fill(fatJet1HasBJetCSVTight, myWeight);
	  histJet1OppositeHemisphereHasBJet->Fill(fatJet1OppositeHemisphereHasBJet, myWeight);
	  histJet1Rho->Fill(fatJet1Rho, myWeight);

	  if(fatJet1PNetXbb > 0.95){
	    histMET_Pass->Fill(MET, myWeight);
	    histJet1Mass_Pass->Fill(fatJet1MassSD, myWeight);
	    histJet1Pt_Pass->Fill(fatJet1Pt, myWeight);
	    histJet1Eta_Pass->Fill(fatJet1Eta, myWeight);
	    histJet1Phi_Pass->Fill(fatJet1Phi, myWeight);
	    histJet1DDB_Pass->Fill(fatJet1DDBTagger, myWeight);
	    histJet1PNetXbb_Pass->Fill(fatJet1PNetXbb, myWeight);
	    histJet1Tau3OverTau2_Pass->Fill(fatJet1Tau3OverTau2, myWeight);
	    histJet1n2b1_Pass->Fill(fatJet1_n2b1,myWeight);
	    histJet1HasMuon_Pass->Fill(fatJet1HasMuon, myWeight);
	    histJet1HasElectron_Pass->Fill(fatJet1HasElectron, myWeight);
	    histJet1HasBJetCSVLoose_Pass->Fill(fatJet1HasBJetCSVLoose, myWeight);
	    histJet1HasBJetCSVMedium_Pass->Fill(fatJet1HasBJetCSVMedium, myWeight);
	    histJet1HasBJetCSVTight_Pass->Fill(fatJet1HasBJetCSVTight, myWeight);
	    histJet1OppositeHemisphereHasBJet_Pass->Fill(fatJet1OppositeHemisphereHasBJet, myWeight);
	    histJet1Rho_Pass->Fill(fatJet1Rho, myWeight);

	  }
	  else{
	    histMET_Fail->Fill(MET, myWeight);
	    histJet1Mass_Fail->Fill(fatJet1MassSD, myWeight);
	    histJet1Pt_Fail->Fill(fatJet1Pt, myWeight);
	    histJet1Eta_Fail->Fill(fatJet1Eta, myWeight);
	    histJet1Phi_Fail->Fill(fatJet1Phi, myWeight);
	    histJet1DDB_Fail->Fill(fatJet1DDBTagger, myWeight);
	    histJet1PNetXbb_Fail->Fill(fatJet1PNetXbb, myWeight);
	    histJet1Tau3OverTau2_Fail->Fill(fatJet1Tau3OverTau2, myWeight);
	    histJet1n2b1_Fail->Fill(fatJet1_n2b1,myWeight);
	    histJet1HasMuon_Fail->Fill(fatJet1HasMuon, myWeight);
	    histJet1HasElectron_Fail->Fill(fatJet1HasElectron, myWeight);
	    histJet1HasBJetCSVLoose_Fail->Fill(fatJet1HasBJetCSVLoose, myWeight);
	    histJet1HasBJetCSVMedium_Fail->Fill(fatJet1HasBJetCSVMedium, myWeight);
	    histJet1HasBJetCSVTight_Fail->Fill(fatJet1HasBJetCSVTight, myWeight);
	    histJet1OppositeHemisphereHasBJet_Fail->Fill(fatJet1OppositeHemisphereHasBJet, myWeight);
	    histJet1Rho_Fail->Fill(fatJet1Rho, myWeight);
	  }
  
   }
}
