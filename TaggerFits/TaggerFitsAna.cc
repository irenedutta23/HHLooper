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
    cerr << "Please give 4 arguments " << "runList " << " " << "outputFileName" << " " << "dataset" << "data type"<<endl;
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
    const char *isData        = argv[4];
    TaggerFitsAna tfit(inputFileList, outFileName, data, isData);
    cout << "dataset " << data << " " << endl;
    
    map<TString,double> proc_scale;
    double lumi = 41.529*1000.;
    double lumi_18 = 59.74*1000.;
    double lumi_16 = 35.9*1000.;
   
    TString procname   = argv[3];
    if(*isData=='F') cout <<"process name: "<<procname<<" scale: "<<proc_scale[procname]<<endl;
   
  if(*isData=='T') tfit.EventLoop(data, isData, 1.0);
  else tfit.EventLoop(data, isData, scalef);
  return 0;
}

void TaggerFitsAna::EventLoop(const char *data,const char *isData, double scale)
{  if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if(jentry%5000 ==0) cout<<jentry<<endl;
   }
}
