#include "THStack.h"
#include "TF1.h"
#include "TString.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include "tdrstyle.C"
#include "CMS_lumi.C"
using std::cout;
using std::endl;
void decorate(TH1*h,const char* xtitle, const char* ytitle, const char* title,
	      int linecolor, int linewidth, int markercolor, int markerstyle, int tofill) ;
void histDraw(TPad *p,TH1F *h1a,TH1F *h1b, THStack *hs,TH1F *vFrame );
void histDraw(TPad *p,TH1F *h1a,TGraphAsymmErrors *h1b, THStack *hs,TH1F *vFrame );
void TaggerFits()
{

  //  gROOT->LoadMacro("tdrstyle.C");
  setTDRStyle();

  //  gROOT->LoadMacro("CMS_lumi.C");

  writeExtraText = true;       // if extra text
  extraText  = "Preliminary";  // default extra text is "Preliminary"
  //lumi_sqrtS = "13 TeV";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
  //lumi_8TeV  = "19.1 fb^{-1}"; // default is "19.7 fb^{-1}"
  //lumi_7TeV  = "4.9 fb^{-1}";  // default is "5.1 fb^{-1}"
  //lumi_13TeV  = "35.9  fb^{-1}";
  lumi_13TeV  = "59.74  fb^{-1}";
  //lumi_13TeV  = "41.529  fb^{-1}";
  int iPeriod = 4;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses lumi_sqrtS)
  int iPos=11; //11 = default: left-aligned
  //int iPos=33; // right-aligned
  //int iPos=22; // center-aligned
  //int iPos=0;//out of frame (in exceptional cases)
  //First declare the file names

  // references for T, B, L, R
  /* float T = 0.08*H_ref;
     float B = 0.12*H_ref; 
     float L = 0.12*W_ref;
     float R = 0.04*W_ref;
  */

  TString workdir = "2020-12-23/";
  TString f[77];
  f[0]=workdir+"Data_2018.root";
  f[1]=workdir+"QCD_2018.root";
  f[2]=workdir+"ttbar_2018.root";
  f[3]=workdir+"wqq_2018.root";
  f[4]=workdir+"zqq_2018.root";
  f[5]=workdir+"top_2018.root";
  
  //Declare other constants, strings that you might need here.

  //Now let us define the plot name to overlay
  TString plotname[33];
  plotname[0]= "h_Jet1Mass";
  plotname[1] = "h_Jet1Pt";
  plotname[2] = "h_Jet1Eta";
  plotname[3] = "h_Jet1Phi";
  plotname[4] ="h_Jet1DDB";
  plotname[5] = "h_Jet1PNetXbb";
  plotname[6]= "h_Jet1Tau3OverTau2";
  plotname[7] = "h_Jet1n2b1";
  plotname[8] = "h_Jet1HasMuon";
  plotname[9] = "h_Jet1HasElectron";
  plotname[10] = "h_Jet1HasBJetCSVLoose";
  plotname[11] = "h_Jet1HasBJetCSVMedium";
  plotname[12] = "h_Jet1HasBJetCSVTight";
  plotname[13] = "h_Jet1OppositeHemisphereHasBJet";
  plotname[14] = "h_Jet1Rho";
  plotname[15] = "h_MET";
 
  gStyle->SetOptStat(0);
  //Also give fancy name for the axis titles
  TString xtitle= "P_{T} (GeV)";
  TString xtitle1= "#eta";
  TString xtitle2= "#phi";
  TString xtitle3= "H_{T} (GeV)";
  TString ytitle = "Number of Events"; // Or "Events"
  
  //Now let us open the files
  TFile *file[77];
  for(int i=0;i<6;i++) file[i] = new TFile(f[i]);
  
  //change directory if required.
  TDirectory *directory_md0[77];
  TDirectory *Main_Dir[77][77];
  TString dir_name[77];
  dir_name[0]="Base";
  dir_name[1]="Pass";
  dir_name[2]="Fail";
  for(int k=0;k<3;k++){
  for(int i=0;i<6;i++) Main_Dir[i][k] = (TDirectory*) file[i]->Get(dir_name[k]);
  }

  TString name_unc;
  double sum=0.0;
  char hist_sum[20];
  TLegend *lg[32];
  char name[100];
  THStack *hs[32][32];
  TCanvas *c[32][32];
  TH1F *vFrame[32][32];
  TPad *pad1[32][32];
  TPad *pad2[32][32];
  TPad *pad3[32][32];
  TH1F *hd[32][32];
  TLine *l[32][32];
  //Now open the respective histograms from the files
  TH1F *hist[32][32][32];
  TH1F *hdata[32][32];
  TH1F *h[32][32], *h_err[32][32], *h_temp,*h_total[32];
  TH2F *h_covar[32];
  TGraphAsymmErrors *h3[32];
  TAxis *xaxis_temp;
  int nbins_temp;
  float xlow;
  float xup;
  TLatex *l_chi[32];
  double bin_err, stat_err,jec_err,jer_err;
  int nbins;
  int color[5] = {91,73,53,46,73};

    //Data
    //h3[j] = (TH1F*)file[j]->Get(plotname[0]);
    //decorate(h3[j],"",ytitle,"",kBlack,2,kBlack,8,0);
  for(int d=0;d<3;d++){
  for(int k=0;k<16;k++){
    cout<<k<<endl;
    hdata[d][k] = (TH1F*)Main_Dir[0][d]->Get(plotname[k]);
    hdata[d][k]->SetMarkerColor(kBlack);
    hdata[d][k]->SetLineColor(kBlack);
    //hdata[k]->SetHistLineStyle(0);
    hdata[d][k]->SetLineWidth(2);
    hdata[d][k]->SetMarkerStyle(8);
    hdata[d][k]->SetMarkerSize(1.3);
  }}
  
  for(int j=1;j<6;j++){
    for(int d=0;d<3;d++){
      for(int k=0;k<16;k++){
	//dy amc 01j
	hist[j-1][d][k] = (TH1F*)Main_Dir[j][d]->Get(plotname[k]);
	sprintf(name,"h%i%i%i",j,d,k);
	if(j==1){
	  hist[j-1][d][k]->Scale(0.85); //QCD
	  //cout<<j-1<<k<<endl;
	  h[d][k] = (TH1F*)hist[j-1][d][k]->Clone(name); //For MC // doing it here so that I can adjust the color of h[j]
	}
	decorate(hist[j-1][d][k],"",ytitle,"",color[j-1],2,color[j-1],20,1);
      }
    }
  }
    
  TFile *file_out = TFile::Open("HHTo4BPlots_Tagger.root", "UPDATE");
  file_out->cd();
  TString process[6] = {"Data_2018","QCD_2018","ttbar_2018","wqq_2018","zqq_2018","top_2018"};
  for(int i=0; i<6; i++) {
    for(int d=0;d<3;d++){
      for(int j=0;j<16;j++){
	TString out_name = plotname[j]+"_"+process[i]+"_"+dir_name[d];
	if(i==0) file_out->WriteTObject(hdata[d][j], out_name, "WriteDelete");
	else file_out->WriteTObject(hist[i-1][d][j], out_name, "WriteDelete");
	
      } 
    }
  }
  file_out->Close();
  delete file_out;       


  cout<<"got all hists\n";

  TLegend *lg1 = new TLegend(0.35,0.75,0.85,0.9);
  lg1->SetTextFont(132);
  lg1->SetBorderSize(0);
  lg1->SetNColumns(3);
  gStyle->SetLegendTextSize(0.04);
  lg1->AddEntry(hdata[0][0],"Data","p");
  
  lg1->AddEntry(hist[0][0][0],"QCD","f");
  lg1->AddEntry(hist[1][0][0],"Top","f");
  lg1->AddEntry(hist[2][0][0],"W+jets","f");
  lg1->AddEntry(hist[3][0][0],"Z+jets","l");
  
  for(int d=0;d<3;d++){
  for(int j=0;j<16;j++){
    sprintf(name,"hs%i%i",j,d);
    hs[d][j] = new THStack(name,name); // Stack of all backgrounds
    //hs[d][j]->Add(hist[3][d][j]);
    hs[d][j]->Add(hist[2][d][j]);
    hs[d][j]->Add(hist[4][d][j]);
    hs[d][j]->Add(hist[1][d][j]);
    hs[d][j]->Add(hist[0][d][j]);
  
    

    sprintf(name,"c%i_%i",j,d);
    c[d][j] = new TCanvas(name,name,800,800);
    sprintf(name,"pad1%i_%i",j,d);
    pad1[d][j]= new TPad(name, name, 0, 0.3, 1.0, 1.0);
    pad1[d][j]->SetBottomMargin(0.05); // Upper and lower plot are joined
    //pad1[d][j]->SetGridx();         // Vertical grid
    pad1[d][j]->Draw();             // Draw the Upper pad: pad1
    pad1[d][j]->cd();               // pad1 becomes the current pad
    pad1[d][j]->SetLogy();
    
    vFrame[d][j] = pad1[d][j]->DrawFrame(0, 0.3, 1.0, 700);
    histDraw(pad1[d][j],hist[0][d][j],hdata[d][j],hs[d][j],vFrame[d][j]);
    
    hdata[d][j]->Draw("PEsames");  //data 
    hist[3][d][j]->Draw("sames");    
    h[d][j]->Sumw2();
    cout<<j<<endl;

    
    h[d][j]->Add(hist[1][d][j]);
    h[d][j]->Add(hist[2][d][j]);
     h[d][j]->Add(hist[4][d][j]);
    //h[d][j]->Add(hist[3][d][j]); //exclude Z(bb)
    
    sprintf(name,"h_err%i",j);
    h_err[d][j] = (TH1F*)h[d][j]->Clone(name);
    int n=h_err[d][j]->GetNbinsX();
    for (Int_t i=1; i<=n; i++){
      h_err[d][j]->SetBinContent(i,1.0);
    }
    //h[d][j]->SetFillColorAlpha(40,0.35);
    //h[d][j]->SetMarkerStyle(20);
    h[d][j]->SetFillColor(12);
    //h[d][j]->SetLineWidth(20);
    h[d][j]->SetFillStyle(3018);
    //h[d][j]->Draw("e2 same");
   
    vFrame[d][j]->GetYaxis()->SetTitle(ytitle);
    vFrame[d][j]->GetYaxis()->SetTitleSize(0.06);
    vFrame[d][j]->GetYaxis()->SetTitleOffset(0.8);
    vFrame[d][j]->GetYaxis()->SetLabelSize(0.05);
    vFrame[d][j]->GetXaxis()->SetLabelOffset(999);
    vFrame[d][j]->GetXaxis()->SetLabelSize(0);
    lg1->Draw();
    CMS_lumi( pad1[d][j], iPeriod, iPos );
    gPad->RedrawAxis();

    
    c[d][j]->cd();          // Go back to the main canvas before defining pad2

    sprintf(name,"pad2%i%i",j,d);
    pad2[d][j] = new TPad(name, name , 0, 0.0, 1.0, 0.3);
    pad2[d][j]->SetTopMargin(0.1);
    pad2[d][j]->SetBottomMargin(0.25);
    pad2[d][j]->SetGridx(); // vertical grid
    pad2[d][j]->SetGridy(); 
    pad2[d][j]->Draw();
    pad2[d][j]->cd();

   

    sprintf(name,"htemp%i%i",j,d);
  
    sprintf(name,"hd%i%i",j,d);
    hd[d][j] = (TH1F*)hdata[d][j]->Clone(name); // For data *************
    
    hd[d][j]->SetLineColor(kBlack);
    hd[d][j]->Sumw2();
    hd[d][j]->SetStats(0);// No statistics on lower plot
    hd[d][j]->SetMinimum(0.);  // Define Y ..
    hd[d][j]->SetMaximum(2.); // .. range
    //hd[d][j]->Divide(h[d][j]); //Data/MC
    hd[d][j]->Add(h[d][j], -1); // Data -Mc (except Z(bb))
    n=hdata[d][j]->GetNbinsX();
    for (Int_t i=1; i<=n; i++) {
      if(h[d][j]->GetBinContent(i)!=0){
	//bin_err =(hdata[d][j]->GetErrorYhigh(i)+hdata[d][j]->GetErrorYlow(i))/(2*h[d][j]->GetBinContent(i));
	bin_err=hdata[d][j]->GetBinError(i)/h[d][j]->GetBinContent(i);
	//cout<<hdata[d][j]->GetErrorYhigh(i)<<"," <<hdata[d][j]->GetErrorYlow(i)<<endl;
	//cout<<hdata[d][j]->GetBinError(i)<<","<<h[d][j]->GetBinContent(i)<<","<<bin_err<<endl;
	//cout<<bin_err<<endl;
	//hd[d][j]->SetBinError(i,bin_err);
      }
    }
    n=hd[d][j]->GetNbinsX();
    double chi_2=0.0;
    int bin_ct=0;
    for (Int_t i=1; i<=n; i++) {
      if(hd[d][j]->GetBinContent(i)!=0){
	chi_2+=pow((hd[d][j]->GetBinContent(i)-1),2.);
	bin_ct+=1;
      }
    }
    chi_2/=bin_ct;
    
    //cout<<"Chi_2/ndof:"<<chi_2<<endl;
    
   
    hd[d][j]->SetMarkerStyle(21);
    hd[d][j]->Draw("ep");
    for (Int_t i=1; i<=n; i++) {
      if(h[d][j]->GetBinContent(i)!=0){
	bin_err =h[d][j]->GetBinError(i)/h[d][j]->GetBinContent(i);
	cout<<i<<","<<bin_err<<h[d][j]->GetBinContent(i)<<endl;
	h_err[d][j]->SetBinError(i,bin_err);
      }
    }
    h_err[d][j]->SetFillColor(40);
    //h[d][j]->SetLineWidth(20);
    h_err[d][j]->SetFillStyle(3001);
    //h_err[d][j]->Draw("e2 same");
    hd[d][j]->Draw("ep same");
    l_chi[j] = new TLatex(0.8,0.85,Form("#chi^{2}/ndof : %.3g",chi_2));
    l_chi[j]->SetNDC();
    l_chi[j]->SetTextSize(0.08);
    //l_chi[d][j]->Draw("same");
    
    hd[d][j]->SetTitle(""); // Remove the ratio title

    // Y axis ratio plot settings
    
    hd[d][j]->GetYaxis()->SetTitle("Data-Bkg");
    hd[d][j]->GetYaxis()->SetNdivisions(5);
    hd[d][j]->GetYaxis()->SetTitleSize(30);
    hd[d][j]->GetYaxis()->SetTitleFont(43);
    hd[d][j]->GetYaxis()->SetTitleOffset(1.);
    hd[d][j]->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    hd[d][j]->GetYaxis()->SetLabelSize(25);

    // X axis ratio plot settings
    hd[d][j]->GetXaxis()->SetTitleSize(30);
    hd[d][j]->GetXaxis()->SetTitleFont(43);
    hd[d][j]->GetXaxis()->SetTitleOffset(2.5);
    hd[d][j]->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    hd[d][j]->GetXaxis()->SetLabelSize(25);
    h_err[d][j]->GetYaxis()->SetTitle("Data/Bkg");
    h_err[d][j]->GetYaxis()->SetNdivisions(5);
    h_err[d][j]->GetYaxis()->SetTitleSize(30);
    h_err[d][j]->GetYaxis()->SetTitleFont(43);
    h_err[d][j]->GetYaxis()->SetTitleOffset(1.);
    h_err[d][j]->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    h_err[d][j]->GetYaxis()->SetLabelSize(25);

    // X axis ratio plot settings
    h_err[d][j]->GetXaxis()->SetTitleSize(30);
    h_err[d][j]->GetXaxis()->SetTitleFont(43);
    h_err[d][j]->GetXaxis()->SetTitleOffset(2.5);
    h_err[d][j]->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    h_err[d][j]->GetXaxis()->SetLabelSize(25);
    h_err[d][j]->GetXaxis()->SetLabelOffset(0.02);
    h_err[d][j]->SetStats(0);// No statistics on lower plot
    h_err[d][j]->SetMinimum(0.);  // Define Y ..
    h_err[d][j]->SetMaximum(2.); // .. range
   sprintf(name,"l%i",j);
   
  }

 
 
  }
  for(int d=0;d<3;d++){
  for(int j=0;j<16;j++){
    
    pad2[d][j]->cd();
    double xmax= hist[0][d][j]->GetXaxis()->GetXmax();
    double xmin = hist[0][d][j]->GetXaxis()->GetXmin();
    l[d][j]= new TLine(xmin,1.0,xmax,1.0);
    l[d][j]->SetLineColor(kBlack); l[d][j]->SetLineWidth(2); l[d][j]->SetLineStyle(1);
    l[d][j]->Draw("same");
   
   
    //c[d][j]->cd();
    //c[d][j]->Update();
    //c[d][j]->RedrawAxis();
    //c[d][j]->GetFrame()->Draw();
    /* name_unc=plotname[d][j]+".pdf";
    c[d][j]->SaveAs(name_unc);
    name_unc=plotname[d][j]+".png";
    c[d][j]->SaveAs(name_unc);
    name_unc=plotname[d][j]+".C";
    c[d][j]->SaveAs(name_unc);
    */

  }
  }
  for(int d=0;d<3;d++){
  for(int j=0;j<16;j++){
    //c[d][j]->SaveAs(name_unc);
    name_unc=workdir+plotname[j]+"_"+dir_name[d]+".png";
    c[d][j]->SaveAs(name_unc);
    //name_unc=plotname[d][j]+".pdf";
    //c[d][j]->SaveAs(name_unc);
    //name_unc=plotname[d][j]+".C";
    //c[d][j]->SaveAs(name_unc);
    delete c[d][j];
  }
  }

}


//Define the function decorate for histograms
void decorate(TH1*h,const char* xtitle, const char* ytitle, const char* title,
	      int linecolor, int linewidth, int markercolor, int markerstyle, int tofill) {

  h->GetXaxis()->SetTitle(xtitle);
  h->GetYaxis()->SetTitle(ytitle);

  h->SetLineColor(linecolor); 
  h->SetLineWidth(linewidth);
  h->SetMarkerColor(markercolor);
  h->SetMarkerStyle(markerstyle);
  if(tofill==1) h->SetFillColor(markercolor);
  
  h->SetMarkerSize(1.3);
  h->SetTitle(title);
  h->SetMinimum(0.00001);
}
//Define the function decorate for legends
void decorate(TLegend *g, float textSize, TString legendheader)
{
  g->SetTextSize(textSize);
  g->SetFillStyle(1);
  g->SetFillColor(100);
  g->SetBorderSize(1);
  g->SetLineColor(10);
  //Usually legends should not have headers
  //g->SetHeader(legendheader);
}
/*
void histDraw(TPad *p,TH1F *h1a,TH1F *h1b, THStack *hs,TH1F *vFrame )
{
  
  double maxCont = -1.0;
  maxCont=h1a->GetBinContent(h1a->GetMaximumBin());
  if(maxCont<h1b->GetBinContent(h1b->GetMaximumBin()))maxCont=h1b->GetBinContent(h1b->GetMaximumBin());
  maxCont*=100.0;
  double maxRange = h1a->GetXaxis()->GetXmax();
  double minRange=h1a->GetXaxis()->GetXmin();
  //if(h1a->GetBinCenter(h1a->GetNbinsX()) > maxRange) maxRange =  h1a->GetBinCenter(h1a->GetNbinsX());
  vFrame = p->DrawFrame(minRange, 0.01, maxRange, maxCont);
  hs->Draw("same hist");
  }*/

void histDraw(TPad *p,TH1F *h1a, TH1F *h1b, THStack *hs,TH1F *vFrame )
{
  
  double maxCont = -1.0;
  maxCont=h1a->GetBinContent(h1a->GetMaximumBin());
  if(maxCont<h1b->GetMaximum())maxCont=h1b->GetMaximum();
  maxCont*=1000.0;
  double maxRange = h1a->GetXaxis()->GetXmax();
  double minRange=h1a->GetXaxis()->GetXmin();
  //if(h1a->GetBinCenter(h1a->GetNbinsX()) > maxRange) maxRange =  h1a->GetBinCenter(h1a->GetNbinsX());
  vFrame = p->DrawFrame(minRange, 0.01, maxRange, maxCont);
  hs->Draw("same hist");
}
// Here are a couple of other utility functions

// For a given histogram hst, return the number of entries between bin_lo and bin_hi
/*
float get_nevents(TH1F *hst, float bin_lo, float bin_hi)
{
  int bin_width = hst->GetBinWidth(1);
  int ibin_begin = 1;
  float nevents = 0;
  while ( hst->GetBinCenter(ibin_begin) < bin_lo )
    ibin_begin++;
  int ibin_end = ibin_begin;
  while ( hst->GetBinCenter(ibin_end) < bin_hi )
    ibin_end++;
  for(int i=ibin_begin; i<ibin_end; i++)
    nevents += hst->GetBinContent(i);

  return nevents;
}
// Partner function for above, returning the error for the above nevents
float get_nevents_err(TH1F *hst, float bin_lo, float bin_hi)
{
  int bin_width = hst->GetBinWidth(1);
  int ibin_begin = 1;
  float nevents_err = 0;
  while ( hst->GetBinCenter(ibin_begin) < bin_lo )
    ibin_begin++;
  int ibin_end = ibin_begin;
  while ( hst->GetBinCenter(ibin_end) < bin_hi )
    ibin_end++;
  for(int i=ibin_begin; i<ibin_end; i++)
    nevents_err += pow(hst->GetBinError(i),2);
  nevents_err = sqrt(nevents_err);

  return nevents;
}
*/
void h12ascii (TH1* h)
{
   Int_t n = h->GetNbinsX();
   
   for (Int_t i=1; i<=n; i++) {
      printf("%g %g\n",
             h->GetBinLowEdge(i)+h->GetBinWidth(i)/2,
             h->GetBinContent(i));
   }
}
