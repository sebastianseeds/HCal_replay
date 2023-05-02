//sseeds 4.17.23 - script to plot SF results on TMultigraph (all kinematics)

#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <TChain.h>
#include <TSystem.h>
#include <TStopwatch.h>
#include <unistd.h>
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TCut.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TLegend.h"
#include "TVector3.h"
#include "TGraphErrors.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TF1.h"

const Int_t nkine = 6; //total kinematics
const Int_t kIdx[6] = {4,7,11,14,8,9}; //indexing kinematics for processing
// const Double_t kDIdx[6] = {4.2,7.2,11.2,14.2,8.2,9.2}; //indexing kinematics for processing
// const Double_t kDIdxmc[6] = {4.,7.,11.,14.,8.,9.}; //indexing kinematics for processing
const Double_t kDIdx[6] = {1.2,2.2,3.2,4.2,5.2,6.2}; //indexing kinematics for processing
const Double_t kDIdxmc[6] = {1.,2.,3.,4.,5.,6.}; //indexing kinematics for processing
// const Double_t kDIdx[6] = {1.65,5.25,7.25,4.05,2.45,2.55}; //nucleon central KE
// const Double_t kDIdxmc[6] = {1.6,5.2,7.2,4.0,2.4,2.5}; //nucleon KE
// const Double_t kDIdx[6] = {0.106,0.376,0.527,0.283,0.166,0.162}; //via manual fits to E spectra
// const Double_t kDIdxmc[6] = {0.125,0.419,0.581,0.320,0.199,0.196}; //via manual fits to E spectra
const Int_t pIdx[6] = {4,5,0,3,2,1};
const Int_t pIdxmc[6] = {2,1,5,4,3,0};
const Double_t kXerr[6] = {0.,0.,0.,0.,0.,0.}; //via manual fits to E spectra
const Double_t kXerrmc[6] = {0.,0.,0.,0.,0.,0.}; //via manual fits to E spectra
// const Double_t kXerr[6] = {0.069,0.165,0.221,0.103,0.085,0.060}; //via manual fits to E spectra
// const Double_t kXerrmc[6] = {0.045,0.103,0.114,0.078,0.060,0.055}; //via manual fits to E spectra

void plotSF_MCdirectComparison( ){
  
  gStyle->SetOptFit();
  gStyle->SetEndErrorSize(0);
  TCanvas *c1 = new TCanvas("c1","HCal Sampling Fraction MC Comparison",1600,1200);
  c1->SetGrid();

  auto mg = new TMultiGraph();

  TFile *fmc[nkine];
  TH1D *hmc[nkine];
  TF1 *fitmc[nkine];
  Double_t fitmcp1[nkine];
  Double_t fitmcp2[nkine];

  for( Int_t i=0; i<nkine; i++ ){

    fmc[i] = TFile::Open(Form("../MC/output/sfracMC_sbs%d.root",kIdx[i]));
    hmc[i]= (TH1D*)fmc[i]->Get("hSampFrac_clus");
    hmc[i]->SetLineWidth(0);
    
    hmc[i]->Fit("gaus","0","",0.04,0.1);
    fitmc[i] = hmc[i]->GetFunction("gaus");
    fitmcp1[i] = fitmc[i]->GetParameter(1)*100.; //convert to percent
    fitmcp2[i] = fitmc[i]->GetParameter(2)*100.; //convert to percent

  }

  auto gr1 = new TGraphErrors(nkine,kDIdxmc,fitmcp1,kXerrmc,fitmcp2);
  gr1->SetTitle("MC");
  gr1->SetMarkerColor(kAzure-6);
  gr1->SetMarkerStyle(33);
  gr1->SetMarkerSize(2);
  gr1->SetLineColor(kAzure-6);
  gr1->SetLineWidth(2);
  mg->Add(gr1);


  //Get MC samp frac all kine
  TFile *f[nkine];
  TH1D *h[nkine];
  TF1 *fit[nkine];
  Double_t fitp1[nkine];
  Double_t fitp2[nkine];

  for( Int_t i=0; i<nkine; i++ ){

    f[i] = TFile::Open(Form("qreplay_sbs%d.root",kIdx[i]));
    h[i]= (TH1D*)f[i]->Get("hSampFrac");
    h[i]->SetLineWidth(0);
    
    h[i]->Fit("gaus","0","",0.04,0.1);
    fit[i] = h[i]->GetFunction("gaus");
    fitp1[i] = fit[i]->GetParameter(1)*100.; //convert to percent
    fitp2[i] = fit[i]->GetParameter(2)*100.; //convert to percent

  }

  auto gr2 = new TGraphErrors(nkine,kDIdx,fitp1,kXerr,fitp2);
  gr2->SetTitle("Data (Calibrated)");
  gr2->SetMarkerColor(kOrange-1);
  gr2->SetMarkerStyle(33);
  gr2->SetMarkerSize(2);
  gr2->SetLineColor(kOrange-1);
  gr2->SetLineWidth(2);
  mg->Add(gr2);

  mg->SetTitle("HCal Sampling Fraction");
  //mg->GetXaxis()->SetTitle("HCal E (GeV)");
  //mg->GetXaxis()->SetTitle("Kinematic");
  mg->GetYaxis()->SetTitle("HCal Cluster E / Hadron KE (%)");
  mg->Draw("AP");

  mg->GetYaxis()->SetRangeUser(0.,15.);
  //mg->GetXaxis()->SetNdivisions(12,3,15,kTRUE);
  mg->GetXaxis()->SetNdivisions(6,"I");
  
  //cout << mg->GetNBins() << endl;

  for( Int_t i=0; i<nkine; i++ ){
    mg->GetXaxis()->SetBinLabel( i*17.5+7, Form("SBS-%d",kIdx[i]));
  }
  mg->GetXaxis()->LabelsOption("h");
  mg->GetXaxis()->SetLabelSize(0.05);

  c1->Modified();

  c1->BuildLegend();


}
