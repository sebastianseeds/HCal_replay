//sseeds 4.17.23 - script to plot SF histos (all kinematics)

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
const Int_t pIdx[6] = {4,5,0,3,2,1};
const Int_t pIdxmc[6] = {2,1,5,4,3,0};

void plotAsame_SF( ){
  
  TCanvas *c1 = new TCanvas("c1","G4SBS Samp Frac Results",1600,1200);
  c1->Divide(2,1);

  TFile *fmc[nkine];
  TH1D *hmc[nkine];
  TF1 *fitmc[nkine];
  Double_t fitmcp1[nkine];

  auto hsmc = new THStack("hsmc","HCal Sampling Fraction (MC)");

  for( Int_t i=0; i<nkine; i++ ){

    fmc[i] = TFile::Open(Form("../MC/output/sfracMC_sbs%d.root",kIdx[i]));
    hmc[i]= (TH1D*)fmc[i]->Get("hSampFrac_clus");
    //hmc[i]->GetXaxis()->SetRange(0,0.17);
    hmc[i]->SetLineWidth(0);
    
    //hmc[i]->GetXaxis()->SetRangeUser(0.,0.2);

    hmc[i]->Fit("gaus","0","",0.04,0.1);
    fitmc[i] = hmc[i]->GetFunction("gaus");
    fitmcp1[i] = fitmc[i]->GetParameter(1)*100.;

    hmc[i]->SetTitle(Form("SBS%d, peak %0.2f%%",kIdx[i],fitmcp1[i]));

  }

  c1->cd(1);

  for( Int_t i=0; i<nkine; i++ ){
    Int_t j = pIdxmc[i];

    hsmc->Add(hmc[j]);
  }

  hsmc->Draw("pfc nostack");

  hsmc->GetXaxis()->SetTitle("E_{dep} / KE_{had}");
  hsmc->GetXaxis()->SetTitleOffset(1.4);

  c1->Modified();
  
  gPad->BuildLegend(0.45,0.65,0.89,0.89,"");


  gStyle->SetOptStat(0);
  //gStyle->SetPalette(kOcean);
  gStyle->SetPalette(kDeepSea);
  //Get MC samp frac all kine
  TFile *f[nkine];
  TH1D *h[nkine];
  TF1 *fit[nkine];
  Double_t fitp1[nkine];

  auto hs = new THStack("hs","HCal Sampling Fraction (Calibrated Data)");

  for( Int_t i=0; i<nkine; i++ ){

    f[i] = TFile::Open(Form("qreplay_sbs%d.root",kIdx[i]));
    h[i]= (TH1D*)f[i]->Get("hSampFrac");
    //h[i]->GetXaxis()->SetRange(0,0.17);
    h[i]->SetLineWidth(0);
    
    //h[i]->GetXaxis()->SetRangeUser(0.,0.2);

    h[i]->Fit("gaus","0","",0.04,0.1);
    fit[i] = h[i]->GetFunction("gaus");
    fitp1[i] = fit[i]->GetParameter(1)*100.;

    h[i]->SetTitle(Form("SBS%d, peak %0.2f%%",kIdx[i],fitp1[i]));


  }

  c1->cd(2);
  

  for( Int_t i=0; i<nkine; i++ ){
    Int_t j = pIdx[i];

    hs->Add(h[j]);
  }
  
  hs->Draw("pfc nostack");

  hs->GetXaxis()->SetTitle("E_{dep} / KE_{had}");
  hs->GetXaxis()->SetTitleOffset(1.4);

  c1->Modified();
  
  gPad->BuildLegend(0.45,0.65,0.89,0.89,"");

  //c1->SaveAs("/u/home/sseeds/Plots/MC_SF_allkine.png");

}
