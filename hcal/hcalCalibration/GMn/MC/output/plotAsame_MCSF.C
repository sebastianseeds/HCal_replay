//sseeds 3.9.23 - script to plot SF histos (all kinematics)

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
const Int_t pIdx[6] = {2,1,5,4,3,0};

void plotAsame_MCSF( ){
  
  TCanvas *c1 = new TCanvas("c1","G4SBS Samp Frac Results",1600,1200);

  gStyle->SetOptStat(0);
  //gStyle->SetPalette(kOcean);
  gStyle->SetPalette(kCividis);
  //Get MC samp frac all kine
  TFile *fmc[nkine];
  TH1D *hmc[nkine];
  TF1 *fitmc[nkine];
  Double_t fitp1[nkine];

  auto hs = new THStack("hs","HCal Sampling Fraction (MC)");

  for( Int_t i=0; i<nkine; i++ ){

    fmc[i] = TFile::Open(Form("sfracMC_sbs%d.root",kIdx[i]));
    hmc[i]= (TH1D*)fmc[i]->Get("hSampFrac_clus");
    hmc[i]->GetXaxis()->SetRange(0,80);
    hmc[i]->SetLineWidth(0);
    
    hmc[i]->Fit("gaus","0","",0.04,0.1);
    fitmc[i] = hmc[i]->GetFunction("gaus");
    fitp1[i] = fitmc[i]->GetParameter(1);

    hmc[i]->SetTitle(Form("SBS%d, %0.2f%%",kIdx[i],fitp1[i]*100.));


  }

  c1->cd();
  

  for( Int_t i=0; i<nkine; i++ ){
    Int_t j = pIdx[i];

    hs->Add(hmc[j]);
  }
  
  hs->Draw("pfc nostack");

  hs->GetXaxis()->SetTitle("E_{deposited} / (E_{beam} - E_{e'})");
  hs->GetXaxis()->SetTitleOffset(1.4);

  c1->Modified();
  
  gPad->BuildLegend(0.45,0.65,0.89,0.89,"");

  //c1->SaveAs("/u/home/sseeds/Plots/MC_SF_allkine.png");

}
