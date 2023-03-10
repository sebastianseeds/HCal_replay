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

void plotA_MCSF( ){
  
  TCanvas *c1 = new TCanvas("c1","G4SBS Samp Frac Results",1600,1200);
  c1->Divide(3,2);

  gStyle->SetOptStat(0);
  //Get MC samp frac all kine
  TFile *fmc[nkine];
  TH1D *hmc[nkine];
  TF1 *fitmc[nkine];
  Double_t fitp1[nkine];
  for( Int_t i=0; i<nkine; i++ ){
    c1->cd(i+1);

    fmc[i] = TFile::Open(Form("sfracMC_sbs%d.root",kIdx[i]));
    hmc[i]= (TH1D*)fmc[i]->Get("hSampFrac_clus");
    hmc[i]->GetXaxis()->SetRange(0,80);
    hmc[i]->SetLineWidth(3);
    hmc[i]->SetLineColor(i+40);
    
    hmc[i]->Fit("gaus","","",0.04,0.1);
    fitmc[i] = hmc[i]->GetFunction("gaus");
    fitp1[i] = fitmc[i]->GetParameter(1);

    //hmc[i]->Draw("hist same");

    // if( i==0 ){
    //   hmc[i]->Draw("hist");
    // }else{
    //   hmc[i]->Draw("hist same");
    // }

    hmc[i]->Draw("hist");
  }

  //Add a legend
  c1->cd(3);
  auto legend = new TLegend(0.45,0.65,0.89,0.89);
  legend->SetTextSize(0.03);
  legend->SetHeader("HCal MC Samp Frac Results");
  for( Int_t i=0; i<nkine; i++ ){

    legend->AddEntry(hmc[i],Form("SBS%d, mean:%f",kIdx[i],fitp1[i]),"l");
  }
  legend->Draw();

  c1->SaveAs("/u/home/sseeds/Plots/MC_SF_allkine.png");

}
