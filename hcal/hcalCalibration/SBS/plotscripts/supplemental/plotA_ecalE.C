//sseeds 3.9.23 - script to plot all SF histo overlays (before calibration / after calibration)

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
const Double_t fit1l[6] = {0.04,0.15,0.25,0.12,0.06,0.06};
const Double_t fit1u[6] = {0.1,0.35,0.45,0.25,0.13,0.13};
const Double_t fit2l[6] = {0.06,0.25,0.35,0.2,0.11,0.11};
const Double_t fit2u[6] = {0.16,0.52,0.75,0.38,0.23,0.23};

void plotA_ecalE( ){
  
  TCanvas *c1 = new TCanvas("c1","HCal Primary Cluster E All kine",1600,1200);
  c1->Divide(3,2);

  for( Int_t i=0; i<nkine; i++ ){

    c1->cd(i+1);

    gStyle->SetOptStat(0);
    //Get pre-ecal energy dist
    TFile *f1 = TFile::Open(Form("eCalEOut_sbs%d.root",kIdx[i]));
    TH1D *h1= (TH1D*)f1->Get("hHCALe");
    //h1->GetXaxis()->SetRange(0,400);
    h1->SetLineWidth(3);
    h1->SetLineColor(kRed);
    //h1->Draw("hist");
    //Get post-ecal energy dist
    TFile *f2 = TFile::Open(Form("qreplay_sbs%d.root",kIdx[i]));
    TH1D *h2 = (TH1D*)f2->Get("hHCALe");
    //h2->GetXaxis()->SetRange(0,400);
    h2->SetLineWidth(2);
    h2->SetLineColor(kBlack);
    //h2->Draw("hist same");

    TF1 *fit1;
    h1->Fit("gaus","","",fit1l[i],fit1u[i]);
    fit1 = h1->GetFunction("gaus");
    Double_t fit1m = fit1->GetParameter(1);

    TF1 *fit2;
    h2->Fit("gaus","","",fit2l[i],fit2u[i]);
    fit2 = h2->GetFunction("gaus");
    Double_t fit2m = fit2->GetParameter(1);
    Double_t fit2s = fit2->GetParameter(2);
    
    h1->Draw("hist");

    h2->Draw("hist same");

    Double_t res = fit2s/fit2m*100;

    //Add a legend
    auto legend = new TLegend(0.43,0.7,0.89,0.89);
    legend->SetTextSize(0.03);
    legend->SetHeader(Form("HCal PCluster E SBS%d",kIdx[i]));
    legend->AddEntry(h1,"Pre-Cal","l");
    legend->AddEntry(h2,Form("Post-Cal, E Res:%d%%",(Int_t)res),"l");
    legend->Draw();

  }

  c1->SaveAs("/u/home/sseeds/Plots/ecalE_allkine.png");

}
