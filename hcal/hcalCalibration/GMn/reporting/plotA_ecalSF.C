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

void plotA_ecalSF( ){
  
  TCanvas *c1 = new TCanvas("c1","HCal eCal Samp Frac Results All kine",1600,1200);
  c1->Divide(3,2);

  for( Int_t i=0; i<nkine; i++ ){

    c1->cd(i+1);

    gStyle->SetOptStat(0);
    //Get pre-ecal samp frac
    TFile *f1 = TFile::Open(Form("eCalEOut_sbs%d.root",kIdx[i]));
    TH1D *h1= (TH1D*)f1->Get("hSampFrac");
    h1->GetXaxis()->SetRange(0,80);
    h1->SetLineWidth(3);
    h1->SetLineColor(kRed);
    h1->Draw("hist");
    //Get post-ecal samp frac
    TFile *f2 = TFile::Open(Form("qreplay_sbs%d.root",kIdx[i]));
    TH1D *h2 = (TH1D*)f2->Get("hSampFrac");
    h2->GetXaxis()->SetRange(0,80);
    h2->SetLineWidth(2);
    h2->SetLineColor(kBlack);
    h2->Draw("hist same");

    TF1 *fit1;
    h1->Fit("gaus","","",0.02,0.07);
    fit1 = h1->GetFunction("gaus");
    Double_t fit1m = fit1->GetParameter(1);

    TF1 *fit2;
    h2->Fit("gaus","","",0.05,0.09);
    fit2 = h2->GetFunction("gaus");
    Double_t fit2m = fit2->GetParameter(1);

    //Add a legend
    auto legend = new TLegend(0.43,0.7,0.89,0.89);
    legend->SetTextSize(0.03);
    legend->SetHeader(Form("HCal eCal Samp Frac SBS%d",kIdx[i]));
    legend->AddEntry(h1,Form("Pre-Cal, mean:%f",fit1m),"l");
    legend->AddEntry(h2,Form("Post-Cal, mean:%f",fit2m),"l");
    legend->Draw();

  }

  c1->SaveAs("/u/home/sseeds/Plots/ecalSF_allkine.png");

}
