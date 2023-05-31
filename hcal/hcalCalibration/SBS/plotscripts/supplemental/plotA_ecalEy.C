//sseeds 3.9.23 - script to plot all E distributions by Y after calibrations and make tgrapherrors

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
#include "TMultiGraph.h"

const Int_t nkine = 6; //total kinematics
const Int_t kIdx[6] = {4,7,11,14,8,9}; //indexing kinematics for processing
const Int_t kNcols = 12; // Total number of HCal rows
const Int_t tfitmin = 30; //min elements per row to fit
const Double_t fitl[6] = {0.02,0.15,0.22,0.12,0.06,0.06};
const Double_t fitu[6] = {0.18,0.55,0.78,0.48,0.28,0.28};

void plotA_ecalEy( ){
  
  TCanvas *c1 = new TCanvas("c1","HCal Primary Cluster E All kine",2200,1200);
  c1->Divide(3,2);

  for( Int_t i=0; i<nkine; i++ ){

    c1->cd(i+1);

    gStyle->SetOptStat(0);

    TFile *f1 = TFile::Open(Form("qreplay_sbs%d.root",kIdx[i]));
    TH2D *h1 = (TH2D*)f1->Get("hHCALeY");
    
    //Make arrays for Ey tgraphs
    Double_t cellerr[kNcols] = {0.};
    Double_t Eycell[kNcols] = {0.};
    Double_t Eycval[kNcols] = {0.};
    Double_t Eycerr[kNcols] = {0.};
    TH1D *Eycellslice[kNcols];

    //Dispersive X counts across rows in HCal
    for( Int_t r=0; r<kNcols; r++ ){
      
      Eycell[r] = r;

      Eycellslice[r] = h1->ProjectionY(Form("Eycellslice_%d",r+1),r+1,r+1);
      Int_t sliceN = Eycellslice[r]->GetEntries();
      if( sliceN<tfitmin ){
	continue;
      }

      TF1 *fit1;
      Eycellslice[r]->Fit("gaus","Q","",fitl[i],fitu[i]);
      fit1 = Eycellslice[r]->GetFunction("gaus");
      Eycval[r] = fit1->GetParameter(1);
      Eycerr[r] = fit1->GetParameter(2);

    }

    // //Make graphs with errors for reporting
    TGraphErrors *gEy = new TGraphErrors( kNcols, Eycell, Eycval, cellerr, Eycerr );

    gEy->SetMinimum(0.0001);
    gEy->SetMaximum(1.);
    gEy->SetTitle(Form("HCal E vs X, SBS%d",kIdx[i]));
    gEy->GetXaxis()->SetTitle("Row");
    gEy->GetYaxis()->SetTitle("E (GeV)");
    gEy->SetMarkerStyle(20); // idx 20 Circles, idx 21 Boxes
    for( Int_t r=0; r<kNcols; r++ ){
      cout << Eycval[r] << endl;
      Double_t idx = r/kNcols;
      if( Eycval[r]==0. ) gEy->RemovePoint(idx);
    }

    gEy->Draw("AP");


  }

  c1->SaveAs("/u/home/sseeds/Plots/ecalEy_allkine.png");

}
