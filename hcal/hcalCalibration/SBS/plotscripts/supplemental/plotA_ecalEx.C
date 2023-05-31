//sseeds 3.9.23 - script to plot all E distributions by X after calibrations and make tgrapherrors

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
const Int_t kNrows = 24; // Total number of HCal rows
const Int_t tfitmin = 30; //min elements per row to fit
const Double_t fitl[6] = {0.02,0.15,0.22,0.12,0.06,0.06};
const Double_t fitu[6] = {0.18,0.55,0.78,0.48,0.28,0.28};

void plotA_ecalEx( ){
  
  TCanvas *c1 = new TCanvas("c1","HCal Primary Cluster E All kine",2200,1200);
  c1->Divide(3,2);

  for( Int_t i=0; i<nkine; i++ ){

    c1->cd(i+1);

    gStyle->SetOptStat(0);
    gStyle->SetPalette(53);

    TFile *f1 = TFile::Open(Form("qreplay_sbs%d.root",kIdx[i]));
    TH2D *h1 = (TH2D*)f1->Get("hHCALeX");
    
    //Make arrays for Ex tgraphs
    Double_t cellerr[kNrows] = {0.};
    Double_t Excell[kNrows] = {0.};
    Double_t Excval[kNrows] = {0.};
    Double_t Excerr[kNrows] = {0.};
    TH1D *Excellslice[kNrows];

    //Dispersive X counts across rows in HCal
    for( Int_t r=0; r<kNrows; r++ ){
      
      Excell[r] = r;

      Excellslice[r] = h1->ProjectionY(Form("Excellslice_%d",r+1),r+1,r+1);
      Int_t sliceN = Excellslice[r]->GetEntries();
      if( sliceN<tfitmin ){
	continue;
      }

      TF1 *fit1;
      Excellslice[r]->Fit("gaus","Q","",fitl[i],fitu[i]);
      fit1 = Excellslice[r]->GetFunction("gaus");
      Excval[r] = fit1->GetParameter(1);
      Excerr[r] = fit1->GetParameter(2);

    }

    // //Make graphs with errors for reporting
    TGraphErrors *gEx = new TGraphErrors( kNrows, Excell, Excval, cellerr, Excerr );

    gEx->SetMinimum(0.);
    gEx->SetMaximum(1.);
    gEx->SetTitle(Form("HCal E vs X, SBS%d",kIdx[i]));
    gEx->GetXaxis()->SetTitle("Row");
    gEx->GetYaxis()->SetTitle("E (GeV)");
    gEx->SetMarkerStyle(20); // idx 20 Circles, idx 21 Boxes
    for( Int_t r=0; r<kNrows; r++ ){
      cout << Excval[r] << endl;
      Double_t idx = r/kNrows;
      if( Excval[r]==0 ) gEx->RemovePoint(idx);
    }

    gEx->Draw("AP");


  }

  c1->SaveAs("/u/home/sseeds/Plots/ecalEx_allkine.png");

}
