//sseeds 4.24.23 - Script to extract sigma from titles on plots and put into tgraph

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

const Int_t kNcell = 288; //total kinematics
const Int_t kNrows = 24;
const Int_t kNcols = 12;

void plotAtdcsigma( Int_t kine = 8 ){
  
  gStyle->SetOptFit();
  gStyle->SetEndErrorSize(0);
  TCanvas *c1 = new TCanvas("c1","HCal TDC Timing Resolution by Block",1600,1200);
  c1->SetGrid();

  //auto mg = new TMultiGraph();

  TH1D *hmc[kNcell];
  TF1 *fit[kNcell];
  Double_t fitp1[kNcell];
  Double_t fitp2[kNcell];

  auto hsmc = new THStack("hsmc","HCal TDC Resolution");
  TFile *fmc = TFile::Open(Form("allc_replay_sbs%d.root",kine));;
  
  TH1D *hh = (TH1D*)fmc->Get("hpp_p");

  TH1D *hall = new TH1D("hall","All TDC, Safety Margin Cut", 100, 0, 3);

  Double_t fitl;
  Double_t fith;
  Double_t fits = 0.5;
  std::string title;
  std::string line;
  Double_t sigs[kNcell]={-1.};
  Double_t xidx[kNcell]={0.};
  for( Int_t i=1; i<kNcell; i++ ){

    Int_t row = i/kNcols;
    Int_t col = i%kNcols;

    xidx[i] = i;

    hmc[i]= (TH1D*)fmc->Get(Form("tcellslice_%d",i));

    
    title = hmc[i]->GetTitle();

    istringstream iss( title );
    
    Int_t mem = 0;
    while( iss >> line ){
      if( i==165 || i==265 || i==277 || i==288 )
	line = "0000000";
      if( mem==1 ){
	line.erase(0,6);
	sigs[i] = abs(std::stof(line));
	if( sigs[i]<1 || sigs[i]>2.4 )
	  sigs[i]=0.;
	//cout << i << " " << line << " " << sigs[i] << endl;
      }
      mem++;
    }

    if( row>1 && row<22 && col>1 && col<11 )
      hall->Fill(sigs[i]);

  }


  auto gr1 = new TGraph(kNcell,xidx,sigs);
  gr1->SetTitle("HCal TDC Sigma vs Block");
  gr1->SetMarkerColor(kViolet+4);
  gr1->SetMarkerStyle(33);
  gr1->SetMarkerSize(2);
  gr1->SetLineColor(kViolet+4);
  gr1->SetLineWidth(2);

  gr1->Draw("AP");
  
  TCanvas *c2 = new TCanvas("c2","HCal TDC Timing Resolution",1600,1200);
  c2->SetGrid();

  hall->Draw();
  

}
