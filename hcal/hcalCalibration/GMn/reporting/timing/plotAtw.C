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

const Int_t kNcell = 288; //total cells
const Int_t kNcols = 12; //total columns

void plotAtw( Int_t kine = 8 ){
  
  gStyle->SetOptFit();

  TCanvas *TWt_top = new TCanvas("TWt_top","TWt_top",1600,1200);
  TCanvas *TWt_bot = new TCanvas("TWt_bot","TWt_bot",1600,1200);

  TWt_top->Divide(12,12);
  TWt_bot->Divide(12,12);

  TH2D *hmc[kNcell];

  TFile *fmc = TFile::Open(Form("allc_replay_sbs%d.root",kine));;
  
  for( Int_t c=0; c<kNcell; c++ ){

    Int_t col = c%kNcols;

    //Index through the canvas
    TWt_top->cd(c+1);
    if( c>=144 ){
      TWt_bot->cd(c-143);
      TWt_bot->SetGrid();
      gStyle->SetOptStat(0);
    }

    hmc[c]= (TH2D*)fmc->Get(Form("htdcVe_bl%d",c));

    hmc[c]->Draw();

  }



}
