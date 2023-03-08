//SSeeds 2.27.23 - Simple macro that overlays two similar histograms post e calibration. NOT WORKING

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


void plottwo( Int_t kine=4 ){
  
  // Declare paths to root analysis files
  string f1path = Form("eCalEOut_sbs%d.txt",kine);
  string f2path = Form("qreplay_sbs%d.txt",kine);
  
  // Declare outfile
  TFile *fout = new TFile( Form("plottwo_sbs%d.root", kine ), "RECREATE" );

  //TChain *C = new TChain("T");
  //TChain *D = new TChain("T");

  TFile *f1 = TFile.Open(f1path);
  TH1D *h1 = (TH1D*)f1->Get("hSampFrac");
  h1->SetLineColor(kRed);
  h1->SetLineWidth(2);
  h1->Draw();

  TFile *f2 = TFile.Open(f2path);
  TH1D *h2 = (TH1D*)f2->Get("hSampFrac");
  h2->SetLineColor(kBlue);
  h2->SetLineWidth(3);
  h2->Draw("same");

  fout->Write();

}

