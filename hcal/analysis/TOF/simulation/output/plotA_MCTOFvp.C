//sseeds 3.9.23 - script to plot TOF vs p histos (all kinematics, both protons an neutrons)

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
const Double_t TOFfitp_p[6][4] = {{73.7816, -35.8806, 12.90060, -1.60996},
				  {52.7549, -2.05791, 0.273797, -0.01309},
				  {54.6380, -1.97939, 0.231943, -0.00934},
				  {65.2843, -10.6005, 2.252520, -0.16723},
				  {60.6177, -18.2422, 5.123620, -0.49684},
				  //{71.7263, -32.1323, 41.57460, -1.18023}}; //Pol3 TOF fit params from proton MC data with SBS fields {30,85,100,70,70,70}
				  {63.9328, -22.3211, 6.673740, -0.68472}}; //Pol3 TOF fit params from neutron MC data with SBS fields {30,85,100,70,70,70}

const Double_t TOFfitp_n[6][4] = {{62.1158, -20.8086, 6.464510, -0.71342},
				  {55.8218, -3.84358, 0.617647, -0.03516},
				  {51.4030, -0.51868, 0.016076, 0.001009},
				  {64.1014, -9.30317, 1.818670, -0.12324},
				  {56.1597, -13.4154, 3.449660, -0.30995},
				  {63.9328, -22.3211, 6.673740, -0.68472}}; //Pol3 TOF fit params from neutron MC data with SBS fields {30,85,100,70,70,70}

void plotA_MCTOFvp( ){
  
  TCanvas *c1 = new TCanvas("c1","G4SBS TOF vs p, Protons",2200,1200);
  c1->Divide(3,2);
  TCanvas *c2 = new TCanvas("c2","G4SBS TOF vs p, Neutrons",2200,1200);
  c2->Divide(3,2);

  TFile *fout = new TFile( "aTOFfits.root", "RECREATE" );

  gStyle->SetOptStat(0);
  gStyle->SetPalette(53);

  //Get TOF v p for protons and send to canvas 1
  TFile *fmc[nkine];
  TH2D *hmc_p[nkine];
  for( Int_t i=0; i<nkine; i++ ){
    c1->cd(i+1);
    c1->SetGrid();

    TF1 *f = new TF1("f","pol3",1.5,9);
    f->SetParameters(TOFfitp_p[i][0],TOFfitp_p[i][1],TOFfitp_p[i][2],TOFfitp_p[i][3]);
    f->SetLineColor(kBlue);

    cout << TOFfitp_p[i][0] << " " << TOFfitp_p[i][1] << " " << TOFfitp_p[i][2] << " " << TOFfitp_p[i][3] << endl;

    fmc[i] = TFile::Open(Form("simTOFout_sbs%d.root",kIdx[i]));
    hmc_p[i]= (TH2D*)fmc[i]->Get("TOF_vs_pp");
    hmc_p[i]->SetTitle(Form("Proton Time of Flight vs Momentum, SBS%d",kIdx[i]));
    hmc_p[i]->GetXaxis()->SetTitle("p (GeV)");
    hmc_p[i]->GetYaxis()->SetTitle("Time of Flight (ns)");
    hmc_p[i]->Draw("hist colz");
    f->Draw("same");
    c1->Update();
  }

  // c1->Update();
  // c1->Modified();
  // c1->Write();
  cout << endl << endl;


  //Get TOF v p for neutrons and send to canvas 2
  TH2D *hmc_n[nkine];
  for( Int_t i=0; i<nkine; i++ ){

    c2->cd(i+1);
    c2->SetGrid();

    TF1 *g = new TF1("g","pol3",1.5,9);
    g->SetParameters(TOFfitp_n[i][0],TOFfitp_n[i][1],TOFfitp_n[i][2],TOFfitp_n[i][3]);
    g->SetLineColor(kBlue);

    cout << TOFfitp_n[i][0] << " " << TOFfitp_n[i][1] << " " << TOFfitp_n[i][2] << " " << TOFfitp_n[i][3] << endl;

    
    fmc[i] = TFile::Open(Form("simTOFout_sbs%d.root",kIdx[i]));
    hmc_n[i]= (TH2D*)fmc[i]->Get("TOF_vs_pn");
    hmc_n[i]->SetTitle(Form("Neutron Time of Flight vs Momentum, SBS%d",kIdx[i]));
    hmc_n[i]->GetXaxis()->SetTitle("p (GeV)");
    hmc_n[i]->GetYaxis()->SetTitle("Time of Flight (ns)");
    hmc_n[i]->Draw("hist colz");
    g->Draw("same");
    c2->Update();
  }

  // c2->Update();
  // c2->Modified();
  // c2->Write();

  fout->Write();

  c1->SaveAs("/u/home/sseeds/Plots/MC_TOFvp_proton_allkine.png");
  c2->SaveAs("/u/home/sseeds/Plots/MC_TOFvp_neutron_allkine.png");

}
