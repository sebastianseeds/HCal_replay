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

void plotA_MCTOFvp( ){
  
  TCanvas *c1 = new TCanvas("c1","G4SBS TOF vs p, Protons",1600,1200);
  c1->Divide(3,2);
  TCanvas *c2 = new TCanvas("c2","G4SBS TOF vs p, Neutrons",1600,1200);
  c2->Divide(3,2);

  gStyle->SetOptStat(0);

  //Get TOF v p for protons and send to canvas 1
  TFile *fmc[nkine];
  TH2D *hmc_p[nkine];
  //TF1 *fitmc[nkine];
  //Double_t fitp1[nkine];
  for( Int_t i=0; i<nkine; i++ ){
    c1->cd(i+1);

    fmc[i] = TFile::Open(Form("simTOFout_sbs%d.root",kIdx[i]));
    hmc_p[i]= (TH2D*)fmc[i]->Get("TOF_vs_pp");
    //hmc[i]->GetXaxis()->SetRange(0,80);
    //hmc[i]->SetLineWidth(3);
    //hmc[i]->SetLineColor(i+40);
    
    //hmc[i]->Fit("gaus","","",0.04,0.1);
    //fitmc[i] = hmc[i]->GetFunction("gaus");
    //fitp1[i] = fitmc[i]->GetParameter(1);

    //hmc[i]->Draw("hist same");

    // if( i==0 ){
    //   hmc[i]->Draw("hist");
    // }else{
    //   hmc[i]->Draw("hist same");
    // }

    hmc_p[i]->Draw("hist");
  }

  //Add a legend
  c1->cd(3);
  auto legend = new TLegend(0.45,0.65,0.89,0.89);
  legend->SetTextSize(0.03);
  legend->SetHeader("HCal MC TOF vs p Results");
  for( Int_t i=0; i<nkine; i++ ){

    legend->AddEntry(hmc[i],Form("SBS%d",kIdx[i]),"l");
  }
  legend->Draw();

  //Get TOF v p for neutrons and send to canvas 2
  TH2D *hmc_n[nkine];
  for( Int_t i=0; i<nkine; i++ ){
    c2->cd(i+1);

    fmc[i] = TFile::Open(Form("simTOFout_sbs%d.root",kIdx[i]));
    hmc_n[i]= (TH2D*)fmc[i]->Get("TOF_vs_pn");
    hmc_n[i]->Draw("hist");
  }

  //Add a legend
  c1->cd(3);
  auto legend = new TLegend(0.45,0.65,0.89,0.89);
  legend->SetTextSize(0.03);
  legend->SetHeader("HCal MC TOF vs p Results");
  for( Int_t i=0; i<nkine; i++ ){

    legend->AddEntry(hmc[i],Form("SBS%d",kIdx[i]),"l");
  }
  legend->Draw();

  c1->SaveAs("/u/home/sseeds/Plots/MC_TOFvp_proton_allkine.pdf");

}
