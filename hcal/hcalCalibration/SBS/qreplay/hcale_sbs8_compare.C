//sseeds 8.2.23 script to compare output of qreplay hcal energy distributions against MC

#include <vector>
#include <iostream>
#include <algorithm>
#include "TCut.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TTreeFormula.h"
#include "TChain.h"
#include "TLatex.h"
#include "TH1D.h"
#include "TStopwatch.h"
#include "../include/sbs.h"

const Double_t E_hard_ulim = 0.9;

//Compares the energy spectra from qreplay_standalone.C root file and mc data
void hcale_sbs8_compare( const char *experiment = "gmn", Int_t config = 9, Int_t pass = 1, Double_t FWHM = 0.03){

  //hardcoded for my environment for now
  std::string outdir_path = "/lustre19/expphy/volatile/halla/sbs/seeds";

  //set up relative paths
  std::string plotdir = Form("../quality_plots/%s/conf%d%s/",experiment,config,"_comp");
  std::string qr_path = outdir_path + Form("/hcal_calibrations/qreplay/qreplay_from_gmn_8_1_to_%s_%d_%d.root",experiment,config,pass); //only looking at sbs8 compared with each other kinematic's full data set (LH2+LD2)
  std::string digmc_path = outdir_path + Form("/hcal_calibrations/MC/hcalE_idx_mc_dig_%s_conf%d.root",experiment,config);
  std::string energy_dmc_path = plotdir + "qr_energy_dmc.png";
  std::string sf_dmc_path = plotdir + "qr_sf_dmc.png";

  std::string fout_path = outdir_path + Form("/hcal_calibrations/qreplay/ecomp_from_gmn_8_1_to_%s_%d_%d.root",experiment,config,pass);

  TFile *fout = new TFile( fout_path.c_str(), "RECREATE" );

  TCanvas *c1 = new TCanvas("c1","mc/data energy",1600,1200);
  c1->cd();
  gStyle->SetOptStat(0);

  TFile *qr_file = new TFile(qr_path.c_str(), "READ"); // Open the ROOT file in read mode
  TH1D *hE = (TH1D*)qr_file->Get("hE_new");
  Int_t hE_entries = hE->GetEntries();

  //get max value of fit to peak for hE
  Int_t hE_binMax = hE->GetMaximumBin();
  Double_t hE_binCenter = hE->GetBinCenter( hE_binMax );
  Double_t hE_fitLL = hE_binCenter - FWHM;
  Double_t hE_fitUL = hE_binCenter + FWHM;
  
  // Fit a Gaussian to this projection
  TF1 *hE_gfit = new TF1("hE_gfit", "gaus");
  hE->Fit("hE_gfit", "QN", "", hE_fitLL, hE_fitUL ); // "Q" for quiet mode, "N" for no draw
  Double_t hE_fit_max = hE_gfit->GetMaximum(hE_fitLL,hE_fitUL);

  //cout << "hE fit max " << hE_fit_max << endl;

  TFile *digmc_file = new TFile(digmc_path.c_str(), "READ"); // Open the ROOT file in read mode
  TH1D *hEmc = (TH1D*)digmc_file->Get("hE_mc");
  Int_t hEmc_entries = hEmc->GetEntries();
  TH1D *hEmc_nocut =  (TH1D*)digmc_file->Get("hE_mc_nocut");
  Int_t hEmc_nocut_entries = hEmc_nocut->GetEntries();

  //get max value of fit to peak for hEmc
  Double_t hEmc_max = hEmc->GetMaximum();

  Int_t hEmc_binMax = hEmc->GetMaximumBin();
  Double_t hEmc_binCenter = hEmc->GetBinCenter( hEmc_binMax );
  Double_t hEmc_fitLL = hEmc_binCenter - FWHM;
  Double_t hEmc_fitUL = hEmc_binCenter + FWHM;
  
  // Fit a Gaussian to this projection
  TF1 *hEmc_gfit = new TF1("hEmc_gfit", "gaus");
  hEmc->Fit("hEmc_gfit", "QN", "", hEmc_fitLL, hEmc_fitUL ); // "Q" for quiet mode, "N" for no draw
  Double_t hEmc_fit_max = hEmc_gfit->GetMaximum(hEmc_fitLL,hEmc_fitUL);

  //get max value of fit to peak for hEmc_nocut
  Double_t hEmc_nocut_max = hEmc_nocut->GetMaximum();

  Int_t hEmc_nocut_binMax = hEmc_nocut->GetMaximumBin();
  Double_t hEmc_nocut_binCenter = hEmc_nocut->GetBinCenter( hEmc_nocut_binMax );
  Double_t hEmc_nocut_fitLL = hEmc_nocut_binCenter - FWHM;
  Double_t hEmc_nocut_fitUL = hEmc_nocut_binCenter + FWHM;
  
  // Fit a Gaussian to this projection
  TF1 *hEmc_nocut_gfit = new TF1("hEmc_nocut_gfit", "gaus");
  hEmc_nocut->Fit("hEmc_nocut_gfit", "QN", "", hEmc_nocut_fitLL, hEmc_nocut_fitUL ); // "Q" for quiet mode, "N" for no draw
  Double_t hEmc_nocut_fit_max = hEmc_nocut_gfit->GetMaximum(hEmc_nocut_fitLL,hEmc_nocut_fitUL);

  // Obtain scale factor and scale other histograms
  Double_t escale_factor = hEmc_fit_max / hE_fit_max;
  hE->Scale(escale_factor);
  Double_t escale_nocut_factor = hEmc_fit_max / hEmc_nocut_fit_max;
  hEmc_nocut->Scale(escale_nocut_factor);
  Double_t hE_max = hE->GetMaximum();

  //cout << hEmc_max << ":" << hE_max << endl;

  Double_t eyrangescale = hEmc_max*1.1;
  if( hEmc_max < hE_max )
    eyrangescale = hE_max*1.1;

  hEmc->SetLineWidth(1);
  hEmc->SetLineColor(kGreen);
  hEmc->SetFillStyle(3004);
  hEmc->SetFillColor(kGreen);
  hEmc->GetXaxis()->SetRangeUser(0,E_hard_ulim);
  hEmc->GetYaxis()->SetRangeUser(0,eyrangescale);
  hEmc->GetXaxis()->SetTitle("GeV");
  hEmc->Draw("hist");

  hEmc_nocut->SetLineWidth(1);
  hEmc_nocut->SetLineColor(kCyan);
  hEmc_nocut->SetFillStyle(3005);
  hEmc_nocut->SetFillColor(kCyan);
  hEmc_nocut->GetXaxis()->SetRangeUser(0,E_hard_ulim);
  hEmc_nocut->GetYaxis()->SetRangeUser(0,eyrangescale);
  hEmc_nocut->GetXaxis()->SetTitle("GeV");
  hEmc_nocut->Draw("hist same");
    
  hE->SetLineWidth(2);
  hE->SetLineColor(kBlack);
  hE->GetXaxis()->SetRangeUser(0,E_hard_ulim);
  hE->GetXaxis()->SetTitle("GeV");
  hE->Draw("hist same");
  
  //Add a legend
  auto elegend = new TLegend(0.43,0.7,0.89,0.89);
  elegend->SetTextSize(0.03);
  elegend->SetHeader(Form("SBS%d, E",config));
  elegend->AddEntry(hEmc,Form("MC, %d entries",hEmc_entries),"l");
  elegend->AddEntry(hEmc_nocut,Form("MC no cuts (scaled), %d entries",hEmc_nocut_entries),"l");
  elegend->AddEntry(hE,Form("Data (scaled), %d entries",hE_entries),"l");
  elegend->Draw();

  c1->SaveAs(energy_dmc_path.c_str());
  //c1->Write();

  //Sampling fraction
  TCanvas *c2 = new TCanvas("c2","mc/data sampling fraction",1600,1200);
  c2->cd();
  gStyle->SetOptStat(0);

  //get mc sf
  //TFile *fSFmc = TFile::Open(digmc_path.c_str());
  TH1D *hSFmc = (TH1D*)digmc_file->Get("hSF_mc");
  Int_t hSFmc_entries = hSFmc->GetEntries();
  TH1D *hSFmc_nocut = (TH1D*)digmc_file->Get("hSF_mc_nocut");
  Int_t hSFmc_nocut_entries = hSFmc_nocut->GetEntries();

  //get max value of hSFmc histogram
  Double_t hSFmc_max = hSFmc->GetMaximum();

  Int_t hSFmc_binMax = hSFmc->GetMaximumBin();
  Double_t hSFmc_binCenter = hSFmc->GetBinCenter( hSFmc_binMax );
  Double_t hSFmc_fitLL = hSFmc_binCenter - FWHM;
  Double_t hSFmc_fitUL = hSFmc_binCenter + FWHM;
  
  // Fit a Gaussian to this projection
  TF1 *hSFmc_gfit = new TF1("hSFmc_gfit", "gaus");
  hSFmc->Fit("hSFmc_gfit", "QN", "", hSFmc_fitLL, hSFmc_fitUL ); // "Q" for quiet mode, "N" for no draw
  Double_t hSFmc_fit_max = hSFmc_gfit->GetMaximum(hSFmc_fitLL,hSFmc_fitUL);

  //get max value of hSFmc_nocut histogram
  Double_t hSFmc_nocut_max = hSFmc_nocut->GetMaximum();

  Int_t hSFmc_nocut_binMax = hSFmc_nocut->GetMaximumBin();
  Double_t hSFmc_nocut_binCenter = hSFmc_nocut->GetBinCenter( hSFmc_nocut_binMax );
  Double_t hSFmc_nocut_fitLL = hSFmc_nocut_binCenter - FWHM;
  Double_t hSFmc_nocut_fitUL = hSFmc_nocut_binCenter + FWHM;
  
  // Fit a Gaussian to this projection
  TF1 *hSFmc_nocut_gfit = new TF1("hSFmc_nocut_gfit", "gaus");
  hSFmc_nocut->Fit("hSFmc_nocut_gfit", "QN", "", hSFmc_nocut_fitLL, hSFmc_nocut_fitUL ); // "Q" for quiet mode, "N" for no draw
  Double_t hSFmc_nocut_fit_max = hSFmc_nocut_gfit->GetMaximum(hSFmc_nocut_fitLL,hSFmc_nocut_fitUL);

  //get data sf
  //TFile *qr_file = new TFile(qr_path.c_str(), "READ"); // Open the ROOT file in read mode
  TH1D *hSF = (TH1D*)qr_file->Get("hSF_new");
  Int_t hSF_entries = hSF->GetEntries();

  //get max value of fit
  Int_t hSF_binMax = hSF->GetMaximumBin();
  Double_t hSF_binCenter = hSF->GetBinCenter( hSF_binMax );
  Double_t hSF_fitLL = hSF_binCenter - FWHM;
  Double_t hSF_fitUL = hSF_binCenter + FWHM;
  
  // Fit a Gaussian to this projection
  TF1 *hSF_gfit = new TF1("hSF_gfit", "gaus");
  hSF->Fit("hSF_gfit", "QN", "", hSF_fitLL, hSF_fitUL ); // "Q" for quiet mode, "N" for no draw
  Double_t hSF_fit_max = hSF_gfit->GetMaximum(hSF_fitLL,hSF_fitUL);

  //get scale factor
  Double_t sfscale_factor = hSFmc_fit_max / hSF_fit_max;
  hSF->Scale(sfscale_factor);
  Double_t sfscale_nocut_factor = hSFmc_fit_max / hSFmc_nocut_fit_max;
  hSFmc_nocut->Scale(sfscale_nocut_factor);
  Double_t hSF_max = hSF->GetMaximum();

  Double_t sfyrangescale = hSFmc_max*1.1;
  if( hSFmc_max < hSF_max )
    sfyrangescale = hSF_max*1.1;

  hSFmc->SetLineWidth(1);
  hSFmc->SetLineColor(kOrange);
  hSFmc->SetFillStyle(3004);
  hSFmc->SetFillColor(kOrange);
  hSFmc->GetXaxis()->SetRangeUser(0,0.5);
  hSFmc->GetYaxis()->SetRangeUser(0,sfyrangescale);
  hSFmc->GetXaxis()->SetTitle("E_{clus}/KE");
  hSFmc->Draw("hist");

  hSFmc_nocut->SetLineWidth(1);
  hSFmc_nocut->SetLineColor(kRed);
  hSFmc_nocut->SetFillStyle(3005);
  hSFmc_nocut->SetFillColor(kRed);
  hSFmc_nocut->GetXaxis()->SetRangeUser(0,0.5);
  hSFmc_nocut->GetYaxis()->SetRangeUser(0,sfyrangescale);
  hSFmc_nocut->GetXaxis()->SetTitle("E_{clus}/KE");
  hSFmc_nocut->Draw("hist same");

  hSF->SetLineWidth(2);
  hSF->SetLineColor(kBlack);
  hSF->GetXaxis()->SetRangeUser(0,0.5);
  hSF->GetXaxis()->SetTitle("E_{clus}/KE");
  hSF->Draw("hist same");

  //Add a legend
  auto sflegend = new TLegend(0.43,0.7,0.89,0.89);
  sflegend->SetTextSize(0.03);
  sflegend->SetHeader(Form("SBS%d, SF",config));
  sflegend->AddEntry(hSFmc,Form("MC, %d entries",hSFmc_entries),"l");
  sflegend->AddEntry(hSFmc_nocut,Form("MC no cuts (scaled), %d entries",hSFmc_nocut_entries),"l");
  sflegend->AddEntry(hSF,Form("Data (scaled), %d entries",hSF_entries),"l");
  //sflegend->AddEntry(hSFmc,"MC","l");
  //sflegend->AddEntry(hSF,"Data (Scaled)","l");
  sflegend->Draw();

  c2->SaveAs(sf_dmc_path.c_str());
  //c2->Write();

  fout->Write();

}
