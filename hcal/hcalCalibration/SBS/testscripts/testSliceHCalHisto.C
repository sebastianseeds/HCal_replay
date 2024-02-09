//SSeeds Test script to verify the functionality of util.C function sliceHCalIDHisto.C

#include <vector>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <string>
#include "TCut.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TTreeFormula.h"
#include "TChain.h"
#include "TLatex.h"
#include "TH1D.h"
#include "TStopwatch.h"
#include "../include/sbs.h"

bool verb = true;

///////////////////////////////////////////////////////
//Manual overide/tune information
bool channel_override = true;
bool channel_override_tdc = false;
Int_t corrChan = 264; //Channel to correct
Int_t corrChan_tdc = -1; //Channel to correct
//Double_t corrChanOffset = -249.948; //sbs4 (-249.448, -249.948), sbs7 -239.234, sbs 11 (-271.841,-267.471,-271.841), sbs14 (-243.765,-243.765), sbs8 -253.499, sbs9 (-250.951)
Double_t corrChanOffset_tdc = 0.;
Double_t binXmax = 1e38;
Double_t binXmin = -1e38;
Double_t binXmax_tdc = 1e38;
Double_t binXmin_tdc = -1e38;
///////////////////////////////////////////////////////

//Main. If no arguments are passed to run_b and run_e, alignment over all runs will proceed.
void testSliceHCalHisto( const char *experiment = "gmn", 
			 int config = 9, 
			 int pass = 1, 
			 int run_b = 0, 
			 int run_e = 0, 
			 int run_exclude_b = 0, 
			 int run_exclude_e = 0, 
			 bool best_clus = false, 
			 double adct_target = 0.,
			 double tdc_target = 0.,
			 double corrChanOffset = -250.951){

  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  // Get the date
  std::string date = util::getDate();

  // Standardize pass 0/1, differences handled during histogram creation
  if( pass==0 )
    pass=1;

  // set style
  gStyle->SetOptFit(0);
  gStyle->SetEndErrorSize(0);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(55);

  std::string out_dir = gSystem->Getenv("OUT_DIR");
  std::string in_path = out_dir + Form("/hcal_calibrations/pass%d/diagnostics/timingdiag_%s_conf%d_pass%d_%d_to_%d_exclude_%d_to_%d.root",pass,experiment,config,pass,run_b,run_e,run_exclude_b,run_exclude_e);

  std::ifstream file(in_path.c_str());
  if (!file.good()){
    cout << "ERROR: No file located at " << in_path << endl;
    return;
  }

  //Path to outfile
  std::string out_path = "slice_out.root";

  // create output analysis file
  TFile *fout = new TFile( out_path.c_str(), "RECREATE" );

  // Get information from .csv files
  std::string struct_dir = Form("../config/%s_p%d/",experiment,pass);

  vector<caldiag> cut;
  std::string targ = "lh2";
  Int_t mag = 70;
  util::ReadDiagnosticCutList(struct_dir,experiment,config,targ,mag,verb,cut); //needed only for min ev
  Int_t min_epc = 300; //TEST ONLY
  Double_t atime_fwhm = cut[0].atime_sig;
  

  // open file, get histograms
  TFile *f1 = TFile::Open(in_path.c_str());

  TH2D *hatimeid = (TH2D*)f1->Get("hatimeid");
  TH2D *hatimerow = (TH2D*)f1->Get("hatimerow");
  TH2D *htimeid = (TH2D*)f1->Get("htimeid");

  /////////////////////////////////
  //Fits to tight elastic cut data
    
  vector<Double_t> c1cell;
  vector<Double_t> c1mean;
  vector<Double_t> c1err;
    
  //Clone the histogram which corresponds to the cluster method chosen
  TH2D *hadctAlign = (TH2D*)hatimeid->Clone("hadctAlign");

  Int_t c1binsX = hadctAlign->GetNbinsX();

  TCanvas *c1v1 = new TCanvas("c1v1","ATime Diff fits tight cuts HCal Top Half",1600,1200);
  TCanvas *c1v2 = new TCanvas("c1v2","ATime Diff fits tight cuts HCal Bottom Half",1600,1200);
   
  util::sliceHCalIDHisto(hadctAlign, c1binsX, atime_fwhm, min_epc, c1v1, c1v2, c1cell, c1mean, c1err, binXmin, binXmax );
  c1v1->Write();
  c1v2->Write();

  // TCanvas *c1 = new TCanvas("c1","ATime Diff vs ID Tight Cut",1600,1200);
  // c1->Divide(1,2);
  // c1->SetGrid();
  // c1->cd(1);

  // // Convert vectors to arrays
  // Double_t* c1x = &c1cell[0];  
  // Double_t* c1y = &c1mean[0];
  // Double_t* c1ey = &c1err[0];
  
  // //hadctAlign->GetYaxis()->SetRangeUser(-10.,10.);
  // hadctAlign->GetYaxis()->SetTitle("ns");
  // hadctAlign->GetXaxis()->SetTitle("");
  // hadctAlign->Draw("colz");
    
  // //Make graphs with errors for reporting
  // TGraphErrors *c1g = new TGraphErrors( c1cell.size(), c1x, c1y, 0, c1ey );
    
  // c1g->SetTitle(Form("HCal adctime vs id, SBS%d",config));
  // c1g->SetMarkerStyle(21); // idx 20 Circles, idx 21 Boxes
  // c1g->SetMarkerColor(kBlack);
  // c1g->SetMarkerSize(0.6);
  // c1g->SetLineColor(kBlack);
  // c1g->SetLineWidth(1);
  // c1g->Draw("same P");

  // c1->Update();
  // c1->cd(2);

  // //Clone the histogram which corresponds to the cluster method chosen
  // TH2D *hadctAlignRow = (TH2D*)hatimerow->Clone("hadctAlignRow");

  // vector<Double_t> c1cell_row;
  // vector<Double_t> c1mean_row;
  // vector<Double_t> c1err_row;
    
  // Int_t c1binsX_row = hadctAlignRow->GetNbinsX();

  // util::sliceHisto(hadctAlignRow, c1binsX_row, atime_fwhm, min_epc, c1cell_row, c1mean_row, c1err_row, binXmin, binXmax );
    
  // // Convert vectors to arrays
  // Double_t* c1x_row = &c1cell_row[0];  
  // Double_t* c1y_row = &c1mean_row[0];
  // Double_t* c1ey_row = &c1err_row[0];
  
  // //hadctAlignRow->GetYaxis()->SetRangeUser(-10.,10.);
  // hadctAlignRow->GetYaxis()->SetTitle("ns");
  // hadctAlignRow->GetXaxis()->SetTitle("");
  // hadctAlignRow->Draw("colz");
    
  // //Make graphs with errors for reporting
  // TGraphErrors *c1g_row = new TGraphErrors( c1cell_row.size(), c1x_row, c1y_row, 0, c1ey_row );
    
  // c1g_row->SetTitle(Form("HCal adctime vs row, SBS%d",config));
  // c1g_row->SetMarkerStyle(21); // idx 20 Circles, idx 21 Boxes
  // c1g_row->SetMarkerColor(kBlack);
  // c1g_row->SetMarkerSize(0.6);
  // c1g_row->SetLineColor(kBlack);
  // c1g_row->SetLineWidth(1);
  // c1g_row->Draw("same P");

  // c1->Update();
  // c1->Draw();
  // c1->Write();
 
  fout->Write();
  fout->Close();
  st->Stop();

  cout << endl << "ADCt Analysis and fits written to " << out_path << endl << endl;

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

}
