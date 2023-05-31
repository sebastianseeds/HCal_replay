//sseeds 5.9.23 - Script to plot before/after timing alignment comparisons for both TDC and ADCt for a given experiment and configuration

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
#include "../../include/hcal.h"

const Double_t Y_range_factor = 20; //plot +/- about average value for comparison. Adjust depending on variance across HCal ID

//Main. Pass exp{gmn,gen,genrp,gep,etc}, config{<int>configuration number}, type{tdc,adtc}, setN{<int> calibration subset number}
void plot_talign_comparison( const char *exp = "gmn", Int_t config = 4, const char *type = "tdc", Int_t setN = 0, Int_t pass=0 ){
  //build the canvas
  gStyle->SetOptFit();
  gStyle->SetEndErrorSize(0);
  TCanvas *c1 = new TCanvas("c1",Form("HCal %s Alignment Comparison, %s Config %d",type,exp,config),1600,1200);
  c1->SetGrid();

  //declare multigraph to store both output tgrapherrors
  auto mg = new TMultiGraph();
  Double_t x,y;
  Double_t avg = 0;

  std::string outdir_path = gSystem->Getenv("OUT_DIR");
  std::string align_path = outdir_path + Form("/hcal_calibrations/pass%d/timing/%salign_%s_conf%d_qr0.root",pass,type,exp,config);
  std::string align_qr_path = outdir_path + Form("/hcal_calibrations/pass%d/timing/%salign_%s_conf%d_qr1.root",pass,type,exp,config);

  //open files and get tgraphs
  // TFile *f = TFile::Open(Form("%salign_%s_conf%d_qr0.root",type,exp,config));
  // TFile *fqr = TFile::Open(Form("%salign_%s_conf%d_qr1.root",type,exp,config));
  TFile *f = TFile::Open(align_path.c_str());
  TFile *fqr = TFile::Open(align_qr_path.c_str());
  TGraphErrors *g = (TGraphErrors*)f->Get(Form("g%s_c_s%d",type,setN));
  for( Int_t i=0; i<hcal::maxHCalChan; i++ ){ //remove points where not enough data to fit (y=0) and get Y average
    x=-1;
    y=-1;
    g->GetPoint(i,x,y);
    if( y==0 ) g->RemovePoint( i );
    avg += y;

    if( x==hcal::maxHCalChan-1 ) break; //avoid adding points to mg
  }
  avg /= hcal::maxHCalChan;
  //options for first tgraph errors
  g->SetTitle("Unaligned");
  g->SetMarkerColor(kGray);
  g->SetMarkerStyle(33);
  g->SetMarkerSize(1);
  g->SetLineColor(kGray);
  g->SetLineWidth(1);
  mg->Add(g);

  TGraphErrors *gqr = (TGraphErrors*)fqr->Get(Form("g%s_c_s%d",type,setN));
  for( Int_t i=0; i<hcal::maxHCalChan; i++ ){ //remove points and offset values by 0.5 in X for clarity
    x=-1;
    y=-1;
    gqr->GetPoint(i,x,y);
    if( y!=0 ){
      gqr->SetPoint( i,x+0.5,y );
    }else{
      gqr->RemovePoint( i );
    }

    if( x==hcal::maxHCalChan-1 ) break; //avoid adding points to mg
  }
  //options for second tgrapherrors
  gqr->SetTitle("Aligned");
  gqr->SetMarkerColor(kBlack);
  gqr->SetMarkerStyle(34);
  gqr->SetMarkerSize(1);
  gqr->SetLineColor(kBlack);
  gqr->SetLineWidth(1);
  mg->Add(gqr);

  //options for mg
  mg->SetTitle(Form("HCal %s Timing Alignment Comparison, %s Config %d",type,exp,config));
  mg->GetXaxis()->SetTitle("HCal Channel");
  mg->GetYaxis()->SetTitle(Form("HCal %s (ns)",type));
  mg->Draw("AP");

  mg->GetYaxis()->SetRangeUser(avg-Y_range_factor,avg+Y_range_factor);

  c1->Modified();

  c1->BuildLegend(); //must happen after mg is drawn

}
