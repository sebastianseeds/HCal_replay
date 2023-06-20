//sseeds 4.19.23 - script to plot E results on TMultigraph (all kinematics)

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
const Int_t kNcell = 288;
const Int_t kNrows = 24;
const Int_t kNcols = 12;
const Int_t tfitmin = 50;
const Double_t fitl = 0.;
const Double_t fitu = 0.2;

void uniformity_check( Int_t set = 0, const char *experiment = "gmn", Int_t config = 4, Int_t pass = 0 ){
  
  gStyle->SetOptFit();
  gStyle->SetEndErrorSize(0);
  // TCanvas *c1 = new TCanvas("c1","HCal Energy MC Comparison",1600,1200);
  // c1->SetGrid();

  auto mg = new TMultiGraph();

  std::string outdir_path = gSystem->Getenv("OUT_DIR");
  std::string hcal_energy_qr1_path = outdir_path + Form("/hcal_calibrations/pass%d/energy/ecal_class_test_%s_conf%d_qr1_pass%d.root",pass,experiment,config,pass);

  TFile *fsf = TFile::Open(hcal_energy_qr1_path.c_str());
  TH2D *hsfvid = (TH2D*)fsf->Get(Form("hEvID_set%d",set));
  TH2D *hsfvrow = (TH2D*)fsf->Get(Form("hEvX_set%d",set));
  TH2D *hsfvcol = (TH2D*)fsf->Get(Form("hEvY_set%d",set));

  //Make arrays for row tgraphs
  Double_t cellerr[kNrows] = {0.};
  Double_t sfrowcell[kNrows] = {0.};
  Double_t sfrowcval[kNrows] = {0.};
  Double_t sfrowcerr[kNrows] = {0.};
  TH1D *sfrowslice[kNrows];
  
  //Dispersive X counts across rows in HCal
  for( Int_t r=0; r<kNrows; r++ ){
    
    sfrowcell[r] = r+0.5;
    
    sfrowslice[r] = hsfvrow->ProjectionY(Form("sfrowslice_%d",r+1),r+1,r+1);
    Int_t sliceN = sfrowslice[r]->GetEntries();
    if( sliceN<tfitmin ){
      continue;
    }
    
    TF1 *fit1;
    sfrowslice[r]->Fit("gaus","Q","",fitl,fitu);
    fit1 = sfrowslice[r]->GetFunction("gaus");
    sfrowcval[r] = fit1->GetParameter(1);
    sfrowcerr[r] = fit1->GetParameter(2);
    
  }

  //Make arrays for col tgraphs
  Double_t sfcolcell[kNcols] = {0.};
  Double_t sfcolcval[kNcols] = {0.};
  Double_t sfcolcerr[kNcols] = {0.};
  TH1D *sfcolslice[kNcols];
  
  //Dispersive X counts across cols in HCal
  for( Int_t r=0; r<kNcols; r++ ){
    
    sfcolcell[r] = r+0.5;
    
    sfcolslice[r] = hsfvcol->ProjectionY(Form("sfcolslice_%d",r+1),r+1,r+1);
    Int_t sliceN = sfcolslice[r]->GetEntries();
    if( sliceN<tfitmin ){
      continue;
    }
    
    TF1 *fit2;
    sfcolslice[r]->Fit("gaus","Q","",fitl,fitu);
    fit2 = sfcolslice[r]->GetFunction("gaus");
    sfcolcval[r] = fit2->GetParameter(1);
    sfcolcerr[r] = fit2->GetParameter(2);
    
  }

  TCanvas *c1 = new TCanvas("c1","HCal SF vs Row/Col",2200,1200);
  c1->Divide(1,2);
  c1->SetGrid();
  gStyle->SetOptStat(0);
  gStyle->SetPalette(53);

  c1->cd(1);
  hsfvrow->Draw("colz");

  // //Make graphs with errors for reporting
  TGraphErrors *gsfrow = new TGraphErrors( kNrows, sfrowcell, sfrowcval, cellerr, sfrowcerr );

  gsfrow->SetMinimum(0.);
  gsfrow->SetMaximum(1.);
  gsfrow->SetTitle("HCal SF vs Row, SBS4");
  gsfrow->GetXaxis()->SetTitle("Row");
  gsfrow->GetYaxis()->SetTitle("Sampling Fraction (%)");
  gsfrow->SetMarkerStyle(33); // idx 20 Circles, idx 21 Boxes
  gsfrow->SetMarkerColor(kMagenta);
  gsfrow->SetMarkerSize(3);
  gsfrow->SetLineWidth(3);
  gsfrow->SetLineColor(kMagenta);
  for( Int_t r=0; r<kNrows; r++ ){
    cout << sfrowcval[r] << endl;
    Double_t idx = r/kNrows;
    if( sfrowcval[r]==0 ) gsfrow->RemovePoint(idx);
  }

  gsfrow->Draw("P");

  c1->cd(2);
  hsfvcol->Draw("colz");

  // //Make graphs with errors for reporting
  TGraphErrors *gsfcol = new TGraphErrors( kNcols, sfcolcell, sfcolcval, cellerr, sfcolcerr );

  gsfcol->SetMinimum(0.);
  gsfcol->SetMaximum(1.);
  gsfcol->SetTitle("HCal SF vs Col, SBS4");
  gsfcol->GetXaxis()->SetTitle("Col");
  gsfcol->GetYaxis()->SetTitle("Sampling Fraction (%)");
  gsfcol->SetMarkerStyle(33); // idx 20 Circles, idx 21 Boxes
  gsfcol->SetMarkerColor(kMagenta);
  gsfcol->SetMarkerSize(3);
  gsfcol->SetLineWidth(3);
  gsfcol->SetLineColor(kMagenta);
  for( Int_t r=0; r<kNcols; r++ ){
    cout << sfcolcval[r] << endl;
    Double_t idx = r/kNcols;
    if( sfcolcval[r]==0 ) gsfcol->RemovePoint(idx);
  }

  gsfcol->Draw("P");




















  // TF1 *fitmc[nkine];
  // Double_t fitmcp1[nkine];
  // Double_t fitmcp2[nkine];

  // Double_t fitl;
  // Double_t fith;
  // Double_t fits = 0.5;
  // for( Int_t i=0; i<nkine; i++ ){
  //   fitl = fmean[i]-fits;
  //   fith = fmean[i]+fits;

  //   fmc[i] = TFile::Open(Form("../MC/output/MChcalE_sbs%d.root",kIdx[i]));
  //   hmc[i]= (TH1D*)fmc[i]->Get("hHCALe_clus");
  //   hmc[i]->SetLineWidth(0);
    
  //   hmc[i]->Fit("gaus","0","",fitl,fith);
  //   fitmc[i] = hmc[i]->GetFunction("gaus");
  //   fitmcp1[i] = fitmc[i]->GetParameter(1);
  //   fitmcp2[i] = fitmc[i]->GetParameter(2);

  // }

  // auto gr1 = new TGraphErrors(nkine,kDIdxmc,fitmcp1,kXerrmc,fitmcp2);
  // gr1->SetTitle("MC");
  // gr1->SetMarkerColor(kViolet+4);
  // gr1->SetMarkerStyle(33);
  // gr1->SetMarkerSize(2);
  // gr1->SetLineColor(kViolet+4);
  // gr1->SetLineWidth(2);
  // mg->Add(gr1);


  // //Get MC samp frac all kine
  // TFile *f[nkine];
  // TH1D *h[nkine];
  // TF1 *fit[nkine];
  // Double_t fitp1[nkine];
  // Double_t fitp2[nkine];

  // for( Int_t i=0; i<nkine; i++ ){
  //   fitl = fmean[i]-fits;
  //   fith = fmean[i]+fits;

  //   f[i] = TFile::Open(Form("qreplay_sbs%d.root",kIdx[i]));
  //   h[i]= (TH1D*)f[i]->Get("hHCALe");
  //   h[i]->SetLineWidth(0);
    
  //   h[i]->Fit("gaus","0","",fitl,fith);
  //   fit[i] = h[i]->GetFunction("gaus");
  //   fitp1[i] = fit[i]->GetParameter(1);
  //   fitp2[i] = fit[i]->GetParameter(2);

  // }

  // auto gr2 = new TGraphErrors(nkine,kDIdx,fitp1,kXerr,fitp2);
  // gr2->SetTitle("Data (Calibrated)");
  // gr2->SetMarkerColor(kAzure);
  // gr2->SetMarkerStyle(33);
  // gr2->SetMarkerSize(2);
  // gr2->SetLineColor(kAzure);
  // gr2->SetLineWidth(2);
  // mg->Add(gr2);

  // mg->SetTitle("HCal Cluster Energy");
  // //mg->GetXaxis()->SetTitle("HCal E (GeV)");
  // //mg->GetXaxis()->SetTitle("Kinematic");
  // mg->GetYaxis()->SetTitle("HCal Cluster E (GeV)");
  // mg->Draw("AP");

  // mg->GetYaxis()->SetRangeUser(0.,0.8);
  // //mg->GetXaxis()->SetNdivisions(12,3,15,kTRUE);
  // mg->GetXaxis()->SetNdivisions(6,"I");
  
  // //cout << mg->GetNBins() << endl;

  // for( Int_t i=0; i<nkine; i++ ){
  //   mg->GetXaxis()->SetBinLabel( i*17.5+7, Form("SBS-%d",kIdx[i]));
  // }
  // mg->GetXaxis()->LabelsOption("h");
  // mg->GetXaxis()->SetLabelSize(0.05);

  // c1->Modified();

  // c1->BuildLegend();


}
