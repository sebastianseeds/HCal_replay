//sseeds 5.23.23 - Script to generate supplemental quality plots for timing alignments and energy calibrations

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

const Double_t hcalE_sigma_factor = 0.05; //Rough correspondance between energy spectra width and beam energy (beamE * hcalE_sigma_factor = hcalE distribution sigma. Used for fit limits only.
const Double_t Y_range_factor = 20; //plot +/- about average value for comparison. Adjust depending on variance across HCal ID
const Double_t SF_hard_llim = 0.02; //Hard lower limit on hcal SF for fitting
const Double_t SF_hard_ulim = 0.30; //Hard lower limit on hcal SF for fitting
const Double_t E_hard_llim = 0.02; //Hard lower limit on hcal E for fitting
const Double_t E_hard_ulim = 0.30; //Hard lower limit on hcal E for fitting

//Gets the maximum bin on a user passed range
Int_t get_max_bin_on_range( TH1D *h, Double_t xlow, Double_t xhigh ){

  Int_t bin1 = h->FindBin(xlow);  // convert xlow to bin number
  Int_t bin2 = h->FindBin(xhigh); // convert xhigh to bin number

  // Initial maximum and its position
  Double_t maxVal = h->GetBinContent(bin1);
  Int_t maxBin = bin1;

  // Loop over bins in the range to find the maximum
  for (Int_t i = bin1; i <= bin2; i++) {
    Double_t binContent = h->GetBinContent(i);
    if (binContent > maxVal) {
      maxVal = binContent;
      maxBin = i;
    }
  }

  return maxBin;
}

//Generates multigraph with timing alignment comparison
void generate_mg( std::string qr0path, std::string qr1path, std::string type, TMultiGraph *mg, Int_t set_N, Double_t avg ){

  //declare necessary variables
  Double_t x,y;
  avg = 0.;
  
  //open files and get tgraphs
  TFile *f = TFile::Open(qr0path.c_str());
  TFile *fqr = TFile::Open(qr1path.c_str()); 
  
  TGraphErrors *g = (TGraphErrors*)f->Get(Form("g%s_c_s%d",type.c_str(),set_N));
  for( Int_t i=0; i<hcal::maxHCalChan; i++ ){ //remove points where not enough data to fit (y=0) and get Y average
    x=-1;
    y=-1;
    g->GetPoint(i,x,y);
    if( y==0 ) g->RemovePoint( i );
    avg += y;

    if( x==hcal::maxHCalChan-1 ) break; //avoid adding points to mg
  }
  avg /= hcal::maxHCalChan;
  g->SetTitle("Unaligned");
  g->SetMarkerColor(kGray);
  g->SetMarkerStyle(33);
  g->SetMarkerSize(1);
  g->SetLineColor(kGray);
  g->SetLineWidth(1);
  mg->Add(g);

  TGraphErrors *gqr = (TGraphErrors*)fqr->Get(Form("g%s_c_s%d",type.c_str(),set_N));
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

}

//Main. Script to generate all necessary plots for calibration quality checks. See overleaf for quality check template
void quality_check( const char *experiment = "gmn", Int_t config = 11, Int_t pass = 1 ){
  
  // set style
  gStyle->SetOptFit(0);
  gStyle->SetEndErrorSize(0);
  gStyle->SetOptStat(0);

  // Get beam energy by configuration
  SBSconfig config_parameters(experiment,config);    
  Double_t beamE = config_parameters.GetEbeam();

  //////////
  //ADCt alignment plots

  // Get ADCt alignment comparison canvas
  TCanvas *c1 = new TCanvas("c1","adct comparison",1600,1200);

  //get paths to calibration output files and type
  std::string outdir_path = gSystem->Getenv("OUT_DIR");
  std::string adct_type = "adct";
  std::string adct_param = "sbs.hcal.adc.timeoffset";
  std::string adct_align_qr0_path = outdir_path + Form("/hcal_calibrations/pass%d/timing/%salign_class_%s_conf%d_qr0_pass%d.root",pass,adct_type.c_str(),experiment,config,pass);
  std::string adct_align_qr1_path = outdir_path + Form("/hcal_calibrations/pass%d/timing/%salign_class_%s_conf%d_qr1_pass%d.root",pass,adct_type.c_str(),experiment,config,pass);
  std::string adct_param_path = Form("../timing/parameters/adctoffsets_class_%s_conf%d_pass%d.txt",experiment,config,pass);

  //determine how many calibration sets exist
  Int_t adct_align_Nsets = util::countSets( adct_param_path, adct_param );
  
  // divide the canvas
  c1->Divide(1,adct_align_Nsets);

  //Generate adct multigraphs by calibration set
  for( Int_t set=0; set<adct_align_Nsets; set++ ){
    
    c1->cd(set+1);

    TMultiGraph *adctmg = new TMultiGraph();
    Double_t adct_avg;

    generate_mg( adct_align_qr0_path, adct_align_qr1_path, adct_type, adctmg, set, adct_avg );

    //options for mg
    adctmg->SetTitle(Form("HCal %s Timing Alignment Comparison, %s Config %d, Calibration Set %d",adct_type.c_str(),experiment,config,set));
    adctmg->GetXaxis()->SetTitle("HCal Channel");
    adctmg->GetYaxis()->SetTitle(Form("HCal %s (ns)",adct_type.c_str()));
    adctmg->Draw("AP");
    
    adctmg->GetYaxis()->SetRangeUser(adct_avg-Y_range_factor,adct_avg+Y_range_factor);
    
    c1->Modified();
    
    c1->BuildLegend(); //must happen after mg is drawn

  }

  //////////
  //TDC plots

  // Get TDC alignment canvas
  TCanvas *c2 = new TCanvas("c2","tdc alignment",1600,1200);

  //get paths to calibration output files and type
  std::string tdc_type = "tdc";
  std::string tdc_param = "sbs.hcal.tdc.offset";
  std::string tdc_align_qr0_path = outdir_path + Form("/hcal_calibrations/pass%d/timing/%salign_class_%s_conf%d_qr0_pass%d.root",pass,tdc_type.c_str(),experiment,config,pass);
  std::string tdc_align_qr1_path = outdir_path + Form("/hcal_calibrations/pass%d/timing/%salign_class_%s_conf%d_qr1_pass%d.root",pass,tdc_type.c_str(),experiment,config,pass);
  std::string tdc_param_path = Form("../timing/parameters/tdcoffsets_class_%s_conf%d_pass%d.txt",experiment,config,pass);

  //determine how many calibration sets exist
  Int_t tdc_align_Nsets = util::countSets( tdc_param_path, tdc_param );
  
  // divide the canvas
  c2->Divide(1,tdc_align_Nsets); 

  //Generate tdc multigraphs by calibration set
  for( Int_t set=0; set<tdc_align_Nsets; set++ ){
    
    c2->cd(set+1);

    TMultiGraph *tdcmg = new TMultiGraph();
    Double_t tdc_avg;

    generate_mg( tdc_align_qr0_path, tdc_align_qr1_path, tdc_type, tdcmg, set, tdc_avg );

    //options for mg
    tdcmg->SetTitle(Form("HCal %s Timing Alignment Comparison, %s Config %d, Calibration Set %d",tdc_type.c_str(),experiment,config,set));
    tdcmg->GetXaxis()->SetTitle("HCal Channel");
    tdcmg->GetYaxis()->SetTitle(Form("HCal %s (ns)",tdc_type.c_str()));
    tdcmg->Draw("AP");
    
    tdcmg->GetYaxis()->SetRangeUser(tdc_avg-Y_range_factor,tdc_avg+Y_range_factor);
    
    c2->Modified();
    
    c2->BuildLegend(); //must happen after mg is drawn

  }

  //////////
  //Energy Plots

  // Get sampling fraction canvas
  TCanvas *c3 = new TCanvas("c3","sampling fraction",1600,1200);
  
  //get paths to calibration output files and type
  //std::string tdc_type = "tdc";
  std::string adc_gain_param = "sbs.hcal.adc.gain";
  std::string hcal_energy_qr0_path = outdir_path + Form("/hcal_calibrations/pass%d/energy/ecal_class_%s_conf%d_qr0_pass%d.root",pass,experiment,config,pass);
  std::string hcal_energy_qr1_path = outdir_path + Form("/hcal_calibrations/pass%d/energy/ecal_class_%s_conf%d_qr1_pass%d.root",pass,experiment,config,pass);
  std::string adc_gain_path = Form("../energy/parameters/adcgaincoeff_%s_conf%d_pass%d.txt",experiment,config,pass);

  //determine how many calibration sets exist
  Int_t hcal_energy_Nsets = util::countSets( adc_gain_path, adc_gain_param );

  // divide the canvas
  c3->Divide(1,hcal_energy_Nsets);

  //Generate adct multigraphs by calibration set
  Double_t obs_fwhm = beamE*hcalE_sigma_factor;
  for( Int_t set=0; set<hcal_energy_Nsets; set++ ){
    
    c3->cd(set+1);

    //Get pre-ecal samp frac
    TFile *f1 = TFile::Open(hcal_energy_qr0_path.c_str());
    TH1D *h1= (TH1D*)f1->Get(Form("hSF_set%d",set));
    Double_t x1min = h1->GetXaxis()->GetXmin();
    Double_t x1max = h1->GetXaxis()->GetXmax();
    Int_t x1totalbins = h1->GetNbinsX();
    Int_t x1maxX = (SF_hard_ulim-x1min)*x1totalbins/(x1max-x1min);
    h1->GetXaxis()->SetRange(0,x1maxX);
    h1->SetLineWidth(3);
    h1->SetLineColor(kRed);
    h1->Draw("hist");

    //Get reasonable limits on fit
    Int_t x1binmax = h1->GetMaximumBin();
    Double_t x1binmaxX = x1min+x1binmax*(x1max-x1min)/x1totalbins;
    Double_t x1llim = x1binmaxX-obs_fwhm/2;
    Double_t x1ulim = x1binmaxX+obs_fwhm/2;
    if( x1llim < SF_hard_llim ){ //if the lower limit is below minimum, recalculate
      x1binmax = get_max_bin_on_range( h1, SF_hard_llim, x1max );
      x1binmaxX = x1min+x1binmax*(x1max-x1min)/x1totalbins;
      x1llim = SF_hard_llim;
      x1ulim = x1binmaxX + ( x1binmaxX - x1llim );
    }

    //Get post-ecal samp frac
    TFile *f2 = TFile::Open(hcal_energy_qr1_path.c_str());
    TH1D *h2 = (TH1D*)f2->Get(Form("hSF_set%d",set));
    Double_t x2min = h2->GetXaxis()->GetXmin();
    Double_t x2max = h2->GetXaxis()->GetXmax();
    Int_t x2totalbins = h2->GetNbinsX();
    h2->SetLineWidth(2);
    h2->SetLineColor(kBlack);
    h2->Draw("hist same");

    //Get reasonable limits on fit
    Int_t x2binmax = h2->GetMaximumBin();
    Double_t x2binmaxX = x2min+x2binmax*(x2max-x2min)/x2totalbins;
    Double_t x2llim = x2binmaxX-obs_fwhm/2;
    Double_t x2ulim = x2binmaxX+obs_fwhm/2;
    if( x2llim < SF_hard_llim ){ //if the lower limit is below minimum, recalculate
      x2binmax = get_max_bin_on_range( h2, SF_hard_llim, x2max );
      x2binmaxX = x2min+x2binmax*(x2max-x2min)/x2totalbins;
      x2llim = SF_hard_llim;
      x2ulim = x2binmaxX + ( x2binmaxX - x2llim );
    }

    //Fit the pre-cal SF distribution with gaussian
    TF1 *fit1;
    h1->Fit("gaus","","",x1llim,x1ulim);
    fit1 = h1->GetFunction("gaus");
    Double_t fit1m = fit1->GetParameter(1);

    //Fit the pre-cal SF distribution with gaussian
    TF1 *fit2;
    h2->Fit("gaus","","",x2llim,x2ulim);
    fit2 = h2->GetFunction("gaus");
    Double_t fit2m = fit2->GetParameter(1);

    //Add a legend
    auto legend = new TLegend(0.43,0.7,0.89,0.89);
    legend->SetTextSize(0.03);
    legend->SetHeader(Form("HCal eCal Samp Frac SBS%d, Cal Set%d",config,set));
    legend->AddEntry(h1,Form("Pre-Cal, mean:%f",fit1m),"l");
    legend->AddEntry(h2,Form("Post-Cal, mean:%f",fit2m),"l");
    legend->Draw();

  }

  // Get energy spectra canvas
  TCanvas *c4 = new TCanvas("c4","energy",1600,1200);

  // divide the canvas
  c4->Divide(1,hcal_energy_Nsets);

  //Generate adct multigraphs by calibration set
  for( Int_t set=0; set<hcal_energy_Nsets; set++ ){
    
    c4->cd(set+1);

    //Get pre-ecal samp frac
    TFile *f1 = TFile::Open(hcal_energy_qr0_path.c_str());
    TH1D *h1= (TH1D*)f1->Get(Form("hE_set%d",set));    
    Double_t x1min = h1->GetXaxis()->GetXmin();
    Double_t x1max = h1->GetXaxis()->GetXmax();
    Int_t x1totalbins = h1->GetNbinsX();
    Int_t x1maxX = (E_hard_ulim-x1min)*x1totalbins/(x1max-x1min);
    h1->GetXaxis()->SetRange(0,x1maxX);
    h1->SetLineWidth(3);
    h1->SetLineColor(kRed);
    h1->Draw("hist");

    //Get reasonable limits on fit
    Int_t x1binmax = h1->GetMaximumBin();
    Double_t x1binmaxX = x1min+x1binmax*(x1max-x1min)/x1totalbins;
    Double_t x1llim = x1binmaxX-obs_fwhm/2;
    Double_t x1ulim = x1binmaxX+obs_fwhm/2;
    if( x1llim < E_hard_llim ){ //if the lower limit is below minimum, recalculate
      x1binmax = get_max_bin_on_range( h1, E_hard_llim, x1max );
      x1binmaxX = x1min+x1binmax*(x1max-x1min)/x1totalbins;
      x1llim = E_hard_llim;
      x1ulim = x1binmaxX + ( x1binmaxX - x1llim );
    }

    //Get post-ecal samp frac
    TFile *f2 = TFile::Open(hcal_energy_qr1_path.c_str());
    TH1D *h2 = (TH1D*)f2->Get(Form("hE_set%d",set));
    Double_t x2min = h2->GetXaxis()->GetXmin();
    Double_t x2max = h2->GetXaxis()->GetXmax();
    Int_t x2totalbins = h2->GetNbinsX();
    h2->SetLineWidth(2);
    h2->SetLineColor(kBlack);
    h2->Draw("hist same");

    //Get reasonable limits on fit
    Int_t x2binmax = h2->GetMaximumBin();
    Double_t x2binmaxX = x2min+x2binmax*(x2max-x2min)/x2totalbins;
    Double_t x2llim = x2binmaxX-obs_fwhm/2;
    Double_t x2ulim = x2binmaxX+obs_fwhm/2;
    if( x2llim < E_hard_llim ){ //if the lower limit is below minimum, recalculate
      x2binmax = get_max_bin_on_range( h2, E_hard_llim, x2max );
      x2binmaxX = x2min+x2binmax*(x2max-x2min)/x2totalbins;
      x2llim = E_hard_llim;
      x2ulim = x2binmaxX + ( x2binmaxX - x2llim );
    }

    //Fit the pre-cal E distribution with gaussian
    TF1 *fit1;
    h1->Fit("gaus","","",x1llim,x1ulim);
    fit1 = h1->GetFunction("gaus");
    Double_t fit1m = fit1->GetParameter(1);

    //Fit the pre-cal E distribution with gaussian
    TF1 *fit2;
    h2->Fit("gaus","","",x2llim,x2ulim);
    fit2 = h2->GetFunction("gaus");
    Double_t fit2m = fit2->GetParameter(1);

    //Add a legend
    auto legend = new TLegend(0.43,0.7,0.89,0.89);
    legend->SetTextSize(0.03);
    legend->SetHeader(Form("HCal eCal Energy Spectrum SBS%d, Cal Set%d",config,set));
    legend->AddEntry(h1,Form("Pre-Cal, mean:%f",fit1m),"l");
    legend->AddEntry(h2,Form("Post-Cal, mean:%f",fit2m),"l");
    legend->Draw();

  }

  // Get monte carlo sampling fraction canvas
  TCanvas *c5 = new TCanvas("c5","MC sampling fraction",1600,1200);
  c5->cd();

  //get path to mc output files
  std::string mc_path = outdir_path + Form("/hcal_calibrations/MC/hcalE_mc_class_%s_conf%d.root",experiment,config);

  TFile *fSFmc = TFile::Open(mc_path.c_str());
  TH1D *hSFmc = (TH1D*)fSFmc->Get("hSampFrac_clus");

  hSFmc->SetLineWidth(3);
  hSFmc->Draw("hist");

  //get reasonable limits on fit
  Double_t SFmcmin = hSFmc->GetXaxis()->GetXmin();
  Double_t SFmcmax = hSFmc->GetXaxis()->GetXmax();
  Int_t SFmctotalbins = hSFmc->GetNbinsX();
  Int_t SFmcbinmax = hSFmc->GetMaximumBin();
  Double_t SFmcbinmaxX = SFmcmin+SFmcbinmax*(SFmcmax-SFmcmin)/SFmctotalbins;
  Double_t SFmcllim = SFmcbinmaxX-obs_fwhm/2;
  Double_t SFmculim = SFmcbinmaxX+obs_fwhm/2;
  if( SFmcllim < E_hard_llim ){ //if the lower limit is below minimum, recalculate
    SFmcbinmax = get_max_bin_on_range( hSFmc, SF_hard_llim, SFmcmax );
    SFmcbinmaxX = SFmcmin+SFmcbinmax*(SFmcmax-SFmcmin)/SFmctotalbins;
    SFmcllim = SF_hard_llim;
    SFmculim = SFmcbinmaxX + ( SFmcbinmaxX - SFmcllim );
  }


  hSFmc->Fit("gaus","","",SFmcllim,SFmculim);
  TF1 *fitSFmc = hSFmc->GetFunction("gaus");
  Double_t fitSFmcp1 = fitSFmc->GetParameter(1)*100.;
  
  hSFmc->SetTitle(Form("MC Samp Frac SBS%d, peak %0.2f%%",config,fitSFmcp1));

  // Get monte carlo energy spectrum canvas
  TCanvas *c6 = new TCanvas("c6","MC energy spectrum",1600,1200);
  c6->cd();

  TFile *fEmc = TFile::Open(mc_path.c_str());
  TH1D *hEmc = (TH1D*)fEmc->Get("hHCALe_clus");

  hEmc->SetLineWidth(3);
  hEmc->Draw("hist");

  //get reasonable limits on fit
  Double_t Emcmin = hEmc->GetXaxis()->GetXmin();
  Double_t Emcmax = hEmc->GetXaxis()->GetXmax();
  Int_t Emctotalbins = hEmc->GetNbinsX();
  Int_t Emcbinmax = hEmc->GetMaximumBin();
  Double_t EmcbinmaxX = Emcmin+Emcbinmax*(Emcmax-Emcmin)/Emctotalbins;
  Double_t Emcllim = EmcbinmaxX-obs_fwhm/2;
  Double_t Emculim = EmcbinmaxX+obs_fwhm/2;
  if( Emcllim < E_hard_llim ){ //if the lower limit is below minimum, recalculate
    Emcbinmax = get_max_bin_on_range( hEmc, E_hard_llim, Emcmax );
    EmcbinmaxX = Emcmin+Emcbinmax*(Emcmax-Emcmin)/Emctotalbins;
    Emcllim = E_hard_llim;
    Emculim = EmcbinmaxX + ( EmcbinmaxX - Emcllim );
  }

  hEmc->Fit("gaus","0","",Emcllim,Emculim);
  TF1 *fitEmc = hEmc->GetFunction("gaus");
  Double_t fitEmcp1 = fitEmc->GetParameter(1);
  
  hEmc->SetTitle(Form("MC Energy Spectrum SBS%d, peak %0.2f",config,fitEmcp1));

}
