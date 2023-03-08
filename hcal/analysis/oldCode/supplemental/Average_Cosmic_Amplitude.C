#include <TH2F.h>
#include <TChain.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <TSystem.h>
#include "hcal.h"
#include <vector>
#include <TStopwatch.h>
#include <TMath.h>
#include <TF1.h>
#include<TStyle.h>
#include<TList.h>
#include<TChain.h>
#include <TDatime.h>
using namespace std;

Int_t gCurrentEntry = 0;
Int_t DISP_MIN_SAMPLE = 0;
Int_t DISP_MAX_SAMPLE = 30;
const Int_t DISP_FADC_SAMPLES = (DISP_MAX_SAMPLE-DISP_MIN_SAMPLE);
Int_t bins = 100;
Int_t r,c,n,idx;
Int_t nevt;
Double_t cmu_hv = 1500, jlab_hv = 1600, cmu_exp = 10.5, jlab_exp = 8.; //Note voltages are actually negative but it's easier to plot if positive.
Double_t target_RAU = 61.425;
Double_t adc_fire_avg = 10.;                  //Minimum fADC average RAU pedestal standard deviations above which fADC 'fired'.
Double_t adc_fire_int = 3.;                  //Minimum fADC integral pedestal standard deviations above which fADC 'fired'.
Double_t adc_int_fire = 100.;
const Int_t kNrows = 24;
const Int_t kNcols = 12;
const Int_t channels = kNrows * kNcols;
Int_t limit_evts = 0;                     //0-> analyze all events. 1-> only analyze events up until max_evts.
Int_t loop_max = 0;                       //Dummy variable set equal to max_evts or nevt for the loop.
Int_t max_evts = 50000;                   //Maximum number of events to analyze if limit_evts = 1.
Int_t use_int = 0;                        //Create plots with the fADC integral too (normally just look at fADC maz RAU). No pedestal subtraction currently.
Double_t f1_res = 0.112;//0.112
Double_t fadc_res = 4.;
TChain *T = 0;
std::string user_input;
Int_t skip = 3, skip1 = 1;                            //Gives number of lines to skip at top of data file. 
Int_t nlines = 0, nlines1 = 0;                          //Counts number of lines in the data file. 
Int_t ncols, ncols1;                               //Set how many columns of data we have in the data file.
char str[1000], str1[1000];                           //Variable to read lines of the data file.
Float_t mod[channels], mod1[channels], hv[channels], ped_avg[channels], ped_int[channels], ped_avg_std_dev[channels], ped_int_std_dev[channels], threshold[channels], avg_amp_tdc[channels], avg_amp_adc[channels], avg_amp_adc_tdc[channels], avg_amp_adc_int[channels], target_voltage_tdc[channels], target_voltage_adc[channels], target_voltage_adc_tdc[channels], target_voltage_adc_int[channels];
Int_t ntdc_hits[channels], nadc_hits[channels], nadc_tdc_hits[channels], ntdc_hits_vert[channels], nadc_int_hits[channels], nadc_hits_vert[channels], nadc_tdc_hits_vert[channels], nadc_int_hits_vert[channels];
Float_t mod_temp, mod1_temp, hv_temp, ped_avg_temp, ped_int_temp, ped_avg_std_dev_temp, ped_int_std_dev_temp, threshold_temp;

//Create Gaussian to fit the timing resolution.
Double_t fit_gaus(Double_t *X,Double_t *par) 
{
  Double_t fitval = par[0]*TMath::Exp(-0.5*pow(((X[0]-par[1])/par[2]),2));
  return fitval;
}

//Create skewed normal distribution to fit the timing resolution.
Double_t fit_skew(Double_t *X,Double_t *par) 
{
  //Double_t fitval = (2./par[1])*(1./pow(2.*TMath::Pi(),0.5))*TMath::Exp(-0.5*pow((X[0]-par[1])/par[0],2.)) * 0.5*(1+TMath::Erf(   par[2]*((X[0]-par[1])/par[0])/pow(2.,0.5)   ));
  //Double_t fitval = (2./par[1])*TMath::Exp(-0.5*pow((X[0]-par[1])/par[0],2.)) * 0.5*(1+TMath::Erf(   par[2]*((X[0]-par[1])/par[0])   ));
  Double_t fitval = par[0]*TMath::Exp(-0.5*pow((X[0]-par[1])/par[2],2.)) * 0.5*(1+TMath::Erf(   par[3]*((X[0]-par[1])/par[2])   ));

  return fitval;
}

TH1F *histos[kNrows][kNcols];

TH1F* MakeHisto(Int_t row, Int_t col, Int_t bins)
{
  TH1F *h = new TH1F(TString::Format("h%02d%02d",row,col),
		     TString::Format("%d-%d",row+1,col+1),bins,DISP_MIN_SAMPLE,DISP_MAX_SAMPLE);
  h->SetStats(0);
  h->SetLineWidth(2);
  return h;
}

//Create Landau to fit the fADC pulses.
Double_t fit_landau(Double_t *X, Double_t *par) 
{
  Double_t landau = par[0]*TMath::Landau(X[0],par[1],par[2])+par[3];
  return landau;
}

//Create a Weibull distribution to try to fit the average cosmic amplitudes.
Double_t fit_weibull(Double_t *X, Double_t *par) 
{
  Double_t weibull = par[2]*(par[0]/par[1])*pow((X[0]-par[3])/par[1],par[0]-1)*exp(-pow((X[0]-par[3])/par[1],par[0]));
  return weibull;
}

//Create a lognormal distribution to try to fit the average cosmic amplitudes.
/*
Double_t fit_lognorm(Double_t *X, Double_t *par) 
{
  //Double_t lognorl = (   1/( par[0]*X[0]*pow(2*TMath::Pi(),0.5) ) * exp(-)   ) ;
  return lognorm;
}
*/

void Average_Cosmic_Amplitude(Int_t run = 1221)
//Int_t Test_Display(Int_t run = 290)
{
  //Define a new stopwatch.
  TStopwatch *st=new TStopwatch();
  st->Start(kTRUE);

  gROOT->SetBatch(kTRUE);//Always run in batch mode otherwise 100s of plot will appear.

  FILE *fp;
  fp = fopen("/home/daq/test_fadc/Pedestals_Updated.txt","r");
 //Read in data.
  while (1) 
    {
      //Skips the first skip lines of the file. 
      if (nlines < skip)
	{
	  fgets(str,1000,fp);
	  nlines++;
	}
      //Reads the two columns of data into x and y.
      else
	{
	  //Read in the number of columns of data in your data file. 
	  ncols = fscanf(fp,"%f %f %f %f %f", &mod_temp, &ped_int_temp, &ped_int_std_dev_temp, &ped_avg_temp, &ped_avg_std_dev_temp);
	  
	  if (ncols < 0) break;    
	  
	  mod[nlines-skip] = mod_temp;
	  ped_avg[nlines-skip] = ped_avg_temp;
	  ped_int[nlines-skip] = ped_int_temp;
	  ped_avg_std_dev[nlines-skip] = ped_avg_std_dev_temp;
	  ped_int_std_dev[nlines-skip] = ped_int_std_dev_temp;
	  //Threshold is defined as 1/4 of the average peak height minus the pedestal during a real event. Then add back in the pedestal value since the Landau fit I use doesn't subtract the pedestal.
	  //threshold[nlines-skip] = (threshold_temp-pedestal[nlines-skip])/4.0+pedestal[nlines-skip];

	  nlines++;
	}
    }

  FILE *fp1;
  fp1 = fopen(Form("/home/daq/test_fadc/Voltage_Scans/HV_Settings/HV_Run_%d.txt",run),"r");
 //Read in data.
  while (1) 
    {
      //Skips the first skip lines of the file. 
      if (nlines1 < skip1)
	{
	  fgets(str1,1000,fp1);
	  nlines1++;
	}
      //Reads the two columns of data into x and y.
      else
	{
	  //Read in the number of columns of data in your data file. 
	  ncols1 = fscanf(fp1,"%f %f", &mod1_temp, &hv_temp);
	  
	  if (ncols1 < 0) break;    
	  
	  mod1[nlines1-skip1] = mod1_temp;
	  hv[nlines1-skip1] = hv_temp;

	  nlines1++;
	}
    }

  if(!T) 
    { 
      T = new TChain("T");
 
      //============  Reading the Rootfile =======================//
      
      const TString rootfilePath = "/home/daq/test_fadc/rootfiles/";
      std::ostringstream str;
      str << rootfilePath<<"fadc_f1tdc_"<<run;
      TString basename = str.str().c_str();
      TString rootfile = basename + ".root";
      cout<<basename<<endl;
      
      Long_t split=0;
      char* file = 0;
      
      //====adding splits rootfiles =======================//
      
      Long_t u=0;
      while ( !gSystem->AccessPathName(rootfile.Data()) ) 
	{
	  T->Add(rootfile.Data());
	  cout << "ROOT file " << rootfile << " added to TChain." << endl;
	  u++;
	  rootfile = basename + "_" + u + ".root";
	}
      
      if(!T->GetEntries())
	{
	  cerr<< "No root file was found" << endl;
	  return;
	}
      //==finish adding splits rootfiles=====================//

      T->SetBranchStatus("*",0);
      T->SetBranchStatus("sbs.hcal.*",1);
      T->SetBranchAddress("sbs.hcal.nsamps",hcalt::nsamps);       //Number of samples for given row-col
      T->SetBranchAddress("sbs.hcal.a",hcalt::a);                 //Raw ADC amplitudes
      T->SetBranchAddress("sbs.hcal.tdc",hcalt::tdc);             //Raw TDC value
      //T->SetBranchAddress("sbs.hcal.ledbit",&hcalt::ledbit);
      //T->SetBranchAddress("sbs.hcal.ledcount",&hcalt::ledcount);
      T->SetBranchAddress("sbs.hcal.samps",hcalt::samps);         //RAW ADC samples
      T->SetBranchAddress("sbs.hcal.samps_idx",hcalt::samps_idx); //Index in samples vector for given row-col module 
      T->SetBranchAddress("sbs.hcal.row",hcalt::row);             //Row for block in data vectors
      T->SetBranchAddress("sbs.hcal.col",hcalt::col);             //Col for block in data vectors
      T->SetBranchStatus("Ndata.sbs.hcal.row",1);
      T->SetBranchAddress("Ndata.sbs.hcal.row",&hcalt::ndata);    //???
      std::cerr << "Opened up tree with nentries=" << T->GetEntries() << std::endl;
      }

  for(Int_t r = 0; r < kNrows; r++) 
    {
      for(Int_t c = 0; c < kNcols; c++) 
	{
	  histos[r][c] = MakeHisto(r,c,DISP_FADC_SAMPLES);
	}
    }

  T->GetEntry(gCurrentEntry);
  cout << "Displaying event " << gCurrentEntry << endl;
  cout<<"ndata = "<<hcalt::ndata<<endl;                                  //Number of channels.
  cout<<"row[] = "<<hcalt::row[16]<<"   col[] = "<<hcalt::col[16]<<endl;
  cout<<"a[] = "<<hcalt::a[1]<<endl;
  cout<<"tdc[] = "<<hcalt::tdc[1]<<endl;

  n = hcalt::nsamps[0];
  DISP_MAX_SAMPLE = hcalt::nsamps[0];
  cout<<"n = "<<n<<endl;

  idx = hcalt::samps_idx[0];
  cout<<"idx = "<<idx<<endl;

  nevt = T->GetEntries();
  cout<<"Number of events = "<<nevt<<endl;

  Float_t peak[kNrows][kNcols];
  Double_t adc[kNrows][kNcols];
  Double_t tdc[kNrows][kNcols];

  TF1 **func_landau = new TF1*[channels];
  for(Int_t i = 0; i<hcalt::ndata; i++)
    {
      func_landau[i] = new TF1("func_landau",fit_landau, DISP_MIN_SAMPLE, DISP_MAX_SAMPLE, 4);
    }

  Double_t max_min = 0, max_max = 1000, max_bins = 201, int_min = 5000, int_max = 20000, int_bins = 1500;

  //Store the max value of the fADC pulse.
  TH1F **hmax_tdc = new TH1F*[channels];
  for(Int_t i = 0; i<channels; i++)
    {
      hmax_tdc[i] = new TH1F(Form("hmax%d_tdc",i),Form("Maximum fADC Values for Module %d (Vertical TDC)",i),max_bins,max_min,max_max);
      hmax_tdc[i]->GetXaxis()->SetTitle("Maximum fADC Value");
      hmax_tdc[i]->GetXaxis()->CenterTitle();
      hmax_tdc[i]->GetYaxis()->SetTitle("Occurrences");
      hmax_tdc[i]->GetYaxis()->CenterTitle();
    }

  TH1F **hmax_adc = new TH1F*[channels];
  for(Int_t i = 0; i<channels; i++)
    {
      hmax_adc[i] = new TH1F(Form("hmax%d_adc",i),Form("Maximum fADC Values for Module %d (Vertical ADC)",i),max_bins,max_min,max_max);
      hmax_adc[i]->GetXaxis()->SetTitle("Maximum fADC Value");
      hmax_adc[i]->GetXaxis()->CenterTitle();
      hmax_adc[i]->GetYaxis()->SetTitle("Occurrences");
      hmax_adc[i]->GetYaxis()->CenterTitle();
    }

  TH1F **hadc_int = new TH1F*[channels];
  for(Int_t i = 0; i<channels; i++)
    {
      hadc_int[i] = new TH1F(Form("hadc_int%d",i),Form("Integral of Average fADC Values for Module %d (Vertical ADC)",i),int_bins,int_min,int_max);
      hadc_int[i]->GetXaxis()->SetTitle("Integral of Average fADC Value");
      hadc_int[i]->GetXaxis()->CenterTitle();
      hadc_int[i]->GetYaxis()->SetTitle("Occurrences");
      hadc_int[i]->GetYaxis()->CenterTitle();
    }

  TH1F **hmax_adc_tdc = new TH1F*[channels];
  for(Int_t i = 0; i<channels; i++)
    {
      hmax_adc_tdc[i] = new TH1F(Form("hmax%d_adc_tdc",i),Form("Maximum fADC Values for Module %d (Vertical ADC and TDC Cuts)",i),max_bins,max_min,max_max);
      hmax_adc_tdc[i]->GetXaxis()->SetTitle("Maximum fADC Value");
      hmax_adc_tdc[i]->GetXaxis()->CenterTitle();
      hmax_adc_tdc[i]->GetYaxis()->SetTitle("Occurrences");
      hmax_adc_tdc[i]->GetYaxis()->CenterTitle();
    }

  //Store the integral of the fADC pulse.
  TH1F **hint = new TH1F*[channels];
  for(Int_t i = 0; i<channels; i++)
    {
      hint[i] = new TH1F(Form("hint%d",i),Form("Integral of fADC Values for Module %d",i),int_bins,int_min,int_max);
      hint[i]->GetXaxis()->SetTitle("Integral of fADC Value");
      hint[i]->GetXaxis()->CenterTitle();
      hint[i]->GetYaxis()->SetTitle("Occurrences");
      hint[i]->GetYaxis()->CenterTitle();
    }

  TH1F **hintegrals = new TH1F*[channels];
  for(Int_t i = 0; i<channels; i++)
    {
      hintegrals[i] = new TH1F(Form("hintegrals%d",i),Form("Integral of fADC Values for Module %d",i),1000,0,30000);
      hintegrals[i]->GetXaxis()->SetTitle("Integral of Landau Fit of fADC Value");
      hintegrals[i]->GetXaxis()->CenterTitle();
      hintegrals[i]->GetYaxis()->SetTitle("Occurrences");
      hintegrals[i]->GetYaxis()->CenterTitle();
    }

  TH2F **hintvpeak = new TH2F*[channels];
  for(Int_t i = 0; i<channels; i++)
    {
      hintvpeak[i] = new TH2F(Form("hintvpeak%d",i),Form("Integral of fADC Values Versus Peak fADC Value for Module %d",i),1000,0,30000,1000,0,10000);
      hintvpeak[i]->GetXaxis()->SetTitle("Integral of Landau Fit of fADC Value");
      hintvpeak[i]->GetXaxis()->CenterTitle();
      hintvpeak[i]->GetYaxis()->SetTitle("Peak fADC Value");
      hintvpeak[i]->GetYaxis()->CenterTitle();
    }

  gStyle->SetOptStat(0);

  TH1F *hntdc_hits = new TH1F(Form("hntdc_hits"),Form("Number of TDC Hits per PMT Module"),channels,0,channels);
  hntdc_hits->GetXaxis()->SetTitle("Module Number (from 0)");
  hntdc_hits->GetXaxis()->CenterTitle();
  hntdc_hits->GetYaxis()->SetTitle("Occurrences");
  hntdc_hits->GetYaxis()->CenterTitle();
  TH1F *hnadc_hits = new TH1F(Form("hnadc_hits"),Form("Number of ADC Hits per PMT Module"),channels,0,channels);
  hnadc_hits->GetXaxis()->SetTitle("Module Number (from 0)");
  hnadc_hits->GetXaxis()->CenterTitle();
  hnadc_hits->GetYaxis()->SetTitle("Occurrences");
  hnadc_hits->GetYaxis()->CenterTitle();
  TH1F *hnadc_tdc_hits = new TH1F(Form("hnadc_tdc_hits"),Form("Number of ADC w/ TDC Hits per PMT Module"),channels,0,channels);
  hnadc_tdc_hits->GetXaxis()->SetTitle("Module Number (from 0)");
  hnadc_tdc_hits->GetXaxis()->CenterTitle();
  hnadc_tdc_hits->GetYaxis()->SetTitle("Occurrences");
  hnadc_tdc_hits->GetYaxis()->CenterTitle();
  TH1F *hntdc_hits_vert = new TH1F(Form("hntdc_hits_vert"),Form("Number of Vertical TDC Hits per PMT Module"),channels,0,channels);
  hntdc_hits_vert->GetXaxis()->SetTitle("Module Number (from 0)");
  hntdc_hits_vert->GetXaxis()->CenterTitle();
  hntdc_hits_vert->GetYaxis()->SetTitle("Occurrences");
  hntdc_hits_vert->GetYaxis()->CenterTitle();
  TH1F *hnadc_hits_vert = new TH1F(Form("hnadc_hits_vert"),Form("Number of Vertical ADC Hits per PMT Module"),channels,0,channels);
  hnadc_hits_vert->GetXaxis()->SetTitle("Module Number (from 0)");
  hnadc_hits_vert->GetXaxis()->CenterTitle();
  hnadc_hits_vert->GetYaxis()->SetTitle("Occurrences");
  hnadc_hits_vert->GetYaxis()->CenterTitle();
  TH1F *hnadc_tdc_hits_vert = new TH1F(Form("hnadc_tdc_hits_vert"),Form("Number of Vertical ADC w/ TDC Hits per PMT Module"),channels,0,channels);
  hnadc_tdc_hits_vert->GetXaxis()->SetTitle("Module Number (from 0)");
  hnadc_tdc_hits_vert->GetXaxis()->CenterTitle();
  hnadc_tdc_hits_vert->GetYaxis()->SetTitle("Voltage (-V)");
  hnadc_tdc_hits_vert->GetYaxis()->CenterTitle();
  TH1F *hhv = new TH1F(Form("hhv"),Form("High Voltage Setting per PMT Module"),channels,0,channels);
  hhv->GetXaxis()->SetTitle("Module Number (from 0)");
  hhv->GetXaxis()->CenterTitle();
  hhv->GetYaxis()->SetTitle("Occurrences");
  hhv->GetYaxis()->CenterTitle();
  TH1F *hadc_amp_avg = new TH1F(Form("hadc_amp_avg"),Form("Average Vertical fADC Cosmic Signal (RAU) per PMT Module"),channels,0,channels);
  hadc_amp_avg->GetXaxis()->SetTitle("Module Number (from 0)");
  hadc_amp_avg->GetXaxis()->CenterTitle();
  hadc_amp_avg->GetYaxis()->SetTitle("Average fADC Cosmic Signal (RAU)");
  hadc_amp_avg->GetYaxis()->CenterTitle();

  Int_t hits[channels] = {};
  Int_t no_tdc[channels] = {};

  //Loop over all events.
  if(limit_evts==1)
    {
      loop_max = max_evts;
    }
  else
    {
      loop_max = nevt;
    }

  for(Int_t i=0; i<loop_max ;i++)
    {
      T->GetEntry(gCurrentEntry);

      //Set arrays to zero;
      for(r  = 0; r < kNrows; r++) 
	{
	  for(c  = 0; c < kNcols; c++) 
	    {
	      peak[r][c] = 0.0;
	      adc[r][c] = 0.0;
	      tdc[r][c] = 0.0;
	    }
	}

      Double_t fadc_int[channels] = {};
      //Loop over all channels and fill the arrays. Do this seperately so these can be cut on later.
      for(Int_t j=0; j<hcalt::ndata; j++)
	{
	  r = hcalt::row[j]-1;
	  c = hcalt::col[j]-1;
	  idx = hcalt::samps_idx[j];
	  n = hcalt::nsamps[j];
	  adc[r][c] = hcalt::a[j];
	  tdc[r][c] = hcalt::tdc[j];
	  for(Int_t s = DISP_MIN_SAMPLE; s < DISP_MAX_SAMPLE && s < n; s++) 
	    {
	      histos[r][c]->SetBinContent(s+1-DISP_MIN_SAMPLE,hcalt::samps[idx+s]);
	    }
	}

      for(Int_t j=0; j<hcalt::ndata; j++)
	{
	  r = hcalt::row[j]-1;
	  c = hcalt::col[j]-1;
	  idx = hcalt::samps_idx[j];
	  n = hcalt::nsamps[j];

	  //Record if each PMT had the tdc fire.
	  if(tdc[r][c]!=0)
	    {
	      ntdc_hits[j] = ntdc_hits[j] + 1;
	    }

	  //Record if each PMT had the adc fire.
	  if(histos[r][c]->GetBinContent(histos[r][c]->GetMaximumBin())>(ped_avg[j]+adc_fire_avg*ped_avg_std_dev[j]))
	    {
	      nadc_hits[j] = nadc_hits[j] + 1;
	    }

	  //Record if each PMT had the adc and tdc fire.
	  if(histos[r][c]->GetBinContent(histos[r][c]->GetMaximumBin())>(ped_avg[j]+adc_fire_avg*ped_avg_std_dev[j]) && tdc[r][c]!=0)
	    {
	      nadc_tdc_hits[j] = nadc_tdc_hits[j] + 1;
	    }

	  if( adc[r][c]>(ped_int[j]+adc_fire_int*ped_int_std_dev[j]) )
	    {
	      nadc_int_hits[j] = nadc_int_hits[j] + 1;
	    }
	  
	  //cout<<Form("Event = %d, PMT = %d = (row %d,col %d). tdc(r-4,c)=%f, tdc(r-3,c)=%f, tdc(r-2,c)=%f, tdc(r-1,c)=%f, tdc(r,c)=%f, tdc(r+1,c)=%f, tdc(r+2,c)=%f, tdc(r+3,c)=%f, tdc(r+4,c)=%f",gCurrentEntry,j,r,c,tdc[r-4][c],tdc[r-3][c],tdc[r-2][c],tdc[r-1][c],tdc[r][c],tdc[r+1][c],tdc[r+2][c],tdc[r+3][c],tdc[r+4][c])<<endl;
	  //cout<<Form("Event = %d, PMT = %d = (row %d,col %d), tdc[r][c] = %f, tdc[%d][%d]=%f",gCurrentEntry,j,r,c,tdc[r][c],r-1,c,tdc[r-1][c])<<endl;

	  //Fill the histogram with the peak value of the fADC with a vertical requiement.
	  //Require 3 PMT TDCs vertically fire in a row.
	  //if(   (  ((tdc[r][c]!=0 && tdc[r-1][c]!=0 && tdc[r+1][c]!=0) || (tdc[r][c]!=0 && tdc[r-1][c]!=0 && tdc[r-2][c]!=0) || (tdc[r][c]!=0 && tdc[r+1][c]!=0 && tdc[r+2][c]!=0))  )   )
	  
	  if(   
	     ( (r-1)>=0 && (r+1)<=(kNrows-1) && tdc[r][c]!=0 && tdc[r-1][c]!=0 && tdc[r+1][c]!=0 ) 
	     || 
	     ( (r-1)>=0 && (r-2)>=0 && tdc[r][c]!=0 && tdc[r-1][c]!=0 && tdc[r-2][c]!=0 ) 
	     || 
	     ( (r+1)<=(kNrows-1) && (r+2)<=(kNrows-1) && tdc[r][c]!=0 && tdc[r+1][c]!=0 && tdc[r+2][c]!=0 )   
	    )
	  
	  //Require 3 PMT ADCs vertically fire in a row using integral of ADC (not finished).
	  //if(   (  ((adc[r][c]>(ped_avg[j]+adc_fire) && adc[r-1][c]>(ped_avg[j]+adc_fire) && adc[r+1][c]>(ped_avg[j]+adc_fire)) || (adc[r][c]>(ped_avg[j]+adc_fire) && adc[r-1][c]>(ped_avg[j]+adc_fire) && adc[r-2][c]>(ped_avg[j]+adc_fire)) || (adc[r][c]>(ped_avg[j]+adc_fire) && adc[r+1][c]>(ped_avg[j]+adc_fire) && adc[r+2][c]>(ped_avg[j]+adc_fire)))  )   )
	  
	    //Require 3 PMTs vertically fire in a row and the surrounding 6 do not fire.
	    //if(   (  ((tdc[r][c]!=0 && tdc[r-1][c]!=0 && tdc[r+1][c]!=0 && tdc[r][c-1]==0 && tdc[r-1][c-1]==0 && tdc[r+1][c-1]==0 && tdc[r][c+1]==0 && tdc[r-1][c+1]==0 && tdc[r+1][c+1]==0) || (tdc[r][c]!=0 && tdc[r-1][c]!=0 && tdc[r-2][c]!=0 && tdc[r][c-1]==0 && tdc[r-1][c-1]==0 && tdc[r-2][c-1]==0 && tdc[r][c+1]==0 && tdc[r-1][c+1]==0 && tdc[r-2][c+1]==0) || (tdc[r][c]!=0 && tdc[r+1][c]!=0 && tdc[r+2][c]!=0 && tdc[r][c-1]==0 && tdc[r+1][c-1]==0 && tdc[r+2][c-1]==0 && tdc[r][c+1]==0 && tdc[r+1][c+1]==0 && tdc[r+2][c+1]==0))  )   )
	  //Require 4 PMTs vertically fire in a row.
	  //if(   (  ((tdc[r+3][c]!=0 && tdc[r+2][c]!=0 && tdc[r+1][c]!=0 && tdc[r][c]!=0) || (tdc[r+2][c]!=0 && tdc[r+1][c]!=0 && tdc[r][c]!=0 && tdc[r-1][c]!=0 ) || (tdc[r+1][c]!=0 && tdc[r][c]!=0 && tdc[r-1][c]!=0 && tdc[r-2][c]!=0) || (tdc[r][c]!=0 && tdc[r-1][c]!=0 && tdc[r-2][c]!=0 && tdc[r-3][c]!=0))  )   )
	    //Require 5 PMTs vertically fire in a row.
	    //if(   (tdc[r+4][c]!=0 && tdc[r+3][c]!=0 && tdc[r+2][c]!=0 && tdc[r+1][c]!=0 && tdc[r][c]!=0 && (r+1)<kNrows) || (tdc[r+3][c]!=0 && tdc[r+2][c]!=0 && tdc[r+1][c]!=0 && tdc[r][c]!=0 && tdc[r-1][c]!=0 && (r+1)<kNrows && (r-1)>-1) || (tdc[r+2][c]!=0 && tdc[r+1][c]!=0 && tdc[r][c]!=0 && tdc[r-1][c]!=0 && tdc[r-2][c]!=0 && (r+1)<kNrows && (r-1)>-1) || (tdc[r+1][c]!=0 && tdc[r][c]!=0 && tdc[r-1][c]!=0 && tdc[r-2][c]!=0 && tdc[r-3][c]!=0 && (r+1)<kNrows && (r-1)>-1) || (tdc[r][c]!=0 && tdc[r-1][c]!=0 && tdc[r-2][c]!=0 && tdc[r-3][c]!=0 && tdc[r-4][c]!=0 && (r-1)>-1)   )
	    {
	      hmax_tdc[j]->Fill(histos[r][c]->GetBinContent(histos[r][c]->GetMaximumBin())-ped_avg[j]); 
	      ntdc_hits_vert[j] = ntdc_hits_vert[j] + 1;
	    }

	  //Require 3 PMT ADCs vertically fire in a row.	  
	  if(   
	     ( (r-1)>=0 && (r+1)<=(kNrows-1) && histos[r][c]->GetBinContent(histos[r][c]->GetMaximumBin())>(ped_avg[j]+adc_fire_avg*ped_avg_std_dev[j]) && histos[r-1][c]->GetBinContent(histos[r-1][c]->GetMaximumBin())>(ped_avg[j-12]+adc_fire_avg*ped_avg_std_dev[j-12]) && histos[r+1][c]->GetBinContent(histos[r+1][c]->GetMaximumBin())>(ped_avg[j+12]+adc_fire_avg*ped_avg_std_dev[j+12]) ) 
	     ||
	     ( (r-1)>=0 && (r-2)>=0 && histos[r][c]->GetBinContent(histos[r][c]->GetMaximumBin())>(ped_avg[j]+adc_fire_avg*ped_avg_std_dev[j]) && histos[r-1][c]->GetBinContent(histos[r-1][c]->GetMaximumBin())>(ped_avg[j-12]+adc_fire_avg*ped_avg_std_dev[j-12]) && histos[r-2][c]->GetBinContent(histos[r-2][c]->GetMaximumBin())>(ped_avg[j-24]+adc_fire_avg*ped_avg_std_dev[j-24]) )
	     ||
	     ( (r+1)<=(kNrows-1) && (r+2)<=(kNrows-1) && histos[r][c]->GetBinContent(histos[r][c]->GetMaximumBin())>(ped_avg[j]+adc_fire_avg*ped_avg_std_dev[j]) && histos[r+1][c]->GetBinContent(histos[r+1][c]->GetMaximumBin())>(ped_avg[j+12]+adc_fire_avg*ped_avg_std_dev[j+12]) && histos[r+2][c]->GetBinContent(histos[r+2][c]->GetMaximumBin())>(ped_avg[j+24]+adc_fire_avg*ped_avg_std_dev[j+24]) )
	    )
	    {
	      hmax_adc[j]->Fill(histos[r][c]->GetBinContent(histos[r][c]->GetMaximumBin())-ped_avg[j]); 
	      nadc_hits_vert[j] = nadc_hits_vert[j] + 1;
	    }

	  //Require 1 ADC above and below to fire ADC (from integral).	  
	  if(   
	     ( (r-1)>=0 && (r+1)<=(kNrows-1) && adc[r][c]>0 && adc[r-1][c]>(ped_int[j-12]+adc_fire_avg*ped_int_std_dev[j-12]) && adc[r+1][c]>(ped_int[j+12]+adc_fire_avg*ped_int_std_dev[j+12]) ) 
	     ||
	     ( (r-1)>=0 && (r-2)>=0 && adc[r][c]>0 && adc[r-1][c]>(ped_int[j-12]+adc_fire_avg*ped_int_std_dev[j-12]) && adc[r-2][c]>(ped_int[j-24]+adc_fire_avg*ped_int_std_dev[j-24]) )
	     ||
	     ( (r+1)<=(kNrows-1) && (r+2)<=(kNrows-1) && adc[r][c]>0 && adc[r+1][c]>(ped_int[j+12]+adc_fire_avg*ped_int_std_dev[j+12]) && adc[r+2][c]>(ped_int[j+24]+adc_fire_avg*ped_int_std_dev[j+24]) )
	    )
	    {
	      hadc_int[j]->Fill(adc[r][c]); 
	      nadc_int_hits_vert[j] = nadc_int_hits_vert[j] + 1;
	    }

	  //Require 3 PMT ADCs and TDCs vertically fire in a row.
	  if(   
	     ( (r-1)>=0 && (r+1)<=(kNrows-1) && histos[r][c]->GetBinContent(histos[r][c]->GetMaximumBin())>(ped_avg[j]+adc_fire_avg*ped_avg_std_dev[j]) && histos[r-1][c]->GetBinContent(histos[r-1][c]->GetMaximumBin())>(ped_avg[j-12]+adc_fire_avg*ped_avg_std_dev[j-12]) && histos[r+1][c]->GetBinContent(histos[r+1][c]->GetMaximumBin())>(ped_avg[j+12]+adc_fire_avg*ped_avg_std_dev[j+12]) && tdc[r][c]!=0 && tdc[r-1][c]!=0 && tdc[r+1][c]!=0 ) 
	     ||
	     ( (r-1)>=0 && (r-2)>=0 && histos[r][c]->GetBinContent(histos[r][c]->GetMaximumBin())>(ped_avg[j]+adc_fire_avg*ped_avg_std_dev[j]) && histos[r-1][c]->GetBinContent(histos[r-1][c]->GetMaximumBin())>(ped_avg[j-12]+adc_fire_avg*ped_avg_std_dev[j-12]) && histos[r-2][c]->GetBinContent(histos[r-2][c]->GetMaximumBin())>(ped_avg[j-24]+adc_fire_avg*ped_avg_std_dev[j-24]) && tdc[r][c]!=0 && tdc[r-1][c]!=0 && tdc[r-2][c]!=0 )
	     ||
	     ( (r+1)<=(kNrows-1) && (r+2)<=(kNrows-1) && histos[r][c]->GetBinContent(histos[r][c]->GetMaximumBin())>(ped_avg[j]+adc_fire_avg*ped_avg_std_dev[j]) && histos[r+1][c]->GetBinContent(histos[r+1][c]->GetMaximumBin())>(ped_avg[j+12]+adc_fire_avg*ped_avg_std_dev[j+12]) && histos[r+2][c]->GetBinContent(histos[r+2][c]->GetMaximumBin())>(ped_avg[j+24]+adc_fire_avg*ped_avg_std_dev[j+24]) && tdc[r][c]!=0 && tdc[r+1][c]!=0 && tdc[r+2][c]!=0 )
	    )
	    {
	      hmax_adc_tdc[j]->Fill(histos[r][c]->GetBinContent(histos[r][c]->GetMaximumBin())-ped_avg[j]); 
	      nadc_tdc_hits_vert[j] = nadc_tdc_hits_vert[j] + 1;
	    }

	}
      
      if(i%5000==0)
	{
	  cout<<i<<" events processed. "<<((double)i/(double)loop_max)*100.<<" % complete."<<endl;
	}
      gCurrentEntry++;
    }

  //Open file to store results histograms etc.
  TFile* file = new TFile(Form("/home/daq/test_fadc/Voltage_Scans/Cosmics/Average_Cosmic_Amplitudes_run%d.root",run),"RECREATE");

  TCanvas* cSingle=new TCanvas("cSingle");
  cSingle->SetGrid();
  histos[5][6]->Draw();

  //Create a date object.
  TDatime time;

  //Open file to write output to.
  std::ofstream output1 (Form("/home/daq/test_fadc/Voltage_Scans/Cosmics/Calibration_Results_Run%d.txt",run), std::ofstream::out);
  output1<<"These pedestals were generated from run "<<run<<"."<<endl;
  output1<<"They were generated at "<<time.GetHour()<<":"<<time.GetMinute()<<" on "<<time.GetMonth()<<"/"<<time.GetDay()<<"/"<<time.GetYear()<<"."<<endl;
  output1<<"PMT Module   Average Cosmic RAU (tdc cut)   Target Voltage (tdc cut)   Average Cosmic RAU (adc cut)   Target Voltage (adc cut)   Average Cosmic RAU (adc & tdc cut)   Target Voltage (adc & tdc cut)   TDC Hits   ADC Hits   TDC/ADC Hits   ADC & TDC Hits   Vert TDC Hits   Vert ADC Hits   Vertical TDC/ADC Hits   Vert ADC & TDC Hits"<<endl;

  TF1 **func_gaus_fit = new TF1*[channels];
  TF1 **func_gaus_fit_max = new TF1*[channels];
  //TF1 **func_weibull_max = new TF1*[channels];
  TF1 **func_skew_max_tdc = new TF1*[channels];
  TF1 **func_skew_max_adc = new TF1*[channels];
  TF1 **func_skew_max_adc_tdc = new TF1*[channels];
  TF1 **func_adc_int = new TF1*[channels];
  TF1 **func_gaus_fit_int = new TF1*[channels];
  for(Int_t i = 0; i<hcalt::ndata; i++)
    {
      //Assign HCal rows and cols (counting from 1) to each PMT module.
      r=(i/12)+1;
      c=i-12*(r-1)+1;

      func_skew_max_tdc[i] = new TF1(Form("func_skew_max_tdc%d",i),fit_skew, max_min, max_max, 4);
      func_skew_max_tdc[i]->SetLineColor(4);
      func_skew_max_tdc[i]->SetNpx(1000);

      func_skew_max_tdc[i]->SetRange(hmax_tdc[i]->GetXaxis()->GetBinCenter(hmax_tdc[i]->GetMaximumBin())-0.5*hmax_tdc[i]->GetStdDev(),hmax_tdc[i]->GetXaxis()->GetBinCenter(hmax_tdc[i]->GetMaximumBin())+1.3*hmax_tdc[i]->GetStdDev());
      func_skew_max_tdc[i]->SetParameter(0,hmax_tdc[i]->GetBinContent(hmax_tdc[i]->GetMaximumBin()));
      //func_skew_max_tdc[i]->SetParLimits(0,0,10000);
      func_skew_max_tdc[i]->SetParameter(1,hmax_tdc[i]->GetXaxis()->GetBinCenter(hmax_tdc[i]->GetMaximumBin()));
      //func_skew_max_tdc[i]->SetParLimits(1,0,8000);
      func_skew_max_tdc[i]->SetParameter(2,hmax_tdc[i]->GetStdDev());
      //func_skew_max_tdc[i]->SetParameter(2,100);
      func_skew_max_tdc[i]->SetParameter(3,2);
      //cout<<Form("Fit parameters for PMT %d (vertical tdc cut)",i)<<endl;
      hmax_tdc[i]->Fit(func_skew_max_tdc[i],"rq");
      hmax_tdc[i]->Write();

      avg_amp_tdc[i] = func_skew_max_tdc[i]->GetMaximumX();

      func_skew_max_adc[i] = new TF1(Form("func_skew_max_adc%d",i),fit_skew, max_min, max_max, 4);
      func_skew_max_adc[i]->SetLineColor(4);
      func_skew_max_adc[i]->SetNpx(1000);

      Double_t hmean_adc = hmax_adc[i]->GetMean();
      Double_t hstd_dev_adc = hmax_adc[i]->GetStdDev();
      Int_t range_min_bin = (int) ( (hmean_adc-1.*hstd_dev_adc - max_min)*(max_bins/(max_max-max_min)) );
      if(range_min_bin<1)
	{
	  range_min_bin=1;
	}
      Int_t range_max_bin = (int) ( (hmean_adc+1.*hstd_dev_adc - max_min)*(max_bins/(max_max-max_min)) );
      if(range_max_bin>max_bins)
	{
	  range_max_bin = max_bins;
	}
      hmax_adc[i]->GetXaxis()->SetRange(range_min_bin,range_max_bin);//Range is in bins!
      Int_t est_amp_adc = hmax_adc[i]->GetBinContent(hmax_adc[i]->GetMaximumBin());
      Int_t est_mean_adc = hmax_adc[i]->GetXaxis()->GetBinCenter(hmax_adc[i]->GetMaximumBin());
      Double_t est_std_dev_adc = hmax_adc[i]->GetStdDev();
      Int_t func_min_bin_adc = (int) (est_mean_adc - est_std_dev_adc);
      Int_t func_max_bin_adc = (int) (est_mean_adc + est_std_dev_adc);
      cout<<Form("PMT %d: histo mean = %f, histo std dev = %f, range min = %f = bin %d, range max = %f = bin %d, est_amp_adc = %d, est_mean_adc = %d, est_std_dev_adc = %f, func_min_bin_adc = %d, func_max_bin_adc = %d",i,hmean_adc,hstd_dev_adc,hmean_adc-1.*hstd_dev_adc,range_min_bin,hmean_adc+1.*hstd_dev_adc,range_max_bin,est_amp_adc,est_mean_adc,est_std_dev_adc,func_min_bin_adc,func_max_bin_adc)<<endl;

      func_skew_max_adc[i]->SetRange(func_min_bin_adc,func_max_bin_adc);
      //func_skew_max_adc[i]->SetRange(hmax_adc[i]->GetXaxis()->GetBinCenter(hmax_adc[i]->GetMaximumBin())-0.5*hmax_adc[i]->GetStdDev(),hmax_adc[i]->GetXaxis()->GetBinCenter(hmax_adc[i]->GetMaximumBin())+1.3*hmax_adc[i]->GetStdDev());
      func_skew_max_adc[i]->SetParameter(0,est_amp_adc);
      //func_skew_max_adc[i]->SetParLimits(0,0,10000);
      func_skew_max_adc[i]->SetParameter(1,est_mean_adc);
      //func_skew_max_adc[i]->SetParLimits(1,0,8000);
      func_skew_max_adc[i]->SetParameter(2,est_std_dev_adc);
      //func_skew_max_adc[i]->SetParameter(2,100);
      func_skew_max_adc[i]->SetParameter(3,2);
      hmax_adc[i]->GetXaxis()->SetRange(0,max_bins);
      cout<<Form("Fit parameters for PMT %d (vertical adc cut)",i)<<endl;
      hmax_adc[i]->Fit(func_skew_max_adc[i],"rq");
      hmax_adc[i]->Write();

      avg_amp_adc[i] = func_skew_max_adc[i]->GetMaximumX();

      func_skew_max_adc_tdc[i] = new TF1(Form("func_skew_max_adc_tdc%d",i),fit_skew, max_min, max_max, 4);
      func_skew_max_adc_tdc[i]->SetLineColor(4);
      func_skew_max_adc_tdc[i]->SetNpx(1000);

      func_skew_max_adc_tdc[i]->SetRange(hmax_adc_tdc[i]->GetXaxis()->GetBinCenter(hmax_adc_tdc[i]->GetMaximumBin())-0.5*hmax_adc_tdc[i]->GetStdDev(),hmax_adc_tdc[i]->GetXaxis()->GetBinCenter(hmax_adc_tdc[i]->GetMaximumBin())+1.3*hmax_adc_tdc[i]->GetStdDev());
      func_skew_max_adc_tdc[i]->SetParameter(0,hmax_adc_tdc[i]->GetBinContent(hmax_adc_tdc[i]->GetMaximumBin()));
      //func_skew_max_adc_tdc[i]->SetParLimits(0,0,10000);
      func_skew_max_adc_tdc[i]->SetParameter(1,hmax_adc_tdc[i]->GetXaxis()->GetBinCenter(hmax_adc_tdc[i]->GetMaximumBin()));
      //func_skew_max_adc_tdc[i]->SetParLimits(1,0,8000);
      func_skew_max_adc_tdc[i]->SetParameter(2,hmax_adc_tdc[i]->GetStdDev());
      //func_skew_max_adc_tdc[i]->SetParameter(2,100);
      func_skew_max_adc_tdc[i]->SetParameter(3,2);
      //cout<<Form("Fit parameters for PMT %d (vertical adc and tdc cut)",i)<<endl;
      hmax_adc_tdc[i]->Fit(func_skew_max_adc_tdc[i],"rq");
      hmax_adc_tdc[i]->Write();

      avg_amp_adc_tdc[i] = func_skew_max_adc_tdc[i]->GetMaximumX();

      func_adc_int[i] = new TF1(Form("func_adc_int%d",i),fit_gaus, max_min, max_max, 3);
      func_adc_int[i]->SetLineColor(4);
      func_adc_int[i]->SetNpx(1000);

      //func_adc_int[i]->SetRange(hadc_int[i]->GetXaxis()->GetBinCenter(hadc_int[i]->GetMaximumBin())-0.5*hadc_int[i]->GetStdDev(),hadc_int[i]->GetXaxis()->GetBinCenter(hadc_int[i]->GetMaximumBin())+1.3*hadc_int[i]->GetStdDev());
      func_adc_int[i]->SetRange(ped_int[i]+4*ped_int_std_dev[i],hadc_int[i]->GetMean()+1.5*hadc_int[i]->GetStdDev());
      Int_t func_min_bin = (int) (ped_int[i]+4*ped_int_std_dev[i]-int_min)/((int_max-int_min)/int_bins);
      Int_t func_min_bin_value = hadc_int[i]->GetBinContent(func_min_bin);
      Int_t func_max_bin = (int) (hadc_int[i]->GetMean()+1.5*hadc_int[i]->GetStdDev()-int_min)/((int_max-int_min)/int_bins);
      Int_t func_max_bin_value = hadc_int[i]->GetBinContent(func_max_bin);
      hadc_int[i]->GetXaxis()->SetRange(func_min_bin,func_max_bin);     //Range is in bins!
      Int_t est_amp = hadc_int[i]->GetBinContent(hadc_int[i]->GetMaximumBin());
      func_adc_int[i]->SetParameter(0,est_amp);
      //cout<<Form("Estimated amplitude = %d",est_amp)<<endl;
      //func_adc_int[i]->SetParameter(0,hadc_int[i]->GetBinContent(hadc_int[i]->GetMaximumBin(ped_int[i]+4*ped_int_std_dev[i],hadc_int[i]->GetMean()+1.5*hadc_int[i]->GetStdDev()),0));
      //cout<<Form("Bin at start of fit func = %.2f   value in that bin = %.2f",(ped_int[i]+4*ped_int_std_dev[i]-int_min)/((int_max-int_min)/int_bins),hadc_int[i]->GetBinContent((int) (ped_int[i]+4*ped_int_std_dev[i]-int_min)/((int_max-int_min)/int_bins)) )<<endl;
      //func_adc_int[i]->SetParLimits(0,0,10000);
      func_adc_int[i]->SetParameter(1,hadc_int[i]->GetXaxis()->GetBinCenter(hadc_int[i]->GetMaximumBin()));
      //func_adc_int[i]->SetParLimits(1,0,8000);
      func_adc_int[i]->SetParameter(2,hadc_int[i]->GetStdDev());
      //func_adc_int[i]->SetParameter(2,100);
      func_adc_int[i]->SetParameter(3,2);
      //cout<<Form("Fit parameters for PMT %d (vertical adc integral cuts)",i)<<endl;
      hadc_int[i]->GetXaxis()->SetRange(0,int_bins);
      hadc_int[i]->Fit(func_adc_int[i],"rq");
      hadc_int[i]->Write();

      avg_amp_adc_int[i] = func_adc_int[i]->GetMaximumX();

      //Calculate the target voltage to set all PMT signals to equal strength for a similar signal.
      if(c==5||c==6||c==7||c==8)
	{
	  target_voltage_tdc[i] = hv[i]/pow(avg_amp_tdc[i]/target_RAU,1./jlab_exp);
	  target_voltage_adc[i] = hv[i]/pow(avg_amp_adc[i]/target_RAU,1./jlab_exp);
	}
      else
	{
	  target_voltage_tdc[i] = hv[i]/pow(avg_amp_tdc[i]/target_RAU,1./cmu_exp);
	  target_voltage_adc[i] = hv[i]/pow(avg_amp_adc[i]/target_RAU,1./cmu_exp);
	}

      cout<<Form("Average cosmic amplitude for PMT %d (tdc cut) is %.2f.",i,avg_amp_tdc[i])<<endl;
      output1<<Form("%d     %.2f    %.0f",i,avg_amp_tdc[i],-target_voltage_tdc[i]);

      cout<<Form("Average cosmic amplitude for PMT %d (adc cut) is %.2f.",i,avg_amp_adc[i])<<endl;
      output1<<Form("    %.2f    %.0f",avg_amp_adc[i],-target_voltage_adc[i]);

      cout<<Form("Average cosmic amplitude for PMT %d (adc and tdc cut) is %.2f.",i,avg_amp_adc_tdc[i])<<endl;
      output1<<Form("    %.2f    %.0f",avg_amp_adc_tdc[i],-target_voltage_adc_tdc[i]);

      output1<<Form("    %d    %d    %.2f    %d    %d    %d    %.2f    %d",ntdc_hits[i],nadc_hits[i],(double)ntdc_hits[i]/(double)nadc_hits[i],nadc_tdc_hits[i],ntdc_hits_vert[i],nadc_hits_vert[i],(double)ntdc_hits_vert[i]/(double)nadc_hits_vert[i],nadc_tdc_hits_vert[i])<<endl;

      if(use_int==1)
	{
	  func_gaus_fit_int[i] = new TF1("func_gaus_fit_int",fit_gaus, int_min, int_max, 3);
	  func_gaus_fit_int[i]->SetLineColor(2);
	  func_gaus_fit_int[i]->SetNpx(1000);
	  
	  //func_gaus_fit_int[i]->SetRange(hint[i]->GetXaxis()->GetBinCenter(hint[i]->GetMaximumBin())-20,hint[i]->GetXaxis()->GetBinCenter(hint[i]->GetMaximumBin())+20);
	  func_gaus_fit_int[i]->SetRange(hint[i]->GetXaxis()->GetBinCenter(hint[i]->GetMaximumBin())-50,hint[i]->GetXaxis()->GetBinCenter(hint[i]->GetMaximumBin())+50);
	  func_gaus_fit_int[i]->SetParameter(0,hint[i]->GetBinContent(hint[i]->GetMaximumBin()));
      //func_gaus_fit_int[i]->SetParLimits(0,0,10000);
	  func_gaus_fit_int[i]->SetParameter(1,hint[i]->GetXaxis()->GetBinCenter(hint[i]->GetMaximumBin()));
	  //func_gaus_fit_int[i]->SetParLimits(1,0,8000);
	  func_gaus_fit_int[i]->SetParameter(2,hint[i]->GetStdDev());
	  hint[i]->Fit(func_gaus_fit_int[i],"rq");
	  hint[i]->Write();
	}
    }

  for(Int_t i=0;i<channels;i++)
    {
      hntdc_hits->SetBinContent(i+1,ntdc_hits[i]);
      hnadc_hits->SetBinContent(i+1,nadc_hits[i]);
      hnadc_tdc_hits->SetBinContent(i+1,nadc_tdc_hits[i]);
      hntdc_hits_vert->SetBinContent(i+1,ntdc_hits_vert[i]);
      hnadc_hits_vert->SetBinContent(i+1,nadc_hits_vert[i]);
      hnadc_tdc_hits_vert->SetBinContent(i+1,nadc_tdc_hits_vert[i]);
      hhv->SetBinContent(i+1,hv[i]*-1);
      hadc_amp_avg->SetBinContent(i+1,avg_amp_adc[i]);
      if(ntdc_hits[i]==0)
	{
	  cout<<Form("The TDC for PMT module %d (from zero) never fired!",i)<<endl;
	}
      if(nadc_hits[i]==0)
	{
	  cout<<Form("The ADC for PMT module %d (from zero) never fired!",i)<<endl;
	}
    }

  hntdc_hits->Write();
  hnadc_hits->Write();
  hnadc_tdc_hits->Write();
  hntdc_hits_vert->Write();
  hnadc_hits_vert->Write();
  hnadc_tdc_hits_vert->Write();
  hhv->Write();
  hadc_amp_avg->Write();

  file->Close();

  //Stop and print out stopwatch information.
  st->Stop();
  cout<<"CPU time = "<<st->CpuTime()<<" s = "<<st->CpuTime()/60.<<" min   Real time = "<<st->RealTime()<<" s = "<<st->RealTime()/60.<<" min"<<endl;
}
