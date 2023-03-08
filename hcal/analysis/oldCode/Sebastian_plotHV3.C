//Version 3. Plotting max ADC, HV, and NPE data - sseeds 3.30.21

#include <TH2F.h>
#include <TChain.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <TSystem.h>
#include <TStopwatch.h>
#include <ctime>
#include "hcal.h"

const int kNrows = 24; //Number of pmt rows
const int kNcols = 12; //Number of pmt columns
const int kNLED = 6; //Number of led channels = 5 leds + off = 6 channels
const int runT = 10; //Total number of runs
const int chan = 144; //Number of active channels

TF1 *t1;
TF1 *t2;
TF1 *t3;
TF1 *t4;
TF1 *t5;

TFile *LEDpulse_HVcal = new TFile("outFiles/aggLEDPlotsV3.root","RECREATE");

double maxMean[kNrows][kNcols][kNLED][runT];
double maxSD[kNrows][kNcols][kNLED][runT];
double HVset[kNrows][kNcols][runT];
double NPE_f1[kNrows][kNcols][kNLED][runT]; //For use plotting NPE v Run
double NPE_f2[kNrows][kNcols][runT][kNLED]; //For use plotting NPE v LED
double maxRes[kNrows][kNcols][kNLED][runT];
double rChiSqrMax[kNrows][kNcols][kNLED];
double alpha[kNrows][kNcols];
double avgNPE[kNrows*kNcols];
double SDNPE[kNrows*kNcols];
double channel[kNrows*kNcols];

double LED[5] = {1,2,3,4,5}; //Array for plotting

TH1F *maxVar[kNLED];
TH1F *maxAlpha[kNrows][kNcols];
TH1F *NPEl1[kNrows][kNcols];
TH1F *maxAlphaVsChannel;
TH1F *maxSdVsChannel;
TH1F *NPEVsChannel;
TH1F *jlabAlphas;
TH1F *cmuAlphas;

double powFit(double *X,double *p) //Single parameter fits - GOING FOR THE GOLD!
{
  double fitval = p[0]*pow(X[0]/p[1],p[2]); //X normalized to lowest HV value
  return fitval;
}

void LEDplots_v3(int run1 = 1288, int run2 = 1289, int run3 = 1290, int run4 = 1291, int run5 = 1292, int run6 = 1293, int run7 = 1294, int run8 = 1297, int run9 = 1300, int run10 = 1301){ 

  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start(kTRUE);
  
  cout << "Beginning loop over runs " << run1 << ", " << run2 << ", " << run3 << ", " << run4 << ", " << run5 << ", " << run6 << ", " << run7 << ", " << run8 << ", " << run9 << ", and  " << run10 << ".." << endl;  
  int runs[runT]={run1,run2,run3,run4,run5,run6,run7,run8,run9,run10};

  //Create array for TGraphErrors over channels
  for(int c=0; c<(kNrows*kNcols); c++){
    channel[c]=c;
  }

  cout << "HV settings and LED intensities hardcoded from HCAL wiki for now.." << endl;

  int runHV_jlab[runT] = {1500,1550,1600,1650,1700,1750,1800,1850,1900,1950};
  int runHV_cmu[runT] = {1400,1425,1450,1475,1500,1525,1550,1575,1600,1625};

  //Loop over all runs and fill master array
  cout << "Filling master array from analysis histograms.." << endl;
  for(int i=0; i<runT; i++){

    TFile *f = TFile::Open(Form("outFiles/ledSpectFile_run%d.root",runs[i]));  //Read in file
    
    for(int r=12; r<kNrows; r++){
      for(int c=0; c<kNcols; c++){

	TH1F *h3 = (TH1F*)f->Get(Form("Pedestal Spect R%d C%d",r,c));
	
	for(int l=1; l<kNLED; l++){ //Ignoring led = 0, as no led is on for this setting
	  	  
	  TH1F *h2 = (TH1F*)f->Get(Form("Max ADC Spect R%d C%d LED%d",r,c,l));  //Save all integrated spectra to master histo

	  if(h2){
	    h2->Fit("gaus","Q");
	    t2=h2->GetFunction("gaus");
	    maxMean[r][c][l][i] = t2->GetParameter(1);
	    maxSD[r][c][l][i] = t2->GetParameter(2);
	    if(h3){
	      h3->Fit("gaus","Q");
	      t3=h3->GetFunction("gaus");
	      NPE_f1[r][c][l][i]=pow(t2->GetParameter(1)/pow(pow(t2->GetParameter(2),2.0)-pow(t3->GetParameter(2),2.0),0.5),2.0); //Number of photoelectrons calculated from max ADC currently
	    }
	  }	  
	}	
      }
    }

    f->Close();
    cout << "Ending run " << runs[i] << ".." << endl;
  }

  ofstream outFile;
  outFile.open("setFiles/alphas.txt");
  time_t now = time(0);
  char *dt = ctime(&now);
  outFile << "#Alpha parameters obtained from LEDplots_v2 on " << dt << "#" << endl;

  //Now to draw TGraphs for gain curves
  cout << "Using master array to draw gain curves and to obtain fit parameters.." << endl;
  for(int r=12; r<kNrows; r++){
    for(int c=0; c<kNcols; c++){

      //Use maxMean to draw graphs
      maxAlpha[r][c] = new TH1F(Form("R%dC%dMaxAlpha",r,c),Form("R%dC%dMaxAlpha",r,c),50,0,0);
      for(int l=1; l<kNLED; l++){
	double HV[runT];
	double Max[runT];
	double SD[runT];
	for(int i=0; i<runT; i++){
	  if(c>3&&c<8){ //jlab tubes
	    HV[i]=runHV_jlab[i];
	    HVset[r][c][i]=runHV_jlab[i];
	  }else{ //cmu tubes
	    HV[i]=runHV_cmu[i];
	    HVset[r][c][i]=runHV_cmu[i];  
	  }	  
	  Max[i]=maxMean[r][c][l][i];
	  SD[i]=maxSD[r][c][l][i]/sqrt(5000);  //Error is (std dev)/sqrt(N) (or is it the std dev of the pedestal?)
	}
	TGraph *G1 = new TGraphErrors(runT,HV,Max,0,SD);	
	TF1 *pFit = new TF1("pFit",powFit,1300,2100,3);

	if(c>3&&c<8){ //jlab tubes

	  pFit->FixParameter(0,maxMean[r][c][l][0]);
	  pFit->FixParameter(1,runHV_jlab[0]);
	  pFit->SetParameter(2,5);
	  pFit->SetParLimits(2,1,10);
	  
	}else{ //cmu tubes

	  pFit->FixParameter(0,maxMean[r][c][l][0]);
	  pFit->FixParameter(1,runHV_cmu[0]);
	  pFit->SetParameter(2,22);
	  pFit->SetParLimits(2,1,30);

	}

	//Fit the gain Curve and obtain reduced chi square from the fit
	G1->Fit("pFit","Q");
	rChiSqrMax[r][c][l] = pFit->GetChisquare()/(runT-1);

	//Draw the TGraph and write to file
	G1->SetTitle(Form("MAX R%d-C%d LED%d: ADCref=%f HVref=%f alpha=%f rChiSqr=%f",r,c,l,pFit->GetParameter(0),pFit->GetParameter(1),pFit->GetParameter(2),rChiSqrMax[r][c][l]));
	G1->GetYaxis()->SetTitle("Average fADC Amplitude (RAU)");  //1 RAU approx 0.5 mV
	G1->GetYaxis()->CenterTitle();
	G1->GetXaxis()->SetTitle("PMT HV Setting (-V)");
	G1->GetXaxis()->CenterTitle();
	G1->SetMinimum(0.0);
	G1->SetLineColor(kWhite);
	G1->SetMarkerStyle(8);
	G1->SetMarkerSize(1);
	LEDpulse_HVcal->cd();

	if(l==1) G1->Write(Form("MAX R%d-C%d LED%d",r,c,l));

	maxAlpha[r][c]->Fill(pFit->GetParameter(2));

	if(l==1) alpha[r][c]=pFit->GetParameter(2); //Obtain exponent from LED 1 only
	
	G1->Draw("AP");

	for(int i=0; i<runT; i++){
	  maxRes[r][c][l][i] = Max[i]-pFit->Eval(HV[i]);
	}
      }      
    }
  }

  //Define another array to handle the NPE vs HV TGraph
  for(int r=12; r<kNrows; r++){
    for(int c=0; c<kNcols; c++){
      NPEl1[r][c] = new TH1F(Form("R%dC%d NPE LED1",r,c),Form("R%dC%d NPE LED1",r,c),50,0,0);
      for(int run=0; run<runT; run++){	
	for(int l=1; l<kNLED; l++){
	  NPE_f2[r][c][run][l-1]=NPE_f1[r][c][l][run];
	  if(l==1) NPEl1[r][c]->Fill(NPE_f1[r][c][l][run]);
	}
      }
    }
  }
  
  //Draw NPE v HV graphs to look for plateau region
  for(int r=12; r<kNrows; r++){
    for(int c=0; c<kNcols; c++){
      
      for(int l=1; l<kNLED; l++){
	TGraph *G3 = new TGraph(runT, HVset[r][c], NPE_f1[r][c][l]);
	G3->SetTitle(Form("R%d-C%d LED %d NPE v HV",r,c,l));
	G3->GetYaxis()->SetTitle("NPE");
	G3->GetYaxis()->CenterTitle();
	G3->GetXaxis()->SetTitle("HV Setting");
	G3->GetXaxis()->CenterTitle();
	G3->SetMinimum(0.0);
	G3->SetLineColor(kWhite);
	G3->SetMarkerStyle(8);
	G3->SetMarkerSize(1);
	LEDpulse_HVcal->cd();
	G3->Write(Form("R%d-C%d LED %d NPE v HV",r,c,l));
	G3->Draw("AP");

	TGraph *G4 = new TGraph(runT, HVset[r][c], maxRes[r][c][l]);
	G4->SetTitle(Form("R%d-C%d LED %d Max Residuals",r,c,l));
	G4->GetYaxis()->SetTitle("Residuals");
	G4->GetYaxis()->CenterTitle();
	G4->GetXaxis()->SetTitle("HV Setting");
	G4->GetXaxis()->CenterTitle();
	G4->SetLineColor(kWhite);
	G4->SetMarkerStyle(8);
	G4->SetMarkerSize(1);
	LEDpulse_HVcal->cd();
	G4->Write(Form("R%d-C%d LED %d Max Residuals",r,c,l));
	G4->Draw("AP");
      }

      for(int run=0; run<runT; run++){
	TGraph *G5 = new TGraph(kNLED, LED, NPE_f2[r][c][run]);
	G5->SetTitle(Form("R%d-C%d run %d NPE v LED",r,c,runs[run]));
	G5->GetYaxis()->SetTitle("NPE");
	G5->GetYaxis()->CenterTitle();
	G5->GetXaxis()->SetTitle("LED");
	G5->GetXaxis()->CenterTitle();
	G5->SetLineColor(kWhite);
	G5->SetMarkerStyle(8);
	G5->SetMarkerSize(1);
	LEDpulse_HVcal->cd();
	G5->Write(Form("R%d-C%d run %d NPE v LED",r,c,runs[run]));
	G5->Draw("AP");

      }
    }
  }

  // General analysis histograms
  
  maxAlphaVsChannel = new TH1F("Average Max Alpha vs Channel","Average Max Alpha vs Channel",kNrows*kNcols,0,kNrows*kNcols-1);
  maxSdVsChannel = new TH1F("Max SD vs Channel","Max SD vs Channel",kNrows*kNcols,0,kNrows*kNcols-1);
  
  jlabAlphas = new TH1F("Jlab Alphas","Jlab Alphas",40,0,15);
  cmuAlphas = new TH1F("CMU Alphas","CMU Alphas",40,15,30);
 
  for(int l=0; l<kNLED; l++){
    maxVar[l] = new TH1F(Form("MaxADC Chi Sqr vs Channel, LED %d",l),Form("MaxADC Chi Sqr vs Channel, LED %d",l),kNrows*kNcols,0,kNrows*kNcols-1);
  } 

  double cmuSum=0.0;
  double jlabSum=0.0;

  for(int c=0; c<144; c++){ //Temporary until LED runs can be taken over both halves of HCal
    avgNPE[c]=0.0;
    SDNPE[c]=0.0;
  }
  
  for(int c=144; c<(kNrows*kNcols); c++){
    
    int rval = floor(c/kNcols);
    int cval = c%kNcols;

    avgNPE[c]=NPEl1[rval][cval]->GetMean();
    SDNPE[c]=NPEl1[rval][cval]->GetRMS();

    //cout << avgNPE[c] << "  " << SDNPE[c] << endl;
        
    maxAlphaVsChannel->SetBinContent(c,alpha[rval][cval]);

    if(cval>3&&cval<8){
      jlabSum+=alpha[rval][cval];
    }else{
      cmuSum+=alpha[rval][cval];
    }
    
    maxSdVsChannel->SetBinContent(c,maxAlpha[rval][cval]->GetRMS());
    for(int l=0; l<kNLED; l++){
      maxVar[l]->SetBinContent(c,rChiSqrMax[rval][cval][l]);
    }
  }

  //Send alpha parameters to file. 0-144 are temporary averages pending LED runs on RHS.
  for(int c=0; c<144; c++){
    int cval = c%kNcols;
    if(cval>3&&cval<8){
      outFile << jlabSum/48 << endl; //For now just take an average of all alphas for jlab tubes
    }else{
      outFile << cmuSum/96 << endl; //For now just take an average of all alphas for cmu tubes
    }
  }

  for(int c=144; c<(kNrows*kNcols); c++){
    int rval = floor(c/kNcols);
    int cval = c%kNcols;
    outFile << alpha[rval][cval] << endl;
    
    if(cval>3&&cval<8){
      jlabAlphas->Fill(alpha[rval][cval]);
    }else{
      cmuAlphas->Fill(alpha[rval][cval]);
    }
    
  }

  TGraph *NPEvChan = new TGraphErrors(kNrows*kNcols,channel,avgNPE,0,SDNPE);	
  NPEvChan->SetTitle("Average NPE vs Channel (Error SD)");
  NPEvChan->GetYaxis()->SetTitle("<NPE>");
  NPEvChan->GetYaxis()->CenterTitle();
  NPEvChan->GetXaxis()->SetTitle("Channel");
  NPEvChan->GetXaxis()->CenterTitle();
  //NPEvChan->SetLineColorAlpha(kRed,0.0);
  NPEvChan->SetMarkerStyle(8);
  NPEvChan->SetMarkerSize(1);
  LEDpulse_HVcal->cd();
  NPEvChan->Write("NPEvChan");
  NPEvChan->Draw("AP");
  
  LEDpulse_HVcal->cd();
  maxAlphaVsChannel->SetTitle("Average Max Alpha vs Channel");
  maxAlphaVsChannel->Write("maxAlphaVsChannel");
  maxAlphaVsChannel->Draw("AP");
  
  LEDpulse_HVcal->cd();
  maxSdVsChannel->SetTitle("Max SD vs Channel");
  maxSdVsChannel->Write("maxSdVsChannel");
  maxSdVsChannel->Draw("AP");

  LEDpulse_HVcal->cd();
  jlabAlphas->SetTitle("Jlab Tube Alpha Values, LED 1");
  jlabAlphas->Write("jlabAlphas");
  jlabAlphas->Draw("AP");
  
  LEDpulse_HVcal->cd();
  cmuAlphas->SetTitle("CMU Tube Alpha Values, LED 1");
  cmuAlphas->Write("cmuAlphas");
  cmuAlphas->Draw("AP");

  for(int l=1; l<kNLED; l++){
    LEDpulse_HVcal->cd();
    maxVar[l]->Write(Form("maxVar LED%d",l));
    maxVar[l]->Draw("AP");
  }
  
  cout << "Completed loop over runs. Histograms written to LEDpulse_HVcal.root. Fit parameters written to file PMTCalibrationData.txt. Alphas written to alphas.txt." << endl;

  st->Stop();

  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

}
