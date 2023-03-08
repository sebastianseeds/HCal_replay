//First attempt at plotting ave int charge vs HV setting for LED pulse data - sseeds 12.10.20

#include <TH2F.h>
#include <TChain.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <TSystem.h>
#include "hcal.h"

const int kNrows = 24; //Number of pmt rows
const int kNcols = 12; //Number of pmt columns
const int runT = 5; //Five total runs read in by function
double chiSqr = 0.0;

TFile *F1 = new TFile("outFiles/cosmicPresPlots.root","RECREATE");

double pmtHV[kNrows][kNcols][runT];
double avgADC[kNrows][kNcols][runT];
double avgTDC[kNrows][kNcols][runT];

TH1F *ChannelADC[runT];
TH1F *ChannelTDC[runT];

TH1F *ADCvChannel[runT];
TH1F *TDCvChannel[runT];

TH1F *ChannelPed[runT];

TGraph *ADCvHV[kNrows][kNcols];
TGraph *TDCvHV[kNrows][kNcols];

double powFit(double *x,double *par){ //"Single parameter" fits
  double ADC_norm = par[0]; //Currently fixed to ADC value of lowest HV point
  double HV_norm = par[1]; //Currently fixed to lowest HV value
  double alpha = par[2];
  double HV = x[0];
  return ADC_norm*pow(HV/HV_norm,alpha);
}

//Landau fit
double landauFit(double *x, double *par) {
  double amp = par[0];
  double mpv = par[1];
  double sigma = par[2];
  double offset = par[3];
  double ADC = x[0];
  return amp * TMath::Landau(ADC,mpv,sigma) + offset;
}


void cosCalPlots(int run1 = 1239, int run2 = 1263, int run3 = 1265, int run4 = 1269, int run5 = 1271){ 
  
  cout << "Beginning loop over runs " << run1 << ", " << run2 << ", " << run3 << ", " << run4 << ", and " << run5 << ".." << endl;  
  int runs[runT]={run1,run2,run3,run4,run5};
  
  cout << "Reading in HV settings per pmt for each run.." << endl;

  for(int i=0; i<runT; i++){
    
    ifstream file(Form("setFiles/HV_run%d.txt",runs[i]));

    cout << "Getting HV settings for each pmt for run " << runs[i] << "." << endl;
    
    int n1=0;
    double d1;
    
    int rval, cval;
    string line;
    
    while (getline(file,line)){
      if(line.at(0) == '#'){
	continue;
      }
      
      stringstream ss(line);
      ss >> d1;

      //cout << d1 << endl;
      
      rval = floor(n1/kNcols);
      cval = n1 % kNcols;
      
      pmtHV[rval][cval][i] = -d1; 
      
      cout << "run = " << runs[i] << ", row = " << rval << ", col = " << cval << ", array val = " << pmtHV[rval][cval][i] << endl;
      
      n1++;
    }
  }

  //Loop over runs and populate ADC and TDC histograms
  cout<< "Filling histograms with average ADC and TDC values by run.." << endl;
  for(int i=0; i<runT; i++){

    cout << "Filling histos from run " << runs[i] << ".." << endl;
    
    TFile *f = TFile::Open(Form("outFiles/cosmicHistograms_run%d.root",runs[i]));  //Read in file

    TH1F *h0 = (TH1F*)f->Get("pedvChannel");
    TH1F *h1 = (TH1F*)f->Get("ADCvChannel"); //Get ADC Histogram
    TH1F *h2 = (TH1F*)f->Get("TDCvChannel"); //Get TDC Histogram

    TF1 *lFit = new TF1("lFit",powFit,1200,2100,3); //Define fit function

    ChannelPed[i] = new TH1F(Form("Avg Pedestal, All Channels, Run %d",runs[i]),Form("Avg Pedestal, All Channels, Run %d",runs[i]),50,0,300);
    
    ChannelADC[i] = new TH1F(Form("Avg ADC, All Channels, Run %d",runs[i]),Form("Avg ADC, All Channels, Run %d",runs[i]),50,0,200);
    ChannelTDC[i] = new TH1F(Form("Avg TDC, All Channels, Run %d",runs[i]),Form("Avg TDC, All Channels, Run %d",runs[i]),50,0,200);

    ADCvChannel[i] = new TH1F(Form("ADC vs Channel Run %d",runs[i]),Form("ADC vs Channel Run %d",runs[i]),kNrows*kNcols,0,288);
    TDCvChannel[i] = new TH1F(Form("TDC vs Channel Run %d",runs[i]),Form("TDC vs Channel Run %d",runs[i]),kNrows*kNcols,0,288);

    for(int j=0;j<h1->GetEntries();j++){
      
      int row = floor(j/kNcols);
      int col = (j%kNcols);
      
      TH1F *h3 = (TH1F*)f->Get(Form("Max ADC Spect R%d C%d",row,col));
      TH1F *h4 = (TH1F*)f->Get(Form("Max ADC Spect TDC R%d C%d",row,col));
    
      if(h1->GetBinContent(j+1)<10||h3->GetMean()<10||h1->GetBinContent(j+1)<150){
	F1->cd();
	h3->Write(Form("A ChannelADC R%d C%d Run%d",row,col,runs[i]));
	h3->Draw("AP");
	//F1->cd();
	//h4->Write(Form("A ChannelTDC R%d C%d Run%d",row,col,runs[i]));
	//h4->Draw("AP");    
      }

      ChannelPed[i]->Fill(h0->GetBinContent(j+1));

      //if(h1->GetBinContent(j+1)<250&&h1->GetBinContent(j+1)>20){

      if(h1->GetBinContent(j+1)>20){
	
	ChannelADC[i]->Fill(h1->GetBinContent(j+1));
	ChannelTDC[i]->Fill(h2->GetBinContent(j+1));

	ADCvChannel[i]->SetBinContent(j+1,h1->GetBinContent(j+1));
	TDCvChannel[i]->SetBinContent(j+1,h2->GetBinContent(j+1));

	avgADC[row][col][i]=h1->GetBinContent(j+1);
	avgTDC[row][col][i]=h2->GetBinContent(j+1);
      }
      
      if(h1->GetBinContent(j+1)<20){
		
	ChannelADC[i]->Fill(h3->GetMean());
	ChannelTDC[i]->Fill(h2->GetMean());
	
	ADCvChannel[i]->SetBinContent(j+1,h3->GetMean());
	TDCvChannel[i]->SetBinContent(j+1,h4->GetMean());

	avgADC[row][col][i-1]=h3->GetMean();
	avgTDC[row][col][i-1]=h4->GetMean();
      }
      /*
      if(h1->GetBinContent(j+1)>250){

	ChannelADC[i]->Fill(h3->GetMaximumBin()*10+5);
	ChannelTDC[i]->Fill(h2->GetMaximumBin()*10+5);
	
	ADCvChannel[i]->SetBinContent(j+1,h3->GetMaximumBin()*10+5);
	TDCvChannel[i]->SetBinContent(j+1,h4->GetMaximumBin()*10+5);

	avgADC[row][col][i]=h3->GetMaximumBin()*10+5;
	avgTDC[row][col][i]=h4->GetMaximumBin()*10+5;
      }
      */
    }
  }
  
  //Draw gain curve tgraphs by pmt
  for(int r=0; r<kNrows; r++){
    for(int c=0; c<kNcols; c++){
      

      TF1 *pFit = new TF1("pFit",powFit,1200,2100,3); //Define fit function
	
      //Set limits for fit based on pmt type (jlab or cmu)
      /*
      if(c>3&&c<8){
	pFit->SetParameters(0,10,8,0);
	pFit->SetParLimits(1,0,1500);
	pFit->SetParLimits(2,1,100);
      }else{
	pFit->SetParameters(0,800,10,10);
      }
      */

      if(c>3&&c<8){ //jlab tubes
	
	pFit->FixParameter(0,avgADC[r][c][0]);
	pFit->FixParameter(1,pmtHV[r][c][0]);
	pFit->SetParameter(2,5);
	pFit->SetParLimits(2,1,10);
	
      }else{ //cmu tubes

	pFit->FixParameter(0,avgADC[r][c][0]);
	pFit->FixParameter(1,pmtHV[r][c][0]);
	pFit->SetParameter(2,22);
	pFit->SetParLimits(2,1,30);
	
      }

      //Draw ADC graph
      ADCvHV[r][c] = new TGraph(runT,pmtHV[r][c],avgADC[r][c]); 
      ADCvHV[r][c]->Fit("pFit","Q");
      ADCvHV[r][c]->SetTitle(Form("Avg ADC vs HV R%d-C%d: coeff=%g shift=%f alpha=%f offset=%f",r,c,pFit->GetParameter(0),pFit->GetParameter(1),pFit->GetParameter(2),pFit->GetParameter(3)));
      ADCvHV[r][c]->GetYaxis()->SetTitle("Avg Max ADC (RAU)");  
      ADCvHV[r][c]->GetYaxis()->CenterTitle();
      ADCvHV[r][c]->GetXaxis()->SetTitle("PMT HV Setting (-V)");
      ADCvHV[r][c]->GetXaxis()->CenterTitle();
      ADCvHV[r][c]->SetMinimum(0.0);
      ADCvHV[r][c]->SetLineColor(kWhite);
      ADCvHV[r][c]->SetMarkerStyle(8);
      ADCvHV[r][c]->SetMarkerSize(1);
      F1->cd();
      ADCvHV[r][c]->Write(Form("Avg ADC R%d-C%d: coeff=%g shift=%f alpha=%f offset=%f",r,c,pFit->GetParameter(0),pFit->GetParameter(1),pFit->GetParameter(2),pFit->GetParameter(3)));
      ADCvHV[r][c]->Draw("AP");
    
      /*
      //Draw TDC graph
      TDCvHV[r][c] = new TGraph(runT,pmtHV[r][c],avgTDC[r][c]); 
      TDCvHV[r][c]->Fit("pFit","Q");
      TDCvHV[r][c]->SetTitle(Form("Avg TDC vs HV R%d-C%d: coeff=%g shift=%f alpha=%f offset=%f",r,c,pFit->GetParameter(0),pFit->GetParameter(1),pFit->GetParameter(2),pFit->GetParameter(3)));
      TDCvHV[r][c]->GetYaxis()->SetTitle("Avg Max ADC (RAU)"); 
      TDCvHV[r][c]->GetYaxis()->CenterTitle();
      TDCvHV[r][c]->GetXaxis()->SetTitle("PMT HV Setting (-V)");
      TDCvHV[r][c]->GetXaxis()->CenterTitle();
      TDCvHV[r][c]->SetMinimum(0.0);
      TDCvHV[r][c]->SetLineColor(kWhite);
      TDCvHV[r][c]->SetMarkerStyle(8);
      TDCvHV[r][c]->SetMarkerSize(1);
      F1->cd();
      TDCvHV[r][c]->Write(Form("Avg TDC R%d-C%d: coeff=%g shift=%f alpha=%f offset=%f",r,c,pFit->GetParameter(0),pFit->GetParameter(1),pFit->GetParameter(2),pFit->GetParameter(3)));
      TDCvHV[r][c]->Draw("AP");
      */
    }
  }
  
  //Draw ADC and TDC Histograms
  for(int i=0; i<runT; i++){
    ChannelADC[i]->SetTitle(Form("Avg ADC, All Channels, Run %d",runs[i]));
    ChannelADC[i]->GetYaxis()->SetTitle("Counts");  
    ChannelADC[i]->GetYaxis()->CenterTitle();
    ChannelADC[i]->GetXaxis()->SetTitle("Avg Max ADC (RAU)");
    ChannelADC[i]->GetXaxis()->CenterTitle();
    ChannelADC[i]->SetMinimum(0.0);
    ChannelADC[i]->SetMaximum(150.0);

    F1->cd();
    ChannelADC[i]->Write(Form("ChannelADC Run%d",runs[i]));
    ChannelADC[i]->Draw("AP");

    ChannelPed[i]->SetTitle(Form("Avg Pedestal, All Channels, Run %d",runs[i]));
    ChannelPed[i]->GetYaxis()->SetTitle("Counts");  
    ChannelPed[i]->GetYaxis()->CenterTitle();
    ChannelPed[i]->GetXaxis()->SetTitle("Avg Pedestal (RAU)");
    ChannelPed[i]->GetXaxis()->CenterTitle();
    ChannelPed[i]->SetMinimum(0.0);
    ChannelPed[i]->SetMaximum(150.0);

    F1->cd();
    ChannelPed[i]->Write(Form("ChannelPed Run%d",runs[i]));
    ChannelPed[i]->Draw("AP");

    ChannelTDC[i]->SetTitle(Form("Avg TDC, All Channels, Run %d",runs[i]));
    ChannelTDC[i]->GetYaxis()->SetTitle("Counts");  
    ChannelTDC[i]->GetYaxis()->CenterTitle();
    ChannelTDC[i]->GetXaxis()->SetTitle("Avg Max TDC (RAU)");
    ChannelTDC[i]->GetXaxis()->CenterTitle();
    ChannelTDC[i]->SetMinimum(0.0);
    ChannelTDC[i]->SetMaximum(150.0);
    
    F1->cd();
    ChannelTDC[i]->Write(Form("ChannelTDC Run%d",runs[i]));
    ChannelTDC[i]->Draw("AP");
    
    ADCvChannel[i]->SetTitle(Form("Avg ADC v Channels, Run %d",runs[i]));
    ADCvChannel[i]->SetMaximum(350.0);
    F1->cd();
    ADCvChannel[i]->Write(Form("ADCvChannel Run%d",runs[i]));
    ADCvChannel[i]->Draw("AP");
    
    TDCvChannel[i]->SetTitle(Form("Avg TDC v Channels, Run %d",runs[i]));
    TDCvChannel[i]->SetMaximum(350.0);
    F1->cd();
    TDCvChannel[i]->Write(Form("TDCvChannel Run%d",runs[i]));
    TDCvChannel[i]->Draw("AP");

  }

  
  cout << "Completed loop over runs. Histograms written to cosmicPresPlots.root. Fit parameters written to file cosPMTCalibrationData.txt." << endl;
}
