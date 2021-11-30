//SSeeds February 2021 - Script to output target HV settings from gain equation

#include <TH2F.h>
#include <TChain.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <TSystem.h>
#include <TStopwatch.h>
#include "hcal.h"

const int kNrows = 24; //Number of pmt rows
const int kNcols = 12; //Number of pmt columns
const int kNsamps = 30; //Total fADC samples taken per trigger
const int bitL = 4095; //Maximum ADC out value corresponding to 12 bits

int gCurrentEntry = -1;
TChain *T = 0;

double ampF = 1.5; //Dynamic range amplifier (3.0 V) / dynamic range fADC (2.0 V)
double attC = 2.697; //Average over all channels via measurements
double maxEDep = 700; //Maximum energy deposited per pmt in GMn at highest Q^2 (in MeV)
double ADCTarg;
double beamEDep;
double maxADC[kNrows][kNcols];
double pmtAlpha[kNrows][kNcols];
double pmtHVSet[kNrows][kNcols];
double pmtHVTarg[kNrows][kNcols];
double pedestals[kNrows][kNcols];
double pedSigma[kNrows][kNcols];

bool gSaturated[kNrows][kNcols];

TH1F *histos[kNrows][kNcols];
TH1F *PMTMaxSpec[kNrows][kNcols];
TH1F *pedSpec[kNrows][kNcols];

TF1 *f1;

//Function declarations
void processMaxADC(int,int);
void goodHistoTest(TH1F*, double, int, int);
void getPedestal(int);
void processEvent(int);


////////////////////////////////
//Create several fit functions// 
////////////////////////////////

//Gaussian fit (may just use built-in gaussian fit when convenient)
double gausFit(double *x, double *par) {
  double amp = par[0];
  double loc = par[1];
  double sigma = par[2];
  double ADC = x[0];
  return amp * TMath::Exp(-0.5 * pow(((ADC - loc) / sigma), 2));
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

//Poisson fit
double poissonFit(double *x, double *par) {
  double N = par[0];
  double mu = par[1];
  double mu1 = par[2];
  double ADC = x[0];
  return N * TMath::Poisson(ADC / mu1, mu);
}

//Power-law fit
double powFit(double *x, double *par) {
  double amp = par[0];
  double ex = par[1];
  double loc = par[2];
  double offset = par[3];
  double ADC = x[0];
  return amp * pow(ADC,ex); //No shift, no normalization
  //return amp*pow(ADC/1200.0,ex); //X normalized to lowest HV value
  //return amp*pow((ADC-loc)/1200,ex)+offset; //X normalized to lowest HV value, y-offset
}

//Skewed-normal fit
double skewFit(double *x, double *par) {
  double amp = par[0];
  double loc = par[1];
  double scale = par[2];
  double shape = par[3];
  double ADC = x[0];
  return amp * exp(-0.5 * pow((ADC - loc) / scale, 2.0)) * 0.5 * (1 + erf(shape * ((ADC - loc) / scale)));
}


TH1F* MakeHisto(int row, int col, int bins){
  TH1F *h = new TH1F(Form("h%02d%02d",row,col),Form("%d-%d",row+1,col+1),bins,0,kNsamps);
  h->SetStats(0);
  h->SetLineWidth(2);
  return h;
}

//Rows and Columns start at zero and go to kNrows-1 and kNcols-1
void processEvent(int entry = -1){
  if(entry == -1) {
    gCurrentEntry++;
  } else {
    gCurrentEntry = entry;
  }

  if(gCurrentEntry<0) {
    gCurrentEntry = 0;
  }

  T->GetEntry(gCurrentEntry);
 
  int r,c,idx,n,sub;
  
  // Clear old histograms, just in case modules are not in the tree
  for(r = 0; r < kNrows; r++) {
    for(c = 0; c < kNcols; c++) {
      histos[r][c]->Reset("ICES M");
      gSaturated[r][c] = false;
    }
  }

  // Create arrays for fADC and TDC values
  float peak[kNrows][kNcols];
  double adc[kNrows][kNcols];
  double tdc[kNrows][kNcols];
  for(r  = 0; r < kNrows; r++) {
    for(c  = 0; c < kNcols; c++) {
      peak[r][c] = 0.0;
      adc[r][c] = 0.0;
      tdc[r][c] = 0.0;
    }
  }

  // Fill fADC histograms by event
  for(int m = 0; m < hcalt::ndata; m++) {
    r = hcalt::row[m]-1;
    c = hcalt::col[m]-1;
    if(r>=kNrows || c>=kNcols) continue;
    idx = hcalt::samps_idx[m];
    n = hcalt::nsamps[m];
    adc[r][c] = hcalt::a[m];
    tdc[r][c] = hcalt::tdc[m];
    bool processed = false;
    for(int s=0; s<kNsamps && s<n; s++) {
      processed = true;
      histos[r][c]->SetBinContent(s+1,hcalt::samps[idx+s]);
      if(peak[r][c]<hcalt::samps[idx+s])
        peak[r][c]=hcalt::samps[idx+s];
      if(peak[r][c]>4095) {
        gSaturated[r][c] = true;
      }
    }
    if(!processed) {
      cerr << "Skipping empty module: " << m << ".." << endl;
      for(int s = 0; s<kNsamps; s++) {
        histos[r][c]->SetBinContent(s+1,-404);
      }
    }
  }

  // Test ADC samples for signal on event window by channel
  for(r = 0; r < kNrows; r++) {
    for(c = 0; c < kNcols; c++) {
      histos[r][c]->SetTitle(Form("R%d C%d (ADC=%g,TDC=%g)",r+1,c+1,adc[r][c],tdc[r][c]));
      if(!gSaturated[r][c]){
	goodHistoTest(histos[r][c],tdc[r][c],r,c);
      }
    }
  }

  // Fill max ADC spectra histograms
  for(r = 0; r < kNrows; r++) {
    for(c = 0; c < kNcols; c++) {
      PMTMaxSpec[r][c]->Fill(histos[r][c]->GetMaximum());
    }	
  }
}

// Acquire pedestal from first 4 bins of each histogram for each module
void getPedestal(int entry=-1){ 
  if(entry == -1) {
    gCurrentEntry++;
  } else {
    gCurrentEntry = entry;
  }
  
  if(gCurrentEntry<0) {
    gCurrentEntry = 0;
  }
  
  T->GetEntry(gCurrentEntry);
  
  int r,c,idx,n,sub;
  
  // Clear old histograms, just in case modules are not in the tree
  for(r=0; r<kNrows; r++) {
    for(c=0; c<kNcols; c++) {
      histos[r][c]->Reset("ICES M");
    }
  }

  // Create arrays for fADC and TDC values
  double adc[kNrows][kNcols];
  double tdc[kNrows][kNcols];
  for(r=0; r<kNrows; r++) {
    for(c=0; c<kNcols; c++) {
      adc[r][c] = 0.0;
      tdc[r][c] = 0.0;
    }
  }
  
  // Fill fADC histograms by event
  for(int m=0; m<hcalt::ndata; m++) {
    r = hcalt::row[m]-1;
    c = hcalt::col[m]-1;
    if(r>=kNrows || c>=kNcols) continue;
    idx = hcalt::samps_idx[m];
    n = hcalt::nsamps[m];
    adc[r][c] = hcalt::a[m];
    tdc[r][c] = hcalt::tdc[m];
    bool processed = false;
    for(int s=0; s<kNsamps && s<n; s++) {
      processed = true;
      histos[r][c]->SetBinContent(s+1,hcalt::samps[idx+s]);      
    }
  
    if(!processed) {
      cerr << "Skipping empty module: " << m << ".." << endl;
      for(int s=0; s<kNsamps; s++) {
        histos[r][c]->SetBinContent(s+1,-404);
      }
    }
  }
  
  // Fill pedestal histograms by channel
  for(r=0; r<kNrows; r++) {
    for(c=0; c<kNcols; c++) {      
      if(tdc[r][c]==0){ //eliminating adc pulses from coincident tdc measurement
	for(int j=0; j<kNsamps; j++) {
	  pedSpec[r][c]->Fill(histos[r][c]->GetBinContent(j+1)); //Fill pedestal histogram by channel
	}
      }
    }
  }
}

// Pedestal subtract then verify that max value is greater than twice sigma of pedestals and lies w/in bounds of histo
void goodHistoTest(TH1F *testHisto, double tdcVal, int row, int col){
  
  // Get pedestal value for the run and pmt 
  double pedVal = pedestals[row][col];
  
  // Subtract pedestal value from each bin
  for(int b=1; b<=testHisto->GetNbinsX(); b++) {
    testHisto->SetBinContent(b,testHisto->GetBinContent(b)-pedVal);
  }

  // Apply software cuts: Signal must be greater than 6sigma pedestal, peak must not be on edges of fADC samp window, TDC must cotrigger
  if(testHisto->GetMaximum()>(6.0*pedSigma[row][col]) && testHisto->GetMaximumBin()>2 && testHisto->GetMaximumBin()<(kNsamps-2) && tdcVal!=0){
    PMTMaxSpec[row][col]->Fill(testHisto->GetMaximum()); //Fill max fADC histograms from run by channel
  }
}

// Get average max ADC by channel for the run
void processMaxADC(int run, int event)
{
  
  //Build spectrum histograms. Empirical limits.
  cout << "Building spectrum histograms.." << endl;
  for(int r=0; r<kNrows; r++){
    for(int c=0; c<kNcols; c++){
      
      PMTMaxSpec[r][c] = new TH1F(Form("Max ADC Spect R%d C%d",r,c),Form("Max ADC Spect R%d C%d",r,c),100,0,1000);
      PMTMaxSpec[r][c]->GetXaxis()->SetTitle("RAU");
      PMTMaxSpec[r][c]->GetXaxis()->CenterTitle();
      pedSpec[r][c] = new TH1F(Form("Pedestal Spect R%d C%d",r,c),Form("Pedestal Spect R%d C%d",r,c),200,120,320);
      pedSpec[r][c]->GetXaxis()->SetTitle("<RAU>");
      pedSpec[r][c]->GetXaxis()->CenterTitle();
      
    }
  }

  cout << "Reading in run data.." << endl;
  
  if(!T) { 
    T = new TChain("T");

    //Read in data file. If multiple files exist: format = 'fadc_f1tdc_%d_n.root' where n runs from 0 to n-1 files.
    T->Add(TString::Format("rootFiles/fadc_f1tdc_%d.root",run));
    /*
    T->Add(TString::Format("rootFiles/fadc_f1tdc_%d_0.root",run));    
    T->Add(TString::Format("rootFiles/fadc_f1tdc_%d_1.root",run));
    T->Add(TString::Format("rootFiles/fadc_f1tdc_%d_2.root",run));
    T->Add(TString::Format("rootFiles/fadc_f1tdc_%d_3.root",run));
    T->Add(TString::Format("rootFiles/fadc_f1tdc_%d_4.root",run));
    T->Add(TString::Format("rootFiles/fadc_f1tdc_%d_5.root",run));
    T->Add(TString::Format("rootFiles/fadc_f1tdc_%d_6.root",run));
    T->Add(TString::Format("rootFiles/fadc_f1tdc_%d_7.root",run));
    T->Add(TString::Format("rootFiles/fadc_f1tdc_%d_8.root",run));
    T->Add(TString::Format("rootFiles/fadc_f1tdc_%d_9.root",run));
    T->Add(TString::Format("rootFiles/fadc_f1tdc_%d_10.root",run));
    T->Add(TString::Format("rootFiles/fadc_f1tdc_%d_11.root",run));
    T->Add(TString::Format("rootFiles/fadc_f1tdc_%d_12.root",run));
    */
    T->SetBranchStatus("*",0);
    T->SetBranchStatus("sbs.hcal.*",1);
    T->SetBranchAddress("sbs.hcal.nsamps",hcalt::nsamps);
    T->SetBranchAddress("sbs.hcal.a",hcalt::a);
    T->SetBranchAddress("sbs.hcal.tdc",hcalt::tdc);
    T->SetBranchAddress("sbs.hcal.samps",hcalt::samps);
    T->SetBranchAddress("sbs.hcal.samps_idx",hcalt::samps_idx);
    T->SetBranchAddress("sbs.hcal.row",hcalt::row);
    T->SetBranchAddress("sbs.hcal.col",hcalt::col);
    T->SetBranchStatus("Ndata.sbs.hcal.row",1);
    T->SetBranchAddress("Ndata.sbs.hcal.row",&hcalt::ndata);
    for(int r = 0; r < kNrows; r++) {
      for(int c = 0; c < kNcols; c++) {
        histos[r][c] = MakeHisto(r,c,kNsamps);
        gSaturated[r][c] = false;
      }
    }
  }
  
  gCurrentEntry = event;

  cout << "Total events to process " << T->GetEntries() << ".." << endl;

  cout << "Filling pedestal histograms.." << endl;
  for (int i=gCurrentEntry; i<T->GetEntries(); i++){
    getPedestal(gCurrentEntry);
    gCurrentEntry++;
    
    // Keep count of events processed for monitoring
    if (gCurrentEntry%10000 == 0){
      cout << "Current pedestal event: " << gCurrentEntry << endl;
    }
  }
  
  cout << "Fitting pedestal histos and extracting mean/sigma.." << endl;
  for(int r=0; r<kNrows; r++){
    for(int c=0; c<kNcols; c++){  
      pedSpec[r][c]->Fit("gaus","Q");
      f1=pedSpec[r][c]->GetFunction("gaus");
      pedestals[r][c]=f1->GetParameter(1);
      pedSigma[r][c]=f1->GetParameter(2);
    }
  }

  gCurrentEntry = event; // Resetting entries for signal analysis
  
  cout << "Processing signals.." << endl;
  for (int i = gCurrentEntry; i < T->GetEntries(); i++){
    processEvent(gCurrentEntry);
    gCurrentEntry++;
    
    // Keep count of events processed for monitoring
    if (gCurrentEntry%10000 == 0){
      cout << "Current event: " << gCurrentEntry << endl;
    }
  }

  TF1 *landauFitMax[kNrows][kNcols] = {};
  
  for(int r=0; r<kNrows; r++){
    for(int c=0; c<kNcols; c++){

      // Create landau fit function per r and c
      landauFitMax[r][c] = new TF1(Form("landauFitMax r%d c%d",r,c),landauFit, 0, 1000, 4);
      landauFitMax[r][c]->SetNpx(1000);

      // Set parameters for landau fit from spectrum histograms
      landauFitMax[r][c]->SetRange(PMTMaxSpec[r][c]->GetXaxis()->GetBinCenter(PMTMaxSpec[r][c]->GetMaximumBin())-0.8*PMTMaxSpec[r][c]->GetStdDev(),PMTMaxSpec[r][c]->GetXaxis()->GetBinCenter(PMTMaxSpec[r][c]->GetMaximumBin())+1.3*PMTMaxSpec[r][c]->GetStdDev());
      landauFitMax[r][c]->SetParameter(0,PMTMaxSpec[r][c]->GetBinContent(PMTMaxSpec[r][c]->GetMaximumBin()));
      landauFitMax[r][c]->SetParameter(1,PMTMaxSpec[r][c]->GetXaxis()->GetBinCenter(PMTMaxSpec[r][c]->GetMaximumBin()));
      landauFitMax[r][c]->SetParameter(2,PMTMaxSpec[r][c]->GetStdDev());
      landauFitMax[r][c]->SetParameter(3,0);

      // Fit max fADC histograms to obtain mean by channel
      if(PMTMaxSpec[r][c]->GetEntries()>0){
	PMTMaxSpec[r][c]->Fit(landauFitMax[r][c],"+RQ");
	if(landauFitMax[r][c]->GetParameter(1)>0&&landauFitMax[r][c]->GetParameter(1)<2000){ // Apply sanity checks to fit
	  maxADC[r][c] = landauFitMax[r][c]->GetParameter(1);
	}else{
	  maxADC[r][c] = PMTMaxSpec[r][c]->GetMean(); // Get mean of histogram if fit fails
	}
      }
    } 
  }
}
  
// Main
int HCal_HV_Settings(double ED, int run){
  
  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start(kTRUE);

  // Create file to hold target HV settings
  ofstream outFile;
  outFile.open(Form("HCal_HV_Settings_Run%d.txt",run));
  outFile << "#Target HV settings for run " << run << "." << endl;
  outFile << "#" << endl;
  outFile << "#Row Col targetHV" << endl;
  
  beamEDep = ED; // User input: average energy deposited in HCal pmt from beam (in MeV)
  
  cout << "Reading in alpha parameters.." << endl;

  // Read in from pmt alpha parameter file
  ifstream fileA("setFiles/PMT_alpha.txt");

  int n1 = 0;
  double d1;
  int r, c;
  string line1;
    
  while (getline(fileA,line1))
    {
      if(line1.at(0) == '#'){
	continue;
      }
      /*
      if(line1.at(0) == '!'){
	break;
      }
      */
      stringstream ss(line1);
      ss >> d1;
      
      r = floor(n1/kNcols);
      c = n1 % kNcols;
      
      pmtAlpha[r][c] = d1;
      
      cout << "PMT at row = " << r << ", col = " << c << ": alpha parameter val = " << pmtAlpha[r][c] << "." << endl;

      n1++;
    }

  cout << "Reading in HV settings.." << endl;
  
  //Read in from pmt HV setting file
  ifstream fileB("setFiles/PMT_HV.txt");

  cout << "Getting HV settings for each pmt for run." << endl;

  n1 = 0;
  string line2;
  
  while (getline(fileB,line2))
    {
      if(line2.at(0) == '#'){
	continue;
      }
      /*
      if(line1.at(0) == '!'){
	break;
      }
      */
      stringstream ss(line2);
      ss >> d1;
      
      r = floor(n1/kNcols);
      c = n1 % kNcols;
      
      pmtHVSet[r][c] = -d1; //Assuming negative HV inputs
      
      cout << "PMT at row = " << r << ", col = " << c << ": HV setting = " << pmtHVSet[r][c] << "." << endl;

      n1++;
    }

  cout << "Calculating ADC targets.." << endl;

  ADCTarg = (bitL*ampF*attC)*(beamEDep/maxEDep); //Calculate target ADC from beamEDep

  cout << "Analyzing data from run " << run << ".." << endl;

  processMaxADC(run,-1);

  cout << "Completed loop over run " << run << " data. Writing target HV settings to file.." << endl;
  
  for(int r=0; r<kNrows; r++){
    for(int c=0; c<kNcols; c++){
      pmtHVTarg[r][c] = pmtHVSet[r][c]/pow(maxADC[r][c]/ADCTarg,1.0/pmtAlpha[r][c]);

      outFile << r << "  " << c << "  " << pmtHVTarg[r][c] << endl; // Write target HV settings to file
    }
  }

  st->Stop();

  cout << "Analysis complete. Target HV settings written to file HCal_HV_Settings_Run" << run << ".txt." << endl << endl;
  
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

  return 0;
}

