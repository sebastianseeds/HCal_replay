//SSeeds 2020 - Cosmic PMT HV calibration code for use with SBS HCAL experiments. 24 rows x 12 cols HCAL modules. Relevant information can be found here: hallaweb.jlab.org/wiki/index.php/HCAL_Cosmics.

#include <TH2F.h>
#include <TChain.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <TSystem.h>
#include <TStopwatch.h>
#include "hcal.h"

const int kNrows = 24;
const int kNcols = 12;

const int DISP_MIN_SAMPLE = 0.0;
const int DISP_MAX_SAMPLE = 50.0;
const int DISP_FADC_SAMPLES = (DISP_MAX_SAMPLE-DISP_MIN_SAMPLE);

int gCurrentEntry = -1;
TChain *T = 0;

//Augment fn's for cosmic Calibration
void goodHistoTest(TH1F*, double, int, int);
void getPedestal(int);

TH1F *ADCvChannel;
TH1F *TDCvChannel;
TH1F *intADCvChannel;
TH1F *intTDCvChannel;
TH1F *NEVvChannel;
TH1F *pedvChannel;
TH1F *PPvCh;
TH1F *PPSDvCh;
TH1F *histos[kNrows][kNcols];
TH1F *pedSpec[kNrows][kNcols];
TH1F *PMTIntSpec[kNrows][kNcols]; //=new TH1F("","",20,0,5000); //Empirical limits
TH1F *PMTMaxSpec[kNrows][kNcols]; //=new TH1F("","",20,0,1000); //Empirical limits
TH1F *PMTIntSpecTDC[kNrows][kNcols]; //=new TH1F("","",20,0,5000); //Empirical limits, cut on tdc coincidence
TH1F *PMTMaxSpecTDC[kNrows][kNcols]; //=new TH1F("","",20,0,1000); //Empirical limits, cut on tdc coincidence
TH1F *PMTMaxPP[kNrows][kNcols];
TF1 *f1;
TF1 *f2;
TF1 *f3;

int signalTotal = 1;
TFile *HistosFile = new TFile("outFiles/HistosFile.root","RECREATE");  //File for checking histogram fits and integrity
TFile *HistosFilePedestal = new TFile("outFiles/HistosFilePedestal.root","RECREATE");  //File for checking histogram fits and integrity
double PPmean[kNrows][kNcols];
double PPsd[kNrows][kNcols];
double pedestals[kNrows][kNcols];
double pedSigma[kNrows][kNcols];
bool gSaturated[kNrows][kNcols];
bool gPulse[kNrows+4][kNcols+4]; //Needs to be larger for false-valued buffers (2 per side)
bool gPulseTDC[kNrows+4][kNcols+4]; //Needs to be larger for false-valued buffers (2 per side)
int DRAWNHISTO = 0; //Counts up per histo drawn
int histNote = 0;

double targetRAU = 61.425; //4095*3.0/2.0*0.5*14/700 or (digMaxADC)*(dynRanAmp/dynRanADC)*(cosEdep/maxGMnEdep)
int runN[9] = {989,988,987,984,983,981,980,978,1235}; //All run numbers

double pmtHV[kNrows][kNcols];
//double cmuExp = 10.5; //Left for comparison to Scott's results using these assumptions
//double jlabExp = 8.0; //Left for comparison to Scott's results using these assumptions
double alphas[kNrows][kNcols];
double targetHV[kNrows][kNcols];
double targetHVTDC[kNrows][kNcols];

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

//Skewed-landau fit
double skewLanFit(double *x, double *par) {
  double amp = par[0];
  double loc = par[1];
  double scale = par[2];
  double shape = par[3];
  double ADC = x[0];
  return amp * exp(-0.5 * pow((ADC - loc) / scale, 2.0)) * (1 + TMath::Landau(ADC,loc,scale));
}


TH1F* MakeHisto(int row, int col, int bins){
  TH1F *h = new TH1F(Form("h%02d%02d",row,col),Form("%d-%d",row+1,col+1),bins,DISP_MIN_SAMPLE,DISP_MAX_SAMPLE);
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
      gPulse[r+2][c+2] = false;
      gPulseTDC[r+2][c+2] = false;
    }
  }
  
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
  
  for(int m = 0; m < hcalt::ndata; m++) {
    r = hcalt::row[m]-1;
    c = hcalt::col[m]-1;
    if(r < 0 || c < 0) {
      std::cerr << "Why is row negative? Or col?" << std::endl;
      continue;
    }
    
    if(r>= kNrows || c >= kNcols) continue;

    //if(r==6||r==9) continue; //Temporary as fADC is in shop-------------------------------------------*
    
    idx = hcalt::samps_idx[m];
    n = hcalt::nsamps[m];
    adc[r][c] = hcalt::a[m];
    tdc[r][c] = hcalt::tdc[m];
    bool processed = false;
    for(int s = DISP_MIN_SAMPLE; s < DISP_MAX_SAMPLE && s < n; s++) {
      processed = true;
      histos[r][c]->SetBinContent(s+1-DISP_MIN_SAMPLE,hcalt::samps[idx+s]);
      if(peak[r][c]<hcalt::samps[idx+s])
        peak[r][c]=hcalt::samps[idx+s];
      if(peak[r][c]>4095) {
        gSaturated[r][c] = true;
      }
    }
    if(!processed) {
      std::cerr << "Skipping empty module: " << m << std::endl;
      for(int s = 0;  s < DISP_FADC_SAMPLES; s++) {
        histos[r][c]->SetBinContent(s+1,-404);
      }
    }
  }
  
  for(r = 0; r < kNrows; r++) {
    for(c = 0; c < kNcols; c++) {
      histos[r][c]->SetTitle(TString::Format("%d-%d (ADC=%g,TDC=%g)",r+1,c+1,adc[r][c],tdc[r][c]));
      if(!gSaturated[r][c]){
	goodHistoTest(histos[r][c],tdc[r][c],r,c);
      }
    }
  }

  //Vertical line test - remember to shift all values up by one due to intialization of gPulse/gPulseTDC and 'false' buffer. Only building functionality for straight vertical line test.

  for(r = 0; r < kNrows; r++) {
    for(c = 0; c < kNcols; c++) {
      if(gPulse[r+2][c+2]==true){
	if((gPulse[r][c+2]==true&&gPulse[r+1][c+2]==true)||(gPulse[r+1][c+2]==true&&gPulse[r+3][c+2]==true)||(gPulse[r+3][c+2]==true&&gPulse[r+4][c+2]==true)){ //Checks if two pulses exist above, below, or one above and below every given pulse to ensure a track exists
	  continue;
	}else{
	  gPulse[r+2][c+2]=false;
	}
      }
    }
  }
  
  //Now, if pulse passes verticality test, pedestal subtract and fill spectra histograms
  for(r = 0; r < kNrows; r++) {
    for(c = 0; c < kNcols; c++) {
      //if(gPulse[r+1][c+1]==true){

      //if(r==6||r==9) continue; //Temporary as fADC is in shop-------------------------------------------*
      
      if(gPulse[r+2][c+2]==true){ 
	
	PMTIntSpec[r][c]->Fill(histos[r][c]->Integral(1,DISP_MAX_SAMPLE)); //Integral from total bin content
	PMTMaxSpec[r][c]->Fill(histos[r][c]->GetMaximum());
	PMTMaxPP[r][c]->Fill((double)(histos[r][c]->GetMaximumBin()));
	
	//Only fill TDC if both ADC and TDC pulses present
	if(gPulseTDC[r+2][c+2]==true){
	  PMTIntSpecTDC[r][c]->Fill(histos[r][c]->Integral(1,DISP_MAX_SAMPLE)); //Integral from total bin content
	  PMTMaxSpecTDC[r][c]->Fill(histos[r][c]->GetMaximum());
	}	
      }	
    }
  }
}

//Acquire pedestal from first 4 bins of each histogram for each module
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
  for(r = 0; r < kNrows; r++) {
    for(c = 0; c < kNcols; c++) {
      histos[r][c]->Reset("ICES M");
    }
  }
  
  double adc[kNrows][kNcols];
  double tdc[kNrows][kNcols];
  for(r  = 0; r < kNrows; r++) {
    for(c  = 0; c < kNcols; c++) {
      adc[r][c] = 0.0;
      tdc[r][c] = 0.0;
    }
  }
  
  
  for(int m = 0; m < hcalt::ndata; m++) {
    r = hcalt::row[m]-1;
    c = hcalt::col[m]-1;
    if(r < 0 || c < 0) {
      std::cerr << "Why is row negative? Or col?" << std::endl;
      continue;
    }
    
    if(r>= kNrows || c >= kNcols) continue;

    //if(r==6||r==9) continue; //Temporary as fADC is in shop-------------------------------------------*
    
    idx = hcalt::samps_idx[m];
    n = hcalt::nsamps[m];
    adc[r][c] = hcalt::a[m];
    tdc[r][c] = hcalt::tdc[m];
    bool processed = false;
    for(int s = DISP_MIN_SAMPLE; s < DISP_MAX_SAMPLE && s < n; s++) {
      processed = true;
      histos[r][c]->SetBinContent(s+1-DISP_MIN_SAMPLE,hcalt::samps[idx+s]);
    }
  
    if(!processed) {
      cerr << "Skipping empty module: " << m << endl;
      for(int s = 0;  s < DISP_FADC_SAMPLES; s++) {
        histos[r][c]->SetBinContent(s+1,0);
      }
    }
  }
  
  for(r = 0; r < kNrows; r++) {
    for(c = 0; c < kNcols; c++) {
            
      if(tdc[r][c]==0){ //eliminating adc pulses from coincident tdc measurement in pedestal calculation
	
	for(int b=0 ; b<DISP_MAX_SAMPLE ; b++) {

	  //cout << 12*r+c << "  " << histos[r][c]->GetBinContent(b+1) << endl;
	  
	  pedSpec[r][c]->Fill(histos[r][c]->GetBinContent(b+1));
	}
      }
    }
  }
}

//Pedestal subtract then verify that max value is greater than twice sigma of pedestals and lies w/in bounds of histo
void goodHistoTest(TH1F *testHisto, double tdcVal, int row, int col){
  
  //Get pedestal value for the run and pmt 
  double pedVal = pedestals[row][col];
  
  //Subtract pedestal value from each bin
  for(int b=1 ; b<=testHisto->GetNbinsX() ; b++) {
    testHisto->SetBinContent(b,testHisto->GetBinContent(b)-pedVal);
  }
 
  if(testHisto->GetMaximum()>(6.0*pedSigma[row][col]) && testHisto->GetMaximumBin()>2 && testHisto->GetMaximumBin()<28){  //Cut out noise by ensuring all pulses are greater than twice the sigma of the ped and that the max lies in the center of the histogram
    
    if(tdcVal!=0) gPulseTDC[row+2][col+2]=true;
    
    gPulse[row+2][col+2]=true;
    
    if (signalTotal % 5000 == 0){
      cout << "Writing a reference histogram to file.." << endl;
      cout << "Pedestal for same histogram = " << pedestals[row][col] << "." << endl;
      cout << "Pedestal sigma for same histogram = " << pedSigma[row][col] << "." << endl;
      cout << "Maximum value for same histogram = " << testHisto->GetMaximum() << "." << endl;
      cout << "TDC value for same histogram = " << tdcVal << "." << endl << endl;
      
      HistosFile->cd();
      testHisto->SetTitle(Form("Sample R%d C%d",row,col));
      testHisto->GetYaxis()->SetTitle("RAU");
      testHisto->GetYaxis()->CenterTitle();
      testHisto->GetXaxis()->SetTitle("ADC Sample");
      testHisto->GetXaxis()->CenterTitle();
      testHisto->Write(Form("sampleHisto%d_r%d_c%d",signalTotal,row,col));
      testHisto->Draw("AP");
      DRAWNHISTO++;
    }
    
    signalTotal++;
    
  }
}

int cosCal_l2(int run = 989, int event = -1){

  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start(kTRUE);
  
  //Declare outfile
  TFile *cosAggFile = new TFile(Form("outFiles/cosmicHistograms_run%d.root",run),"RECREATE");
  
  //Build spectrum histograms. Empirical limits.
  cout << "Building spectrum histograms.." << endl;
  for(int r=0; r<kNrows; r++){
    for(int c=0; c<kNcols; c++){  
      PMTIntSpec[r][c] = new TH1F(Form("Int ADC Spect R%d C%d",r,c),Form("Int ADC Spect R%d C%d",r,c),100,0,4000);
      PMTIntSpec[r][c]->GetXaxis()->SetTitle("sRAU");
      PMTIntSpec[r][c]->GetXaxis()->CenterTitle();
      
      PMTMaxSpec[r][c] = new TH1F(Form("Max ADC Spect R%d C%d",r,c),Form("Max ADC Spect R%d C%d",r,c),100,0,1000);
      PMTMaxSpec[r][c]->GetXaxis()->SetTitle("RAU");
      PMTMaxSpec[r][c]->GetXaxis()->CenterTitle();
      
      PMTIntSpecTDC[r][c] = new TH1F(Form("Int TDC Spect R%d C%d",r,c),Form("Int TDC Spect R%d C%d",r,c),100,0,4000);
      PMTIntSpecTDC[r][c]->GetXaxis()->SetTitle("sRAU");
      PMTIntSpecTDC[r][c]->GetXaxis()->CenterTitle();
      
      PMTMaxSpecTDC[r][c] = new TH1F(Form("Max TDC Spect R%d C%d",r,c),Form("Max TDC Spect R%d C%d",r,c),100,0,1000);
      PMTMaxSpecTDC[r][c]->GetXaxis()->SetTitle("RAU");
      PMTMaxSpecTDC[r][c]->GetXaxis()->CenterTitle();
      
      pedSpec[r][c] = new TH1F(Form("Pedestal Spect R%d C%d",r,c),Form("Pedestal Spect R%d C%d",r,c),200,120,320);
      pedSpec[r][c]->GetXaxis()->SetTitle("<RAU>");
      pedSpec[r][c]->GetXaxis()->CenterTitle();

      PMTMaxPP[r][c] = new TH1F(Form("Max ADC Peak Position R%d C%d",r,c),Form("Max ADC Peak Position R%d C%d",r,c),50,0,50);
      PMTMaxPP[r][c]->GetXaxis()->SetTitle("ADC Sample");
      PMTMaxPP[r][c]->GetXaxis()->CenterTitle();
      
    }
  }

  //Build extra analysis histograms
  ADCvChannel = new TH1F("ADCvChannel", "ADCvChannel", kNcols*kNrows, 0, kNcols*kNrows-1);
  TDCvChannel = new TH1F("TDCvChannel", "TDCvChannel", kNcols*kNrows, 0, kNcols*kNrows-1);
  intADCvChannel = new TH1F("intADCvChannel", "intADCvChannel", kNcols*kNrows, 0, kNcols*kNrows-1);
  intTDCvChannel = new TH1F("intTDCvChannel", "intTDCvChannel", kNcols*kNrows, 0, kNcols*kNrows-1);
  NEVvChannel = new TH1F("NEVvChannel", "NEVvChannel", kNcols*kNrows, 0, kNcols*kNrows-1);
  pedvChannel = new TH1F("pedvChannel", "pedvChannel", kNcols*kNrows, 0, kNcols*kNrows-1);
  PPvCh = new TH1F("peakPosvChannel", "peakPosvChannel", kNcols*kNrows, 0, kNcols*kNrows-1);
  PPSDvCh = new TH1F("peakPosSDvChannel", "peakPosSDvChannel", kNcols*kNrows, 0, kNcols*kNrows-1);
    
  if(!T) { 
    T = new TChain("T");
    T->Add(TString::Format("rootFiles/cosmic/fadc_f1tdc_%d_*.root",run));
    /*
    T->Add(TString::Format("rootFiles/cosmic/fadc_f1tdc_%d_0.root",run));    
    T->Add(TString::Format("rootFiles/cosmic/fadc_f1tdc_%d_1.root",run));
    T->Add(TString::Format("rootFiles/cosmic/fadc_f1tdc_%d_2.root",run));
    T->Add(TString::Format("rootFiles/cosmic/fadc_f1tdc_%d_3.root",run));
    T->Add(TString::Format("rootFiles/cosmic/fadc_f1tdc_%d_4.root",run));
    T->Add(TString::Format("rootFiles/cosmic/fadc_f1tdc_%d_5.root",run));
    T->Add(TString::Format("rootFiles/cosmic/fadc_f1tdc_%d_6.root",run));
    T->Add(TString::Format("rootFiles/cosmic/fadc_f1tdc_%d_7.root",run));
    
    T->Add(TString::Format("rootFiles/cosmic/fadc_f1tdc_%d_8.root",run));
    T->Add(TString::Format("rootFiles/cosmic/fadc_f1tdc_%d_9.root",run));
    T->Add(TString::Format("rootFiles/cosmic/fadc_f1tdc_%d_10.root",run));
    T->Add(TString::Format("rootFiles/cosmic/fadc_f1tdc_%d_11.root",run));
    T->Add(TString::Format("rootFiles/cosmic/fadc_f1tdc_%d_12.root",run));
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
    cout << "Opened up tree with nentries=" << T->GetEntries() << endl;
    for(int r = 0; r < kNrows; r++) {
      for(int c = 0; c < kNcols; c++) {
        histos[r][c] = MakeHisto(r,c,DISP_FADC_SAMPLES);
        gSaturated[r][c] = false;
      }
    }
  }

  //Set appropriate HV and alphas for run. HV settings from HCAL wiki. Alphas from LED analysis. Must have accompanying text file, one double per line by module number. Assuming negative voltage inputs.

  ifstream file(Form("setFiles/HV_run%d.txt",run));

  cout << "Getting HV settings for each pmt for run " << run << "." << endl;

  int n1=0;
  double d1;
  
  int rval, cval;
  string line;
    
  while(getline(file,line)){
    if(line.at(0) == '#'){
      continue;
    }
    
    stringstream ss(line);
    ss >> d1;
    
    rval = floor(n1/kNcols);
    cval = n1 % kNcols;
    
    pmtHV[rval][cval] = -d1; 
    
    //cout << "row = " << rval << ", col = " << cval << ", HV val = " << pmtHV[rval][cval] << endl;
    
    n1++;
  }

  ifstream file2("setFiles/alphas.txt");
  
  n1=0;
  string line2;

  while(getline(file2,line2)){
    if(line2.at(0)=='#'){
      continue;
    }

    stringstream ss(line2);
    ss >> d1;

    rval = floor(n1/kNcols);
    cval = n1 % kNcols;

    alphas[rval][cval] = d1;
    
    //cout << "row = " << rval << ", col = " << cval << ", alphas val = " << alphas[rval][cval] << endl;

    n1++;
  }
  
  //Set default values for pulse check histograms. Arrays contain false-valued buffer on all sides
  for(int r = 0; r < kNrows+4; r++) {
    for(int c = 0; c < kNcols+4; c++) {
      gPulse[r][c] = false;
      gPulseTDC[r][c] = false;
    }
  }
  
  gCurrentEntry = event;

  cout << "Total events to process " << T->GetEntries() << endl;

  cout << "Filling pedestal histograms.." << endl;
  for (int i = gCurrentEntry; i < T->GetEntries(); i++){ //Not sure about GetEntries() here
    getPedestal(gCurrentEntry);
    gCurrentEntry++;
    
    //Keep count of events processed for monitoring
    if (gCurrentEntry%10000 == 0){
      cout << "Current pedestal event: " << gCurrentEntry << endl;
    }
  }
  
  for(int r=0; r<kNrows; r++){
    for(int c=0; c<kNcols; c++){  
      HistosFilePedestal->cd();
      pedSpec[r][c]->SetTitle(Form("Pedestal R%d C%d",r,c));
      pedSpec[r][c]->Write(Form("pedHisto_r%d_c%d",r,c));
      pedSpec[r][c]->Draw("AP");
    }
  }
  
  cout << "Fitting pedestal histos and extracting mean/sigma.." << endl;
  for(int r=0; r<kNrows; r++){
    for(int c=0; c<kNcols; c++){

      //cout << r << "  " << c << endl;
      
      //if(r==6||r==9) continue; //Temporary as fADC is in shop-------------------------------------------*

      pedSpec[r][c]->Fit("gaus","Q");
      f1=pedSpec[r][c]->GetFunction("gaus");
      pedestals[r][c]=f1->GetParameter(1);
      pedSigma[r][c]=f1->GetParameter(2);      
    }
  }

  gCurrentEntry = event; //Resetting entries for pulse analysis
  
  cout << "Processing pedestal subtracted pulses.." << endl;
  for (int i = gCurrentEntry; i < T->GetEntries(); i++){ //Not sure about GetEntries() here
    processEvent(gCurrentEntry);
    gCurrentEntry++;
    
    //Keep count of events processed for monitoring
    if (gCurrentEntry%10000 == 0){
      cout << "Current event: " << gCurrentEntry << endl;
    }
  }
  
  cout << "Finished loop over run " << run << "." << endl;
  cout << "Total good signals = " << signalTotal << "." << endl;

  cout << "Creating file to hold fit parameters.." << endl;
  ofstream outFile;
  outFile.open(Form("outFiles/cosCalOut_run%d.txt",run));
  outFile << "#Target HV settings for run " << run << "." << endl;
  outFile << "#Row Col targetHV targetHV_TDCCut" << endl;

  cout << "Writing spectrum histograms and calibration constants to file.." << endl;

  TF1 *landauFitMax[kNrows][kNcols] = {};
  TF1 *landauFitMaxTDC[kNrows][kNcols] = {};
  TF1 *landauFitInt[kNrows][kNcols] = {};
  TF1 *landauFitIntTDC[kNrows][kNcols] = {};

  double fitMaxX;
  double TDCfitMaxX;
  double fitIntX;
  double TDCfitIntX;
  
  for(int r=0; r<kNrows; r++){
    for(int c=0; c<kNcols; c++){

      //if(r==6||r==9) continue; //Temporary as fADC is in shop-------------------------------------------*
      
      //Create landau fit function per r and c
      landauFitMax[r][c] = new TF1(Form("landauFitMax r%d c%d",r,c), landauFit, 0, 1000, 4);
      landauFitMax[r][c]->SetLineColor(4);
      landauFitMax[r][c]->SetNpx(1000);
      landauFitMaxTDC[r][c] = new TF1(Form("landauFitMaxTDC r%d c%d",r,c), landauFit, 0, 1000, 4);
      landauFitMaxTDC[r][c]->SetLineColor(4);
      landauFitMaxTDC[r][c]->SetNpx(1000);
      landauFitInt[r][c] = new TF1(Form("landauFitInt r%d c%d",r,c), landauFit, 0, 4000, 4);
      landauFitInt[r][c]->SetLineColor(4);
      landauFitInt[r][c]->SetNpx(1000);
      landauFitIntTDC[r][c] = new TF1(Form("landauFitIntTDC r%d c%d",r,c), landauFit, 0, 4000, 4);
      landauFitIntTDC[r][c]->SetLineColor(4);
      landauFitIntTDC[r][c]->SetNpx(1000);

      //Set parameters for landau fit from spectrum histograms
      landauFitMax[r][c]->SetRange(PMTMaxSpec[r][c]->GetXaxis()->GetBinCenter(PMTMaxSpec[r][c]->GetMaximumBin())-0.8*PMTMaxSpec[r][c]->GetStdDev(),PMTMaxSpec[r][c]->GetXaxis()->GetBinCenter(PMTMaxSpec[r][c]->GetMaximumBin())+1.3*PMTMaxSpec[r][c]->GetStdDev());
      landauFitMax[r][c]->SetParameter(0,PMTMaxSpec[r][c]->GetBinContent(PMTMaxSpec[r][c]->GetMaximumBin()));

      //landauFitMax[r][c]->SetParameter(1,PMTMaxSpec[r][c]->GetMean()); //Try this to avoid fitting pedestal
      //landauFitMax[r][c]->SetParameter(1,PMTMaxSpec[r][c]->GetXaxis()->GetBinCenter(PMTMaxSpec[r][c]->GetMaximumBin()));
      landauFitMax[r][c]->SetParameter(1,PMTMaxSpec[r][c]->GetXaxis()->GetBinCenter(PMTMaxSpecTDC[r][c]->GetMaximumBin()));
      //landauFitMax[r][c]->SetParameter(1,PMTMaxSpecTDC[r][c]->GetMaximumBin());

      landauFitMax[r][c]->SetParameter(2,PMTMaxSpec[r][c]->GetStdDev());
      landauFitMax[r][c]->SetParameter(3,0);
      landauFitMaxTDC[r][c]->SetRange(PMTMaxSpecTDC[r][c]->GetXaxis()->GetBinCenter(PMTMaxSpecTDC[r][c]->GetMaximumBin())-0.8*PMTMaxSpecTDC[r][c]->GetStdDev(),PMTMaxSpecTDC[r][c]->GetXaxis()->GetBinCenter(PMTMaxSpecTDC[r][c]->GetMaximumBin())+1.3*PMTMaxSpecTDC[r][c]->GetStdDev());
      landauFitMaxTDC[r][c]->SetParameter(0,PMTMaxSpecTDC[r][c]->GetBinContent(PMTMaxSpecTDC[r][c]->GetMaximumBin()));
      landauFitMaxTDC[r][c]->SetParameter(1,PMTMaxSpecTDC[r][c]->GetMean()); //Try this to avoid fitting pedestal
      //landauFitMaxTDC[r][c]->SetParameter(1,PMTMaxSpecTDC[r][c]->GetXaxis()->GetBinCenter(PMTMaxSpecTDC[r][c]->GetMaximumBin()));

      landauFitMaxTDC[r][c]->SetParameter(2,PMTMaxSpecTDC[r][c]->GetStdDev());
      landauFitMaxTDC[r][c]->SetParameter(3,0);
      landauFitInt[r][c]->SetRange(PMTIntSpec[r][c]->GetXaxis()->GetBinCenter(PMTIntSpec[r][c]->GetMaximumBin())-0.8*PMTIntSpec[r][c]->GetStdDev(),PMTIntSpec[r][c]->GetXaxis()->GetBinCenter(PMTIntSpec[r][c]->GetMaximumBin())+1.3*PMTIntSpec[r][c]->GetStdDev());
      landauFitInt[r][c]->SetParameter(0,PMTIntSpec[r][c]->GetBinContent(PMTIntSpec[r][c]->GetMaximumBin()));


      //landauFitInt[r][c]->SetParameter(1,PMTIntSpec[r][c]->GetXaxis()->GetBinCenter(PMTIntSpec[r][c]->GetMaximumBin()));
      landauFitInt[r][c]->SetParameter(1,PMTIntSpec[r][c]->GetXaxis()->GetBinCenter(PMTIntSpecTDC[r][c]->GetMaximumBin()));
      //landauFitInt[r][c]->SetParameter(1,PMTIntSpecTDC[r][c]->GetMaximumBin());

      landauFitInt[r][c]->SetParameter(2,PMTIntSpec[r][c]->GetStdDev());
      landauFitInt[r][c]->SetParameter(3,0);
      landauFitIntTDC[r][c]->SetRange(PMTIntSpecTDC[r][c]->GetXaxis()->GetBinCenter(PMTIntSpecTDC[r][c]->GetMaximumBin())-0.8*PMTIntSpecTDC[r][c]->GetStdDev(),PMTIntSpecTDC[r][c]->GetXaxis()->GetBinCenter(PMTIntSpecTDC[r][c]->GetMaximumBin())+1.3*PMTIntSpecTDC[r][c]->GetStdDev());
      landauFitIntTDC[r][c]->SetParameter(0,PMTIntSpecTDC[r][c]->GetBinContent(PMTIntSpecTDC[r][c]->GetMaximumBin()));
      landauFitIntTDC[r][c]->SetParameter(1,PMTIntSpecTDC[r][c]->GetXaxis()->GetBinCenter(PMTIntSpecTDC[r][c]->GetMaximumBin()));
      landauFitIntTDC[r][c]->SetParameter(2,PMTIntSpecTDC[r][c]->GetStdDev());
      landauFitIntTDC[r][c]->SetParameter(3,0);

      pedSpec[r][c]->Fit("gaus","Q");
      f1=pedSpec[r][c]->GetFunction("gaus");
      
      cosAggFile->cd();
      pedSpec[r][c]->SetTitle(Form("Pedestal Spect R%d C%d",r,c));
      pedSpec[r][c]->GetXaxis()->SetTitle("RAU");
      pedSpec[r][c]->GetXaxis()->CenterTitle();
      pedSpec[r][c]->Write(Form("Pedestal Spect R%d C%d",r,c));
      pedSpec[r][c]->Draw("AP");
 
      //Spectra without TDC cut
      if( PMTIntSpec[r][c]->GetEntries() > 0){

	PMTIntSpec[r][c]->Fit(landauFitInt[r][c],"+RQ");
	if(landauFitInt[r][c]->GetParameter(1)>0&&landauFitInt[r][c]->GetParameter(1)<200){
	  fitIntX = landauFitInt[r][c]->GetParameter(1);
	  landauFitInt[r][c]->SetLineColor(kBlue);
	}else{
	  fitIntX = PMTIntSpec[r][c]->GetMean();
	  landauFitInt[r][c]->SetLineColor(kRed);
	}
	cosAggFile->cd();
	PMTIntSpec[r][c]->SetTitle(Form("Int ADC Spect R%d C%d MaxX%f",r,c,fitIntX));
	PMTIntSpec[r][c]->Write(Form("Int ADC Spect R%d C%d",r,c));
	PMTIntSpec[r][c]->Draw("AP");
      }
      
      if( PMTMaxSpec[r][c]->GetEntries() > 0){
	PMTMaxSpec[r][c]->Fit(landauFitMax[r][c],"+RQ");
	if(landauFitMax[r][c]->GetParameter(1)>0&&landauFitMax[r][c]->GetParameter(1)<2000){
	  fitMaxX = landauFitMax[r][c]->GetParameter(1);
	  landauFitMax[r][c]->SetLineColor(kBlue);
	}else{

	  fitMaxX = PMTMaxSpec[r][c]->GetMean();
	  landauFitMax[r][c]->SetLineColor(kRed);
	}
	cosAggFile->cd();
	PMTMaxSpec[r][c]->SetTitle(Form("Max ADC Spect R%d C%d MaxX%f",r,c,fitMaxX));
	PMTMaxSpec[r][c]->Write(Form("Max ADC Spect R%d C%d",r,c));
	PMTMaxSpec[r][c]->Draw("AP");
      }

      if(PMTMaxPP[r][c]->GetEntries()>0){
	PMTMaxPP[r][c]->Fit("gaus","Q");
	f3=PMTMaxPP[r][c]->GetFunction("gaus");
	PMTMaxPP[r][c]->SetTitle(Form("Max ADC Peak Pos R%d C%d Mean%f SD%f",r,c,f3->GetParameter(1),f3->GetParameter(2)));
	PPmean[r][c] = PMTMaxPP[r][c]->GetMean();
	PPsd[r][c] = PMTMaxPP[r][c]->GetRMS();
	cosAggFile->cd();
	PMTMaxPP[r][c]->Write(Form("PP R%d C%d",r,c));
	PMTMaxPP[r][c]->Draw("AP");
      }
	            
      targetHV[r][c] = pmtHV[r][c]/pow(fitMaxX/targetRAU,1.0/alphas[r][c]);
      
      //cout << targetHV[r][c] << endl;

      //Spectra with TDC cut
      if( PMTIntSpecTDC[r][c]->GetEntries() > 0){
	PMTIntSpecTDC[r][c]->Fit(landauFitIntTDC[r][c],"+RQ");
	if(landauFitIntTDC[r][c]->GetParameter(1)>0&&landauFitIntTDC[r][c]->GetParameter(1)<2000){
	  TDCfitIntX = landauFitIntTDC[r][c]->GetParameter(1);
	  landauFitIntTDC[r][c]->SetLineColor(kBlue);
	}else{
	  TDCfitIntX = PMTIntSpecTDC[r][c]->GetMean();
	  landauFitIntTDC[r][c]->SetLineColor(kRed);
	}
	cosAggFile->cd();
	PMTIntSpecTDC[r][c]->SetTitle(Form("Int TDC Spect R%d C%d MaxX%f",r,c,TDCfitIntX));
	PMTIntSpecTDC[r][c]->Write(Form("Int TDC Spect R%d C%d",r,c));
	PMTIntSpecTDC[r][c]->Draw("AP");
      }
      
      if( PMTMaxSpecTDC[r][c]->GetEntries() > 0){
	PMTMaxSpecTDC[r][c]->Fit(landauFitMaxTDC[r][c],"+RQ");
	if(landauFitMaxTDC[r][c]->GetParameter(1)>0&&landauFitMaxTDC[r][c]->GetParameter(1)<400){
	  TDCfitMaxX = landauFitMaxTDC[r][c]->GetParameter(1);
	  landauFitMaxTDC[r][c]->SetLineColor(kBlue);
	}else{
	  TDCfitMaxX = PMTMaxSpecTDC[r][c]->GetMean();
	  landauFitMaxTDC[r][c]->SetLineColor(kRed);
	}
	if(TDCfitMaxX>2000)cout << "TDCfitMaxX>2000 on " << r << "  " << c << endl;
	cosAggFile->cd();
	PMTMaxSpecTDC[r][c]->SetTitle(Form("Max TDC Spect R%d C%d MaxX%f",r,c,TDCfitMaxX));
	PMTMaxSpecTDC[r][c]->Write(Form("Max TDC Spect R%d C%d",r,c));
	PMTMaxSpecTDC[r][c]->Draw("AP");
      }

      targetHVTDC[r][c] = pmtHV[r][c]/pow(TDCfitMaxX/targetRAU,1.0/alphas[r][c]);
      
      //Fill Max ADC vs Channel Histogram
      ADCvChannel->SetBinContent(kNcols*r+c+1,fitMaxX);
      if(fitMaxX > 200) cout << "Very high average ADC value on r " << r << ", c " << c << "." << endl;
      TDCvChannel->SetBinContent(kNcols*r+c+1,TDCfitMaxX);      
      intADCvChannel->SetBinContent(kNcols*r+c+1,fitIntX);
      intTDCvChannel->SetBinContent(kNcols*r+c+1,TDCfitIntX);

      //Fill Number of events per module histogram
      NEVvChannel->SetBinContent(kNcols*r+c+1,PMTMaxSpec[r][c]->GetEntries());

      //Fill average pedestal value per module histograms
      pedvChannel->SetBinContent(kNcols*r+c+1,pedestals[r][c]);

      //Fill peak position and peak position sigma per module histograms
      PPvCh->SetBinContent(kNcols*r+c+1,PPmean[r][c]);
      PPSDvCh->SetBinContent(kNcols*r+c+1,PPsd[r][c]);

      cout << "For row " << r << " and col " << c << ", target HV with landau fit is " << targetHV[r][c] << ", and target HV (TDC cut) with landau fit is " << targetHVTDC[r][c] << "." << endl;
      outFile << r << "  " << c << "  " << targetHV[r][c] << "  " << targetHVTDC[r][c] << "  " << endl;
      
    }
  }

  //Draw analysis histograms

  cosAggFile->cd();
  ADCvChannel->SetTitle("Average Max ADC Val vs Channel");
  ADCvChannel->Write("ADCvChannel");
  ADCvChannel->Draw("AP");

  cosAggFile->cd();
  TDCvChannel->SetTitle("Average Max ADC Val with TDC Cut vs Channel");
  TDCvChannel->Write("TDCvChannel");
  TDCvChannel->Draw("AP");

  cosAggFile->cd();
  intADCvChannel->SetTitle("Average Int ADC Val vs Channel");
  intADCvChannel->Write("intADCvChannel");
  intADCvChannel->Draw("AP");

  cosAggFile->cd();
  intTDCvChannel->SetTitle("Average Int ADC Val with TDC Cut vs Channel");
  intTDCvChannel->Write("intTDCvChannel");
  intTDCvChannel->Draw("AP");
  
  cosAggFile->cd();
  PPvCh->SetTitle("Average ADC Sample Peak Position vs Channel");
  PPvCh->Write("PPvCh");
  PPvCh->Draw("AP");

  cosAggFile->cd();
  PPSDvCh->SetTitle("Average ADC Sample Peak Pos Stand Dev vs Channel");
  PPSDvCh->Write("PPSDvCh");
  PPSDvCh->Draw("AP");
  
  cosAggFile->cd();
  NEVvChannel->SetTitle("Number of Event Pulses vs Channel");
  NEVvChannel->Write("NEVvChannel");
  NEVvChannel->Draw("AP");

  cosAggFile->cd();
  pedvChannel->SetTitle("Avg Pedestal vs Channel");
  pedvChannel->Write("pedvChannel");
  pedvChannel->Draw("AP");

  outFile << endl << endl;
  
  for(int r=0; r<kNrows; r++){
    for(int c=0; c<kNcols; c++){
      
      outFile << "Pedestal for row " << r << " and col " << c << " = " << pedestals[r][c] << endl;

    }
  }
  
  cout << "ADC spectrum histograms written to file cosSpectFile_run" << run << ".root" << endl;
  cout << "Sample histograms drawn to file HistosFile.root" << endl;
  cout << "Target HV settings written to cosCalOut_run" << run << ".txt" << endl;

  st->Stop();

  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

  
  return 0;
}
  
