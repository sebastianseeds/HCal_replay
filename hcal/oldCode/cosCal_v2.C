#include <TH2F.h>
#include <TChain.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <TSystem.h>
#include "hcal.h"

const int kNrows = 12;
const int kNcols = 12;
const float HV[2]={1500,1600};
const int HVpts = 1;

const int DISP_MIN_SAMPLE = 0.0;
const int DISP_MAX_SAMPLE = 20.0;
const int DISP_FADC_SAMPLES = (DISP_MAX_SAMPLE-DISP_MIN_SAMPLE);

int gCurrentEntry = -1;
TChain *T = 0;

//Augment fn's for cosmic Calibration
void goodHistoTest(TH1F*, double, int, int);
void getPedestal(int);

TH1F *histos[kNrows][kNcols];
TH1F *testHistos[kNrows][kNcols];
TH1F *pedSpec[kNrows][kNcols];
TH1F *PMTIntSpec[kNrows][kNcols]; //=new TH1F("","",20,0,5000); //Empirical limits
TH1F *PMTMaxSpec[kNrows][kNcols]; //=new TH1F("","",20,0,1000); //Empirical limits
TH1F *PMTIntSpecTDC[kNrows][kNcols]; //=new TH1F("","",20,0,5000); //Empirical limits, cut on tdc coincidence
TH1F *PMTMaxSpecTDC[kNrows][kNcols]; //=new TH1F("","",20,0,1000); //Empirical limits, cut on tdc coincidence
TF1 *f1;
TF1 *f2;
TF1 *f3;
TF1 *f4;
TF1 *f5;

int signalTotal = 1;
TFile *HistosFile = new TFile("HistosFile.root","RECREATE");  //File for checking histogram fits and integrity
double pedestals[kNrows][kNcols];
double pedSigma[kNrows][kNcols];
bool gSaturated[kNrows][kNcols];
bool gPulse[kNrows+2][kNcols+2]; //Needs to be larger for false-valued buffers
bool gPulseTDC[kNrows+2][kNcols+2]; //Needs to be larger for false-valued buffers
int DRAWNHISTO = 0; //Counts up per histo drawn
int rejectNote = 0;
int passNote = 0;
int rejectNoteTDC = 0;
int passNoteTDC = 0;
int histNote = 0;
int r0c0check = 0;

double targetRAU = 61.425; //Taken from Scott B.
int runN[9] = {989,988,987,984,983,981,980,978,1235}; //All run numbers
double jlab_lHV[9] = {1600.0,1650.0,1750.0,1850.0,2000.0,1900.0,1800.0,1700.0,1600.0}; //For all runs, left half jlab pmt HV
double cmu_lHV[9] = {1500.0,1550.0,1650.0,1750.0,1900.0,1800.0,1700.0,1600.0,1500.0}; //For all runs, left half cmu pmt HV
double jlab_rHV = 1600.0;
double cmu_rHV = 1500.0;
double jlab_lHVset, cmu_lHVset;
double cmuExp = 10.5; //Exp parameter per gain equation on data sheet
double jlabExp = 8.0; //Exp parameter per gain equation on data sheet
double targetHV[kNrows][kNcols];
double targetHVTDC[kNrows][kNcols];

double powFit(double *X,double *p) 
{
  double fitval = p[0]*pow(X[0],p[1]);
  //double fitval = p[0]*pow((X[0]-p[1])/1200,p[2]) + p[3]; //X normalized to lowest HV value
  return fitval;
}

//Create skewed normal distribution to fit the timing resolution.
double skewFit(double *X,double *p) 
{
  double fitval = p[0]*exp(-0.5*pow((X[0]-p[1])/p[2],2.0))*0.5*(1+erf(p[3]*((X[0]-p[1])/p[2])));
  return fitval;
}


TH1F* MakeHisto(int row, int col, int bins){
  TH1F *h = new TH1F(TString::Format("h%02d%02d",row,col),
      TString::Format("%d-%d",row+1,col+1),bins,DISP_MIN_SAMPLE,DISP_MAX_SAMPLE);
  h->SetStats(0);
  h->SetLineWidth(2);
  return h;
}

//Rows and Columns start at zero and go to kNrows-1 and kNcols-1
void displayEvent(int entry = -1){
  if(entry == -1) {
    gCurrentEntry++;
  } else {
    gCurrentEntry = entry;
  }

  if(gCurrentEntry<0) {
    gCurrentEntry = 0;
  }

  T->GetEntry(gCurrentEntry);

  //histNote=0;
  
  int r,c,idx,n,sub;
  // Clear old histograms, just in case modules are not in the tree
  for(r = 0; r < kNrows; r++) {
    for(c = 0; c < kNcols; c++) {
      histos[r][c]->Reset("ICES M");
      gSaturated[r][c] = false;
      gPulse[r+1][c+1] = false;
      gPulseTDC[r+1][c+1] = false;
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
    
    if(r>= kNrows || c >= kNcols)
      continue;
    idx = hcalt::samps_idx[m];
    n = hcalt::nsamps[m];
    adc[r][c] = hcalt::a[m];
    tdc[r][c] = hcalt::tdc[m];
    bool displayed = false;
    for(int s = DISP_MIN_SAMPLE; s < DISP_MAX_SAMPLE && s < n; s++) {
      displayed = true;
      histos[r][c]->SetBinContent(s+1-DISP_MIN_SAMPLE,hcalt::samps[idx+s]);
      //testHistos[r][c]->SetBinContent(s+1-DISP_MIN_SAMPLE,hcalt::samps[idx+s]);
      if(peak[r][c]<hcalt::samps[idx+s])
        peak[r][c]=hcalt::samps[idx+s];
      if(peak[r][c]>4095) {
        gSaturated[r][c] = true;
      }
    }
    if(!displayed) {
      std::cerr << "Skipping empty module: " << m << std::endl;
      for(int s = 0;  s < DISP_FADC_SAMPLES; s++) {
        histos[r][c]->SetBinContent(s+1,-404);
        //testHistos[r][c]->SetBinContent(s+1,-404);
      }
    }
  }
  
  for(r = 0; r < kNrows; r++) {
    for(c = 0; c < kNcols; c++) {
      histos[r][c]->SetTitle(TString::Format("%d-%d (ADC=%g,TDC=%g)",r+1,c+1,adc[r][c],tdc[r][c]));
      if(!gSaturated[r][c]){
	//testHistos[r][c]=histos[r][c];
	//goodHistoTest(testHistos[r][c],tdc[r][c],r,c);
	goodHistoTest(histos[r][c],tdc[r][c],r,c);
	//checkHistos++;
      }
    }
  }

  //Vertical line test - remember to shift all values up by one due to intialization of gPulse/gPulseTDC and 'false' buffer
  for(r = 0; r < kNrows; r++) {
    for(c = 0; c < kNcols; c++) {
      if(gPulse[r+1][c+1]==true){
	if(gPulse[r+2][c]==false&&gPulse[r+2][c+1]==false&&gPulse[r+2][c+2]==false){
	  if(gPulse[r][c]==false&&gPulse[r][c+1]==false&&gPulse[r][c+2]==false){
	    gPulse[r+1][c+1]=false;

	    if(r==0&&c==0&&tdc[r][c]!=0) cout << histos[r][c]->GetBinContent(histos[r][c]->GetMaximumBin()) << "!!!!!!!" << endl;
	    
	    if(rejectNote<10) cout << "Pulse r" << r << " c" << c << " REJECTED on verticality cut." << endl;
	    rejectNote++;
	  }
	}
      }
      if(gPulse[r+1][c+1]==true&&passNote<10){
	cout << "Pulse r" << r << " c" << c << " passes verticality cut." << endl;
	passNote++;
      }
      if(gPulseTDC[r+1][c+1]==true){
	if(gPulseTDC[r+2][c]==false&&gPulseTDC[r+2][c+1]==false&&gPulseTDC[r+2][c+2]==false){
	  if(gPulseTDC[r][c]==false&&gPulseTDC[r][c+1]==false&&gPulseTDC[r][c+2]==false){
	    gPulseTDC[r+1][c+1]=false;
	    if(rejectNoteTDC<10) cout << "TDC cut pulse r" << r << " c" << c << " REJECTED on verticality cut." << endl;
	    rejectNoteTDC++;
	  }
	}
      }
      if(gPulse[r+1][c+1]==true&&passNoteTDC<10){
	cout << "TDC cut pulse r" << r << " c" << c << " passes verticality cut." << endl;
	passNoteTDC++;
      }
    }
  }
  
  //Now, if pulse passes verticality test, pedestal subtract and fill spectra histograms
  for(r = 0; r < kNrows; r++) {
    for(c = 0; c < kNcols; c++) {

      if(gPulse[r+1][c+1]==true&&r==0&&c==0) cout << "Pulse on PMT r0 c0!!!" << endl;

      if(gPulse[r+1][c+1]==true){

	//double pedVal = pedestals[r][c];
  
	//Subtract pedestal value from each bin
	//for(int b=1 ; b<=histos[r][c]->GetNbinsX() ; b++) {
	//histos[r][c]->SetBinContent(b,histos[r][c]->GetBinContent(b)-pedVal);
	//}
	
	//maxy = peak[r][c]; //WRONG, not pedestal subtracted

	if(r==0&&c==0&&histNote<100){
	  cout << "Int value added to histogram = " << histos[r][c]->Integral(1,DISP_MAX_SAMPLE) << "." << endl;
	  cout << "Max value added to histogram = " << histos[r][c]->GetMaximum() << "." << endl;
	  cout << "Max value added to histogram, method 2,  = " << histos[r][c]->GetBinContent(histos[r][c]->GetMaximumBin()) << "." << endl;
	  histNote++;
	}
	
	PMTIntSpec[r][c]->Fill(histos[r][c]->Integral(1,DISP_MAX_SAMPLE)); //Integral from total bin content
	PMTMaxSpec[r][c]->Fill(histos[r][c]->GetMaximum());
	//cout << PMTIntSpec[r][c]->GetEntries() << endl;
	//cout << PMTMaxSpec[r][c]->GetEntries() << endl;
	/*
	if (signalTotal % 10000 == 0){
	  cout << "Writing a reference histogram to file.." << endl;
	  cout << "Pedestal for same histogram = " << pedestals[r][c] << "." << endl;
	  cout << "Pedestal sigma for same histogram = " << pedSigma[r][c] << "." << endl;
	  cout << "Maximum value for same histogram = " << histos[r][c]->GetMaximum() << "." << endl;
	  cout << "TDC value for same histogram = " << tdcVal << "." << endl << endl;
	  
	  HistosFile->cd();
	  testHisto->SetTitle(Form("Sample R%d C%d",r,c));
	  testHisto->Write(Form("sampleHisto%d_r%d_c%d",signalTotal,r,c));
	  testHisto->Draw("AP");
	  DRAWNHISTO++;
	}
	
	signalTotal++;
	*/
	if(gPulseTDC[r+1][c+1]==true){
	  PMTIntSpecTDC[r][c]->Fill(histos[r][c]->Integral(1,DISP_MAX_SAMPLE)); //Integral from total bin content
	  PMTMaxSpecTDC[r][c]->Fill(histos[r][c]->GetMaximum());
	}	
      }	
    }
  }
}

//Open pedestal text file and read out all pedestal values by HV run
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
  for(r  = 0; r < kNrows; r++) {
    for(c  = 0; c < kNcols; c++) {
      adc[r][c] = 0.0;
    }
  }
  
  
  for(int m = 0; m < hcalt::ndata; m++) {
    r = hcalt::row[m]-1;
    c = hcalt::col[m]-1;
    if(r < 0 || c < 0) {
      std::cerr << "Why is row negative? Or col?" << std::endl;
      continue;
    }
    
    if(r>= kNrows || c >= kNcols)
      continue;
    idx = hcalt::samps_idx[m];
    n = hcalt::nsamps[m];
    adc[r][c] = hcalt::a[m];
    bool displayed = false;
    for(int s = DISP_MIN_SAMPLE; s < DISP_MAX_SAMPLE && s < n; s++) {
      displayed = true;
      histos[r][c]->SetBinContent(s+1-DISP_MIN_SAMPLE,hcalt::samps[idx+s]);
    }
  
    if(!displayed) {
      std::cerr << "Skipping empty module: " << m << std::endl;
      for(int s = 0;  s < DISP_FADC_SAMPLES; s++) {
        histos[r][c]->SetBinContent(s+1,0);
      }
    }
  }
  
  for(r = 0; r < kNrows; r++) {
    for(c = 0; c < kNcols; c++) {
      
      //pedSpec[r][c]->Fill(histos[r][c]->Integral()/DISP_MAX_SAMPLE);

      double sum = 0.0;
      
      //Subtract pedestal value from each bin
      for(int b=1 ; b<5 ; b++) {
      sum+=histos[r][c]->GetBinContent(b);
      }

      pedSpec[r][c]->Fill(sum/4);
      
    }
  }
}

//Pedestal subtract then verify that max value is greater than twice sigma of pedestals and lies w/in bounds of histo
void goodHistoTest(TH1F *testHisto, double tdcVal, int row, int col){
  
  //Get pedestal value for the run and pmt 
  double pedVal = pedestals[row][col];
  
  //Subtract pedestal value from each bin and cut off bin 20 (overflow)
  for(int b=1 ; b<=testHisto->GetNbinsX() ; b++) {
    testHisto->SetBinContent(b,testHisto->GetBinContent(b)-pedVal);
    //if(testHisto->GetBinContent(b) < 0) testHisto->SetBinContent(b,0.0); //romoved to check pedestals
  }

  //if(row==0&&col==0&&tdcVal!=0) cout << " PEDSIGMA " << pedSigma[row][col] << ", MAX " << testHisto->GetMaximum() << ", ..WHAT THE HELL?" << endl;

  if (r0c0check < 10 && row==0 && col==0 && tdcVal!=0){
    HistosFile->cd();
    testHisto->SetTitle(Form("Sample R%d C%d",row,col));
    testHisto->Write(Form("r0c0sampleHisto%d_r%d_c%d",r0c0check,row,col));
    testHisto->Draw("AP");
    r0c0check++;
  }
  
  if(testHisto->GetMaximum()>50.0 && testHisto->GetMaximum()>(3.0*pedSigma[row][col]) && testHisto->GetMaximumBin()>3 && testHisto->GetMaximumBin()<17){  //Cut out noise by ensuring all pulses are greater than twice the sigma of the ped and that the max lies in the center of the histogram
    
    if(tdcVal!=0) gPulseTDC[row+1][col+1]=true;

    if(row==0&&col==0) cout << "IS GOING ON?" << endl;
    
    gPulse[row+1][col+1]=true;
    
    if (signalTotal % 10000 == 0){
      cout << "Writing a reference histogram to file.." << endl;
      cout << "Pedestal for same histogram = " << pedestals[row][col] << "." << endl;
      cout << "Pedestal sigma for same histogram = " << pedSigma[row][col] << "." << endl;
      cout << "Maximum value for same histogram = " << testHisto->GetMaximum() << "." << endl;
      cout << "TDC value for same histogram = " << tdcVal << "." << endl << endl;
      
      HistosFile->cd();
      testHisto->SetTitle(Form("Sample R%d C%d",row,col));
      testHisto->Write(Form("sampleHisto%d_r%d_c%d",signalTotal,row,col));
      testHisto->Draw("AP");
      DRAWNHISTO++;
    }
    
    signalTotal++;
    
  }
}

int cosCal_v2(int run = 989, int event = -1)
{

  //Declare outfile
  TFile *cosAggFile = new TFile(Form("cosSpectFile_run%d.root",run),"RECREATE");
  
  //Build spectrum histograms. Empirical limits.
  cout << "Building spectrum histograms.." << endl;
  for(int r=0; r<kNrows; r++){
    for(int c=0; c<kNcols; c++){  
      PMTIntSpec[r][c] = new TH1F(Form("Int ADC Spect R%d C%d",r,c),Form("Int ADC Spect R%d C%d",r,c),50,0,4000);
      PMTIntSpec[r][c]->GetXaxis()->SetTitle("sRAU");
      PMTIntSpec[r][c]->GetXaxis()->CenterTitle();
      PMTMaxSpec[r][c] = new TH1F(Form("Max ADC Spect R%d C%d",r,c),Form("Max ADC Spect R%d C%d",r,c),50,0,1000);
      PMTMaxSpec[r][c]->GetXaxis()->SetTitle("RAU");
      PMTMaxSpec[r][c]->GetXaxis()->CenterTitle();
      PMTIntSpecTDC[r][c] = new TH1F(Form("Int ADC Spect R%d C%d, TDC Cut",r,c),Form("Int ADC Spect R%d C%d",r,c),50,0,4000);
      PMTIntSpecTDC[r][c]->GetXaxis()->SetTitle("sRAU");
      PMTIntSpecTDC[r][c]->GetXaxis()->CenterTitle();
      PMTMaxSpecTDC[r][c] = new TH1F(Form("Max ADC Spect R%d C%d, TDC Cut",r,c),Form("Max ADC Spect R%d C%d",r,c),50,0,1000);
      PMTMaxSpecTDC[r][c]->GetXaxis()->SetTitle("RAU");
      PMTMaxSpecTDC[r][c]->GetXaxis()->CenterTitle();
      pedSpec[r][c] = new TH1F(Form("Pedestal Spect R%d C%d",r,c),Form("Pedestal Spect R%d C%d",r,c),50,0,0);
      pedSpec[r][c]->GetXaxis()->SetTitle("<RAU>");
      pedSpec[r][c]->GetXaxis()->CenterTitle();
    }
  }
    
  if(!T) { 
    T = new TChain("T");
    T->Add(TString::Format("fadc_f1tdc_%d.root",run));
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
	//testHistos[r][c] = MakeHisto(r,c,DISP_FADC_SAMPLES);
        gSaturated[r][c] = false;
      }
    }
  }

  //Set appropriate HVs for run
  for(int i=0; i<9; i++){
    if(run==runN[i]){
      jlab_lHVset = jlab_lHV[i];
      cmu_lHVset = cmu_lHV[i];
    }
  }
  
  //Set default values for pulse check histograms. Arrays contain false-valued buffer on all sides
  for(int r = 0; r < kNrows+2; r++) {
    for(int c = 0; c < kNcols+2; c++) {
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

  cout << "Fitting pedestal histos and extracting mean/sigma.." << endl;
  for(int r=0; r<kNrows; r++){
    for(int c=0; c<kNcols; c++){  
      pedSpec[r][c]->Fit("gaus","Q");
      f1=pedSpec[r][c]->GetFunction("gaus");
      pedestals[r][c]=f1->GetParameter(1);
      pedSigma[r][c]=f1->GetParameter(2);
    }
  }

  gCurrentEntry = event; //Resetting entries for pulse analysis
  
  cout << "Processing pedestal subtracted pulses.." << endl;
  for (int i = gCurrentEntry; i < T->GetEntries(); i++){ //Not sure about GetEntries() here
    displayEvent(gCurrentEntry);
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
  outFile.open(Form("cosCalOut_run%d.txt",run));
  outFile << "#Target HV settings for run " << run << "." << endl;
  outFile << "#Row Col targetHV targetHV_TDCCut" << endl;

  cout << "Writing spectrum histograms and calibration constants to file.." << endl;

  //Now to fit these spectra
  //TF1 **skewFitMax = new TF1*[kNrows][kNcols];
  //TF1 **skewFitMaxTDC = new TF1*[kNrows][kNcols];
  //TF1 **skewFitInt = new TF1*[kNrows][kNcols];
  //TF1 **skewFitIntTDC = new TF1*[kNrows][kNcols];

  TF1 *skewFitMax[kNrows][kNcols] = {};
  TF1 *skewFitMaxTDC[kNrows][kNcols] = {};
  TF1 *skewFitInt[kNrows][kNcols] = {};
  TF1 *skewFitIntTDC[kNrows][kNcols] = {};
  
  double fitMaxX;
  double TDCfitMaxX;
  double fitIntX;
  double TDCfitIntX;
  
  for(int r=0; r<kNrows; r++){
    for(int c=0; c<kNcols; c++){

      skewFitMax[r][c] = new TF1();
      skewFitMaxTDC[r][c] = new TF1();
      skewFitInt[r][c] = new TF1();
      skewFitIntTDC[r][c] = new TF1();
      
      //Create skew fit function per r and c
      skewFitMax[r][c] = new TF1(Form("SkewFitMax r%d c%d",r,c),skewFit, 0, 1000, 4);
      skewFitMax[r][c]->SetLineColor(4);
      skewFitMax[r][c]->SetNpx(1000);
      skewFitMaxTDC[r][c] = new TF1(Form("SkewFitMaxTDC r%d c%d",r,c),skewFit, 0, 1000, 4);
      skewFitMaxTDC[r][c]->SetLineColor(4);
      skewFitMaxTDC[r][c]->SetNpx(1000);
      skewFitInt[r][c] = new TF1(Form("SkewFitInt r%d c%d",r,c),skewFit, 0, 4000, 4);
      skewFitInt[r][c]->SetLineColor(4);
      skewFitInt[r][c]->SetNpx(1000);
      skewFitIntTDC[r][c] = new TF1(Form("SkewFitIntTDC r%d c%d",r,c),skewFit, 0, 4000, 4);
      skewFitIntTDC[r][c]->SetLineColor(4);
      skewFitIntTDC[r][c]->SetNpx(1000);

      //Set parameters for function from spectrum histograms
      skewFitMax[r][c]->SetRange(PMTMaxSpec[r][c]->GetXaxis()->GetBinCenter(PMTMaxSpec[r][c]->GetMaximumBin())-0.8*PMTMaxSpec[r][c]->GetStdDev(),PMTMaxSpec[r][c]->GetXaxis()->GetBinCenter(PMTMaxSpec[r][c]->GetMaximumBin())+1.3*PMTMaxSpec[r][c]->GetStdDev());
      skewFitMax[r][c]->SetParameter(0,PMTMaxSpec[r][c]->GetBinContent(PMTMaxSpec[r][c]->GetMaximumBin()));;
      skewFitMax[r][c]->SetParameter(1,PMTMaxSpec[r][c]->GetXaxis()->GetBinCenter(PMTMaxSpec[r][c]->GetMaximumBin()));
      skewFitMax[r][c]->SetParameter(2,PMTMaxSpec[r][c]->GetStdDev());
      skewFitMax[r][c]->SetParameter(3,2);
      skewFitMaxTDC[r][c]->SetRange(PMTMaxSpecTDC[r][c]->GetXaxis()->GetBinCenter(PMTMaxSpecTDC[r][c]->GetMaximumBin())-0.8*PMTMaxSpecTDC[r][c]->GetStdDev(),PMTMaxSpecTDC[r][c]->GetXaxis()->GetBinCenter(PMTMaxSpecTDC[r][c]->GetMaximumBin())+1.3*PMTMaxSpecTDC[r][c]->GetStdDev());
      skewFitMaxTDC[r][c]->SetParameter(0,PMTMaxSpecTDC[r][c]->GetBinContent(PMTMaxSpecTDC[r][c]->GetMaximumBin()));;
      skewFitMaxTDC[r][c]->SetParameter(1,PMTMaxSpecTDC[r][c]->GetXaxis()->GetBinCenter(PMTMaxSpecTDC[r][c]->GetMaximumBin()));
      skewFitMaxTDC[r][c]->SetParameter(2,PMTMaxSpecTDC[r][c]->GetStdDev());
      skewFitMaxTDC[r][c]->SetParameter(3,2);
      skewFitInt[r][c]->SetRange(PMTIntSpec[r][c]->GetXaxis()->GetBinCenter(PMTIntSpec[r][c]->GetMaximumBin())-0.5*PMTIntSpec[r][c]->GetStdDev(),PMTIntSpec[r][c]->GetXaxis()->GetBinCenter(PMTIntSpec[r][c]->GetMaximumBin())+1.3*PMTIntSpec[r][c]->GetStdDev());
      skewFitInt[r][c]->SetParameter(0,PMTIntSpec[r][c]->GetBinContent(PMTIntSpec[r][c]->GetMaximumBin()));;
      skewFitInt[r][c]->SetParameter(1,PMTIntSpec[r][c]->GetXaxis()->GetBinCenter(PMTIntSpec[r][c]->GetMaximumBin()));
      skewFitInt[r][c]->SetParameter(2,PMTIntSpec[r][c]->GetStdDev());
      skewFitInt[r][c]->SetParameter(3,2);
      skewFitIntTDC[r][c]->SetRange(PMTIntSpecTDC[r][c]->GetXaxis()->GetBinCenter(PMTIntSpecTDC[r][c]->GetMaximumBin())-0.5*PMTIntSpecTDC[r][c]->GetStdDev(),PMTIntSpecTDC[r][c]->GetXaxis()->GetBinCenter(PMTIntSpecTDC[r][c]->GetMaximumBin())+1.3*PMTIntSpecTDC[r][c]->GetStdDev());
      skewFitIntTDC[r][c]->SetParameter(0,PMTIntSpecTDC[r][c]->GetBinContent(PMTIntSpecTDC[r][c]->GetMaximumBin()));;
      skewFitIntTDC[r][c]->SetParameter(1,PMTIntSpecTDC[r][c]->GetXaxis()->GetBinCenter(PMTIntSpecTDC[r][c]->GetMaximumBin()));
      skewFitIntTDC[r][c]->SetParameter(2,PMTIntSpecTDC[r][c]->GetStdDev());
      skewFitIntTDC[r][c]->SetParameter(3,2);
      
      int iHV=0; 
      if(c>3&&c<8) iHV=1;

      pedSpec[r][c]->Fit("gaus","Q");
      f1=pedSpec[r][c]->GetFunction("gaus");
      
      cosAggFile->cd();
      pedSpec[r][c]->Write(Form("Pedestal Spect R%d C%d",r,c));
      pedSpec[r][c]->Draw("AP");


      //Spectra without TDC cut
      if( PMTIntSpec[r][c]->GetEntries() > 0){
      cout << "PMTIntSpec entries for row " << r << ", col " << c << " = " << PMTIntSpec[r][c]->GetEntries() << "." << endl;
      //PMTIntSpec[r][c]->Fit("gaus","Q");
      PMTIntSpec[r][c]->Fit(skewFitInt[r][c],"RQ");
      //f2=PMTIntSpec[r][c]->GetFunction("gaus");
      fitIntX = skewFitInt[r][c]->GetMaximumX();
      cosAggFile->cd();
      PMTIntSpec[r][c]->SetTitle(Form("Int ADC Spect R%d C%d MaxX%f",r,c,fitIntX));
      PMTIntSpec[r][c]->Write(Form("Int ADC Spect R%d C%d",r,c));
      PMTIntSpec[r][c]->Draw("AP");
      }
      
      if( PMTMaxSpec[r][c]->GetEntries() > 0){
      cout << "PMTMaxSpec entries for row " << r << ", col " << c << " = " << PMTMaxSpec[r][c]->GetEntries() << "." << endl;
      //PMTMaxSpec[r][c]->Fit("gaus","Q");
      PMTMaxSpec[r][c]->Fit(skewFitMax[r][c],"RQ");
      //f3=PMTMaxSpec[r][c]->GetFunction("gaus");
      fitMaxX = skewFitMax[r][c]->GetMaximumX();
      cosAggFile->cd();
      //PMTMaxSpec[r][c]->SetTitle(Form("Max ADC Spect R%d C%d NPE%f",r,c,pow(f3->GetParameter(1)/pow(pow(f3->GetParameter(2),2.0)-pow(f1->GetParameter(2),2.0),0.5),2.0)));
      PMTMaxSpec[r][c]->SetTitle(Form("Max ADC Spect R%d C%d MaxX%f",r,c,fitMaxX));
      PMTMaxSpec[r][c]->Write(Form("Max ADC Spect R%d C%d",r,c));
      PMTMaxSpec[r][c]->Draw("AP");
      }

      //Calculate target HV from gain equation on data sheet
      if(iHV==1){
	targetHV[r][c] = jlab_lHVset/pow(fitMaxX/targetRAU,1.0/jlabExp);
      }else{
	targetHV[r][c] = cmu_lHVset/pow(fitMaxX/targetRAU,1.0/cmuExp);
      }
      
      //Spectra with TDC cut
      if( PMTIntSpecTDC[r][c]->GetEntries() > 0){
      cout << "PMTIntSpecTDC entries for row " << r << ", col " << c << " = " << PMTIntSpecTDC[r][c]->GetEntries() << "." << endl;
      PMTIntSpecTDC[r][c]->Fit(skewFitInt[r][c],"RQ");
      TDCfitIntX = skewFitIntTDC[r][c]->GetMaximumX();
      cosAggFile->cd();
      PMTIntSpecTDC[r][c]->SetTitle(Form("Int ADC Spect TDCcut R%d C%d MaxX%f",r,c,TDCfitIntX));
      PMTIntSpecTDC[r][c]->Write(Form("Int ADC SpecTDCt R%d C%d",r,c));
      PMTIntSpecTDC[r][c]->Draw("AP");
      }
      
      if( PMTMaxSpecTDC[r][c]->GetEntries() > 0){
      cout << "PMTMaxSpecTDC entries for row " << r << ", col " << c << " = " << PMTMaxSpecTDC[r][c]->GetEntries() << "." << endl;
      PMTMaxSpecTDC[r][c]->Fit(skewFitMax[r][c],"RQ");
      TDCfitMaxX = skewFitMaxTDC[r][c]->GetMaximumX();
      cosAggFile->cd();
      PMTMaxSpecTDC[r][c]->SetTitle(Form("Max ADC Spect TDCcut R%d C%d MaxX%f",r,c,TDCfitMaxX));
      PMTMaxSpecTDC[r][c]->Write(Form("Max ADC SpecTDCt R%d C%d",r,c));
      PMTMaxSpecTDC[r][c]->Draw("AP");
      }

      //Calculate target HV with TDC cut from gain equation on data sheet
      if(iHV==1){
	targetHVTDC[r][c] = jlab_lHVset/pow(TDCfitMaxX/targetRAU,1.0/jlabExp);
      }else{
	targetHVTDC[r][c] = cmu_lHVset/pow(TDCfitMaxX/targetRAU,1.0/cmuExp);
      }

      cout << "For row " << r << " and col " << c << ", target HV is " << targetHV[r][c] << " and target HV (TDC cut) is " << targetHVTDC[r][c] << "." << endl;
      outFile << r << "  " << c << "  " << targetHV[r][c] << "  " << targetHVTDC[r][c] << endl;
      
      
      //PMTMaxSpec[r][c]->Fit("gaus","Q");
      //f3=PMTMaxSpec[r][c]->GetFunction("gaus");
      
      
      //No need to fit, simply calculate calibration parameters
      /*
	double x[1] = {HV[iHV]};
	
	//double y[1] = {f3->GetParameter(1)}; //Mean of gaussian fit
	//double y[1] = {}; //"Mean" of skewed fit
	double y[1] = {PMTMaxSpec[r][c]->GetMaximumBin()};  //Simple sanity check on output tgraph
	
	double ex[1] = {0};
	//double ey[1] = {f3->GetParameter(2)/pow(PMTMaxSpec[r][c]->GetEntries(),0.5)};
	double ey[1] = {0};
	
	//TGraph *G3 = new TGraphErrors(HVpts,HV[iHV],f3->GetParameter(1),0,(f3->GetParameter(2)/pow(PMTMaxSpec[r][c]->GetEntries(),0.5)));
	TGraph *G3 = new TGraphErrors(HVpts,x,y,ex,ey);
	TF1 *pFit3 = new TF1("pFit3",powFit,0,2000,2);
	
	if(c>3&&c<8){
	pFit3->SetParameters(0,8);
	pFit3->SetParLimits(1,0,100);
	
	}else{
	pFit3->SetParameters(0,10);
	pFit3->SetParLimits(1,0,100);
	
	}
	
	G3->Fit("pFit3","Q");
	G3->SetTitle(Form("Max R%d-C%d: coeff=%g alpha=%f",r,c,pFit3->GetParameter(0),pFit3->GetParameter(1)));
	G3->GetYaxis()->SetTitle("Average fADC Amplitude (RAU)");  //Need actual units here
	G3->GetYaxis()->CenterTitle();
	G3->GetXaxis()->SetTitle("PMT HV Setting (-V)");
	G3->GetXaxis()->CenterTitle();
	//G3->SetMinimum(0.0);
	G3->SetLineColor(kWhite);
	G3->SetMarkerStyle(8);
	G3->SetMarkerSize(1);
	G3->Write(Form("Max R%d-C%d: coeff=%g alpha=%f",r,c,pFit3->GetParameter(0),pFit3->GetParameter(1)));
	G3->Draw("AP");
	
	}
      */
      
      /*
      //Spectra with TDC cut
      if( PMTIntSpecTDC[r][c]->GetEntries() > 0){
      cout << "PMTIntSpecTDC entries for row " << r << ", col " << c << " = " << PMTIntSpecTDC[r][c]->GetEntries() << "." << endl;
      PMTIntSpecTDC[r][c]->Fit("gaus","Q");
      f4=PMTIntSpecTDC[r][c]->GetFunction("gaus");
      cosAggFile->cd();
      PMTIntSpecTDC[r][c]->Write(Form("Int ADC Spect R%d C%d, TDC Cut",r,c));
      PMTIntSpecTDC[r][c]->Draw("AP");
      }
      
      if( PMTMaxSpecTDC[r][c]->GetEntries() > 0){
      cout << "PMTMaxSpecTDC entries for row " << r << ", col " << c << " = " << PMTMaxSpecTDC[r][c]->GetEntries() << "." << endl;
      PMTMaxSpecTDC[r][c]->Fit("gaus","Q");
      f5=PMTMaxSpecTDC[r][c]->GetFunction("gaus");
      cosAggFile->cd();
      PMTMaxSpecTDC[r][c]->SetTitle(Form("Max ADC Spect R%d C%d NPE%f, TDC Cut",r,c,pow(f5->GetParameter(1)/pow(pow(f5->GetParameter(2),2.0)-pow(f1->GetParameter(2),2.0),0.5),2.0)));
      PMTMaxSpecTDC[r][c]->Write(Form("Max ADC Spect R%d C%d, TDC Cut",r,c));
      PMTMaxSpecTDC[r][c]->Draw("AP");
      }
      //outFile << r << "  " << c << "  " << f3->GetParameter(1) << "  " << f2->GetParameter(1) << "  " << pow(f3->GetParameter(1)/pow(pow(f3->GetParameter(2),2.0)-pow(f1->GetParameter(2),2.0),0.5),2.0) << f5->GetParameter(1) << "  " << f4->GetParameter(1) << "  " << pow(f5->GetParameter(1)/pow(pow(f5->GetParameter(2),2.0)-pow(f1->GetParameter(2),2.0),0.5),2.0) << endl;  //Sigma error in quadrature, mean already pedestal subtracted
      
      */
            
    }
  }

  outFile << endl << endl;
  
  for(int r=0; r<kNrows; r++){
    for(int c=0; c<kNcols; c++){
      
      outFile << "Pedestal for row " << r << " and col " << c << " = " << pedestals[r][c] << endl;

    }
  }
  
  cout << "ADC spectrum histograms written to file cosSpectFile_run" << run << ".root" << endl;
  cout << "Sample histograms drawn to file HistosFile.root" << endl;
  cout << "Target HV settings written to cosCalOut_run" << run << ".txt" << endl;
  
  return 0;
}
  
