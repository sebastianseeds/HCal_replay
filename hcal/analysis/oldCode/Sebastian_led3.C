#include <TH2F.h>
#include <TChain.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <TSystem.h>
#include <TStopwatch.h>
#include "hcal.h"

const int rmin = 12;
const int kNrows = 24; //All available rows for HCAL
const int kNcols = 12; //Half of HCAL has 12 columns
const int kNLED = 6; //For 5 LEDs and an off position at zero

const int minSample = 0; //Initialize histogram minimum
const int maxSample = 30; //Initialize histogram maximum - cannot exceed hcalt::nsamps per pmt per event
const int totSample = (maxSample-minSample);

int gCurrentEntry = -1;

TChain *T = 0;
TFile *HistosFile = new TFile("Comparisons/HistosFile.root","RECREATE");  //File for checking histogram fits and integrity

void getPedestal(int);
void pulseSpec(TH1F*, TF1*, double, int, int, int);

TH1F *ADCvChannel[kNLED];
TH1F *intADCvChannel[kNLED];
TH1F *NEVvChannel[kNLED];
TH1F *NPEvChannel[kNLED];
TH1F *pedvChannel;

TH1F *histos[kNrows][kNcols][kNLED];
TH1F *PMTIntSpec[kNrows][kNcols][kNLED]; //=new TH1F("","",20,0,5000); //Empirical limits
TH1F *PMTMaxSpec[kNrows][kNcols][kNLED]; //=new TH1F("","",20,0,1000); //Empirical limits
TH1F *pedSpec[kNrows][kNcols];

double NPE[kNrows][kNcols][kNLED-1];
double LED[kNLED-1];
TGraph *NPEvLED[kNrows][kNcols];

TF1 *f1;
TF1 *f2;
TF1 *f3;
//vector<TH1F*> sampleHistos;

int signalTotal = 1;
int pedTotal = 1;

double pedestals[kNrows][kNcols];
double Vpedestals[kNrows][kNcols];
double pedSigma[kNrows][kNcols];
bool gSaturated[kNrows][kNcols][kNLED];

int DRAWNHISTO = 0; //Counts up per histo drawn
double maxy; //var for landau fit. Needed for maximum value of each histo
int ledNumber; //to convert from ledbit to led on
double pedVal; //Value passed into pedestals histogram
double pedSig; //std dev of pedestals histogram
double sampSum; //var for sum of samples from pedestal bins
double sampSumN; //Normalized sum of samples from ped bins

/*
//Create Gaussian for fitting.
double fit_gaus(double *x,double *p) 
{
  double fitval = p[0]*TMath::Exp(-0.5*pow(((x[0]-p[1])/p[2]),2));
  return fitval;
}
*/

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

int pInc = 0;

int INC = 0;

TH1F* MakeHisto(int row, int col, int led, int bins){
  TH1F *h = new TH1F(TString::Format("h%02d%02d LED%d",row,col,led),TString::Format("%d-%d LED-%d",row+1,col+1,led),bins,minSample,maxSample);
  h->SetStats(0);
  h->SetLineWidth(2);
  return h;
}

void getPedestal(int entry=-1){  //Taken only from the set of data where the LEDs are off.
  if(entry == -1) {
    gCurrentEntry++;
  } else {
    gCurrentEntry = entry;
  }

  if(gCurrentEntry<0) {
    gCurrentEntry = 0;
  }

  T->GetEntry(gCurrentEntry);

  //General vars by event
  int r,c,l,idx,n,sub;
  int ledbit = hcalt::ledbit; 
  ledNumber = (int) (ceil(log2(ledbit+1)));  //Verified map

  if(ledNumber==0){  //When ledNumber == 0, the LEDs are OFF

    // Clear old histograms, just in case modules are not in the tree
    for(r = 0; r < kNrows; r++) {
      for(c = 0; c < kNcols; c++) {
	for(l = 0; l < kNLED; l++){
	  histos[r][c][l]->Reset("ICES M");
	  gSaturated[r][c][l] = false;
	}
      }
    }

    
    //Get rows and columns from Tree for the event and fill histos and gsaturated
    for(int m = 0; m < hcalt::ndata; m++) {
      r = hcalt::row[m]-1;
      c = hcalt::col[m]-1;
      
      if(r < 0 || c < 0) {
        cerr << "Row/Col error: R" << r << " C" << c << "." << endl;
	continue;
      }
      
      if(r >= kNrows || c >= kNcols)
	continue;
      idx = hcalt::samps_idx[m];
      n = hcalt::nsamps[m];  
      bool processed = false;
      for(int s = minSample; s < maxSample && s < n; s++) {
	processed = true;
	histos[r][c][ledNumber]->SetBinContent(s+1-minSample,hcalt::samps[idx+s]);
      }
      
      if(!processed) {
	//cout << "Skipping empty module: " << m << endl;
	for(int s = 0;  s < totSample; s++) {
	  histos[r][c][ledNumber]->SetBinContent(s+1,-404); //CORRECT minsample
	}
      }
    }
    
    
    /*
      sampSum = 0.0;
      sampSumN = 0.0;
      
      for (int i=minSample; i<maxSample; i++){  
      sampSum+=histos[r][c][ledNumber]->GetBinContent(i+1); //Data starts at bin 1, not zero (underflow)
      }
      sampSumN = sampSum/maxSample; //Normalize by number of bins used per event
      
      pedSpec[r][c]->Fill(sampSumN); //Fill histogram with pedestal value from current event
      
      pInc++;
    */
    for(r = rmin; r < kNrows; r++) {
      for(c = 0; c < kNcols; c++) {

	for(int b=0 ; b<totSample ; b++) {
	  pedSpec[r][c]->Fill(histos[r][c][ledNumber]->GetBinContent(b+1));
	}
	
	if (ledNumber==0 && pedTotal % 35000 == 0 && r>11 && gSaturated[r][c][ledNumber]==false){ //Right half is shut off, hence r > 12, and cut all saturated events
	  //cout << "Saving a reference pedestal histogram to write to file.." << endl;
	  //cout << "Integrated value for this histogram = " << histos[r][c][ledNumber]->Integral(1,30) << "." << endl;
	  //cout << "Maximum value for this histogram = " << maxy << "." << endl;
	  //cout << "Pedestal for this histogram = " << pedestals[r][c] << "." << endl;
	  //cout << "LED number for this histogram = " << ledNumber << "." << endl << endl;
	  histos[r][c][ledNumber]->SetTitle(Form("Pedestal Histogram R%d C%d LED%d",r,c,ledNumber));
	  histos[r][c][ledNumber]->Draw();
	  HistosFile->cd();
	  histos[r][c][ledNumber]->SetTitle(Form("Pedestal Histogram R%d C%d LED%d",r,c,ledNumber));
	  histos[r][c][ledNumber]->GetYaxis()->SetTitle("RAU");
	  histos[r][c][ledNumber]->GetYaxis()->CenterTitle();
	  histos[r][c][ledNumber]->GetXaxis()->SetTitle("fADC Channel");
	  histos[r][c][ledNumber]->GetXaxis()->CenterTitle();
	  histos[r][c][ledNumber]->Write(Form("Pedestal Histogram R%d C%d LED%d",r,c,ledNumber));
	  histos[r][c][ledNumber]->Draw("AP");
	  
	}
	pedTotal++;
	/*  //Cannot use until TDC issue is corrected.
	    if(tdc[r][c]==0){ //eliminating adc pulses from coincident tdc measurement in pedestal calculation
	    
	    for(int b=0 ; b<totSample ; b++) {
	    pedSpec[r][c]->Fill(histos[r][c]->GetBinContent(b+1));
	    }
	    }
	*/
      }
    }   
  }
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

  if(gCurrentEntry%1000==0){
  cout << "Looping over entry " << gCurrentEntry << "." << endl;
  }
  
  //General vars by event
  int r,c,l,idx,n,sub;
  int ledbit = hcalt::ledbit;  //Attempt 1 at including ledbit info in output histograms. One ledbit per event, no loops over ledbit.
  
  ledNumber = (int) (ceil(log2(ledbit+1)));  //Verified map
  
  // Clear old histograms, just in case modules are not in the tree
  for(r = 0; r < kNrows; r++) {
    for(c = 0; c < kNcols; c++) {
      for(l = 0; l < kNLED; l++){
	histos[r][c][l]->Reset("ICES M");
	gSaturated[r][c][l] = false;
      }
    }
  }
    
  //Initializing peak, adc, and tdc arrays. No ledbit dimension as redefined by event.
  Float_t peak[kNrows][kNcols];
  Double_t adc[kNrows][kNcols];
  Double_t tdc[kNrows][kNcols];
  for(r  = 0; r < kNrows; r++) {
    for(c  = 0; c < kNcols; c++) {
      peak[r][c] = 0.0;
      adc[r][c] = 0.0;
      tdc[r][c] = 0.0;
    }
  }
    
  //Get rows and columns from Tree for the event and fill histos and gsaturated
  for(int m = 0; m < hcalt::ndata; m++) {
    r = hcalt::row[m]-1;
    c = hcalt::col[m]-1;
    if(r < 0 || c < 0) {
      cout << "Why is row negative? Or col?" << endl;
      continue;
    }
      
    if(r >= kNrows || c >= kNcols)
      continue;
    idx = hcalt::samps_idx[m];
    n = hcalt::nsamps[m];  
    adc[r][c] = hcalt::a[m];
    tdc[r][c] = hcalt::tdc[m];
    bool processed = false;
    for(int s = minSample; s < maxSample && s < n; s++) {
      processed = true;
      histos[r][c][ledNumber]->SetBinContent(s+1-minSample,hcalt::samps[idx+s]); 
      if(peak[r][c]<hcalt::samps[idx+s])
	peak[r][c]=hcalt::samps[idx+s];
      if(peak[r][c]>4095) {
	gSaturated[r][c][ledNumber] = true;
      }
    }
    if(!processed) {
      //cout << "Skipping empty module: " << m << endl;
      for(int s = 0;  s < totSample; s++) {
	histos[r][c][ledNumber]->SetBinContent(s+1,-404);
      }
    }
  }
    
  //Perform basic checks on histograms and fill integrated/max pulse histograms
  for(r = rmin; r < kNrows; r++) {
    for(c = 0; c < kNcols; c++) {
      histos[r][c][ledNumber]->SetTitle(Form("%d-%d LED-%d (ADC=%g,TDC=%g)",r+1,c+1,ledNumber,adc[r][c],tdc[r][c]));

      //Pedestal subtract
      //Get pedestal value and std dev for the event 
      pedVal = pedestals[r][c];
      pedSig = pedSigma[r][c];

      
      //Subtract pedestal value from each bin
      for(int b=1 ; b<=histos[r][c][ledNumber]->GetNbinsX() ; b++) {
	histos[r][c][ledNumber]->SetBinContent(b,histos[r][c][ledNumber]->GetBinContent(b)-pedVal);
      }
	
      //Obtain maximum and integrated value from the histogram
      maxy = histos[r][c][ledNumber]->GetMaximum();

      
      if(ledNumber!=0 && r>11 && gSaturated[r][c][ledNumber]==false){ //Right half is shut off, hence r > 12, and cut all saturated events
	PMTIntSpec[r][c][ledNumber]->Fill(histos[r][c][ledNumber]->Integral(1,30)); //Integral from total bin content
	PMTMaxSpec[r][c][ledNumber]->Fill(maxy);
      }
	
      //Write out regular histograms - LED bit changes every 1000 events.
      if (ledNumber!=0 && signalTotal % 35000 == 0 && r>11 && gSaturated[r][c][ledNumber]==false){ //Right half is shut off, hence r > 12, and cut all saturated events
	//cout << "Saving a reference good fit histogram to write to file.." << endl;
	//cout << "Integrated value for this histogram = " << histos[r][c][ledNumber]->Integral(1,30) << "." << endl;
	//cout << "Maximum value for this histogram = " << maxy << "." << endl;
	//cout << "Pedestal for this histogram = " << pedestals[r][c] << "." << endl;
	//cout << "TDC value for this histogram = " << tdc[r][c] << "." << endl;
	//cout << "LED number for this histogram = " << ledNumber << "." << endl << endl;
	histos[r][c][ledNumber]->SetTitle(Form("Sample Histo R%d C%d LED%d",r,c,ledNumber));
	histos[r][c][ledNumber]->Draw();
	HistosFile->cd();
	histos[r][c][ledNumber]->SetTitle(Form("Sample Event Histogram R%d C%d LED%d",r,c,ledNumber));
	histos[r][c][ledNumber]->GetYaxis()->SetTitle("RAU");
	histos[r][c][ledNumber]->GetYaxis()->CenterTitle();
	histos[r][c][ledNumber]->GetXaxis()->SetTitle("fADC Channel");
	histos[r][c][ledNumber]->GetXaxis()->CenterTitle();
	histos[r][c][ledNumber]->Write(Form("Sample Event Histogram R%d C%d LED%d",r,c,ledNumber));
	histos[r][c][ledNumber]->Draw("AP");
  
      }
      signalTotal++;
    }
  }
}

int Sebastian_led3(int run = 1205, int event = 1000)  //Start after LEDs warm up ~1000 events
{
  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start(kTRUE);
  
  //Declare outfile
  TFile *ledAggFile = new TFile(Form("Comparisons/ledSpectFile_run%d.root",run),"RECREATE");

  //Build spectra histograms. Empirical limits.
  cout << "Building spectrum histograms.." << endl;
  for(int r=0; r<kNrows; r++){
    for(int c=0; c<kNcols; c++){
      for(int l=0; l<kNLED; l++){
	PMTIntSpec[r][c][l] = new TH1F(Form("Int ADC Spect R%d C%d LED%d",r,c,l),Form("Int ADC Spect R%d C%d LED%d",r,c,l),20,0,0);
	PMTIntSpec[r][c][l]->GetXaxis()->SetTitle("sRAU");
	PMTIntSpec[r][c][l]->GetXaxis()->CenterTitle();

	PMTMaxSpec[r][c][l] = new TH1F(Form("Max ADC Spect R%d C%d LED%d",r,c,l),Form("Max ADC Spect R%d C%d LED%d",r,c,l),20,0,0);
	PMTMaxSpec[r][c][l]->GetXaxis()->SetTitle("RAU");
	PMTMaxSpec[r][c][l]->GetXaxis()->CenterTitle();

	if(l==0){
	  PMTIntSpec[r][c][l]->Fill(-404);  //l=0 corresponds to no led on, used only for extracting pedestals
	  PMTMaxSpec[r][c][l]->Fill(-404);  //l=0 corresponds to no led on, used only for extracting pedestals
	}
      }      
      pedSpec[r][c] = new TH1F(Form("Pedestal Spect R%d C%d",r,c),Form("Pedestal Spect R%d C%d",r,c),200,120,320);
      pedSpec[r][c]->GetXaxis()->SetTitle("<RAU>");
      pedSpec[r][c]->GetXaxis()->CenterTitle();
      
    }
  }

  //Build extra analysis histograms and array for NPEvLED TGraph
  for(int l=1; l<kNLED; l++){
    ADCvChannel[l] = new TH1F(Form("ADCvChannel LED %d",l), Form("ADCvChannel LED %d",l), kNcols*kNrows, 0, kNcols*kNrows-1);
    ADCvChannel[l]->GetXaxis()->SetTitle("Channel");
    ADCvChannel[l]->GetXaxis()->CenterTitle();
    ADCvChannel[l]->GetYaxis()->SetTitle("<RAU>");
    ADCvChannel[l]->GetYaxis()->CenterTitle();
    intADCvChannel[l] = new TH1F(Form("intADCvChannel LED %d",l), Form("intADCvChannel LED %d",l), kNcols*kNrows, 0, kNcols*kNrows-1);
    intADCvChannel[l]->GetXaxis()->SetTitle("Channel");
    intADCvChannel[l]->GetXaxis()->CenterTitle();
    intADCvChannel[l]->GetYaxis()->SetTitle("<RAU>");
    intADCvChannel[l]->GetYaxis()->CenterTitle();
    NEVvChannel[l] = new TH1F(Form("NEVvChannel LED %d",l), Form("NEVvChannel LED %d",l), kNcols*kNrows, 0, kNcols*kNrows-1);
    NEVvChannel[l]->GetXaxis()->SetTitle("Channel");
    NEVvChannel[l]->GetXaxis()->CenterTitle();
    NEVvChannel[l]->GetYaxis()->SetTitle("<RAU>");
    NEVvChannel[l]->GetYaxis()->CenterTitle();

    if(l!=0){
      LED[l-1]=l;
    }
  }
  pedvChannel = new TH1F("pedvChannel", "pedvChannel", kNcols*kNrows, 0, kNcols*kNrows-1);
  pedvChannel->GetXaxis()->SetTitle("Channel");
  pedvChannel->GetXaxis()->CenterTitle();
  pedvChannel->GetYaxis()->SetTitle("<RAU>");
  pedvChannel->GetYaxis()->CenterTitle();
    
  if(!T) { 
    T = new TChain("T");

    T->Add(TString::Format("rootFiles/LED/fadc_f1tdc_%d*.root",run));
    

    
    T->SetBranchStatus("*",0);
    T->SetBranchStatus("sbs.hcal.*",1);
    T->SetBranchAddress("sbs.hcal.ledbit",&hcalt::ledbit);
    T->SetBranchAddress("sbs.hcal.ledcount",&hcalt::ledcount);
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
	for(int l = 0; l < kNLED; l++) {
	  histos[r][c][l] = MakeHisto(r,c,l,totSample);
	  gSaturated[r][c][l] = false;
	}
      }
    }
  }
  
  gCurrentEntry = event;

  cout << "Total events to process " << T->GetEntries() << endl;
  
  //Obtain pedestal histograms from all entries for each led
  cout << "Looping over events to populate pedestal histograms from LED bit 00.." << endl;
  for (int i = gCurrentEntry; i < T->GetEntries(); i++){ 
    getPedestal(gCurrentEntry);
    gCurrentEntry++;
    
    //Keep count of events processed for pedestal
    if (gCurrentEntry%1000 == 0){
      cout << "Pedestal -> processing event: " << gCurrentEntry << endl;
    }
  }
  
  //Get mean and std dev from each pedestal spectra histo and put in array
  cout << "Obtaining mean and std dev from pedestal histograms.." << endl;  
  for(int r = rmin; r < kNrows; r++) {
    for(int c = 0; c < kNcols; c++) {
      if(pedSpec[r][c]->GetEntries()==0){
	cout << "Empty pedestal histograms at r c " << r << " " << c << "." << endl;
	continue;
      }else{

	if(r==6&&c==0) pedSpec[r][c]->Draw();
	
	pedSpec[r][c]->Fit("gaus","Q");
	f1=pedSpec[r][c]->GetFunction("gaus");
	pedestals[r][c]=f1->GetParameter(1); //Get the mean
	pedSigma[r][c]=f1->GetParameter(2); //Get std dev

      }
    }
  }
  
  cout << "Total pedestal events processed: " << pInc << "." << endl;
  
  gCurrentEntry = event; //Reset event count for Int/Max calculation
  
  //Obtain average max/int value over histograms for all entries by led at each ledbit
  cout << "Pedestal subtracting and obtaining integrated/max values from event histograms.." << endl;
  for (int i = gCurrentEntry; i < T->GetEntries(); i++){ 
    processEvent(gCurrentEntry);
    gCurrentEntry++;
    
    //Keep count of events processed for max/int avg values
    if (gCurrentEntry%10000 == 0){
      cout << "Max/Int Value -> processing event: " << gCurrentEntry << endl;
    }
  }
  cout << "Finished loop over run " << run << "." << endl;
  cout << "Total good signals = " << signalTotal << "." << endl;

  cout << "Reading in output parameter from Vanessa's code for comparison.." << endl;

  ifstream file1(Form("Comparisons/VPed_%d.txt",run));

  int n1=0;
  double d1;
  int rval, cval;
  string line;

  while(getline(file1,line)){
    if(line.at(0) == '#'){
      continue;
    }
    
    stringstream ss(line);
    ss >> d1;
    
    rval = floor(n1/kNcols);
    cval = n1 % kNcols;
    
    Vpedestals[rval][cval] = d1; 
    
    //cout << "row = " << rval << ", col = " << cval << ", HV val = " << pmtHV[rval][cval] << endl;
    
    n1++;
  }
  
  cout << "Writing spectrum histograms to file and creating file to hold fit parameters.." << endl;
  ofstream outFile;
  outFile.open(Form("Comparisons/NPE_run%d.txt",run));
  outFile << "#All information for run " << run << "." << endl;
  outFile << "#Row Col LED <Max> <Int> NPE" << endl;
  
  for(int r = rmin; r < kNrows; r++){ //Only left half of HCal on for these tests
    for(int c = 0; c < kNcols; c++){
      pedSpec[r][c]->Fit("gaus","Q");
      f1=pedSpec[r][c]->GetFunction("gaus");
      
      ledAggFile->cd();
      pedSpec[r][c]->Write(Form("Pedestal Spect R%d C%d",r,c));
      pedSpec[r][c]->Draw("AP");

      pedvChannel->SetBinContent(kNcols*r+c+1,pedestals[r][c]);
      
      //Get mean of gaussian fit to pedestals histogram
      outFile << "For row " << r << " and col " << c << ", pedestal is " << f1->GetParameter(1) << endl;
      
      for(int l=1; l<kNLED; l++){
	if(PMTIntSpec[r][c][l]->GetEntries()!=0){
	  PMTIntSpec[r][c][l]->Fit("gaus","Q");
	  f2=PMTIntSpec[r][c][l]->GetFunction("gaus");
	  ledAggFile->cd();
	  PMTIntSpec[r][c][l]->Write(Form("Int ADC Spect R%d C%d LED%d",r,c,l));
	  PMTIntSpec[r][c][l]->Draw("AP");
	
	  intADCvChannel[l]->SetBinContent(kNcols*r+c+1,f2->GetParameter(1));
	}

	if(PMTMaxSpec[r][c][l]->GetEntries()!=0){
	  PMTMaxSpec[r][c][l]->Fit("gaus","Q");
	  f3=PMTMaxSpec[r][c][l]->GetFunction("gaus");
	  ledAggFile->cd();
	  PMTMaxSpec[r][c][l]->SetTitle(Form("Max ADC Spect R%d C%d LED%d NPE%f",r,c,l,pow(f3->GetParameter(1)/f3->GetParameter(2),2)));
	  PMTMaxSpec[r][c][l]->Write(Form("Max ADC Spect R%d C%d LED%d",r,c,l));
	  PMTMaxSpec[r][c][l]->Draw("AP");

	  ADCvChannel[l]->SetBinContent(kNcols*r+c+1,f3->GetParameter(1));
	}
	
	NEVvChannel[l]->SetBinContent(kNcols*r+c+1,PMTIntSpec[r][c][l]->GetEntries());
	
	if(PMTMaxSpec[r][c][l]->GetEntries()!=0){

	  outFile << r << "  " << c << "  " << l << "  " << f3->GetParameter(1) << "  " << f2->GetParameter(1) << "  " << pow(f3->GetParameter(1)/pow(pow(f3->GetParameter(2),2.0)-pow(f1->GetParameter(2),2.0),0.5),2.0) << endl;  //Sigma error in quadrature, mean already pedestal subtracted

	  //cout << r << "  " << c << "  " << l << "  " << f3->GetParameter(1) << "  " << f2->GetParameter(1) << "  " << pow(f3->GetParameter(1)/pow(pow(f3->GetParameter(2),2.0)-pow(f1->GetParameter(2),2.0),0.5),2.0) << endl;  //Sigma error in quadrature, mean already pedestal subtracted

	  //if(l==1) cout << "Module: " << r*12+c << ",  Pedestal: " << pedestals[r][c] << endl;

	  if(l==1) cout << "Module: " << r*12+c << ", Pedestal Difference (V-S): " << Vpedestals[r][c]-pedestals[r][c] << endl;
	  
	  NPE[r][c][l-1]=pow(f3->GetParameter(1)/pow(pow(f3->GetParameter(2),2.0)-pow(f1->GetParameter(2),2.0),0.5),2.0);
	  
	}
      }
    }
  }  

  //Draw NPE by pmt histograms
  for(int r=rmin; r<kNrows; r++){
    for(int c=0; c<kNcols; c++){      
      NPEvLED[r][c]=new TGraph(kNLED-1,LED,NPE[r][c]);
      ledAggFile->cd();
      NPEvLED[r][c]->SetTitle(Form("Avg NPE vs LED R%d C%d",r,c));
      NPEvLED[r][c]->GetXaxis()->SetTitle("LED Setting");
      NPEvLED[r][c]->GetXaxis()->CenterTitle();
      NPEvLED[r][c]->GetYaxis()->SetTitle("<NPE>");
      NPEvLED[r][c]->GetYaxis()->CenterTitle();
      NPEvLED[r][c]->SetMinimum(0.0);
      NPEvLED[r][c]->SetLineColor(kWhite);
      NPEvLED[r][c]->SetMarkerStyle(8);
      NPEvLED[r][c]->Draw("P same");
      NPEvLED[r][c]->Write(Form("NPEvLEDr%dc%d",r,c));
    }
  }
      
  //Draw analysis histograms
  for(int l=1; l<kNLED; l++){
    ledAggFile->cd();
    ADCvChannel[l]->SetTitle(Form("Average Max ADC Val vs Channel LED %d",l));
    ADCvChannel[l]->Write("ADCvChannel");
    ADCvChannel[l]->Draw("AP");
    
    ledAggFile->cd();
    intADCvChannel[l]->SetTitle(Form("Average Int ADC Val vs Channel LED %d",l));
    intADCvChannel[l]->Write("intADCvChannel");
    intADCvChannel[l]->Draw("AP");
    
    ledAggFile->cd();
    NEVvChannel[l]->SetTitle(Form("Number of Event Pulses vs Channel LED %d",l));
    NEVvChannel[l]->Write("NEVvChannel");
    NEVvChannel[l]->Draw("AP");
  }
  
  ledAggFile->cd();
  pedvChannel->SetTitle("Avg Pedestal vs Channel");
  pedvChannel->Write("pedvChannel");
  pedvChannel->Draw("AP");

  cout << "Alpha parameters written to alphas.txt" << endl;
  cout << "ADC spectrum histograms written to file ledSpectFile_run" << run << ".root" << endl;
  cout << "NPE (and other parameters) written to file NPE_run" << run << ".txt." << endl;
  cout << "Total sample histograms drawn to file ledSpectFile_run" << run << ".root = " << DRAWNHISTO << "." << endl;

  st->Stop();

  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

  return 0;
}

