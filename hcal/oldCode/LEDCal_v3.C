#include <TH2F.h>
#include <TChain.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <TSystem.h>
#include <TStopwatch.h>
#include "hcal.h"

const int rmin = 0;
const int rmax = 12;
const int kNrows = 24; //All available rows for HCAL
const int kNcols = 12; //Half of HCAL has 12 columns
const int kNLED = 10; //Brute force to accomodate all bits

const int minSample = 0; //Initialize histogram minimum
const int maxSample = 50; //Initialize histogram maximum - cannot exceed hcalt::nsamps per pmt per event
const int totSample = (maxSample-minSample);

int gCurrentEntry = -1;

TChain *T = 0;
//TFile *HistosFile = new TFile(Form("%s/hcal_ledHistosFile.root",getenv("OUT_DIR")),"RECREATE");  //File for checking histogram fits and integrity

TFile *HistosFile = new TFile("hcal_ledHistosFile.root","RECREATE");

void pulseSpec(TH1F*, TF1*, double, int, int, int);

TH1F *ADCvChannel[kNLED];
TH1F *intADCvChannel[kNLED];
//TH1F *ADCTimevChannel[kNLED];
TH1F *NEVvChannel[kNLED];
TH1F *NPEvChannel[kNLED];
TH1F *pedvChannel[kNLED];

TH1F *histos[kNrows][kNcols][kNLED];
//TH1F *ADCTime[kNrows][kNcols][kNLED];
TH1F *PMTIntSpec[kNrows][kNcols][kNLED]; //=new TH1F("","",20,0,5000); //Empirical limits
TH1F *PMTMaxSpec[kNrows][kNcols][kNLED]; //=new TH1F("","",20,0,1000); //Empirical limits
TH1F *pedSpec[kNrows][kNcols][kNLED];

TF1 *f1;
TF1 *f2;
TF1 *f3;
TF1 *f4;
TF1 *f5;

int signalTotal = 1;

int DRAWNHISTO = 0; //Counts up per histo drawn
double maxy; //var for landau fit. Needed for maximum value of each histo
int ledNumber; //to convert from ledbit to led on
double pedVal; //Value passed into pedestals histogram
double pedSig; //std dev of pedestals histogram
double sampSum; //var for sum of samples from pedestal bins
double sampSumN; //Normalized sum of samples from ped bin

//Power-law fit
double powFit(double *x, double *par) {
  double amp = par[0];
  double ex = par[1];
  double loc = par[2];
  double offset = par[3];
  double ADC = x[0];
  return amp * pow(ADC,ex); //No shift, no normalization
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
  int led = hcalt::ledbit;  //Attempt 1 at including ledbit info in output histograms. One ledbit per event, no loops over ledbit.
    
  //Initializing peak, adc, and tdc arrays. No ledbit dimension as redefined by event.
  double peak[kNrows][kNcols]={0.0};
  double adc[kNrows][kNcols]={0.0};
  double amp[kNrows][kNcols]={0.0};
  double adc_p[kNrows][kNcols]={0.0};
  //double adc_time[kNrows][kNcols]={0.0};
  double amp_p[kNrows][kNcols]={0.0};
  double tdc[kNrows][kNcols]={0.0};
  double ped[kNrows][kNcols]={0.0};
  bool saturated[kNrows][kNcols]={false};

  //Get rows and columns from Tree for the event and fill histos and gsaturated
  for(int m = 0; m < hcalt::ndata; m++) {
    r = hcalt::row[m];
    c = hcalt::col[m];
    if(r < 0 || c < 0) {
      cout << "Error: row or col negative." << endl;
      continue;
    }
      
    if(r >= kNrows || c >= kNcols)
      continue;
    adc[r][c] = hcalt::a[m];
    tdc[r][c] = hcalt::tdc[m];
    adc_p[r][c] = hcalt::a_p[m];
    amp[r][c] = hcalt::a_amp[m];
    amp_p[r][c] = hcalt::a_amp_p[m];
    //adc_time[r][c] = hcalt::a_time[m];
    ped[r][c] = hcalt::ped[m];

    if( amp[r][c] > 4094 ) {
      saturated[r][c] = true;
      cout << "Saturation on r " << r << ", col " << c << "." << endl;
      //gSaturated++;
    }

  }
  
  //Perform basic checks on histograms and fill integrated/max pulse histograms
  for(r = rmin; r < kNrows; r++) {
    for(c = 0; c < kNcols; c++) {

      PMTIntSpec[r][c][led]->Fill( adc_p[r][c] );
      PMTMaxSpec[r][c][led]->Fill( amp_p[r][c] );
      //ADCTime[r][c][led]->Fill( adc_time[r][c] );
      pedSpec[r][c][led]->Fill( ped[r][c] );

    }
  }
}

int LEDCal_v3(int run = 1205, int event = 1000)  //Start after LEDs warm up ~1000 events
{
  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start(kTRUE);
  
  //Declare outfile
  TFile *ledAggFile = new TFile(Form("/adaqeb1/data1/hcal_ledSpectFile_%d.root",run),"RECREATE");

  //Build spectra histograms. Empirical limits.
  cout << "Building spectrum histograms.." << endl;
  for(int r=0; r<kNrows; r++){
    for(int c=0; c<kNcols; c++){
      for(int l=0; l<kNLED; l++){
	PMTIntSpec[r][c][l] = new TH1F(Form("Int ADC Spect R%d C%d LED%d",r,c,l),Form("Int ADC Spect R%d C%d LED%d",r,c,l),40,0,0);
	PMTIntSpec[r][c][l]->GetXaxis()->SetTitle("pC");
	PMTIntSpec[r][c][l]->GetXaxis()->CenterTitle();

	pedSpec[r][c][l] = new TH1F(Form("Pedestal Spect R%d C%d LED%d",r,c,l),Form("Pedestal Spect R%d C%d LED%d",r,c,l),40,0,0);
	pedSpec[r][c][l]->GetXaxis()->SetTitle("mV");
	pedSpec[r][c][l]->GetXaxis()->CenterTitle();

	//ADCTime[r][c][l] = new TH1F(Form("ADC Time R%d C%d LED%d",r,c,l),Form("Int ADC Time R%d C%d LED%d",r,c,l),40,0,0);
        //ADCTime[r][c][l]->GetXaxis()->SetTitle("ns");
	//ADCTime[r][c][l]->GetXaxis()->CenterTitle();

	PMTMaxSpec[r][c][l] = new TH1F(Form("Max ADC Spect R%d C%d LED%d",r,c,l),Form("Max ADC Spect R%d C%d LED%d",r,c,l),20,0,0);
	PMTMaxSpec[r][c][l]->GetXaxis()->SetTitle("mV");
	PMTMaxSpec[r][c][l]->GetXaxis()->CenterTitle();

      }     
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

    //ADCTimevChannel[l] = new TH1F(Form("ADCTimevChannel LED %d",l), Form("ADCTimevChannel LED %d",l), kNcols*kNrows, 0, kNcols*kNrows-1);
    //ADCTimevChannel[l]->GetXaxis()->SetTitle("Channel");
    //ADCTimevChannel[l]->GetXaxis()->CenterTitle();
    //ADCTimevChannel[l]->GetYaxis()->SetTitle("<ns>");
    //ADCTimevChannel[l]->GetYaxis()->CenterTitle();

    NEVvChannel[l] = new TH1F(Form("NEVvChannel LED %d",l), Form("NEVvChannel LED %d",l), kNcols*kNrows, 0, kNcols*kNrows-1);
    NEVvChannel[l]->GetXaxis()->SetTitle("Channel");
    NEVvChannel[l]->GetXaxis()->CenterTitle();
    NEVvChannel[l]->GetYaxis()->SetTitle("Events");
    NEVvChannel[l]->GetYaxis()->CenterTitle();

    pedvChannel[l] = new TH1F(Form("pedvChannel LED %d",l), Form("pedvChannel LED %d",l), kNcols*kNrows, 0, kNcols*kNrows-1);
    pedvChannel[l]->GetXaxis()->SetTitle("Channel");
    pedvChannel[l]->GetXaxis()->CenterTitle();
    pedvChannel[l]->GetYaxis()->SetTitle("<mV>");
    pedvChannel[l]->GetYaxis()->CenterTitle();
  }

    
  if(!T) { 
    T = new TChain("T");

    //T->Add(TString::Format("$ROOTFILES/hcal_%d*.root",run));
    //T->Add(Form("%s/e1209019_%d.root",getenv("DATA_DIR"),run));
    //T->Add(Form("~/sbs/Rootfiles/hcal_gmn_%d*.root",run));
    T->Add("/adaqfs/home/a-onl/sbs/Rootfiles/hcal_11358_1000.root");

    T->SetBranchStatus("*",0);
    T->SetBranchStatus("sbs.hcal.*",1);
    T->SetBranchAddress("sbs.hcal.ledbit",&hcalt::ledbit);
    T->SetBranchAddress("sbs.hcal.ledcount",&hcalt::ledcount);
    T->SetBranchAddress("sbs.hcal.nsamps",hcalt::nsamps);
    T->SetBranchAddress("sbs.hcal.a",hcalt::a);
    T->SetBranchAddress("sbs.hcal.a_p",hcalt::a_p);
    T->SetBranchAddress("sbs.hcal.a_amp",hcalt::a_amp);
    T->SetBranchAddress("sbs.hcal.a_amp_p",hcalt::a_amp_p);
    //T->SetBranchAddress("sbs.hcal.a_time",hcalt::a_time);
    T->SetBranchAddress("sbs.hcal.tdc",hcalt::tdc);
    T->SetBranchAddress("sbs.hcal.samps",hcalt::samps);
    T->SetBranchAddress("sbs.hcal.samps_idx",hcalt::samps_idx);
    T->SetBranchAddress("sbs.hcal.adcrow",hcalt::row);
    T->SetBranchAddress("sbs.hcal.adccol",hcalt::col);
    T->SetBranchAddress("sbs.hcal.ped",hcalt::ped);
    T->SetBranchStatus("Ndata.sbs.hcal.adcrow",1);
    T->SetBranchAddress("Ndata.sbs.hcal.adcrow",&hcalt::ndata);
    //cout << "Opened up tree with nentries=" << T->GetEntries() << endl;
  }
  
  if( T->GetEntries()<=0 ){
    cout << "Root file is empty. Check run number entry." << endl;
    return 0;
  }
  
  gCurrentEntry = event;

  cout << "Total events to process " << T->GetEntries() << endl;
  
  //Obtain average max/int value over histograms for all entries by led at each ledbit
  cout << "Building distributions.." << endl;
  for (int i = gCurrentEntry; i < T->GetEntries(); i++){ 
    processEvent(gCurrentEntry);
    gCurrentEntry++;
    
    //Keep count of events processed for max/int avg values
    if (gCurrentEntry%10000 == 0){
      cout << "Max/Int Value -> processing event: " << gCurrentEntry << endl;
    }
  }
  cout << "Finished loop over run " << run << "." << endl;
  //cout << "Total good signals = " << signalTotal << "." << endl;

  cout << "Writing spectrum histograms to file and creating file to hold fit parameters.." << endl;
  ofstream outFile;
  outFile.open(Form("NPE_run%d.txt",run));
  outFile << "#All information for run " << run << "." << endl;
  outFile << "#Row Col LED <Max> <Int> NPE" << endl;
  
  for(int r = rmin; r < kNrows; r++){ //Only left half of HCal on for these tests
    for(int c = 0; c < kNcols; c++){
      //pedSpec[r][c]->Fit("gaus","Q");
      //f1=pedSpec[r][c]->GetFunction("gaus");
      
      ledAggFile->cd();
      //pedSpec[r][c]->Write(Form("Pedestal Spect R%d C%d",r,c));
      //pedSpec[r][c]->Draw("AP");

      //pedvChannel->SetBinContent(kNcols*r+c+1,pedestals[r][c]);
      
      //Get mean of gaussian fit to pedestals histogram
      //outFile << "For row " << r << " and col " << c << ", pedestal is " << f1->GetParameter(1) << endl;
      
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
	/*
	if(ADCTime[r][c][l]->GetEntries()!=0){
	  ADCTime[r][c][l]->Fit("gaus","Q");
	  f4=ADCTime[r][c][l]->GetFunction("gaus");
	  ledAggFile->cd();
	  ADCTime[r][c][l]->Write(Form("ADC Time Spect R%d C%d LED%d",r,c,l));
	  ADCTime[r][c][l]->Draw("AP");
	
	  ADCTimevChannel[l]->SetBinContent(kNcols*r+c+1,f4->GetParameter(1));
	}
	*/
	NEVvChannel[l]->SetBinContent(kNcols*r+c+1,PMTIntSpec[r][c][l]->GetEntries());
	



	if(pedSpec[r][c][l]->GetEntries()!=0){
	  pedSpec[r][c][l]->Fit("gaus","Q");
	  f5=pedSpec[r][c][l]->GetFunction("gaus");
	  ledAggFile->cd();
	  pedSpec[r][c][l]->Write(Form("Pedestal Spect R%d C%d LED%d",r,c,l));
	  pedSpec[r][c][l]->Draw("AP");
	
	  pedvChannel[l]->SetBinContent(kNcols*r+c+1,f5->GetParameter(1));
	}
      }
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
    /*
    ledAggFile->cd();
    ADCTimevChannel[l]->SetTitle(Form("Average ADC Time vs Channel LED %d",l));
    ADCTimevChannel[l]->Write("ADCTimevChannel");
    ADCTimevChannel[l]->Draw("AP");
    */
    ledAggFile->cd();
    NEVvChannel[l]->SetTitle(Form("Number of Event Pulses vs Channel LED %d",l));
    NEVvChannel[l]->Write("NEVvChannel");
    NEVvChannel[l]->Draw("AP");

    ledAggFile->cd();
    pedvChannel[l]->SetTitle(Form("Pedestal vs Channel LED %d",l));
    pedvChannel[l]->Write("pedvChannel");
    pedvChannel[l]->Draw("AP");
  }

  cout << "ADC spectrum histograms written to file ledSpectFile_run" << run << ".root" << endl;
  cout << "NPE (and other parameters) written to file NPE_run" << run << ".txt." << endl;
  cout << "Total sample histograms drawn to file ledSpectFile_run" << run << ".root = " << DRAWNHISTO << "." << endl;

  st->Stop();

  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

  return 0;
}

