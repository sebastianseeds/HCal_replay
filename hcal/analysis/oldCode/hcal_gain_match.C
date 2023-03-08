//SSeeds August 2021 - Cosmic PMT HV calibration code for use with SBS HCAL experiments. 24 rows x 12 cols HCAL modules. Relevant information can be found here: hallaweb.jlab.org/wiki/index.php/HCAL_Cosmics.

#include <TH2F.h>
#include <TChain.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <TSystem.h>
#include <TStopwatch.h>
#include <iomanip>
#include <ctime>
#include "hcal.h"

// Detector parameters and flags
const int kNrows = 24;
const int kNcols = 12;
const int kCanvSize = 100;

const double kTargetpC = 406/50; // Set by looking for signature amp saturation distribution between maxADC values of 975 mv to 1000 mV and taking the rising edge start. 1000 is the minimum maxADC value to saturate where the amp dynamic range is 0-3V and the average attenuation over cables is 1/3. Looking at 14MeV/700MeV for mean energy deposition from muons by peak energy deposition in highest GMn kinematic (E_m/E_gmn = 1/50), so multiply this target (375 mV) by 1/50 to maximize dynamic range of the PMT and avoid saturation on amplifier where (owing to attenuation over long cables) saturation occurs there before the fADC

// Counter to keep track of T tree entry for processing
int gCurrentEntry = -1;
TChain *T = 0;

// Augment fn's for cosmic Calibration and recording
string getDate( );

// Declare necessary histograms and canvas
TH1F *gIntSpec[kNrows][kNcols]; 
TH1F *gMaxSpec[kNrows][kNcols]; 
TH1F *gIntSpecTDC[kNrows][kNcols]; 
TH1F *gMaxSpecTDC[kNrows][kNcols]; 
TH1F *gTDCSpec[kNrows][kNcols];
TH1F *gPedSpec[kNrows][kNcols];
TH1F *gTHVvChannel;
TH1F *gADCvChannel;
TH1F *gAmpvChannel;
TH1F *gNEVvChannel;
TH1F *gPedvChannel;
TH1F *gTDCvChannel;
TH1F *gTDCRMSvChannel;
TH2F *gNEV;

// Declare fit function
TF1 *g1;

// Create marker to track number of signals that pass cuts and saturation
int gSignalTotal = 0;
int gSaturated = 0;

// Declare necessary arrays
double gHV[kNrows][kNcols];
double gAlphas[kNrows][kNcols];

// Make machine-based date function
string getDate(){
  time_t now = time(0);
  tm ltm = *localtime(&now);

  string yyyy = std::to_string(1900 + ltm.tm_year);
  string mm = std::to_string(1 + ltm.tm_mon);
  string dd = std::to_string(ltm.tm_mday);
  string date = mm + '/' + dd + '/' + yyyy;

  return date;
} 

void processEvent(int entry = -1, double cut = 5, int diagPlots = 0 ){
  // Check event increment and increment
  if( entry == -1 ) {
    gCurrentEntry++;
  } else {
    gCurrentEntry = entry;
  }

  if( gCurrentEntry < 0 ) {
    gCurrentEntry = 0;
  }

  // Get the event from the TTree
  T->GetEntry( gCurrentEntry );
  
  int r,c;

  // Declare signal amplitude, adc, and tdc arrays for this event
  double adc[kNrows][kNcols] = {0.0};
  double tdc[kNrows][kNcols] = {0.0};
  double adc_p[kNrows][kNcols] = {0.0};
  double amp[kNrows][kNcols] = {0.0};
  double amp_p[kNrows][kNcols] = {0.0};
  double ped[kNrows][kNcols] = {0.0};
  bool saturated[kNrows][kNcols] = { {false} };

  // Process event with m data
  for( int m = 0; m < hcalt::ndata; m++ ) {
    // Define row and column
    r = hcalt::row[m];
    c = hcalt::col[m];
    if( r < 0 || c < 0 ) {
      cerr << "Error: row or col negative." << endl;
      continue;
    }
    
    if( r >= kNrows || c >= kNcols ) continue;

    // Fill adc and tdc arrays
    adc[r][c] = hcalt::a[m];
    adc_p[r][c] = hcalt::a_p[m];
    amp[r][c] = hcalt::a_amp[m];
    amp_p[r][c] = hcalt::a_amp_p[m];
    tdc[r][c] = hcalt::tdc[m];
    if( diagPlots==1 ) ped[r][c] = hcalt::ped[m];
    // Mark saturated array when amplitude meets max RAU
    if( amp[r][c] > 3900 ) {
      saturated[r][c] = true;
      //cout << "Saturation on r " << r << ", col " << c << endl;
      gSaturated++;
    }
  }

  // Vertical line test - accept only if three hits in vertical track. Reject if hit adjacent in row. Assuming pedestal sigma reasonable for all channels all events (sigma << 5 mV or 10 RAU (observed roughly 1 mV or 2 RAU))
 
  for( r = 0; r < kNrows; r++ ) {
    for( c = 0; c < kNcols; c++ ) {

      if( diagPlots==1 ) gPedSpec[r][c]->Fill( ped[r][c] );

      if( amp_p[r][c]>cut && saturated[r][c]==false ) { // Only examine channels with max > 5 sigma pedestal for hits and no saturation
	if( r==0 && c==0 ) { // Logic for top left channel
	  if( amp_p[r+1][c] > cut && amp_p[r+2][c] > cut ){
	    if( amp_p[r][c+1] < 6 ) {
	      if( tdc[r][c] != 0 ) {
		gIntSpecTDC[r][c]->Fill( adc_p[r][c] );
		gMaxSpecTDC[r][c]->Fill( amp_p[r][c] );
	      }
	      gIntSpec[r][c]->Fill( adc_p[r][c] );
	      gMaxSpec[r][c]->Fill( amp_p[r][c] );
	      gTDCSpec[r][c]->Fill( tdc[r][c] );
	      gSignalTotal++;
	    }
	  }
	}
	else if( r==0 && c==11 ) { // Logic for top right channel
	  if( amp_p[r+1][c] > cut && amp_p[r+2][c] > cut ){
	    if( amp_p[r][c-1] < 6 ) {
	      if( tdc[r][c] != 0 ) {
		gIntSpecTDC[r][c]->Fill( adc_p[r][c] );
		gMaxSpecTDC[r][c]->Fill( amp_p[r][c] );
	      }
	      gIntSpec[r][c]->Fill( adc_p[r][c] );
	      gMaxSpec[r][c]->Fill( amp_p[r][c] );
	      gTDCSpec[r][c]->Fill( tdc[r][c] );
	      gSignalTotal++;
	    }
	  }
	}
	else if( r==23 && c==0 ) { // Logic for bottom left channel
	  if( amp_p[r-1][c] > cut && amp_p[r-2][c] > cut ){
	    if( amp_p[r][c+1] < 6 ) {
	      if( tdc[r][c] != 0 ) {
		gIntSpecTDC[r][c]->Fill( adc_p[r][c] );
		gMaxSpecTDC[r][c]->Fill( amp_p[r][c] );
	      }
	      gIntSpec[r][c]->Fill( adc_p[r][c] );
	      gMaxSpec[r][c]->Fill( amp_p[r][c] );
	      gTDCSpec[r][c]->Fill( tdc[r][c] );
	      gSignalTotal++;
	    }
	  }
	}
	else if( r==23 && c==11 ) { // Logic for bottom right channel
	  if( amp_p[r-1][c] > cut && amp_p[r-2][c] > cut ){
	    if( amp_p[r][c-1] < 6 ) {
	      if( tdc[r][c] != 0 ) {
		gIntSpecTDC[r][c]->Fill( adc_p[r][c] );
		gMaxSpecTDC[r][c]->Fill( amp_p[r][c] );
	      }
	      gIntSpec[r][c]->Fill( adc_p[r][c] );
	      gMaxSpec[r][c]->Fill( amp_p[r][c] );
	      gTDCSpec[r][c]->Fill( tdc[r][c] );
	      gSignalTotal++;
	    }
	  }
	}
	else if( r==0 || r==1) { // Logic for other row 1 or 2 channels
	  if( (amp_p[r+1][c] > cut && amp_p[r+2][c] > cut) ||
	      (r==1 && amp_p[r-1][c] > cut && amp_p[r+1][c] > cut) ) {
	    if( amp_p[r][c-1] < 6 && amp_p[r][c+1] < 6 ) {
	      if( tdc[r][c] != 0 ) {
		gIntSpecTDC[r][c]->Fill( adc_p[r][c] );
		gMaxSpecTDC[r][c]->Fill( amp_p[r][c] );
	      }
	      gIntSpec[r][c]->Fill( adc_p[r][c] );
	      gMaxSpec[r][c]->Fill( amp_p[r][c] );
	      gTDCSpec[r][c]->Fill( tdc[r][c] );
	      gSignalTotal++;
	    }
	  }
	}
	else if( r==23 || r==22 ) { // Logic for other row 22 or 23 channels
	  if( (amp_p[r-1][c] > cut && amp_p[r-2][c] > cut) ||
	      (r==22 && amp_p[r-1][c] > cut && amp_p[r+1][c] > cut) ) {
	    if( amp_p[r][c-1] < 6 && amp_p[r][c+1] < 6 ) {
	      if( tdc[r][c] != 0 ) {
		gIntSpecTDC[r][c]->Fill( adc_p[r][c] );
		gMaxSpecTDC[r][c]->Fill( amp_p[r][c] );
	      }
	      gIntSpec[r][c]->Fill( adc_p[r][c] );
	      gMaxSpec[r][c]->Fill( amp_p[r][c] );
	      gTDCSpec[r][c]->Fill( tdc[r][c] );
	      gSignalTotal++;
	    }
	  }
	}
	else { // Logic for all other channels
	  if( (amp_p[r-2][c] > cut && amp_p[r-1][c] > cut) ||
	      (amp_p[r-1][c] > cut && amp_p[r+1][c] > cut) ||
	      (amp_p[r+1][c] > cut && amp_p[r+2][c] > cut) ) {
	    if( amp_p[r][c-1] < 6 && amp_p[r][c+1] < 6 ) {
	      if( tdc[r][c] != 0 ) {
		gIntSpecTDC[r][c]->Fill( adc_p[r][c] );
		gMaxSpecTDC[r][c]->Fill( amp_p[r][c] );
	      }
	      gIntSpec[r][c]->Fill( adc_p[r][c] );
	      gMaxSpec[r][c]->Fill( amp_p[r][c] );
	      gTDCSpec[r][c]->Fill( tdc[r][c] );
	      gSignalTotal++;
	    }
	  }
	}
      }
    }
  }
}

// Main
int hcal_gain_match(int run = -1, int event = -1){
    
  // Initialize function with user commands
  bool diagPlots = 0;
  bool vertCut = 0;
  bool HVLimit = 0;
  string date = getDate();
  double adcCut = 5;

  cout << "Enter run number for analysis." << endl;
  cin >> run;
  cout << "Cut on vertical tracks? Enter 0 for NO and 1 for YES." << endl;
  cin >> vertCut;
  cout << "Print diagnostic plots? Enter 0 for NO and 1 for YES." << endl;
  cin >> diagPlots;
  cout << "Enter ADC pulse threshold over pedestal (in mV)." << endl;
  cin >> adcCut;
  cout << "Limit HV Targets to within 100V of HV settings? Enter 0 for NO and 1 for YES." << endl;
  cin >> HVLimit;
  
  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );
  
  //Declare outfile
  TFile *cosHistFile = new TFile( Form( "outFiles/gainMatchResults_run%d.root", run ), "RECREATE" );
  ofstream fitData;
  fitData.open( Form( "outFiles/HCal_%d_FitData.txt", run ) );
  fitData << "Run Number: " << run << " Desired Peak Position: " << kTargetpC << endl;
  fitData << "Module " << " " << " Target HV " << " " << " Stat " << " " << " ErrStat " << " " << " Peak Pos " << " " << " ErrPPos " << " " << " Peak Width " << " " << " ErrPWid " << " " << " NGoodEvents " << " " <<  " " << " Flag " << std::endl;

  //Set ADC spectrum histogram limits
  int SpecBins = 150;
  //INT
  double IntSpec_min = 0.0, IntSpec_max = 75.0;
  //MAX
  double MaxSpec_min = 0.0, MaxSpec_max = 300.0;

  //Set TDC spectrum histogram limits
  int TDCBins = 400;
  double TDC_min = -200.0, TDC_max = 200.0;

  //Build spectrum histograms. Empirical limits.
  cout << "Building ADC and TDC spectrum histograms.." << endl;

  for( int r=0; r<kNrows; r++ ){
    for( int c=0; c<kNcols; c++ ){  

      int m = r*12+c;

      gIntSpec[r][c] = new TH1F( Form( "Int ADC Spect R%d C%d PMT%d", r, c, m ), Form( "Int ADC Spect R%d C%d PMT%d", r, c, m ), SpecBins, IntSpec_min, IntSpec_max );
      gIntSpec[r][c]->GetXaxis()->SetTitle( "pC" );
      gIntSpec[r][c]->GetXaxis()->CenterTitle();
      
      gMaxSpec[r][c] = new TH1F( Form( "Max ADC Spect R%d C%d PMT%d", r, c, m ), Form( "Max ADC Spect R%d C%d PMT%d", r, c, m ), SpecBins, MaxSpec_min, MaxSpec_max );
      gMaxSpec[r][c]->GetXaxis()->SetTitle( "mV" );
      gMaxSpec[r][c]->GetXaxis()->CenterTitle();
      
      gIntSpecTDC[r][c] = new TH1F( Form( "Int ADC Spect R%d C%d PMT%d, TDC Cut", r, c, m ), Form( "Int ADC Spect R%d C%d PMT%d, TDC Cut", r, c, m ), SpecBins, IntSpec_min, IntSpec_max );
      gIntSpecTDC[r][c]->GetXaxis()->SetTitle( "pC" );
      gIntSpecTDC[r][c]->GetXaxis()->CenterTitle();
      
      gMaxSpecTDC[r][c] = new TH1F( Form( "Max ADC Spect R%d C%d PMT%d, TDC Cut", r, c, m ), Form( "Max ADC Spect R%d C%d PMT%d, TDC Cut", r, c, m ), SpecBins, MaxSpec_min, MaxSpec_max );
      gMaxSpecTDC[r][c]->GetXaxis()->SetTitle( "mV" );
      gMaxSpecTDC[r][c]->GetXaxis()->CenterTitle();

      gPedSpec[r][c] = new TH1F( Form( "Pedestal Spect R%d C%d PMT%d", r, c, m ), Form( "Pedestal Spect R%d C%d PMT%d", r, c, m ), SpecBins, MaxSpec_min, MaxSpec_max );
      gPedSpec[r][c]->GetXaxis()->SetTitle( "mV" );
      gPedSpec[r][c]->GetXaxis()->CenterTitle();

      gTDCSpec[r][c] = new TH1F( Form( "TDC Spect%d C%d PMT%d", r, c, m ), Form( "Pedestal Spect R%d C%d PMT%d", r, c, m ), TDCBins, TDC_min, TDC_max );
      gTDCSpec[r][c]->GetXaxis()->SetTitle( "t-ref difference (ns)" );
      gTDCSpec[r][c]->GetXaxis()->CenterTitle();
    }
  }

  // Build diagnostic histograms
  if( diagPlots==1 ){
    gADCvChannel = new TH1F( "ADCvChannel", "ADC vs Channel", kNcols*kNrows, 0, kNcols*kNrows );
    gADCvChannel->SetMaximum(100);
    gAmpvChannel = new TH1F( "AmpvChannel", "Amplitude vs Channel", kNcols*kNrows, 0, kNcols*kNrows );
    gNEVvChannel = new TH1F( "NEVvChannel", "Number of events vs Channel", kNcols*kNrows, 0, kNcols*kNrows );
    gNEV = new TH2F( "NEV", "Number of events heatmap", kNrows, 0, kNrows, kNcols, 0, kNcols );
    gPedvChannel = new TH1F( "PedvChannel", "Pedestal vs Channel", kNcols*kNrows, 0, kNcols*kNrows );
    gTHVvChannel = new TH1F( "THVvChannel", "Target HV vs Channel", kNcols*kNrows, 0, kNcols*kNrows );
    gTDCvChannel = new TH1F( "TDCvChannel", "TDC vs Channel", kNcols*kNrows, 0, kNcols*kNrows );
    gTDCRMSvChannel = new TH1F( "TDCRMSvChannel", "TDC RMS vs Channel", kNcols*kNrows, 0, kNcols*kNrows );
  }
  
  // Read in data produced by analyzer in root format
  cout << "Reading replayed root file.." << endl;
  if( !T ) { 
    T = new TChain( "T" );
    T->Add( Form( "/volatile/halla/sbs/seeds/rootfiles/hcal_%d*.root", run ) );
    T->SetBranchStatus( "*", 0 );
    T->SetBranchStatus( "sbs.hcal.*", 1 );
    T->SetBranchAddress( "sbs.hcal.a", hcalt::a );
    T->SetBranchAddress( "sbs.hcal.a_amp", hcalt::a_amp );
    T->SetBranchAddress( "sbs.hcal.a_p", hcalt::a_p );
    T->SetBranchAddress( "sbs.hcal.a_amp_p", hcalt::a_amp_p );
    T->SetBranchAddress( "sbs.hcal.tdc", hcalt::tdc );
    T->SetBranchAddress( "sbs.hcal.ped", hcalt::ped );
    T->SetBranchAddress( "sbs.hcal.adcrow", hcalt::row );
    T->SetBranchAddress( "sbs.hcal.adccol", hcalt::col );
    T->SetBranchStatus( "Ndata.sbs.hcal.adcrow", 1 );
    T->SetBranchAddress( "Ndata.sbs.hcal.adcrow", &hcalt::ndata );
    cout << "Opened up tree with nentries=" << T->GetEntries() << endl;
  }

  // Set appropriate HV and alphas for run. HV settings from HCAL wiki. Alphas from LED analysis. Must have accompanying text file, one double per line by module number. Assuming negative voltage inputs.

  ifstream HVFile( Form( "setFiles/HV_run%d.txt", run ) );

  if( !HVFile ){
    cerr << "No HV settings from run present -> setFiles/HV_run" << run << ".txt expected." << endl;
    return 0;
  }
  
  cout << "Getting HV settings for each pmt for run " << run << "." << endl;

  int n1=0;
  double d1;
  
  int rval, cval;
  string line;
    
  while( getline( HVFile, line ) ){
    if( line.at(0) == '#' ) {
      continue;
    }
    
    stringstream ss( line );
    ss >> d1;
    
    rval = floor( n1/kNcols );
    cval = n1 % kNcols;
    
    gHV[rval][cval] = -d1; 
        
    n1++;
  }

  ifstream alphaFile( "setFiles/alphas.txt" );

  if( !alphaFile ){
    cerr << "No PMT alphas file present -> setFiles/alphas.txt expected." << endl;
    return 0;
  }
  
  n1=0;
  string line2;

  while( getline( alphaFile, line2 ) ){
    if( line2.at( 0 )=='#' ) {
      continue;
    }

    stringstream ss( line2 );
    ss >> d1;

    rval = floor( n1/kNcols );
    cval = n1 % kNcols;

    gAlphas[rval][cval] = d1;
    
    n1++;
  }

  gCurrentEntry = event;

  cout << "Total events to process " << T->GetEntries() << ".." << endl;

  cout << "Processing hits by event.." << endl;

  for ( int i = gCurrentEntry; i < T->GetEntries(); i++ ){ 
    processEvent( gCurrentEntry, adcCut, diagPlots );
    gCurrentEntry++;
    
    // Keep count of events processed for monitoring
    if ( gCurrentEntry%10000 == 0 ){
      cout << "Current event: " << gCurrentEntry << endl;
    }
  }

  // Fit ADC spectra with landau to extract mpv for HV calibration
  cout << "Writing spectrum histograms and calibration constants to file.." << endl;
  ofstream reportFile;
  reportFile.open( "outFiles/HVReport.txt" ); // Create text file to hold target HVs
  ofstream outFile;
  outFile.open( Form("outFiles/HVTargets_run%d.txt", run ) ); // Create text file to hold target HVs
  ofstream HVsetFile;
  HVsetFile.open( "setFiles/HV_run####.txt" );
  time_t now = time(0); 
  char *dt = ctime(&now);
  reportFile << "#Target HV settings for run " << run << " from hcal_gain_match.C on " << dt << "#" << endl;
  reportFile << "#Row Col targetHVUnlimited targetHV setHV Error" << endl;

  outFile << "Module" << '\t' << " HV" << endl;

  // Implement two-fit scheme
  TF1 *landauInt[kNrows][kNcols] = {};
  TF1 *landauMax[kNrows][kNcols] = {};

  // Fit parameter array for built-in landau
  double parsInt[3];
  double parErrInt[3];
  double parsMax[3];
  double parErrMax[3];

  // Array to store target high voltage per module as it's calculated
  double targetHV[kNrows][kNcols];
  double targetHVUnlimited[kNrows][kNcols];

  // Fit all spectra histograms
  for(int r=0; r<kNrows; r++){
    for(int c=0; c<kNcols; c++){

      // Create fit functions for each module
      // INT
      landauInt[r][c] = new TF1(Form("landauInt_r%d_c%d",r,c),"landau");
      landauInt[r][c]->SetLineColor(2);
      landauInt[r][c]->SetNpx(1000);
      // MAX
      landauMax[r][c] = new TF1(Form("landauMax_r%d_c%d",r,c),"landau");;
      landauMax[r][c]->SetLineColor(4);
      landauMax[r][c]->SetNpx(1000);

      // First fit
      // INT
      int maxBinInt = gIntSpec[r][c]->GetMaximumBin();
      double maxBinCenterInt = gIntSpec[r][c]->GetXaxis()->GetBinCenter( maxBinInt );
      double maxCountInt = gIntSpec[r][c]->GetMaximum();
      double binWidthInt = gIntSpec[r][c]->GetBinWidth( maxBinInt );
      double stdDevInt = gIntSpec[r][c]->GetStdDev();
      // MAX
      int maxBinMax = gMaxSpec[r][c]->GetMaximumBin();
      double maxBinCenterMax = gMaxSpec[r][c]->GetXaxis()->GetBinCenter( maxBinMax );
      double maxCountMax = gMaxSpec[r][c]->GetMaximum();
      double binWidthMax = gMaxSpec[r][c]->GetBinWidth( maxBinMax );
      double stdDevMax = gMaxSpec[r][c]->GetStdDev();
      // INT
      landauInt[r][c]->SetParameters( maxCountInt, maxBinCenterInt, stdDevInt );
      landauInt[r][c]->SetRange( IntSpec_min, IntSpec_max );
      gIntSpec[r][c]->Fit(landauInt[r][c],"No+RQ");
      landauInt[r][c]->GetParameters( parsInt );
      // MAX
      landauMax[r][c]->SetParameters( maxCountMax, maxBinCenterMax, stdDevMax );
      landauMax[r][c]->SetRange( MaxSpec_min, MaxSpec_max );
      gMaxSpec[r][c]->Fit(landauMax[r][c],"No+RQ");
      landauMax[r][c]->GetParameters( parsMax );

  
      // Second fit with tailored range
      // INT
      int lowBinCentInt = ( maxBinInt*binWidthInt ) - ( 4.5*parsInt[2] );
      int upBinCentInt = ( maxBinInt*binWidthInt ) + ( 4.5*parsInt[2] );
      landauInt[r][c]->SetParameters( parsInt[0], parsInt[1], parsInt[2] );
      landauInt[r][c]->SetRange( lowBinCentInt, upBinCentInt );
      gIntSpec[r][c]->Fit( landauInt[r][c], "+RQ" );
      landauInt[r][c]->GetParameters( parsInt );
      for ( int i=0; i<3; i++ ) parErrInt[i] = landauInt[r][c]->GetParError( i );
      // MAX
      int lowBinCentMax = ( maxBinMax*binWidthMax ) - ( 4.5*parsMax[2] );
      int upBinCentMax = ( maxBinMax*binWidthMax ) + ( 4.5*parsMax[2] );
      landauMax[r][c]->SetParameters( parsMax[0], parsMax[1], parsMax[2] );
      landauMax[r][c]->SetRange( lowBinCentMax, upBinCentMax );
      gMaxSpec[r][c]->Fit( landauMax[r][c], "+RQ" );
      landauMax[r][c]->GetParameters( parsMax );
      for ( int i=0; i<3; i++ ) parErrMax[i] = landauMax[r][c]->GetParError( i ); 

      // Count up the good events
      // INT
      double goodEvInt = gIntSpec[r][c]->Integral( gIntSpec[r][c]->FindFixBin( lowBinCentInt ), gIntSpec[r][c]->FindFixBin( upBinCentInt ), "" );
      // MAX
      double goodEvMax = gMaxSpec[r][c]->Integral( gMaxSpec[r][c]->FindFixBin( lowBinCentMax ), gMaxSpec[r][c]->FindFixBin( upBinCentMax ), "" );

      // Calculate target HV
      targetHV[r][c] = gHV[r][c]/pow(parsInt[1]/kTargetpC,1.0/gAlphas[r][c]);
      targetHVUnlimited[r][c] = targetHV[r][c];
      

      // Checking on goodness of fit for max spectra
      string Flag = "Good"; // Flag to keep track of fit status
      if( gMaxSpec[r][c]->GetEntries() == 0 ){
	cout << "**Module " << r << " " << c << " histo empty." << endl;
	Flag = "No_Data";
	for( int i=0; i<4; i++ ) { parsMax[i] = 0.0; parErrMax[i] = 0.0; }
	if( HVLimit == 1 ) targetHV[r][c] = gHV[r][c];
	goodEvMax = 0.0;
      }else if( parsMax[2] > 60.0 ){
	cout << "**Module " << r << " " << c << " fit too wide." << endl;
	cout << "parsMax[2] = " << parsMax[2] << endl;
	Flag = "Wide";
	for( int i=0; i<4; i++ ) { parsMax[i] = 0.0; parErrMax[i] = 0.0; }
	if( HVLimit == 1 ) targetHV[r][c] = gHV[r][c];
	goodEvMax = 0.0;
      }else if( parsMax[2] < 1.0 ){
	cout << "**Module " << r << " " << c << " fit too narrow." << endl;
	cout << "parsMax[2] = " << parsMax[2] << endl;
	Flag = "Narrow";
	for( int i=0; i<4; i++ ) { parsMax[i] = 0.0; parErrMax[i] = 0.0; }
	if( HVLimit == 1 ) targetHV[r][c] = gHV[r][c];
	goodEvMax = 0.0;
      }else if( parErrMax[1] > 20.0 || parErrMax[2] > 20.0 ){
	cout << "**Module " << r << " " << c << " error bar too high." << endl;
	cout << "parErrMax[1] = " << parErrMax[1] << " parErrMax[2] = " << parErrMax[2] << endl;
	Flag = "Big_error";
	for( int i=0; i<4; i++ ) { parsMax[i] = 0.0; parErrMax[i] = 0.0; }
	if( HVLimit == 1 ) targetHV[r][c] = gHV[r][c];
	goodEvMax = 0.0;
      }else if( targetHV[r][c]-gHV[r][c] > 100 ){
	if( HVLimit==1 ){
	  targetHV[r][c] = gHV[r][c]+100;
	  Flag = "Big_UP";
	}
      }else if( targetHV[r][c]-gHV[r][c] < -100 ){
	if( HVLimit==1 ){
	  targetHV[r][c] = gHV[r][c]-100;
	  Flag = "Big_DOWN";
	}
      }else if( parsMax[2] < 5.0 ){
	cout << "**Warning: Module " << r << " " << c << " appears narrow." << endl;
      }else if( targetHV[r][c] > 2500 ){
	targetHV[r][c] = 2500;
	Flag = "HV_MAX";
      }

      if( Flag != "Good" )
	cout << "Problem fit to module " << r << " " << c << " -> " << Flag << endl;

      // Write target HV and error flags to file
      cout << "Target HV for r=" << r << ", c=" << c << ", target pC=" << kTargetpC << ", HV setting=" << gHV[r][c] << ", mean ADC=" << parsInt[1] << ", and PMT alpha=" << gAlphas[r][c] << " is: " << targetHV[r][c] << "." << endl;
      reportFile << r << " " << c << " " << -targetHVUnlimited[r][c] << " " << -targetHV[r][c] << " " << -gHV[r][c] << " " << Flag << endl; 

      outFile << r*12+c+1 << '\t' << -targetHV[r][c] << endl;
      HVsetFile << -targetHV[r][c] << endl;

      // Write fitted histograms to file

      cout << "Writing fitted histograms to file.." << endl;
      
      int m = r*12+c;

      cosHistFile->cd();
      // INT
      gIntSpec[r][c]->SetTitle( Form( "Integrated Spectrum R%d C%d ADC%d FitMean %f", r, c, m, parsInt[1] ) );
      gIntSpec[r][c]->Write( Form( "ADCSpect_PMT%d", m ) );
      gIntSpec[r][c]->Draw( "AP" );
      // MAX
      gMaxSpec[r][c]->SetTitle( Form( "Max Spectrum R%d C%d FitMean %f", r, c, parsMax[1] ) );
      gMaxSpec[r][c]->Write( Form( "AmpSpect_PMT%d", m ) );
      gMaxSpec[r][c]->Draw( "AP" );
      // TDC
      gTDCSpec[r][c]->SetTitle( Form( "TDC Spectrum R%d C%d Mean %f", r, c, gTDCSpec[r][c]->GetMean() ) );
      gTDCSpec[r][c]->Write( Form( "TDCSpect_PMT%d", m ) );
      gTDCSpec[r][c]->Draw( "AP" );
      

      if( diagPlots==1 ){
	//cout << "Writing diagnostic plots to file.." << endl;
	gADCvChannel->SetBinContent( kNcols*r+c+1, parsInt[1] );
	gAmpvChannel->SetBinContent( kNcols*r+c+1, parsMax[1] );
	gPedvChannel->SetBinContent( kNcols*r+c+1, gPedSpec[r][c]->GetMean() );
	gNEVvChannel->SetBinContent( kNcols*r+c+1, gIntSpec[r][c]->GetEntries() );
	gTDCvChannel->SetBinContent( kNcols*r+c+1, gTDCSpec[r][c]->GetMean() );
	gTDCRMSvChannel->SetBinContent( kNcols*r+c+1, gTDCSpec[r][c]->GetRMS() );

	gNEV->SetBinContent( r+1, c+1, gIntSpec[r][c]->GetEntries() );
	gTHVvChannel->SetBinContent( kNcols*r+c, targetHV[r][c] );

	gPedSpec[r][c]->SetTitle( Form( "Pedestal Spectrum R%d C%d PMT%d", r, c, m ) );
	gPedSpec[r][c]->Write( Form( "PedSpect_PMT%d", m ) );
	gPedSpec[r][c]->Draw( "AP" );
      }
    }
  }
    
  if( diagPlots==1 ){
    cosHistFile->cd();
    gADCvChannel->SetTitle( "Average ADC vs Channel" );
    gADCvChannel->Write( "ADCvChannel" );
    gADCvChannel->Draw( "AP" );

    gAmpvChannel->SetTitle( "Average Amplitude vs Channel" );
    gAmpvChannel->Write( "AmpvChannel" );
    gAmpvChannel->Draw( "AP" );

    gPedvChannel->SetTitle( "Average Pedestal vs Channel" );
    gPedvChannel->Write( "PedvChannel" );
    gPedvChannel->Draw( "AP" );

    gNEVvChannel->SetTitle( "Number of Hits vs Channel" );
    gNEVvChannel->Write( "NEVvChannel" );
    gNEVvChannel->Draw( "AP" );

    gNEV->SetTitle( "Number of Hits Heatmap" );
    gNEV->Write( "NEV" );
    gNEV->Draw( "colz" );

    gTHVvChannel->SetTitle( "Target HV vs Channel" );
    gTHVvChannel->Write( "THVvChannel" );
    gTHVvChannel->Draw( "AP" );

    gTDCvChannel->SetTitle( "TDC Mean vs Channel" );
    gTDCvChannel->Write( "TDCvChannel" );
    gTDCvChannel->Draw( "AP" );

    gTDCRMSvChannel->SetTitle( "TDC RMS vs Channel" );
    gTDCRMSvChannel->Write( "TDCRMSvChannel" );
    gTDCRMSvChannel->Draw( "AP" );
  }

  int tab=4;
  TCanvas *c1[tab];

  if( diagPlots==1 ){
    cout << "Writing spectra to single canvas.." << endl;
    int j=0;
    int k=0;
    int m=0;

    for( int i=0; i<tab; i++){
      c1[i] = new TCanvas( Form("c1[%d]",i),"Integrated ADC Spectra (pC)",kCanvSize*kNrows,kCanvSize*kNrows);
      gStyle->SetOptStat(0);
      c1[i]->Divide(kNcols,kNrows/tab);
    }

    for( int r=0; r<kNrows; r++ ){
      for( int c=0; c<kNcols; c++ ){
	m = kNcols*r+c;
	k = m/72;
	j = m-72*k+1;	  

	cout << m << " " << j << " " << k << endl;

	c1[k]->cd(j);
	gIntSpec[r][c]->SetTitle(Form("R%d C%d",r,c));
	gIntSpec[r][c]->Draw();
      }
    }
    
    for( int i=0; i<tab; i++){
      gStyle->SetImageScaling(3.);
      c1[i]->SaveAs(Form("outFiles/Run%d_Q%d.png",run,i));
    }
    


  }

  reportFile.close();
  outFile.close();
  HVsetFile.close();

  // Post analysis reporting
  cout << "Finished loop over run " << run << "." << endl;
  cout << "Total good signals = " << gSignalTotal << "." << endl;
  cout << "Total saturated events = " << gSaturated << "." << endl;
  cout << "Target HV settings written to HVTargets_run" << run << ".txt." << endl;

  st->Stop();

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

  return 0;
}
  
