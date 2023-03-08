//SSeeds June 2021 - Cosmic PMT HV calibration code version 4 for use with SBS HCAL experiments. 24 rows x 12 cols HCAL modules. Relevant information can be found here: hallaweb.jlab.org/wiki/index.php/HCAL_Cosmics.

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

const int kMinSample = 0.0;
const int kMaxSample = 50.0;
const int kTotSample = ( kMaxSample-kMinSample ); // Should be the total fADC window size with each samp = 4ns
const double kTargetRAU = 78.355; // Previous value of 61.425*(2.55/2) or kTargetRAU_prev*(attLongCable_newMeas/attLongCable_oldMeas)
const int kMaxSteps = 10; // Corresponding to sample delta of 20 RAU (1 sigma LED) over 200 RAU (cosmics RAU range)
const double kMaxStep = 20; // Sample delta of 20 RAU (1 sigma LED on either side)

// Counter to keep track of T tree entry for processing
int gCurrentEntry = -1;
TChain *T = 0;

// Augment fn's for cosmic Calibration and recording
void goodHistoTest( TH1F*, double, int, int );
void getPedestal( int );
void diagnosticPlots( int, string );
string getDate( );

// Declare necessary histograms
TH1F *gHistos[kNrows][kNcols];
TH1F *gPedSpec[kNrows][kNcols];
TH1F *gRetToBase[kNrows][kNcols][kMaxSteps];
TH1F *gPMTIntSpec[kNrows][kNcols]; 
TH1F *gPMTMaxSpec[kNrows][kNcols]; 
TH1F *gPMTIntSpecTDC[kNrows][kNcols]; 
TH1F *gPMTMaxSpecTDC[kNrows][kNcols]; 

// Declare fit function
TF1 *g1;

// Create marker to track number of signals that pass cuts
int gSignalTotal = 1;

// Create files to keep temporary histograms for integrity checks
TFile *HistosFile = new TFile( "CosHistosFile.root", "RECREATE" );  // File for checking histogram fits and integrity
TFile *HistosFilePedestal = new TFile( "CosHistosFilePedestal.root", "RECREATE" );  // File for checking histogram fits and integrity

// Declare necessary arrays
double gPedestals[kNrows][kNcols];
double gPedSigma[kNrows][kNcols];
bool gSaturated[kNrows][kNcols];
double gPMTHV[kNrows][kNcols];
double gAlphas[kNrows][kNcols];
bool gPulse[kNrows+4][kNcols+4]; // Needs to be larger for false-valued buffers
bool gPulseTDC[kNrows+4][kNcols+4];
bool gVert[kNrows+4][kNcols+4];

double gMasterHisto[kNrows][kNcols];
int nevents = 0;

// Declare storage vectors
vector<double> gModule;
vector<double> gPeakPosInt, gPeakPosErrInt, gPeakPosMax, gPeakPosErrMax;
vector<double> gRMSInt, gRMSErrInt, gRMSMax, gRMSErrMax;
vector<double> gGoodEventMax, gGoodEventInt, gTargetHVOut;

// Create generic histogram function
TH1F* MakeHisto(int row, int col, int bins){
  TH1F *h = new TH1F(Form("h%02d%02d",row,col),Form("%d-%d",row+1,col+1),bins,kMinSample,kMaxSample);
  h->SetStats(0);
  h->SetLineWidth(2);
  return h;
}

// Make machine-based date function
string getDate(){
  time_t now = time(0);
  tm ltm = *localtime(&now);

  string yyyy = std::to_string(1900 + ltm.tm_year);
  string mm = std::to_string(1 + ltm.tm_mon);
  string dd = std::to_string(ltm.tm_mday);
  string date = mm + '/' + dd + '/' + yyyy;

  return date;
} // getDate

// Create diagnostic plots upon command
void diagnosticPlots( int run, string date ){
    char CName[9], CTitle[100];
    TCanvas *CGr[4];
    TGraph *Gr[4];

    CGr[0] = new TCanvas( "c1SH", "pPosvsModule", 100, 10, 700, 500 );
    CGr[1] = new TCanvas( "c2SH", "pRMSvsModule", 100, 10, 700, 500 );
    CGr[2] = new TCanvas( "c3SH", "#EvInPeakvsModule", 100, 10, 700, 500 );
    CGr[3] = new TCanvas( "c4SH", "HVcorrvsModule", 100, 10, 700, 500 );
 
    int totalModules = kNrows*kNcols;
    double xErr[totalModules]; // Holds x-error
    for( int i=0; i<totalModules; i++ ) xErr[i] = 0.0 ;
	 
    Gr[0] = new TGraphErrors( totalModules, &(gModule[0]), &(gPeakPosMax[0]), xErr, &(gPeakPosErrMax[0]) );
    Gr[1] = new TGraphErrors( totalModules, &(gModule[0]), &(gRMSMax[0]), xErr, &(gRMSErrMax[0]) );
    Gr[2] = new TGraph( totalModules, &(gModule[0]), &(gGoodEventMax[0]) );
    Gr[3] = new TGraph( totalModules, &(gModule[0]), &(gTargetHVOut[0]) );

    for( int i = 0; i < 4; i++ ){
      CGr[i]->cd();
      gPad->SetGridy();
      Gr[i]->SetLineColor( 2 );
      Gr[i]->SetLineWidth( 2 );
      Gr[i]->SetMarkerColor( 1 );
      Gr[i]->SetMarkerStyle( 20 );
      if( i == 0 ){
    	Gr[i]->SetTitle( Form( "Run# %d | Peak Position vs Module | %s", run, date.c_str() ) );
    	Gr[i]->GetXaxis()->SetTitle( "Module Number");
    	Gr[i]->GetYaxis()->SetTitle( "Peak Position (RAU)");
    	Gr[i]->GetYaxis()->SetRangeUser( 0.0, 200.0 );
      } else if( i == 1 ){
    	Gr[i]->SetTitle( Form( "Run# %d | Peak RMS vs Module | %s", run, date.c_str() ) ); 
    	Gr[i]->GetXaxis()->SetTitle( "Module Number" );
    	Gr[i]->GetYaxis()->SetTitle( "Peak RMS (RAU)" );
    	Gr[i]->GetYaxis()->SetTitleOffset( 1.4 );
    	Gr[i]->GetYaxis()->SetRangeUser( 0.0, 50.0 );
      }else if( i == 2 ){
    	Gr[i]->SetTitle( Form( "Run# %d | N of Events in Peak (fitted region) vs Module | %s", run, date.c_str() ) );
    	Gr[i]->GetXaxis()->SetTitle( "Module Number" );
    	Gr[i]->GetYaxis()->SetTitle( "N of Events in Peak (fitted region)" );
    	Gr[i]->GetYaxis()->SetTitleOffset( 1.4 );
      }else if( i == 3 ){
    	Gr[i]->SetTitle( Form( "Run# %d | HV Correction Factor vs Module | %s", run, date.c_str() ) );
    	Gr[i]->GetXaxis()->SetTitle( "Module Number" );
    	Gr[i]->GetYaxis()->SetTitle( "HV Correction Factor" );
    	Gr[i]->GetYaxis()->SetTitleOffset( 1.4 );
	Gr[i]->GetYaxis()->SetRangeUser( 0.4, 1.6 );
      }   
      Gr[i]->Draw( "AP" );
      CGr[i]->Write();
    }  

    CGr[0]->SaveAs( Form( "outFiles/HCal_Cosmic_%d.pdf[", run) );
    for( int i=0; i<4; i++ ) CGr[i]->SaveAs( Form( "outFiles/HCal_Cosmic_%d.pdf", run) );
    CGr[3]->SaveAs( Form( "outFiles/HCal_Cosmic_%d.pdf]", run) );
}

void processEvent(int entry = -1){
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
  
  int r,c,idx,n,sub;
  // Clear old signal histograms, just in case modules are not in the tree
  for( r = 0; r < kNrows; r++ ) {
    for( c = 0; c < kNcols; c++ ) {
      gHistos[r][c]->Reset( "ICES M" );
      gSaturated[r][c] = false;
      gPulse[r+2][c+2] = false;
      gPulseTDC[r+2][c+2] = false;
      gVert[r+2][c+2] = false;
    }
  }

  // Reset signal peak, adc, and tdc arrays
  float peak[kNrows][kNcols];
  double adc[kNrows][kNcols];
  double tdc[kNrows][kNcols];
  for( r = 0; r < kNrows; r++ ) {
    for( c = 0; c < kNcols; c++ ) {
      peak[r][c] = 0.0;
      adc[r][c] = 0.0;
      tdc[r][c] = 0.0;
    }
  }

  // Process event with m data
  for( int m = 0; m < hcalt::ndata; m++ ) {
    // Define row and column
    r = hcalt::row[m];
    c = hcalt::col[m];
    if( r < 0 || c < 0 ) {
      cerr << "Error: R or C negative." << endl;
      continue;
    }
    
    if( r >= kNrows || c >= kNcols ) continue;

    // Define index, number of samples, fill adc and tdc arrays, and switch processed marker for error reporting
    idx = hcalt::samps_idx[m];
    n = hcalt::nsamps[m];
    adc[r][c] = hcalt::a[m];
    tdc[r][c] = hcalt::tdc[m];
    bool processed = false;

    // Fill signal histogram from samps and mark saturated array if applicable
    for( int s = kMinSample; s < kMaxSample && s < n; s++ ) {
      processed = true;
      gHistos[r][c]->SetBinContent( s+1-kMinSample, hcalt::samps[idx+s] );

      //double lastSample = gMasterHisto[r][c]->GetBinContent( s+1-kMinSample );
      gMasterHisto[r][c]->Fill( s+1-kMinSample, hcalt::samps[idx+s] );

      if( peak[r][c] < hcalt::samps[idx+s] )
        peak[r][c] = hcalt::samps[idx+s];
      if( peak[r][c] > 4095 ) {
        gSaturated[r][c] = true;
      }
    }

    nevents++;
    // Report error if module is empty
    if(!processed) {
      cerr << "Skipping empty module: " << m << ".." << endl;
      for( int s = kMinSample; s < kTotSample; s++ ) {
        gHistos[r][c]->SetBinContent( s+1 - kMinSample, -404 );
      }
    }
  }

  // Pass unsaturated signal histos to goodHistoTest to see if a signal pulse exists there. gPulse array value  marked true if test passed.
  for( r = 0; r < kNrows; r++ ) {
    for( c = 0; c < kNcols; c++ ) {
      gHistos[r][c]->SetTitle( Form( "%d-%d (ADC=%g,TDC=%g)", r+1, c+1, adc[r][c], tdc[r][c] ) );
      if( !gSaturated[r][c] ){
	goodHistoTest( gHistos[r][c], tdc[r][c], r, c );
      }
    }
  }

  // Vertical line test - remember to shift all values up by one due to intialization of gPulse/gPulseTDC and 'false' buffer. Only building functionality for straight vertical line test.
  for( r = 0; r < kNrows; r++ ) {
    for( c = 0; c < kNcols; c++ ) {
      if( gPulse[r+2][c+2] == true ){
	if( ( gPulse[r][c+2] == true && gPulse[r+1][c+2] == true ) ||
	    ( gPulse[r+1][c+2] == true && gPulse[r+3][c+2] == true) ||
	    ( gPulse[r+3][c+2] == true && gPulse[r+4][c+2] ==true ) ) { // Checks if two pulses exist above, below, or one above and below every given pulse to ensure a track exists.
	  if( gPulse[r+2][c+1] == false && gPulse[r+2][c+3] == false ) // Check on each side of track for pulses to exclude events which clip the corner of the module.
	    gVert[r+2][c+2]=true; //Flip the bool to determine a good cosmic pulse.
	} else { // Else excluded. Diagonal tracks also excluded.
	  continue;
	}
      }
    }
  }

  // Now, if pulse passes verticality test, pedestal subtract and fill spectra histograms
  for( r = 0; r < kNrows; r++ ) {
    for( c = 0; c < kNcols; c++ ) {
      if( gVert[r+2][c+2] == true ) { 	
	// Get return to baseline parameter (in bins (ADC samples, 4ns per sample))
	double baseline = 0.0;
	double maxValue = 0.0;
	for( int s = 0; s < 4; s++ ){ // Get baseline from first four bins 
	  baseline += gHistos[r][c]->GetBinContent(s+1-kMinSample);  
	}
	baseline /= 4.0; // Get average value of the baseline per event
	maxValue = gHistos[r][c]->GetBinContent( gHistos[r][c]->GetMaximumBin() );
	for( int i = 0; i < kMaxSteps; i++){
	  if( maxValue > i*kMaxStep && maxValue < (i+1)*kMaxStep){
	    for( int s = gHistos[r][c]->GetMaximumBin(); s < kTotSample; s++ ){
	      double binContent = gHistos[r][c]->GetBinContent( s );
	      if( binContent < baseline+3*gPedSigma[r][c] ){
		gRetToBase[r][c][i]->Fill( s-gHistos[r][c]->GetMaximumBin() );
		break;
	      }else if( s == kTotSample-1 ){
		gRetToBase[r][c][i]->Fill( 50 ); // If there is never a return, use 50
	      }else{
		continue;
	      }
	    }
	  }
	}
	
	gPMTIntSpec[r][c]->Fill( gHistos[r][c]->Integral( 1, kMaxSample ) ); // Integral from total bin content
	gPMTMaxSpec[r][c]->Fill( gHistos[r][c]->GetMaximum() );
	
	// Only fill TDC if both ADC and TDC pulses present
	if( gPulseTDC[r+2][c+2] == true ){
	  gPMTIntSpecTDC[r][c]->Fill( gHistos[r][c]->Integral( 1, kMaxSample ) ); // Integral from total bin content
	  gPMTMaxSpecTDC[r][c]->Fill( gHistos[r][c]->GetMaximum() );
	}	
      }	
    }
  }
}

// Acquire pedestal for each module on events where no cosmic passed through said module (TDC==0). See comments for processEvent().
void getPedestal( int entry=-1 ){ 
  if( entry == -1 ) {
    gCurrentEntry++;
  } else {
    gCurrentEntry = entry;
  }
  
  if( gCurrentEntry < 0 ) {
    gCurrentEntry = 0;
  }
  
  T->GetEntry( gCurrentEntry );
  
  int r,c,idx,n,sub;
  
  for( r = 0; r < kNrows; r++ ) {
    for( c = 0; c < kNcols; c++ ) {
      gHistos[r][c]->Reset( "ICES M" );
    }
  }

  double adc[kNrows][kNcols];
  double tdc[kNrows][kNcols];
  for( r = 0; r < kNrows; r++ ) {
    for( c = 0; c < kNcols; c++ ) {
      adc[r][c] = 0.0;
      tdc[r][c] = 0.0;
    }
  }
  
  for( int m = 0; m < hcalt::ndata; m++ ) {
    r = hcalt::row[m];
    c = hcalt::col[m];
    if( r < 0 || c < 0 ) {
      cerr << "Error: R or C negative." << endl;
      continue;
    }
    
    if( r >= kNrows || c >= kNcols ) continue;
    
    idx = hcalt::samps_idx[m];
    n = hcalt::nsamps[m];
    adc[r][c] = hcalt::a[m];
    tdc[r][c] = hcalt::tdc[m];
    bool processed = false;
    for( int s = kMinSample; s < kMaxSample && s < n; s++ ) {
      processed = true;
      gHistos[r][c]->SetBinContent( s+1-kMinSample, hcalt::samps[idx+s] );

      //cout << s << " " << idx << " " << hcalt::samps[idx+s] << endl;
    }
  
    if( !processed ) {
      cerr << "Skipping empty module: " << m << endl;
      for( int s = kMinSample;  s < kTotSample; s++ ) {
        gHistos[r][c]->SetBinContent( s+1-kMinSample, 0 );
      }
    }
  }
  
  // With a histogram populated for each module on this event, fill pedestal histogram per module over all events with contents
  for( r = 0; r < kNrows; r++ ) {
    for( c = 0; c < kNcols; c++ ) {
      //if( tdc[r][c] == 0 ){ // Eliminating adc pulses from coincident tdc measurement in pedestal calculation	
	for( int b=0 ; b<kMaxSample ; b++) {	  
	  gPedSpec[r][c]->Fill( gHistos[r][c]->GetBinContent(b) );


	  //cout << gHistos[r][c]->GetBinContent(b) << endl;;

	}
	// }
    }
  }
}

// Pedestal subtract then verify that max value is greater than twice sigma of pedestals and lies w/in bounds of histo
void goodHistoTest( TH1F *testHisto, double tdcVal, int row, int col ){
  
  // Get pedestal value for the run and pmt 
  double pedVal = gPedestals[row][col];
  
  // Subtract pedestal value from each bin
  for( int b=1 ; b<=testHisto->GetNbinsX() ; b++ ) {
    testHisto->SetBinContent( b, testHisto->GetBinContent(b)-pedVal );
  }

  // Primary signal cut. All defined signals are greater than twice the sigma of the pedestal and the maximum value is not in the first or last bin
  if( testHisto->GetMaximum() > (6.0*gPedSigma[row][col]) &&
      testHisto->GetMaximumBin() > 2 &&
      testHisto->GetMaximumBin() < 28){

    // Switch if both ADC signal and TDC signal exists
    if( tdcVal != 0 ) gPulseTDC[row+2][col+2]=true;

    // Switch if ADC signal exists
    gPulse[row+2][col+2]=true;

    // Write out sample event histograms for independent verification
    if ( gSignalTotal % 50000 == 0 ){
      cout << "Writing a reference histogram to file.." << endl;
      cout << "Pedestal for same histogram = " << gPedestals[row][col] << "." << endl;
      cout << "Pedestal sigma for same histogram = " << gPedSigma[row][col] << "." << endl;
      cout << "Maximum value for same histogram = " << testHisto->GetMaximum() << "." << endl;
      cout << "TDC value for same histogram = " << tdcVal << "." << endl << endl;
      
      HistosFile->cd();
      testHisto->SetTitle( Form( "Sample R%d C%d", row, col ) );
      testHisto->GetYaxis()->SetTitle( "RAU" );
      testHisto->GetYaxis()->CenterTitle();
      testHisto->GetXaxis()->SetTitle( "ADC Sample" );
      testHisto->GetXaxis()->CenterTitle();
      testHisto->Write( Form( "sampleHisto%d_r%d_c%d", gSignalTotal, row, col ) );
      testHisto->Draw( "AP" );
    }
    gSignalTotal++; 
  }
}

// Main
int cosCal_v5(int run = 1725, int event = -1){

  // Initialize function with user commands
  bool diagPlots = 0;
  bool vertCut = 0;
  string date = getDate();

  cout << "Enter run number for analysis." << endl;
  cin >> run;
  cout << "Cut on vertical tracks? Enter 0 for NO and 1 for YES." << endl;
  cin >> vertCut;
  cout << "Print diagnostic plots? Enter 0 for NO and 1 for YES." << endl;
  cin >> diagPlots;
  
  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );
  
  //Declare outfile
  //TFile *cosHistFile = new TFile( Form( "outFiles/cosHistv4_run%d.root", run ), "RECREATE" );
  TFile *cosHistFile = new TFile( Form( "cosHistv4_DEBUG_run%d.root", run ), "RECREATE" ); //DEBUGGING
  ofstream fitData;
  fitData.open( Form( "outFiles/HCal_%d_FitResults.txt", run ) );
  fitData << "Run Number: " << run << " Desired Peak Position: " << kTargetRAU << endl;
  fitData << "Module " << " " << " Target HV " << " " << " Stat " << " " << " ErrStat " << " " << " Peak Pos " << " " << " ErrPPos " << " " << " Peak Width " << " " << " ErrPWid " << " " << " NGoodEvents " << " " <<  " " << " Flag " << std::endl;

  //Set spectrum histogram limits
  int PMTSpecBins = 350;
  //INT
  double PMTIntSpec_min = -30.0, PMTIntSpec_max = 670.0;
  //MAX
  double PMTMaxSpec_min = 0.0, PMTMaxSpec_max = 1050.0;
  
  //Build spectrum histograms. Empirical limits.
  cout << "Building ADC and TDC spectrum histograms.." << endl;
  for( int r=0; r<kNrows; r++ ){
    for( int c=0; c<kNcols; c++ ){  
      gPMTIntSpec[r][c] = new TH1F( Form( "Int ADC Spect R%d C%d", r, c ), Form( "Int ADC Spect R%d C%d", r, c ), PMTSpecBins, PMTIntSpec_min, PMTIntSpec_max );
      gPMTIntSpec[r][c]->GetXaxis()->SetTitle( "sRAU" );
      gPMTIntSpec[r][c]->GetXaxis()->CenterTitle();
      
      gPMTMaxSpec[r][c] = new TH1F( Form( "Max ADC Spect R%d C%d", r, c ), Form( "Max ADC Spect R%d C%d", r, c ), PMTSpecBins, PMTMaxSpec_min, PMTMaxSpec_max );
      gPMTMaxSpec[r][c]->GetXaxis()->SetTitle( "RAU" );
      gPMTMaxSpec[r][c]->GetXaxis()->CenterTitle();
      
      gPMTIntSpecTDC[r][c] = new TH1F( Form( "Int ADC Spect R%d C%d, TDC Cut", r, c ), Form( "Int ADC Spect R%d C%d", r, c ), PMTSpecBins, PMTIntSpec_min, PMTIntSpec_max );
      gPMTIntSpecTDC[r][c]->GetXaxis()->SetTitle( "sRAU" );
      gPMTIntSpecTDC[r][c]->GetXaxis()->CenterTitle();
      
      gPMTMaxSpecTDC[r][c] = new TH1F( Form( "Max ADC Spect R%d C%d, TDC Cut", r, c ), Form( "Max ADC Spect R%d C%d", r, c ), PMTSpecBins, PMTMaxSpec_min, PMTMaxSpec_max );
      gPMTMaxSpecTDC[r][c]->GetXaxis()->SetTitle( "RAU" );
      gPMTMaxSpecTDC[r][c]->GetXaxis()->CenterTitle();

      gPedSpec[r][c] = new TH1F( Form( "Pedestal Spect R%d C%d", r, c ), Form( "Pedestal Spect R%d C%d", r, c ), 200, -20, 100 ); // Empirical limits
      gPedSpec[r][c]->GetXaxis()->SetTitle( "<RAU>" );
      gPedSpec[r][c]->GetXaxis()->CenterTitle();;

      gMasterHisto[r][c] = new TH1F(Form("h%02d%02d",row,col),Form("%d-%d",row+1,col+1),kTotSample,kMinSample,kMaxSample);
      gMasterHisto[r][c]->GetXaxis()->SetTitle( "<RAU>" );
      gMasterHisto[r][c]->GetXaxis()->CenterTitle();;


      for( int i = 0; i < kMaxSteps; i++){
	gRetToBase[r][c][i] = new TH1F( Form( "Return to Baseline R%d C%d MaxADC(20*%d)", r, c, i ), Form( "Return to Baseline R%d C%d MaxADC(20*%d)", r, c, i ), kTotSample, kMinSample, kMaxSample );
	gRetToBase[r][c][i]->GetXaxis()->SetTitle( "Samples" );
	gRetToBase[r][c][i]->GetXaxis()->CenterTitle();
      }
    }
  }
    
  // Read in data produced by analyzer in root format
  cout << "Reading raw data from analyzer.." << endl;
  if( !T ) { 
    T = new TChain( "T" );
    //T->Add( Form( "$ROOTFILES/hcal_%d_*.root", run ) );
    T->Add( Form( "$ROOTFILES/hcal_%d*.root", run ) );
    
    //4T->Add( Form( "$ROOTFILES/hcal_%d*.root", run ) );
    
    //T->Add( Form( "rootFiles/cosmic/fadc_f1tdc_%d_1.root", run ) ); //DEBUGGING
    T->SetBranchStatus( "*", 0 );
    T->SetBranchStatus( "sbs.hcal.*", 1 );
    T->SetBranchAddress( "sbs.hcal.nsamps", hcalt::nsamps );
    T->SetBranchAddress( "sbs.hcal.a", hcalt::a );
    T->SetBranchAddress( "sbs.hcal.tdc", hcalt::tdc );
    T->SetBranchAddress( "sbs.hcal.samps", hcalt::samps );
    T->SetBranchAddress( "sbs.hcal.samps_idx", hcalt::samps_idx );
    T->SetBranchAddress( "sbs.hcal.adcrow", hcalt::row );
    T->SetBranchAddress( "sbs.hcal.adccol", hcalt::col );
    T->SetBranchStatus( "Ndata.sbs.hcal.adcrow", 1 );
    T->SetBranchAddress( "Ndata.sbs.hcal.adcrow", &hcalt::ndata );
    cout << "Opened up tree with nentries=" << T->GetEntries() << endl;
    for( int r = 0; r < kNrows; r++ ) {
      for( int c = 0; c < kNcols; c++ ) {
        gHistos[r][c] = MakeHisto( r, c, kTotSample );
        //gMasterHisto[r][c] = MakeHisto( r, c, kTotSample );
	gSaturated[r][c] = false;
      }
    }
  }

  // Set appropriate HV and alphas for run. HV settings from HCAL wiki. Alphas from LED analysis. Must have accompanying text file, one double per line by module number. Assuming negative voltage inputs.

  ifstream file( Form( "setFiles/HV_run%d.txt", run ) );

  cout << "Getting HV settings for each pmt for run " << run << "." << endl;

  int n1=0;
  double d1;
  
  int rval, cval;
  string line;
    
  while( getline( file, line ) ){
    if( line.at(0) == '#' ) {
      continue;
    }
    
    stringstream ss( line );
    ss >> d1;
    
    rval = floor( n1/kNcols );
    cval = n1 % kNcols;
    
    gPMTHV[rval][cval] = -d1; 
        
    n1++;
  }

  ifstream file2( "setFiles/alphas.txt" );
  
  n1=0;
  string line2;

  while( getline( file2, line2 ) ){
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
  
  // Set default values for pulse check bool arrays. Arrays contain false-valued buffer for possible future diagonal track implementation.
  cout << "Resetting particle track test flags.." << endl;
  for( int r = 0; r < kNrows+4; r++ ){
    for( int c = 0; c < kNcols+4; c++ ){
      gPulse[r][c] = false;
      gPulseTDC[r][c] = false;
      gVert[r][c] = false;
    }
  }
  
  gCurrentEntry = event;

  cout << "Total events to process " << T->GetEntries() << endl;

  // Fill pedestal histograms from data where f1TDC did not fire
  cout << "Filling pedestal histograms.." << endl;
  for ( int i = gCurrentEntry; i < T->GetEntries(); i++ ){ 
    getPedestal(gCurrentEntry);
    gCurrentEntry++;
    
    // Keep count of events processed for monitoring
    if ( gCurrentEntry%10000 == 0 ){
      cout << "Current pedestal event: " << gCurrentEntry << endl;
    }
  }

  /*
  // Write pedestal histograms to file
  cout << "Writing pedestal histograms to file.." << endl;
  for( int r=0; r<kNrows; r++ ){
    for( int c=0; c<kNcols; c++ ){  
      HistosFilePedestal->cd();
      gPedSpec[r][c]->SetTitle( Form( "Pedestal R%d C%d", r, c ) );
      gPedSpec[r][c]->Write( Form( "pedHisto_r%d_c%d", r, c ) );
      gPedSpec[r][c]->Draw( "AP" );
    }
  }
  */

  for( int r=0; r<kNrows; r++ ){
    for( int c=0; c<kNcols; c++ ){
      cout << r << " " << c << " " << gPedSpec[r][c]->GetEntries() << " " << gPedSpec[r][c] << endl;
    }
  }


  // Fitting pedestal histograms and getting mean to subtract event pulses by module
  cout << "Processing pedestals and saving to file.." << endl;
  for( int r=0; r<kNrows; r++ ){
    for( int c=0; c<kNcols; c++ ){

      int m = kNcols*r+c;

      gPedSpec[r][c]->Fit( "gaus", "Q" );
      g1=gPedSpec[r][c]->GetFunction( "gaus" );
      gPedestals[r][c]=g1->GetParameter( 1 );
      gPedSigma[r][c]=g1->GetParameter( 2 );
      HistosFilePedestal->cd();
      gPedSpec[r][c]->SetTitle(Form( "Pedestal R%d C%d M%d", r, c, m ) );
      gPedSpec[r][c]->Write( Form( "pedHisto_r%d_c%d_m%d", r, c, m ) );
      gPedSpec[r][c]->Draw( "AP" );
    }
  }

  gCurrentEntry = event; // Resetting entries for pulse analysis

  // Pedestal subtract each pulse by module and populate fADC and f1TDC spectra histograms
  cout << "Pedestal subtracting and processing signals.." << endl;
  for ( int i = gCurrentEntry; i < T->GetEntries(); i++ ){ 
    processEvent( gCurrentEntry );
    gCurrentEntry++;
    
    // Keep count of events processed for monitoring
    if ( gCurrentEntry%10000 == 0 ){
      cout << "Current event: " << gCurrentEntry << endl;
    }
  }

  // Reset vectors
  gModule.clear();
  gPeakPosInt.clear();
  gPeakPosErrInt.clear();
  gPeakPosMax.clear();
  gPeakPosErrMax.clear();
  gRMSInt.clear();
  gRMSErrInt.clear();
  gRMSMax.clear();
  gRMSErrMax.clear();
  gGoodEventMax.clear();
  gGoodEventInt.clear();
  gTargetHVOut.clear();

  // Fit ADC spectra with landau to extract mpv for HV calibration
  cout << "Writing spectrum histograms and calibration constants to file.." << endl;
  ofstream outFile;
  outFile.open( Form( "outFiles/HVTargetsV4_run%d.txt", run ) ); // Create text file to hold target HVs
  time_t now = time(0); 
  char *dt = ctime(&now);
  outFile << "#Target HV settings for run " << run << " from cosCal_v4 on " << dt << "#" << endl;
  outFile << "#Row Col targetHV" << endl;

  // Implement P. Datta's two-fit scheme
  TF1 *landauInt[kNrows][kNcols] = {};
  TF1 *landauMax[kNrows][kNcols] = {};

  // Fit parameter array for built-in landau
  double parsInt[3];
  double parErrInt[3];
  double parsMax[3];
  double parErrMax[3];

  // Array to store target high voltage per module as it's calculated
  double targetHV[kNrows][kNcols];

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
      int maxBinInt = gPMTIntSpec[r][c]->GetMaximumBin();
      double maxBinCenterInt = gPMTIntSpec[r][c]->GetXaxis()->GetBinCenter( maxBinInt );
      double maxCountInt = gPMTIntSpec[r][c]->GetMaximum();
      double binWidthInt = gPMTIntSpec[r][c]->GetBinWidth( maxBinInt );
      double stdDevInt = gPMTIntSpec[r][c]->GetStdDev();
      // MAX
      int maxBinMax = gPMTMaxSpec[r][c]->GetMaximumBin();
      double maxBinCenterMax = gPMTMaxSpec[r][c]->GetXaxis()->GetBinCenter( maxBinMax );
      double maxCountMax = gPMTMaxSpec[r][c]->GetMaximum();
      double binWidthMax = gPMTMaxSpec[r][c]->GetBinWidth( maxBinMax );
      double stdDevMax = gPMTMaxSpec[r][c]->GetStdDev();
      // INT
      landauInt[r][c]->SetParameters( maxCountInt, maxBinCenterInt, stdDevInt );
      landauInt[r][c]->SetRange( PMTIntSpec_min, PMTIntSpec_max );
      gPMTIntSpec[r][c]->Fit(landauInt[r][c],"No+RQ");
      landauInt[r][c]->GetParameters( parsInt );
      // MAX
      landauMax[r][c]->SetParameters( maxCountMax, maxBinCenterMax, stdDevMax );
      landauMax[r][c]->SetRange( PMTMaxSpec_min, PMTMaxSpec_max );
      gPMTMaxSpec[r][c]->Fit(landauMax[r][c],"No+RQ");
      landauMax[r][c]->GetParameters( parsMax );

  
      // Second fit with tailored range
      // INT
      int lowBinCentInt = ( maxBinInt*binWidthInt ) - ( 4.5*parsInt[2] );
      int upBinCentInt = ( maxBinInt*binWidthInt ) + ( 4.5*parsInt[2] );
      landauInt[r][c]->SetParameters( parsInt[0], parsInt[1], parsInt[2] );
      landauInt[r][c]->SetRange( lowBinCentInt, upBinCentInt );
      gPMTIntSpec[r][c]->Fit( landauInt[r][c], "+RQ" );
      landauInt[r][c]->GetParameters( parsInt );
      for ( int i=0; i<3; i++ ) parErrInt[i] = landauInt[r][c]->GetParError( i );
      // MAX
      int lowBinCentMax = ( maxBinMax*binWidthMax ) - ( 4.5*parsMax[2] );
      int upBinCentMax = ( maxBinMax*binWidthMax ) + ( 4.5*parsMax[2] );
      landauMax[r][c]->SetParameters( parsMax[0], parsMax[1], parsMax[2] );
      landauMax[r][c]->SetRange( lowBinCentMax, upBinCentMax );
      gPMTMaxSpec[r][c]->Fit( landauMax[r][c], "+RQ" );
      landauMax[r][c]->GetParameters( parsMax );
      for ( int i=0; i<3; i++ ) parErrMax[i] = landauMax[r][c]->GetParError( i ); 

      // Count up the good events
      // INT
      double goodEvInt = gPMTIntSpec[r][c]->Integral( gPMTIntSpec[r][c]->FindFixBin( lowBinCentInt ), gPMTIntSpec[r][c]->FindFixBin( upBinCentInt ), "" );
      // MAX
      double goodEvMax = gPMTMaxSpec[r][c]->Integral( gPMTMaxSpec[r][c]->FindFixBin( lowBinCentMax ), gPMTMaxSpec[r][c]->FindFixBin( upBinCentMax ), "" );

      // Calculate target HV
      targetHV[r][c] = gPMTHV[r][c]/pow(parsMax[1]/kTargetRAU,1.0/gAlphas[r][c]);
      
      // Checking on goodness of fit for max spectra
      string Flag = "Good"; // Flag to keep track of fit status
      if( gPMTMaxSpec[r][c]->GetEntries() == 0 ){
	cout << "**Module " << r << " " << c << " histo empty." << endl;
	Flag = "No_Data";
	for( int i=0; i<4; i++ ) { parsMax[i] = 0.0; parErrMax[i] = 0.0; }
	targetHV[r][c] = 1.0;
	goodEvMax = 0.0;
      }else if( parsMax[2] > 60.0 ){
	cout << "**Module " << r << " " << c << " fit too wide." << endl;
	cout << "parsMax[2] = " << parsMax[2] << endl;
	Flag = "Wide";
	for( int i=0; i<4; i++ ) { parsMax[i] = 0.0; parErrMax[i] = 0.0; }
	targetHV[r][c] = 1.0;
	goodEvMax = 0.0;
      }else if( parsMax[2] < 1.0 ){
	cout << "**Module " << r << " " << c << " fit too narrow." << endl;
	cout << "parsMax[2] = " << parsMax[2] << endl;
	Flag = "Narrow";
	for( int i=0; i<4; i++ ) { parsMax[i] = 0.0; parErrMax[i] = 0.0; }
	targetHV[r][c] = 1.0;
	goodEvMax = 0.0;
      }else if( parErrMax[1] > 20.0 || parErrMax[2] > 20.0 ){
	cout << "**Module " << r << " " << c << " error bar too high." << endl;
	cout << "parErrMax[1] = " << parErrMax[1] << " parErrMax[2] = " << parErrMax[2] << endl;
	Flag = "Big_error";
	for( int i=0; i<4; i++ ) { parsMax[i] = 0.0; parErrMax[i] = 0.0; }
	targetHV[r][c] = 1.0;
	goodEvMax = 0.0;
      }else if( parsMax[2] < 5.0 ){
	cout << "**Warning: Module " << r << " " << c << " appears narrow." << endl;
      }

      // Record fit parameters and target HV by module
      gModule.push_back( kNcols*r+c );
      gTargetHVOut.push_back( targetHV[r][c] );
      // INT
      gPeakPosInt.push_back( parsInt[1] );
      gPeakPosErrInt.push_back( parErrInt[1] );
      gRMSInt.push_back( parsInt[2] );
      gRMSErrInt.push_back( parErrInt[2] );
      gGoodEventInt.push_back( goodEvInt );
      // MAX
      gPeakPosMax.push_back( parsMax[1] );
      gPeakPosErrMax.push_back( parErrMax[1] );
      gRMSMax.push_back( parsMax[2] );
      gRMSErrMax.push_back( parErrMax[2] );
      gGoodEventMax.push_back( goodEvMax );

      
      // Write all the important fit parameters in a text file
      fitData.setf( ios::fixed );
      fitData.setf( ios::showpoint );
      fitData.precision( 2 );
      fitData.width( 5 ); fitData << kNcols*r+c; 
      fitData.width( 12 ); fitData << parsMax[0];
      fitData.width( 12 ); fitData << parErrMax[0];
      fitData.width( 12 ); fitData << parsMax[1]; 
      fitData.width( 12 ); fitData << parErrMax[1];
      fitData.width( 12 ); fitData << parsMax[2]; 
      fitData.width( 12 ); fitData << parErrMax[2];
      fitData.width( 12 ); fitData << goodEvMax;
      fitData.width( 12 ); fitData.precision( 4 ); fitData << targetHV[r][c];
      fitData.width( 12 ); fitData << Flag << endl;

      if( Flag != "Good" )
	cout << "Problem fit to module " << r << " " << c << " -> " << Flag << endl;

      // Write fitted histograms to file
      cosHistFile->cd();
      // INT
      gPMTIntSpec[r][c]->SetTitle( Form( "Int Spect R%d C%d FitMean%f", r, c, parsInt[1] ) );
      gPMTIntSpec[r][c]->Write( Form( "Int Spect R%d C%d", r, c ) );
      gPMTIntSpec[r][c]->Draw( "AP" );
      // MAX
      gPMTMaxSpec[r][c]->SetTitle( Form( "Max Spect R%d C%d FitMean%f", r, c, parsMax[1] ) );
      gPMTMaxSpec[r][c]->Write( Form( "Max Spect R%d C%d", r, c ) );
      gPMTMaxSpec[r][c]->Draw( "AP" );


      gMasterHisto[r][c]->SetTitle( Form( "Master Histo R%d C%d", r, c) );
      gMasterHisto[r][c]->Write( Form( "Master Histo R%d C%d", r, c ) );
      gMasterHisto[r][c]->Draw( "AP" );

      

      // Draw return to baseline histograms
      for( int i = 0; i < kMaxSteps; i++){
	gRetToBase[r][c][i]->SetTitle( Form( "Return to Baseline R%d C%d MaxADC(%d*20)", r, c, i ) );
	gRetToBase[r][c][i]->Write(Form("RTB_R%d_C%d_DSamp(%d*20)", r, c, i ) );
	gRetToBase[r][c][i]->Draw("AP");
      }
    }
  }
    
  // Post analysis reporting
  cout << "Finished loop over run " << run << "." << endl;
  cout << "Total good signals = " << gSignalTotal << "." << endl;
  cout << "Sample histograms drawn to file HistosFile.root." << endl;
  cout << "Target HV settings written to HVTargetsV3_run" << run << ".txt." << endl;
  if( diagPlots ){
    cout << "Printing diagnostic plots.." << endl;
    diagnosticPlots( run, date );
  }

  st->Stop();

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

  return 0;
}
  
