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

const double kTargetRAU = 78.355; // Previous value of 61.425*(2.55/2) or kTargetRAU_prev*(attLongCable_newMeas/attLongCable_oldMeas)
const int kMaxSteps = 10; // Corresponding to sample delta of 20 RAU (1 sigma LED) over 200 RAU (cosmics RAU range)
const double kMaxStep = 20; // Sample delta of 20 RAU (1 sigma LED on either side)

// Counter to keep track of T tree entry for processing
int gCurrentEntry = -1;
TChain *T = 0;

// Augment fn's for cosmic Calibration and recording
//void goodHistoTest( TH1F*, double, int, int );
//void getPedestal( int );
void diagnosticPlots( int, string );
string getDate( );

// Declare necessary histograms
//TH1F *gHistos[kNrows][kNcols];
//TH1F *gPedSpec[kNrows][kNcols];
//TH1F *gRetToBase[kNrows][kNcols][kMaxSteps];
TH1F *gPMTIntSpec[kNrows][kNcols]; 
TH1F *gPMTMaxSpec[kNrows][kNcols]; 
TH1F *gPMTIntSpecTDC[kNrows][kNcols]; 
TH1F *gPMTMaxSpecTDC[kNrows][kNcols]; 

// Declare fit function
TF1 *g1;

// Create marker to track number of signals that pass cuts
int gSignalTotal = 0;

// Create files to keep temporary histograms for integrity checks
//TFile *HistosFile = new TFile( "outFiles/cosHistosFile.root", "RECREATE" );  // File for checking histogram fits and integrity
//TFile *HistosFilePedestal = new TFile( "outFiles/CosHistosFilePedestal.root", "RECREATE" );  // File for checking histogram fits and integrity

// Declare necessary arrays
//double gPedestals[kNrows][kNcols];
//double gPedSigma[kNrows][kNcols];
//bool gSaturated[kNrows][kNcols];
double gPMTHV[kNrows][kNcols];
double gAlphas[kNrows][kNcols];
//bool gPulse[kNrows+4][kNcols+4]; // Needs to be larger for false-valued buffers
//bool gPulseTDC[kNrows+4][kNcols+4];
//bool gVert[kNrows+4][kNcols+4];

// Declare storage vectors
//vector<double> gModule;
//vector<double> gPeakPosInt, gPeakPosErrInt, gPeakPosMax, gPeakPosErrMax;
//vector<double> gRMSInt, gRMSErrInt, gRMSMax, gRMSErrMax;
//vector<double> gGoodEventMax, gGoodEventInt, gTargetHVOut;

/* Not yet sure if this is necessary
// Create generic histogram function
TH1F* MakeHisto(int row, int col, int bins){
  TH1F *h = new TH1F(Form("h%02d%02d",row,col),Form("%d-%d",row+1,col+1),bins,kMinSample,kMaxSample);
  h->SetStats(0);
  h->SetLineWidth(2);
  return h;
}
*/

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
/*
// Not yet necessary - will parse
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

*/

void processEvent(int entry = -1, double cut = 5){
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
  
  //int r,c,idx,n,sub;

  int r,c;
  
  /*
  // Clear old signal histograms, just in case modules are not in the tree
  for( r = 0; r < kNrows; r++ ) {
    for( c = 0; c < kNcols; c++ ) {
      //gHistos[r][c]->Reset( "ICES M" );
      gSaturated[r][c] = false;
      //gPulse[r+2][c+2] = false;
      //gPulseTDC[r+2][c+2] = false;
      //gVert[r+2][c+2] = false;
    }
  }
  */
  // Declare signal peak, adc, and tdc arrays for this event
  //float peak[kNrows][kNcols];
  double adc[kNrows][kNcols] = {0.0};
  double tdc[kNrows][kNcols] = {0.0};
  double adc_p[kNrows][kNcols] = {0.0};
  //double ped[kNrows][kNcols];
  double amp[kNrows][kNcols] = {0.0};
  double amp_p[kNrows][kNcols] = {0.0};
  bool saturated[kNrows][kNcols] = { {false} };
  /*
  for( r = 0; r < kNrows; r++ ) {
    for( c = 0; c < kNcols; c++ ) {
      //peak[r][c] = 0.0;
      adc[r][c] = 0.0;
      tdc[r][c] = 0.0;
    }
  }
  */
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

    // Define index, number of samples, fill adc and tdc arrays, and switch processed marker for error reporting
    //idx = hcalt::samps_idx[m];
    //n = hcalt::nsamps[m];
    adc[r][c] = hcalt::a[m];
    adc_p[r][c] = hcalt::a_p[m];
    amp[r][c] = hcalt::amp[m];
    amp_p[r][c] = hcalt::amp_p[m];
    //ped[r][c] = hcalt::ped[m];
    tdc[r][c] = hcalt::tdc[m];
    
    //bool processed = false;

    /*
    // Fill signal histogram from samps and mark saturated array if applicable
    for( int s = kMinSample; s < kMaxSample && s < n; s++ ) {
      processed = true;

      gHistos[r][c]->SetBinContent( s+1-kMinSample, hcalt::samps[idx+s] ); //No need to fill gHistos
      if( peak[r][c] < hcalt::samps[idx+s] )
        peak[r][c] = hcalt::samps[idx+s];
      if( peak[r][c] > 4095 ) {
        gSaturated[r][c] = true;
      }

    }
    */

    // Mark saturated array when amplitude meets max RAU
    if( amp[r][c] > 4094 ) {
      saturated[r][c] = true;
    }
    /*
    // Report error if module is empty
    if(!processed) {
      cerr << "Missing data on event " << entry << ", module " << m << ".." << endl;
      
      for( int s = kMinSample; s < kTotSample; s++ ) {
        gHistos[r][c]->SetBinContent( s+1 - kMinSample, -404 );
      }
      
      }
    */
  }
  /*
  // Pass unsaturated signal histos to goodHistoTest to see if a signal pulse exists there. gPulse array value  marked true if test passed.
  for( r = 0; r < kNrows; r++ ) {
    for( c = 0; c < kNcols; c++ ) {
      //gHistos[r][c]->SetTitle( Form( "%d-%d (ADC=%g,TDC=%g)", r+1, c+1, adc[r][c], tdc[r][c] ) );
      if( !gSaturated[r][c] ){
	goodHistoTest( gHistos[r][c], tdc[r][c], r, c );
      }
    }
  }
  */

  // Assuming pedestal sigma reasonable for all channels all events (sigma << 5 mV or 10 RAU (observed roughly 1 mV or 2 RAU))
  /*
  // Vertical line test - remember to shift all values up by one due to intialization of gPulse/gPulseTDC and 'false' buffer. Only building functionality for straight vertical line test.
  */
  // Vertical line test - accept only if three hits in vertical track. Reject if hit adjacent in row.
  for( r = 0; r < kNrows; r++ ) {
    for( c = 0; c < kNcols; c++ ) {
      if( amp_p[r][c]>cut && saturated[r][c]==false ) { // Only examine channels with max > 5 sigma pedestal for hits and no saturation
	if( r==0 && c==0 ) { // Logic for top left channel
	  if( amp_p[r+1][c] > cut && amp_p[r+2][c] > cut ){
	    if( amp_p[r][c+1] < 6 ) {
	      if( tdc[r][c] != 0 ) {
		gPMTIntSpecTDC[r][c]->Fill( adc_p[r][c] );
		gPMTMaxSpecTDC[r][c]->Fill( amp_p[r][c] );
	      }
	      gPMTIntSpec[r][c]->Fill( adc_p[r][c] );
	      gPMTMaxSpec[r][c]->Fill( amp_p[r][c] );
	      gSignalTotal++;
	    }
	  }
	}
	else if( r==0 && c==11 ) { // Logic for top right channel
	  if( amp_p[r+1][c] > cut && amp_p[r+2][c] > cut ){
	    if( amp_p[r][c-1] < 6 ) {
	      if( tdc[r][c] != 0 ) {
		gPMTIntSpecTDC[r][c]->Fill( adc_p[r][c] );
		gPMTMaxSpecTDC[r][c]->Fill( amp_p[r][c] );
	      }
	      gPMTIntSpec[r][c]->Fill( adc_p[r][c] );
	      gPMTMaxSpec[r][c]->Fill( amp_p[r][c] );
	      gSignalTotal++;
	    }
	  }
	}
	else if( r==23 && c==0 ) { // Logic for bottom left channel
	  if( amp_p[r-1][c] > cut && amp_p[r-2][c] > cut ){
	    if( amp_p[r][c+1] < 6 ) {
	      if( tdc[r][c] != 0 ) {
		gPMTIntSpecTDC[r][c]->Fill( adc_p[r][c] );
		gPMTMaxSpecTDC[r][c]->Fill( amp_p[r][c] );
	      }
	      gPMTIntSpec[r][c]->Fill( adc_p[r][c] );
	      gPMTMaxSpec[r][c]->Fill( amp_p[r][c] );
	      gSignalTotal++;
	    }
	  }
	}
	else if( r==23 && c==11 ) { // Logic for bottom right channel
	  if( amp_p[r-1][c] > cut && amp_p[r-2][c] > cut ){
	    if( amp_p[r][c-1] < 6 ) {
	      if( tdc[r][c] != 0 ) {
		gPMTIntSpecTDC[r][c]->Fill( adc_p[r][c] );
		gPMTMaxSpecTDC[r][c]->Fill( amp_p[r][c] );
	      }
	      gPMTIntSpec[r][c]->Fill( adc_p[r][c] );
	      gPMTMaxSpec[r][c]->Fill( amp_p[r][c] );
	      gSignalTotal++;
	    }
	  }
	}
	else if( r==0 || r==1) { // Logic for other row 1 or 2 channels
	  if( (amp_p[r+1][c] > cut && amp_p[r+2][c] > cut) ||
	      (r==1 && amp_p[r-1][c] > cut && amp_p[r+1][c] > cut) ) {
	    if( amp_p[r][c-1] < 6 && amp_p[r][c+1] < 6 ) {
	      if( tdc[r][c] != 0 ) {
		gPMTIntSpecTDC[r][c]->Fill( adc_p[r][c] );
		gPMTMaxSpecTDC[r][c]->Fill( amp_p[r][c] );
	      }
	      gPMTIntSpec[r][c]->Fill( adc_p[r][c] );
	      gPMTMaxSpec[r][c]->Fill( amp_p[r][c] );
	      gSignalTotal++;
	    }
	  }
	}
	else if( r==23 || r==22 ) { // Logic for other row 22 or 23 channels
	  if( (amp_p[r-1][c] > cut && amp_p[r-2][c] > cut) ||
	      (r==22 && amp_p[r-1][c] > cut && amp_p[r+1][c] > cut) ) {
	    if( amp_p[r][c-1] < 6 && amp_p[r][c+1] < 6 ) {
	      if( tdc[r][c] != 0 ) {
		gPMTIntSpecTDC[r][c]->Fill( adc_p[r][c] );
		gPMTMaxSpecTDC[r][c]->Fill( amp_p[r][c] );
	      }
	      gPMTIntSpec[r][c]->Fill( adc_p[r][c] );
	      gPMTMaxSpec[r][c]->Fill( amp_p[r][c] );
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
		gPMTIntSpecTDC[r][c]->Fill( adc_p[r][c] );
		gPMTMaxSpecTDC[r][c]->Fill( amp_p[r][c] );
	      }
	      gPMTIntSpec[r][c]->Fill( adc_p[r][c] );
	      gPMTMaxSpec[r][c]->Fill( amp_p[r][c] );
	      gSignalTotal++;
	    }
	  }
	}
	/*
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
	*/


      }
    }
  }

/*
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
*/
/*
//SBS-Offline now processes this - no further analysis needed
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
	//cout << gHistos[r][c]->GetBinContent(b+1) << endl;
	gPedSpec[r][c]->Fill( gHistos[r][c]->GetBinContent(b+1) );
      }
      //}
    }
  }
}
*/

/*
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
*/
}

// Main
int hcal_gain_match(int run = -1, int event = -1){
    
  // Initialize function with user commands
  bool diagPlots = 0;
  bool vertCut = 0;
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
  
  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );
  
  //Declare outfile
  TFile *cosHistFile = new TFile( Form( "outFiles/gainMatchResults_run%d.root", run ), "RECREATE" );
  ofstream fitData;
  fitData.open( Form( "outFiles/HCal_%d_FitData.txt", run ) );
  fitData << "Run Number: " << run << " Desired Peak Position: " << kTargetRAU << endl;
  fitData << "Module " << " " << " Target HV " << " " << " Stat " << " " << " ErrStat " << " " << " Peak Pos " << " " << " ErrPPos " << " " << " Peak Width " << " " << " ErrPWid " << " " << " NGoodEvents " << " " <<  " " << " Flag " << std::endl;

  //Set spectrum histogram limits
  int PMTSpecBins = 100;
  //INT
  double PMTIntSpec_min = 0.0, PMTIntSpec_max = 4000.0;
  //MAX
  double PMTMaxSpec_min = 0.0, PMTMaxSpec_max = 1000.0;
  
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
      /*
      gPedSpec[r][c] = new TH1F( Form( "Pedestal Spect R%d C%d", r, c ), Form( "Pedestal Spect R%d C%d", r, c ), 100, -5, 5 ); // Empirical limits
      gPedSpec[r][c]->GetXaxis()->SetTitle( "<RAU>" );
      gPedSpec[r][c]->GetXaxis()->CenterTitle();
      */
      /*
      for( int i = 0; i < kMaxSteps; i++){
	gRetToBase[r][c][i] = new TH1F( Form( "Return to Baseline R%d C%d MaxADC(20*%d)", r, c, i ), Form( "Return to Baseline R%d C%d MaxADC(20*%d)", r, c, i ), kTotSample, kMinSample, kMaxSample );
	gRetToBase[r][c][i]->GetXaxis()->SetTitle( "Samples" );
	gRetToBase[r][c][i]->GetXaxis()->CenterTitle();
      }
      */

    }
  }
  
  // Read in data produced by analyzer in root format
  cout << "Reading replayed root file.." << endl;
  if( !T ) { 
    T = new TChain( "T" );
    T->Add( Form( "/volatile/halla/sbs/seeds/rootFiles/hcal_%d_*.root", run ) );
    T->SetBranchStatus( "*", 0 );
    T->SetBranchStatus( "sbs.hcal.*", 1 );
    //T->SetBranchAddress( "sbs.hcal.nsamps", hcalt::nsamps );
    T->SetBranchAddress( "sbs.hcal.a", hcalt::a );
    T->SetBranchAddress( "sbs.hcal.a_amp", hcalt::amp );
    T->SetBranchAddress( "sbs.hcal.a_p", hcalt::a_p );
    T->SetBranchAddress( "sbs.hcal.a_amp_p", hcalt::amp );
    //T->SetBranchAddress( "sbs.hcal.ped", hcalt::ped );
    T->SetBranchAddress( "sbs.hcal.tdc", hcalt::tdc );
    //T->SetBranchAddress( "sbs.hcal.samps", hcalt::samps );
    //T->SetBranchAddress( "sbs.hcal.samps_idx", hcalt::samps_idx );
    T->SetBranchAddress( "sbs.hcal.adcrow", hcalt::row );
    T->SetBranchAddress( "sbs.hcal.adccol", hcalt::col );
    T->SetBranchStatus( "Ndata.sbs.hcal.adcrow", 1 );
    T->SetBranchAddress( "Ndata.sbs.hcal.adcrow", &hcalt::ndata );
    cout << "Opened up tree with nentries=" << T->GetEntries() << endl;
    /*
      for( int r = 0; r < kNrows; r++ ) {
      for( int c = 0; c < kNcols; c++ ) {
      gHistos[r][c] = MakeHisto( r, c, kTotSample );
      gSaturated[r][c] = false;
      }
    */
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
    
    gPMTHV[rval][cval] = -d1; 
        
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
  /*
  // Set default values for pulse check bool arrays. Arrays contain false-valued buffer for possible future diagonal track implementation.
  cout << "Declaring particle track test flags.." << endl;
  for( int r = 0; r < kNrows+4; r++ ){
    for( int c = 0; c < kNcols+4; c++ ){
      gPulse[r][c] = false;
      gPulseTDC[r][c] = false;
      gVert[r][c] = false;
    }
  }
  */
  gCurrentEntry = event;

  cout << "Total events to process " << T->GetEntries() << ".." << endl;
  /*
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

  // Fitting pedestal histograms and getting mean to subtract event pulses by module
  cout << "Processing pedestals and saving to file.." << endl;
  for( int r=0; r<kNrows; r++ ){
    for( int c=0; c<kNcols; c++ ){
      gPedSpec[r][c]->Fit( "gaus", "Q" );
      g1=gPedSpec[r][c]->GetFunction( "gaus" );
      gPedestals[r][c]=g1->GetParameter( 1 );
      gPedSigma[r][c]=g1->GetParameter( 2 );
      HistosFilePedestal->cd();
      gPedSpec[r][c]->SetTitle(Form( "Pedestal R%d C%d", r, c ) );
      gPedSpec[r][c]->Write( Form( "pedHisto_r%d_c%d", r, c ) );
      gPedSpec[r][c]->Draw( "AP" );
    }
  }

  gCurrentEntry = event; // Resetting entries for pulse analysis
  */

  // Pedestal subtract each pulse by module and populate fADC and f1TDC spectra histograms
  //cout << "Pedestal subtracting and processing signals.." << endl;

  cout << "Processing hits by event.." << endl;

  for ( int i = gCurrentEntry; i < T->GetEntries(); i++ ){ 
    processEvent( gCurrentEntry, adcCut );
    gCurrentEntry++;
    
    // Keep count of events processed for monitoring
    if ( gCurrentEntry%10000 == 0 ){
      cout << "Current event: " << gCurrentEntry << endl;
    }
  }
  /*
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
  */
  // Fit ADC spectra with landau to extract mpv for HV calibration
  cout << "Writing spectrum histograms and calibration constants to file.." << endl;
  ofstream outFile;
  outFile.open( Form( "outFiles/HVTargets_run%d.txt", run ) ); // Create text file to hold target HVs
  time_t now = time(0); 
  char *dt = ctime(&now);
  outFile << "#Target HV settings for run " << run << " from hcal_gain_match.C on " << dt << "#" << endl;
  outFile << "#Row Col targetHV" << endl;

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
      outFile << r << " " << c << " " << targetHV[r][c] << endl; 

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

      /*
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
      */

      /*
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
      */

      if( Flag != "Good" )
	cout << "Problem fit to module " << r << " " << c << " -> " << Flag << endl;

      // Write fitted histograms to file

      cout << "Writing fitted histograms to file.." << endl;
      
      cosHistFile->cd();
      // INT
      gPMTIntSpec[r][c]->SetTitle( Form( "Int Spect R%d C%d FitMean %f", r, c, parsInt[1] ) );
      gPMTIntSpec[r][c]->Write( Form( "Int Spect R%d C%d", r, c ) );
      gPMTIntSpec[r][c]->Draw( "AP" );
      // MAX
      gPMTMaxSpec[r][c]->SetTitle( Form( "Max Spect R%d C%d FitMean %f", r, c, parsMax[1] ) );
      gPMTMaxSpec[r][c]->Write( Form( "Max Spect R%d C%d", r, c ) );
      gPMTMaxSpec[r][c]->Draw( "AP" );
      /*
      // Draw return to baseline histograms
      for( int i = 0; i < kMaxSteps; i++){
	gRetToBase[r][c][i]->SetTitle( Form( "Return to Baseline R%d C%d MaxADC(%d*20)", r, c, i ) );
	gRetToBase[r][c][i]->Write(Form("RTB_R%d_C%d_DSamp(%d*20)", r, c, i ) );
	gRetToBase[r][c][i]->Draw("AP");
      }
      */
    }
  }
    
  // Post analysis reporting
  cout << "Finished loop over run " << run << "." << endl;
  cout << "Total good signals = " << gSignalTotal << "." << endl;
  //cout << "Sample histograms drawn to file HistosFile.root." << endl;
  cout << "Target HV settings written to HVTargets_run" << run << ".txt." << endl;
  if( diagPlots ){
    //cout << "Printing diagnostic plots.." << endl;
    //diagnosticPlots( run, date );

    cout << "No diagnostic plots yet.. Stay tuned!" << endl;
    
  }

  st->Stop();

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

  return 0;
}
  
