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

const int kMinSample = 0;
const int kMaxSample = 165;
const int kTotSample = ( kMaxSample-kMinSample ); // Should be the total fADC window size with each samp = 4ns
const double kTargetRAU = 78.355; // Previous value of 61.425*(2.55/2) or kTargetRAU_prev*(attLongCable_newMeas/attLongCable_oldMeas)
const int kMaxSteps = 10; // Corresponding to sample delta of 20 RAU (1 sigma LED) over 200 RAU (cosmics RAU range)
const double kMaxStep = 20; // Sample delta of 20 RAU (1 sigma LED on either side)
const double PI = TMath::Pi();
const double M_e = 0.00051;
const double M_p = 0.938272;
const double M_n = 0.939565;
const double c_light = 299792458.0;

// Counter to keep track of T tree entry for processing
int gCurrentEntry = -1;
TChain *T = 0;

// Declare necessary histograms
TH1F *gHistos[kNrows][kNcols];

// Create files to keep temporary histograms for integrity checks
TFile *HistosFile = new TFile( "elasWaveformFile.root", "RECREATE" );  // File for checking histogram fits and integrity

// Declare necessary arrays
bool gPulse[kNrows][kNcols]; 
bool gPulseTDC[kNrows][kNcols];

int nevents = 0;
int elasticEvents=0;

TH1F *gMasterHisto[kNrows][kNcols];

// Create generic histogram function
TH1F* MakeHisto(int row, int col, int bins){
  TH1F *h = new TH1F(Form("h%02d%02d",row,col),Form("%d-%d",row+1,col+1),bins,kMinSample,kMaxSample);
  h->SetStats(0);
  h->SetLineWidth(2);
  return h;
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
    }
  }

  // Reset signal peak, adc, and tdc arrays
  double adc[kNrows][kNcols];
  double tdc[kNrows][kNcols];
  for( r = 0; r < kNrows; r++ ) {
    for( c = 0; c < kNcols; c++ ) {
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

    }

    gMasterHisto[r][c]->Add( gHistos[r][c] );
    nevents++;
    // Report error if module is empty
    if(!processed) {
      cerr << "Skipping empty module: " << m << ".." << endl;
      for( int s = kMinSample; s < kTotSample; s++ ) {
        gHistos[r][c]->SetBinContent( s+1 - kMinSample, -404 );
      }
    }
    if( elasticEvents<100 ) {
      HistosFile->cd();
      gHistos[r][c]->Write(Form("ElasticWaveform r%d c%d event%d", r, c, gCurrentEntry));
      elasticEvents++;
    }
    
  }
}
// Main
int elasWaveform(int run = 11311, int event = -1){

  cout << "Enter run number for analysis." << endl;
  cin >> run;
  
  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );
   
  // Read in data produced by analyzer in root format
  cout << "Reading raw data from analyzer.." << endl;
  if( !T ) { 
    T = new TChain( "T" );
    T->Add( Form( "/adaqfs/home/a-onl/sbs/Rootfiles/hcal_gmn_%d*.root", run ) );
    T->SetBranchStatus( "*", 0 );
    T->SetBranchStatus( "sbs.hcal.*", 1 );
    T->SetBranchAddress( "sbs.hcal.a_p", hcalt::a_p );
    T->SetBranchAddress( "sbs.hcal.a_amp", hcalt::a_amp );
    T->SetBranchAddress( "sbs.hcal.a_c", hcalt::a_c );
    T->SetBranchAddress( "sbs.hcal.adcrow", hcalt::row );
    T->SetBranchAddress( "sbs.hcal.adccol", hcalt::col );
    T->SetBranchAddress( "sbs.hcal.a_time", hcalt::a_time );
    T->SetBranchAddress( "bb.tr.chi2", hcalt::BBtr_chi2 );
    T->SetBranchAddress( "bb.tr.n", hcalt::BBtr_n );
    T->SetBranchAddress( "bb.tr.px", hcalt::BBtr_px );
    T->SetBranchAddress( "bb.tr.py", hcalt::BBtr_py );
    T->SetBranchAddress( "bb.tr.pz", hcalt::BBtr_pz );
    T->SetBranchAddress( "bb.tr.p", hcalt::BBtr_p );
    T->SetBranchAddress("sbs.hcal.clus.id",hcalt::cid);
    T->SetBranchAddress("sbs.hcal.clus.row",hcalt::crow);
    T->SetBranchAddress("sbs.hcal.clus.col",hcalt::ccol);
    T->SetBranchAddress("sbs.hcal.clus.atime",hcalt::catime); //a_time for main block in each cluster
    T->SetBranchAddress("sbs.hcal.clus_blk.atime",hcalt::cblkatime); //a_time for each block in main cluster

    T->SetBranchAddress("sbs.hcal.clus_blk.row",hcalt::cbrow);
    T->SetBranchAddress("sbs.hcal.clus_blk.col",hcalt::cbcol);
    T->SetBranchAddress("sbs.hcal.clus_blk.atime",hcalt::cblkatime);
    T->SetBranchAddress("sbs.hcal.clus.e",hcalt::ce);
    T->SetBranchAddress("sbs.hcal.clus.eblk",hcalt::ceblk); // Highest energy block
    T->SetBranchAddress("sbs.hcal.clus.nblk",hcalt::cnblk);
    //T->SetBranchAddress("sbs.hcal.nclus",hcalt::nclus);
    T->SetBranchAddress("sbs.hcal.nblk",hcalt::nblk);
    T->SetBranchAddress("sbs.hcal.clus_blk.id",hcalt::cblkid);
    T->SetBranchAddress("sbs.hcal.clus_blk.e",hcalt::cblke); // Array of block energies
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
      }
    }
  }

  Long64_t Nevents = T->GetEntries();

  for(int r=0; r<kNrows; r++){
    for(int c=0; c<kNcols; c++){
      gMasterHisto[r][c] = new TH1F(Form("Master histo r%d c%d",r,c),Form("%d-%d",r+1,c+1),kTotSample,kMinSample,kMaxSample);
    }
  }

  // Need to get electron energy and W for cuts on elastic events
  long mevent = 0;
  
  // e' momentum correction factor
  double eCorr = 1.04;
  double E_e = 1.92;

  ////GET PHYSICS////
  while( T->GetEntry( mevent++ ) ){ 
    
    if( mevent%1000 == 0 ) cout << "Loop: " << mevent << "/" << Nevents << endl;
    
    // Sort all tracks by lowest Chi^2 (highest confidence track)
    int track_tot = int(hcalt::BBtr_n[0]);
    int track = 0;
    double min = 1000.0; //Set arbitrarily high chi^2 to minimize with loop over tracks

    for( int elem=0; elem<track_tot; elem++ ){            
      if( hcalt::BBtr_chi2[elem] < min )
	track = elem;      
    }
    
    TVector3 lVep( hcalt::BBtr_px[track], hcalt::BBtr_py[track], hcalt::BBtr_pz[track] ); // Construct the scattered electron momentum vector from the tree (GEM reconstructed)
    
    double E_ep = sqrt( pow(M_e,2) + pow(hcalt::BBtr_p[track],2) ); // Obtain the scattered electron energy

    double p_ep = eCorr*hcalt::BBtr_p[track]; // Obtain the magnitude of scattered electron momentum with correction factor

    double theta_e = acos( ( hcalt::BBtr_pz[track])/p_ep ) * 180.0 / PI; // Obtain the electron's scattering angle in degrees
    
    double phi_e = TMath::ASin( hcalt::BBtr_py[track]/( p_ep*TMath::Sin( theta_e ) ) );
    
    double Q2 = 2*E_e*E_ep*( 1-(hcalt::BBtr_pz[track]/p_ep) ); // Obtain Q2 from beam energy, outgoing electron energy, and momenta

    double nu = E_e-E_ep; // Obtain energy transfer

    double W2 = pow( M_p,2 )+2*M_p*nu-Q2; // Obtain W2 from Q2 and nu

    double W = 0;
    if( W2>0 ){
      W = sqrt( W2 ); // Obtain W for real events by throwing out negative W2 solns
    }
    
    // Use the electron kinematics to predict the proton momentum assuming elastic scattering on free proton at rest:

    double E_pp = nu+M_p; // Get energy of the proton

    double p_p = sqrt( pow( E_pp,2 )-W2 ); // Magnitude of the scattered proton

    //double KE_p = pow( p_p,2 )/(2*M_p);
    double KE_p = nu; //For elastics

    // Each component in the form p_p_z = p_p*cos(phi)
    double phi_p = TMath::ACos( ( E_e-hcalt::BBtr_pz[track])/p_p ) * 180.0 / PI;

    // Get the projection of each event onto the face of HCal - consider only X (column) components to avoid calculations with 48D48 field map. Should still be useful

    ////MAKE CUT ON MOMENTUM////
    // Liberal cut on W (W=Mp for elastics)
    if( W<1.1&&W>0.3 ){
      
      gCurrentEntry = event; // Resetting entries for pulse analysis
      
      for ( int i = gCurrentEntry; i < T->GetEntries(); i++ ){ 
	processEvent( gCurrentEntry );
	gCurrentEntry++;
	
	if ( gCurrentEntry%10000 == 0 ){
	  cout << "Current event: " << gCurrentEntry << endl;
	}
      }
    }
  }
  st->Stop();
  
  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;
  
  return 0;
  
}  
