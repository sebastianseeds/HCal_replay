//SSeeds 2.27.23 - Pass2 Update - Calibration code which employs best current cuts on elastic events to obtain ADC time and TDC offset parameters (ns). Applies timewalk corrections. Now loops over all data and cuts on both dx and dy with array of chains.
//NOTE: should be used only after gain parameters are available from ecal.C
//3.13.23 Update - added iteration 2 fuctionality which enables TOF corrections and internal cluster timing to output files, otherwise identical to iteration 1. Note that MC fits for TOF extrapolate well to low momentum, but NOT for high momentum.

#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <TChain.h>
#include <TSystem.h>
#include <TStopwatch.h>
#include <unistd.h>
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TCut.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TLegend.h"
#include "TVector3.h"
#include "TGraphErrors.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TF1.h"

//Detector constants
const Int_t ifac = 3; // Inclusion factor, number of sigma to keep around cut peaks
const Int_t kNcell = 288; // Total number of HCal modules
const Int_t kNrows = 24; // Total number of HCal rows
const Int_t kNcols = 12; // Total number of HCal columns
const Int_t maxTracks = 1000; // Reasonable limit on tracks to be stored per event
const Int_t uni_N = 400; // Total number of bins used to measure detection uniformity (hSampFrac histos)
const Int_t xN = 48; //2*kNrows, total number of dispersive bins detection uni
const Int_t yN = 24; //2*kNcols, total number of transverse bins detection uni
const Double_t hcalheight = -0.2897;
const Double_t TDCCalib = 0.112;
const Double_t TDC_target = -75.; //Target tdc values for aligning block times
const Double_t ADCt_target = 50.; //Target adct values for aligning block times

//Target constants
const Double_t dEdx_tgt=0.00574; //According to NIST ESTAR, the collisional stopping power of hydrogen is about 5.74 MeV*cm2/g at 2 GeV energy
const Double_t dEdx_Al = 0.0021; //According to NIST ESTAR, the collisional stopping power of Aluminum is about 2.1 MeV*cm2/g between 1-4 GeV
const Double_t uwallthick_LH2 = 0.0145; //cm
const Double_t dwallthick_LH2 = 0.015; //cm
const Double_t cellthick_LH2 = 0.02; //cm, this is a guess;
const Double_t Alshieldthick = 2.54/8.0; //= 1/8 inch * 2.54 cm/inch  
const Double_t dxlim_l = 4.0;
const Double_t dxlim_u = 3.0;
const Double_t l_tgt = 0.15; // Length of the target (m)
const Double_t rho_tgt = 0.0723; // Density of target (g/cc)
const Double_t rho_Al = 2.7; // Density of aluminum windows (g/cc)
const Double_t celldiameter = 1.6*2.54; //cm, right now this is a guess
const Double_t Ztgt = 1.0;
const Double_t Atgt = 1.0;
const Double_t Mmol_tgt = 1.008; //g/mol
const UInt_t second = 1000000; 
const Int_t nkine = 6; // Total number of kinematic points
const Int_t maxfset = 5;
const Int_t nfset_lh2[6] = {3,1,2,2,4,1}; //Number of LH2 magnetic field settings in order SBS{4,7,11,14,8,9}
const Int_t nfset_ld2[6] = {3,1,2,1,4,1}; //Number of LD2 magnetic field settings in order SBS{4,7,11,14,8,9}
const Int_t fset_lh2[6][4] = {{0,30,50,-1},
			      {85,-1,-1,-1},
			      {0,100,-1,-1},
			      {0,70,-1,-1},
			      {0,50,70,100},
			      {70,-1,-1,-1}}; //All field settings for all kinematics LH2. -1 indicates no setting
const Int_t fset_ld2[6][4] = {{0,30,50,-1},
			      {85,-1,-1,-1},
			      {0,100,-1,-1},
			      {70,-1,-1,-1},
			      {0,50,70,100},
			      {70,-1,-1,-1}}; //All field settings for all kinematics LH2. -1 indicates no setting
const Double_t TOFfitp_p[6][4] = {{73.7816, -35.8806, 12.90060, -1.60996},
				  {52.7549, -2.05791, 0.273797, -0.01309},
				  {54.6380, -1.97939, 0.231943, -0.00934},
				  {65.2843, -10.6005, 2.252520, -0.16723},
				  {60.6177, -18.2422, 5.123620, -0.49684},
				  {71.7263, -32.1323, 41.57460, -1.18023}}; //Pol3 TOF fit params from proton MC data with SBS fields {30,85,100,70,70,70}
const Double_t TOFfitp_n[6][4] = {{62.1158, -20.8086, 6.464510, -0.71342},
				  {55.8218, -3.84358, 06176470, -0.03516},
				  {51.4030, -0.51868, 0.016076, 0.001009},
				  {64.1014, -9.30317, 1.818670, -0.12324},
				  {56.1597, -13.4154, 3.449660, -0.30995},
				  {63.9328, -22.3211, 6.673740, -0.68472}}; //Pol3 TOF fit params from neutron MC data with SBS fields {30,85,100,70,70,70} 
const Double_t pulim[6] = {2.8,6.7,8.9,5.4,3.9,3.7}; //p upper limit on fits, beyond this fits diverge
const Int_t tofB[6] = {30,85,100,70,70,70};

//Math/other
const Double_t PI = TMath::Pi();
const Double_t M_e = 0.00051;
const Double_t M_p = 0.938272;
const Int_t tFitMin = 120;

// Get today's date
string getDate(){
  time_t now = time(0);
  tm ltm = *localtime(&now);
  
  string yyyy = to_string(1900 + ltm.tm_year);
  string mm = to_string(1 + ltm.tm_mon);
  string dd = to_string(ltm.tm_mday);
  string date = mm + '_' + dd + '_' + yyyy;
  
  return date;
}

// Dying exponential fit for timewalk corrections
Double_t TW_fit(Double_t *x, Double_t *par){
  Double_t amp = par[0];
  Double_t str = par[1];
  Double_t offset = par[2];
  return amp*exp(-str*x[0])+offset;
}

//Keep this around just in case it's needed later
Double_t Gfit(Double_t *x, Double_t *par){
  Double_t amp = par[0];
  Double_t offset = par[1];
  Double_t sigma = par[2];
  return amp*exp(-0.5*pow((x[0]-offset)/sigma,2.));
}

//Main. Calibration one kinematic at a time. Loads all configs for both lh2 and ld2.
void tcal( Int_t kine=-1, Int_t iter=0 ){

  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );
  
  // Get the date
  string date = getDate();

  // Declare Chain for many root files
  // First obtain index for nfset constants
  Int_t kIdx;
  if( kine == 4 ) kIdx=0;
  if( kine == 7 ) kIdx=1;
  if( kine == 11 ) kIdx=2;
  if( kine == 14 ) kIdx=3;
  if( kine == 8 ) kIdx=4;
  if( kine == 9 ) kIdx=5;

  //Declare bool indicating iteration
  bool qreplay = iter>0;

  bool tofreplay = iter==2;

  //Set up arrays to hold gain parameters for energy corrections
  Double_t gOldConst_pass0[kNcell] = {0.};
  Double_t gConst_iter1[kNcell] = {0.}; 

  // Declare arrays to hold old offset parameters
  Double_t oldADCtoffsets[kNcell]={0.};
  Double_t oldTDCoffsets[kNcell]={0.};
  Double_t otdcP0[kNcell]={0.};
  Double_t otdcP1[kNcell]={0.};
  Double_t otdcP2[kNcell]={0.};
  Double_t oadctP0[kNcell]={0.};
  Double_t oadctP1[kNcell]={0.};
  Double_t oadctP2[kNcell]={0.};
  Double_t calADCtoffsets[kNcell]={0.};
  Double_t calTDCoffsets[kNcell]={0.};
  bool pass0 = false;

  //Get gain parameters from ecal.C for use with timewalk corrections
  string inConstPath = Form("coefficients/coeff_sbs%d.txt",kine);

  cout << endl << "Loading gain coefficents from new calibration: " << inConstPath << ".." << endl;
  ifstream inConstFile( inConstPath );
  if( !inConstFile ){
    cerr << endl << "ERROR: No input constant file present -> path to coeff_sbs<kine>.txt expected." << endl;
    return 0;
  }

  Int_t n1=0;
  Double_t d1=0;
  string line;
  bool skip_line = true;

  while( getline( inConstFile, line ) ){
      
    if( n1==( kNcell ) ) break;
      
    TString Tline = (TString)line;
      
    if( Tline.BeginsWith("sbs.hcal.adc.gain") && skip_line==true ){	
      skip_line = false;
      continue;
    }
      
    if( skip_line==false ){
      istringstream iss( line );
      while( iss >> d1 ){
	gConst_iter1[n1] = d1;
	n1++;
      }
    }
    
  }

  cout << endl << "ADC gain params from ecal.C : " << endl;
    
  for( Int_t r=0; r<kNrows; r++){
    for( Int_t c=0; c<kNcols; c++){
      Int_t i = r*kNcols+c;
      cout << gConst_iter1[i] << " ";
    }
    cout << endl;
  }

  usleep( 3*second ); //Give some time for review of step 
  
  cout << endl;
  
  if( kine==4 || kine==7 ){ //end if iteration 1 (quasi-replay)
    
    pass0=true;
    
    // Reading ADC gain parameters from database 
    string inOldConstPath = "/u/home/sbs-gmn/pass0/SBS_REPLAY/SBS-replay/DB/db_sbs.hcal.dat";
    
    cout << "Loading previous gain coefficients from file for pass0 kinematics: " << inOldConstPath << ".." << endl;
    ifstream inOldConstFile( inOldConstPath );
    if( !inOldConstFile ){
      cerr << endl << "ERROR: No input constant file present -> path to db_sbs.hcal.dat expected." << endl;
      return 0;
    }
    
    Int_t n1=0;
    Double_t d1=0;
    string line;
    bool skip_line = true;
    bool skip_one_line = true;
    bool skip_first_instance = true;
    
    while( getline( inOldConstFile, line ) ){
      
      if( n1==( kNcell ) ) break;
      
      TString Tline = (TString)line;
      
      if( Tline.BeginsWith("sbs.hcal.adc.gain") && skip_first_instance==true ){
	skip_first_instance = false;
	continue;
      }
      
      if( Tline.BeginsWith("sbs.hcal.adc.gain") && skip_line==true && skip_first_instance==false ) skip_line = false;
      
      if( skip_line==false && skip_one_line==true ){
	skip_one_line = false;
	continue;
      }
      
      if( skip_line==false && skip_one_line==false ){
	istringstream iss( line );
	while( iss >> d1 ){
	  gOldConst_pass0[n1] = d1;
	  n1++;
	}
      }
    }
    
    cout << endl << "Old ADC gain params from pass0 db file: " << endl;
    
    for( Int_t r=0; r<kNrows; r++){
      for( Int_t c=0; c<kNcols; c++){
	Int_t i = r*kNcols+c;
	cout << gOldConst_pass0[i] << " ";
      }
      cout << endl;
    }
    
    cout << endl;
  }//end if pass0

  // Paths to read constants
  string RtdctwP0path = Form("parameters/tdctwP0_sbs%d.txt",kine);
  string RtdctwP1path = Form("parameters/tdctwP1_sbs%d.txt",kine);
  string RtdctwP2path = Form("parameters/tdctwP2_sbs%d.txt",kine);
  string RadcttwP0path = Form("parameters/adcttwP0_sbs%d.txt",kine);
  string RadcttwP1path = Form("parameters/adcttwP1_sbs%d.txt",kine);
  string RadcttwP2path = Form("parameters/adcttwP2_sbs%d.txt",kine);
  string tdcoffpath = Form("parameters/tdcoffsets_sbs%d.txt",kine);
  string adctoffpath = Form("parameters/adctoffsets_sbs%d.txt",kine);

  if( qreplay ){

    // Reading fit parameters if quasi-replay iteration (1)
    cout << "Loading fit parameters.." << endl;

    //Get first parameter (all channels)
    ifstream tdctwP0file( RtdctwP0path );
    if( !tdctwP0file ){
      cerr << endl << "ERROR: No input constant file present -> path to tdctwP0.txt expected." << endl;
      return 0;
    }
  
    Int_t n1=0;
    Double_t d1;
    string line1;
    
    while( getline( tdctwP0file, line1 ) ){
      if( line1.at(0) == '#' ) {
	continue;
      }
      stringstream ss( line1 );
      ss >> d1;     
      otdcP0[n1] = d1;      
      n1++;
    }

    //Get second parameter (all channels)
    ifstream tdctwP1file( RtdctwP1path );
    if( !tdctwP1file ){
      cerr << endl << "ERROR: No input constant file present -> path to tdctwP1.txt expected." << endl;
      return 0;
    }

    n1=0;
    d1=0;
    string line2;
    
    while( getline( tdctwP1file, line2 ) ){
      if( line2.at(0) == '#' ) {
	continue;
      }
      stringstream ss( line2 );
      ss >> d1;    
      otdcP1[n1] = d1; 
      n1++;
    }

    //Get third parameter (all channels)
    ifstream tdctwP2file( RtdctwP2path );
    if( !tdctwP2file ){
      cerr << endl << "ERROR: No input constant file present -> path to tdctwP2.txt expected." << endl;
      return 0;
    }

    n1=0;
    d1=0;
    string line3;
    
    while( getline( tdctwP2file, line3 ) ){
      if( line3.at(0) == '#' ) {
	continue;
      }
      stringstream ss( line3 );
      ss >> d1;     
      otdcP2[n1] = d1;     
      n1++;
    }

    //Get first parameter (all channels)
    ifstream adcttwP0file( RadcttwP0path );
    if( !adcttwP0file ){
      cerr << endl << "ERROR: No input constant file present -> path to adcttwP0.txt expected." << endl;
      return 0;
    }
  
    n1=0;
    d1;
    string line4;
    
    while( getline( adcttwP0file, line4 ) ){
      if( line4.at(0) == '#' ) {
	continue;
      }
      stringstream ss( line4 );
      ss >> d1;     
      oadctP0[n1] = d1;      
      n1++;
    }

    //Get second parameter (all channels)
    ifstream adcttwP1file( RadcttwP1path );
    if( !adcttwP1file ){
      cerr << endl << "ERROR: No input constant file present -> path to adcttwP1.txt expected." << endl;
      return 0;
    }

    n1=0;
    d1=0;
    string line5;
    
    while( getline( adcttwP1file, line5 ) ){
      if( line5.at(0) == '#' ) {
	continue;
      }
      stringstream ss( line5 );
      ss >> d1;    
      oadctP1[n1] = d1; 
      n1++;
    }

    //Get third parameter (all channels)
    ifstream adcttwP2file( RadcttwP2path );
    if( !adcttwP2file ){
      cerr << endl << "ERROR: No input constant file present -> path to adcttwP2.txt expected." << endl;
      return 0;
    }

    n1=0;
    d1=0;
    string line6;
    
    while( getline( adcttwP2file, line6 ) ){
      if( line6.at(0) == '#' ) {
	continue;
      }
      stringstream ss( line6 );
      ss >> d1;     
      oadctP2[n1] = d1;     
      n1++;
    }

    //Get new TDC offsets (all channels)
    ifstream newtdcfile( tdcoffpath );
    if( !newtdcfile ){
      cerr << endl << "ERROR: No input constant file present -> path to tdcoffsets.txt expected." << endl;
      return 0;
    }

    n1=0;
    d1=0;
    string line7;
    
    while( getline( newtdcfile, line7 ) ){
      if( line7.at(0) == '#' ) {
	continue;
      }
      stringstream ss( line7 );
      ss >> d1;     
      calTDCoffsets[n1] = d1;     
      n1++;
    }

    //Get new ADCt offsets (all channels)
    ifstream newadctfile( adctoffpath );
    if( !newadctfile ){
      cerr << endl << "ERROR: No input constant file present -> path to adctoffsets.txt expected." << endl;
      return 0;
    }

    n1=0;
    d1=0;
    string line8;
    
    while( getline( newadctfile, line8 ) ){
      if( line8.at(0) == '#' ) {
	continue;
      }
      stringstream ss( line8 );
      ss >> d1;     
      calADCtoffsets[n1] = d1;     
      n1++;
    }


    //Print loaded fit parameters to console
    cout << endl << endl << "TDC Timewalk Fit P0 vs Channels: " << endl;
    
    for( Int_t r=0; r<kNrows; r++){
      for( Int_t c=0; c<kNcols; c++){
	Int_t i = r*kNcols+c;
	cout << otdcP0[i] << " ";
      }
      cout << endl;
    }
    
    usleep( 2*second ); //Give some time for review of step

    cout << endl << endl << "TDC Timewalk Fit P1 vs Channels: " << endl;
    
    for( Int_t r=0; r<kNrows; r++){
      for( Int_t c=0; c<kNcols; c++){
	Int_t i = r*kNcols+c;
	cout << otdcP1[i] << " ";
      }
      cout << endl;
    }
    
    usleep( 2*second ); //Give some time for review of step

    cout << endl << endl << "TDC Timewalk Fit P2 vs Channels: " << endl;
    
    for( Int_t r=0; r<kNrows; r++){
      for( Int_t c=0; c<kNcols; c++){
	Int_t i = r*kNcols+c;
	cout << otdcP2[i] << " ";
      }
      cout << endl;
    }

    //Print loaded fit parameters to console
    cout << endl << endl << "ADCt Timewalk Fit P0 vs Channels: " << endl;
    
    for( Int_t r=0; r<kNrows; r++){
      for( Int_t c=0; c<kNcols; c++){
	Int_t i = r*kNcols+c;
	cout << oadctP0[i] << " ";
      }
      cout << endl;
    }
    
    usleep( 2*second ); //Give some time for review of step

    cout << endl << endl << "ADCt Timewalk Fit P1 vs Channels: " << endl;
    
    for( Int_t r=0; r<kNrows; r++){
      for( Int_t c=0; c<kNcols; c++){
	Int_t i = r*kNcols+c;
	cout << oadctP1[i] << " ";
      }
      cout << endl;
    }
    
    usleep( 2*second ); //Give some time for review of step

    cout << endl << endl << "ADCt Timewalk Fit P2 vs Channels: " << endl;
    
    for( Int_t r=0; r<kNrows; r++){
      for( Int_t c=0; c<kNcols; c++){
	Int_t i = r*kNcols+c;
	cout << oadctP2[i] << " ";
      }
      cout << endl;
    }
    
    usleep( 2*second ); //Give some time for review of step

    cout << endl << endl << "TDC offsets vs Channels: " << endl;
    
    for( Int_t r=0; r<kNrows; r++){
      for( Int_t c=0; c<kNcols; c++){
	Int_t i = r*kNcols+c;
	cout << calTDCoffsets[i] << " ";
      }
      cout << endl;
    }
   
    usleep( 2*second ); //Give some time for review of step

    cout << endl << endl << "ADCt offsets vs Channels: " << endl;
    
    for( Int_t r=0; r<kNrows; r++){
      for( Int_t c=0; c<kNcols; c++){
	Int_t i = r*kNcols+c;
	cout << calADCtoffsets[i] << " ";
      }
      cout << endl;
    }
  }

  // Path for previous ADCt and TDC constants
  string inOldConstPath = "/u/home/sbs-gmn/pass0/SBS_REPLAY/SBS-replay/DB/db_sbs.hcal.dat";

  // Reading ADC and TDC timing offsets from database
  cout << endl << "Loading previous offsets from database file: " << inOldConstPath << ".." << endl;
  ifstream inOldConstFile2( inOldConstPath );
  if( !inOldConstFile2 ){
    cerr << endl << "ERROR: No input constant file present -> path to db_sbs.hcal.dat expected." << endl;
    return 0;
  }

  Int_t n0=0;
  Double_t d0=0;
  string line0;
  skip_line = true;
  bool skip_one_line = true;
  bool pass_first_cond = false;
  bool pass_all_cond = false;

  while( getline( inOldConstFile2, line0 ) ){

    if( pass_first_cond && n0==( kNcell ) ) break;

    TString Tline = (TString)line0;

    if( Tline.BeginsWith("sbs.hcal.adc.timeoffset") && skip_line==true ) skip_line = false;
    if( Tline.BeginsWith("sbs.hcal.tdc.offset") && skip_line==true ) skip_line = false;

    if( skip_line==false && skip_one_line==true ){
      skip_one_line = false;
      continue;
    }

    if( n0==( kNcell ) && pass_first_cond==false ){
      skip_line = true;
      skip_one_line = true;
      n0=0;
      pass_first_cond = true;
    }

    if( skip_line==false && skip_one_line==false ){
      istringstream iss( line0 );
      while( iss >> d0 ){
	if( pass_first_cond==false ){
	  oldADCtoffsets[n0] = d0;
	}else{
	  oldTDCoffsets[n0] = d0;
	}
	n0++;
      }
    }
  }

  cout << endl << endl << "Old TDC offsets: " << endl;

  for( Int_t r=0; r<kNrows; r++){
    for( Int_t c=0; c<kNcols; c++){
      Int_t i = r*kNcols+c;
      cout << oldTDCoffsets[i] << " ";
    }
    cout << endl;
  }

  cout << endl << endl << "Old ADC time offsets: " << endl;

  for( Int_t r=0; r<kNrows; r++){
    for( Int_t c=0; c<kNcols; c++){
      Int_t i = r*kNcols+c;
      cout << oldADCtoffsets[i] << " ";
    }
    cout << endl;
  }

  cout << endl << endl << "Database parameters loaded." << endl;

  usleep( 3*second ); //Give some time for review of step

  //Set up multiple chains and elists for hydrogen and deuterium
  TChain *Ch[nfset_lh2[kIdx]];
  TChain *Cd[nfset_ld2[kIdx]];

  TEventList *elist_h[nfset_lh2[kIdx]];
  Long64_t Nevents_h[nfset_lh2[kIdx]];
  Long64_t NTevents_h[nfset_lh2[kIdx]];

  TEventList *elist_d[nfset_ld2[kIdx]];
  Long64_t Nevents_d[nfset_ld2[kIdx]];
  Long64_t NTevents_d[nfset_ld2[kIdx]];

  TCut globalcut_h[nfset_lh2[kIdx]];
  TCut globalcut_d[nfset_ld2[kIdx]];

  // Declare general detector parameters - lh2
  Double_t BBtr_p_h[nfset_lh2[kIdx]][maxTracks], BBtr_px_h[nfset_lh2[kIdx]][maxTracks], BBtr_py_h[nfset_lh2[kIdx]][maxTracks], BBtr_pz_h[nfset_lh2[kIdx]][maxTracks];
  Double_t BBtr_vz_h[nfset_lh2[kIdx]][maxTracks], BBtr_chi2_h[nfset_lh2[kIdx]][maxTracks];
  Double_t BBtr_n_h[nfset_lh2[kIdx]], BBps_x_h[nfset_lh2[kIdx]], BBps_y_h[nfset_lh2[kIdx]], BBps_e_h[nfset_lh2[kIdx]], BBsh_x_h[nfset_lh2[kIdx]], BBsh_y_h[nfset_lh2[kIdx]], BBsh_e_h[nfset_lh2[kIdx]];

  Double_t HCALx_h[nfset_lh2[kIdx]], HCALy_h[nfset_lh2[kIdx]], HCALe_h[nfset_lh2[kIdx]];
  Double_t crow_h[nfset_lh2[kIdx]], ccol_h[nfset_lh2[kIdx]], nblk_h[nfset_lh2[kIdx]];
  Double_t cblkid_h[nfset_lh2[kIdx]][kNcell], cblke_h[nfset_lh2[kIdx]][kNcell], cblkatime_h[nfset_lh2[kIdx]][kNcell], cblktime_h[nfset_lh2[kIdx]][kNcell];
  Double_t cblkagain_h[nfset_lh2[kIdx]][kNcell];
  Double_t ekineW2_h[nfset_lh2[kIdx]];
  Double_t HODOtmean_h[nfset_lh2[kIdx]];

  // Declare general detector parameters - ld2
  Double_t BBtr_p_d[nfset_ld2[kIdx]][maxTracks], BBtr_px_d[nfset_ld2[kIdx]][maxTracks], BBtr_py_d[nfset_ld2[kIdx]][maxTracks], BBtr_pz_d[nfset_ld2[kIdx]][maxTracks];
  Double_t BBtr_vz_d[nfset_ld2[kIdx]][maxTracks], BBtr_chi2_d[nfset_ld2[kIdx]][maxTracks];
  Double_t BBtr_n_d[nfset_ld2[kIdx]], BBps_x_d[nfset_ld2[kIdx]], BBps_y_d[nfset_ld2[kIdx]], BBps_e_d[nfset_ld2[kIdx]], BBsh_x_d[nfset_ld2[kIdx]], BBsh_y_d[nfset_ld2[kIdx]], BBsh_e_d[nfset_ld2[kIdx]];

  Double_t HCALx_d[nfset_ld2[kIdx]], HCALy_d[nfset_ld2[kIdx]], HCALe_d[nfset_ld2[kIdx]];
  Double_t crow_d[nfset_ld2[kIdx]], ccol_d[nfset_ld2[kIdx]], nblk_d[nfset_ld2[kIdx]];
  Double_t cblkid_d[nfset_ld2[kIdx]][kNcell], cblke_d[nfset_ld2[kIdx]][kNcell], cblkatime_d[nfset_ld2[kIdx]][kNcell], cblktime_d[nfset_ld2[kIdx]][kNcell];
  Double_t cblkagain_d[nfset_ld2[kIdx]][kNcell];
  Double_t ekineW2_d[nfset_ld2[kIdx]];
  Double_t HODOtmean_d[nfset_ld2[kIdx]];

  // Declare file name paths for writing constants
  string tdcconstPath;
  string adctconstPath;
  string tdctwP0path;
  string tdctwP1path;
  string tdctwP2path;
  string adcttwP0path;
  string adcttwP1path;
  string adcttwP2path;
  string errPath;
  if( qreplay ){
    tdcconstPath = Form("parameters/tdcoffsets_sbs%d_qreplay.txt",kine);
    adctconstPath = Form("parameters/adctoffsets_sbs%d_qreplay.txt",kine);
    tdctwP0path = Form("parameters/tdctwP0_sbs%d_qreplay.txt",kine);
    tdctwP1path = Form("parameters/tdctwP1_sbs%d_qreplay.txt",kine);
    tdctwP2path = Form("parameters/tdctwP2_sbs%d_qreplay.txt",kine);
    adcttwP0path = Form("parameters/adcttwP0_sbs%d_qreplay.txt",kine);
    adcttwP1path = Form("parameters/adcttwP1_sbs%d_qreplay.txt",kine);
    adcttwP2path = Form("parameters/adcttwP2_sbs%d_qreplay.txt",kine);
    errPath = Form("reporting/term/tcalReport_sbs%d_qreplay.txt",kine);
  }else{
    tdcconstPath = Form("parameters/tdcoffsets_sbs%d.txt",kine);
    adctconstPath = Form("parameters/adctoffsets_sbs%d.txt",kine);
    tdctwP0path = Form("parameters/tdctwP0_sbs%d.txt",kine);
    tdctwP1path = Form("parameters/tdctwP1_sbs%d.txt",kine);
    tdctwP2path = Form("parameters/tdctwP2_sbs%d.txt",kine);
    adcttwP0path = Form("parameters/adcttwP0_sbs%d.txt",kine);
    adcttwP1path = Form("parameters/adcttwP1_sbs%d.txt",kine);
    adcttwP2path = Form("parameters/adcttwP2_sbs%d.txt",kine);
    errPath = Form("reporting/term/tcalReport_sbs%d.txt",kine);
  }

  //Declare report outfile
  ofstream report;

  report.open( errPath );
  report << "#Error and performance report from SBS-" << kine << " obtained " << date.c_str() << endl << endl;

  string configfilename_h[nfset_lh2[kIdx]];
  for( Int_t h=0; h<nfset_lh2[kIdx]; h++){
    Int_t field = fset_lh2[kIdx][h];
    if(field==-1) cout << "Error. configfilename array out of bounds." << endl;
    configfilename_h[h] = Form("../config/GMn/SBS%d/secal_lh2_sbs%d_f%d.cfg",kine,kine,field);
  }
  string configfilename_d[nfset_ld2[kIdx]];
  for( Int_t d=0; d<nfset_ld2[kIdx]; d++){
    Int_t field = fset_ld2[kIdx][d];
    if(field==-1) cout << "Error. configfilename array out of bounds." << endl;
    configfilename_d[d] = Form("../config/GMn/SBS%d/secal_ld2_sbs%d_f%d.cfg",kine,kine,field);
  }

  // Declare general physics/fit parameters
  Int_t minEventPerCell = 30; // Minimum number of scattered p in cell required to calibrate
  Int_t maxEventPerCell = 4000; // Maximum number of scattered p events to contribute
  Double_t HCal_divx = 0.15875; // Transverse width in x and y per cell
  Double_t HCal_divy = 0.15494;
  Double_t HCal_Xi = -2.355005; // Distance from beam center to top of HCal in m from MC database
  Double_t HCal_Xf = 1.454995; // Distance from beam center to bottom of HCal in m from MC database
  Double_t HCal_Yi = -0.92964; // Distance from beam center to opposite-beam side of HCal in m from MC database
  Double_t HCal_Yf = 0.92964;
  Double_t aTDCsig = 3.0; //approximate TDC sigma until proven otherwise (ns)
  Double_t aADCtsig = 3.0; //approximate TDC sigma until proven otherwise (ns)
  //Double_t highDelta = 0.1; // Minimum M(i,j)/b(i) factor allowed 

  // Declare general physics parameters to be modified by input config file - LH2
  //Double_t E_e_h[nfset_lh2[kIdx]] = {0.}; // Energy of beam (incoming electrons from accelerator)  
  Double_t E_e_h[maxfset] = {-1.}; // Energy of beam (incoming electrons from accelerator)
  Double_t HCal_d_h[maxfset] = {-1.}; // Distance to HCal from scattering chamber for comm1
  Double_t HCal_th_h[maxfset] = {-1.}; // Angle that the center of HCal is at
  Double_t BB_th_h[maxfset] = {-1.}; // Angle that the center of BBCal is at
  Double_t W2_mean_h[maxfset] = {-1.}; // Mean of W at current kinematic
  Double_t W2_sig_h[maxfset] = {-1.}; // Width of W at current kinematic
  Double_t dx0_n_h[maxfset] = {-1.}; // Position of neutron spot, x-x_expected
  Double_t dx0_p_h[maxfset] = {-1.}; // Position of proton spot, x-x_expected
  Double_t dy0_h[maxfset] = {-1.}; // Position of hadron spot, y-y_expected
  Double_t dx_sig_n_h[maxfset] = {-1.}; // Max spread of neutron spot, x-x_expected
  Double_t dx_sig_p_h[maxfset] = {-1.}; // Max spread of proton spot, x-x_expected
  Double_t dy_sig_h[maxfset] = {-1.}; // Max spread of hadron spot, y-y_expected
  Double_t atime0_h[maxfset] = {-1.}; // Expected location in ADC time of signal
  Double_t atime_sig_h[maxfset] = {-1.}; // 1 sig of atime distribution
  Int_t useAlshield_h[maxfset] = {0};

  // Declare general physics parameters to be modified by input config file - LD2
  //Double_t E_e_d[nfset_ld2[kIdx]] = {0.}; // Energy of beam (incoming electrons from accelerator)
  Double_t E_e_d[maxfset] = {-1.}; // Energy of beam (incoming electrons from accelerator)
  Double_t HCal_d_d[maxfset] = {-1.}; // Distance to HCal from scattering chamber for comm1
  Double_t HCal_th_d[maxfset] = {-1.}; // Angle that the center of HCal is at
  Double_t BB_th_d[maxfset] = {-1.}; // Angle that the center of BBCal is at
  Double_t W2_mean_d[maxfset] = {-1.}; // Mean of W at current kinematic
  Double_t W2_sig_d[maxfset] = {-1.}; // Width of W at current kinematic
  Double_t dx0_n_d[maxfset] = {-1.}; // Position of neutron spot, x-x_expected
  Double_t dx0_p_d[maxfset] = {-1.}; // Position of proton spot, x-x_expected
  Double_t dy0_d[maxfset] = {-1.}; // Position of hadron spot, y-y_expected
  Double_t dx_sig_n_d[maxfset] = {-1.}; // Max spread of neutron spot, x-x_expected
  Double_t dx_sig_p_d[maxfset] = {-1.}; // Max spread of proton spot, x-x_expected
  Double_t dy_sig_d[maxfset] = {-1.}; // Max spread of hadron spot, y-y_expected
  Double_t atime0_d[maxfset] = {-1.}; // Expected location in ADC time of signal
  Double_t atime_sig_d[maxfset] = {-1.}; // 1 sig of atime distribution
  Int_t useAlshield_d[maxfset] = {0};

  Int_t elasYield = 0; // Keep track of total elastics analyzed
  Int_t badtimeblkclus = 0; // Keep track of multi-blk clusters with out of time blks
  Int_t multblkclus = 0; // Keep track of total multi-block clusters (nblk>1)

  //For position reconstruction
  Double_t HCal_Xmin = HCal_Xi-HCal_divx/2;
  Double_t HCal_Xmax = HCal_Xf+HCal_divx/2;
  Double_t HCal_Ymin = HCal_Yi-HCal_divy/2;
  Double_t HCal_Ymax = HCal_Yf+HCal_divy/2;

  // Reading all relevant config files for lh2 and setting up chain
  for( Int_t f=0; f<nfset_lh2[kIdx]; f++ ){
    Ch[f] = new TChain("T");

    ifstream configfile(configfilename_h[f]);
    TString currentline;
    while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endlist") ){

      if( !currentline.BeginsWith("#") ){
	Ch[f]->Add(currentline);
	cout << "Loading file: " << currentline << ".." << endl;
      }    
    }
    //TCut globalcut = "";
    while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endcut") ){
      if( !currentline.BeginsWith("#") ){
	globalcut_h[f] += currentline;
      }    
    }
    while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("#") ){
      TObjArray *tokens = currentline.Tokenize(" ");
      Int_t ntokens = tokens->GetEntries();
      if( ntokens>1 ){
	TString skey = ( (TObjString*)(*tokens)[0] )->GetString();
	if( skey == "E_e" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  E_e_h[f] = sval.Atof();
	  cout << "Loading beam energy: " << E_e_h[f] << endl;
	}
	if( skey == "HCal_d" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  HCal_d_h[f] = sval.Atof();
	  cout << "Loading HCal distance: " << HCal_d_h[f] << endl;
	}
	if( skey == "HCal_th" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  HCal_th_h[f] = sval.Atof() * TMath::DegToRad();	
	  cout << "Loading HCal angle: " << HCal_th_h[f] << endl;
	}
	if( skey == "BB_th" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  BB_th_h[f] = sval.Atof() * TMath::DegToRad();	
	  cout << "Loading BBCal angle: " << BB_th_h[f] << endl;
	}
	if( skey == "W2_mean" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  W2_mean_h[f] = sval.Atof();
	  cout << "Loading W2 mean cut: " << W2_mean_h[f] << endl;
	}
	if( skey == "W2_sig" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  W2_sig_h[f] = sval.Atof();
	  cout << "Loading W2 sigma cut: " << W2_sig_h[f] << endl;
	}
	if( skey == "dx0_n" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  dx0_n_h[f] = sval.Atof();
	  cout << "Loading x position of neutron spot: " << dx0_n_h[f] << endl;
	}
	if( skey == "dx0_p" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  dx0_p_h[f] = sval.Atof();
	  cout << "Loading y position of proton spot: " << dx0_p_h[f] << endl;
	}
	if( skey == "dy0" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  dy0_h[f] = sval.Atof();
	  cout << "Loading y position of both hadron spots: " << dy0_h[f] << endl;
	}
	if( skey == "dx_sig_n" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  dx_sig_n_h[f] = sval.Atof();
	  cout << "Loading x sigma of neutron spot: " << dx_sig_n_h[f] << endl;
	}
	if( skey == "dx_sig_p" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  dx_sig_p_h[f] = sval.Atof();
	  cout << "Loading x sigma of proton spot: " << dx_sig_p_h[f] << endl;
	}
	if( skey == "dy_sig" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  dy_sig_h[f] = sval.Atof();
	  cout << "Loading y sigma of both hadron spots: " << dy_sig_h[f] << endl;
	}
	if( skey == "atime0" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  atime0_h[f] = sval.Atof();
	  cout << "Loading ADC time mean: " << atime0_h[f] << endl;
	}
	if( skey == "atime_sig" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  atime_sig_h[f] = sval.Atof();
	  cout << "Loading ADC time sigma: " << atime_sig_h[f] << endl;
	}
	if( skey == "useAlshield" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  useAlshield_h[f] = sval.Atoi();
	  cout << "Loading Aluminum absorber option: " << useAlshield_h[f] << endl;
	}
      }
      delete tokens;
    }

    NTevents_h[f] = Ch[f]->GetEntries();
    elist_h[f] = new TEventList(Form("elist_lh2_sbs%d_f%d",kine,fset_lh2[kIdx][f]),Form("Elastic Event List, LH2, SBS%d, Mag%d",kine,fset_lh2[kIdx][f]));
    Ch[f]->Draw(Form(">>elist_lh2_sbs%d_f%d",kine,fset_lh2[kIdx][f]),globalcut_h[f]);
    Nevents_h[f] = elist_h[f]->GetN();

    Ch[f]->SetBranchStatus( "*", 0 );
    Ch[f]->SetBranchStatus( "sbs.hcal.x", 1 );
    Ch[f]->SetBranchStatus( "sbs.hcal.y", 1 );
    Ch[f]->SetBranchStatus( "sbs.hcal.e", 1 );
    Ch[f]->SetBranchStatus( "sbs.hcal.rowblk", 1 );
    Ch[f]->SetBranchStatus( "sbs.hcal.colblk", 1 );
    Ch[f]->SetBranchStatus( "sbs.hcal.nblk", 1 );
    Ch[f]->SetBranchStatus( "sbs.hcal.clus_blk.id", 1 );
    Ch[f]->SetBranchStatus( "sbs.hcal.clus_blk.e", 1 );
    Ch[f]->SetBranchStatus( "sbs.hcal.clus_blk.tdctime", 1 );
    Ch[f]->SetBranchStatus( "sbs.hcal.clus_blk.atime", 1 );
    Ch[f]->SetBranchStatus( "sbs.hcal.clus_blk.again", 1 );
    Ch[f]->SetBranchStatus( "bb.tr.chi2", 1 );
    Ch[f]->SetBranchStatus( "bb.tr.n", 1 );
    Ch[f]->SetBranchStatus( "bb.tr.px", 1 );
    Ch[f]->SetBranchStatus( "bb.tr.py", 1 );
    Ch[f]->SetBranchStatus( "bb.tr.pz", 1 );    
    Ch[f]->SetBranchStatus( "bb.tr.vz", 1 );
    Ch[f]->SetBranchStatus( "bb.tr.p", 1 );
    Ch[f]->SetBranchStatus( "bb.ps.e", 1 );
    Ch[f]->SetBranchStatus( "bb.ps.x", 1 );
    Ch[f]->SetBranchStatus( "bb.ps.y", 1 );
    Ch[f]->SetBranchStatus( "bb.sh.e", 1 );
    Ch[f]->SetBranchStatus( "bb.sh.x", 1 );
    Ch[f]->SetBranchStatus( "bb.sh.y", 1 );
    Ch[f]->SetBranchStatus( "bb.hodotdc.clus.tmean", 1 );
    Ch[f]->SetBranchStatus( "e.kine.W2", 1 );

    Ch[f]->SetBranchAddress( "sbs.hcal.x", &HCALx_h[f] );
    Ch[f]->SetBranchAddress( "sbs.hcal.y", &HCALy_h[f] );
    Ch[f]->SetBranchAddress( "sbs.hcal.e", &HCALe_h[f] );
    Ch[f]->SetBranchAddress( "sbs.hcal.rowblk", &crow_h[f] );
    Ch[f]->SetBranchAddress( "sbs.hcal.colblk", &ccol_h[f] );
    Ch[f]->SetBranchAddress( "sbs.hcal.nblk", &nblk_h[f] ); // Total number of blocks in highest E clus
    Ch[f]->SetBranchAddress( "sbs.hcal.clus_blk.id", cblkid_h[f] ); // kNcell-1 index for each block
    Ch[f]->SetBranchAddress( "sbs.hcal.clus_blk.e", cblke_h[f] ); // Array of block energies
    Ch[f]->SetBranchAddress( "sbs.hcal.clus_blk.tdctime", cblktime_h[f] ); // Array of block TDC times
    Ch[f]->SetBranchAddress( "sbs.hcal.clus_blk.atime", cblkatime_h[f] ); // Array of block ADC times
    Ch[f]->SetBranchAddress( "sbs.hcal.clus_blk.again", cblkagain_h[f] ); // Array of block ADC gain coeff
    Ch[f]->SetBranchAddress( "bb.tr.n", &BBtr_n_h[f] );
    Ch[f]->SetBranchAddress( "bb.tr.chi2", BBtr_chi2_h[f] );
    Ch[f]->SetBranchAddress( "bb.tr.px", BBtr_px_h[f] );
    Ch[f]->SetBranchAddress( "bb.tr.py", BBtr_py_h[f] );
    Ch[f]->SetBranchAddress( "bb.tr.pz", BBtr_pz_h[f] );
    Ch[f]->SetBranchAddress( "bb.tr.vz", BBtr_vz_h[f] );
    Ch[f]->SetBranchAddress( "bb.tr.p", BBtr_p_h[f] );
    Ch[f]->SetBranchAddress( "bb.ps.e", &BBps_e_h[f] );
    Ch[f]->SetBranchAddress( "bb.ps.x", &BBps_x_h[f] );
    Ch[f]->SetBranchAddress( "bb.ps.y", &BBps_y_h[f] );
    Ch[f]->SetBranchAddress( "bb.sh.e", &BBsh_e_h[f] );
    Ch[f]->SetBranchAddress( "bb.sh.x", &BBsh_x_h[f] );
    Ch[f]->SetBranchAddress( "bb.sh.y", &BBsh_y_h[f] ); 
    Ch[f]->SetBranchAddress( "bb.hodotdc.clus.tmean", &HODOtmean_h[f] );
    Ch[f]->SetBranchAddress( "e.kine.W2", &ekineW2_h[f] );


  }//end config file for over lh2 field settings

 // Reading all relevant config files for ld2
  for( Int_t f=0; f<nfset_ld2[kIdx]; f++ ){
    Cd[f] = new TChain("T");

    ifstream configfile(configfilename_d[f]);
    TString currentline;
    while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endlist") ){
      if( !currentline.BeginsWith("#") ){
	Cd[f]->Add(currentline);
	cout << "Loading file: " << currentline << ".." << endl;
      }    
    }
    //TCut globalcut = "";
    while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endcut") ){
      if( !currentline.BeginsWith("#") ){
	globalcut_d[f] += currentline;
      }    
    }
    while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("#") ){
      TObjArray *tokens = currentline.Tokenize(" ");
      Int_t ntokens = tokens->GetEntries();
      if( ntokens>1 ){
	TString skey = ( (TObjString*)(*tokens)[0] )->GetString();
	if( skey == "E_e" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  E_e_d[f] = sval.Atof();
	  cout << "Loading beam energy: " << E_e_d[f] << endl;
	}
	if( skey == "HCal_d" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  HCal_d_d[f] = sval.Atof();
	  cout << "Loading HCal distance: " << HCal_d_d[f] << endl;
	}
	if( skey == "HCal_th" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  HCal_th_d[f] = sval.Atof() * TMath::DegToRad();	
	  cout << "Loading HCal angle: " << HCal_th_d[f] << endl;
	}
	if( skey == "BB_th" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  BB_th_d[f] = sval.Atof() * TMath::DegToRad();	
	  cout << "Loading BBCal angle: " << BB_th_d[f] << endl;
	}
	if( skey == "W2_mean" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  W2_mean_d[f] = sval.Atof();
	  cout << "Loading W2 mean cut: " << W2_mean_d[f] << endl;
	}
	if( skey == "W2_sig" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  W2_sig_d[f] = sval.Atof();
	  cout << "Loading W2 sigma cut: " << W2_sig_d[f] << endl;
	}
	if( skey == "dx0_n" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  dx0_n_d[f] = sval.Atof();
	  cout << "Loading x position of neutron spot: " << dx0_n_d[f] << endl;
	}
	if( skey == "dx0_p" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  dx0_p_d[f] = sval.Atof();
	  cout << "Loading y position of proton spot: " << dx0_p_d[f] << endl;
	}
	if( skey == "dy0" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  dy0_d[f] = sval.Atof();
	  cout << "Loading y position of both hadron spots: " << dy0_d[f] << endl;
	}
	if( skey == "dx_sig_n" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  dx_sig_n_d[f] = sval.Atof();
	  cout << "Loading x sigma of neutron spot: " << dx_sig_n_d[f] << endl;
	}
	if( skey == "dx_sig_p" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  dx_sig_p_d[f] = sval.Atof();
	  cout << "Loading x sigma of proton spot: " << dx_sig_p_d[f] << endl;
	}
	if( skey == "dy_sig" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  dy_sig_d[f] = sval.Atof();
	  cout << "Loading y sigma of both hadron spots: " << dy_sig_d[f] << endl;
	}
	if( skey == "atime0" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  atime0_d[f] = sval.Atof();
	  cout << "Loading ADC time mean: " << atime0_d[f] << endl;
	}
	if( skey == "atime_sig" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  atime_sig_d[f] = sval.Atof();
	  cout << "Loading ADC time sigma: " << atime_sig_d[f] << endl;
	}
	if( skey == "useAlshield" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  useAlshield_d[f] = sval.Atoi();
	  cout << "Loading Aluminum absorber option: " << useAlshield_d[f] << endl;
	}
      }
      delete tokens;
    }

    NTevents_d[f] = Cd[f]->GetEntries();
    elist_d[f] = new TEventList(Form("elist_ld2_sbs%d_f%d",kine,fset_ld2[kIdx][f]),Form("Elastic Event List, LD2, SBS%d, Mag%d",kine,fset_ld2[kIdx][f]));
    Cd[f]->Draw(Form(">>elist_ld2_sbs%d_f%d",kine,fset_ld2[kIdx][f]),globalcut_d[f]);
    Nevents_d[f] = elist_d[f]->GetN();

    Cd[f]->SetBranchStatus( "*", 0 );
    Cd[f]->SetBranchStatus( "sbs.hcal.x", 1 );
    Cd[f]->SetBranchStatus( "sbs.hcal.y", 1 );
    Cd[f]->SetBranchStatus( "sbs.hcal.e", 1 );
    Cd[f]->SetBranchStatus( "sbs.hcal.rowblk", 1 );
    Cd[f]->SetBranchStatus( "sbs.hcal.colblk", 1 );
    Cd[f]->SetBranchStatus( "sbs.hcal.nblk", 1 );
    Cd[f]->SetBranchStatus( "sbs.hcal.clus_blk.id", 1 );
    Cd[f]->SetBranchStatus( "sbs.hcal.clus_blk.e", 1 );
    Cd[f]->SetBranchStatus( "sbs.hcal.clus_blk.tdctime", 1 );
    Cd[f]->SetBranchStatus( "sbs.hcal.clus_blk.atime", 1 );
    Cd[f]->SetBranchStatus( "sbs.hcal.clus_blk.again", 1 );
    Cd[f]->SetBranchStatus( "bb.tr.chi2", 1 );
    Cd[f]->SetBranchStatus( "bb.tr.n", 1 );
    Cd[f]->SetBranchStatus( "bb.tr.px", 1 );
    Cd[f]->SetBranchStatus( "bb.tr.py", 1 );
    Cd[f]->SetBranchStatus( "bb.tr.pz", 1 );    
    Cd[f]->SetBranchStatus( "bb.tr.vz", 1 );
    Cd[f]->SetBranchStatus( "bb.tr.p", 1 );
    Cd[f]->SetBranchStatus( "bb.ps.e", 1 );
    Cd[f]->SetBranchStatus( "bb.ps.x", 1 );
    Cd[f]->SetBranchStatus( "bb.ps.y", 1 );
    Cd[f]->SetBranchStatus( "bb.sh.e", 1 );
    Cd[f]->SetBranchStatus( "bb.sh.x", 1 );
    Cd[f]->SetBranchStatus( "bb.sh.y", 1 );
    Cd[f]->SetBranchStatus( "bb.hodotdc.clus.tmean", 1 );
    Cd[f]->SetBranchStatus( "e.kine.W2", 1 );

    Cd[f]->SetBranchAddress( "sbs.hcal.x", &HCALx_d[f] );
    Cd[f]->SetBranchAddress( "sbs.hcal.y", &HCALy_d[f] );
    Cd[f]->SetBranchAddress( "sbs.hcal.e", &HCALe_d[f] );
    Cd[f]->SetBranchAddress( "sbs.hcal.rowblk", &crow_d[f] );
    Cd[f]->SetBranchAddress( "sbs.hcal.colblk", &ccol_d[f] );
    Cd[f]->SetBranchAddress( "sbs.hcal.nblk", &nblk_d[f] ); // Total number of blocks in highest E clus
    Cd[f]->SetBranchAddress( "sbs.hcal.clus_blk.id", cblkid_d[f] ); // kNcell-1 index for each block
    Cd[f]->SetBranchAddress( "sbs.hcal.clus_blk.e", cblke_d[f] ); // Array of block energies
    Cd[f]->SetBranchAddress( "sbs.hcal.clus_blk.tdctime", cblktime_d[f] ); // Array of block TDC times
    Cd[f]->SetBranchAddress( "sbs.hcal.clus_blk.atime", cblkatime_d[f] ); // Array of block ADC times
    Cd[f]->SetBranchAddress( "sbs.hcal.clus_blk.again", cblkagain_d[f] ); // Array of block ADC gain coeff
    Cd[f]->SetBranchAddress( "bb.tr.n", &BBtr_n_d[f] );
    Cd[f]->SetBranchAddress( "bb.tr.chi2", BBtr_chi2_d[f] );
    Cd[f]->SetBranchAddress( "bb.tr.px", BBtr_px_d[f] );
    Cd[f]->SetBranchAddress( "bb.tr.py", BBtr_py_d[f] );
    Cd[f]->SetBranchAddress( "bb.tr.pz", BBtr_pz_d[f] );
    Cd[f]->SetBranchAddress( "bb.tr.vz", BBtr_vz_d[f] );
    Cd[f]->SetBranchAddress( "bb.tr.p", BBtr_p_d[f] );
    Cd[f]->SetBranchAddress( "bb.ps.e", &BBps_e_d[f] );
    Cd[f]->SetBranchAddress( "bb.ps.x", &BBps_x_d[f] );
    Cd[f]->SetBranchAddress( "bb.ps.y", &BBps_y_d[f] );
    Cd[f]->SetBranchAddress( "bb.sh.e", &BBsh_e_d[f] );
    Cd[f]->SetBranchAddress( "bb.sh.x", &BBsh_x_d[f] );
    Cd[f]->SetBranchAddress( "bb.sh.y", &BBsh_y_d[f] ); 
    Cd[f]->SetBranchAddress( "bb.hodotdc.clus.tmean", &HODOtmean_d[f] );
    Cd[f]->SetBranchAddress( "e.kine.W2", &ekineW2_d[f] );

  }//end config file for over ld2 field settings

  cout << endl << endl << "Setup parameters loaded and chains linked." << endl;

  // Create stopwatch to track processing time
  TStopwatch *sw = new TStopwatch();
  sw->Start();
  
  // Declare outfile
  TFile *fout;
  if( qreplay ){
    fout = new TFile( Form("reporting/timing/qtreplay_sbs%d.root", kine ), "RECREATE" );
  }else{
    fout = new TFile( Form("reporting/timing/tCalOut_sbs%d.root", kine ), "RECREATE" );
  }

  // Initialize vectors and arrays
  Double_t TDCoffsets[kNcell] = {0.0};
  Double_t TDCdeviations[kNcell] = {0.0};
  Double_t ADCtoffsets[kNcell] = {0.0};
  Double_t TDCsig[kNcell] = {0.0};
  Double_t ADCtsig[kNcell] = {0.0};
  Double_t TvseCorr[kNcell] = {0.0};

  Int_t lh2Events = 0;
  for( Int_t f=0; f<nfset_lh2[kIdx]; f++ ){
    lh2Events+=Nevents_h[f];
  }
  Int_t ld2Events = 0;
  for( Int_t f=0; f<nfset_ld2[kIdx]; f++ ){
    ld2Events+=Nevents_d[f];
  }
  
  cout << "Opened many trees with " << lh2Events << " LH2 events passing global cuts." << endl;
  cout << "Opened many trees with " << ld2Events << " LD2 events passing global cuts." << endl << endl;
  cout << "Total available events to calibrate = " << lh2Events + ld2Events << "." << endl << endl;

  report << "Opened many trees with " << lh2Events << " LH2 events passing global cuts." << endl;
  report << "Opened many trees with " << ld2Events << " LD2 events passing global cuts." << endl << endl;
  report << "Total available events to calibrate = " << lh2Events + ld2Events << "." << endl << endl;

  //Declare matrices for chi-square min calibration scheme and keep track of calibrated events
  Int_t TNEV_h = 0;
  Int_t TNEV_d = 0;
  Int_t NEV[kNcell] = {0};
  
  //Declare limits on histograms to avoid fitting mistakes
  Double_t tllim = -150.;
  Double_t tulim = 0.;
  Double_t ttotb = 2*abs(tllim-tulim);
  Double_t allim = 0.;
  Double_t aulim = 150.;
  Double_t atotb = 2*abs(allim-aulim);

  //Declare physics histograms for match to MC
  TH1D *hpp_p = new TH1D( "hpp_p", "Expected proton momentum from BB; GeV", 100, 0, 10 );
  TH1D *hpp_n = new TH1D( "hpp_n", "Expected neutron momentum from BB; GeV", 100, 0, 10 );

  //Declare diagnostic histograms (as sparse as possible)
  TH2D *hdx_mag_h = new TH2D( "hdx_mag_h", "Delta X vs Field Setting (LH2); field (percent); x_{HCAL} - x_{exp} (m)", 21, 0, 105, 250, -4, 3 );
  TH2D *hdy_mag_h = new TH2D( "hdy_mag_h", "Delta Y vs Field Setting (LH2); field (percent); y_{HCAL} - y_{exp} (m)", 21, 0, 105, 250, -1.25, 1.25 );
  TH2D *hdx_mag_d = new TH2D( "hdx_mag_d", "Delta X vs Field Setting (LD2); field (percent); x_{HCAL} - x_{exp} (m)", 21, 0, 105, 250, -4, 3 );
  TH2D *hdy_mag_d = new TH2D( "hdy_mag_d", "Delta Y vs Field Setting (LD2); field (percent); y_{HCAL} - y_{exp} (m)", 21, 0, 105, 250, -1.25, 1.25 );
  TH1D *htdc = new TH1D("htdc","TDC (All Channels); ns",ttotb,tllim,tulim);
  TH1D *htdc_corr = new TH1D("htdc_corr","TDC TWalk/Hodo Corrected TW (All Channels); ns",ttotb,tllim,tulim);
  TH1D *htdc_allcorr = new TH1D("htdc_allcorr","TDC TWalk/Hodo Corrected TW/TOF (All Channels); ns",ttotb,tllim,tulim);
  TH1D *hadct = new TH1D("hadct","ADCt (All Channels); ns",atotb,allim,aulim);
  TH1D *hadct_corr = new TH1D("hadct_corr","ADCt TWalk/Hodo Corrected (All Channels); ns",atotb,allim,aulim);
  TH2D *ht_ID = new TH2D("ht_ID","TDC;Channel;TDC_{HCAL} (ns)",288,0,288,ttotb,tllim,tulim);
  TH2D *htp_ID = new TH2D("htp_ID","TDC Primary Block;Channel;TDC_{HCAL} (ns)",288,0,288,ttotb,tllim,tulim);
  TH2D *htpDiff_ID = new TH2D("htpDiff_ID","TDC Primary Block - TDC hodo;Channel;TDC_{HCAL}-TDC_{HODO} (ns)",288,0,288,ttotb,tllim,tulim);
  TH2D *htpCorr_ID = new TH2D("htpCorr_ID","TDC Primary Block, Corrected;Channel;TDC_{HCAL}-TDC_{HODO}-Corr_{TW} (ns)",288,0,288,ttotb,tllim,tulim);
  TH2D *htDiff_ID = new TH2D("htDiff_ID",";Channel;TDC_{HCAL}-TDC_{HODO} (ns)",288,0,288,ttotb,tllim,tulim);
  TH2D *htCorr_ID = new TH2D("htCorr_ID",";Channel;TDC_{HCAL}-TDC_{HODO} w/TWalk (ns)",288,0,288,ttotb,tllim,tulim);
  TH2D *htAllCorr_ID = new TH2D("htAllCorr_ID",";Channel;TDC_{HCAL}-TDC_{HODO} w/TWalk w/TOF (ns)",288,0,288,ttotb,tllim,tulim);
  TH2D *ha_ID = new TH2D("ha_ID","ADCt;Channel;ADCtime_{HCAL} (ns)",288,0,288,atotb,allim,aulim);
  TH2D *hap_ID = new TH2D("hap_ID","ADCt Primary Block;Channel;ADCtime_{HCAL} (ns)",288,0,288,atotb,allim,aulim);
  TH2D *hapDiff_ID = new TH2D("hapDiff_ID","ADCt Primary Block - TDC hodo;Channel;ADCtime_{HCAL} - TDC_{HODO} (ns)",288,0,288,atotb,allim,aulim);
  TH2D *hapCorr_ID = new TH2D("hapCorr_ID","ADCt Primary Block, Corrected;Channel;ADCtime_{HCAL}-TDC_{HODO}-Corr_{TW} (ns)",288,0,288,atotb,allim,aulim);
  TH2D *haDiff_ID = new TH2D("haDiff_ID",";Channel;ADCtime_{HCAL}-TDC_{HODO} (ns)",288,0,288,atotb,allim,aulim);
  TH2D *haCorr_ID = new TH2D("haCorr_ID",";Channel;ADCtime_{HCAL}-TDC_{HODO}-Corr_{TW} (ns)",288,0,288,atotb,allim,aulim);

  //Declare self timing histograms
  TH1D *hclusmean = new TH1D( "hclusmean", "HCal Cluster Mean Time (Primary Block Reftime, nblk>1)", 1000, -50, 50 );
  TH2D *hclusmeanID = new TH2D( "hclusmeanID", "HCal Cluster Mean Time, (Primary Block Reftime/ID, nblk>1), vs ID", 288, 0, 288, 1000, -50, 50 );
  TH1D *hclusdiff = new TH1D( "hclusdiff", "HCal Cluster Seed - Block Diff, (Primary Block Reftime, nblk>1)", 1000, -50, 50 );  
  TH2D *hclusdiffID = new TH2D( "hclusdiffID", "HCal Cluster Seed - Block Diff, (Primary Block Reftime, nblk>1), vs ID", 288, 0, 288, 1000, -50, 50 );

  //TW fit diagnostic histos
  TH1D *htdcP0 = new TH1D("htdcP0","TDC P0, P0exp(P1*x)+P2",60,0,30);
  TH1D *htdcP1 = new TH1D("htdcP1","TDC P1, P0exp(P1*x)+P2",50,0,0.5);
  TH1D *htdcP2 = new TH1D("htdcP2","TDC P2, P0exp(P1*x)+P2",250,-200,50);
  TH1D *hadctP0 = new TH1D("hadctP0","ADCt P0, P0exp(P1*x)+P2",60,0,10);
  TH1D *hadctP1 = new TH1D("hadctP1","ADCt P1, P0exp(P1*x)+P2",50,0,0.5);
  TH1D *hadctP2 = new TH1D("hadctP2","ADCt P2, P0exp(P1*x)+P2",250,-50,200);

  //Construct many histograms to extract timewalk fit parameters
  TH2D *htdcVe[kNcell];
  TH2D *hadctVe[kNcell];
  for( Int_t i=0; i<kNcell; i++ ){
    htdcVe[i] = new TH2D(Form("htdcVe_bl%d",i),Form(";E_{bl%d} (GeV);TDC_{HCAL} (ns)",i),1000,0.0,500,250,-200,50);
    hadctVe[i] = new TH2D(Form("hadctVe_bl%d",i),Form(";E_{bl%d} (GeV);ADCt_{HCAL} (ns)",i),1000,0.0,500,250,-50,200);
  }

  //Loop over events
  cout << "Looping over hydrogen data.." << endl;
  report << "Looping over hydrogen data.." << endl;
  if(pass0) gROOT->ProcessLine( "gErrorIgnoreLevel = 6001" ); //Suppress error output to avoid undetectable sbs.hcal.clus_blk.again for pass0

  //Loop over all hydrogen data
  for( Int_t f=0; f<nfset_lh2[kIdx]; f++ ){
    
    Int_t hfieldset = fset_lh2[kIdx][f]; //Check B field setting
    bool tofready = tofreplay && hfieldset==tofB[kIdx]; //Check if TOF corrections can be applied for these data

    //Declare energy loss parameters for beam going through the target
    Double_t pBeam = E_e_h[f]/(1.+E_e_h[f]/M_p*(1.-cos(BB_th_h[f])));
    Double_t Eloss_outgoing = celldiameter/2.0/sin(BB_th_h[f]) * rho_tgt * dEdx_tgt; //Mean energy loss of the beam prior to the scattering, approximately 1 MeV, could correct further with raster position (likely negligable)
    
    for( Long64_t nevent = 1; nevent <Nevents_h[f]; nevent++){

      if ( nevent%10000==0 ) cout << "LH2 kinematic " << kine << " at field " << hfieldset << "% , entry: " << nevent << "/" << Nevents_h[f] << ". Total events gathered for calibration: " << TNEV_h << " \r";
      cout.flush();

      if ( nevent%100000==0 ) report << "LH2 kinematic " << kine << " at field " << hfieldset << "% , entry: " << nevent << "/" << Nevents_h[f] << ". Total events gathered for calibration: " << TNEV_h << endl;
      
      Ch[f]->GetEntry( elist_h[f]->GetEntry( nevent ) ); 

      Double_t A[kNcell] = {0.0}; // Array to keep track of ADC values per cell. Outscope on each ev
      Double_t A_oneblock[kNcell] = {0.0}; // Array to keep track of ADC values per cell for one block clusters only. Outscope on each ev
      
      //Correct the beam energy with energy loss in target using vertex position
      Double_t Eloss = (BBtr_vz_h[f][0]+l_tgt/2.0) * rho_tgt * dEdx_tgt + uwallthick_LH2 * rho_Al * dEdx_Al; //approximately 3 MeV
      Double_t E_corr = E_e_h[f] - Eloss;
      Double_t p_corr = BBtr_p_h[f][0] - Eloss_outgoing; //Neglecting the mass of e'
      Double_t etheta = acos( BBtr_pz_h[f][0]/BBtr_p_h[f][0] );
      Double_t ephi = atan2( BBtr_py_h[f][0], BBtr_px_h[f][0] );

      TVector3 vertex(0,0,BBtr_vz_h[f][0]); // z location of vertex in hall coordinates
      TLorentzVector Pbeam(0,0,E_corr,E_corr); //Mass of e negligable
      TLorentzVector kprime(BBtr_px_h[f][0],BBtr_py_h[f][0],BBtr_pz_h[f][0],BBtr_p_h[f][0]);
      TLorentzVector Ptarg(0,0,0,M_p); // assume proton for both LH2 and LD2 - can refine where useful with exclusive LH2 data set. Likely better to refine first with dxdy spot cuts on both protons and neutrons.
      TLorentzVector q = Pbeam - kprime;
      TLorentzVector PgammaN = Ptarg + q; //(-px, -py, ebeam - pz, Mp + ebeam - p)

      Double_t pel = E_corr/(1.+E_corr/M_p*(1.-cos(etheta)));
      Double_t nu = E_corr - BBtr_p_h[f][0];
      Double_t pp = sqrt(pow(nu,2)+2.*M_p*nu); //momentum of proton
      Double_t pn = sqrt(pow(nu,2)+2.*M_n*nu); //momentum of neutron
      Double_t phinucleon = ephi + PI; //assume coplanarity
      Double_t thetanucleon = acos( (E_corr - BBtr_pz_h[f][0])/pp ); //use elastic constraint on nucleon kinematic
      TVector3 pNhat( sin(thetanucleon)*cos(phinucleon),sin(thetanucleon)*sin(phinucleon),cos(thetanucleon));

      //Define HCal coordinate system
      TVector3 HCAL_zaxis(sin(-HCal_th_h[f]),0,cos(-HCal_th_h[f]));
      TVector3 HCAL_xaxis(0,-1,0);
      TVector3 HCAL_yaxis = HCAL_zaxis.Cross(HCAL_xaxis).Unit();	
      TVector3 HCAL_origin = HCal_d_h[f] * HCAL_zaxis + hcalheight * HCAL_xaxis;

      //Define intersection points for hadron vector
      Double_t sintersect = ( HCAL_origin - vertex ).Dot( HCAL_zaxis ) / ( pNhat.Dot( HCAL_zaxis ) );
      TVector3 HCAL_intersect = vertex + sintersect * pNhat;
      Double_t yexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_yaxis );
      Double_t xexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_xaxis );
      Double_t E_ep = sqrt( pow(M_e,2) + pow(BBtr_p_h[f][0],2) ); // Obtain the scattered electron energy
      Double_t p_ep = BBtr_p_h[f][0];
      Double_t Q2 = 2*E_corr*E_ep*( 1-(BBtr_pz_h[f][0]/p_ep) ); // Obtain Q2 from beam energy, outgoing electron energy, and momenta
	
      Double_t W = PgammaN.M();
      Double_t W2 = ekineW2_h[f];

      //Use the electron kinematics to predict the proton momentum assuming elastic scattering on free proton at rest:
      Double_t E_pp = nu+M_p; // Get energy of the proton
      Double_t Enucleon = sqrt(pow(pp,2)+pow(M_p,2)); // Check on E_pp, same
      Double_t KE_p = nu; //For elastics
      Double_t dx = HCALx_h[f] - xexpect_HCAL;
      Double_t dy = HCALy_h[f] - yexpect_HCAL;
      
      /////////////////////////////////////////////////
      //Primary W2 cut on elastics and HCal active area
      if( fabs(W2-W2_mean_h[f])>W2_sig_h[f] ) continue; //Observed mean W2 cut on elastic peak
      elasYield++; //Events that pass the above cuts constitute elastics
      /////////////////////////////////////////////////

      //Calculate/declare new variables for analysis
      Double_t SFrac = HCALe_h[f]/KE_p;
      //Double_t E_exp = KE_p*sampFrac[kIdx];
      Int_t rec_row = ( HCALx_h[f] - HCal_Xmin )/HCal_divx;
      Int_t rec_col = ( HCALy_h[f] - HCal_Ymin )/HCal_divy;
      Int_t rec_cell = rec_row*kNcols + rec_col;

      hdx_mag_h->Fill( hfieldset, dx );
      hdy_mag_h->Fill( hfieldset, dy );
      
      //////////////////////////////////////////////////////////////////
      //Cut on dx and dy.
      bool pass_y = abs(dy-dy0_h[f])<ifac*dy_sig_h[f];
      if( !pass_y ) continue;
      bool pass_p = abs(dx-dx0_p_h[f])<ifac*dx_sig_p_h[f];
      bool pass_n = abs(dx-dx0_n_h[f])<ifac*dx_sig_n_h[f]; //Redundant for LH2
      bool isproton = pass_p && !pass_n;
      bool isneutron = pass_n && !pass_p;
      bool isamb = pass_p && pass_n;
      if( !pass_p && !pass_n ) continue; //Cut on both n and p spots for each event, cannot know which apriori
      //////////////////////////////////////////////////////////////////
      if( isproton ) hpp_p->Fill( pp );
      if( isneutron ) hpp_n->Fill( pn );
      

      //Hodo cluster mean time
      Double_t hodot = HODOtmean_h[f];

      //Fill some histograms with only the primary block timing
      Int_t pblkid = int(cblkid_h[f][0])-1;
      Double_t pTDC = cblktime_h[f][0];
      Double_t pADCt = cblkatime_h[f][0];
      Double_t pTWcorr = 0.; //primary block timewalk correction
      Double_t pADCtcorr = 0.;
      Double_t pblke = 0.;
      if( pass0 ){
	pblke = cblke_h[f][0]/gOldConst_pass0[pblkid]*gConst_iter1[pblkid];
      }else{
	pblke = cblke_h[f][0]/cblkagain_h[f][0]*gConst_iter1[pblkid];
      }
      if( qreplay ){
	pTWcorr = otdcP0[pblkid]*exp(-otdcP1[pblkid]*pblke);
	pTDC = cblktime_h[f][0] + oldTDCoffsets[pblkid]*TDCCalib - calTDCoffsets[pblkid]*TDCCalib;
      }
      htp_ID->Fill( pblkid, pTDC );
      htpDiff_ID->Fill( pblkid, pTDC-hodot );
      htpCorr_ID->Fill( pblkid, pTDC-hodot-pTWcorr );
      
      if( qreplay ){
	pADCtcorr = oadctP0[pblkid]*exp(-oadctP1[pblkid]*pblke);
	pADCt = cblkatime_h[f][0] + oldADCtoffsets[pblkid] - calADCtoffsets[pblkid];
      }
      hap_ID->Fill( pblkid, pADCt );
      hapDiff_ID->Fill( pblkid, pADCt-hodot );
      hapCorr_ID->Fill( pblkid, pADCt-hodot-pADCtcorr );

      //Loop over primary cluster
      Double_t blkseed = 0.;
      Double_t tavg = 0;
      Int_t nblk = (int)nblk_h[f];
      for( Int_t blk = 0; blk<nblk; blk++ ){
	Int_t blkid = int(cblkid_h[f][blk])-1; //-1 necessary since sbs.hcal.clus_blk.id ranges from 1 to 288 (different than just about everything else)
	Double_t blkatime = cblkatime_h[f][blk];
	Double_t blktime = cblktime_h[f][blk];
	Double_t blke = 0.;
	if( pass0 ){
	  blke = cblke_h[f][blk]/gOldConst_pass0[blkid]*gConst_iter1[blkid];
	}else{
	  blke = cblke_h[f][blk]/cblkagain_h[f][blk]*gConst_iter1[blkid];
	}
	
	Double_t TOFcorr = 0; //time of flight correction
	Double_t TWcorr = 0.; //timewalk correction
	Double_t ADCtcorr = 0.;
	if( qreplay ){
	  TWcorr = otdcP0[blkid]*exp(-otdcP1[blkid]*blke); //timewalk correction applied
	  if( isproton && tofready){ //Check if tofreplay, in proton dxdy, not ambiguous, if field set matches current MC fits
	    Double_t protmom = pp;
	    if( pp>pulim[kIdx] ) protmom = pulim[kIdx]; //do not pass corrections for momenta > fit limit in MC
	    TOFcorr = TOFfitp_p[kIdx][0]+TOFfitp_p[kIdx][1]*protmom+TOFfitp_p[kIdx][2]*pow(protmom,2)+TOFfitp_p[kIdx][3]*pow(protmom,3); //proton correction, only overwrite if data matches currently available MC fit results
	  }
	  if( isneutron && tofready){ //Check if tofreplay, in neutron dxdy, not ambiguous, if field set matches current MC fits
	    Double neutmom = pn;
	    if( pn>pulim[kIdx] ) neutmom = pulim[kIdx]; //do not pass corrections for momenta > fit limit in MC
	    TOFcorr = TOFfitp_n[kIdx][0]+TOFfitp_n[kIdx][1]*neutmom+TOFfitp_n[kIdx][2]*pow(neutmom,2)+TOFfitp_n[kIdx][3]*pow(neutmom,3); //proton correction, only overwrite if data matches currently available MC fit results
	  }
	  ADCtcorr = oadctP0[blkid]*exp(-oadctP1[blkid]*blke);
	  Double_t blkatime_new = blkatime + oldADCtoffsets[blkid] - calADCtoffsets[blkid]; //reverse replay "+", then qreplay "-"
	  Double_t blktime_new = blktime + oldTDCoffsets[blkid]*TDCCalib - calTDCoffsets[blkid]*TDCCalib; //reverse replay "+", then qreplay "-", but first note that both offset values are in tdc units and need conversion to ns with tdccalib
	  if( blkatime_new>0 ){
	    hadct->Fill( blkatime_new );
	    hadct_corr->Fill( blkatime_new - hodot - ADCtcorr );
	    ha_ID->Fill( blkid, blkatime_new );
	    haDiff_ID->Fill( blkid, blkatime_new - hodot );
	    haCorr_ID->Fill( blkid, blkatime_new - hodot - ADCtcorr );
	    hadctVe[blkid]->Fill( blke*1000, blkatime_new ); //Convert to MeV for clarity in plots
	  }
	  if( blktime_new<10000 && blktime_new>-500){ //ensure that good tdc times exist on tree
	    htdc->Fill( blktime_new );
	    htdc_corr->Fill( blktime_new - hodot - TWcorr );
	    ht_ID->Fill( blkid, blktime_new );
	    htDiff_ID->Fill( blkid, blktime_new - hodot );
	    htCorr_ID->Fill( blkid, blktime_new - hodot - TWcorr );
	    htdcVe[blkid]->Fill( blke*1000, blktime_new ); //Convert to MeV for clarity in plots
	    if( tofready ){
	      htdc_allcorr->Fill( blktime_new - hodot - TWcorr - TOFcorr );
	      htAllCorr_ID->Fill( blkid, blktime_new - hodot - TWcorr - TOFcorr);
	    }

	    //Get self timing histograms
	    if( blk==0 ){ //Just set the reference time with the time of the primary block
	      blkseed = blktime_new;
	    }else if( blkseed != 0. ){
	      
	      Double_t tcdiff = blktime_new-blkseed;
	      
	      tavg += tcdiff;
	      hclusdiffID->Fill( blkid, tcdiff );
	      hclusdiff->Fill( tcdiff );
	    }
	    
	  }
	}else{
	  if( blkatime>0 ){
	    hadct->Fill( blkatime );
	    hadct_corr->Fill( blkatime - hodot - ADCtcorr );
	    ha_ID->Fill( blkid, blkatime );
	    haDiff_ID->Fill( blkid, blkatime - hodot );
	    haCorr_ID->Fill( blkid, blkatime - hodot - ADCtcorr );
	    hadctVe[blkid]->Fill( blke*1000, blkatime ); //Convert to MeV for clarity in plots
	  }
	  if( blktime<10000 && blktime>-500){ //ensure that good tdc times exist on tree
	    htdc->Fill( blktime );
	    htdc_corr->Fill( blktime - hodot - TWcorr );
	    ht_ID->Fill( blkid, blktime );
	    htDiff_ID->Fill( blkid, blktime - hodot );
	    htCorr_ID->Fill( blkid, blktime - hodot - TWcorr );
	    htdcVe[blkid]->Fill( blke*1000, blktime ); //Convert to MeV for clarity in plots
	    
	    //Get self timing histograms (no recovery on bad tdc pblks yet)
	    if( blk==0 ){ //Just set the reference time with the time of the primary block
	      blkseed = blktime;
	    }else if( blkseed != 0. ){
	      
	      Double_t tcdiff = blktime-blkseed;
	      
	      tavg += tcdiff;
	      hclusdiffID->Fill( blkid, tcdiff );
	      hclusdiff->Fill( tcdiff );
	    }
	  }
	}
	
	NEV[blkid]++;
	TNEV_h++;

      } //loop over blocks in primary cluster

      //Finish getting self timing mean time
      if( nblk<2 ) continue;
      tavg /= (nblk-1);
      if( tavg!=0 ){
	hclusmeanID->Fill( pblkid, tavg );
	hclusmean->Fill( tavg );
      }
      
    } //loop over hydrogen events
  } //loop for lh2

  if( !qreplay ){
    Int_t cell = 0;

    cout << endl << "Number of events available for calibration from hydrogen alone: " << endl;
    report << endl << "Number of events available for calibration from hydrogen alone: " << endl;

    for( Int_t r = 0; r<kNrows; r++){
      for( Int_t c = 0; c<kNcols; c++){
	cout << NEV[cell] << "  ";
	report << NEV[cell] << "  ";
	cell++;
      }
      cout << endl;
      report << endl;
    }
  
    cout << endl;
    report << endl;

    usleep( 3*second ); //Give some time for review of step
  }

  cout << "Looping over deuterium data.." << endl;
  report << "Looping over deuterium data.." << endl;

  //Loop over all deuterium data
  for( Int_t f=0; f<nfset_ld2[kIdx]; f++ ){

    Int_t dfieldset = fset_ld2[kIdx][f];

    //Declare energy loss parameters for beam going through the target
    Double_t pBeam = E_e_d[f]/(1.+E_e_d[f]/M_p*(1.-cos(BB_th_d[f])));
    Double_t Eloss_outgoing = celldiameter/2.0/sin(BB_th_d[f]) * rho_tgt * dEdx_tgt; //Mean energy loss of the beam prior to the scattering, approximately 1 MeV, could correct further with raster position (likely negligable)

    for( Long64_t nevent = 1; nevent <Nevents_d[f]; nevent++){

      if ( nevent%10000==0 ) cout << "LD2 kinematic " << kine << " at field " << dfieldset << "%, entry: " << nevent << "/" << Nevents_d[f] << ". Total events gathered for calibration: " << TNEV_d << " \r";
      cout.flush();
      
      if ( nevent%100000==0 ) report << "LD2 kinematic " << kine << " at field " << dfieldset << "%, entry: " << nevent << "/" << Nevents_d[f] << ". Total events gathered for calibration: " << TNEV_d << endl;

      Cd[f]->GetEntry( elist_d[f]->GetEntry( nevent ) ); 

      Double_t A[kNcell] = {0.0}; // Array to keep track of ADC values per cell. Outscope on each ev
      Double_t A_oneblock[kNcell] = {0.0}; // Array to keep track of ADC values per cell for one block clusters only. Outscope on each ev
         
      //Correct the beam energy with energy loss in target using vertex position
      Double_t Eloss = (BBtr_vz_d[f][0]+l_tgt/2.0) * rho_tgt * dEdx_tgt + uwallthick_LH2 * rho_Al * dEdx_Al; //approximately 3 MeV
      Double_t E_corr = E_e_d[f] - Eloss;
      Double_t p_corr = BBtr_p_d[f][0] - Eloss_outgoing; //Neglecting the mass of e'
      Double_t etheta = acos( BBtr_pz_d[f][0]/BBtr_p_d[f][0] );
      Double_t ephi = atan2( BBtr_py_d[f][0], BBtr_px_d[f][0] );

      TVector3 vertex(0,0,BBtr_vz_d[f][0]); // z location of vertex in hall coordinates
      TLorentzVector Pbeam(0,0,E_corr,E_corr); //Mass of e negligable
      TLorentzVector kprime(BBtr_px_d[f][0],BBtr_py_d[f][0],BBtr_pz_d[f][0],BBtr_p_d[f][0]);
      TLorentzVector Ptarg(0,0,0,M_p); // assume proton for both LH2 and LD2 - can refine where useful with exclusive LH2 data set. Likely better to refine first with dxdy spot cuts on both protons and neutrons.
      TLorentzVector q = Pbeam - kprime;
      TLorentzVector PgammaN = Ptarg + q; //(-px, -py, ebeam - pz, Mp + ebeam - p)

      Double_t pel = E_corr/(1.+E_corr/M_p*(1.-cos(etheta)));
      Double_t nu = E_corr - BBtr_p_d[f][0];
      Double_t pp = sqrt(pow(nu,2)+2.*M_p*nu); //momentum of proton
      Double_t pn = sqrt(pow(nu,2)+2.*M_n*nu); //momentum of neutron
      Double_t phinucleon = ephi + PI; //assume coplanarity
      Double_t thetanucleon = acos( (E_corr - BBtr_pz_d[f][0])/pp ); //use elastic constraint on nucleon kinematics
	
      TVector3 pNhat( sin(thetanucleon)*cos(phinucleon),sin(thetanucleon)*sin(phinucleon),cos(thetanucleon));
      //Define HCal coordinate system
      TVector3 HCAL_zaxis(sin(-HCal_th_d[f]),0,cos(-HCal_th_d[f]));
      TVector3 HCAL_xaxis(0,-1,0);
      TVector3 HCAL_yaxis = HCAL_zaxis.Cross(HCAL_xaxis).Unit();	
      TVector3 HCAL_origin = HCal_d_d[f] * HCAL_zaxis + hcalheight * HCAL_xaxis;

      //Define intersection points for hadron vector
      Double_t sintersect = ( HCAL_origin - vertex ).Dot( HCAL_zaxis ) / ( pNhat.Dot( HCAL_zaxis ) );
      TVector3 HCAL_intersect = vertex + sintersect * pNhat;
      Double_t yexpect_dCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_yaxis );
      Double_t xexpect_dCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_xaxis );
      Double_t E_ep = sqrt( pow(M_e,2) + pow(BBtr_p_d[f][0],2) ); // Obtain the scattered electron energy	
      Double_t p_ep = BBtr_p_d[f][0];
      Double_t Q2 = 2*E_corr*E_ep*( 1-(BBtr_pz_d[f][0]/p_ep) ); // Obtain Q2 from beam energy, outgoing electron energy, and momenta
	
      Double_t W = PgammaN.M();
      Double_t W2 = ekineW2_d[f];

      //Use the electron kinematics to predict the proton momentum assuming elastic scattering on free proton at rest:
      Double_t E_pp = nu+M_p; // Get energy of the proton
      Double_t Enucleon = sqrt(pow(pp,2)+pow(M_p,2)); // Check on E_pp, same
      Double_t KE_p = nu; //For elastics
      Double_t dx = HCALx_d[f] - xexpect_dCAL;
      Double_t dy = HCALy_d[f] - yexpect_dCAL;
      
      /////////////////////////////////////////////////
      //Primary W2 cut on elastics and HCal active area
      if( fabs(W2-W2_mean_d[f])>W2_sig_d[f] ) continue; //Observed mean W2 cut on elastic peak
      elasYield++; //Events that pass the above cuts constitute elastics
      /////////////////////////////////////////////////

      //Calculate/declare new variables for analysis
      Double_t SFrac = HCALe_d[f]/KE_p;
      //Double_t E_exp = KE_p*sampFrac[kIdx];
      Int_t rec_row = ( HCALx_d[f] - HCal_Xmin )/HCal_divx;
      Int_t rec_col = ( HCALy_d[f] - HCal_Ymin )/HCal_divy;
      Int_t rec_cell = rec_row*kNcols + rec_col;
      
      hdx_mag_d->Fill( dfieldset, dx );
      hdy_mag_d->Fill( dfieldset, dy );

      //////////////////////////////////////////////////////////////////
      //Cut on dx and dy.
      bool pass_y = abs(dy-dy0_d[f])<ifac*dy_sig_d[f];
      if( !pass_y ) continue;
      bool pass_p = abs(dx-dx0_p_d[f])<ifac*dx_sig_p_d[f];
      bool pass_n = abs(dx-dx0_n_d[f])<ifac*dx_sig_n_d[f];
      bool isproton = pass_p && !pass_n;
      bool isneutron = pass_n && !pass_p;
      bool isamb = pass_p && pass_n;
      if( !pass_p && !pass_n ) continue; //Cut on both n and p spots for each event, cannot know which apriori      
      //////////////////////////////////////////////////////////////////
      if( isproton ) hpp_p->Fill( pp );
      if( isneutron ) hpp_n->Fill( pn );

      //Hodo cluster mean time
      Double_t hodot = HODOtmean_d[f];

      //Fill some histograms with only the primary block timing
      Int_t pblkid = int(cblkid_d[f][0])-1;
      Double_t pTDC = cblktime_d[f][0];
      Double_t pADCt = cblkatime_d[f][0];
      Double_t pTWcorr = 0.;
      Double_t pADCtcorr = 0.;
      Double_t pblke = 0.;
      if( pass0 ){
	pblke = cblke_d[f][0]/gOldConst_pass0[pblkid]*gConst_iter1[pblkid];
      }else{
	pblke = cblke_d[f][0]/cblkagain_d[f][0]*gConst_iter1[pblkid];
      }
      if( qreplay ){
	pTWcorr = otdcP0[pblkid]*exp(-otdcP1[pblkid]*pblke);
	pTDC = cblktime_d[f][0] + oldTDCoffsets[pblkid]*TDCCalib - calTDCoffsets[pblkid]*TDCCalib;
      }
      htp_ID->Fill( pblkid, pTDC );
      htpDiff_ID->Fill( pblkid, pTDC-hodot );
      htpCorr_ID->Fill( pblkid, pTDC-hodot-pTWcorr );
      
      if( qreplay ){
	pADCtcorr = oadctP0[pblkid]*exp(-oadctP1[pblkid]*pblke);
	pADCt = cblkatime_d[f][0] + oldADCtoffsets[pblkid] - calADCtoffsets[pblkid];
      }
      hap_ID->Fill( pblkid, pADCt );
      hapDiff_ID->Fill( pblkid, pADCt-hodot );
      hapCorr_ID->Fill( pblkid, pADCt-hodot-pADCtcorr );

      //Loop over primary cluster
      Double_t blkseed = 0.;
      Double_t tavg = 0;
      Int_t nblk = (int)nblk_d[f];
      for( Int_t blk = 0; blk<nblk; blk++ ){
	Int_t blkid = int(cblkid_d[f][blk])-1; //-1 necessary since sbs.hcal.clus_blk.id ranges from 1 to 288 (different than just about everything else)
	Double_t blkatime = cblkatime_d[f][blk];
	Double_t blktime = cblktime_d[f][blk];
	Double_t blke = 0.;
	if( pass0 ){
	  blke = cblke_d[f][blk]/gOldConst_pass0[blkid]*gConst_iter1[blkid];
	}else{
	  blke = cblke_d[f][blk]/cblkagain_d[f][blk]*gConst_iter1[blkid];
	}

	Double_t TOFcorr = 0.; //time of flight correction
	Double_t TWcorr = 0.; //timewalk correction
	Double_t ADCtcorr = 0.;
	if( qreplay ){
	  TWcorr = otdcP0[blkid]*exp(-otdcP1[blkid]*blke);
	  if( isproton && tofready){ //Check if tofreplay, in proton dxdy, not ambiguous, if field set matches current MC fits
	    Double_t protmom = pp;
	    if( pp>pulim[kIdx] ) protmom = pulim[kIdx]; //do not pass corrections for momenta > fit limit in MC
	    TOFcorr = TOFfitp_p[kIdx][0]+TOFfitp_p[kIdx][1]*protmom+TOFfitp_p[kIdx][2]*pow(protmom,2)+TOFfitp_p[kIdx][3]*pow(protmom,3); //proton correction, only overwrite if data matches currently available MC fit results
	  }
	  if( isneutron && tofready){ //Check if tofreplay, in neutron dxdy, not ambiguous, if field set matches current MC fits
	    Double neutmom = pn;
	    if( pn>pulim[kIdx] ) neutmom = pulim[kIdx]; //do not pass corrections for momenta > fit limit in MC
	    TOFcorr = TOFfitp_n[kIdx][0]+TOFfitp_n[kIdx][1]*neutmom+TOFfitp_n[kIdx][2]*pow(neutmom,2)+TOFfitp_n[kIdx][3]*pow(neutmom,3); //proton correction, only overwrite if data matches currently available MC fit results
	  }
	  ADCtcorr = oadctP0[blkid]*exp(-oadctP1[blkid]*blke);
	  Double_t blkatime_new = blkatime + oldADCtoffsets[blkid] - calADCtoffsets[blkid]; //reverse replay "+", then qreplay "-"
	  Double_t blktime_new = blktime + oldTDCoffsets[blkid]*TDCCalib - calTDCoffsets[blkid]*TDCCalib; //reverse replay "+", then qreplay "-", but first note that both offset values are in tdc units and need conversion to ns with tdccalib
	  if( blkatime_new>0 ){
	    hadct->Fill( blkatime_new );
	    hadct_corr->Fill( blkatime_new - hodot - ADCtcorr );
	    ha_ID->Fill( blkid, blkatime_new );
	    haDiff_ID->Fill( blkid, blkatime_new - hodot );
	    haCorr_ID->Fill( blkid, blkatime_new - hodot - ADCtcorr );
	    hadctVe[blkid]->Fill( blke*1000, blkatime_new ); //Convert to MeV for clarity in plots
	  }
	  if( blktime_new<10000 && blktime_new>-500){ //ensure that good tdc times exist on tree
	    htdc->Fill( blktime_new );
	    htdc_corr->Fill( blktime_new - hodot - TWcorr );
	    ht_ID->Fill( blkid, blktime_new );
	    htDiff_ID->Fill( blkid, blktime_new - hodot );
	    htCorr_ID->Fill( blkid, blktime_new - hodot - TWcorr );
	    htdcVe[blkid]->Fill( blke*1000, blktime_new ); //Convert to MeV for clarity in plots
	    if( tofready ){
	      htdc_allcorr->Fill( blktime_new - hodot - TWcorr - TOFcorr );
	      htAllCorr_ID->Fill( blkid, blktime_new - hodot - TWcorr - TOFcorr);
	    }

	    //Get self timing histograms
	    if( blk==0 ){ //Just set the reference time with the time of the primary block
	      blkseed = blktime_new;
	    }else if( blkseed != 0.){
	      Double_t tcdiff = blktime_new-blkseed;
	      
	      tavg += tcdiff;
	      hclusdiffID->Fill( blkid, tcdiff );
	      hclusdiff->Fill( tcdiff );
	    }

	  }	  
	}else{
	  if( blkatime>0 ){
	    hadct->Fill( blkatime );
	    hadct_corr->Fill( blkatime - hodot - ADCtcorr );
	    ha_ID->Fill( blkid, blkatime );
	    haDiff_ID->Fill( blkid, blkatime - hodot );
	    haCorr_ID->Fill( blkid, blkatime - hodot - ADCtcorr );
	    hadctVe[blkid]->Fill( blke*1000, blkatime ); //Convert to MeV for clarity in plots
	  }
	  if( blktime<10000 && blktime>-500){ //ensure that good tdc times exist on tree
	    htdc->Fill( blktime );
	    htdc_corr->Fill( blktime - hodot - TWcorr );
	    ht_ID->Fill( blkid, blktime );
	    htDiff_ID->Fill( blkid, blktime - hodot );
	    htCorr_ID->Fill( blkid, blktime - hodot - TWcorr );
	    htdcVe[blkid]->Fill( blke, blktime ); //Convert to MeV for clarity in plots
	  
	    //Get self timing histograms
	    if( blk==0 ){ //Just set the reference time with the time of the primary block
	      blkseed = blktime;
	    }else if( blkseed != 0. ){
	      
	      Double_t tcdiff = blktime-blkseed;

	      tavg += tcdiff;
	      hclusdiffID->Fill( blkid, tcdiff );
	      hclusdiff->Fill( tcdiff );
	    }
	  }
	}

	NEV[blkid]++;
	TNEV_d++;
	
      } //loop over blocks in primary cluster

      //Finish getting self timing mean time
      if( nblk<2 ) continue;
      tavg /= (nblk-1);
      if( tavg!=0 ){
	hclusmeanID->Fill( pblkid, tavg );
	hclusmean->Fill( tavg );
      }

    } //loop over deuterium events
  } //loop for ld2

  cout << endl << "Fitting TDC/ADCt/TW distributions" << endl;
  report << endl << "Fitting TDC/ADCt/TW distributions" << endl;

  //Make arrays for TW fits
  Double_t TDCvseP0[kNcell] = {0.0};
  Double_t TDCvseP1[kNcell] = {0.0};
  Double_t TDCvseP2[kNcell] = {0.0};
  Double_t ADCtvseP0[kNcell] = {0.0};
  Double_t ADCtvseP1[kNcell] = {0.0};
  Double_t ADCtvseP2[kNcell] = {0.0};

  //Set fitting parameters for gaussian fits
  Double_t asetpar[3];
  Double_t tsetpar[3];

  //No error on cell location
  Double_t cellerr[kNcell] = {0.};

  //Make arrays for tdc tgraphs
  Double_t tcell[kNcell] = {0.};
  Double_t tcval[kNcell] = {0.};
  Double_t tcvalw[kNcell] = {0.};
  Double_t tcerr[kNcell] = {0.};
  TH1D *tcellslice[kNcell];

  //Make arrays for adc tgraphs
  Double_t acell[kNcell] = {0.};
  Double_t acval[kNcell] = {0.};
  Double_t acvalw[kNcell] = {0.};
  Double_t acerr[kNcell] = {0.};
  TH1D *acellslice[kNcell];

  //Get averages
  Double_t tcval_avg = 0.;
  Int_t tcval_Ng = 0;
  Double_t acval_avg = 0.;
  Int_t acval_Ng = 0;
  
  //Fits for timewalk corrections
  for(Int_t c=0; c<kNcell; c++){
    //Fit the TDC vs E plots
    TF1 *fittdcTW = new TF1( "fittdcTW", TW_fit, 0, 300, 3 );
    fittdcTW->SetParameters(14,0.04,-77);
    fittdcTW->SetParLimits(0,2,26);
    fittdcTW->SetParLimits(1,0.01,0.2);
    fittdcTW->SetParLimits(2,-200,50);

    if( htdcVe[c]->GetEntries()>tFitMin ){
      htdcVe[c]->Fit("fittdcTW","Q","",5,300);
      TDCvseP0[c] = fittdcTW->GetParameter(0);
      TDCvseP1[c] = fittdcTW->GetParameter(1);
      TDCvseP2[c] = fittdcTW->GetParameter(2);
      htdcVe[c]->SetTitle(Form("P0:%f P1:%f P2:%f",TDCvseP0[c],TDCvseP1[c],TDCvseP2[c]));
    }
    htdcP0->Fill(TDCvseP0[c]);
    htdcP1->Fill(TDCvseP1[c]);
    htdcP2->Fill(TDCvseP2[c]);
  }

  //Separate loops to improve fits (avoid parameter recall)
  for(Int_t c=0; c<kNcell; c++){
    //Fit the ADCt vs E plots
    TF1 *fitadctTW = new TF1( "fitadctTW", TW_fit, 0, 300, 3 );
    fitadctTW->SetParameters(10,0.01,75);
    fitadctTW->SetParLimits(0,2,9);
    fitadctTW->SetParLimits(1,0.008,0.05);
    fitadctTW->SetParLimits(2,-50,200);

    if( hadctVe[c]->GetEntries()>tFitMin ){
      hadctVe[c]->Fit("fitadctTW","Q","",5,300);
      ADCtvseP0[c] = fitadctTW->GetParameter(0);
      ADCtvseP1[c] = fitadctTW->GetParameter(1);
      ADCtvseP2[c] = fitadctTW->GetParameter(2);
      hadctVe[c]->SetTitle(Form("P0:%f P1:%f P2:%f",ADCtvseP0[c],ADCtvseP1[c],ADCtvseP2[c]));
    }
    hadctP0->Fill(ADCtvseP0[c]);
    hadctP1->Fill(ADCtvseP1[c]);
    hadctP2->Fill(ADCtvseP2[c]);
  }

  //Fits for TDC time
  TCanvas *TDC_top = new TCanvas("TDC_top","TDC_top",1600,1200);
  TCanvas *TDC_bot = new TCanvas("TDC_bot","TDC_bot",1600,1200);

  TDC_top->Divide(12,12);
  TDC_bot->Divide(12,12);

  gStyle->SetOptStat(0);

  for(Int_t c=0; c<kNcell; c++){

    //Index through the canvas
    TDC_top->cd(c+1);
    if( c>=144 ){
      TDC_bot->cd(c-143);
      gStyle->SetOptStat(0);
    }

    //Get slices from htDiff_ID and fit for mean vals
    Double_t tfitl = 0.;
    Double_t tfith = 0.;
    tcell[c] = c;
    tcellslice[c] = htpDiff_ID->ProjectionY(Form("tcellslice_%d",c+1),c+1,c+1); //Trying htpDiff_ID from htDiff_ID
    tcval[c] = oldTDCoffsets[c]; //will overwrite if fit is good.
    tcvalw[c] = 0.; //leave as zero to evaluate fit values only

    Int_t sliceN = tcellslice[c]->GetEntries();
    if( sliceN<tFitMin ){
      tcellslice[c]->Draw();
      continue;
    }

    Double_t arimean = tcellslice[c]->GetMean();
    tsetpar[0] = sliceN;
    tsetpar[1] = arimean;
    tsetpar[2] = aTDCsig;
    tfitl = arimean - 4*aTDCsig;
    tfith = arimean + 4*aTDCsig;
    TF1 *gausfit = new TF1("gausfit",Gfit,tfitl,tfith,3);
    gausfit->SetLineWidth(4);
    gausfit->SetParameter(0,tsetpar[0]);
    gausfit->SetParameter(1,tsetpar[1]);
    gausfit->SetParLimits(1,tfitl,tfith);
    gausfit->SetParameter(2,tsetpar[2]);

    tcellslice[c]->Fit("gausfit","RBM");
    tcellslice[c]->Draw();

    tcval[c] = gausfit->GetParameter(1);
    tcerr[c] = gausfit->GetParameter(2);
    tcvalw[c] = gausfit->GetParameter(1);
    tcellslice[c]->SetTitle(Form("Mean:%f Sigma:%f",tcval[c],tcerr[c]));    

    tcval_avg += tcval[c];
    tcval_Ng++;

  }    
  TDC_top->Write();
  TDC_bot->Write();

  tcval_avg /= tcval_Ng;

  //Fits for ADC time 
  TCanvas *ADCt_top = new TCanvas("ADCt_top","ADCt_top",1600,1200);
  TCanvas *ADCt_bot = new TCanvas("ADCt_bot","ADCt_bot",1600,1200);

  ADCt_top->Divide(12,12);
  ADCt_bot->Divide(12,12);

  gStyle->SetOptStat(0);

  for(Int_t c=0; c<kNcell; c++){

    Int_t col = c%kNcols;

    //Index through the canvas
    ADCt_top->cd(c+1);
    if( c>=144 ){
      ADCt_bot->cd(c-143);
      gStyle->SetOptStat(0);
    }

    //Get slices from haDiff_ID and fit for mean vals
    Double_t afitl = 0.; //50. empirically
    Double_t afith = 0.; //75.
    acell[c] = c;
    acellslice[c] = haDiff_ID->ProjectionY(Form("acellslice_%d",c+1),c+1,c+1); //Trying hapDiff_ID from haDiff_ID
    acval[c] = oldADCtoffsets[c]; //will overwrite if fit is good.
    acvalw[c] = 0.; //will overwrite if fit is good.

    Int_t sliceN = acellslice[c]->GetEntries();
    if( sliceN<tFitMin ){
      acellslice[c]->Draw();
      continue;
    }
    Double_t arimean = acellslice[c]->GetMean();
    asetpar[0] = sliceN;
    asetpar[1] = arimean;
    asetpar[2] = aADCtsig;
    afitl = arimean - 5*aADCtsig;
    if( col==11 ) afitl=40.;
    afith = arimean + 4*aADCtsig;
    TF1 *gausfit = new TF1("gausfit",Gfit,afitl,afith,3);
    gausfit->SetLineWidth(4);
    gausfit->SetParameter(0,asetpar[0]);
    gausfit->SetParameter(1,asetpar[1]);
    gausfit->SetParLimits(1,afitl,afith);
    gausfit->SetParameter(2,asetpar[2]);

    acellslice[c]->Fit("gausfit","RBM");
    acellslice[c]->Draw();

    acval[c] = gausfit->GetParameter(1);
    acerr[c] = gausfit->GetParameter(2);
    acvalw[c] = gausfit->GetParameter(1);
    acellslice[c]->SetTitle(Form("Mean:%f Sigma:%f",acval[c],acerr[c]));    

    acval_avg += acval[c];
    acval_Ng++;

  }
  ADCt_top->Write();
  ADCt_bot->Write();

  acval_avg /= acval_Ng;

  
  //Make graphs with errors for reporting
  TGraphErrors *gtdc_c = new TGraphErrors( kNcell, tcell, tcvalw, cellerr, tcerr );
  gtdc_c->GetXaxis()->SetLimits(-10,290);  
  gtdc_c->GetYaxis()->SetLimits(tllim,tulim);
  gtdc_c->SetTitle("TDC_{hcal}-TDCmean_{hodo} vs Cell");
  gtdc_c->GetXaxis()->SetTitle("Cell");
  gtdc_c->GetYaxis()->SetTitle("TDC_{HCAL}-TDCMEAN_{HODO}");
  gtdc_c->SetMarkerStyle(20); // idx 20 Circles, idx 21 Boxes
  gtdc_c->Draw();
  gtdc_c->Write("gtdc_c");


  TGraphErrors *gadct_c = new TGraphErrors( kNcell, acell, acvalw, cellerr, acerr );
  gadct_c->GetXaxis()->SetLimits(-10,290);  
  gadct_c->GetYaxis()->SetLimits(allim,aulim);
  gadct_c->SetTitle("ADCt_{hcal}-TDCmean_{hodo} vs Cell");
  gadct_c->GetXaxis()->SetTitle("Cell");
  gadct_c->GetYaxis()->SetTitle("ADCt_{HCAL}-TDCMEAN_{HODO}");
  gadct_c->SetMarkerStyle(20); // idx 20 Circles, idx 21 Boxes
  gadct_c->Draw();
  gadct_c->Write("gadct_c");

  // Create parameter files
  ofstream tdcP0;
  ofstream tdcP1;
  ofstream tdcP2;
  ofstream adctP0;
  ofstream adctP1;
  ofstream adctP2;
  ofstream tdcoff;
  ofstream adctoff;

  //Write to outfiles if !qreplay
  if( !qreplay ){
    //Write first TDCvE parameter
    tdcP0.open( tdctwP0path );
    tdcP0 << "#HCal tdc vs E fit parameter P0 obtained " << date.c_str() << endl;
    tdcP0 << "#Exponential fit for i PMTs -> y = P0_i*exp(-P1_i*x)+P2_i, P0_i = " << endl;
    report << "#HCal tdc vs E fit parameter P0 obtained " << date.c_str() << endl;
    report << "#Exponential fit for i PMTs -> y = P0_i*exp(-P1_i*x)+P2_i, P0_i = " << endl;

    for( Int_t i=0; i<kNcell; i++ ){   
      tdcP0 << TDCvseP0[i] << endl;
      report << TDCvseP0[i] << endl;
    }

    tdcP0.close();
    report << endl << endl;

    //Write second TDCvE parameter
    tdcP1.open( tdctwP1path );
    tdcP1 << "#HCal tdc vs E fit parameter P1 obtained " << date.c_str() << endl;
    tdcP1 << "#Exponential fit for i PMTs -> y = P0_i*exp(-P1_i*x)+P2_i, P1_i = " << endl;
    report << "#HCal tdc vs E fit parameter P1 obtained " << date.c_str() << endl;
    report << "#Exponential fit for i PMTs -> y = P0_i*exp(-P1_i*x)+P2_i, P1_i = " << endl;

    for( Int_t i=0; i<kNcell; i++ ){   
      tdcP1 << TDCvseP1[i] << endl;
      report << TDCvseP1[i] << endl;
    }

    tdcP1.close();
    report << endl << endl;

    //Write third TDCvE parameter
    tdcP2.open( tdctwP2path );
    tdcP2 << "#HCal tdc vs E fit parameter P2 obtained " << date.c_str() << endl;
    tdcP2 << "#Exponential fit for i PMTs -> y = P0_i*exp(-P1_i*x)+P2_i, P2_i = " << endl;
    report << "#HCal tdc vs E fit parameter P2 obtained " << date.c_str() << endl;
    report << "#Exponential fit for i PMTs -> y = P0_i*exp(-P1_i*x)+P2_i, P2_i = " << endl;

    for( Int_t i=0; i<kNcell; i++ ){   
      tdcP2 << TDCvseP2[i] << endl;
      report << TDCvseP2[i] << endl;
    }

    tdcP2.close();
    report << endl << endl;

    //Write first ADCtvE parameter
    adctP0.open( adcttwP0path );
    adctP0 << "#HCal adct vs E fit parameter P0 obtained " << date.c_str() << endl;
    adctP0 << "#Exponential fit for i PMTs -> y = P0_i*exp(-P1_i*x)+P2_i, P0_i = " << endl;
    report << "#HCal adct vs E fit parameter P0 obtained " << date.c_str() << endl;
    report << "#Exponential fit for i PMTs -> y = P0_i*exp(-P1_i*x)+P2_i, P0_i = " << endl;

    for( Int_t i=0; i<kNcell; i++ ){   
      adctP0 << ADCtvseP0[i] << endl;
      report << ADCtvseP0[i] << endl;
    }

    adctP0.close();
    report << endl << endl;

    //Write second ADCtvE parameter
    adctP1.open( adcttwP1path );
    adctP1 << "#HCal adct vs E fit parameter P1 obtained " << date.c_str() << endl;
    adctP1 << "#Exponential fit for i PMTs -> y = P0_i*exp(-P1_i*x)+P2_i, P1_i = " << endl;
    report << "#HCal adct vs E fit parameter P1 obtained " << date.c_str() << endl;
    report << "#Exponential fit for i PMTs -> y = P0_i*exp(-P1_i*x)+P2_i, P1_i = " << endl;

    for( Int_t i=0; i<kNcell; i++ ){   
      adctP1 << ADCtvseP1[i] << endl;
      report << ADCtvseP1[i] << endl;
    }

    adctP1.close();
    report << endl << endl;

    //Write third ADCtvE parameter
    adctP2.open( adcttwP2path );
    adctP2 << "#HCal adct vs E fit parameter P2 obtained " << date.c_str() << endl;
    adctP2 << "#Exponential fit for i PMTs -> y = P0_i*exp(-P1_i*x)+P2_i, P2_i = " << endl;
    report << "#HCal adct vs E fit parameter P2 obtained " << date.c_str() << endl;
    report << "#Exponential fit for i PMTs -> y = P0_i*exp(-P1_i*x)+P2_i, P2_i = " << endl;

    for( Int_t i=0; i<kNcell; i++ ){   
      adctP2 << ADCtvseP2[i] << endl;
      report << ADCtvseP2[i] << endl;
    }

    adctP2.close();
    report << endl << endl;

    //Write TDC offsets
    tdcoff.open( tdcconstPath );
    tdcoff << "#HCal TDC offsets obtained " << date.c_str() << endl;
    tdcoff << "#Offsets obtained from fits over TDC distributions = " << endl;
    report << "#HCal TDC offsets obtained " << date.c_str() << endl;
    report << "#Offsets obtained from fits over TDC distributions = " << endl;
  
    for( Int_t i=0; i<kNcell; i++ ){ //check this!!
      //bool firstpeak = tcval[i]
      if( oldTDCoffsets[i]>0 && tcval[i]==0 ){
	tdcoff << tcval_avg << endl;
	report << tcval_avg << " with adjustment" << endl;
      }else{
	tdcoff << tcval[i]/TDCCalib + oldTDCoffsets[i] - TDC_target/TDCCalib << endl;
	report << tcval[i]/TDCCalib + oldTDCoffsets[i] - TDC_target/TDCCalib << endl;
      }
    }  

    tdcoff.close();

    report << endl << endl;

    //Write ADCt offsets
    adctoff.open( adctconstPath );
    adctoff << "#HCal ADCT offsets obtained " << date.c_str() << endl;
    adctoff << "#Offsets obtained from fits over ADCT distributions = " << endl;
    report << "#HCal ADCT offsets obtained " << date.c_str() << endl;
    report << "#Offsets obtained from fits over ADCT distributions = " << endl;
  
    for( Int_t i=0; i<kNcell; i++ ){   
      adctoff << acval[i] + oldADCtoffsets[i] - ADCt_target << endl;
      report << acval[i] + oldADCtoffsets[i] - ADCt_target << endl;
    }
  
    adctoff.close();
  }

  fout->Write();

  cout << endl << endl << "Elastic yield for analyzed runs: " << elasYield << ". Total event bins available for calibration = LH2 + LD2: " << TNEV_h << " + " << TNEV_d << " = " << TNEV_h + TNEV_d << ". Total number of events analyzed: " << lh2Events + ld2Events << "." << endl << endl;
  report << endl << endl << "Elastic yield for analyzed runs: " << elasYield << ". Total event bins available for calibration = LH2 + LD2: " << TNEV_h << " + " << TNEV_d << " = " << TNEV_h + TNEV_d << ". Total number of events analyzed: " << lh2Events + ld2Events << "." << endl << endl;
  
  if( qreplay ){
    cout << "Diagnostic iteration complete. Histograms written to file." << endl;
    report << "Diagnostic iteration complete. Histograms written to file." << endl << endl;
  }else{
    cout << "Timing analysis complete. Constants and histograms written to file." << endl;
    report << "Timing analysis complete. Constants and histograms written to file." << endl << endl;
  }

  st->Stop();

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;
  report << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

  report.close();

}
