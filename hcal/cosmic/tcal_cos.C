//SSeeds 9.12.22 Script to get TDC spectral data for use in alignment and timing resolution
#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <TChain.h>
#include <TSystem.h>
#include <TStopwatch.h>
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

const Int_t kNcell = 288; // Total number of HCal modules
const Int_t kNrows = 24; // Total number of HCal rows
const Int_t kNcols = 12; // Total number of HCal columns
const Int_t maxTdcChan = 10; // Set to accomodate original 5 TDCTrig channels with buffer
const Int_t maxTracks = 1000; // Reasonable limit on tracks to be stored per event
const Double_t hcalheight = 0.365; //m The height of the center of HCAL above beam

const Double_t PI = TMath::Pi();
const Double_t M_e = 0.00051;
const Double_t M_p = 0.938272;

const Double_t TDCCalib = 0.112;

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

// Inputs: run (run number to be analyzed), fresh (toggle to check previous ADC/TDC offsets or not, 1=no)
void tcal_cos( Int_t run, const char *outputfilename="tcal_cosmicOut.root", int fresh = 0, bool GMnOption = false ){
  
  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  // Get the date
  string date = getDate();
  
  // Declare Chain for many root files
  TChain *C = new TChain("T");

  // Declare general physics parameters to be modified by input config file
  Double_t fresh = 1; // Toggle for input fit params. 0: no input
  Double_t test = 1; // Keep track of DB read for analysis. 0: ../SBS_REPLAY/SBS_Replay/DB   1: ../seeds/SBS-replay/DB
  Double_t diag = 0; // Keep track of diagnostic canvas pdf option printed at end
  Double_t kine = 8; // Keep track of kinematic calibrating from
  Double_t spotcut = 0; 
  Double_t tFitMin = 30; // Minimum number of entries per channel to calibrate ADC/TDC time
  Double_t t_trig = 510; // Mean tdc trig value (HCAL - BB) 
  Double_t E_e = 1.92; // Energy of beam (incoming electrons from accelerator)
  Int_t TDCmiss = 0; // Keep track of TDC misses
  Double_t HCal_d = 14.5; // Distance to HCal from scattering chamber for comm1
  Double_t HCal_th = 35.0; // Angle that the center of HCal is at  
  Double_t W_mean = 0.93; // Mean of W at current kinematic
  Double_t W_sig = 0.039; // Width of W at current kinematic
  Double_t dx0 = 0.9; // Position of proton spot, x-x_expected
  Double_t dy0 = 0.62; // Position of proton spot, y-y_expected
  Double_t dx_sig = 0.09; // Max spread of proton spot, x-x_expected
  Double_t dy_sig = 0.15; // Max spread of proton spot, y-y_expected
  Long64_t elasYield = 0; // Keep track of total elastics analyzed

  // Declare arrays to hold old offset parameters
  Double_t oldADCtoffsets[kNcell]={0.};
  Double_t oldTDCoffsets[kNcell]={0.};
  Double_t oldP0[kNcell]={0.};
  Double_t oldP1[kNcell]={0.};
  Double_t oldP2[kNcell]={0.};

  cout << endl;

  if( GMnOption==true ){
    C->Add( Form("/lustre19/expphy/volatile/halla/sbs/seeds/rootfiles/hcal_general_%d*", run) );
  }else{
    C->Add( Form("/lustre19/expphy/volatile/halla/sbs/seeds/rootfiles/hcal_GEn_%d*", run) );
  }

  /*
  // Reading config file
  ifstream configfile(configfilename);
  TString currentline;
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      if(!currentline) cout << "WARNING: No file exists at " << currentline << "." << endl;
      C->Add(currentline);
    }    
  }
  TCut globalcut = "";
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endcut") ){
    if( !currentline.BeginsWith("#") ){
      globalcut += currentline;
    }    
  }
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("#") ){
    TObjArray *tokens = currentline.Tokenize(" ");
    Int_t ntokens = tokens->GetEntries();
    if( ntokens>1 ){
      TString skey = ( (TObjString*)(*tokens)[0] )->GetString();
      if( skey == "fresh" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	fresh = sval.Atof();
	cout << "Loading fresh setting: " << fresh << endl;
      }
      if( skey == "test" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	test = sval.Atof();
	cout << "Loading test setting: " << test << endl;
      }
      if( skey == "diag" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	diag = sval.Atof();
	cout << "Loading diagnostic setting: " << diag << endl;
      }
      if( skey == "kine" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	kine = sval.Atof();
	cout << "Loading kinematic setting: " << kine << endl;
      }
      if( skey == "tFitMin" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	tFitMin = sval.Atof();
	cout << "Loading timing fit min entries: " << tFitMin << endl;
      }
      if( skey == "t_trig" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	t_trig = sval.Atof();
	cout << "Loading mean timing difference BB/HCal trigger: " << t_trig << endl;
      }
      if( skey == "E_e" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	E_e = sval.Atof();
	cout << "Loading beam energy: " << E_e << endl;
      }
      if( skey == "HCal_d" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	HCal_d = sval.Atof();
	cout << "Loading HCal distance: " << HCal_d << endl;
      }
      if( skey == "HCal_th" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	HCal_th = sval.Atof();	
	cout << "Loading HCal angle: " << HCal_th << endl;
      }
      if( skey == "W_mean" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        W_mean = sval.Atof();
	cout << "Loading W mean cut: " << W_mean << endl;
      }
      if( skey == "W_sig" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        W_sig = sval.Atof();
	cout << "Loading W sigma cut: " << W_sig << endl;
      }
      if( skey == "dx0" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	dx0 = sval.Atof();
	cout << "Loading x position of spot: " << dx0 << endl;
      }
      if( skey == "dy0" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	dy0 = sval.Atof();
	cout << "Loading y position of spot: " << dy0 << endl;
      }
      if( skey == "dx_sig" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	dx_sig = sval.Atof();
	cout << "Loading x sigma of spot: " << dx_sig << endl;
      }
      if( skey == "dy_sig" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	dy_sig = sval.Atof();
	cout << "Loading y sigma of spot: " << dy_sig << endl;
      }
    }
    delete tokens;
  }
  */

  // Path for fit parameters
  string P0path = "params/P0.txt";
  string P1path = "params/P1.txt";
  string P2path = "params/P2.txt";
  string offpath = "params/offsets.txt";

  /*
  // Path for previous ADCt and TDC constants
  string inConstPath;
  if( test==0 ){
    inConstPath = "/w/halla-scshelf2102/sbs/SBS_REPLAY/SBS-replay/DB/db_sbs.hcal.dat";
  }else if( test==1 ){
    inConstPath = "/w/halla-scshelf2102/sbs/seeds/SBS-replay/DB/db_sbs.hcal.dat";
  }else{
    cout << "ERROR: No database file where specified" << endl;
    return 0;
  }

  // Reading ADC and TDC timing offsets from database
  cout << endl << endl << "Loading previous offsets from database file: " << inConstPath << ".." << endl;
  ifstream inConstFile( inConstPath );
  if( !inConstFile ){
    cerr << endl << "ERROR: No input constant file present -> path to db_sbs.hcal.dat expected." << endl;
    return 0;
  }

  Int_t n0=0;
  Double_t d0=0;
  string line0;
  bool skip_line = true;
  bool skip_one_line = true;
  bool pass_first_cond = false;
  bool pass_all_cond = false;

  while( getline( inConstFile, line0 ) ){

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
  */

  // Reading fit parameters if fresh=0
  if( fresh==0 ){
    cout << "Loading fit parameters.." << endl;

    //Get first parameter (all channels)
    ifstream P0file( P0path );
    if( !P0file ){
      cerr << endl << "ERROR: No input constant file present -> path to P0.txt expected." << endl;
      return 0;
    }
  
    Int_t n1=0;
    Double_t d1;
    string line1;
    
    while( getline( P0file, line1 ) ){
      if( line1.at(0) == '#' ) {
	continue;
      }
      stringstream ss( line1 );
      ss >> d1;     
      oldP0[n1] = d1;      
      n1++;
    }

    //Get second parameter (all channels)
    ifstream P1file( P1path );
    if( !P1file ){
      cerr << endl << "ERROR: No input constant file present -> path to P1.txt expected." << endl;
      return 0;
    }

    n1=0;
    d1=0;
    string line2;
    
    while( getline( P1file, line2 ) ){
      if( line2.at(0) == '#' ) {
	continue;
      }
      stringstream ss( line2 );
      ss >> d1;    
      oldP1[n1] = d1; 
      n1++;
    }

    //Get third parameter (all channels)
    ifstream P2file( P2path );
    if( !P2file ){
      cerr << endl << "ERROR: No input constant file present -> path to P2.txt expected." << endl;
      return 0;
    }

    n1=0;
    d1=0;
    string line3;
    
    while( getline( P2file, line3 ) ){
      if( line3.at(0) == '#' ) {
	continue;
      }
      stringstream ss( line3 );
      ss >> d1;     
      oldP2[n1] = d1;     
      n1++;
    }

    //Get offsets (all channels)
    ifstream offfile( offpath );
    if( !offfile ){
      cerr << endl << "ERROR: No input constant file present -> path to P2.txt expected." << endl;
      return 0;
    }

    n1=0;
    d1=0;
    string line4;
    
    while( getline( P2file, line4 ) ){
      if( line4.at(0) == '#' ) {
	continue;
      }
      stringstream ss( line4 );
      ss >> d1;     
      oldTDCoffsets[n1] = d1;     
      n1++;
    }

    //Print old fit parameters to console
    cout << endl << endl << "Timewalk Fit P0 vs Channels: " << endl;
    
    for( Int_t r=0; r<kNrows; r++){
      for( Int_t c=0; c<kNcols; c++){
	Int_t i = r*kNcols+c;
	cout << oldP0[i] << " ";
      }
      cout << endl;
    }
    
    cout << endl << endl << "Timewalk Fit P1 vs Channels: " << endl;
    
    for( Int_t r=0; r<kNrows; r++){
      for( Int_t c=0; c<kNcols; c++){
	Int_t i = r*kNcols+c;
	cout << oldP1[i] << " ";
      }
      cout << endl;
    }

   cout << endl << endl << "Timewalk Fit P2 vs Channels: " << endl;
    
    for( Int_t r=0; r<kNrows; r++){
      for( Int_t c=0; c<kNcols; c++){
	Int_t i = r*kNcols+c;
	cout << oldP2[i] << " ";
      }
      cout << endl;
    }

    cout << endl << endl << "Old TDC Offsets vs Channels: " << endl;
    
    for( Int_t r=0; r<kNrows; r++){
      for( Int_t c=0; c<kNcols; c++){
	Int_t i = r*kNcols+c;
	cout << oldTDCoffsets[i] << " ";
      }
      cout << endl;
    }
  }

  cout << endl << endl << "Fit parameters loaded." << endl;
  
  TEventList *elist = new TEventList("elist","Elastic Event List");
  C->Draw(">>elist",globalcut);

  cout << endl << "Event list populated with cut placed on elastics." << endl;

  // Declare general detector and physics parameters
  Double_t TDCT_id[maxTdcChan], TDCT_tdc[maxTdcChan]; 
  Int_t TDCTndata;
  UInt_t runI=1; 
  UInt_t runN=1;
  UInt_t TBits;
  ULong64_t runT;
  Double_t kineW2;

  // BB params
  Double_t BBtr_p[maxTracks], BBtr_px[maxTracks], BBtr_py[maxTracks], BBtr_pz[maxTracks];
  Double_t BBtr_vz[maxTracks], BBtr_chi2[maxTracks];
  Double_t BBtr_n, BBps_x, BBps_y, BBps_e, BBsh_x, BBsh_y, BBsh_e;

  // HCal params
  Double_t HCALx, HCALy, HCALe;
  Double_t HCALblktdc[kNcell], HCALblka[kNcell];
  Double_t HCALtdc[kNcell], HCALatime[kNcell], HCALa_p[kNcell];
  Double_t crow, ccol, nblk;
  Double_t cblkid[kNcell], cblke[kNcell];
  Double_t ce[kNcell], catime[kNcell], cid[kNcell], cnblk[kNcell], ctdc[kNcell];
  Double_t nclus;

  // Hodo params
  Double_t HODOtmean;
  Double_t HODOnclus;
  
  // Declare root tree variables and set values to memory locations in root file
  // Turn off all branches
  C->SetBranchStatus( "*", 0 );
  // Turn on specific branches for analysis
  C->SetBranchStatus( "sbs.hcal.x", 1 );
  C->SetBranchStatus( "sbs.hcal.y", 1 );
  C->SetBranchStatus( "sbs.hcal.e", 1 );
  C->SetBranchStatus( "sbs.hcal.rowblk", 1 );
  C->SetBranchStatus( "sbs.hcal.colblk", 1 );
  C->SetBranchStatus( "sbs.hcal.nblk", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.tdctime", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.atime", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.id", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.e", 1 );
  C->SetBranchStatus( "sbs.hcal.clus.e", 1 );
  C->SetBranchStatus( "sbs.hcal.clus.atime", 1 );
  C->SetBranchStatus( "sbs.hcal.clus.id", 1 );
  C->SetBranchStatus( "sbs.hcal.clus.nblk", 1 );
  C->SetBranchStatus( "sbs.hcal.clus.tdctime", 1 );
  C->SetBranchStatus( "sbs.hcal.nclus", 1 );  
  C->SetBranchStatus( "bb.tr.chi2", 1 );
  C->SetBranchStatus( "bb.tr.n", 1 );
  C->SetBranchStatus( "bb.tr.px", 1 );
  C->SetBranchStatus( "bb.tr.py", 1 );
  C->SetBranchStatus( "bb.tr.pz", 1 );    
  C->SetBranchStatus( "bb.tr.vz", 1 );
  C->SetBranchStatus( "bb.tr.p", 1 );
  C->SetBranchStatus( "bb.ps.e", 1 );
  C->SetBranchStatus( "bb.ps.x", 1 );
  C->SetBranchStatus( "bb.ps.y", 1 );
  C->SetBranchStatus( "bb.sh.e", 1 );
  C->SetBranchStatus( "bb.sh.x", 1 );
  C->SetBranchStatus( "bb.sh.y", 1 );
  C->SetBranchStatus( "bb.tdctrig.tdc", 1 );
  C->SetBranchStatus( "bb.tdctrig.tdcelemID", 1 );
  C->SetBranchStatus( "Ndata.bb.tdctrig.tdcelemID", 1 );
  C->SetBranchStatus( "bb.hodotdc.clus.tmean", 1 );
  C->SetBranchStatus( "bb.hodotdc.nclus", 1 );
  C->SetBranchStatus( "fEvtHdr.fRun", 1 );
  C->SetBranchStatus( "fEvtHdr.fEvtTime", 1 );
  C->SetBranchStatus( "fEvtHdr.fEvtNum", 1 );
  C->SetBranchStatus( "fEvtHdr.fTrigBits", 1 );
  C->SetBranchStatus( "e.kine.W2", 1 );

  // Map the branches to variables
  C->SetBranchAddress( "sbs.hcal.x", &HCALx );
  C->SetBranchAddress( "sbs.hcal.y", &HCALy );
  C->SetBranchAddress( "sbs.hcal.e", &HCALe );
  C->SetBranchAddress( "sbs.hcal.rowblk", &crow );
  C->SetBranchAddress( "sbs.hcal.colblk", &ccol );
  C->SetBranchAddress( "sbs.hcal.nblk", &nblk );
  C->SetBranchAddress( "sbs.hcal.clus_blk.id", cblkid );
  C->SetBranchAddress( "sbs.hcal.clus_blk.e", cblke );
  C->SetBranchAddress( "sbs.hcal.clus_blk.tdctime", HCALblktdc );
  C->SetBranchAddress( "sbs.hcal.clus_blk.atime", HCALblka );
  C->SetBranchAddress( "sbs.hcal.clus.e", ce );
  C->SetBranchAddress( "sbs.hcal.clus.atime", catime );
  C->SetBranchAddress( "sbs.hcal.clus.id", cid );
  C->SetBranchAddress( "sbs.hcal.clus.nblk", cnblk );
  C->SetBranchAddress( "sbs.hcal.clus.tdctime", ctdc );
  C->SetBranchAddress( "sbs.hcal.nclus", &nclus );
  C->SetBranchAddress( "bb.tr.chi2", BBtr_chi2 );
  C->SetBranchAddress( "bb.tr.px", BBtr_px );
  C->SetBranchAddress( "bb.tr.py", BBtr_py );
  C->SetBranchAddress( "bb.tr.pz", BBtr_pz );
  C->SetBranchAddress( "bb.tr.vz", BBtr_vz );
  C->SetBranchAddress( "bb.tr.p", BBtr_p );
  C->SetBranchAddress( "bb.ps.e", &BBps_e );
  C->SetBranchAddress( "bb.ps.x", &BBps_x );
  C->SetBranchAddress( "bb.ps.y", &BBps_y );
  C->SetBranchAddress( "bb.sh.e", &BBsh_e );
  C->SetBranchAddress( "bb.sh.x", &BBsh_x );
  C->SetBranchAddress( "bb.sh.y", &BBsh_y );
  C->SetBranchAddress( "bb.tdctrig.tdcelemID", TDCT_id );
  C->SetBranchAddress( "bb.tdctrig.tdc", TDCT_tdc );
  C->SetBranchAddress( "Ndata.bb.tdctrig.tdcelemID", &TDCTndata );  
  C->SetBranchAddress( "bb.hodotdc.clus.tmean", &HODOtmean );
  C->SetBranchAddress( "bb.hodotdc.nclus", &HODOnclus );
  C->SetBranchAddress( "fEvtHdr.fRun", &runI );
  C->SetBranchAddress( "fEvtHdr.fEvtTime", &runT );
  C->SetBranchAddress( "fEvtHdr.fEvtNum", &runN );
  C->SetBranchAddress( "fEvtHdr.fTrigBits", &TBits );
  C->SetBranchAddress( "e.kine.W2", &kineW2 );

  cout << "Tree variables linked." << endl;

  // Define a clock to check macro processing time
  TStopwatch *sw = new TStopwatch();
  sw->Start();
  
  // Declare outfile
  TFile *fout = new TFile( outputfilename, "RECREATE" );
  
  // Initialize vectors and arrays
  Double_t TDCoffsets[kNcell] = {0.0};
  Double_t ADCtoffsets[kNcell] = {0.0};
  Double_t TDCsig[kNcell] = {0.0};
  Double_t ADCtsig[kNcell] = {0.0};
  Double_t TvseP0[kNcell] = {0.0};
  Double_t TvseP1[kNcell] = {0.0};
  Double_t TvseP2[kNcell] = {0.0};
  Double_t TvseCorr[kNcell] = {0.0};

  // Initialize histograms
  TH1D *htdcAll = new TH1D("htdcAll","Aggregate TDC (All Channels); ns",250,-200,50);
  TH1D *ht_HCAL_HODO_All = new TH1D("ht_HCAL_HODO_All","TDC (HCal - Hodo, All Channels); ns",250,-200,50);  
  TH2D *ht_HCAL_HODO = new TH2D("ht_HCAL_HODO",";TDC_{HCAL} (ns);TDC_{HODO} (ns)", 150, -150, 0, 150, -15, 15);  
  TH2D *ht_HCAL_HODO_corr = new TH2D("ht_HCAL_HODO_corr",";TDC_{HCAL} (ns);TDC_{HODO} (ns)", 150, -150, 0, 150, -15, 15);  
  TH2D *htDiff_ID = new TH2D("htDiff_ID",";Channel;TDC_{HCAL}-TDC_{HODO} (ns)",288,0,288,300,-150,0);
  TH2D *haDiff_ID = new TH2D("haDiff_ID",";Channel;ADCtime_{HCAL}-TDC_{HODO} (ns)",288,0,288,300,0,150);
  TH1D *hdx = new TH1D( "hdx", "HCal X - X Expected", 400, -2, 2 );
  TH1D *hdy = new TH1D( "hdy", "HCal Y - Y Expected", 400, -2, 2 );
  TH2D *hdxdy_HCAL = new TH2D("hdxdy_HCAL",";y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)", 250, -5.0, 5.0, 250, -10, 10 );
  TH1D *hHCALe = new TH1D( "hHCALe","E HCal Cluster E", 500, 0., 500 );

  // Physics and trigger histograms
  //TH1D *el_ev = new TH1D( "el_ev",";Event", 50000, 0, 50000);
  TH1D *hDiff = new TH1D( "hDiff","HCal time - BBCal time (ns)", 1300, -500, 800 );
  TH1D *hW = new TH1D( "W", "W; GeV", 250, 0.3, 1.5 );
  TH1D *hNBlk = new TH1D( "hNBlk", "Number of Blocks in Primary Cluster; Number", 25, 0, 25 );
  TH1D *hNClus = new TH1D( "hNClus", "Number of Clusters; Number", 25, 0, 25 );
  TH1D *hW_cut = new TH1D( "W_cut", "W Elastic Cut; GeV", 250, 0.3, 1.5 );
  TH1D *hQ2 = new TH1D( "Q2", "Q2; GeV", 250, 0.5, 3.0 );
  TH1D *hE_ep = new TH1D( "hE_ep","E_ep; GeV", 500, 0.0, E_e*1.5 ); 
  TH1D *hE_pp = new TH1D( "hE_pp", "E_pp; GeV", 500, 0.0, E_e*1.5 );
  TH1D *hKE_p = new TH1D( "hKE_p", "KE_pp; GeV", 500, 0.0, E_e*1.5 );
  TH1D *hTBBt = new TH1D( "hTBBt","BBCal Trigger (L1A) (ns)", 500, 340, 390 ); //trigger for experiment
  TH1D *hTRF = new TH1D( "hTRF","RF Signature (ns)", 1700, -10, 160 );
  TH1D *hTRFmod = new TH1D( "hTRFmod","RF Signature Modulo 4ns (ns)", 200, -10, 10 );
  TH1D *hTEDTM = new TH1D( "hTEDTM","Electronic Dead Time Monitor Signal (ns)", 2000, -1, 2 );
  TH1D *hTBBlo = new TH1D( "hTBBlo","BBCal Lo Trigger (ns)", 199, 1, 200 );
  TH1D *hTBBhiv = new TH1D( "hTBBhiv","BBCal Lo (ns)", 2010, -10, 2000 );
  TH1D *hTHCAL = new TH1D( "hTHCAL","HCal time (ns)", 1000, 800, 900 );
  TH1D *hTHODO = new TH1D( "hTHODO","Hodo mean TDC time (ns)", 2000, -1000, 1000 );
  TH1D *hTHCALvRF = new TH1D( "hTHCALvRF","HCal time - RF Signature Modulo 4ns (ns)", 1000, 800, 900 );
  TH1D *hHODOnclus = new TH1D( "hHODOnclus","Number of Hodoscope Clusters", 50, 0., 50. );

  // Histograms by cell
  TH1D *htdc[kNcell];
  TH1D *htdcDiff[kNcell];
  TH1D *hcorr[kNcell];
  TH1D *htdc_corr[kNcell];
  TH1D *htdcDiff_corr[kNcell];
  TH1D *he[kNcell];
  TH1D *hatime[kNcell];
  TH2D *htdcVe[kNcell];
  
  for( Int_t i=0; i<kNcell; i++ ){
    htdc[i] = new TH1D(Form("htdc_bl%d",i),";TDC_{HCAL} (ns)",250,-200,50);
    htdcDiff[i] = new TH1D(Form("htdcDiff_bl%d",i),";TDC_{HCAL}-TDC_{HODO} (ns)",250,-200,50);
    hcorr[i] = new TH1D(Form("hcorr_bl%d",i),"TDC Timewalk Correction; ns",200,0,20);
    htdc_corr[i] = new TH1D(Form("htdc_corr_bl%d",i),";TDC_{HCAL} (ns)",250,-200,50);
    htdcDiff_corr[i] = new TH1D(Form("htdcDiff_corr_bl%d",i),";TDC_{HCAL}-TDC_{HODO} (ns)",250,-200,50);
    he[i] = new TH1D(Form("he_bl%d",i),"E; GeV",1000,0,100);
    hatime[i] = new TH1D(Form("hatime_bl%d",i),";ADCt_{HCAL} (ns)",180,0,180);
    htdcVe[i] = new TH2D(Form("htdcVe_bl%d",i),Form(";E_{bl%d} (GeV);TDC_{HCAL} (ns)",i),1000,0.0,500,250,-200,50);
  }

  cout << "Variables and histograms defined." << endl;

  // Set long int to keep track of total entries
  Long64_t Nevents = elist->GetN();
  UInt_t run_number = 0;

  cout << endl << "All parameters loaded and initialization complete." << endl << endl;
  cout << "Opened up TChain with nentries: " << C->GetEntries() << "." << endl << endl;

  //Loop over events
  cout << "Main loop over all data commencing.." << endl;
  Double_t progress = 0.;
  Double_t timekeeper = 0., timeremains = 0.;
  while(progress<1.0){
    Int_t barwidth = 70;
    for(Long64_t nevent = 0; nevent<Nevents; nevent++){
      
      //Create a progress bar
      cout << "[";
      Int_t pos = barwidth * progress;
      for(Int_t i=0; i<barwidth; ++i){
	if(i<pos) cout << "=";
	else if(i==pos) cout << ">";
	else cout << " ";
      }
      
      //Calculate remaining time 
      sw->Stop();
      timekeeper += sw->RealTime();
      if( nevent%100000 == 0 && nevent!=0 ) 
	timeremains = timekeeper*( Double_t(Nevents)/Double_t(nevent) - 1. ); 
      sw->Reset();
      sw->Continue();
      
      progress = (double)((nevent+1.)/Nevents);
      cout << "] " << int(progress*100.) << "%, elastic events: " << elasYield << "/" << nevent << ", run: " << run_number << ", time remaining: " << int(timeremains/60.) << "m \r";
      cout.flush();

      if( nevent == Nevents ){
	cout << "Hit limit on event " << nevent << endl;
	break;
      }
      
      C->GetEntry( elist->GetEntry( nevent ) ); 

      if( run_number!=runI ){
	run_number=runI;
	//cout << "Now analyzing run " << run_number << "." << endl;
      }
      
      //BEGIN coincidence timing cut
      //Cut on BBCal and HCal trigger coincidence and plot all other trigger timing. All ref to tdcs need index adjustment (ihit+1).
      Double_t bbcal_time=0., hcal_time=0., RF_time=0.; 
      Double_t bbcalLO_time=0., bbcalHIveto_time=0., edtm_time=0.;
      for(Int_t ihit=0; ihit<TDCTndata; ihit++){
	//cout << "TDCT_id[ihit]:" << TDCT_id[ihit] << ", TDCT_tdc:" << TDCT_tdc[ihit] << endl;
	if(TDCT_id[ihit]==5) bbcal_time=TDCT_tdc[ihit];
	if(TDCT_id[ihit]==4) RF_time=TDCT_tdc[ihit];
	if(TDCT_id[ihit]==3) edtm_time=TDCT_tdc[ihit];
	if(TDCT_id[ihit]==2) bbcalLO_time=TDCT_tdc[ihit];
	if(TDCT_id[ihit]==1) bbcalHIveto_time=TDCT_tdc[ihit];
	if(TDCT_id[ihit]==0) hcal_time=TDCT_tdc[ihit];
      }
      Double_t diff = hcal_time - bbcal_time; 
      Double_t RFmod = std::fmod(RF_time,4); //RF Signature measures the bunch crossing of beam electrons at 4ns intervals
      hTBBt->Fill( bbcal_time );
      hTRF->Fill( RF_time );
      hTRF->Fill( RFmod );
      hTEDTM->Fill( edtm_time );
      hTBBlo->Fill( bbcalLO_time );
      hTBBhiv->Fill( bbcalHIveto_time );
      hTHCAL->Fill( hcal_time );
      hTHCALvRF->Fill( hcal_time - RFmod );
      hDiff->Fill( diff );

      if( fabs( diff-t_trig )>30 ) continue;
      //cout << "hcal_time:" << hcal_time << ", BBCal time:" << bbcal_time << endl;

      //END coincidence timing cut

      //BEGIN elastic projections/cuts
      Double_t etheta = acos( BBtr_pz[0]/BBtr_p[0] );
      Double_t ephi = atan2( BBtr_py[0], BBtr_px[0] );
      TVector3 vertex(0,0,BBtr_vz[0]); // z location of vertex in hall coordinates
      TLorentzVector Pbeam(0,0,E_e,E_e); //Mass of e negligable
      TLorentzVector kprime(BBtr_px[0],BBtr_py[0],BBtr_pz[0],BBtr_p[0]);
      TLorentzVector Ptarg(0,0,0,M_p);
      TLorentzVector q = Pbeam - kprime;
      TLorentzVector PgammaN = Ptarg + q; //(-px, -py, ebeam - pz, Mp + ebeam - p)
      Double_t pel = E_e/(1.+E_e/M_p*(1.-cos(etheta)));
      Double_t nu = E_e - BBtr_p[0];
      Double_t pp = sqrt(pow(nu,2)+2.*M_p*nu); //momentum of the proton
      Double_t phinucleon = ephi + PI; //assume coplanarity
      Double_t thetanucleon = acos( (E_e - BBtr_p[0]*cos(etheta))/pp ); //use elastic constraint on nucleon kinematics
      TVector3 pNhat( sin(thetanucleon)*cos(phinucleon),sin(thetanucleon)*sin(phinucleon),cos(thetanucleon));

      //Define HCal coordinate system
      TVector3 HCAL_zaxis(sin(-HCal_th),0,cos(-HCal_th));
      TVector3 HCAL_xaxis(0,-1,0);
      TVector3 HCAL_yaxis = HCAL_zaxis.Cross(HCAL_xaxis).Unit();	
      TVector3 HCAL_origin = HCal_d * HCAL_zaxis + hcalheight * HCAL_xaxis;

      //Define intersection points for hadron vector
      Double_t sintersect = ( HCAL_origin - vertex ).Dot( HCAL_zaxis ) / ( pNhat.Dot( HCAL_zaxis ) );
      TVector3 HCAL_intersect = vertex + sintersect * pNhat;
      Double_t yexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_yaxis );
      Double_t xexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_xaxis );
      Double_t E_ep = sqrt( pow(M_e,2) + pow(BBtr_p[0],2) ); // Obtain the scattered electron energy
      Double_t p_ep = BBtr_p[0];
      Double_t Q2 = 2*E_e*E_ep*( 1-(BBtr_pz[0]/p_ep) ); // Obtain Q2 from beam energy, outgoing electron energy, and momenta
      //Double_t W = PgammaN.M();
      Double_t W = sqrt( kineW2 );
      //Use the electron kinematics to predict the proton momentum assuming elastic scattering on free proton at rest:
      Double_t E_pp = nu+M_p; // Get total energy of the proton
      //Double_t Enucleon = sqrt(pow(pp,2)+pow(M_p,2)); // Check on E_pp, same
      Double_t KE_p = nu; //For elastics
      hKE_p->Fill( KE_p );   
      hE_ep->Fill( E_ep );	
      hE_pp->Fill( E_pp );	
      hW->Fill( W );
      hQ2->Fill( p_ep );	

      //Fill some HCal histos
      hHCALe->Fill( HCALe*1000 ); //GeV to MeV
      hHODOnclus->Fill( HODOnclus );

      Double_t xDiff = HCALx - xexpect_HCAL;
      Double_t yDiff = HCALy - yexpect_HCAL;

      //Check how often the tdc failed to register for an otherwise good event
      if( HCALe>0.02 && HCALblktdc[0]<-400 && HODOnclus<10 ) TDCmiss++;
	
      //Tight cut: W elastic peak, BBCal/HCal timing coincidence, and expected position of scattered nucleon in HCal
      bool tightelastic = fabs( W-W_mean )<W_sig && fabs( diff+t_trig )<30 && (fabs(xDiff-dx0)<dx_sig || fabs(yDiff-dy0)<dy_sig); 

      //Loose cut: W elastic peak, BBCal/HCal timing coincidence
      bool elastic = fabs( W-W_mean )<W_sig && fabs( diff-t_trig )<30;

      if( elastic ){
	hdxdy_HCAL->Fill( yDiff, xDiff );
	hdy->Fill( yDiff );
	hdx->Fill( xDiff );
      }

      if( elastic==false ) continue;
      //if( tightelastic==false ) continue;
      elasYield++;

      hW_cut->Fill( W );
      hNBlk->Fill( nblk );
      hNClus->Fill( nclus );
      //el_ev->Fill( nevent );

      //Get timing offsets for tdc and adc from primary block in cluster

      //Fill histograms. Note that this assumes that atime and tdc are fixed size arrays!
      for( Int_t i=0; i<nblk; i++ ){

	Int_t ID = (Int_t)cblkid[i]-1;
	Double_t blkE = cblke[i]*1000; //Convert to MeV at this setting

	//cout << i << " " << ID << endl;

	Double_t TDCcorr = 0;
	//if( fresh==0 ) TDCcorr = oldP0[i]*exp(-oldP1[i]*HCALa_p[i])-oldP2[i];
	if( fresh==0 ) TDCcorr = oldP0[ID]*exp(-oldP1[ID]*cblke[i]);
	
	hcorr[i]->Fill( TDCcorr );

	if( HCALblka[i]>0 ){
	  hatime[ID]->Fill( HCALblka[i] );
	  haDiff_ID->Fill( ID, HCALblka[i]-HODOtmean );
	  he[ID]->Fill( blkE );
	}
	if( HCALblktdc[i]<10000 ){
	  htdc[ID]->Fill( HCALblktdc[i] );
	  htdcDiff[ID]->Fill( HCALblktdc[i]-HODOtmean );
	  htdc_corr[ID]->Fill( HCALblktdc[i]-TDCcorr );
	  htdcDiff_corr[ID]->Fill( HCALblktdc[i]-HODOtmean-TDCcorr );
	  htdcVe[ID]->Fill( blkE, HCALblktdc[i] );
	  htdcAll->Fill( HCALblktdc[i] );
	  ht_HCAL_HODO_All->Fill( HCALblktdc[i]-HODOtmean );
	  ht_HCAL_HODO->Fill( HCALblktdc[i], HODOtmean );
	  ht_HCAL_HODO_corr->Fill( HCALblktdc[i]-TDCcorr, HODOtmean );
	  htDiff_ID->Fill( ID, HCALblktdc[i]-HODOtmean );

	}
      }
    }
  }
  
  cout << endl << "Fitting TDC/ADCt distributions";

  //Fit timing data and extract offsets

  Double_t posErr[kNcell] = {0.};
  Double_t X[kNcell];
  Double_t Xval[kNcell];
  Double_t Xerr[kNcell];
  Double_t TDCerr[kNcell];
  Double_t ADCterr[kNcell];

  for(Int_t i=0; i<kNcell; i++){

    X[i] = i;

    if( i==244 || i==264 ) continue;

    Int_t r = (i)/kNcols;
    Int_t c = (i)%kNcols;

    TF1 *f1;
    TF1 *f2;
    TF1 *f3;
    TF1 *fitTW = new TF1( "fitTW", TW_fit, 0, 50, 3 );
    fitTW->SetParameters(20,0.2,-80);
    fitTW->SetParLimits(0,0,30);
    fitTW->SetParLimits(1,0,1);

    //cout << i << " " << htdcDiff_corr[i]->GetEntries() << endl;

    if( htdcDiff_corr[i]->GetEntries()>tFitMin ){
      htdcDiff_corr[i]->Fit("gaus","Q","",-200,50);
      f1=htdcDiff_corr[i]->GetFunction("gaus");
      TDCoffsets[i] = f1->GetParameter(1);
      TDCsig[i] = f1->GetParameter(2);
      TDCerr[i] = TDCsig[i]/sqrt(htdc_corr[i]->GetEntries());
      //htdcDiff_corr[i]->SetName(Form("htdcDiff_corr_bl%d_r%d_c%d",i,r,c));
      htdcDiff_corr[i]->SetTitle(Form("TDC Diff Corr, Mean:%f Sig:%f",f1->GetParameter(1),f1->GetParameter(2)));
    }

    if( hatime[i]->GetEntries()>tFitMin ){
      hatime[i]->Fit("gaus","Q","",0,150);
      f2=hatime[i]->GetFunction("gaus");
      ADCtoffsets[i] = f2->GetParameter(1);
      ADCtsig[i] = f2->GetParameter(2);
      ADCterr[i] = ADCtsig[i]/sqrt(hatime[i]->GetEntries());
      hatime[i]->SetName(Form("hatime_bl%d_r%d_c%d",i,r,c));
      hatime[i]->SetTitle(Form("ADCt, Mean:%f Sig:%f",f2->GetParameter(1),f2->GetParameter(2)));
    }
    
    if( htdc[i]->GetEntries()>tFitMin ){
      htdc[i]->Fit("gaus","Q","",-200,50);
      f3=htdc[i]->GetFunction("gaus");
      TDCoffsets[i] = f3->GetParameter(1);
      TDCsig[i] = f3->GetParameter(2);
      TDCerr[i] = TDCsig[i]/sqrt(htdc[i]->GetEntries());
      htdc[i]->SetName(Form("htdc_bl%d_r%d_c%d",i,r,c));
      htdc[i]->SetTitle(Form("TDC, Mean:%f Sig:%f",f3->GetParameter(1),f3->GetParameter(2)));
    }

    if( htdcVe[i]->GetEntries()>tFitMin ){
      htdcVe[i]->Fit("fitTW","Q","",5,300);
      TvseP0[i] = fitTW->GetParameter(0);
      TvseP1[i] = fitTW->GetParameter(1);
      TvseP2[i] = fitTW->GetParameter(2);
      htdcVe[i]->SetTitle(Form("P0:%f P1:%f P2:%f",TvseP0[i],TvseP1[i],TvseP2[i]));

      //cout << ".";
    }
  }

  TGraphErrors *cTDC = new TGraphErrors( kNcell, X, TDCoffsets, posErr, TDCsig );
  cTDC->GetXaxis()->SetLimits(0,kNcell);  
  cTDC->GetYaxis()->SetLimits(-200,50);
  cTDC->SetTitle("TDC Time");
  cTDC->GetXaxis()->SetTitle("Channel");
  cTDC->GetYaxis()->SetTitle("TDC Time (ns)");
  cTDC->SetMarkerStyle(20); // idx 20 Circles, idx 21 Boxes
  cTDC->Write("cTDC_off");

  
  TGraphErrors *cADCt = new TGraphErrors( kNcell, X, ADCtoffsets, posErr, ADCtsig );
  cADCt->GetXaxis()->SetLimits(0,kNcell);  
  cADCt->GetYaxis()->SetLimits(0,150);
  cADCt->SetTitle("ADC Time");
  cADCt->GetXaxis()->SetTitle("Channel");
  cADCt->GetYaxis()->SetTitle("ADC Time (ns)");
  cADCt->SetMarkerStyle(20); // idx 20 Circles, idx 21 Boxes
  cADCt->Write("cADCt_off");

  TGraphErrors *cTvseP0 = new TGraphErrors( kNcell, X, TvseP0, posErr, posErr );
  cTvseP0->GetXaxis()->SetLimits(0,kNcell);  
  //cTvsE->GetYaxis()->SetLimits(-0.5,0.5);
  cTvseP0->SetTitle("Timewalk Fit Constant P0");
  cTvseP0->GetXaxis()->SetTitle("Channel");
  cTvseP0->GetYaxis()->SetTitle("ns/GeV");
  cTvseP0->SetMarkerStyle(20); // idx 20 Circles, idx 21 Boxes
  cTvseP0->Write("cTvseP0");
  
  // Create fit parameter files
  ofstream P0;
  ofstream P1;
  ofstream P2;
  ofstream off;
  //Write to outfiles if fresh = 0

  if( fresh==1 ){
    //Write first parameter
    P0.open( P0path );
    P0 << "#HCal cosmic fit parameters obtained " << date.c_str() << endl;
    P0 << "#Exponential fit for i PMTs -> y = P0_i*exp(-P1_i*x)+P2_i, P0_i = " << endl;

    for( Int_t i=0; i<kNcell; i++ ){   
      P0 << TvseP0[i] << endl;
    }

    P0.close();

    //Write second parameter
    P1.open( P1path );
    P1 << "#HCal cosmic fit parameters obtained " << date.c_str() << endl;
    P1 << "#Exponential fit for i PMTs -> y = P0_i*exp(-P1_i*x)+P2_i, P1_i = " << endl;

    for( Int_t i=0; i<kNcell; i++ ){   
      P1 << TvseP1[i] << endl;
    }

    P1.close();

    //Write third parameter
    P2.open( P2path );
    P2 << "#HCal cosmic fit parameters obtained " << date.c_str() << endl;
    P2 << "#Exponential fit for i PMTs -> y = P0_i*exp(-P1_i*x)+P2_i, P2_i = " << endl;

    for( Int_t i=0; i<kNcell; i++ ){   
      P2 << TvseP2[i] << endl;
    }

    P2.close();
  }

  //Write TDC offsets
  off.open( offpath );
  off << "#HCal TDC offsets obtained " << date.c_str() << endl;
  off << "#Offsets obtained from fits over cosmic TDC distributions = " << endl;
  
  for( Int_t i=0; i<kNcell; i++ ){   
    off << TDCoffsets[i]/TDCCalib + oldTDCoffsets[i] << endl;
  }
  
  off.close();
  /*
  //Make canvas to hold all fits for comparison to HCal geometry
  TCanvas *TDC_top = new TCanvas("TDC_top","TDC_top",1600,1200);
  TCanvas *TDC_bot = new TCanvas("TDC_bot","TDC_bot",1600,1200);
  TCanvas *ADCt_top = new TCanvas("ADCt_top","ADCt_top",1600,1200);
  TCanvas *ADCt_bot = new TCanvas("ADCt_bot","ADCt_bot",1600,1200);

  TDC_top->Divide(12,12);
  TDC_bot->Divide(12,12);
  ADCt_top->Divide(12,12);
  ADCt_bot->Divide(12,12);

  gStyle->SetOptStat(0);
  for(Int_t i=0; i<kNcell; i++){
    TDC_top->cd(i+1);
    if( i>=144 ){
      TDC_bot->cd(i-143);
      gStyle->SetOptStat(0);
    }
    if(htdcDiff_corr[i]->GetEntries()<tFitMin){
      htdcDiff_corr[i]->SetAxisColor(2);
    }else{
      htdcDiff_corr[i]->SetAxisColor(1);
    }
    htdcDiff_corr[i]->Draw();
  }

  gStyle->SetOptStat(0);
  for(Int_t i=0; i<kNcell; i++){
    ADCt_top->cd(i+1);
    if( i>=144 ){
      ADCt_bot->cd(i-143);
      gStyle->SetOptStat(0);
    }
    if(hatime[i]->GetEntries()<tFitMin){
      hatime[i]->SetAxisColor(2);
    }else{
      hatime[i]->SetAxisColor(1);
    }
    hatime[i]->Draw();
  }
  */
  //Write out diagnostic histos and print to console
  fout->Write();

  cout << endl << endl << "TDC Offsets: " << endl << endl;
  
  Int_t cell = 0;
  for( Int_t r=0; r<kNrows; r++ ){
    for( Int_t c=0; c<kNcols; c++ ){
      cout << TDCoffsets[cell]/TDCCalib + oldTDCoffsets[cell] << " ";
      cell++;
    }
    cout << endl;
  }

  cout << endl << endl << "Total events analyzed: " << Nevents << "." << endl << endl;

  cout << "Timing offset analysis complete." << endl << endl;

  st->Stop();

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

}
