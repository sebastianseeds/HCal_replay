//SSeeds 9.7.22 - Post-production - Code to parse several runs back by making cuts on elastics and produce a single root file for further analysis. Configured to read in many settings of SBS and write HCAL-projected x and y to tree. 
#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <TSystem.h>
#include <TProfile.h>
#include <stdio.h>
#include <stdlib.h>
#include "TChain.h"
#include "TSystem.h"
#include "TStopwatch.h"
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
#include "TLine.h"
#include "TEllipse.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "pelas.h"
#include "TTreeFormula.h"
//#include "THaRun.h"
//#include "Event.h"
//#include "includes.h"


//Declare Functions
vector <Double_t> getDBParam( const char *file="", const char *param="", Int_t skip_lines = 0, Int_t chan = 0 );
string getDate();

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

vector<Double_t> getDBParam( string file="", string param="", Int_t skip_lines = 0, Int_t chan = 0 ){

  vector<Double_t> vec;

  cout << "Loading constants from database file: " << file << ".." << endl;
  ifstream inConstFile( file );
  if( !inConstFile ){
    cerr << endl << "ERROR: No input constant file present at " << file << "." << endl;
    return vec;
  }

  Int_t n0 = 0;
  Double_t d0 = 0;
  string line0;
  bool start_read = false;

  while( getline( inConstFile, line0 ) ){

    if( n0 == chan ) continue;

    TString Tline = (TString)line0;
    if( Tline.BeginsWith( param ) ) start_read = true;
    if( start_read == true && skip_lines > 0 ) skip_lines--;
    if( start_read == true && skip_lines == 0 ){
      istringstream iss( line0 );
      while( iss >> d0 ){
	vec.push_back(d0);
	n0++;
      }
    }
  }
  return vec;
}

void pelas_sh( const char *configfilename="spelas_sh.cfg", const char *outputfilename="pelas_shOut.root" ){
  
  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  // Get the date
  string date = getDate();

  // Declare Chain for many root files
  TChain *C = new TChain("T");

  // Declare general physics parameters to be modified by input config file
  Int_t kine = 8; // Keep track of kinematic being analyzed
  Double_t E_e = 1.916; // Energy of beam (incoming electrons from accelerator)
  Double_t t_trig = 510; // Time difference between BBCal Hi (L1A) and HCal Trig (ns)
  Double_t BB_d = 1.7988; // Distance to bigbite spectrometer from target chamber (m)
  Double_t BB_th = 36.; // Angle BB spectrometer makes with exit beamline
  Double_t HCal_d = 14.5; // Distance to HCal from scattering chamber for comm1
  Double_t HCal_th = 35.; // Angle HCal center makes with exit beamline  
  Double_t W_mean = 0.93; // Mean of W at current kinematic
  Double_t W_sig = 0.039; // Width of W at current kinematic
  Int_t magSet = 30;  // SBS Magnetic field strength (percent)
  Int_t useAlshield = 0; //Use 1/8" al shield on scattering chamber exit? 0:no
  string tar = "LH2";
  vector<TString> log;

  // Reading config file
  ifstream configfile(configfilename);
  TString currentline;
  cout << endl << "Chaining the following runs: " << endl;
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      if(!currentline) cout << "WARNING: No file exists at " << currentline << "." << endl;
      C->Add(currentline);
      if(currentline){
	log.push_back(currentline);
	cout << currentline << " ..check" << endl;
      }
    }    
  }

  TCut globalcut = "";
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endcut") ){
    if( !currentline.BeginsWith("#") ){
      globalcut += currentline;
    }    
  }

  cout << endl << "Loading the following parameters from " << configfilename << ":" << endl;

  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("#") ){
    TObjArray *tokens = currentline.Tokenize(" ");
    Int_t ntokens = tokens->GetEntries();
    if( ntokens>1 ){
      TString skey = ( (TObjString*)(*tokens)[0] )->GetString();
      if( skey == "kine" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	kine = sval.Atoi();
	cout << "Kinematic Setting: " << kine << endl;
      }
      if( skey == "E_e" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	E_e = sval.Atof();
	cout << "Beam Energy: " << E_e << endl;
      }
      if( skey == "t_trig" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	t_trig = sval.Atof();
	cout << "BBCal/HCal Trigger Difference: " << E_e << endl;
      }
      if( skey == "BB_d" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	BB_d = sval.Atof();
	cout << "BigBite Spectrometer distance: " << BB_d << endl;
      }
      if( skey == "BB_th" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	BB_th = sval.Atof() * TMath::DegToRad();
	cout << "BigBite Spectrometer angle (rad): " << BB_th << endl;
      }
      if( skey == "HCal_d" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	HCal_d = sval.Atof();
	cout << "HCal distance: " << HCal_d << endl;
      }
      if( skey == "HCal_th" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	HCal_th = sval.Atof() * TMath::DegToRad();	
	cout << "HCal angle (rad): " << HCal_th << endl;
      }
      if( skey == "W_mean" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        W_mean = sval.Atof();
	cout << "W mean cut: " << W_mean << endl;
      }
      if( skey == "W_sig" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        W_sig = sval.Atof();
	cout << "W sigma cut: " << W_sig << endl;
      } 
      if( skey == "useAlshield" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	useAlshield = sval.Atoi();
	cout << "Aluminum absorber option: " << useAlshield << endl;
      } 
      if( skey == "tar" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	tar = sval;
	cout << "target: " << tar << endl;
      }
    }
    delete tokens;
  }
  
  cout << "Setup parameters loaded." << endl;

  // Paths for DB constants - should update for sbs-gmn
  string HCalConstPath = "/w/halla-scshelf2102/sbs/SBS_REPLAY/SBS-replay/DB/db_sbs.hcal.dat";
  string BBCalShConstPath = "/w/halla-scshelf2102/sbs/SBS_REPLAY/SBS-replay/DB/db_bb.sh.dat";
  string BBCalPSConstPath = "/w/halla-scshelf2102/sbs/SBS_REPLAY/SBS-replay/DB/db_bb.ps.dat";
  string GEMConstPath = "/w/halla-scshelf2102/sbs/SBS_REPLAY/SBS-replay/DB/db_bb.gem.dat";
  string HodoConstPath = "/w/halla-scshelf2102/sbs/SBS_REPLAY/SBS-replay/DB/db_bb.hodotdc.dat";
  string GRINCHConstPath = "/w/halla-scshelf2102/sbs/SBS_REPLAY/SBS-replay/DB/db_bb.grinch_tdc.dat";

  ///////////////
  ////HCal params
  //ADC time offsets
  string HCal_ADCtoff = "sbs.hcal.adc.timeoffset";
  vector<Double_t> HCal_ADCtoffsets = getDBParam( HCalConstPath, HCal_ADCtoff, 1, maxHCalChan );
  //TDC offsets
  string HCal_TDCoff = "sbs.hcal.tdc.offset";
  vector<Double_t> HCal_TDCoffsets = getDBParam( HCalConstPath, HCal_TDCoff, 1, maxHCalChan );
  //ADC gain coefficients
  string HCal_ADCg = "-------[ 2021-10-24 04:30:00 ]";
  vector<Double_t> HCal_ADCgains_211024043000 = getDBParam( HCalConstPath, HCal_ADCg, 2, maxHCalChan );

  ////////////////
  ////BBCal params
  //ADC gain coefficients
  //Shower 2022-02-01 22:00:00
  string BBCalSh_ADCg = "-------[ 2022-02-01 22:00:00 ]";
  vector<Double_t> BBCalSh_ADCgains_220201220000 = getDBParam( BBCalShConstPath, BBCalSh_ADCg, 2, maxBBCalShChan );
  //Preshower 2022-02-01 22:00:00
  string BBCalPS_ADCg = "-------[ 2022-02-01 22:00:00 ]";
  vector<Double_t> BBCalPS_ADCgains_220201220000 = getDBParam( BBCalPSConstPath, BBCalPS_ADCg, 2, maxBBCalPSChan );
  
  /*
  cout << endl << endl << "Test offsets: " << endl;

  Int_t kNrows = 26;
  Int_t kNcols = 2;

  for( Int_t r=0; r<kNrows; r++){
    for( Int_t c=0; c<kNcols; c++){
      Int_t i = r*kNcols+c;
      cout << BBCalPS_ADCgains_220201220000[i] << " ";
    }
    cout << endl;
  }
  */

  cout << endl << "Database parameters loaded." << endl;

  TTreeFormula *GlobalCut = new TTreeFormula( "GlobalCut", globalcut, C );

  cout << "Elastic global cut obtained." << endl;

  //Set simple structures for complex branches
  typedef struct { 
    Double_t Ep;
    Double_t KEp;
    Double_t dx;
    Double_t dy;
    Double_t W2;
    Double_t etheta;
    Double_t ephi;
    Double_t ebeam;
    Double_t Q2recon;
    Double_t W2recon;
    Double_t nurecon;
    Double_t Eprecon;
    Double_t Epelastic;
    Double_t Epincident;
    Double_t Hx_exp;
    Double_t Hy_exp;
  } PHYSICS;
  
  PHYSICS physics;

  typedef struct {
    UInt_t run;
    UInt_t trigger;
    UInt_t event;
  } RUNDATA;

  RUNDATA rundata;

  typedef struct {
    Double_t event;
    Double_t charge;
    Double_t clock;
  } SCALAR;

  SCALAR scalar;

  /*
  typedef struct {
    Int_t chan;
    Int_t rows;
    Int_t cols;
    Double_t adcoff
  } HCALDET;

  HCALDET hcal_det;
  */

  typedef struct {    
    Double_t vz;
    Double_t vy;
    Double_t vx;
    Double_t pz;
    Double_t py;
    Double_t px;
    Double_t p;
    Int_t n;
  } TRACK;

  TRACK track;
  
  typedef struct {
    Double_t x;
    Double_t y;
    Double_t e;
    Int_t nblk;
    Int_t nclus;
    Double_t cblk_e[ maxHCalChan ];
    Double_t cblk_tdc[ maxHCalChan ];
    Double_t cblk_adct[ maxHCalChan ];
    Double_t catime[ maxClusters ];
    Double_t ce[ maxClusters ];
    Double_t ceblk[ maxClusters ];
    Double_t cid[ maxClusters ];
    Double_t cnblk[ maxClusters ];
    Double_t ctdc[ maxClusters ];
    Double_t cx[ maxClusters ];
    Double_t cy[ maxClusters ];
  } HCAL;

  HCAL hcal;

  typedef struct {
    Double_t x;
    Double_t y;
    Double_t e;
    Int_t nblk;
    Int_t nclus;
    Double_t cblk_e[ maxBBCalPSChan ];
    Double_t cblk_ec[ maxBBCalPSChan ];
    Double_t cblk_atime[ maxBBCalPSChan ];
  } BBCALPS;

  BBCALPS bbcalps;

  typedef struct {
    Double_t x;
    Double_t y;
    Double_t e;
    Int_t nblk;
    Int_t nclus;
    Double_t cblk_e[ maxBBCalShChan ];
    Double_t cblk_ec[ maxBBCalShChan ];
    Double_t cblk_atime[ maxBBCalShChan ];
  } BBCALSH;

  BBCALSH bbcalsh;
  
  //TSsbs tree struct for SCALAR data extraction and interpolation
  typedef struct{
    UInt_t rnumber;
    Double_t enumber;
    Double_t aggcharge;
    Double_t clock;
  } SCALARAGG;

  vector<SCALARAGG> Sagg;
  
  /////////////////////////////////
  ////Declare input tree parameters
  //Event and physics
  UInt_t runI; 
  UInt_t runN, TBits; 
  //ULong64_t runT;

  Double_t kineW2;
  //BigBite Tracks
  Double_t BBtr_p[maxTracks], BBtr_px[maxTracks], BBtr_py[maxTracks], BBtr_pz[maxTracks];
  Double_t BBtr_vx[maxTracks], BBtr_vy[maxTracks], BBtr_vz[maxTracks], BBtr_chi2[maxTracks];
  Double_t BBfp_x[maxTracks], BBfp_y[maxTracks], BBfp_th[maxTracks], BBfp_ph[maxTracks];
  Double_t BBtgt_x[maxTracks], BBtgt_y[maxTracks], BBtgt_th[maxTracks], BBtgt_ph[maxTracks];
  //BBCal
  Double_t BBtr_n, BBps_x, BBps_y, BBps_e, BBps_nblk, BBps_nclus, BBsh_x, BBsh_y, BBsh_e, BBsh_nblk, BBsh_nclus, BB_eOp;
  Double_t BBps_cblk_e[maxBBCalPSChan], BBps_cblk_ec[maxBBCalPSChan], BBps_cblk_atime[maxBBCalPSChan], BBps_cblk_id[maxBBCalPSChan];
  Double_t BBsh_cblk_e[maxBBCalShChan], BBsh_cblk_ec[maxBBCalShChan], BBsh_cblk_atime[maxBBCalShChan], BBsh_cblk_id[maxBBCalShChan];
  //GEM
  Double_t BBgem_nplanes;
  //Trigger TDC
  Double_t TDCT_id[maxTDCTrigChan], TDCT_tdc[maxTDCTrigChan]; 
  Int_t TDCTndata;
  //HCal
  Double_t HCAL_x, HCAL_y, HCAL_e, HCAL_nblk, HCAL_nclus;
  Double_t HCAL_cblk_tdctime[maxHCalChan], HCAL_cblk_atime[maxHCalChan], HCAL_cblk_id[maxHCalChan], HCAL_cblk_e[maxHCalChan]; 
  Double_t HCAL_ce[maxClusters], HCAL_ceblk[maxClusters], HCAL_catime[maxClusters], HCAL_cid[maxClusters], HCAL_cnblk[maxClusters], HCAL_ctdc[maxClusters], HCAL_cx[maxClusters], HCAL_cy[maxClusters];

  // Declare root tree variables and set values to memory locations in root file
  C->SetBranchStatus( "*", 0 );  //Switch everything off first to save processing time
  // HCal
  C->SetBranchStatus( "sbs.hcal.x", 1 ); //Switch the essential ones on individually
  C->SetBranchStatus( "sbs.hcal.y", 1 );
  C->SetBranchStatus( "sbs.hcal.e", 1 );
  C->SetBranchStatus( "sbs.hcal.nblk", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.tdctime", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.atime", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.id", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.e", 1 );
  C->SetBranchStatus( "sbs.hcal.clus.e", 1 );
  C->SetBranchStatus( "sbs.hcal.clus.eblk", 1 );
  C->SetBranchStatus( "sbs.hcal.clus.atime", 1 );
  C->SetBranchStatus( "sbs.hcal.clus.id", 1 );
  C->SetBranchStatus( "sbs.hcal.clus.nblk", 1 );
  C->SetBranchStatus( "sbs.hcal.clus.tdctime", 1 );
  C->SetBranchStatus( "sbs.hcal.clus.x", 1 );
  C->SetBranchStatus( "sbs.hcal.clus.y", 1 );
  C->SetBranchStatus( "sbs.hcal.nclus", 1 );
  // BB track
  C->SetBranchStatus( "bb.tr.chi2", 1 );
  C->SetBranchStatus( "bb.tr.n", 1 );
  C->SetBranchStatus( "bb.tr.px", 1 );
  C->SetBranchStatus( "bb.tr.py", 1 );
  C->SetBranchStatus( "bb.tr.pz", 1 );    
  C->SetBranchStatus( "bb.tr.p", 1 );
  C->SetBranchStatus( "bb.tr.vx", 1 );
  C->SetBranchStatus( "bb.tr.vy", 1 );
  C->SetBranchStatus( "bb.tr.vz", 1 );
  C->SetBranchStatus( "bb.tr.r_x", 1 );
  C->SetBranchStatus( "bb.tr.r_y", 1 );
  C->SetBranchStatus( "bb.tr.r_th", 1 );
  C->SetBranchStatus( "bb.tr.r_ph", 1 );
  C->SetBranchStatus( "bb.tr.tg_x", 1 );
  C->SetBranchStatus( "bb.tr.tg_y", 1 );
  C->SetBranchStatus( "bb.tr.tg_th", 1 );
  C->SetBranchStatus( "bb.tr.tg_ph", 1 );
  // BBCal shower preshower
  C->SetBranchStatus( "bb.ps.e", 1 );
  C->SetBranchStatus( "bb.ps.x", 1 );
  C->SetBranchStatus( "bb.ps.y", 1 );
  C->SetBranchStatus( "bb.sh.e", 1 );
  C->SetBranchStatus( "bb.sh.x", 1 );
  C->SetBranchStatus( "bb.sh.y", 1 );
  C->SetBranchStatus( "bb.ps.clus_blk.e", 1 );
  C->SetBranchStatus( "bb.ps.clus_blk.e_c", 1 );
  C->SetBranchStatus( "bb.ps.clus_blk.atime", 1 );
  C->SetBranchStatus( "bb.ps.clus_blk.id", 1 );
  C->SetBranchStatus( "bb.ps.nblk", 1 );
  C->SetBranchStatus( "bb.ps.nclus", 1 );
  C->SetBranchStatus( "bb.sh.clus_blk.e", 1 );
  C->SetBranchStatus( "bb.sh.clus_blk.e_c", 1 );
  C->SetBranchStatus( "bb.sh.clus_blk.atime", 1 );
  C->SetBranchStatus( "bb.sh.clus_blk.id", 1 );
  C->SetBranchStatus( "bb.sh.nblk", 1 );
  C->SetBranchStatus( "bb.sh.nclus", 1 );
  C->SetBranchStatus( "bb.etot_over_p", 1 );
  // GEM
  C->SetBranchStatus( "bb.gem.track.nhits", 1 );
  // Trigger TDC
  C->SetBranchStatus( "bb.tdctrig.tdc", 1 );
  C->SetBranchStatus( "bb.tdctrig.tdcelemID", 1 );
  C->SetBranchStatus( "Ndata.bb.tdctrig.tdcelemID", 1 );
  // Event and physics
  C->SetBranchStatus( "Event_Branch", 1 );
  C->SetBranchStatus( "e.kine.W2", 1 );
  C->SetBranchStatus( "fEvtHdr.fRun", 1 );
  //C->SetBranchStatus( "fEvtHdr.fEvtTime", 1 );
  C->SetBranchStatus( "fEvtHdr.fEvtNum", 1 );
  C->SetBranchStatus( "fEvtHdr.fTrigBits", 1 );

  // Set the variable addresses
  // HCal
  C->SetBranchAddress( "sbs.hcal.x", &HCAL_x );
  C->SetBranchAddress( "sbs.hcal.y", &HCAL_y );
  C->SetBranchAddress( "sbs.hcal.e", &HCAL_e );
  C->SetBranchAddress( "sbs.hcal.nblk", &HCAL_nblk );
  C->SetBranchAddress( "sbs.hcal.clus_blk.tdctime", HCAL_cblk_tdctime );
  C->SetBranchAddress( "sbs.hcal.clus_blk.atime", HCAL_cblk_atime );
  C->SetBranchAddress( "sbs.hcal.clus_blk.id", HCAL_cblk_id );
  C->SetBranchAddress( "sbs.hcal.clus_blk.e", HCAL_cblk_e );
  C->SetBranchAddress( "sbs.hcal.clus.e", HCAL_ce );
  C->SetBranchAddress( "sbs.hcal.clus.eblk", HCAL_ceblk );
  C->SetBranchAddress( "sbs.hcal.clus.atime", HCAL_catime );
  C->SetBranchAddress( "sbs.hcal.clus.id", HCAL_cid );
  C->SetBranchAddress( "sbs.hcal.clus.nblk", HCAL_cnblk );
  C->SetBranchAddress( "sbs.hcal.clus.tdctime", HCAL_ctdc );
  C->SetBranchAddress( "sbs.hcal.clus.x", HCAL_cx );
  C->SetBranchAddress( "sbs.hcal.clus.y", HCAL_cy );
  C->SetBranchAddress( "sbs.hcal.nclus", &HCAL_nclus );
  // BB track
  C->SetBranchAddress( "bb.tr.chi2", BBtr_chi2 );
  C->SetBranchAddress( "bb.tr.n", &BBtr_n );
  C->SetBranchAddress( "bb.tr.px", BBtr_px );
  C->SetBranchAddress( "bb.tr.py", BBtr_py );
  C->SetBranchAddress( "bb.tr.pz", BBtr_pz );
  C->SetBranchAddress( "bb.tr.p", BBtr_p );
  C->SetBranchAddress( "bb.tr.vx", BBtr_vx );
  C->SetBranchAddress( "bb.tr.vy", BBtr_vy );
  C->SetBranchAddress( "bb.tr.vz", BBtr_vz );
  C->SetBranchAddress( "bb.tr.r_x", BBfp_x );
  C->SetBranchAddress( "bb.tr.r_y", BBfp_y );
  C->SetBranchAddress( "bb.tr.r_th", BBfp_th );
  C->SetBranchAddress( "bb.tr.r_ph", BBfp_ph );
  C->SetBranchAddress( "bb.tr.tg_x", BBtgt_x );
  C->SetBranchAddress( "bb.tr.tg_y", BBtgt_y );
  C->SetBranchAddress( "bb.tr.tg_th", BBtgt_th );
  C->SetBranchAddress( "bb.tr.tg_ph", BBtgt_ph );
  // BBCal shower preshower
  C->SetBranchAddress( "bb.ps.e", &BBps_e );
  C->SetBranchAddress( "bb.ps.x", &BBps_x );
  C->SetBranchAddress( "bb.ps.y", &BBps_y );
  C->SetBranchAddress( "bb.sh.e", &BBsh_e );
  C->SetBranchAddress( "bb.sh.x", &BBsh_x );
  C->SetBranchAddress( "bb.sh.y", &BBsh_y );
  C->SetBranchAddress( "bb.ps.clus_blk.e", BBps_cblk_e );
  C->SetBranchAddress( "bb.ps.clus_blk.e_c", BBps_cblk_ec );
  C->SetBranchAddress( "bb.ps.clus_blk.atime", BBps_cblk_atime );
  C->SetBranchAddress( "bb.ps.clus_blk.id", BBps_cblk_id );
  C->SetBranchAddress( "bb.ps.nblk", &BBps_nblk );
  C->SetBranchAddress( "bb.ps.nclus", &BBps_nclus );
  C->SetBranchAddress( "bb.sh.clus_blk.e", BBsh_cblk_e );
  C->SetBranchAddress( "bb.sh.clus_blk.e_c", BBsh_cblk_ec );
  C->SetBranchAddress( "bb.sh.clus_blk.atime", BBsh_cblk_atime );
  C->SetBranchAddress( "bb.sh.clus_blk.id", BBsh_cblk_id );
  C->SetBranchAddress( "bb.sh.nblk", &BBsh_nblk );
  C->SetBranchAddress( "bb.sh.nclus", &BBsh_nclus );
  C->SetBranchAddress( "bb.etot_over_p", &BB_eOp );
  // GEM
  C->SetBranchAddress( "bb.gem.track.nhits", &BBgem_nplanes );
  // Trigger TDC
  C->SetBranchAddress( "bb.tdctrig.tdcelemID", TDCT_id );
  C->SetBranchAddress( "bb.tdctrig.tdc", TDCT_tdc );
  C->SetBranchAddress( "Ndata.bb.tdctrig.tdcelemID", &TDCTndata );
  // Event and physics
  C->SetBranchAddress( "e.kine.W2", &kineW2 );
  C->SetMakeClass(1);
  C->SetBranchAddress( "fEvtHdr.fRun", &runI );
  //C->SetBranchAddress( "fEvtHdr.fEvtTime", &runT );
  C->SetBranchAddress( "fEvtHdr.fEvtNum", &runN );
  C->SetBranchAddress( "fEvtHdr.fTrigBits", &TBits );

  cout << "Tree variables linked." << endl;

  // Declare outfiles
  //outputfilename = Form( "/lustre19/expphy/volatile/halla/sbs/seeds/GMN_parsed/e1209019_parsed_short_SBS%d_tar%s_%s.root", kine, tar.c_str(), date.c_str() );
  outputfilename = "test_out.root";
  TFile *fout = new TFile( outputfilename, "RECREATE" );
  string logpath = Form( "/lustre19/expphy/volatile/halla/sbs/seeds/GMN_parsed/logs/parsedLog_%s.txt", date.c_str() );

  // Initialize misc. variables
  Int_t elasYield = 0; // Keep track of total elastics analyzed
  Double_t pBeam = E_e/(1.+E_e/M_p*(1.-cos(BB_th)));  
  Double_t Eloss_outgoing = celldiameter/2.0/sin(BB_th) * rho_tgt * dEdx_tgt; //Approximately 1 MeV
  if( useAlshield != 0 ) Eloss_outgoing += Alshieldthick * rho_Al * dEdx_Al;

  // Declare Histograms
  TH1D *h_W2 = new TH1D("h_W2",";W2 (GeV^2);",250,0,2);
  TH1D *h_W = new TH1D("h_W",";W (GeV);",250,0,2);
  TH2D *h_charge = new TH2D("h_charge","Accumulated Charge vs CODA Run; Run Number; Coulombs", 10000, 10000, 20000, 1000, 0, 10000 ); 
  //TH2D *h_test = new TH2D("h_test","Scalar Clock vs Events", 6000, 0, 6000000, 100000, 0, 100000 );

  // Create output tree
  TTree *AT = new TTree("AT","Tree containing aggregate variables and analysis parameters with simple elastic cuts.");

  ////////////////
  ////Set branches
  //Physics
  AT->Branch( "physics", &physics, "Ep/D:KEp/D:dx/D:dy/D:W/D:etheta/D:ephi/D:ebeam/D:Q2recon/D:W2recon/D:nurecon/D:Eprecon/D:Epelastic/D:Epincident/D:Hx_exp/D:Hy_exp/D" );
  //Run Data
  AT->Branch( "rundata", &rundata, "run/i:trigger/i:event/i" );
  //Scalar
  AT->Branch( "scalar", &scalar, "event/D:charge/D:clock/D" );
  //Tracks
  AT->Branch( "track", &track, "vz/D:vy/D:vx/D:pz/D:py/D:px/D:p/D:n/I" );
  //HCal
  AT->Branch( "hcal", &hcal, Form("x/D:y/D:e/D:nblk/I:nclus/I:cblk_e[%d]/D:cblk_tdc[%d]/D:cblk_adct[%d]/D:catime[%d]/D:ce[%d]/D:ceblk[%d]/D:cid[%d]/D:cnblk[%d]/D:ctdc[%d]/D:cx[%d]/D:cy[%d]/D",maxHCalChan,maxHCalChan,maxHCalChan,maxClusters,maxClusters,maxClusters,maxClusters,maxClusters,maxClusters,maxClusters,maxClusters) );
  //BBCal Preshower
  AT->Branch( "bbcalps", &bbcalps, Form("x/D:y/D:e/D:nblk/I:nclus/I:cblk_e[%d]/D:cblk_ec[%d]/D:cblk_atime[%d]/D",maxBBCalPSChan,maxBBCalPSChan,maxBBCalPSChan) );
  //BBCal Shower
  AT->Branch( "bbcalsh", &bbcalsh, Form("x/D:y/D:e/D:nblk/I:nclus/I:cblk_e[%d]/D:cblk_ec[%d]/D:cblk_atime[%d]/D",maxBBCalShChan,maxBBCalShChan,maxBBCalShChan) );

  //Get charge and bcm data for all runs
  //getCharge( log, Sagg);

  //Create chain for epics/scalar tree
  TChain *TS = new TChain("TSsbs");
  //Chain together all files
  for( Int_t l=0; l<log.size(); l++){
    TS->Add( log[l] );
  }
  //Declare dummies and counters for needed vars
  UInt_t Srunnum = 0;
  Double_t Sevnum = 0., Sdnew = 0., Sclk = 0., clkLast = 0., clkOff = 0., clkTCoff = 0., clkDiff = 0.;
  Int_t DumTN = 0, DumCurTN = 0, DumEvLast = 0, DumEvOff = 0, DumRun = 0, evDiff = 0, evTCoff = 0;
  Int_t fileN = -1;
  Long64_t CurEv = 0;
  Long64_t NumEv = TS->GetEntries();
  bool treechange = false;
  bool begin = true;
  Double_t chargeLast = 0;
  Double_t testevent = 0;
  //Switch on and set addresses for branches
  TS->SetBranchStatus( "*", 0 );
  TS->SetBranchStatus( "evNumber", 1 );
  TS->SetBranchStatus( "sbs.bcm.dnew.cnt", 1 );
  TS->SetBranchStatus( "sbs.104kHz_CLK.cnt", 1 );
  TS->SetMakeClass(1);
  TS->SetBranchAddress( "evNumber", &Sevnum );
  TS->SetBranchAddress( "sbs.bcm.dnew.cnt", &Sdnew );
  TS->SetBranchAddress( "sbs.104kHz_CLK.cnt", &Sclk );
  //Loop over events (dt=2s for epics tree)
  while( TS->GetEntry( CurEv++ ) ){

    DumCurTN = TS->GetTreeNumber();
    //Get the current run number if the tree changes
    if( CurEv == 1 || DumCurTN != DumTN ){

      DumEvOff += DumEvLast; //Get the last event number from the previous tree
      clkOff += clkLast; //Get the last corrected clock value from the previous segment

      //cout << "Segment Change!" << endl;
      treechange = true;
      fileN++;
      DumTN = DumCurTN;
      TString testline = (string)log[ fileN ]; 
      TFile *f1 = TFile::Open( testline );
      THaRun *run;
      f1->GetObject("Run_Data",run);
      Srunnum = run->GetNumber();
 
      if(begin==true){ //If it is the first event of the chain, record the run number for comparison
	DumRun = run->GetNumber();
	begin=false;
      }
      if(DumRun != run->GetNumber()){ //If the run number changes, reset both clock and event counter and write charge to histogram
	DumEvOff = 0;
	clkOff = 0.;
	h_charge->Fill(DumRun,chargeLast);
	DumRun = run->GetNumber();
      }
    }
    
    if( treechange==false ){
      evDiff = Sevnum - DumEvLast; //Get the difference between the last two event numbers
      clkDiff = Sclk/103700 - clkLast; //Same for clock
      clkTCoff = 0;
      evTCoff = 0;
    }else{
      evTCoff = evDiff-Sevnum; //When the segment changes, get event number correction to offset
      if( Sevnum==1 ) evTCoff = 0; //Correct for the first event in the chain
      clkTCoff = clkDiff-Sclk/103700; //Same for clock

    }

    DumEvLast = Sevnum;
    clkLast = Sclk/103700;
    chargeLast = Sdnew/3318;

    //Fill the scalar data structure with charge and clock corrections applied
    SCALARAGG thisSCALAR = { Srunnum, Sevnum+DumEvOff+evTCoff, Sdnew/3318, Sclk/103700+clkOff+clkTCoff };
    Sagg.push_back( thisSCALAR );

    testevent = Sevnum+DumEvOff+evTCoff;

    treechange = false;
  }

  cout << DumRun << " " << chargeLast << endl;
  h_charge->Fill(DumRun,chargeLast);
  cout << "Total events counted on scalar tree: " << testevent << "." << endl;

  cout << "Scalar and accumulated charge data extracted." << endl;

  // Set long int to keep track of total entries in T tree
  Long64_t Nevents = C->GetEntries();

  // Set variables to keep track of event sequence
  Int_t treenum = 0, currenttreenum = 0, currentrun = 0;
  UInt_t runnumb = 0;
  Int_t SCALARidx = 0, EvLast = 0, EvOff = 0;
  Double_t SCALAR_event = 0, SCALAR_charge = 0, SCALAR_clock = 0;
  SCALARAGG curSCALAR;
  Int_t dummy = 0;
  TString message = "";
  begin=true;
  DumRun=0;
  evTCoff=0;
  testevent=0;
  treechange=false;

  // Set variable to keep track of charge
  Double_t accCharge = 0.;

  cout << endl << "Initialization complete." << endl << endl;
  cout << "Opened up tree with nentries: " << C->GetEntries() << "." << endl << endl;

  //Loop over events
  TStopwatch *sw = new TStopwatch();
  sw->Start( kTRUE );
  cout << "Main loop over all data commencing.." << endl;
  Double_t progress = 0.;
  Double_t timekeeper = 0., timeremains = 0.;
  while(progress<1.0){
    Int_t barwidth = 50;
    for(Long64_t nevent = 1; nevent<Nevents; nevent++){
      
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
	timeremains = timekeeper*( double(Nevents)/double(nevent) - 1. ); 
      sw->Reset();
      sw->Continue();
      progress = (double)((nevent+1.)/Nevents);
      cout << "] " << int(progress*100.) << "%, elastic events: " << elasYield << "/" << Nevents << ", " << message << ", time remaining: " << int(timeremains/60.) << "m \r";
      cout.flush();
      
      //Get each event for analysis
      C->GetEntry( nevent );


      ////////
      ////Cuts
      //Use TTreeFormula to apply globalcut to speed up processing and cut on track qty=1
      currenttreenum = C->GetTreeNumber();
      if( nevent == 1 || currenttreenum != treenum ){

	EvOff += EvLast;

	//cout << "Segment Change!" << endl;
	treechange = true;
	treenum = currenttreenum; 
	GlobalCut->UpdateFormulaLeaves();
	message = Form("analyzing run %d", runI);
	
	if(begin==true){ //If it is the first event of the chain, record the run number for comparison
	  DumRun = runI;
	  begin=false;
	}
	if(DumRun != runI){ //If the run number changes, reset both clock and event counter
	  EvOff = 0;
	}

      }

      if( treechange==false ){
	evDiff = runN - EvLast; //Get the difference between the last two event numbers
	evTCoff = 0;
	//cout << evDiff << " " << runN << " " << EvLast;
      }else{
	evTCoff = evDiff-runN; //When the segment changes, get event number correction to offset
	if( runN==1 ) evTCoff = 0; //Correct for the first event in the chain
	treechange = false;
      }

      EvLast = runN;
      
      //Return SCALAR data for current index
      curSCALAR = Sagg[ SCALARidx ];
      //Check if T tree event is greater than TSsbs tree event and increment as necessary to continue interpolation
      while( runN+EvOff+evTCoff > (UInt_t)curSCALAR.enumber ){
	SCALARidx++;
	curSCALAR = Sagg[ SCALARidx ];
      }

      //Set SCALAR variables for T tree block where T_tree_event < TSsbs_tree_event
      SCALAR_event = curSCALAR.enumber;
      SCALAR_charge = curSCALAR.aggcharge;
      SCALAR_clock = curSCALAR.clock;      

      //Accumulate charge on current CODA run before any cuts are applied
      if( accCharge < Sdnew/chargeConvert ) accCharge = Sdnew/chargeConvert;

      bool failedglobal = GlobalCut->EvalInstance(0) == 0;
      if( int(BBtr_n) != 1 || failedglobal ) continue;
       
      //BBCal and HCal trigger coincidence
      Double_t bbcal_time=0., hcal_time=0.;
      for(Int_t ihit=0; ihit<TDCTndata; ihit++){
	if(TDCT_id[ihit]==5) bbcal_time=TDCT_tdc[ihit];
	if(TDCT_id[ihit]==0) hcal_time=TDCT_tdc[ihit];
      }
      Double_t diff = hcal_time - bbcal_time;
      //Cut on BB trigger/HCal trigger coincidence time
      if( fabs(diff-510.)>10. ) continue;

      ///////////////
      ////Corrections
      //Correct the beam energy with energy loss in target using vertex position
      Double_t Eloss = (BBtr_vz[0]+l_tgt/2.0) * rho_tgt * dEdx_tgt + uwallthick_LH2 * rho_Al * dEdx_Al; //approximately 3 MeV
      Double_t E_corr = E_e - Eloss;

      //cout << runI << endl;
      //TBranch *bla = C->GetBranch("fEvtHdr.fRun");
      //TBranch *bla = C->GetBranch("sbs.hcal.e");
      //cout << kineW2 << endl;

      //////////////////////////////
      ////Start Physics Calculations
      Double_t p_corr = BBtr_p[0] - Eloss_outgoing; //Neglecting the mass of e'
      Double_t etheta = acos( BBtr_pz[0]/BBtr_p[0] ); //Will use the uncorrected track momentum to reconstruct e' theta
      Double_t ephi = atan2( BBtr_py[0], BBtr_px[0] );   
      TVector3 vertex( 0, 0, BBtr_vz[0] ); // z location of vertex in hall coordinates
      TLorentzVector Pbeam( 0, 0, E_corr, E_corr ); //Mass of e negligable
      TLorentzVector kprime( BBtr_px[0], BBtr_py[0], BBtr_pz[0], BBtr_p[0] );
      TLorentzVector Ptarg( 0, 0, 0, M_p );   
      TLorentzVector q = Pbeam - kprime;
      TLorentzVector PgammaN = Ptarg + q; //(-px, -py, ebeam - pz, Mp + ebeam - p)     
      Double_t pel = E_corr/( 1.+E_corr/M_p*( 1.-cos(etheta) ) );
      Double_t nu = E_corr - BBtr_p[0];
      Double_t pp = sqrt( pow(nu,2)+2.*M_p*nu );
      Double_t phinucleon = ephi + PI; //assume coplanarity
      Double_t thetanucleon = acos( (E_corr - BBtr_pz[0])/pp ); //use elastic constraint on nucleon kinematics     
      TVector3 pNhat( sin(thetanucleon)*cos(phinucleon), sin(thetanucleon)*sin(phinucleon), cos(thetanucleon) );
      
      //Define HCal coordinate system
      TVector3 HCAL_zaxis( sin(-HCal_th), 0, cos(-HCal_th) );
      TVector3 HCAL_xaxis( 0, -1, 0 );
      TVector3 HCAL_yaxis = HCAL_zaxis.Cross(HCAL_xaxis).Unit();      
      TVector3 HCAL_origin = HCal_d * HCAL_zaxis + HCALHeight * HCAL_xaxis;      
      Double_t sintersect = ( HCAL_origin - vertex ).Dot( HCAL_zaxis ) / (pNhat.Dot( HCAL_zaxis ) ); 
      TVector3 HCAL_intersect = vertex + sintersect * pNhat;   

      //Calculate quantities of interest
      Double_t yexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_yaxis );
      Double_t xexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_xaxis );  
      Double_t E_ep = sqrt( pow(M_e,2) + pow(BBtr_p[0],2) ); // Obtain the scattered electron energy 
      Double_t p_ep = BBtr_p[0];     
      Double_t Q2 = 2*E_corr*E_ep*( 1-(BBtr_pz[0]/p_ep) ); // Obtain Q2 from beam energy, outgoing electron energy, and momenta  
      Double_t W2 = kineW2; //Get invarient mass transfer W from the four-momentum of the scattered nucleon
      //Double_t W = PgammaN.M();
      //Double_t W = sqrt(kineW2);
      Double_t E_pp = nu+M_p; // Get energy of the proton
      Double_t Enucleon = sqrt(pow(pp,2)+pow(M_p,2)); // Check on E_pp, same     
      Double_t KE_p = nu; // For elastics
      Double_t dx = HCAL_x - xexpect_HCAL;
      Double_t dy = HCAL_y - yexpect_HCAL;  
      Double_t pelastic = E_corr/(1.+(E_corr/M_p)*(1.0-cos(etheta))); 	
      Double_t precon = BBtr_p[0] + Eloss_outgoing; //reconstructed momentum, corrected for mean energy loss exiting the target (later we'll also want to include Al shielding on scattering chamber)
      Double_t nu_recon = E_corr - precon;
      Double_t Q2recon = 2.0*E_corr*precon*(1.0-cos(etheta));
      Double_t W2recon = pow(M_p,2) + 2.0*M_p*nu_recon - Q2recon;

      //Obtain good clusters in HCal
      

      //Fill physics variables
      physics.Ep = E_pp;
      physics.KEp = KE_p;
      physics.dx = dx;
      physics.dy = dy;
      physics.W2 = W2;
      physics.etheta = etheta;
      physics.ephi = ephi;
      physics.ebeam = E_corr;
      physics.Q2recon = Q2recon;
      physics.W2recon = W2recon;
      physics.nurecon = nu_recon;
      physics.Eprecon = precon;
      physics.Epelastic = pelastic;
      physics.Epincident = pelastic - Eloss_outgoing;
      physics.Hx_exp = xexpect_HCAL;
      physics.Hy_exp = yexpect_HCAL;
 
      //Fill run data variables
      rundata.run = runI;
      rundata.trigger = TBits;
      rundata.event = runN+EvOff+evTCoff;

      testevent = runN+EvOff+evTCoff;

      //cout << runN+EvOff << endl;

      //Fill scalar variables
      scalar.event = SCALAR_event;
      scalar.charge = SCALAR_charge;
      scalar.clock = SCALAR_clock;

      //Fill HCal detector variables. TODO - Inefficient, don't need to save per event
      /*
      T_BB_d = BB_d;
      T_BB_th = BB_th;
      T_HCal_d = HCal_d;
      T_HCal_th = HCal_th;
      */

      //Fill track variables
      track.vz = BBtr_vz[0];
      track.vy = BBtr_vy[0];
      track.vx = BBtr_vx[0];
      track.pz = BBtr_pz[0];
      track.py = BBtr_py[0];
      track.px = BBtr_px[0];
      track.p = BBtr_p[0];
      track.n = (Int_t)BBtr_n;

      //Fill HCal variables from primary cluster
      hcal.x = HCAL_x;
      hcal.y = HCAL_y;
      hcal.e = HCAL_e;
      hcal.nblk = (Int_t)HCAL_nblk;
      hcal.nclus = (Int_t)HCAL_nclus;
      for( Int_t b=0; b<HCAL_nblk; b++ ){	
	Int_t ID = HCAL_cblk_id[b];
	
	hcal.cblk_e[ ID ] = HCAL_cblk_e[ b ];
	hcal.cblk_tdc[ ID ] = HCAL_cblk_tdctime[ b ];
	hcal.cblk_adct[ ID ] = HCAL_cblk_atime[ b ];
      }
      for( Int_t c=0; c<HCAL_nclus; c++ ){
	hcal.catime[ c ] = HCAL_catime[ c ];
	hcal.ce[ c ] = HCAL_ce[ c ];
	hcal.ceblk[ c ] = HCAL_ceblk[ c ];
	hcal.cid[ c ] = HCAL_cid[ c ];
	hcal.cnblk[ c ] = HCAL_cnblk[ c ];
	hcal.ctdc[ c ] = HCAL_ctdc[ c ];
	hcal.cx[ c ] = HCAL_cx[ c ];
	hcal.cy[ c ] = HCAL_cy[ c ];
      }
	
      //Fill BBCal PS
      bbcalps.x = BBps_x;
      bbcalps.y = BBps_y;
      bbcalps.e = BBps_e;
      bbcalps.nblk = (Int_t)BBps_nblk;
      bbcalps.nclus = (Int_t)BBps_nclus;
      for( Int_t b=0; b<BBps_nblk; b++ ){	
	Int_t ID = BBps_cblk_id[b];
	
	bbcalps.cblk_e[ ID ] = BBps_cblk_e[ b ];
	bbcalps.cblk_ec[ ID ] = BBps_cblk_ec[ b ];
	bbcalps.cblk_atime[ ID ] = BBps_cblk_atime[ b ];
      }
      //Fill BBCal Sh
      bbcalsh.x = BBsh_x;
      bbcalsh.y = BBsh_y;
      bbcalsh.e = BBsh_e;
      bbcalsh.nblk = (Int_t)BBsh_nblk;
      bbcalsh.nclus = (Int_t)BBsh_nclus;
      for( Int_t b=0; b<BBsh_nblk; b++ ){	
	Int_t ID = BBsh_cblk_id[b];
	
	bbcalsh.cblk_e[ ID ] = BBsh_cblk_e[ b ];
	bbcalsh.cblk_ec[ ID ] = BBsh_cblk_ec[ b ];
	bbcalsh.cblk_atime[ ID ] = BBsh_cblk_atime[ b ];
      }

      elasYield++;
      AT->Fill();

    }
  }

  cout << endl;
  
  cout << "Total events counted on T tree: " << testevent << "." << endl;

  ofstream logfile;
  logfile.open( logpath );
  logfile << "#Log of runs corresponding to root file with same date" << endl;
  for( Int_t i=0; i<log.size(); i++ ){
    logfile << log[i] << endl;
  }

  logfile.close();

  cout << endl << endl << "Elastic yield for analyzed runs: " << elasYield << ". Total events analyzed: " << Nevents << "." << endl << endl;
  cout << "New root file generated with elastic cuts placed at " << outputfilename << endl;

  fout->Write();

  st->Stop();
  
  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;
  

  return;
}
