//SSeeds 2.1.23 - Pass2 Update - Calibration code which employs best current cuts on elastic events to obtain ADC gain calibration parameters (pC/GeV). Adds beam energy loss to target and structure updates to easily run on batch farm.

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

//Detector constants
const Int_t kNcell = 288; // Total number of HCal modules
const Int_t kNrows = 24; // Total number of HCal rows
const Int_t kNcols = 12; // Total number of HCal columns
const Int_t maxTracks = 1000; // Reasonable limit on tracks to be stored per event
const Int_t maxTdcChan = 10; // Set to accomodate original 5 TDCTrig channels with buffer
const Int_t uni_N = 400; // Total number of bins used to measure detection uniformity (hSampFrac histos)
const Int_t xN = 48; //2*kNrows, total number of dispersive bins detection uni
const Int_t yN = 24; //2*kNcols, total number of transverse bins detection uni
//Double_t hcalheight = 0.365; //m The height of the center of HCAL above beam
const Double_t hcalheight = -0.2897;
const Double_t sampFrac = 0.0795; //HCal sampling frac (0.06588 GeV/0.8286 GeV) = 0.0795 = 7.95% -> (MC E_dep per proton) / (fit to data KE_p)
//Target constants
const double dEdx_tgt=0.00574; //According to NIST ESTAR, the collisional stopping power of hydrogen is about 5.74 MeV*cm2/g at 2 GeV energy
const double dEdx_Al = 0.0021; //According to NIST ESTAR, the collisional stopping power of Aluminum is about 2.1 MeV*cm2/g between 1-4 GeV
const double uwallthick_LH2 = 0.0145; //cm
const double dwallthick_LH2 = 0.015; //cm
const double cellthick_LH2 = 0.02; //cm, this is a guess;
const double Alshieldthick = 2.54/8.0; //= 1/8 inch * 2.54 cm/inch  
const double dxlim_l = 4.0;
const double dxlim_u = 3.0;
const double l_tgt = 0.15; // Length of the target (m)
const double rho_tgt = 0.0723; // Density of target (g/cc)
const double rho_Al = 2.7; // Density of aluminum windows (g/cc)
const double celldiameter = 1.6*2.54; //cm, right now this is a guess
const double Ztgt = 1.0;
const double Atgt = 1.0;
const double Mmol_tgt = 1.008; //g/mol
//Math
const Double_t PI = TMath::Pi();
const Double_t M_e = 0.00051;
const Double_t M_p = 0.938272;


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

void ecal_hcal( const char *configfilename = "secal_hcal.cfg" ){
//void ecal_hcal( Int_t kine=9, Int_t iter=1 ){
  
  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  // Get the date
  string date = getDate();
  
  // Declare Chain for many root files
  TChain *C = new TChain("T");

  // Declare file name paths
  //const char *configfilename = Form("../config/SBS%d/secal_sbs%d.cfg",kine,kine);
  string constPath = "/w/halla-scshelf2102/sbs/seeds/HCal_replay/hcal/energyCalibration/outfiles/newRatio.txt";
  string ratioPath = "/w/halla-scshelf2102/sbs/seeds/HCal_replay/hcal/energyCalibration/outfiles/ratio.txt";
  string oneBlockPath = "/w/halla-scshelf2102/sbs/seeds/HCal_replay/hcal/energyCalibration/outfiles/oneBlock.txt";

  // Declare general physics/fit parameters
  Int_t minEventPerCell = 500; // Minimum number of scattered p in cell required to calibrate
  Int_t maxEventPerCell = 4000; // Maximum number of scattered p events to contribute
  Double_t HCal_divx = 0.15875; // Transverse width in x and y per cell
  Double_t HCal_divy = 0.15494;
  Double_t HCal_Xi = -2.165; // Distance from beam center to top of HCal in cad (m)
  Double_t HCal_Xf = 1.435; // Distance from beam center to bottom of HCal in cad (m)
  Double_t HCal_Yi = -0.9; // Distance from beam center to opposite-beam side of HCal in cad (m)
  Double_t HCal_Yf = 0.9; // Distance from beam center to beam side of HCal in cad (m)
  Double_t highDelta = 0.1; // Minimum M(i,j)/b(i) factor allowed 
  Double_t oldGain[kNcell] = {0.0};
  Double_t oldRatio[kNcell] = {0.0};
  Double_t maxDiff = 50.;

  // Declare general physics parameters to be modified by input config file
  Int_t kine = -1000; // Keep track of kinematic calibrating from
  Double_t cointime = -1000.; // Expected coincidence time difference in 1190 TDC between BBCal HI and HCal trigger
  Double_t E_e = -1000.; // Energy of beam (incoming electrons from accelerator)
  Double_t HCal_d = -1000.; // Distance to HCal from scattering chamber for comm1
  Double_t HCal_th = -1000.; // Angle that the center of HCal is at
  Double_t BB_th = -1000.; // Angle that the center of BBCal is at
  Double_t W2_mean = -1000.; // Mean of W at current kinematic
  Double_t W2_sig = -1000.; // Width of W at current kinematic
  Double_t dx0_n = -1000.; // Position of neutron spot, x-x_expected
  Double_t dx0_p = -1000.; // Position of proton spot, x-x_expected
  Double_t dy0 = -1000.; // Position of hadron spot, y-y_expected
  Double_t dx_sig_n = -1000.; // Max spread of neutron spot, x-x_expected
  Double_t dx_sig_p = -1000.; // Max spread of proton spot, x-x_expected
  Double_t dy_sig = -1000.; // Max spread of hadron spot, y-y_expected
  Double_t atime0 = -1000.; // Expected location in ADC time of signal
  Double_t atime_sig = -1000.; // 1 sig of atime distribution
  Int_t useAlshield = -1000;

  Int_t elasYield = 0; // Keep track of total elastics analyzed
  Double_t cluspC_min = 90.; // Check the minimum value of pC measurements
  Double_t cluspC_max = 0.; // Check the maximum value of pC measurements
  Int_t badtimeblk = 0; // Keep track of cluster blks out of time with primary blk
  Int_t badtimeblkclus = 0; // Keep track of multi-blk clusters with out of time blks
  Int_t multblkclus = 0; // Keep track of total multi-block clusters (nblk>1)

  //For position reconstruction
  Double_t HCal_Xmin = HCal_Xi-HCal_divx/2;
  Double_t HCal_Xmax = HCal_Xf+HCal_divx/2;
  Double_t HCal_Ymin = HCal_Yi-HCal_divy/2;
  Double_t HCal_Ymax = HCal_Yf-HCal_divy/2;

  // Reading config file
  ifstream configfile(configfilename);
  TString currentline;
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      C->Add(currentline);
      cout << "Loading file: " << currentline << ".." << endl;
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
      if( skey == "kine" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	kine = sval.Atoi();
	cout << "Loading kinematic setting: " << kine << endl;
      }
      if( skey == "cointime" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	cointime = sval.Atof();
	cout << "Loading L1A/HCal coincidence time: " << cointime << endl;
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
	HCal_th = sval.Atof() * TMath::DegToRad();	
	cout << "Loading HCal angle: " << HCal_th << endl;
      }
      if( skey == "BB_th" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	BB_th = sval.Atof() * TMath::DegToRad();	
	cout << "Loading BBCal angle: " << BB_th << endl;
      }
      if( skey == "W2_mean" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        W2_mean = sval.Atof();
	cout << "Loading W2 mean cut: " << W2_mean << endl;
      }
      if( skey == "W2_sig" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        W2_sig = sval.Atof();
	cout << "Loading W2 sigma cut: " << W2_sig << endl;
      }
      if( skey == "dx0_n" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	dx0_n = sval.Atof();
	cout << "Loading x position of neutron spot: " << dx0_n << endl;
      }
      if( skey == "dx0_p" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	dx0_p = sval.Atof();
	cout << "Loading y position of proton spot: " << dx0_p << endl;
      }
      if( skey == "dy0" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	dy0 = sval.Atof();
	cout << "Loading y position of both hadron spots: " << dy0 << endl;
      }
      if( skey == "dx_sig_n" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	dx_sig_n = sval.Atof();
	cout << "Loading x sigma of neutron spot: " << dx_sig_n << endl;
      }
      if( skey == "dx_sig_p" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	dx_sig_p = sval.Atof();
	cout << "Loading x sigma of proton spot: " << dx_sig_p << endl;
      }
      if( skey == "dy_sig" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	dy_sig = sval.Atof();
	cout << "Loading y sigma of both hadron spots: " << dy_sig << endl;
      }
      if( skey == "atime0" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	atime0 = sval.Atof();
	cout << "Loading ADC time mean: " << atime0 << endl;
      }
      if( skey == "atime_sig" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	atime_sig = sval.Atof();
	cout << "Loading ADC time sigma: " << atime_sig << endl;
      }
      if( skey == "useAlshield" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	useAlshield = sval.Atoi();
	cout << "Loading Aluminum absorber option: " << useAlshield << endl;
      }
    }
    delete tokens;
  }

  if( kine == -1000 ||
      cointime == -1000 ||
      E_e == -1000 || 
      HCal_d == -1000 ||
      HCal_th == -1000 ||
      BB_th == -1000 ||
      W2_mean == -1000 ||
      W2_sig == -1000 ||
      dx0_n == -1000 ||
      dx0_p == -1000 ||
      dy0 == -1000 ||
      dx_sig_n == -1000 ||
      dx_sig_p == -1000 ||
      dy_sig == -1000 ||
      atime0 == -1000 ||
      atime_sig == -1000 ||
      useAlshield == -1000 ){
    cout << "Error: Setup parameters not fully loaded. Check config file." << endl;
    return;
  }

  // Path for previous ADC gain constants
  string inConstPath;
  Int_t test = 1; //TODO: will remove on batch farm reconfigure for multiple iterations
  if( test==0 ){
    inConstPath = "/w/halla-scshelf2102/sbs/SBS_REPLAY/SBS-replay/DB/db_sbs.hcal.dat";
  }else if( test==1 ){
    inConstPath = "/w/halla-scshelf2102/sbs/seeds/SBS-replay/DB/db_sbs.hcal.dat";
  }else{
    cout << "ERROR: No database file where specified" << endl;
    return 0;
  }
  
  // Reading ADC gain parameters from database
  Double_t gOldConst[kNcell];
  
  cout << "Loading previous gain coefficients from file: " << inConstPath << ".." << endl;
  ifstream inConstFile( inConstPath );
  if( !inConstFile ){
    cerr << endl << "ERROR: No input constant file present -> path to db_sbs.hcal.dat expected." << endl;
    return 0;
  }

  Int_t n1=0;
  Double_t d1=0;
  string line;
  bool skip_line = true;
  bool skip_one_line = true;
  bool skip_first_instance = true;
  
  while( getline( inConstFile, line ) ){

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
	gOldConst[n1] = d1;
	n1++;
      }
    }
  }
  
  cout << endl << endl << "Old ADC gain params: " << endl;

  for( Int_t r=0; r<kNrows; r++){
    for( Int_t c=0; c<kNcols; c++){
      Int_t i = r*kNcols+c;
      cout << gOldConst[i] << " ";
    }
    cout << endl;
  }
  
  cout << endl << endl << "Setup parameters loaded." << endl;
  
  TEventList *elist = new TEventList("elist","Elastic Event List");
  C->Draw(">>elist",globalcut);
  
  cout << "Event list populated with cut placed on elastics." << endl;

  // Declare general detector parameters
  Double_t BBtr_p[maxTracks], BBtr_px[maxTracks], BBtr_py[maxTracks], BBtr_pz[maxTracks];
  Double_t BBtr_vz[maxTracks], BBtr_chi2[maxTracks];
  Double_t BBtr_n, BBps_x, BBps_y, BBps_e, BBsh_x, BBsh_y, BBsh_e;

  Double_t TDCT_id[maxTdcChan], TDCT_tdc[maxTdcChan]; 
  Int_t TDCTndata;

  Double_t HCALx, HCALy, HCALe;
  Double_t crow, ccol, nblk;
  Double_t cblkid[kNcell], cblke[kNcell], cblkatime[kNcell];
  Double_t ekineW2;

  // Declare root tree variables and set values to memory locations in root file
  C->SetBranchStatus( "*", 0 );
  C->SetBranchStatus( "sbs.hcal.x", 1 );
  C->SetBranchStatus( "sbs.hcal.y", 1 );
  C->SetBranchStatus( "sbs.hcal.e", 1 );
  C->SetBranchStatus( "sbs.hcal.rowblk", 1 );
  C->SetBranchStatus( "sbs.hcal.colblk", 1 );
  C->SetBranchStatus( "sbs.hcal.nblk", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.id", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.e", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.atime", 1 );
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
  C->SetBranchStatus( "e.kine.W2", 1 );

  C->SetBranchAddress( "sbs.hcal.x", &HCALx );
  C->SetBranchAddress( "sbs.hcal.y", &HCALy );
  C->SetBranchAddress( "sbs.hcal.e", &HCALe );
  C->SetBranchAddress( "sbs.hcal.rowblk", &crow );
  C->SetBranchAddress( "sbs.hcal.colblk", &ccol );
  C->SetBranchAddress( "sbs.hcal.nblk", &nblk ); // Total number of blocks in highest E clus
  C->SetBranchAddress( "sbs.hcal.clus_blk.id", cblkid ); // kNcell-1 index for each block
  C->SetBranchAddress( "sbs.hcal.clus_blk.e", cblke ); // Array of block energies
  C->SetBranchAddress( "sbs.hcal.clus_blk.atime", cblkatime ); // Array of block ADC times
  C->SetBranchAddress( "bb.tr.n", &BBtr_n );
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
  C->SetBranchAddress( "e.kine.W2", &ekineW2 );

  // Create stopwatch to track processing time
  TStopwatch *sw = new TStopwatch();
  sw->Start();
  
  // Declare outfile
  TFile *fout = new TFile( "eCalEOut.root", "RECREATE" );
  
  // Initialize vectors and arrays
  Double_t gRatio[kNcell] = {0.0}; //Old ratios from first iteration
  Double_t oneBlock[kNcell] = {0.0};
  Double_t GCoeff[kNcell] = {0.0};
  Double_t GCoeff_oneblock[kNcell] = {0.0};
  Double_t GCoeff_divide[kNcell] = {0.0};

  // Initialize histograms
  TH1D *hatime = new TH1D( "hatime", "Cluster adc time, primary block; ns", 220, -50, 170);
  TH1D *hatime_cut = new TH1D( "hatime_cut", "Cluster adc time, primary block, W2 cut; ns", 220, -50, 170);
  TH1D *hatime_diffallblks = new TH1D( "hatime_diffallblks", "Cluster adc time, blk - primary block; ns", 200, -100, 100);
  TH1D *hpC = new TH1D( "hpC", "Cluster element pC", 2000,0,2000);
  TH1D *hDeltaE = new TH1D( "hDeltaE","1.0-Eclus/p_rec", 100, -1.5, 1.5 );
  TH1D *hHCALe = new TH1D( "hHCALe","HCal Cluster E", 400, 0., 1. );
  TH1D *hHCALe_cut = new TH1D( "hHCALe_cut","HCal Cluster E, W2 cut", 400, 0., 1. );
  TH1D *hSampFrac = new TH1D( "hSampFrac","W2 Cut HCal Cluster E / Expected KE", 400, 0., 1. );
  TH1D *hDiff = new TH1D( "hDiff","HCal time - BBCal time (ns)", 1300, -500, 800 );
  TH1D *hEDiffChan = new TH1D( "hEDiffChan","EDiff Events Per Channel", 288, 0, 288 );
  TH1D *hClusE_calc = new TH1D( "hClusE_calc","Best Cluster Energy Calculated from HCal Vars", 100, 0.0, 2.0);
  TH2D *hClusBlk_ClusE = new TH2D( "hClusBlk_ClusE","Cluster Size vs Recon Cluster E", 8, 0, 8, 100, 2.0, 3.0 );

  TH1D *hmclusblkE = new TH1D( "hmclusblkE","HCal multi-blk cluster blk E", 400, 0., 1. );
  TH1D *hmclusblkE_tcut = new TH1D( "hmclusblkE_tcut","HCal multi-blk cluster blk E, atime cut", 400, 0., 1. );
  TH1D *hsclusblkE = new TH1D( "hsclusblkE","HCal single-blk cluster blk E", 400, 0., 1. );

  TH2D *hPAngleCorr = new TH2D( "hPAngCorr","Track p vs Track ang", 100, 30, 60, 100, 0.4, 1.2 );
  TH2D *hPAngleCorr_2 = new TH2D( "hPAngCorr_2","Track p vs Track ang v2", 100, 30, 60, 100, 0.4, 1.2 );
  TH2D *hClusE_vs_X = new TH2D("hClusE_vs_X",";X-pos (m);E_dep (GeV)",500,-3,2,100,0.0,1.0);  
  TH2D *hClusE_vs_Y = new TH2D("hClusE_vs_Y",";Y-pos (m);E_dep (GeV)",200,-1,1,100,0.0,1.0);    
  TH1D *hpp = new TH1D( "hpp", "Elastic Proton Momentum", 600, 0, 6 );
  TH1D *hdx = new TH1D( "hdx","; x_{HCAL} - x_{exp} (m)", 250, -dxlim_l,dxlim_u);
  TH1D *hdx_coincut = new TH1D("hdx_coincut","delta x with coincidence time cut; x_{HCAL} - x_{exp} (m)", 250, -dxlim_l,dxlim_u);  
  TH1D *hdx_dycut = new TH1D("hdx_dycut","delta x with delta y cut; x_{HCAL} - x_{exp} (m)", 250, -dxlim_l,dxlim_u);
  TH1D *hdx_atimecut = new TH1D("hdx_atimecut","delta x with atime cut; x_{HCAL} - x_{exp} (m)", 250, -dxlim_l,dxlim_u);
  TH1D *hdy = new TH1D("hdy","; y_{HCAL} - y_{exp} (m);", 250, -1.25,1.25);
  TH2D *hdxdy = new TH2D("hdxdy",";y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)", 250, -1.25, 1.25, 250, -dxlim_l, dxlim_u );
  TH1D *hW2_calc = new TH1D( "W2_calc", "W2 Calculated; GeV", 250, 0.3, 1.5 );
  TH1D *hW2_tree = new TH1D( "W2_tree", "W2 From Tree; GeV", 250, 0.3, 1.5 );
  TH1D *hW2_coincut = new TH1D( "W2_coincut", "W2 Coin Cut; GeV", 250, 0.3, 1.5 );
  TH1D *hNBlk = new TH1D( "hNBlk", "Number of Blocks in Primary Cluster", 25, 0, 25 );
  TH1D *hNBlk_cut = new TH1D( "hNBlk_cut", "Number of Blocks in Primary Cluster W2 Cut", 25, 0, 25 );
  TH1D *hQ2 = new TH1D( "Q2", "Q2; GeV", 250, 0.5, 3.0 );
  TH1D *hE_ep = new TH1D( "hE_ep","Scattered Electron Energy; GeV", 500, 0.0, E_e*1.5 ); 
  TH1D *hE_pp = new TH1D( "hE_pp", "Scattered Proton Energy; GeV", 500, 0.0, E_e*1.5 );
  TH1D *hKE_p = new TH1D( "hKE_p", "Scattered Proton Kinetic Energy; GeV", 500, 0.0, E_e*1.5 );
  TH2D *hADC = new TH2D( "hADC", "HCal Int_ADC Spectra: W2 Cut", 288, 0, 288, 100., 0., 1. );
  TH2D *coefficients = new TH2D( "coefficients",";Channel ;GeV/mV", 288,0,288,250,0,0.025 );
  TH2D *coefficients_1b = new TH2D( "coefficients_1b",";Channel ;GeV/mV", 288,0,288,250,0,0.025 );
  TH2D *hEDiff_vs_X = new TH2D( "hEDiff_vs_X",";X-pos (m);(E_exp-E_dep)/E_exp (GeV)",xN,HCal_Xi,HCal_Xf,100,-3,3 ); 
  TH2D *hEDiff_vs_Y = new TH2D( "hEDiff_vs_Y",";Y-pos (m);(E_exp-E_dep)/E_exp (GeV)",yN,HCal_Yi,HCal_Yf,100,-3,3 );
  TH2D *hEDiff_vs_block = new TH2D( "hEDiff_vs_block",";Block (cell);(E_exp-E_dep)/E_exp (GeV)",288,0,288,100,-3,3 );
  TH2D *hSampFrac_vs_X = new TH2D( "hSampFrac_vs_X","Sampling Fraction ;X (m) ;HCal_E / Exp_KE", xN, HCal_Xi, HCal_Xf, uni_N, 0., 1. );
  TH2D *hSampFrac_vs_Y = new TH2D( "hSampFrac_vs_Y","Sampling Fraction ;Y (m) ;HCal_E / Exp_KE", yN, HCal_Yi, HCal_Yf, uni_N, 0., 1. );
  TH1D *hE_pp_cut = new TH1D( "hE_pp_cut", "Deposited Elastic Proton Energy; GeV", 500, 0.0, E_e*1.5 );

  // Set long Int_t to keep track of total entries
  Long64_t Nevents = elist->GetN();

  cout << endl << "All parameters loaded." << endl << endl;
  cout << "Opened up tree with nentries: " << Nevents << ".." << endl << endl;

  //Declare matrices for chi-square min calibration scheme and keep track of calibrated events
  TMatrixD Ma(kNcell,kNcell);
  TMatrixD Ma_err(kNcell,kNcell);
  TVectorD ba(kNcell);
  TVectorD bb(kNcell);
  TVectorD ba_err(kNcell);
  Int_t NEV[kNcell] = {0};
  Int_t NEV_oneblock[kNcell] = {0};
  Double_t err[kNcell] = {0.};
  Double_t err_ev[kNcell] = {0.};
  Double_t err_oneblock[kNcell] = {0.};
  Double_t err_ev_oneblock[kNcell] = {0.};
  
  //Declare energy loss parameters for beam going through the target
  Double_t pBeam = E_e/(1.+E_e/M_p*(1.-cos(BB_th)));
  Double_t Eloss_outgoing = celldiameter/2.0/sin(BB_th) * rho_tgt * dEdx_tgt; //Mean energy loss of the beam prior to the scattering, approximately 1 MeV, could correct further with raster position (likely negligable)
  if( useAlshield != 0 ) Eloss_outgoing += Alshieldthick * rho_Al * dEdx_Al;

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
	timeremains = timekeeper*( double(Nevents)/double(nevent) - 1. ); 
      sw->Reset();
      sw->Continue();
      
      progress = (double)((nevent+1.)/Nevents);
      cout << "] " << int(progress*100.) << "%, elastic events: " << elasYield << ", time remaining: " << int(timeremains/60.) << "m \r";
      cout.flush();
      
      C->GetEntry( elist->GetEntry( nevent ) ); 

      Double_t A[kNcell] = {0.0}; // Array to keep track of ADC values per cell. Outscope on each ev
      
      //////////////////////////////////////////////////////
      //Cut on timing first to avoid unnecessary calculation
      //////////////////////////////////////////////////////
      Double_t bbcal_time=0., hcal_time=0.;
      for(Int_t ihit=0; ihit<TDCTndata; ihit++){
	if(TDCT_id[ihit]==5) bbcal_time=TDCT_tdc[ihit];
	if(TDCT_id[ihit]==0) hcal_time=TDCT_tdc[ihit];
      }
      Double_t diff = hcal_time - bbcal_time; 
      hDiff->Fill( diff ); // Fill histogram
      //if( fabs(diff-cointime)>maxDiff ) continue; //Do not include initially
      //if( diff>-355 && diff<490 ) continue;
      //if( diff>570 ) continue;
      
      //Correct the beam energy with energy loss in target using vertex position
      Double_t Eloss = (BBtr_vz[0]+l_tgt/2.0) * rho_tgt * dEdx_tgt + uwallthick_LH2 * rho_Al * dEdx_Al; //approximately 3 MeV
      Double_t E_corr = E_e - Eloss;
      Double_t p_corr = BBtr_p[0] - Eloss_outgoing; //Neglecting the mass of e'

      Double_t etheta = acos( BBtr_pz[0]/BBtr_p[0] );
      Double_t ephi = atan2( BBtr_py[0], BBtr_px[0] );

      TVector3 vertex(0,0,BBtr_vz[0]); // z location of vertex in hall coordinates
      TLorentzVector Pbeam(0,0,E_corr,E_corr); //Mass of e negligable
      TLorentzVector kprime(BBtr_px[0],BBtr_py[0],BBtr_pz[0],BBtr_p[0]);
      TLorentzVector Ptarg(0,0,0,M_p); // assume proton for both LH2 and LD2 - can refine where useful with exclusive LH2 data set. Likely better to refine first with dxdy spot cuts on both protons and neutrons.

      TLorentzVector q = Pbeam - kprime;
      TLorentzVector PgammaN = Ptarg + q; //(-px, -py, ebeam - pz, Mp + ebeam - p)

      Double_t pel = E_corr/(1.+E_corr/M_p*(1.-cos(etheta)));
      Double_t nu = E_corr - BBtr_p[0];
      Double_t pp = sqrt(pow(nu,2)+2.*M_p*nu);
      Double_t phinucleon = ephi + PI; //assume coplanarity
      Double_t thetanucleon = acos( (E_corr - BBtr_pz[0])/pp ); //use elastic constraint on nucleon kinematics
	
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
      Double_t Q2 = 2*E_corr*E_ep*( 1-(BBtr_pz[0]/p_ep) ); // Obtain Q2 from beam energy, outgoing electron energy, and momenta
	
      Double_t W = PgammaN.M();
      Double_t W2 = ekineW2;

      //Use the electron kinematics to predict the proton momentum assuming elastic scattering on free proton at rest:
      Double_t E_pp = nu+M_p; // Get energy of the proton
      Double_t Enucleon = sqrt(pow(pp,2)+pow(M_p,2)); // Check on E_pp, same
      Double_t KE_p = nu; //For elastics
      Double_t dx = HCALx - xexpect_HCAL;
      Double_t dy = HCALy - yexpect_HCAL;
      
      //Fill pre-elastic-cut histograms
      hE_ep->Fill( E_ep );
      hQ2->Fill( Q2 );
      hW2_calc->Fill( pow(W,2) );
      hW2_tree->Fill( W2 );
      hE_pp->Fill( E_pp );
      hKE_p->Fill( KE_p );   
      hHCALe->Fill( HCALe );
      hNBlk->Fill( nblk );
      hatime->Fill( cblkatime[0] );

      //Coincidence timing cut to check effect on W2
      if( fabs(diff-cointime)<40 ) hW2_coincut->Fill( W2 );

      /////////////////////////////////////////////////
      //Primary W2 cut on elastics and HCal active area
      /////////////////////////////////////////////////
      if( fabs(W2-W2_mean)>W2_sig ) continue; //Observed mean W2 cut on elastic peak
      if( crow==0 || crow==23 || ccol==0 || ccol==11 ) continue; //All events with primary cluster element on edge blocks cut
      //Events that pass the above cuts constitute elastics
      elasYield++;

      //Calculate/declare new variables for analysis
      Double_t SFrac = HCALe/KE_p;
      Double_t E_exp = KE_p*sampFrac;
      Int_t rec_row = ( HCALx - HCal_Xmin )/HCal_divx;
      Int_t rec_col = ( HCALy - HCal_Ymin )/HCal_divy;
      Int_t rec_cell = rec_row*kNcols + rec_col;

      //Fill delta histograms
      hpp->Fill( pp );
      hdx->Fill( dx );
      hdy->Fill( dy );
      hdxdy->Fill( dy, dx );
      hatime_cut->Fill( cblkatime[0] );

      //////////////////////////////////////////////////////////////////
      //Cut on dy and atime. Done here to preserve diagnostic histograms
      //////////////////////////////////////////////////////////////////
      if( abs(dy-dy0)<5*dy_sig ) continue;
      if( abs(cblkatime[0]-atime0)>3*atime_sig ) continue;

      /*
      //Does not improve dx. Not used.
      if( diff<-355 || (diff>490 && diff<570) ) {
	hdx_coincut->Fill( dx );
      }
      //Significantly improves dx
      if( abs(dy-dy0)<5*dy_sig ){
	hdx_dycut->Fill( dx );
      }
      //Significantly improves dx
      if( abs(cblkatime[0]-atime0)<3*atime_sig ){
	hdx_atimecut->Fill( dx );
      }
      */

      //Fill histograms after elastics cuts
      hHCALe_cut->Fill( HCALe );
      hSampFrac->Fill( SFrac );
      hSampFrac_vs_X->Fill( HCALx, SFrac );
      hSampFrac_vs_Y->Fill( HCALy, SFrac );
      hNBlk_cut->Fill( nblk );
      hE_pp_cut->Fill( E_exp );
      hPAngleCorr->Fill( ephi, p_ep );
      hPAngleCorr_2->Fill( etheta, p_ep );
      hEDiffChan->Fill( rec_cell );

      // Get energies with simplest scheme from clusters only
      Double_t clusE = 0.0;
      Double_t cluspC = 0.0;
      bool badclus = false;
      for( Int_t blk = 0; blk<(int)nblk; blk++ ){
	Int_t blkid = int(cblkid[blk])-1; //-1 necessary since sbs.hcal.clus_blk.id ranges from 1 to 288 (different than just about everything else)
  
	bool badblk = false;
	if( nblk>1 && blk!=0 ) hatime_diffallblks->Fill(cblkatime[blk]-cblkatime[0]);
	if( abs(cblkatime[blk]-cblkatime[0])>3*atime_sig ){
	  badtimeblk++; //Should be the same cut as primary on timing
	  badclus = true;
	  badblk = true;
	}

	if(nblk==1) hsclusblkE->Fill(cblke[blk]);
	if(nblk>1) hmclusblkE->Fill(cblke[blk]);
	if(nblk>1&&!badblk) hmclusblkE_tcut->Fill(cblke[blk]);

	if( abs(cblkatime[blk]-cblkatime[0])>3*atime_sig ) continue; //Cut blocks that are out of time

	clusE += cblke[blk];
	A[blkid] += cblke[blk];
  	  
	//Error per event (2.1.23 estimate)
	cluspC = cblke[blk]/gOldConst[blkid];

	//pC diagnostics
	if( cluspC<cluspC_min ) cluspC_min = cluspC;
	if( cluspC>cluspC_max ) cluspC_max = cluspC;

	hpC->Fill( cluspC );

	//Error per event is DpC/pC(ev) + DE/E(ev); DpC=0.01(RAU limit)+0.5*pC(HCal Eres); DE=0.05*E(BBCal Eres). Estimate - glossing over elements which have minimal impact such as calculation params of KE_p
	//err_ev[blkid] += pow( ((0.01+0.5*cluspC)/cluspC+(0.05*KE_p)/KE_p) , 2 ); //Getting ready for std. dev of mean
	err_ev[blkid] += pow( ((1.+0.5*cluspC)/cluspC+(0.05*KE_p)/KE_p) , 2 ); //Getting ready for std. dev of mean

	//cout << (0.01+0.5*cluspC)/cluspC << " " << (0.05*KE_p)/KE_p << endl;

	// Simple estimation of the coefficients assuming 100% energy deposition in one block. Will check against the more sophisticated version with chi^2 reduction.
	if(nblk==1) {
	  //err_ev_oneblock[blkid] += pow( ((0.01+0.5*cluspC)/cluspC+(0.05*KE_p)/KE_p) , 2 ); //Probably better approximation of the error since all the KE_p should be deposited in this block.
	  err_ev_oneblock[blkid] += pow( ((1.+0.5*cluspC)/cluspC+(0.05*KE_p)/KE_p) , 2 ); //Probably better approximation of the error since all the KE_p should be deposited in this block.
	  NEV_oneblock[blkid]++;
	}
	NEV[blkid]++;

	hADC->Fill( blkid, cblke[blk] );

      }
      if( badclus ) badtimeblkclus++;
      if( nblk>1 ) multblkclus++;

      //Calculate clusE variables
      Double_t E_dev = (E_exp-clusE)/E_exp; //Energy deviation as percent of expected
	
      //Fill some clusE histograms
      hDeltaE->Fill( 1.0-(clusE/KE_p) );
      hClusBlk_ClusE->Fill( nblk, clusE ); //As function of reconstructed energy
      hClusE_vs_X->Fill( HCALx, clusE );
      hClusE_vs_Y->Fill( HCALy, clusE );
      hClusE_calc->Fill( clusE );
      hEDiff_vs_X->Fill( HCALx, E_dev );
      hEDiff_vs_Y->Fill( HCALy, E_dev );
      hEDiff_vs_block->Fill( rec_cell, E_dev );      
      
      //Build the matrix as simply as possible
      for(Int_t icol = 0; icol<kNcell; icol++){
	ba(icol)+= A[icol];
	if(nblk==1){
	  //ba_err(icol)+=E_exp;
	  //bb(icol)+= A[icol];
	  bb(icol)+= E_exp;
	  //oneBlock[icol]+= A[icol]*A[icol]/E_exp;
	  oneBlock[icol]+= E_exp;
	}
	for(Int_t irow = 0; irow<kNcell; irow++){
	  //Ma(icol,irow) += A[icol]*A[irow]/(KE_p*sampFrac);
	  Ma(icol,irow) += A[icol]*A[irow]/E_exp;
	  //if(nblk==1) Ma_err(icol,irow) += E_exp;
	} 
      }
    }
  }

  //Tests on various HCal cuts
  /*
  TH1D *hdx_coincut_diff = (TH1D*)hdx->Clone();
  hdx_coincut_diff->SetName("hdx_coincut_diff");
  hdx_coincut_diff->Add(hdx_coincut,-1);
  TH1D *hdx_dycut_diff = (TH1D*)hdx->Clone();
  hdx_dycut_diff->SetName("hdx_dycut_diff");
  hdx_dycut_diff->Add(hdx_dycut,-1);
  TH1D *hdx_atimecut_diff = (TH1D*)hdx->Clone();
  hdx_atimecut_diff->SetName("hdx_atimecut_diff");
  hdx_atimecut_diff->Add(hdx_atimecut,-1);
  */

  cout << endl << "Checking data, inverting matrix, and solving for coefficients.." << endl << endl;

  //Reject the bad cells and normalize the oneblock check
  Int_t badcell[kNcell];
  Int_t badcell_oneblock[kNcell];
  Int_t cellBad = 0;
  Double_t y[kNcell] = {0.0}; // For easy TGraphErrors build
  
  for(Int_t i=0; i<kNcell; i++){
    badcell[i] = 0;
    y[i] = i;
    
    //Do not change ADC gain coeff if insufficient events or energy dep in cell
    if( NEV[i] < minEventPerCell || Ma(i,i) < 0.1*ba(i) ){ 

      cellBad = 1;

      Double_t elemRatio = Ma(i,i)/ba(i);

      ba(i) = 1.0;  // Set RHS vector for cell i to 1.0 
      Ma(i,i) = 1.0; // Set diagonal element of Matrix M for cell i to 1.0 
      
      for(Int_t j=0; j<kNcell; j++){
	if( j != i ){
	  
	  Ma(i,j) = 0.0; 
	  Ma(j,i) = 0.0;
	}
      }
      badcell[i] = 1;
      cout << "cell" << i << " bad" << endl;
      cout << "Number of events in bad cell =" << NEV[i] << endl;
      cout << "Matrix element/vector element ratio =" << elemRatio << endl;
    }  

    //Perform same rejection on single block analysis
    badcell_oneblock[i] = 0;
    if( NEV_oneblock[i] < minEventPerCell || oneBlock[i] < 0.1*bb[i] ){
      bb(i) = 1.0;
      oneBlock[i] = 1.0;
      badcell_oneblock[i] = 1;
    }

    //Calculate error per block (element of A)
    err[i] = sqrt(err_ev[i]/NEV[i]);
    err_oneblock[i] = sqrt(err_ev_oneblock[i]/NEV_oneblock[i]);
  }

  if( cellBad==0 ) cout << "No bad cells detected.." << endl << endl;

  //Perform same rejection on single block analysis
  /*
  for(Int_t i=0; i<kNcell; i++){
    badcell_oneblock[i] = 0;
    if( NEV_oneblock[i] < minEventPerCell || oneBlock[i] < 0.1*bb[i] ){
      bb(i) = 1.0;
      oneBlock[i] = 1.0;
      badcell_oneblock[i] = 1;
    }
  }
  */

  //Invert the matrix, solve for ratios
  TMatrixD M_inv = Ma.Invert();
  //TMatrixD M_inv_err = Ma_err.Invert();
  TVectorD Coeff = M_inv*ba; // Stays unmodified for reference
  Double_t oneBlockCoeff[kNcell];

  for( Int_t i=0; i<kNcell; i++ ){
    if(badcell_oneblock[i]==0){
      oneBlockCoeff[i] = gOldConst[i]*(bb[i]/oneBlock[i]);

    }else{
      oneBlockCoeff[i] = gOldConst[i];
    }
  }

  for(Int_t i=0; i<kNcell; i++){
    if(badcell[i]==0){
      GCoeff[i]=gOldConst[i]*Coeff[i]; // The new gain coefficient is the old coefficient multiplied by the gain factor that we just solved for. We will take these coefficents as the input for another set of data to estimate the error per channel.
      GCoeff_divide[i]=GCoeff[i]/oneBlockCoeff[i];

    }else{
      GCoeff[i]=gOldConst[i]; // If the cell is bad, use the old coefficient
      GCoeff_divide[i]=-1.0;
    }
  }

  cout << "Inversion complete. Building histograms and writing out coefficients.." << endl << endl;

  for(Int_t i=0; i<kNcell; i++){
    coefficients->Fill(i,GCoeff[i]);
    coefficients_1b->Fill(i,oneBlockCoeff[i]);
  }

  Double_t yErr[kNcell] = {0.0};
  //Double_t constErr[kNcell] = {0.0}; //Will need to improve error here

  TGraphErrors *ccgraph = new TGraphErrors( kNcell, y, GCoeff, yErr, err ); 
  ccgraph->GetXaxis()->SetLimits(0.0,288.0);  
  ccgraph->GetYaxis()->SetLimits(0.0,0.25);
  ccgraph->SetTitle("Calibration Coefficients");
  ccgraph->GetXaxis()->SetTitle("Channel");
  ccgraph->GetYaxis()->SetTitle("Unitless");
  ccgraph->SetMarkerStyle(20); // idx 20 Circles, idx 21 Boxes
  ccgraph->Write("constants");

  TGraphErrors *ccgraph_oneBlock = new TGraphErrors( kNcell, y, oneBlockCoeff, yErr, err_oneblock ); 
  ccgraph_oneBlock->GetXaxis()->SetLimits(0.0,288.0);  
  ccgraph_oneBlock->GetYaxis()->SetLimits(0.0,0.25);
  ccgraph_oneBlock->SetTitle("Calibration Coefficients One Block");
  ccgraph_oneBlock->GetXaxis()->SetTitle("Channel");
  ccgraph_oneBlock->GetYaxis()->SetTitle("Unitless");
  ccgraph_oneBlock->SetMarkerStyle(21);
  ccgraph_oneBlock->Write("constants_oneblock");

  TGraphErrors *ccgraph_divide = new TGraphErrors( kNcell, y, GCoeff_divide, yErr, yErr ); 
  ccgraph_divide->GetXaxis()->SetLimits(0.0,288.0);  
  ccgraph_divide->SetTitle("Calibration Coefficients / OneBlock Coeff");
  ccgraph_divide->GetXaxis()->SetTitle("Channel");
  ccgraph_divide->GetYaxis()->SetTitle("Unitless");
  ccgraph_divide->SetMarkerStyle(21);
  ccgraph_divide->Write("constants_divide");

  //Write out diagnostic histos and print to console
  //fout->Write();

  cout << "Gain Coefficients: " << endl << endl;

  Int_t cell = 0;

  for( Int_t r = 0; r<kNrows; r++){
    for( Int_t c = 0; c<kNcols; c++){
      if( GCoeff[cell]==1 ) GCoeff[cell] = 0.00175;
      cout << GCoeff[cell] << "  ";
      cell++;
    }
    cout << endl;
  }

  cell = 0;

  for( Int_t i=0; i<kNcell; i++){
    GCoeff[i] = GCoeff[i] - 0.00175; //Cosmic value in mV/GeV for all channels
  }

  TGraphErrors *ccgraph_Cdiff = new TGraphErrors( kNcell, y, GCoeff, yErr, yErr ); 
  ccgraph_Cdiff->GetXaxis()->SetLimits(0.0,288.0);  
  ccgraph_Cdiff->GetYaxis()->SetLimits(0.0,0.025);
  ccgraph_Cdiff->SetTitle("Cosmic Delta Calibration Coefficients");
  ccgraph_Cdiff->GetXaxis()->SetTitle("Channel");
  ccgraph_Cdiff->GetYaxis()->SetTitle("Unitless");
  ccgraph_Cdiff->SetMarkerStyle(20); // idx 20 Circles, idx 21 Boxes
  ccgraph_Cdiff->Write("constants_diff");
  
  //Construct graphs for uniformity check
  Double_t posErr[xN] = {0.};
  TF1 *f1;
  //For X (dispersive)
  Double_t X[xN];
  Double_t Xval[xN];
  Double_t Xerr[xN];
  TH1D *Xslice[xN];
  for( Int_t x=0; x<xN; x++ ){
    X[x] = HCal_Xi + x*HCal_divx;
    Xslice[x] = hSampFrac_vs_X->ProjectionY(Form("Xslice_%d",x+1),x+1,x+1);
    Xslice[x]->Fit("gaus","Q","",0.01,0.22);
    f1=Xslice[x]->GetFunction("gaus");
    if(Xslice[x]->GetEntries()>0){
      Xval[x] = f1->GetParameter(1);
      Xerr[x] = f1->GetParameter(2);
    }
  }
  TGraphErrors *csampFrac_X = new TGraphErrors( xN, X, Xval, posErr, Xerr );
  csampFrac_X->GetXaxis()->SetLimits(HCal_Xi-0.05,HCal_Xf+0.05);  
  csampFrac_X->GetYaxis()->SetLimits(0.0,1.0);
  csampFrac_X->SetTitle("Sampling Fraction - Dispersive X");
  csampFrac_X->GetXaxis()->SetTitle("X (m)");
  csampFrac_X->GetYaxis()->SetTitle("E_{HCAL}/KE_{exp}");
  csampFrac_X->SetMarkerStyle(20); // idx 20 Circles, idx 21 Boxes
  csampFrac_X->Write("csampFrac_X");

  //For Y (transverse)
  Double_t Y[yN];
  Double_t Yval[yN];
  Double_t Yerr[yN];
  TH1D *Yslice[yN];
  for( Int_t x=0; x<yN; x++ ){
    Y[x] = HCal_Yi + x*HCal_divy;
    Yslice[x] = hSampFrac_vs_Y->ProjectionY(Form("Yslice_%d",x+1),x+1,x+1);
    Yslice[x]->Fit("gaus","Q","",0.01,0.22);
    f1=Yslice[x]->GetFunction("gaus");
    if(Yslice[x]->GetEntries()>0){
      Yval[x] = f1->GetParameter(1);
      Yerr[x] = f1->GetParameter(2);
    }
  }
  TGraphErrors *csampFrac_Y = new TGraphErrors( yN, Y, Yval, posErr, Yerr );
  csampFrac_Y->GetXaxis()->SetLimits(HCal_Yi-0.05,HCal_Yf+0.05);  
  csampFrac_Y->GetYaxis()->SetLimits(0.0,1.0);
  csampFrac_Y->SetTitle("Sampling Fraction - Transverse Y");
  csampFrac_Y->GetXaxis()->SetTitle("Y (m)");
  csampFrac_Y->GetYaxis()->SetTitle("E_{HCAL}/KE_{exp}");
  csampFrac_Y->SetMarkerStyle(20); // idx 20 Circles, idx 21 Boxes
  csampFrac_Y->Write("csampFrac_Y");

  //Construct graphs for energy calibration check
  TF1 *f2;
  //For X (dispersive)
  for( Int_t x=0; x<xN; x++ ){
    X[x] = HCal_Xi + x*HCal_divx;
    Xslice[x] = hEDiff_vs_X->ProjectionY(Form("eXslice_%d",x+1),x+1,x+1);
    Xslice[x]->Fit("gaus","Q","",-1.0,1.0);
    f2=Xslice[x]->GetFunction("gaus");
    if(Xslice[x]->GetEntries()>0){
      Xval[x] = f2->GetParameter(1);
      Xerr[x] = f2->GetParameter(2);
    }else{
      Xval[x] = 0.;
      Xerr[x] = 0.;
      cout << "Warning: Not enough statistics at X-pos " << x << " for calibration check." << endl;
    }
  }
  TGraphErrors *cEDiff_X = new TGraphErrors( xN, X, Xval, posErr, Xerr );
  cEDiff_X->GetXaxis()->SetLimits(HCal_Xi-0.05,HCal_Xf+0.05);  
  cEDiff_X->GetYaxis()->SetLimits(0.0,1.0);
  cEDiff_X->SetTitle("Energy Calibration Check - Dispersive X");
  cEDiff_X->GetXaxis()->SetTitle("X (m)");
  cEDiff_X->GetYaxis()->SetTitle("KE_{exp}-E_{HCAL}/KE_{exp}");
  cEDiff_X->SetMarkerStyle(20); // idx 20 Circles, idx 21 Boxes
  cEDiff_X->Write("cEDiff_X");

  //For Y (transverse)
  for( Int_t x=0; x<yN; x++ ){
    Y[x] = HCal_Yi + x*HCal_divy;
    Yslice[x] = hEDiff_vs_Y->ProjectionY(Form("eYslice_%d",x+1),x+1,x+1);
    Yslice[x]->Fit("gaus","Q","",-1.,1.);
    f2=Yslice[x]->GetFunction("gaus");
    if(Yslice[x]->GetEntries()>0){
      Yval[x] = f2->GetParameter(1);
      Yerr[x] = f2->GetParameter(2);
    }else{
      Yval[x] = 0.;
      Yerr[x] = 0.;
      cout << "Warning: Not enough statistics at Y-pos " << x << " for calibration check." << endl;
    }
  }
  TGraphErrors *cEDiff_Y = new TGraphErrors( yN, Y, Yval, posErr, Yerr );
  cEDiff_Y->GetXaxis()->SetLimits(HCal_Yi-0.05,HCal_Yf+0.05);  
  cEDiff_Y->GetYaxis()->SetLimits(0.0,1.0);
  cEDiff_Y->SetTitle("Energy Calibration Check - Transverse Y");
  cEDiff_Y->GetXaxis()->SetTitle("Y (m)");
  cEDiff_Y->GetYaxis()->SetTitle("KE_{exp}-E_{HCAL}/KE_{exp}");
  cEDiff_Y->SetMarkerStyle(20); // idx 20 Circles, idx 21 Boxes
  cEDiff_Y->Write("cEDiff_Y");

  fout->Write();

  cout << endl << "Gain Coefficients diff from cosmic value: " << endl;

  for( Int_t r = 0; r<kNrows; r++){
    for( Int_t c = 0; c<kNcols; c++){
      cout << GCoeff[cell] << "  ";
      cell++;
    }
    cout << endl;
  }

  cell = 0;

  cout << endl << "Number of events available for calibration: " << endl << endl;

  for( Int_t r = 0; r<kNrows; r++){
    for( Int_t c = 0; c<kNcols; c++){
      cout << NEV[cell] << "  ";
      cell++;
    }
    cout << endl;
  }
  cell = 0;

  cout << endl << "One Block:" << endl;

  for( Int_t r = 0; r<kNrows; r++){
    for( Int_t c = 0; c<kNcols; c++){
      cout << oneBlockCoeff[cell] << "  ";
      cell++;
    }
    cout << endl;
  }

  cout << endl << endl;

  //Declare outfiles
  ofstream ratio;
  ofstream GainCoeff;
  ofstream GainCoeff_oneblock;

  //Write to outfiles
  ratio.open( ratioPath );
  ratio << "#HCal Ratios from SBS-" << kine << " obtained " << date.c_str() << endl;
  ratio << "#HCal_ratio = " << endl;

  for( Int_t i=0; i<kNcell; i++ ){   
    ratio << Coeff[i] << endl;
  }

  ratio.close();

  GainCoeff.open( constPath );
  GainCoeff << "#HCal gain coefficients from SBS-" << kine << " obtained " << date.c_str() << endl;
  GainCoeff << "#HCal_gainCoeff = " << endl;

  cell = 0;
  for( Int_t r=0; r<kNrows; r++ ){
    for( Int_t c=0; c<kNcols; c++ ){
      GainCoeff << GCoeff[cell] << "  ";
      cell++;
    }
    GainCoeff << endl;
  }

  GainCoeff.close();

  GainCoeff_oneblock.open( oneBlockPath );
  GainCoeff_oneblock << "#HCal single block gain coefficients from SBS-" << kine << " obtained " << date.c_str() << endl;
  GainCoeff_oneblock << "#HCal_gainCoeff_oneBlock = " << endl;

  cell = 0;
  for( Int_t r=0; r<kNrows; r++ ){
    for( Int_t c=0; c<kNcols; c++ ){
      GainCoeff_oneblock << oneBlockCoeff[cell] << "  ";
      cell++;
    }
    GainCoeff_oneblock << endl;
  }
  cell = 0;

  GainCoeff_oneblock.close();

  cout << endl << endl << "Cluster pC measurement observed MIN: " << cluspC_min << ", MAX: " << cluspC_max << endl;

  cout << endl << endl << "Total blocks out of time with primary block / Multi-block clusters: " << badtimeblkclus << "/" << multblkclus << endl;

  cout << endl << endl << "Elastic yield for analyzed runs: " << elasYield << ". Total events analyzed: " << Nevents << "." << endl << endl;
  

  cout << "Calibration complete and constants written to file." << endl;

  st->Stop();

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

}

//Sample setup file, comment with #
////////////////////////////////////////////////////////
//setup_HCal_Calibration.txt
////////////////////////////////////////////////////////
//#LH2 full field
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11263*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11265*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11267*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11268*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11252*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11251*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11250*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11249*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11248*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11247*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11246*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11242*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11244*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11245*
//#LH2 quarter field
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11301*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11302*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11303*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11304*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11305*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11306*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11307*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11308*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11309*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11310*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11311*
//#LH2 inverted quarter field
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11276*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11277*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11279*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11280*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11281*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11282*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11283*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11284*
//#LH2 Zero Field 
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11332*
//endlist
//endcut
//E_e 1.92
//HCal_d 14.5
//HCal_th 35.0
//opticsCorr 1.05
//W_mean 0.93
//W_sig 0.1
//ScaleFac 1.0
/////////////////////////////////////////////////////////
