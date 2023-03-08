//SSeeds 3.3.23 - Script to run over all data from a given kinematic and plot energy spectra per block with several cuts

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
//Double_t hcalheight = 0.365; //m The height of the center of HCAL above beam
const Double_t hcalheight = -0.2897;
//const Double_t sampFrac[6] = {0.0797,0.0812,0.0812,0.0824,0.0811,0.0817}; //HCal sampling fractions from MC by kinematic. Position dependence confirmed to be negligible. 
//const Double_t sampFrac[6] = {0.0997,0.1012,0.1012,0.1024,0.1011,0.1017}; //HCal sampling fractions tuned by hand 
//Target constants
const Double_t sampFrac[6] = {0.0641,0.0718,0.0734,0.0710,0.0686,0.0682}; //HCal sampling fraction from 3.3.23 study adding cluster constraint and fair comparison to MC
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

//Keep this around just in case it's needed later
Double_t Gfit(Double_t *x, Double_t *par){
  Double_t amp = par[0];
  Double_t offset = par[1];
  Double_t sigma = par[2];
  return amp*exp(-0.5*pow((x[0]-offset)/sigma,2.));
}

//Main. Diagnostic by kinematic. Loads all configs for both lh2 and ld2.
void ecal_blockE_spectra( Int_t kine=-1){

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

  //Get gain parameters from database for pass0 where sbs.hcal.again not on tree
  Double_t gOldConst_pass0[kNcell] = {0.};

  bool pass0 = false;

  cout << endl;

  if( kine==4 || kine==7 ){ //end if iteration 1 (quasi-replay)

    pass0=true;
    
    // Reading ADC gain parameters from database 
    string inConstPath = "/u/home/sbs-gmn/pass0/SBS_REPLAY/SBS-replay/DB/db_sbs.hcal.dat";
    
    cout << "Loading previous gain coefficients from file for pass0 kinematics: " << inConstPath << ".." << endl;
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
	  gOldConst_pass0[n1] = d1;
	  n1++;
	}
      }
    }
    
    cout << endl << endl << "Old ADC gain params from pass0 db file: " << endl;
    
    for( Int_t r=0; r<kNrows; r++){
      for( Int_t c=0; c<kNcols; c++){
	Int_t i = r*kNcols+c;
	cout << gOldConst_pass0[i] << " ";
      }
      cout << endl;
    }
    
    cout << endl << endl;
  }//end if pass0

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
  Double_t cblkid_h[nfset_lh2[kIdx]][kNcell], cblke_h[nfset_lh2[kIdx]][kNcell], cblkatime_h[nfset_lh2[kIdx]][kNcell];
  Double_t cblkagain_h[nfset_lh2[kIdx]][kNcell];
  Double_t ekineW2_h[nfset_lh2[kIdx]];

  // Declare general detector parameters - ld2
  Double_t BBtr_p_d[nfset_ld2[kIdx]][maxTracks], BBtr_px_d[nfset_ld2[kIdx]][maxTracks], BBtr_py_d[nfset_ld2[kIdx]][maxTracks], BBtr_pz_d[nfset_ld2[kIdx]][maxTracks];
  Double_t BBtr_vz_d[nfset_ld2[kIdx]][maxTracks], BBtr_chi2_d[nfset_ld2[kIdx]][maxTracks];
  Double_t BBtr_n_d[nfset_ld2[kIdx]], BBps_x_d[nfset_ld2[kIdx]], BBps_y_d[nfset_ld2[kIdx]], BBps_e_d[nfset_ld2[kIdx]], BBsh_x_d[nfset_ld2[kIdx]], BBsh_y_d[nfset_ld2[kIdx]], BBsh_e_d[nfset_ld2[kIdx]];

  Double_t HCALx_d[nfset_ld2[kIdx]], HCALy_d[nfset_ld2[kIdx]], HCALe_d[nfset_ld2[kIdx]];
  Double_t crow_d[nfset_ld2[kIdx]], ccol_d[nfset_ld2[kIdx]], nblk_d[nfset_ld2[kIdx]];
  Double_t cblkid_d[nfset_ld2[kIdx]][kNcell], cblke_d[nfset_ld2[kIdx]][kNcell], cblkatime_d[nfset_ld2[kIdx]][kNcell];
  Double_t cblkagain_d[nfset_ld2[kIdx]][kNcell];
  Double_t ekineW2_d[nfset_ld2[kIdx]];

  // Declare file name paths
  string configfilename_h[nfset_lh2[kIdx]];
  for( Int_t h=0; h<nfset_lh2[kIdx]; h++){
    Int_t field = fset_lh2[kIdx][h];
    if(field==-1) cout << "Error. configfilename array out of bounds." << endl;
    configfilename_h[h] = Form("../config/GMn/SBS%d/secal_lh2_sbs%d_f%d.cfg",kine,kine,fset_lh2[kIdx][h]);
  }
  string configfilename_d[nfset_ld2[kIdx]];
  for( Int_t d=0; d<nfset_ld2[kIdx]; d++){
    Int_t field = fset_ld2[kIdx][d];
    if(field==-1) cout << "Error. configfilename array out of bounds." << endl;
    configfilename_d[d] = Form("../config/GMn/SBS%d/secal_ld2_sbs%d_f%d.cfg",kine,kine,fset_ld2[kIdx][d]);
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
  Double_t highDelta = 0.1; // Minimum M(i,j)/b(i) factor allowed 

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
    Ch[f]->SetBranchStatus( "e.kine.W2", 1 );

    Ch[f]->SetBranchAddress( "sbs.hcal.x", &HCALx_h[f] );
    Ch[f]->SetBranchAddress( "sbs.hcal.y", &HCALy_h[f] );
    Ch[f]->SetBranchAddress( "sbs.hcal.e", &HCALe_h[f] );
    Ch[f]->SetBranchAddress( "sbs.hcal.rowblk", &crow_h[f] );
    Ch[f]->SetBranchAddress( "sbs.hcal.colblk", &ccol_h[f] );
    Ch[f]->SetBranchAddress( "sbs.hcal.nblk", &nblk_h[f] ); // Total number of blocks in highest E clus
    Ch[f]->SetBranchAddress( "sbs.hcal.clus_blk.id", cblkid_h[f] ); // kNcell-1 index for each block
    Ch[f]->SetBranchAddress( "sbs.hcal.clus_blk.e", cblke_h[f] ); // Array of block energies
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
    Cd[f]->SetBranchStatus( "e.kine.W2", 1 );

    Cd[f]->SetBranchAddress( "sbs.hcal.x", &HCALx_d[f] );
    Cd[f]->SetBranchAddress( "sbs.hcal.y", &HCALy_d[f] );
    Cd[f]->SetBranchAddress( "sbs.hcal.e", &HCALe_d[f] );
    Cd[f]->SetBranchAddress( "sbs.hcal.rowblk", &crow_d[f] );
    Cd[f]->SetBranchAddress( "sbs.hcal.colblk", &ccol_d[f] );
    Cd[f]->SetBranchAddress( "sbs.hcal.nblk", &nblk_d[f] ); // Total number of blocks in highest E clus
    Cd[f]->SetBranchAddress( "sbs.hcal.clus_blk.id", cblkid_d[f] ); // kNcell-1 index for each block
    Cd[f]->SetBranchAddress( "sbs.hcal.clus_blk.e", cblke_d[f] ); // Array of block energies
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
    Cd[f]->SetBranchAddress( "e.kine.W2", &ekineW2_d[f] );

  }//end config file for over ld2 field settings

  cout << endl << endl << "Setup parameters loaded and chains linked." << endl;

  // Create stopwatch to track processing time
  TStopwatch *sw = new TStopwatch();
  sw->Start();
  
  // Declare outfile
  TFile *fout;
  fout = new TFile( Form("output/SBS%d/blockEspec_sbs%d.root", kine, kine ), "RECREATE" );

  // Initialize vectors and arrays
  Double_t gOldConst[kNcell] = {0.0};

  //No histograms needed. Will have to replay data with new database params, then use diagnostic script to assess
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

  //Declare tracking variables
  Int_t TNEV_h = 0;
  Int_t TNEV_d = 0;
  Int_t NEV[kNcell] = {0};
  Double_t err[kNcell] = {0.};
  Int_t deadclus = 0;
  Double_t pCllim = 0.;
  Double_t pCulim = 150.;
  Int_t pCbins = 300;
  
  //Declare diagnostic histograms (as sparse as possible)
  TH2D *hdx_mag_h = new TH2D( "hdx_mag_h", "Delta X vs Field Setting (LH2); field (percent); x_{HCAL} - x_{exp} (m)", 21, 0, 105, 250, -4, 3 );
  TH2D *hdy_mag_h = new TH2D( "hdy_mag_h", "Delta Y vs Field Setting (LH2); field (percent); y_{HCAL} - y_{exp} (m)", 21, 0, 105, 250, -1.25, 1.25 );
  TH2D *hdx_mag_d = new TH2D( "hdx_mag_d", "Delta X vs Field Setting (LD2); field (percent); x_{HCAL} - x_{exp} (m)", 21, 0, 105, 250, -4, 3 );
  TH2D *hdy_mag_d = new TH2D( "hdy_mag_d", "Delta Y vs Field Setting (LD2); field (percent); y_{HCAL} - y_{exp} (m)", 21, 0, 105, 250, -1.25, 1.25 );
  TH1D *hHCALe = new TH1D( "hHCALe","HCal Cluster E", 400, 0., 1. );
  TH2D *hHCALeID = new TH2D( "hHCALeID","HCal Cluster E vs ID", 288, 0, 288, 400, 0., 1. );
  TH1D *hHCALblkse = new TH1D( "hHCALblke","HCal Cluster Seed E", 400, 0., 1. );
  TH2D *hHCALblkeID = new TH2D( "hHCALblkeID","HCal Cluster blk E vs ID", 288, 0, 288, 400, 0., 1. );
  TH1D *hSampFrac = new TH1D( "hSampFrac","W2 Cut HCal Cluster E / Expected KE", 400, 0., 0.3 );
  TH2D *hSampFracID = new TH2D( "hSampFracID","W2 Cut HCal Cluster E / Expected KE vs ID", 288, 0, 288, 400, 0., 0.5 );
  TH1D *hHCALblkpC = new TH1D( "hHCALblkpC","HCal Cluster Element pC", pCbins, pCllim, pCulim );
  TH1D *hHCALblkpC_1b = new TH1D( "hHCALblkpC_1b","HCal Cluster Element pC, Single Block Clusters", pCbins, pCllim, pCulim );
  TH2D *hHCALblkpCID_1b = new TH2D( "hHCALblkpCID_1b","HCal Cluster Element pC, Single Block Clusters vs Cell", 288, 0, 288, pCbins, pCllim, pCulim );
  TH1D *hHCALblkspC = new TH1D( "hHCALblkpC","HCal Cluster Seed pC", pCbins, pCllim, pCulim );
  TH2D *hHCALblkpCID = new TH2D( "hHCALblkpCID","HCal Cluster blk pC vs ID", 288, 0, 288, pCbins, pCllim, pCulim );
  TH2D *hHCALblkpCID_h = new TH2D( "hHCALblkpCID_h","HCal Cluster blk pC vs ID (LH2)", 288, 0, 288, pCbins, pCllim, pCulim );
  TH2D *hHCALblkpCID_d = new TH2D( "hHCALblkpCID_d","HCal Cluster blk pC vs ID (LD2)", 288, 0, 288, pCbins, pCllim, pCulim );
  TH1D *hblkid = new TH1D( "hblkid","HCal Cluster Block ID, ADC time cut", 400, -50, 350 );
  TH1D *hblkid_badatime = new TH1D( "hblkid_badatime","HCal Block ID, Bad Cluster ADCtime", 400, -50, 350 );
  TH1D *hblkc_badatime = new TH1D( "hblkc_badatime","HCal Block Col, Bad Cluster ADCtime", 400, -50, 350 );
  TH1D *hblkr_badatime = new TH1D( "hblkr_badatime","HCal Block Row, Bad Cluster ADCtime", 400, -50, 350 );
  TH1D *hHCALeDev = new TH1D( "hHCALeDev","HCal Cluster E deviation from E expected", 400, -1., 1. );

  TH1D *hsigEratio = new TH1D( "hsigEratio","Sigma to Energy ratio for single block clusters", 100, 0.5, 1.5 );

  //Loop over events
  cout << "Main loop over all data commencing.." << endl;
  if( pass0 ) gROOT->ProcessLine( "gErrorIgnoreLevel = 6001" ); //Suppress error output to avoid undetectable sbs.hcal.clus_blk.again for pass0

  //Loop over all hydrogen data
  for( Int_t f=0; f<nfset_lh2[kIdx]; f++ ){
    
    Int_t hfieldset = fset_lh2[kIdx][f];

    //Declare energy loss parameters for beam going through the target
    Double_t pBeam = E_e_h[f]/(1.+E_e_h[f]/M_p*(1.-cos(BB_th_h[f])));
    Double_t Eloss_outgoing = celldiameter/2.0/sin(BB_th_h[f]) * rho_tgt * dEdx_tgt; //Mean energy loss of the beam prior to the scattering, approximately 1 MeV, could correct further with raster position (likely negligable)
    
    for( Long64_t nevent = 1; nevent <Nevents_h[f]; nevent++){

      if ( nevent%10000==0 ) cout << "LH2 kinematic " << kine << " at field " << hfieldset << "% , entry: " << nevent << "/" << Nevents_h[f] << ". Cut on atime killed clusters: " << deadclus << ". Total events gathered for analysis: " << TNEV_h << " \r";
      cout.flush();
      
      Ch[f]->GetEntry( elist_h[f]->GetEntry( nevent ) ); 

      //////////////////////////////////////////////////////////////////////
      //Make cuts that don't rely on projections first to improve efficiency
      if( crow_h[f]==0 || 
	  crow_h[f]==23 || 
	  ccol_h[f]==0 || 
	  ccol_h[f]==11 ) continue; //All events with primary cluster element on edge blocks cut
      if( abs(cblkatime_h[f][0]-atime0_h[f])>ifac*atime_sig_h[f] ) continue; //All events where adctime outside of reasonable window cut
      //////////////////////////////////////////////////////////////////////
      
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
      Double_t pp = sqrt(pow(nu,2)+2.*M_p*nu);
      Double_t phinucleon = ephi + PI; //assume coplanarity
      Double_t thetanucleon = acos( (E_corr - BBtr_pz_h[f][0])/pp ); //use elastic constraint on nucleon kinematics
	
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
      //Double_t KE_p = Enucleon; //For sanity check on sampling fraction (2.28.23)
      Double_t dx = HCALx_h[f] - xexpect_HCAL;
      Double_t dy = HCALy_h[f] - yexpect_HCAL;
      
      /////////////////////////////////////////////////
      //Primary W2 cut on elastics and HCal active area
      if( fabs(W2-W2_mean_h[f])>W2_sig_h[f] ) continue; //Observed mean W2 cut on elastic peak
      elasYield++; //Events that pass the above cuts constitute elastics
      /////////////////////////////////////////////////

      //Calculate/declare new variables for analysis
      Double_t SFrac = HCALe_h[f]/KE_p;
      Double_t E_exp = KE_p*sampFrac[kIdx];
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
      if( !pass_p && !pass_n ) continue; //Cut on both n and p spots for each event, cannot know which apriori
      //////////////////////////////////////////////////////////////////

      //Write out diagnostic histograms
      Double_t clusE = 0.0;
      Double_t clusblksE = cblke_h[f][0];
      Double_t clusblkspC = 0.0;
      Double_t cluspC = 0.0;
      Double_t blksatime = cblkatime_h[f][0];

      for( Int_t blk = 0; blk<(int)nblk_h[f]; blk++ ){
      	Int_t blkid = int(cblkid_h[f][blk])-1; //-1 necessary since sbs.hcal.clus_blk.id ranges from 1 to 288 (different than just about everything else)
	Double_t blkE = cblke_h[f][blk];
	Double_t blkatime = cblkatime_h[f][blk];
	Double_t blkagain = cblkagain_h[f][blk];
	Int_t blkr = blkid/kNcols;
	Int_t blkc = blkid%kNcols;

	//////////////////////////////
	//Time cut on cluster elements
	if( abs(blkatime-blksatime)>ifac*atime_sig_h[f] ){
	  hblkid_badatime->Fill(blkid);
	  hblkr_badatime->Fill(blkr);
	  hblkc_badatime->Fill(blkc);
	  continue;
	}
	//////////////////////////////

	if( pass0 ){
	  if( blk==0 ) clusblkspC = blkE/gOldConst_pass0[blkid];
	  cluspC += blkE/gOldConst_pass0[blkid]; 
	  hHCALblkpC->Fill( blkE/gOldConst_pass0[blkid] );
	  hHCALblkpCID->Fill( blkid, blkE/gOldConst_pass0[blkid] );
	  hHCALblkpCID_h->Fill( blkid, blkE/gOldConst_pass0[blkid] );
	  if( (int)nblk_h[f]==1 ){
	    hHCALblkpC_1b->Fill( blkE/gOldConst_pass0[blkid] );
	    hHCALblkpCID_1b->Fill( blkid, blkE/gOldConst_pass0[blkid] );
	  }

	}else{
	  if( blk==0 ) clusblkspC = blkE/blkagain;
	  cluspC += blkE/blkagain; 
	  hHCALblkpC->Fill( blkE/blkagain );
	  hHCALblkpCID->Fill( blkid, blkE/blkagain );
	  hHCALblkpCID_h->Fill( blkid, blkE/blkagain );
	  if( (int)nblk_h[f]==1 ){
	    hHCALblkpC_1b->Fill( blkE/blkagain );
	    hHCALblkpCID_1b->Fill( blkid, blkE/blkagain );
	  }

	}
	clusE += blkE;

	hblkid->Fill( blkid );
	hHCALblkeID->Fill( blkid, blkE );

	NEV[blkid]++;
	TNEV_h++;
      }

      hHCALe->Fill( clusE );
      hHCALeID->Fill( int(cblkid_h[f][0])-1, clusE ); //Fill with primary block location
      hHCALblkse->Fill( clusblksE );
      hHCALblkspC->Fill( clusblkspC );
      hSampFrac->Fill( clusE/KE_p );
      hSampFracID->Fill( int(cblkid_h[f][0])-1, clusE/KE_p ); //Fill with primary block location

      if(clusE<0.025){
	deadclus++;
	//cout << "Warning: Cluster block cut on atime killed cluster." << endl;
	continue;
      }

      //Calculate clusE variables
      Double_t E_dev = (E_exp-clusE)/E_exp; //Energy deviation as percent of expected
      hHCALeDev->Fill(E_dev);

    } //loop over events
  } //loop for lh2

  Int_t cell = 0;

  cout << endl << "Number of hydrogen events analyzed by cell: " << endl << endl;

  for( Int_t r = 0; r<kNrows; r++){
    for( Int_t c = 0; c<kNcols; c++){
      cout << NEV[cell] << "  ";
      cell++;
    }
    cout << endl;
  }
  
  cout << endl;

  usleep( 3*second ); //Give some time for review of step
  

  //Loop over all deuterium data
  for( Int_t f=0; f<nfset_ld2[kIdx]; f++ ){

    Int_t dfieldset = fset_ld2[kIdx][f];

    //Declare energy loss parameters for beam going through the target
    Double_t pBeam = E_e_d[f]/(1.+E_e_d[f]/M_p*(1.-cos(BB_th_d[f])));
    Double_t Eloss_outgoing = celldiameter/2.0/sin(BB_th_d[f]) * rho_tgt * dEdx_tgt; //Mean energy loss of the beam prior to the scattering, approximately 1 MeV, could correct further with raster position (likely negligable)

    for( Long64_t nevent = 1; nevent <Nevents_d[f]; nevent++){

      if ( nevent%10000==0 ) cout << "LD2 kinematic " << kine << " at field " << dfieldset << "%, entry: " << nevent << "/" << Nevents_d[f] << ". Total events gathered for analysis: " << TNEV_d << " \r";
      cout.flush();
      
      Cd[f]->GetEntry( elist_d[f]->GetEntry( nevent ) ); 

      //////////////////////////////////////////////////////////////////////
      //Make cuts that don't rely on projections first to improve efficiency
      // if( crow_d[f]==0 || 
      // 	  crow_d[f]==23 || 
      // 	  ccol_d[f]==0 || 
      // 	  ccol_d[f]==11 ) continue; //All events with primary cluster element on edge blocks cut
      if( abs(cblkatime_d[f][0]-atime0_d[f])>3*atime_sig_d[f] ) continue; //All events where adctime outside of reasonable window cut
      //////////////////////////////////////////////////////////////////////

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
      Double_t pp = sqrt(pow(nu,2)+2.*M_p*nu);
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
      //Double_t KE_p = Enucleon; //For sanity check on sampling fraction (2.28.23)
      Double_t dx = HCALx_d[f] - xexpect_dCAL;
      Double_t dy = HCALy_d[f] - yexpect_dCAL;
      
      /////////////////////////////////////////////////
      //Primary W2 cut on elastics and HCal active area
      if( fabs(W2-W2_mean_d[f])>W2_sig_d[f] ) continue; //Observed mean W2 cut on elastic peak
      elasYield++; //Events that pass the above cuts constitute elastics
      /////////////////////////////////////////////////

      //Calculate/declare new variables for analysis
      Double_t SFrac = HCALe_d[f]/KE_p;
      Double_t E_exp = KE_p*sampFrac[kIdx];
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
      if( !pass_p && !pass_n ) continue; //Cut on both n and p spots for each event, cannot know which apriori      
      //////////////////////////////////////////////////////////////////

      //Write out diagnostic histograms
      Double_t clusE = 0.0;
      Double_t clusblksE = cblke_d[f][0];
      Double_t clusblkspC = 0.0;
      Double_t cluspC = 0.0;
      Double_t blksatime = cblkatime_d[f][0];

      for( Int_t blk = 0; blk<(int)nblk_d[f]; blk++ ){
      	Int_t blkid = int(cblkid_d[f][blk])-1; //-1 necessary since sbs.hcal.clus_blk.id ranges from 1 to 288 (different than just about everything else)
	Double_t blkE = cblke_d[f][blk];
	Double_t blkatime = cblkatime_d[f][blk];
	Double_t blkagain = cblkagain_d[f][blk];
	Int_t blkr = blkid/kNcols;
	Int_t blkc = blkid%kNcols;

	//////////////////////////////
	//Time cut on cluster elements
	if( abs(blkatime-blksatime)>ifac*atime_sig_d[f] ){
	  hblkid_badatime->Fill(blkid);
	  hblkr_badatime->Fill(blkr);
	  hblkc_badatime->Fill(blkc);
	  continue;
	}
	//////////////////////////////
	
	if( pass0 ){
	  if( blk==0 ) clusblkspC = blkE/gOldConst_pass0[blkid];
	  cluspC += blkE/gOldConst_pass0[blkid]; 
	  hHCALblkpC->Fill( blkE/gOldConst_pass0[blkid] );
	  hHCALblkpCID->Fill( blkid, blkE/gOldConst_pass0[blkid] );
	  hHCALblkpCID_d->Fill( blkid, blkE/gOldConst_pass0[blkid] );
	  if( (int)nblk_d[f]==1 ){
	    hHCALblkpC_1b->Fill( blkE/gOldConst_pass0[blkid] );
	    hHCALblkpCID_1b->Fill( blkid, blkE/gOldConst_pass0[blkid] );
	  }

	}else{
	  if( blk==0 ) clusblkspC = blkE/blkagain;
	  cluspC += blkE/blkagain; 
	  hHCALblkpC->Fill( blkE/blkagain );
	  hHCALblkpCID->Fill( blkid, blkE/blkagain );
	  hHCALblkpCID_d->Fill( blkid, blkE/blkagain );
	  if( (int)nblk_d[f]==1 ){
	    hHCALblkpC_1b->Fill( blkE/blkagain );
	    hHCALblkpCID_1b->Fill( blkid, blkE/blkagain );
	  }

	}
	clusE += blkE;

	hblkid->Fill( blkid );
	hHCALblkeID->Fill( blkid, blkE );

	NEV[blkid]++;
	TNEV_d++;
      }

      hHCALe->Fill( clusE );
      hHCALeID->Fill( int(cblkid_d[f][0])-1, clusE ); //Fill with primary block location
      hHCALblkse->Fill( clusblksE );
      hHCALblkspC->Fill( clusblkspC );
      hSampFrac->Fill( clusE/KE_p );
      hSampFracID->Fill( int(cblkid_d[f][0])-1, clusE/KE_p ); //Fill with primary block location

      if(clusE<0.025){
	deadclus++;
	//cout << "Warning: Cluster block cut on atime killed cluster." << endl;
	continue;
      }

      //Calculate clusE variables
      Double_t E_dev = (E_exp-clusE)/E_exp; //Energy deviation as percent of expected
      hHCALeDev->Fill(E_dev);

    } //loop over events
  } //loop for ld2

  TCanvas *ADC_top = new TCanvas("ADC_top","ADC_top",1600,1200);
  TCanvas *ADC_bot = new TCanvas("ADC_bot","ADC_bot",1600,1200);
  ADC_top->Divide(12,12);
  ADC_bot->Divide(12,12);
  gStyle->SetOptStat(0);

  //Make arrays for adc tgraphs
  Double_t cellerr[kNcell] = {0.};
  Double_t pCcell[kNcell] = {0.};
  Double_t pCcval[kNcell] = {0.};
  Double_t pCcerr[kNcell] = {0.};
  Int_t fitmin = 600;
  TH1D *pCcellslice[kNcell];
  Double_t pCcval_avg = 0.;
  Double_t pCcerr_avg = 0.;
  Int_t pCcval_Ng = 0;
  Double_t pCsetpar[3];

  for(Int_t c=0; c<kNcell; c++){

    //Index through the canvas
    ADC_top->cd(c+1);
    if( c>=144 ){
      ADC_bot->cd(c-143);
      gStyle->SetOptStat(0);
    }

    //Get slices from haDiff_ID and fit for mean vals
    Double_t pCfitl = 0.;
    Double_t pCfith = 0.;
    Double_t ampl = pCbins/5;
    pCcell[c] = c;
    pCcellslice[c] = hHCALblkpCID_1b->ProjectionY(Form("pCcellslice_%d",c+1),c+1,c+1);
    pCcval[c] = 0.; //will overwrite if fit is good.

    Int_t sliceN = pCcellslice[c]->GetEntries();
    if( sliceN<fitmin ){
      pCcellslice[c]->Draw();
      continue;
    }
    Double_t arimean = pCcellslice[c]->GetMean();
    Double_t aristddev = pCcellslice[c]->GetStdDev();
    pCsetpar[0] = 10; //empirical limits
    pCsetpar[1] = arimean-aristddev;
    pCsetpar[2] = 2.*aristddev;
    pCfitl = 5.; //empirical minimum on pCvID plot
    pCfith = pCsetpar[1] + 3.5*aristddev; 
    TF1 *gausfit = new TF1("gausfit",Gfit,pCfitl,pCfith,3);
    gausfit->SetLineWidth(4);
    gausfit->SetParameter(0,pCsetpar[0]);
    gausfit->SetParLimits(0,1,ampl);
    gausfit->SetParameter(1,pCsetpar[1]);
    gausfit->SetParLimits(1,6,pCfith); //empirical limits
    gausfit->SetParameter(2,pCsetpar[2]);
    gausfit->SetParLimits(2,10,50);

    pCcellslice[c]->Fit("gausfit","RBM");
    pCcellslice[c]->Draw();

    pCcval[c] = gausfit->GetParameter(1);
    pCcerr[c] = gausfit->GetParameter(2);
    pCcellslice[c]->SetTitle(Form("E:%f sigma:%f",pCcval[c],pCcerr[c]));

    hsigEratio->Fill(pCcerr[c]/pCcval[c]);
    pCcval_avg += pCcval[c];
    pCcerr_avg += pCcerr[c];
    pCcval_Ng++;

  }
  pCcval_avg /= pCcval_Ng;
  cout << "pCcval_avg " << pCcval_avg << endl;
  pCcerr_avg /= pCcval_Ng;
  cout << "pCcerr_avg " << pCcerr_avg << endl;

  ADC_top->Write();
  ADC_bot->Write();

  // for(Int_t c=0; c<kNcell; c++){
  //   if( pCcval[c]==0. ){
  //     pCcval[c]=pCcval_avg;
  //     pCcerr[c]=pCcerr_avg;
  //   }
  // }

  TGraphErrors *gpC = new TGraphErrors( kNcell, pCcell, pCcval, cellerr, pCcerr );
  gpC->GetXaxis()->SetLimits(-10,290);  
  gpC->GetYaxis()->SetLimits(pCllim,pCulim);
  gpC->SetTitle("HCal ADC vs Cell");
  gpC->GetXaxis()->SetTitle("Cell");
  gpC->GetYaxis()->SetTitle("pC");
  gpC->SetMarkerStyle(20); // idx 20 Circles, idx 21 Boxes
  gpC->Draw();
  gpC->Write("gpC");

  fout->Write();

  cout << endl << endl;

  cout << "Calibration complete. Constants and histograms written to file." << endl;

  st->Stop();

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

}
