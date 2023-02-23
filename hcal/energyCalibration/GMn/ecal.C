//SSeeds 2.23.23 - Pass2 Update - Calibration code which employs best current cuts on elastic events to obtain ADC gain calibration parameters (pC/GeV). Adds beam energy loss to target and structure updates to easily run on batch farm. Restructured to remove memset and use sbs.hcal.again branch. Now loops over all data and cuts on both dx and dy with array of chains.

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
#include "TString.h"
#include "TF1.h"

//Detector constants
const Int_t kNcell = 288; // Total number of HCal modules
const Int_t kNrows = 24; // Total number of HCal rows
const Int_t kNcols = 12; // Total number of HCal columns
const Int_t maxTracks = 1000; // Reasonable limit on tracks to be stored per event
const Int_t uni_N = 400; // Total number of bins used to measure detection uniformity (hSampFrac histos)
const Int_t xN = 48; //2*kNrows, total number of dispersive bins detection uni
const Int_t yN = 24; //2*kNcols, total number of transverse bins detection uni
//Double_t hcalheight = 0.365; //m The height of the center of HCAL above beam
const Double_t hcalheight = -0.2897;
const Double_t sampFrac = 0.0795; //HCal sampling frac (0.06588 GeV/0.8286 GeV) = 0.0795 = 7.95% -> (MC E_dep per proton) / (fit to data KE_p)
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
const Int_t nkine = 6; // Total number of kinematic points
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

Double_t Gfit(Double_t *x, Double_t *par){
  Double_t amp = par[0];
  Double_t offset = par[1];
  Double_t sigma = par[2];
  return amp*exp(-0.5*pow((x[0]-offset)/sigma,2.));
}

//Passing kine=0 loads the test file
void ecal( Int_t kine=-1, Int_t iter=1 ){
  
  // if( kine==4 || kine==7 ){
  //   cout << "Warning: This script intended for use on kinematics post-pass0." << endl;
  // }

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
    if(fset_lh2[h]==-1) cout << "Error. configfilename array out of bounds." << endl;
    configfilename_h[h] = Form("/w/halla-scshelf2102/sbs/seeds/HCal_replay/hcal/energyCalibration/secal_lh2_sbs%d_f%d.cfg",kine,fset_lh2[h]);
  }
  string configfilename_d[nfset_ld2[kIdx]];
  for( Int_t d=0; d<nfset_ld2[kIdx]; h++){
    if(fset_ld2[d]==-1) cout << "Error. configfilename array out of bounds." << endl;
    configfilename_d[d] = Form("/w/halla-scshelf2102/sbs/seeds/HCal_replay/hcal/energyCalibration/secal_ld2_sbs%d_f%d.cfg",kine,fset_ld2[d]);
  }

  // Declare general physics/fit parameters
  Int_t minEventPerCell = 500; // Minimum number of scattered p in cell required to calibrate
  Int_t maxEventPerCell = 4000; // Maximum number of scattered p events to contribute
  Double_t HCal_divx = 0.15875; // Transverse width in x and y per cell
  Double_t HCal_divy = 0.15494;
  Double_t HCal_Xi = -2.355005; // Distance from beam center to top of HCal in m from MC database
  Double_t HCal_Xf = 1.454995; // Distance from beam center to bottom of HCal in m from MC database
  Double_t HCal_Yi = -0.92964; // Distance from beam center to opposite-beam side of HCal in m from MC database
  Double_t HCal_Yf = 0.92964;
  Double_t highDelta = 0.1; // Minimum M(i,j)/b(i) factor allowed 
  Double_t oldGain[kNcell] = {0.0};
  Double_t oldRatio[kNcell] = {0.0};
  Double_t maxDiff = 50.;

  // Declare general physics parameters to be modified by input config file - LH2
  Double_t E_e_h[nfset_lh2[kIdx]] = {0.}.; // Energy of beam (incoming electrons from accelerator)
  Double_t HCal_d_h[nfset_lh2[kIdx]] = {0.}.; // Distance to HCal from scattering chamber for comm1
  Double_t HCal_th_h[nfset_lh2[kIdx]] = {0.}.; // Angle that the center of HCal is at
  Double_t BB_th_h[nfset_lh2[kIdx]] = {0.}.; // Angle that the center of BBCal is at
  Double_t W2_mean_h[nfset_lh2[kIdx]] = {0.}.; // Mean of W at current kinematic
  Double_t W2_sig_h[nfset_lh2[kIdx]] = {0.}.; // Width of W at current kinematic
  Double_t dx0_n_h[nfset_lh2[kIdx]] = {0.}.; // Position of neutron spot, x-x_expected
  Double_t dx0_p_h[nfset_lh2[kIdx]] = {0.}.; // Position of proton spot, x-x_expected
  Double_t dy0_h[nfset_lh2[kIdx]] = {0.}.; // Position of hadron spot, y-y_expected
  Double_t dx_sig_n_h[nfset_lh2[kIdx]] = {0.}.; // Max spread of neutron spot, x-x_expected
  Double_t dx_sig_p_h[nfset_lh2[kIdx]] = {0.}.; // Max spread of proton spot, x-x_expected
  Double_t dy_sig_h[nfset_lh2[kIdx]] = {0.}.; // Max spread of hadron spot, y-y_expected
  Double_t atime0_h[nfset_lh2[kIdx]] = {0.}.; // Expected location in ADC time of signal
  Double_t atime_sig_h[nfset_lh2[kIdx]] = {0.}.; // 1 sig of atime distribution
  Int_t useAlshield_h[nfset_lh2[kIdx]] = {0.};

  // Declare general physics parameters to be modified by input config file - LD2
  Double_t E_e_d[nfset_ld2[kIdx]] = {0.}.; // Energy of beam (incoming electrons from accelerator)
  Double_t HCal_d_d[nfset_ld2[kIdx]] = {0.}.; // Distance to HCal from scattering chamber for comm1
  Double_t HCal_th_d[nfset_ld2[kIdx]] = {0.}.; // Angle that the center of HCal is at
  Double_t BB_th_d[nfset_ld2[kIdx]] = {0.}.; // Angle that the center of BBCal is at
  Double_t W2_mean_d[nfset_ld2[kIdx]] = {0.}.; // Mean of W at current kinematic
  Double_t W2_sig_d[nfset_ld2[kIdx]] = {0.}.; // Width of W at current kinematic
  Double_t dx0_n_d[nfset_ld2[kIdx]] = {0.}.; // Position of neutron spot, x-x_expected
  Double_t dx0_p_d[nfset_ld2[kIdx]] = {0.}.; // Position of proton spot, x-x_expected
  Double_t dy0_d[nfset_ld2[kIdx]] = {0.}.; // Position of hadron spot, y-y_expected
  Double_t dx_sig_n_d[nfset_ld2[kIdx]] = {0.}.; // Max spread of neutron spot, x-x_expected
  Double_t dx_sig_p_d[nfset_ld2[kIdx]] = {0.}.; // Max spread of proton spot, x-x_expected
  Double_t dy_sig_d[nfset_ld2[kIdx]] = {0.}.; // Max spread of hadron spot, y-y_expected
  Double_t atime0_d[nfset_ld2[kIdx]] = {0.}.; // Expected location in ADC time of signal
  Double_t atime_sig_d[nfset_ld2[kIdx]] = {0.}.; // 1 sig of atime distribution
  Int_t useAlshield_d[nfset_ld2[kIdx]] = {0.};

  Int_t elasYield = 0; // Keep track of total elastics analyzed
  Int_t badtimeblkclus = 0; // Keep track of multi-blk clusters with out of time blks
  Int_t multblkclus = 0; // Keep track of total multi-block clusters (nblk>1)

  //For position reconstruction
  Double_t HCal_Xmin = HCal_Xi-HCal_divx/2;
  Double_t HCal_Xmax = HCal_Xf+HCal_divx/2;
  Double_t HCal_Ymin = HCal_Yi-HCal_divy/2;
  Double_t HCal_Ymax = HCal_Yf+HCal_divy/2;

  // Reading all relevant config files for lh2
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
	  cout << "Loading beam energy: " << E_e << endl;
	}
	if( skey == "HCal_d" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  HCal_d_h[f] = sval.Atof();
	  cout << "Loading HCal distance: " << HCal_d << endl;
	}
	if( skey == "HCal_th" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  HCal_th_h[f] = sval.Atof() * TMath::DegToRad();	
	  cout << "Loading HCal angle: " << HCal_th << endl;
	}
	if( skey == "BB_th" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  BB_th_h[f] = sval.Atof() * TMath::DegToRad();	
	  cout << "Loading BBCal angle: " << BB_th << endl;
	}
	if( skey == "W2_mean" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  W2_mean_h[f] = sval.Atof();
	  cout << "Loading W2 mean cut: " << W2_mean << endl;
	}
	if( skey == "W2_sig" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  W2_sig_h[f] = sval.Atof();
	  cout << "Loading W2 sigma cut: " << W2_sig << endl;
	}
	if( skey == "dx0_n" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  dx0_n_h[f] = sval.Atof();
	  cout << "Loading x position of neutron spot: " << dx0_n << endl;
	}
	if( skey == "dx0_p" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  dx0_p_h[f] = sval.Atof();
	  cout << "Loading y position of proton spot: " << dx0_p << endl;
	}
	if( skey == "dy0" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  dy0_h[f] = sval.Atof();
	  cout << "Loading y position of both hadron spots: " << dy0 << endl;
	}
	if( skey == "dx_sig_n" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  dx_sig_n_h[f] = sval.Atof();
	  cout << "Loading x sigma of neutron spot: " << dx_sig_n << endl;
	}
	if( skey == "dx_sig_p" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  dx_sig_p_h[f] = sval.Atof();
	  cout << "Loading x sigma of proton spot: " << dx_sig_p << endl;
	}
	if( skey == "dy_sig" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  dy_sig_h[f] = sval.Atof();
	  cout << "Loading y sigma of both hadron spots: " << dy_sig << endl;
	}
	if( skey == "atime0" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  atime0_h[f] = sval.Atof();
	  cout << "Loading ADC time mean: " << atime0 << endl;
	}
	if( skey == "atime_sig" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  atime_sig_h[f] = sval.Atof();
	  cout << "Loading ADC time sigma: " << atime_sig << endl;
	}
	if( skey == "useAlshield" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  useAlshield_h[f] = sval.Atoi();
	  cout << "Loading Aluminum absorber option: " << useAlshield << endl;
	}
      }
      delete tokens;
    }

    NTevents_h[f] = Ch[f]->GetEntries();
    elist_h[f] = new TEventList(Form("elist_lh2_sbs%d_f%d",kine,fset_lh2[f]),Form("Elastic Event List, LH2, SBS%d, Mag%d",kine,fset_lh2[f]));
    Ch[f]->Draw(Form(">>elist_lh2_sbs%d_f%d",kine,fset_lh2[f]),globalcut_h[f]);
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
	  cout << "Loading beam energy: " << E_e << endl;
	}
	if( skey == "HCal_d" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  HCal_d_d[f] = sval.Atof();
	  cout << "Loading HCal distance: " << HCal_d << endl;
	}
	if( skey == "HCal_th" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  HCal_th_d[f] = sval.Atof() * TMath::DegToRad();	
	  cout << "Loading HCal angle: " << HCal_th << endl;
	}
	if( skey == "BB_th" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  BB_th_d[f] = sval.Atof() * TMath::DegToRad();	
	  cout << "Loading BBCal angle: " << BB_th << endl;
	}
	if( skey == "W2_mean" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  W2_mean_d[f] = sval.Atof();
	  cout << "Loading W2 mean cut: " << W2_mean << endl;
	}
	if( skey == "W2_sig" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  W2_sig_d[f] = sval.Atof();
	  cout << "Loading W2 sigma cut: " << W2_sig << endl;
	}
	if( skey == "dx0_n" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  dx0_n_d[f] = sval.Atof();
	  cout << "Loading x position of neutron spot: " << dx0_n << endl;
	}
	if( skey == "dx0_p" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  dx0_p_d[f] = sval.Atof();
	  cout << "Loading y position of proton spot: " << dx0_p << endl;
	}
	if( skey == "dy0" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  dy0_d[f] = sval.Atof();
	  cout << "Loading y position of both hadron spots: " << dy0 << endl;
	}
	if( skey == "dx_sig_n" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  dx_sig_n_d[f] = sval.Atof();
	  cout << "Loading x sigma of neutron spot: " << dx_sig_n << endl;
	}
	if( skey == "dx_sig_p" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  dx_sig_p_d[f] = sval.Atof();
	  cout << "Loading x sigma of proton spot: " << dx_sig_p << endl;
	}
	if( skey == "dy_sig" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  dy_sig_d[f] = sval.Atof();
	  cout << "Loading y sigma of both hadron spots: " << dy_sig << endl;
	}
	if( skey == "atime0" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  atime0_d[f] = sval.Atof();
	  cout << "Loading ADC time mean: " << atime0 << endl;
	}
	if( skey == "atime_sig" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  atime_sig_d[f] = sval.Atof();
	  cout << "Loading ADC time sigma: " << atime_sig << endl;
	}
	if( skey == "useAlshield" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  useAlshield_d[f] = sval.Atoi();
	  cout << "Loading Aluminum absorber option: " << useAlshield << endl;
	}
      }
      delete tokens;
    }

    NTevents_d[f] = Cd[f]->GetEntries();
    elist_d[f] = new TEventList(Form("elist_ld2_sbs%d_f%d",kine,fset_ld2[f]),Form("Elastic Event List, LD2, SBS%d, Mag%d",kine,fset_ld2[f]));
    Cd[f]->Draw(Form(">>elist_ld2_sbs%d_f%d",kine,fset_ld2[f]),globalcut_d[f]);
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

  cout << endl << endl << "Setup parameters loaded." << endl;

  // Create stopwatch to track processing time
  TStopwatch *sw = new TStopwatch();
  sw->Start();
  
  // Declare outfile
  TFile *fout = new TFile( Form("eCalEOut_sbs%d_test.root", kine ), "RECREATE" );
  
  // Initialize vectors and arrays
  Double_t GCoeff[kNcell] = {0.0};
  Double_t GCoeff_oneblock[kNcell] = {0.0};
  Double_t GCoeff_divide[kNcell] = {0.0};
  Double_t gOldConst[kNcell] = {0.0};

  // Initialize histograms
  // TH1D *hatime = new TH1D( "hatime", "Cluster adc time, primary block; ns", 220, -50, 170);
  // TH1D *hatime_cut = new TH1D( "hatime_cut", "Cluster adc time, primary block, W2 cut; ns", 220, -50, 170);
  // TH1D *hpC = new TH1D( "hpC", "Cluster element pC", 2000,0,2000);
  // TH2D *hpC_vs_ID = new TH2D( "hpC_vs_ID", "Cluster element pC vs Channel", 300,-10,290,2000,0,2000);
  // TH2D *hpC_vs_ID_oneblock = new TH2D( "hpC_vs_ID_oneblock", "Cluster element pC vs Channel Single Block Clusters", 288,0,288,2000,0,2000);
  // TH2D *hCoeffvID_oneblock = new TH2D( "hCoeffvID_oneblock", "Energy Expected / Integrated ADC; channel; GeV/pC", 288,0,288,100,0,0.03);
  // TH1D *hDeltaE = new TH1D( "hDeltaE","1.0-Eclus/p_rec", 100, -1.5, 1.5 );
  // TH1D *hHCALe = new TH1D( "hHCALe","HCal Cluster E", 400, 0., 1. );
  // TH1D *hHCALe_cut = new TH1D( "hHCALe_cut","HCal Cluster E, W2 cut", 400, 0., 1. );
  // TH1D *hSampFrac = new TH1D( "hSampFrac","W2 Cut HCal Cluster E / Expected KE", 400, 0., 1. );
  // TH1D *hDiff = new TH1D( "hDiff","HCal time - BBCal time (ns)", 1300, -500, 800 );
  // TH1D *hEDiffChan = new TH1D( "hEDiffChan","EDiff Events Per Channel", 288, 0, 288 );
  // TH1D *hClusE_calc = new TH1D( "hClusE_calc","Best Cluster Energy Calculated from HCal Vars", 100, 0.0, 2.0);
  // TH2D *hClusBlk_ClusE = new TH2D( "hClusBlk_ClusE","Cluster Size vs Recon Cluster E", 8, 0, 8, 100, 2.0, 3.0 );

  // TH1D *hmclusblkE = new TH1D( "hmclusblkE","HCal multi-blk cluster blk E", 400, 0., 1. );
  // TH1D *hsclusblkE = new TH1D( "hsclusblkE","HCal single-blk cluster blk E", 400, 0., 1. );

  // TH2D *hPAngleCorr = new TH2D( "hPAngCorr","Track p vs Track ang", 100, 30, 60, 100, 0.4, 1.2 );
  // TH2D *hPAngleCorr_2 = new TH2D( "hPAngCorr_2","Track p vs Track ang v2", 100, 30, 60, 100, 0.4, 1.2 );
  // TH2D *hClusE_vs_X = new TH2D("hClusE_vs_X",";X-pos (m);E_dep (GeV)",500,-3,2,100,0.0,1.0);  
  // TH2D *hClusE_vs_Y = new TH2D("hClusE_vs_Y",";Y-pos (m);E_dep (GeV)",200,-1,1,100,0.0,1.0);    
  // TH1D *hpp = new TH1D( "hpp", "Elastic Proton Momentum", 600, 0, 6 );
  // TH1D *hdx = new TH1D( "hdx","; x_{HCAL} - x_{exp} (m)", 250, -dxlim_l,dxlim_u);
  // TH1D *hdy = new TH1D("hdy","; y_{HCAL} - y_{exp} (m);", 250, -1.25,1.25);
  // TH2D *hdxdy = new TH2D("hdxdy",";y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)", 250, -1.25, 1.25, 250, -dxlim_l, dxlim_u );
  // TH1D *hW2_calc = new TH1D( "W2_calc", "W2 Calculated; GeV", 250, 0.3, 1.5 );
  // TH1D *hW2_tree = new TH1D( "W2_tree", "W2 From Tree; GeV", 250, 0.3, 1.5 );
  // TH1D *hNBlk = new TH1D( "hNBlk", "Number of Blocks in Primary Cluster", 25, 0, 25 );
  // TH1D *hNBlk_cut = new TH1D( "hNBlk_cut", "Number of Blocks in Primary Cluster W2 Cut", 25, 0, 25 );
  // TH1D *hQ2 = new TH1D( "Q2", "Q2; GeV", 250, 0.5, 3.0 );
  // TH1D *hE_ep = new TH1D( "hE_ep","Scattered Electron Energy; GeV", 500, 0.0, E_e*1.5 ); 
  // TH1D *hE_pp = new TH1D( "hE_pp", "Scattered Proton Energy; GeV", 500, 0.0, E_e*1.5 );
  // TH1D *hKE_p = new TH1D( "hKE_p", "Scattered Proton Kinetic Energy; GeV", 500, 0.0, E_e*1.5 );
  // TH2D *hADC = new TH2D( "hADC", "HCal Int_ADC Spectra: W2 Cut", 288, 0, 288, 100., 0., 1. );
  // TH2D *coefficients = new TH2D( "coefficients",";Channel ;GeV/mV", 288,0,288,250,0,0.025 );
  // TH2D *coefficients_1b = new TH2D( "coefficients_1b",";Channel ;GeV/mV", 288,0,288,250,0,0.025 );
  // TH2D *hEDiff_vs_X = new TH2D( "hEDiff_vs_X",";X-pos (m);(E_exp-E_dep)/E_exp (GeV)",xN,HCal_Xi,HCal_Xf,100,-3,3 ); 
  // TH2D *hEDiff_vs_Y = new TH2D( "hEDiff_vs_Y",";Y-pos (m);(E_exp-E_dep)/E_exp (GeV)",yN,HCal_Yi,HCal_Yf,100,-3,3 );
  // TH2D *hEDiff_vs_block = new TH2D( "hEDiff_vs_block",";Block (cell);(E_exp-E_dep)/E_exp (GeV)",288,0,288,100,-3,3 );
  // TH2D *hSampFrac_vs_X = new TH2D( "hSampFrac_vs_X","Sampling Fraction ;X (m) ;HCal_E / Exp_KE", xN, HCal_Xi, HCal_Xf, uni_N, 0., 1. );
  // TH2D *hSampFrac_vs_Y = new TH2D( "hSampFrac_vs_Y","Sampling Fraction ;Y (m) ;HCal_E / Exp_KE", yN, HCal_Yi, HCal_Yf, uni_N, 0., 1. );
  // TH1D *hE_pp_cut = new TH1D( "hE_pp_cut", "Deposited Elastic Proton Energy; GeV", 500, 0.0, E_e*1.5 );

  // Set long Int_t to keep track of total entries
  //Long64_t Nevents = elist->GetN();

  cout << endl << "All parameters loaded." << endl << endl;
  cout << "Opened up tree with nentries: " << Nevents << ".." << endl << endl;

  //Declare matrices for chi-square min calibration scheme and keep track of calibrated events
  TMatrixD Ma(kNcell,kNcell);
  TMatrixD Ma_oneblock(kNcell,kNcell);
  TMatrixD Ma_err(kNcell,kNcell);
  TVectorD ba(kNcell);
  TVectorD ba_oneblock(kNcell);
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
  //if( useAlshield != 0 ) Eloss_outgoing += Alshieldthick * rho_Al * dEdx_Al;

  //Loop over events
  cout << "Main loop over all data commencing.." << endl;

  //Loop over all hydrogen data
  for( Int_t f=0, f<nfset_lh2[kIdx]; f++ ){

    for( Long64_t nevent = 1; nevent <Nevents_h[f]; nevent++){

      if ( nevent%10000==0 ) cout << "LH2 kinematic: " << kine << ", entry: " << nevent << "/" << Nevents[f] << " \r";
      cout.flush();
      
      Ch[f]->GetEntry( elist_h[f]->GetEntry( nevent ) ); 


      //////////////////////////////////////////////////////////////////////
      //Make cuts that don't rely on projections first to improve efficiency
      if( crow_h[f]==0 || 
	  crow_h[f]==23 || 
	  ccol_h[f]==0 || 
	  ccol_h[f]==11 ) continue; //All events with primary cluster element on edge blocks cut
      if( abs(cblkatime_h[f][0]-atime0_h[f])>3*atime_sig_h[f] ) continue; //All events where adctime outside of reasonable window cut
      //////////////////////////////////////////////////////////////////////


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
      Double_t dx = HCALx_h[f] - xexpect_HCAL;
      Double_t dy = HCALy_h[f] - yexpect_HCAL;
      
      /////////////////////////////////////////////////
      //Primary W2 cut on elastics and HCal active area
      if( fabs(W2-W2_mean_h[f])>W2_sig_h[f] ) continue; //Observed mean W2 cut on elastic peak
      elasYield++; //Events that pass the above cuts constitute elastics
      /////////////////////////////////////////////////

      //Calculate/declare new variables for analysis
      Double_t SFrac = HCALe_h[f]/KE_p;
      Double_t E_exp = KE_p*sampFrac;
      Int_t rec_row = ( HCALx_h[f] - HCal_Xmin )/HCal_divx;
      Int_t rec_col = ( HCALy_h[f] - HCal_Ymin )/HCal_divy;
      Int_t rec_cell = rec_row*kNcols + rec_col;
      
      //////////////////////////////////////////////////////////////////
      //Cut on dx and dy.
      if( abs(dy-dy0_h[f])<5*dy_sig_h[f] ) continue;
      if( abs(dx-dx0_h[f])<5*dx_sig_h[f] ) continue;
      //////////////////////////////////////////////////////////////////

      // Get energies with simplest scheme from clusters only
      bool badclus = false;
      Double_t clusE = 0.0;
      Double_t cluspC = 0.0;
      for( Int_t blk = 0; blk<(int)nblk_h[f]; blk++ ){
	Int_t blkid = int(cblkid_h[f][blk])-1; //-1 necessary since sbs.hcal.clus_blk.id ranges from 1 to 288 (different than just about everything else)
	
	//Populate list of old gain coefficients in case dof is bad
	if( gOldConst[blkid]==0.0 ) gOldConst[blkid]=cblkagain_h[f][blk];
	
	//Diagnostic cut for now
	if( abs(cblkatime_h[f][blk]-cblkatime_h[f][0])>3*atime_sig ) {
	  badclus = true;
	  continue; //Cut blocks that are out of time
	}
	
	clusE += cblke_h[f][blk];
	cluspC += cblke_h[f][blk]/cblkagain_h[f][blk]; 
	A[blkid] += cblke_h[f][blk]/cblkagain_h[f][blk];
	
	if(nblk_h[f]==1) {
	  A_oneblock[blkid] += cblke_h[f][blk]/cblkagain_h[f][blk];
	}
	
	err_ev[blkid] += pow( ((1.+cluspC)/cluspC+(0.05*KE_p)/KE_p) , 2 ); //Getting ready for std. dev of mean

	// Simple estimation of the coefficients assuming 100% energy deposition in one block. Will check against the more sophisticated version with chi^2 reduction.
	if(nblk_h[f]==1) {
	 
	  err_ev_oneblock[blkid] += pow( ((1.+0.5*cluspC)/cluspC+(0.05*KE_p)/KE_p) , 2 ); //Probably better approximation of the error since all the KE_p should be deposited in this block.
	  NEV_oneblock[blkid]++;
	}
	NEV[blkid]++;

      }

      if( badclus ) badtimeblkclus++;
      if( nblk_h[f]>1 ) multblkclus++;

      //Calculate clusE variables
      Double_t E_dev = (E_exp-clusE)/E_exp; //Energy deviation as percent of expected
      
      //Build the matrix as simply as possible
      for(Int_t icol = 0; icol<kNcell; icol++){
	ba(icol)+= A[icol];

	for(Int_t irow = 0; irow<kNcell; irow++){
	  Ma(icol,irow) += A[icol]*A[irow]/E_exp;
	  if(nblk_h[f]==1){	    
	    // cout << "do something to correct oneblock" << endl;
	    Ma_oneblock(icol,irow) += A_oneblock[icol]*A_oneblock[irow]/E_exp;
	  } //oneblock
	} //loop over matrix element rows
      } //loop over matrix element cols
    } //loop over events
  } //loop for lh2

  //Loop over all deuterium data
  for( Int_t f=0, f<nfset_ld2[kIdx]; f++ ){

    for( Long64_t nevent = 1; nevent <Nevents_d[f]; nevent++){

      if ( nevent%10000==0 ) cout << "LH2 kinematic: " << kine << ", entry: " << nevent << "/" << Nevents[f] << " \r";
      cout.flush();
      
      Ch[f]->GetEntry( elist_d[f]->GetEntry( nevent ) ); 


      //////////////////////////////////////////////////////////////////////
      //Make cuts that don't rely on projections first to improve efficiency
      if( crow_d[f]==0 || 
	  crow_d[f]==23 || 
	  ccol_d[f]==0 || 
	  ccol_d[f]==11 ) continue; //All events with primary cluster element on edge blocks cut
      if( abs(cblkatime_d[f][0]-atime0_d[f])>3*atime_sig_d[f] ) continue; //All events where adctime outside of reasonable window cut
      //////////////////////////////////////////////////////////////////////


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
      Double_t dx = HCALx_d[f] - xexpect_dCAL;
      Double_t dy = HCALy_d[f] - yexpect_dCAL;
      
      /////////////////////////////////////////////////
      //Primary W2 cut on elastics and HCal active area
      if( fabs(W2-W2_mean_d[f])>W2_sig_d[f] ) continue; //Observed mean W2 cut on elastic peak
      elasYield++; //Events that pass the above cuts constitute elastics
      /////////////////////////////////////////////////

      //Calculate/declare new variables for analysis
      Double_t SFrac = HCALe_d[f]/KE_p;
      Double_t E_exp = KE_p*sampFrac;
      Int_t rec_row = ( HCALx_d[f] - HCal_Xmin )/HCal_divx;
      Int_t rec_col = ( HCALy_d[f] - HCal_Ymin )/HCal_divy;
      Int_t rec_cell = rec_row*kNcols + rec_col;
      
      //////////////////////////////////////////////////////////////////
      //Cut on dx and dy.
      if( abs(dy-dy0_d[f])<5*dy_sig_d[f] ) continue;
      if( abs(dx-dx0_d[f])<5*dx_sig_d[f] ) continue;
      //////////////////////////////////////////////////////////////////

      // Get energies with simplest scheme from clusters only
      bool badclus = false;
      Double_t clusE = 0.0;
      Double_t cluspC = 0.0;
      for( Int_t blk = 0; blk<(int)nblk_d[f]; blk++ ){
	Int_t blkid = int(cblkid_d[f][blk])-1; //-1 necessary since sbs.hcal.clus_blk.id ranges from 1 to 288 (different than just about everything else)
	
	//Populate list of old gain coefficients in case dof is bad
	if( gOldConst[blkid]==0.0 ) gOldConst[blkid]=cblkagain_d[f][blk];
	
	//Diagnostic cut for now
	if( abs(cblkatime_d[f][blk]-cblkatime_d[f][0])>3*atime_sig ) {
	  badclus = true;
	  continue; //Cut blocks that are out of time
	}
	
	clusE += cblke_d[f][blk];
	cluspC += cblke_d[f][blk]/cblkagain_d[f][blk]; 
	A[blkid] += cblke_d[f][blk]/cblkagain_d[f][blk];
	
	if(nblk_d[f]==1) {
	  A_oneblock[blkid] += cblke_d[f][blk]/cblkagain_d[f][blk];
	}
	
	err_ev[blkid] += pow( ((1.+cluspC)/cluspC+(0.05*KE_p)/KE_p) , 2 ); //Getting ready for std. dev of mean

	// Simple estimation of the coefficients assuming 100% energy deposition in one block. Will check against the more sophisticated version with chi^2 reduction.
	if(nblk_d[f]==1) {
	 
	  err_ev_oneblock[blkid] += pow( ((1.+0.5*cluspC)/cluspC+(0.05*KE_p)/KE_p) , 2 ); //Probably better approximation of the error since all the KE_p should be deposited in this block.
	  NEV_oneblock[blkid]++;
	}
	NEV[blkid]++;

      }

      if( badclus ) badtimeblkclus++;
      if( nblk_d[f]>1 ) multblkclus++;

      //Calculate clusE variables
      Double_t E_dev = (E_exp-clusE)/E_exp; //Energy deviation as percent of expected
      
      //Build the matrix as simply as possible
      for(Int_t icol = 0; icol<kNcell; icol++){
	ba(icol)+= A[icol];

	for(Int_t irow = 0; irow<kNcell; irow++){
	  Ma(icol,irow) += A[icol]*A[irow]/E_exp;
	  if(nblk_d[f]==1){	    
	    // cout << "do something to correct oneblock" << endl;
	    Ma_oneblock(icol,irow) += A_oneblock[icol]*A_oneblock[irow]/E_exp;
	  } //oneblock
	} //loop over matrix element rows
      } //loop over matrix element cols
    } //loop over events
  } //loop for ld2

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

    badcell_oneblock[i] = 0;

    //Do not change ADC gain coeff if insufficient events or energy dep in cell, oneblock
    if( NEV[i] < minEventPerCell || Ma_oneblock(i,i) < 0.1*ba_oneblock(i) ){ 

      cellBad = 1;

      Double_t elemRatio = Ma_oneblock(i,i)/ba_oneblock(i);

      ba_oneblock(i) = 1.0;  // Set RHS vector for cell i to 1.0 
      Ma_oneblock(i,i) = 1.0; // Set diagonal element of Matrix M for cell i to 1.0 
      
      for(Int_t j=0; j<kNcell; j++){
	if( j != i ){
	  
	  Ma_oneblock(i,j) = 0.0; 
	  Ma_oneblock(j,i) = 0.0;
	}
      }
      badcell_oneblock[i] = 1;
      cout << "oneblock analysis" << endl;
      cout << "cell" << i << " bad" << endl;
      cout << "Number of events in bad cell =" << NEV[i] << endl;
      cout << "Matrix element/vector element ratio =" << elemRatio << endl;
    }   

    //Calculate error per block (element of A)
    err[i] = sqrt(err_ev[i]/NEV[i]);
    err_oneblock[i] = sqrt(err_ev_oneblock[i]/NEV_oneblock[i]);
  }

  if( cellBad==0 ) cout << "No bad cells detected!" << endl << endl;

  //Invert the matrix, solve for ratios
  TMatrixD M_inv = Ma.Invert();
  TMatrixD M_inv_oneblock = Ma_oneblock.Invert();
  TVectorD Coeff = M_inv*ba; // Stays unmodified for reference
  TVectorD Coeff_oneblock = M_inv_oneblock*ba_oneblock; // Stays unmodified for reference

  for(Int_t i=0; i<kNcell; i++){
    if(badcell_oneblock[i]==0){
      GCoeff_oneblock[i]=Coeff_oneblock[i];
    }else{
      GCoeff_oneblock[i]=gOldConst[i];
    }

    if(badcell[i]==0){
      GCoeff[i]=Coeff[i];
      GCoeff_divide[i]=GCoeff[i]/GCoeff_oneblock[i];

    }else{
      GCoeff[i]=gOldConst[i]; // If the cell is bad, use the old coefficient
      GCoeff_divide[i]=-1.0;
    }
  }

  Double_t yErr[kNcell] = {0.0};

  TGraphErrors *ccgraph = new TGraphErrors( kNcell, y, GCoeff, yErr, err ); 
  ccgraph->GetXaxis()->SetLimits(0.0,288.0);  
  ccgraph->GetYaxis()->SetLimits(0.0,0.25);
  ccgraph->SetTitle("Calibration Coefficients");
  ccgraph->GetXaxis()->SetTitle("Channel");
  ccgraph->GetYaxis()->SetTitle("Unitless");
  ccgraph->SetMarkerStyle(20); // idx 20 Circles, idx 21 Boxes
  ccgraph->Write("constants");

  TGraphErrors *ccgraph_oneBlock = new TGraphErrors( kNcell, y, GCoeff_oneblock, yErr, err_oneblock ); 
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

  fout->Write();

  //Console/txt outs
  Int_t cell = 0;

  cout << endl << "Number of events available for calibration: " << endl << endl;

  for( Int_t r = 0; r<kNrows; r++){
    for( Int_t c = 0; c<kNcols; c++){
      cout << NEV[cell] << "  ";
      cell++;
    }
    cout << endl;
  }
  cell = 0;

  cout << endl << "One Block Std:" << endl;

  for( Int_t r = 0; r<kNrows; r++){
    for( Int_t c = 0; c<kNcols; c++){
      cout << GCoeff_oneblock[cell] << "  ";
      cell++;
    }
    cout << endl;
  }

  cout << endl << endl;

  //Declare outfiles
  ofstream GainCoeff;

  GainCoeff.open( constPath );
  GainCoeff << "#HCal gain coefficients from SBS-" << kine << " obtained " << date.c_str() << endl;
  GainCoeff << "#HCal_gainCoeff = " << endl;

  cell = 0;
  cout << "Gain Coefficients" << endl;
  for( Int_t r=0; r<kNrows; r++ ){
    for( Int_t c=0; c<kNcols; c++ ){
      GainCoeff << GCoeff[cell] << "  ";
      cout << GCoeff[cell] << "  ";
      cell++;
    }
    GainCoeff << endl;
    cout << endl;
  }

  GainCoeff.close();

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
