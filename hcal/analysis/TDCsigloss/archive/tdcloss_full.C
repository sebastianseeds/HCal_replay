//SSeeds 1.19.23 - Short script to run over all data from a given kinematic and check tdc signals with elastic cuts. Configured to include full hcal replays and cosmic data

#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TString.h"
#include <iostream>
#include <fstream>
//#include "hcal.h"
#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <TChain.h>
#include <TSystem.h>
#include <TStopwatch.h>
#include "TMath.h"
//#include "../../include/gmna.h" 

const Int_t maxTracks = 1000; // Reasonable limit on tracks to be stored per event
//Trigger TDC
const Int_t maxTDCTrigChan = 10; // Set to accomodate original 5 TDCTrig channels with buffer
const Double_t tdiffwidecut = 50; // Set at 50ns from nominal 510ns (passed by user) from GMn
//HCal - Note that actual measurement of vertical length is 381.6cm, indicating that the MC figures are correct
const Int_t maxHCalChan = 288; // Total HCal channels
const Int_t maxHCalRows = 24; // Total HCal rows
const Int_t maxHCalCols = 12; // Total HCal cols
const Int_t maxClusters = 10; // Total HCal clusters with information saved
//const Double_t HCALHeight = 0.365; // Height of HCal above beamline in m
const Double_t HCALHeight = -0.2897; // Height of HCal above beamline in m
const Double_t HCalblk_l = 0.15; // Width and height of HCAL blocks in m
const Double_t HCalblk_l_h_MC = 0.15494; // Horizontal length of all HCAL blocks in m from MC database
const Double_t HCalblk_l_v_MC = 0.15875; // Vertical length of all HCAL blocks in m from MC database
const Double_t posHCalXi = -2.165; // Distance from beam center to top of HCal in m
const Double_t posHCalXf = 1.435; // Distance from beam center to bottom of HCal in m
const Double_t posHCalYi = -0.9; // Distance from beam center to opposite-beam side of HCal in m
const Double_t posHCalYf = 0.9; // Distance from beam center to beam side of HCal in m
const Double_t posHCalXi_MC = -2.355005; // Distance from beam center to top of HCal in m from MC database
const Double_t posHCalXf_MC = 1.454995; // Distance from beam center to bottom of HCal in m from MC database
const Double_t posHCalYi_MC = -0.92964; // Distance from beam center to opposite-beam side of HCal in m from MC database
const Double_t posHCalYf_MC = 0.92964; // Distance from beam center to beam side of HCal in m from MC database
const Double_t HCalSampFrac = 0.077;  //Re-evaluated with MC GEn settings using second to outermost shower column for kin2
//Physics/Math
const Double_t PI = TMath::Pi();
const Double_t M_e = 0.00051;
const Double_t M_p = 0.938272;
const Double_t M_n = 0.939565;
const UInt_t us = 1000000; //For conversion to seconds used by reporting time delays
//Static Target/Scattering Chamber Parameters
const Double_t l_tgt = 0.15; // Length of the target (m)
const Double_t rho_tgt = 0.0723; // Density of target (g/cc)
const Double_t rho_Al = 2.7; // Density of aluminum windows (g/cc)
const Double_t celldiameter = 1.6*2.54; //cm, right now this is a guess
const Double_t Ztgt = 1.0;
const Double_t Atgt = 1.0;
const Double_t Mmol_tgt = 1.008; //g/mol
const Double_t dEdx_tgt=0.00574; //According to NIST ESTAR, the collisional stopping power of hydrogen is about 5.74 MeV*cm2/g at 2 GeV energy
const Double_t dEdx_Al = 0.0021; //According to NIST ESTAR, the collisional stopping power of Aluminum is about 2.1 MeV*cm2/g between 1-4 GeV
const Double_t uwallthick_LH2 = 0.0145; //cm
const Double_t dwallthick_LH2 = 0.015; //cm
const Double_t cellthick_LH2 = 0.02; //cm, this is a guess;
const Double_t Alshieldthick = 2.54/8.0; //= 1/8 inch * 2.54 cm/inch 
const Double_t tdiffmax = 20; //ns. Maximum allowed deviation from HCal trig = BBCal trig = t_trig (from config file)

//Use kine=0 for cosmics
void tdcloss_full( Int_t kine = 4 ){

  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  TChain *C = new TChain("T");

  // Declare simple variables and strings to keep track of file extensions
  //TString configfilename = Form( "%s/parse_conf/SBS%d/%s/mag%d/sbs%d%s%dpar.cfg", config_prefix, kine, tar, mag, kine, tar, mag );
  TString outputfilename = Form("tdcloss_full_sbs%d_out.root",kine);
  vector<TString> log;

  // Declare general physics parameters to be modified by input config file
  Double_t E_e = -1.; // Energy of beam (incoming electrons from accelerator)
  Double_t t_trig = -1.; // Time difference between BBCal Hi (L1A) and HCal Trig (ns)
  Double_t BB_d = -1.; // Distance to bigbite spectrometer from target chamber (m)
  Double_t BB_th = -1.; // Angle BB spectrometer makes with exit beamline
  Double_t HCal_d = -1.; // Distance to HCal from scattering chamber for comm1
  Double_t HCal_th = -1.; // Angle HCal center makes with exit beamline  
  Double_t W2_mean = -1.; // Mean of W at current kinematic
  Double_t W2_sig = -1.; // Width of W at current kinematic
  Int_t useAlshield = -1.; //Use 1/8" al shield on scattering chamber exit? 1:yes 0:no
  Int_t debug = 3;

  // Reading config file
  TString configfilename = "stdcloss.cfg";
  if( kine==4 ){
    configfilename = "stdcloss_sbs4.cfg";
  }else if( kine==11 ){
    configfilename = "stdcloss_sbs11.cfg";
  }

  ifstream configfile(configfilename);
  TString currentline;
  cout << endl << "Chaining the following runs: " << endl;
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      if(!currentline) cout << "WARNING: No file exists at " << currentline << "." << endl;
      //C->Add(currentline);
      if(debug==0){
	C->Add(currentline);
	debug=1;
      }else if(debug==3){
	C->Add(currentline);
      }
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
      if( skey == "W2_mean" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        W2_mean = sval.Atof();
	cout << "W2 mean cut: " << W2_mean << endl;
      }
      if( skey == "W2_sig" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        W2_sig = sval.Atof();
	cout << "W2 sigma cut: " << W2_sig << endl;
      } 
      if( skey == "useAlshield" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	useAlshield = sval.Atoi();
	cout << "Aluminum absorber option: " << useAlshield << endl;
      } 
    }
    delete tokens;
  }
  
  if( E_e==-1 || 
      t_trig==-1 || 
      BB_d==-1 || 
      BB_th==-1 || 
      HCal_d==-1 || 
      HCal_th==-1 || 
      W2_mean==-1 || 
      W2_sig==-1 || 
      useAlshield==-1 ){
    cout << "Error: Unable to load setup parameters correctly. Please source setup_gmna.sh and try again. If problem persists, check configuration file." << endl;
    return;
  }

  cout << "Setup parameters loaded." << endl;

  cout << endl;

  // Declare general detector and physics parameters
  Double_t TDCT_id[maxTDCTrigChan], TDCT_tdc[maxTDCTrigChan]; 
  Int_t TDCTndata;

  // BB params
  Double_t BBtr_p[maxTracks], BBtr_px[maxTracks], BBtr_py[maxTracks], BBtr_pz[maxTracks];
  Double_t BBtr_vz[maxTracks], BBtr_chi2[maxTracks];
  Double_t BBtr_n, BBps_x, BBps_y, BBps_e, BBsh_x, BBsh_y, BBsh_e, GEMtr_hits, BBeoverp;

  // HCal params
  Double_t hcal_id[maxHCalChan], hcal_tdc[maxHCalChan], hcal_amp_p[maxHCalChan], hcal_amp[maxHCalChan], hcal_atime[maxHCalChan]; 
  Int_t hcalndata;
  Double_t cid[maxHCalChan], crow[maxHCalRows], ccol[maxHCalCols];
  Double_t ce[maxHCalChan], catime[maxHCalChan], ctdc[maxHCalChan];
  Double_t nclus, nblk;
  Double_t hcalx, hcaly, hcale, kineW2;

  //Switch off all unnecessary branches
  C->SetBranchStatus( "*", 0 );
  C->SetBranchStatus( "sbs.hcal.x", 1 );
  C->SetBranchStatus( "sbs.hcal.y", 1 );
  C->SetBranchStatus( "sbs.hcal.e", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.row", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.col", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.tdctime", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.atime", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.id", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.e", 1 );
  C->SetBranchStatus( "sbs.hcal.nblk", 1 );
  C->SetBranchStatus( "sbs.hcal.nclus", 1 );
  C->SetBranchStatus( "bb.tr.n", 1 );
  C->SetBranchStatus( "bb.tr.p", 1 );
  C->SetBranchStatus( "bb.tr.px", 1 );
  C->SetBranchStatus( "bb.tr.py", 1 );
  C->SetBranchStatus( "bb.tr.pz", 1 );
  C->SetBranchStatus( "bb.tr.vx", 1 );
  C->SetBranchStatus( "bb.tr.vy", 1 );
  C->SetBranchStatus( "bb.tr.vz", 1 );
  C->SetBranchStatus( "bb.ps.e", 1 );
  C->SetBranchStatus( "bb.ps.x", 1 );
  C->SetBranchStatus( "bb.ps.y", 1 );
  C->SetBranchStatus( "bb.sh.e", 1 );
  C->SetBranchStatus( "bb.sh.x", 1 );
  C->SetBranchStatus( "bb.sh.y", 1 );
  C->SetBranchStatus( "bb.tdctrig.tdc", 1 );
  C->SetBranchStatus( "bb.tdctrig.tdcelemID", 1 );
  C->SetBranchStatus( "Ndata.bb.tdctrig.tdcelemID", 1 );
  C->SetBranchStatus( "sbs.hcal.a_amp", 1 );
  C->SetBranchStatus( "sbs.hcal.a_amp_p", 1 );
  C->SetBranchStatus( "sbs.hcal.a_time", 1 );
  C->SetBranchStatus( "sbs.hcal.tdc", 1 );
  C->SetBranchStatus( "sbs.hcal.tdcelemID", 1 );
  C->SetBranchStatus( "Ndata.sbs.hcal.tdcelemID", 1 );
  C->SetBranchStatus( "e.kine.W2", 1 );
  //C->SetBranchStatus( "bb.gem.track.nhits", 1 );
  //C->SetBranchStatus( "bb.etot_over_p", 1 );

  //Link the branches to vars
  C->SetBranchAddress( "sbs.hcal.x", &hcalx );
  C->SetBranchAddress( "sbs.hcal.y", &hcaly );
  C->SetBranchAddress( "sbs.hcal.e", &hcale );
  C->SetBranchAddress( "sbs.hcal.clus_blk.row", crow );
  C->SetBranchAddress( "sbs.hcal.clus_blk.col", ccol );
  C->SetBranchAddress( "sbs.hcal.clus_blk.tdctime", ctdc );
  C->SetBranchAddress( "sbs.hcal.clus_blk.atime", catime );
  C->SetBranchAddress( "sbs.hcal.clus_blk.id", cid );
  C->SetBranchAddress( "sbs.hcal.clus_blk.e", ce );
  C->SetBranchAddress( "sbs.hcal.nblk", &nblk );
  C->SetBranchAddress( "sbs.hcal.nclus", &nclus );
  C->SetBranchAddress( "bb.tr.n", &BBtr_n );
  C->SetBranchAddress( "bb.tr.p", BBtr_p );
  C->SetBranchAddress( "bb.tr.px", BBtr_px );
  C->SetBranchAddress( "bb.tr.py", BBtr_py );
  C->SetBranchAddress( "bb.tr.pz", BBtr_pz );
  C->SetBranchAddress( "bb.tr.vz", BBtr_vz );
  C->SetBranchAddress( "bb.ps.x", &BBps_x );
  C->SetBranchAddress( "bb.ps.y", &BBps_y );
  C->SetBranchAddress( "bb.ps.e", &BBps_e );
  C->SetBranchAddress( "bb.sh.x", &BBsh_x );
  C->SetBranchAddress( "bb.sh.y", &BBsh_y );
  C->SetBranchAddress( "bb.sh.e", &BBsh_e );
  C->SetBranchAddress( "bb.tdctrig.tdc", TDCT_tdc );
  C->SetBranchAddress( "bb.tdctrig.tdcelemID", TDCT_id );
  C->SetBranchAddress( "Ndata.bb.tdctrig.tdcelemID", &TDCTndata );
  C->SetBranchAddress( "sbs.hcal.a_amp", hcal_amp );
  C->SetBranchAddress( "sbs.hcal.a_amp_p", hcal_amp_p );
  C->SetBranchAddress( "sbs.hcal.a_time", hcal_atime );
  C->SetBranchAddress( "sbs.hcal.tdc", hcal_tdc );
  C->SetBranchAddress( "sbs.hcal.tdcelemID", hcal_id );
  C->SetBranchAddress( "Ndata.sbs.hcal.tdcelemID", &hcalndata );
  C->SetBranchAddress( "e.kine.W2", &kineW2 );
  //C->SetBranchAddress( "bb.gem.track.nhits", &GEMtr_hits );
  //C->SetBranchAddress( "bb.etot_over_p", &BBeoverp );

  Int_t nentries = C->GetEntries();
  cout << endl;

  cout << "Opened tree with " << nentries << " entries." << endl;
  
  // Define some global variables
  Double_t emin = 0.0; // Minimum HCal energy to be considered as a hit, added for hard check on 1D distributions
  Double_t efficiency_rel; // Relative detection efficiency of HCAL (elastic events detected by HCal) / (elastic events as defined by BB tracking)
  Int_t hits_elBB = 0; // Count of total elastic hits detected in BB
  Int_t hits_gHCAL = 0; // Count of all elastic events that HCal detected
  Double_t pBeam = E_e/(1.+E_e/M_p*(1.-cos(BB_th))); // Momentum of beam calculated from theta

  //Mean energy loss of the beam prior to the scattering:
  Double_t Eloss_outgoing = celldiameter/2.0/sin(BB_th) * rho_tgt * dEdx_tgt; //Approximately 1 MeV, could correct further with raster position

  // Get elastic region of interest in HCal. No fermi smearing on LH2.
  Double_t dx0_p = 0.01089; // Position of proton spot, x-x_expected
  Double_t dy0_p = -0.02226; // Position of proton spot, y-y_expected
  Double_t dx_sig_p = 0.0878; // Max spread of proton spot, x-x_expected
  Double_t dy_sig_p = 0.1245; // Max spread of proton spot, y-y_expected

  //Keep a wide cut on W for allowed elastics
  //Double_t W2min_elastic = W2_mean - 2.*W2_sig;
  //Double_t W2max_elastic = W2_mean + 2.*W2_sig;

  //Hard-coded to use kinematic variable W2 for now (from parsed data files)
  //Double_t W2_mean = W2_mean;
  //Double_t W2_sig = W2_sig;
  Double_t W2min_elastic = W2_mean - 2.*W2_sig;
  Double_t W2max_elastic = W2_mean + 2.*W2_sig;

  //BigBite tracks
  Double_t ntrack;
  
  //Outfiles
  TFile *fout = new TFile(outputfilename,"RECREATE");
  //ofstream clusblk_bad_tdctime;
  //clusblk_bad_tdctime.open("clusblk_bad_tdctime.txt");

  //Histograms
  TH1D *hHCAL_ttime_nocut = new TH1D("hHCAL_ttime_nocut","Cluster TDC, primary blk, no cut; ns", 1200, -2000, 10000 );  
  //TH1D *hHCAL_ttime_activecut = new TH1D("hHCAL_ttime_activecut","Cluster TDC, primary blk, HCal active area cut; ns", 1200, -2000, 10000 );  
  //TH1D *hHCAL_ttime_gtzblkcut = new TH1D("hHCAL_ttime_gtzblkcut","Cluster TDC, primary blk, nblk>0 cut; ns", 1200, -2000, 10000 );
  //TH1D *hHCAL_ttime_ezblkcut = new TH1D("hHCAL_ttime_ezblkcut","Cluster TDC, primary blk, nblk==0 cut; ns", 1200, -2000, 10000 );  
  TH1D *hHCAL_ttime_BBcut = new TH1D("hHCAL_ttime_BBcut","Cluster TDC, primary blk, elas cut; ns", 1200, -2000, 10000 );
  TH1D *hHCAL_ttime_BBcut_badTDCallblks = new TH1D("hHCAL_ttime_BBcut_badTDCallblks","Cluster TDC, all blks, elas cut, bad clus tdc; ns", 1200, -2000, 10000 );
  TH1D *hHCAL_atime_nocut = new TH1D("hHCAL_atime_nocut","Cluster ADCt, primary blk, no cut; ns", 1000, -500, 500 );  
  TH1D *hHCAL_atime_BBcut = new TH1D("hHCAL_atime_BBcut","Cluster ADCt, primary blk, elas cut; ns", 1000, -500, 500 );  
  TH1D *hHCAL_atime_BBcut_badTDCallblks = new TH1D("hHCAL_atime_BBcut_badTDCallblks","Cluster ADCt, primary blk, elas cut, bad clus tdc; ns", 1000, -500, 500 );  
  TH1D *hHCAL_ttime_badTDCmaxdiffallblks = new TH1D("hHCAL_ttime_badTDCmaxdiffallblks","Cluster TDC, bad cluster tdc, tdc diff pblk - blk; ns", 1200, 0, 12000 );  
  TH1D *hHCAL_ttime_goodTDCmaxdiffallblks = new TH1D("hHCAL_ttime_goodTDCmaxdiffallblks","Cluster TDC, good cluster tdc, tdc diff pblk - blk; ns", 1200, 0, 12000 );  
  TH1D *hHCALF_ttime_nocut = new TH1D("hHCALF_ttime_nocut","All TDC, cluster id cut; ns", 1200, -2000, 10000 );  
  TH1D *hHCALF_ttime_BBcut = new TH1D("hHCALF_ttime_BBcut","All TDC, cluster id cut, elas cut; ns", 1200, -2000, 10000 );
  TH1D *hHCALF_ttime_BBcut_badTDCallblks = new TH1D("hHCALF_ttime_BBcut_badTDCallblks","All TDC, all blks, cluster id cut, elas cut, bad clus tdc; ns", 1200, -2000, 10000 );

  TH1D *hHCAL_diff = new TH1D("hHCAL_diff","TDC; ns", 2000, -1000, 1000 );  
  TH1D *hHCAL_nblk_nocut = new TH1D("hHCAL_nblk_nocut","Number of cluster elements; N", 20, 0, 20 );  
  TH1D *hHCAL_nblk_BBcut = new TH1D("hHCAL_nblk_BBcut","Number of cluster elements; N", 20, 0, 20 );  
  TH1D *hHCAL_time = new TH1D("hHCAL_time","TDC; ns", 1000, -500, 500 );  
  TH1D *hHCAL_avgtime = new TH1D("hHCAL_avgtime","TDC Average on all good tdc blocks; ns", 1000, -500, 500 );
  TH1D *hHCAL_gtime = new TH1D("hHCAL_gtime","TDC on highest E block with data; ns", 1000, -500, 500 );
  TH1D *hHCAL_allblk_gtdc = new TH1D("hHCAL_allblk_gtdc","TDC values on all blocks, good tdc time; ns", 1000, -500, 500 );
  TH1D *hW2_nocut = new TH1D( "hW2_nocut",";W^2 (GeV);", 420, 0.0, 1.4 );
  TH1D *hW2_BBcut = new TH1D( "hW2_BBcut",";W^2 (GeV);", 420, 0.0, 1.4 );
  TH1D *hW2_BBcut_norange = new TH1D( "hW2_BBcut_norange",";W^2 (GeV);", 420, 0.0, 1.4 );
  TH1D *hW2_BBcut_HCalcut = new TH1D( "hW2_BBcut_HCalcut",";W^2 (GeV);", 420, 0.0, 1.4 );
  TH1D *htDiff_BBcut = new TH1D( "htDiff_BBcut",";TDC_{HCal}-TDC_{BBCal} (ns);", 1300, -500, 800 );
  TH1D *htDiff_nocut = new TH1D( "htDiff_nocut",";TDC_{HCal}-TDC_{BBCal} (ns);", 1300, -500, 800 );
  TH2D *hdxdy_HCAL_BBcut = new TH2D("hdxdy_HCAL_BBcut",";y_{HCAL} (m); x_{HCAL} (m)", 250, posHCalYi_MC-2*HCalblk_l_h_MC, posHCalYf_MC+2*HCalblk_l_h_MC, 250, posHCalXi_MC-2*HCalblk_l_v_MC, posHCalXf_MC+2*HCalblk_l_v_MC );
  TH2D *hdxdy_HCAL_nocut = new TH2D("hdxdy_HCAL_nocut",";y_{HCAL} (m); x_{HCAL} (m)", 250, posHCalYi_MC-2*HCalblk_l_h_MC, posHCalYf_MC+2*HCalblk_l_h_MC, 250, posHCalXi_MC-2*HCalblk_l_v_MC, posHCalXf_MC+2*HCalblk_l_v_MC );
  TH1D *hdx_HCAL_nocut = new TH1D("hdx_HCAL_nocut",";x_{HCAL}-x_{expect} (m);", 500, posHCalXi_MC-2*HCalblk_l_v_MC, posHCalXf_MC+2*HCalblk_l_v_MC);
  TH1D *hdy_HCAL_nocut = new TH1D("hdy_HCAL_nocut",";y_{HCAL}-y_{expect} (m);", 500, posHCalYi_MC-2*HCalblk_l_h_MC, posHCalYf_MC+2*HCalblk_l_h_MC);
  TH1D *hdx_HCAL_BBcut = new TH1D("hdx_HCAL_BBcut",";x_{HCAL}-x_{expect} (m);", 500, posHCalXi_MC-2*HCalblk_l_v_MC, posHCalXf_MC+2*HCalblk_l_v_MC);
  TH1D *hdy_HCAL_BBcut = new TH1D("hdy_HCAL_BBcut",";y_{HCAL}-y_{expect} (m);", 500, posHCalYi_MC-2*HCalblk_l_h_MC, posHCalYf_MC+2*HCalblk_l_h_MC);
  TH1D *hdx_HCAL_BBcut_dycut = new TH1D("hdx_HCAL_BBcut_dycut",";x_{HCAL}-x_{expect} (m);", 500, posHCalXi_MC-2*HCalblk_l_v_MC, posHCalXf_MC+2*HCalblk_l_v_MC);

  TH1D *hamp_BBcut = new TH1D("hamp_nocut","Cluster amp, primary blk, atime>0, elas cut; mV",160,0,320);
  TH1D *hamp_p_BBcut = new TH1D("hamp_p_nocut","Cluster amp_p, primary blk, atime>0, elas cut; mV",160,0,320);
  TH1D *hamp_badtdccut = new TH1D("hamp_badtdccut","Cluster amp, primary blk, atime>0, BBcut, bad tdc; mV",160,0,320);
  TH1D *hamp_p_badtdccut = new TH1D("hamp_p_badtdccut","Cluster amp_p, primary blk, atime>0, BBcut, bad tdc; mV",160,0,320);
  TH1D *hblkID = new TH1D("hblkID","Primary Block ID", 300, -10, 290);


  Int_t maxblock = 12;
  TH1D *httime_blkN[maxblock];
  for( Int_t b=1; b<=maxblock; b++ ){
    httime_blkN[b-1] = new TH1D(Form("httime_blkN_%d",b),Form("Cluster TDC, blk %d, BB cut, Active Cut, Expect Cut; ns",b), 1200, -2000, 10000 );  
  }

  //Accounting and diagnostic variables
  Long64_t nevent = 0;
  Long64_t nbadclus = 0;
  Long64_t ngoodtdc = 0;
  Long64_t nexp = 0; //Number of expected events from BB and active area cuts
  Long64_t npdet = 0; //Number of detected events using only the primary blk in each cluster
  Long64_t nrdet = 0; //Number of detected events using highest E blk with tdc data 

  cout << "Processing events.." << endl;

  Double_t progress = 0.;
  
  while( progress<1.0 ){
    
    Int_t barwidth = 70;
    Int_t step = 1;
    
    while( C->GetEntry( nevent++ ) ){
      
      //Correct the beam energy with energy loss in target using vertex position
      Double_t Eloss = (BBtr_vz[0]+l_tgt/2.0) * rho_tgt * dEdx_tgt + uwallthick_LH2 * rho_Al * dEdx_Al; //approximately 3 MeV
      Double_t E_corr = E_e - Eloss;

      //Physics Calculations
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
      TVector3 HCALpos = HCAL_origin + hcalx * HCAL_xaxis + hcaly * HCAL_yaxis;
      Double_t sintersect = ( HCAL_origin - vertex ).Dot( HCAL_zaxis ) / (pNhat.Dot( HCAL_zaxis ) ); 
      TVector3 HCAL_intersect = vertex + sintersect * pNhat;   

      //Calculate quantities of interest
      Double_t yexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_yaxis );
      Double_t xexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_xaxis ); 
      Double_t dx = hcalx - xexpect_HCAL;
      Double_t dy = hcaly - yexpect_HCAL;
      //Double_t dr = sqrt( pow( dx, 2 )+pow( dy, 2 ));
      Double_t E_ep = sqrt( pow(M_e,2) + pow(BBtr_p[0],2) ); // Obtain the scattered electron energy 
      Double_t p_ep = BBtr_p[0];     
      Double_t Q2 = 2*E_corr*E_ep*( 1-(BBtr_pz[0]/p_ep) ); // Obtain Q2 from beam energy, outgoing electron energy, and momenta  
	
      //Double_t W2 = kineW2; //Get invarient mass transfer W from the four-momentum of the scattered nucleon
      Double_t W2 = pow(PgammaN.M(),2);

      //BBCal and HCal trigger coincidence
      Double_t bbcal_time=0., hcal_time=0.;
      for(Int_t ihit=0; ihit<TDCTndata; ihit++){
	if(TDCT_id[ihit]==5) bbcal_time=TDCT_tdc[ihit];
	if(TDCT_id[ihit]==0) hcal_time=TDCT_tdc[ihit];
      }
      Double_t diff = hcal_time - bbcal_time; 

      hHCAL_diff->Fill( diff );

      //Record best tdc time per event
      //Double_t hcal_gtime = 400.;
      Double_t hcal_rtime = ctdc[0];
      Double_t tdcmaxdiff = 0;
      Int_t pblkid = (int)(cid[0]-1); //Must correct for index switch on blocks
      Double_t pblktdc = hcal_tdc[pblkid]; //Correction between clus_blk vars and tdc
      Double_t pblkamp = hcal_amp[pblkid];
      Double_t pblkamp_p = hcal_amp_p[pblkid];
      Double_t pblkatime = hcal_atime[pblkid];

      if( pblktdc > 10000 ) pblktdc = -1500.; //Get value into window for comparison

      //Check blk index
      hblkID->Fill(pblkid);
	
      //Fill histograms without any cuts
      hW2_nocut->Fill( W2 );
      htDiff_nocut->Fill( diff );
      hdxdy_HCAL_nocut->Fill( dy, dx );
      hdx_HCAL_nocut->Fill( dx );
      hdy_HCAL_nocut->Fill( dy );
      hHCAL_nblk_nocut->Fill( nblk );
      hHCAL_atime_nocut->Fill( catime[0] );
      hHCAL_ttime_nocut->Fill( ctdc[0] );
      hHCALF_ttime_nocut->Fill( pblktdc );
      //if( nblk>0 ) hHCAL_ttime_gtzblkcut->Fill( ctdc[0] );
      //if( nblk==0 ) hHCAL_ttime_ezblkcut->Fill( ctdc[0] );

      //perform w2 cut to tighten elastic selection
      bool elastic = abs( W2 - W2_mean ) < W2_sig;
	
      bool exphcalactive = false;
      if(  yexpect_HCAL > (posHCalYf_MC-HCalblk_l_h_MC) &&
	   yexpect_HCAL < (posHCalYi_MC+HCalblk_l_h_MC) &&
	   xexpect_HCAL > (posHCalXf_MC-HCalblk_l_v_MC) &&
	   xexpect_HCAL < (posHCalXi_MC+HCalblk_l_v_MC)){
	  
	exphcalactive=true;
      }
      bool hcalactive = false;
      if( crow[0] > 0 &&
	  crow[0] < 23 &&
	  ccol[0] > 0 &&
	  ccol[0] < 11 ){
	hcalactive = true;
      }

      //Cut all events where no good tdc is expected in hcal and a cluster exists in HCal from fADC (BBcut)
      if( elastic == false || 
	  exphcalactive == false || 
	  hcalactive == false || 
	  kine==0 ) 
	continue;
      
      nexp++;
      
      //Fill histograms with BigBite and HCal active area cuts
      hW2_BBcut_norange->Fill( W2 );
      hW2_BBcut->Fill( W2 );
      htDiff_BBcut->Fill( diff );
      hdxdy_HCAL_BBcut->Fill( dy, dx );
      hdx_HCAL_BBcut->Fill( dx );
      hdy_HCAL_BBcut->Fill( dy );
      hHCAL_nblk_BBcut->Fill( nblk );
      hHCAL_atime_BBcut->Fill( catime[0] );
      hHCAL_ttime_BBcut->Fill( ctdc[0] );
      hHCALF_ttime_BBcut->Fill( pblktdc );
      if( pblkatime>0 ) {
	hamp_BBcut->Fill( pblkamp );
	hamp_p_BBcut->Fill( pblkamp_p );
      }
      

      //Get best possible tdc time from the cluster blocks

      //If the primary element of the cluster has a bad TDC time, record the adc time and write elements to txt
      if( ctdc[0]<-500||ctdc[0]>500 ){
	
	if( pblkatime>0 ) {
	  hamp_badtdccut->Fill( pblkamp );
	  hamp_p_badtdccut->Fill( pblkamp_p );
	}
	
	hcal_rtime = 400.; //Send all bad times to the same bin
	
	for( Int_t b=0; b<nblk; b++ ){
	  hHCAL_ttime_BBcut_badTDCallblks->Fill( ctdc[b] );
	  hHCAL_atime_BBcut_badTDCallblks->Fill( catime[b] );
	  //clusblk_bad_tdctime << ctdc[b] << endl;
	}
	if( nblk>1 ){ //Check if a good tdc time can be recovered
	  for( Int_t b=0; b<nblk; b++ ){
	    if( tdcmaxdiff<abs(ctdc[0]-ctdc[b]) ) tdcmaxdiff=abs(ctdc[0]-ctdc[b]);
	  }
	}
	hHCAL_ttime_badTDCmaxdiffallblks->Fill(tdcmaxdiff);
	
      }else{
	if( nblk>1 ){ //Check if a good tdc time can be recovered
	  for( Int_t b=0; b<nblk; b++ ){
	    if( tdcmaxdiff<abs(ctdc[0]-ctdc[b]) ) tdcmaxdiff=abs(ctdc[0]-ctdc[b]);
	  }
	}
	hHCAL_ttime_goodTDCmaxdiffallblks->Fill(tdcmaxdiff);
	
      }
      
      //Check average tdc on all blocks and compare to tdc time for highest E blk available
      //On sbs4/sbs11 it appears that the highest E blk has the tighter distribution
      vector<Double_t> avgtimes;
      Double_t avgtimes_sum = 0.;
      Double_t hcal_gtime = -2000.;
      //Double_t hcal_avgtime = -499.;
      for( Int_t b=0; b<nblk; b++ ){ //Blocks are sorted in descending energy
	if( ctdc[b]>-500&&ctdc[b]<500 ){
	  if( hcal_gtime == -2000. ) hcal_gtime = ctdc[b];
	  avgtimes.push_back(ctdc[b]);
	  hHCAL_allblk_gtdc->Fill(ctdc[b]);
	}
      }
      for( Int_t a=0; a<avgtimes.size(); a++){
	if( avgtimes.size()<1 ){
	  //hcal_avgtime=-499.;
	  continue;
	}
	avgtimes_sum += avgtimes[a];
      }
      Double_t hcal_avgtime = avgtimes_sum / avgtimes.size();
      if( avgtimes.size()==1 ) hcal_avgtime=-499.;
      
      hHCAL_avgtime->Fill(hcal_avgtime);
      hHCAL_gtime->Fill(hcal_gtime);
      hHCAL_time->Fill(hcal_rtime);
      
      //Fill dx with dy cut
      if( abs(dy-dy0_p)<dy_sig_p ){
	hdx_HCAL_BBcut_dycut->Fill( dx );
	if( abs(dx-dx0_p)<dx_sig_p ){
	  hW2_BBcut_HCalcut->Fill( W2 );
	}
      }
      
      //hHCAL_ttime_activecut->Fill( ctdc[0] );
      
      //Determine distributions at each block
      for( Int_t b=0; b<nblk; b++ ){ //Blocks are sorted in descending energy
	httime_blkN[b]->Fill(ctdc[b]);
      }
      
      //Obtain tdc efficiencies
      if( ctdc[0]<-500||ctdc[0]>500 ){
	for( Int_t b=0; b<nblk; b++ ){ //Blocks are sorted in descending energy
	  if( ctdc[b]>-500&&ctdc[b]<500 ){
	    nrdet++;
	  }
	}
	
      }else{
	npdet++;
      }
      
    }
    
    cout << "[";
    Int_t pos = barwidth*progress;
    for( Int_t i=0; i<barwidth; ++i ){
      if( i<pos ) cout << " ";
      else if( i==pos ){ 
	
	if( step%4==0 ){
	  cout << "0_0";
	}
	if( step%4==1 ){
	  cout << ">_>";
	}
	if( step%4==2 ){
	  cout << "<_<";
	}
	if( step%4==3 ){
	  cout << "-_-";
	}
	
      }
      else cout << " ";
    }
    progress = (Double_t)( ( nevent+1. )/nentries );
    
    cout << "]" << Int_t( progress*100 ) << "%\r";
    cout.flush();
    if( nevent%10000==0 ) step++;
    
    if( nevent >= nentries ) break;


  }
  
  

  
  cout << "TDC efficiency all blks primary blk only: " << (Double_t)npdet/(Double_t)nexp << endl << endl;

  cout << "TDC efficiency all blks highest E available blk: " << (Double_t)nrdet/(Double_t)nexp << endl << endl;



  fout->Write();

  st->Stop();

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

}
