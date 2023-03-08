//SSeeds 5.7.22 - Post-production - Code to parse several runs back by making cuts on elastics and produce a single root file for further analysis. Must take care to attend to all setup file parameters for proper configuration of the the output and one should not mix LH2 with LD2 or different magnet settings when stringing several files together.

#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <map>
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

//Static Detector Parameters
const int maxTracks = 1000; // Reasonable limit on tracks to be stored per event
const int maxTdcChan = 10; // Set to accomodate original 5 TDCTrig channels with buffer
const double hcalheight = 0.365; // Height of HCal above beamline

//Constants
const double PI = TMath::Pi();
const double M_e = 0.00051;
const double M_p = 0.938272;
const double M_n = 0.939565;

//Static Target Parameters
const double l_tgt = 0.15; // Length of the target (m)
const double rho_tgt = 0.0723; // Density of target (g/cc)
const double rho_Al = 2.7; // Density of aluminum windows (g/cc)
const double celldiameter = 1.6*2.54; //cm, right now this is a guess
const double Ztgt = 1.0;
const double Atgt = 1.0;
const double Mmol_tgt = 1.008; //g/mol

//For energy-loss correction to beam energy:
const double dEdx_tgt=0.00574; //According to NIST ESTAR, the collisional stopping power of hydrogen is about 5.74 MeV*cm2/g at 2 GeV energy
const double dEdx_Al = 0.0021; //According to NIST ESTAR, the collisional stopping power of Aluminum is about 2.1 MeV*cm2/g between 1-4 GeV
const double uwallthick_LH2 = 0.0145; //cm
const double dwallthick_LH2 = 0.015; //cm
const double cellthick_LH2 = 0.02; //cm, this is a guess;
const double Alshieldthick = 2.54/8.0; //= 1/8 inch * 2.54 cm/inch  

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

void pelas( const char *configfilename="setup_pelas.cfg", const char *outputfilename="outfiles/pelas_out_11547_0.root" ){
  
  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  // Get the date
  string date = getDate();
  
  // Declare Chain for many root files
  TChain *C = new TChain("T");

  // Declare general physics parameters to be modified by input config file
  int kine = 8; // Keep track of kinematic being analyzed
  double E_e = 1.916; // Energy of beam (incoming electrons from accelerator)
  double BB_d = 1.7988; // Distance to bigbite spectrometer from target chamber (m)
  double BB_th = 36.; // Angle BB spectrometer makes with exit beamline
  double HCal_d = 14.5; // Distance to HCal from scattering chamber for comm1
  double HCal_th = 35.; // Angle HCal center makes with exit beamline  
  double W_mean = 0.93; // Mean of W at current kinematic
  double W_sig = 0.039; // Width of W at current kinematic
  double dx0_n = 0.9; // Position of neutron spot, x-x_expected
  double dy0_n = 0.62; // Position of neutron spot, y-y_expected
  double dx_sig_n = 0.09; // Max spread of neutron spot, x-x_expected
  double dy_sig_n = 0.15; // Max spread of neutron spot, y-y_expected
  double dx0_p = 0.9; // Position of proton spot, x-x_expected
  double dy0_p = 0.62; // Position of proton spot, y-y_expected
  double dx_sig_p = 0.09; // Max spread of proton spot, x-x_expected
  double dy_sig_p = 0.15; // Max spread of proton spot, y-y_expected
  double dx_max = 0.65; // Max distance expected for p-n separation
  int magSet = 30;  // SBS Magnetic field strength (percent)
  int useAlshield = 0; //Use 1/8" aluminum shield on scattering chamber exit?
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
    int ntokens = tokens->GetEntries();
    if( ntokens>1 ){
      TString skey = ( (TObjString*)(*tokens)[0] )->GetString();
      if( skey == "kine" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	kine = sval.Atoi();
	cout << "kinematic setting: " << kine << endl;
      }
      if( skey == "E_e" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	E_e = sval.Atof();
	cout << "beam energy: " << E_e << endl;
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
      if( skey == "dx0_n" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	dx0_n = sval.Atof();
	cout << "x position of neutron spot: " << dx0_n << endl;
      }
      if( skey == "dy0_n" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	dy0_n = sval.Atof();
	cout << "y position of neutron spot: " << dy0_n << endl;
      }
      if( skey == "dx_sig_n" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	dx_sig_n = sval.Atof();
	cout << "x sigma of neutron spot: " << dx_sig_n << endl;
      }
      if( skey == "dy_sig_n" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	dy_sig_n = sval.Atof();
	cout << "y sigma of neutron spot: " << dy_sig_n << endl;
      }
      if( skey == "dx0_p" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	dx0_p = sval.Atof();
	cout << "x position of proton spot: " << dx0_p << endl;
      }
      if( skey == "dy0_p" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	dy0_p = sval.Atof();
	cout << "y position of proton spot: " << dy0_p << endl;
      }
      if( skey == "dx_sig_p" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	dx_sig_p = sval.Atof();
	cout << "x sigma of proton spot: " << dx_sig_p << endl;
      }
      if( skey == "dy_sig_p" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	dy_sig_p = sval.Atof();
	cout << "y sigma of proton spot: " << dy_sig_p << endl;
      }
      if( skey == "dx_max" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	dx_max = sval.Atof();
	cout << "vertical max p-n displacement: " << dx_max << endl;
      }
      if( skey == "magSet" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	magSet = sval.Atoi();
	cout << "SBS magnet field strength: " << magSet << endl;
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
  
  cout << "Setup parameters loaded." << endl << endl;

  TEventList *elist = new TEventList("elist","Elastic Event List");
  C->Draw(">>elist",globalcut);

  cout << "Event list populated with cut placed on elastics." << endl;

  // Declare general tree parameters
  double BBtr_p[maxTracks], BBtr_px[maxTracks], BBtr_py[maxTracks], BBtr_pz[maxTracks];
  double BBtr_vx[maxTracks], BBtr_vy[maxTracks], BBtr_vz[maxTracks], BBtr_chi2[maxTracks];
  double BBfp_x[maxTracks], BBfp_y[maxTracks], BBfp_th[maxTracks], BBfp_ph[maxTracks];
  double BBtgt_x[maxTracks], BBtgt_y[maxTracks], BBtgt_th[maxTracks], BBtgt_ph[maxTracks];
  double BBtr_n, BBps_x, BBps_y, BBps_e, BBsh_x, BBsh_y, BBsh_e;

  double TDCT_id[maxTdcChan], TDCT_tdc[maxTdcChan]; 
  int TDCTndata;

  double HCALx, HCALy, HCALe;

  double BLdnew, BLclk;

  // Declare root tree variables and set values to memory locations in root file
  // Switch them on
  C->SetBranchStatus( "*", 0 );
  // HCal
  C->SetBranchStatus( "sbs.hcal.x", 1 );
  C->SetBranchStatus( "sbs.hcal.y", 1 );
  C->SetBranchStatus( "sbs.hcal.e", 1 );
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
  // Trigger TDC
  C->SetBranchStatus( "bb.tdctrig.tdc", 1 );
  C->SetBranchStatus( "bb.tdctrig.tdcelemID", 1 );
  C->SetBranchStatus( "Ndata.bb.tdctrig.tdcelemID", 1 );
  // Beamline
  C->SetBranchStatus( "sbs.bcm.dnew.cnt", 1 );
  C->SetBranchStatus( "sbs.104kHz_CLK.cnt", 1 );

  // Set the variables
  // HCal
  C->SetBranchAddress( "sbs.hcal.x", &HCALx );
  C->SetBranchAddress( "sbs.hcal.y", &HCALy );
  C->SetBranchAddress( "sbs.hcal.e", &HCALe );
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
  // Trigger TDC
  C->SetBranchAddress( "bb.tdctrig.tdcelemID", TDCT_id );
  C->SetBranchAddress( "bb.tdctrig.tdc", TDCT_tdc );
  C->SetBranchAddress( "Ndata.bb.tdctrig.tdcelemID", &TDCTndata );  
  // Beamline
  C->SetBranchAddress( "sbs.bcm.dnew.cnt", &BLdnew );
  C->SetBranchAddress( "sbs.104kHz_CLK.cnt", &BLclk );

  cout << "Tree variables linked." << endl;

  // Define a clock to check macro processing time
  TStopwatch *sw = new TStopwatch();
  sw->Start();
  
  // Declare outfiles
  outputfilename = Form( "/lustre19/expphy/volatile/halla/sbs/seeds/GMN_parsed/e1209019_parsed_SBS%d_tar%s_%s.root", magSet, tar.c_str(), date.c_str() );
  TFile *fout = new TFile( outputfilename, "RECREATE" );
  string logpath = Form( "/lustre19/expphy/volatile/halla/sbs/seeds/GMN_parsed/logs/parsedLog_%s.txt", date.c_str() );

  // Initialize misc. variables
  int elasYield = 0; // Keep track of total elastics analyzed
  double pBeam = E_e/(1.+E_e/M_p*(1.-cos(BB_th)));  

  //Mean energy loss of the beam prior to the scattering:
  double Eloss_outgoing = celldiameter/2.0/sin(BB_th) * rho_tgt * dEdx_tgt; //Approximately 1 MeV, could correct further with raster position
  if( useAlshield != 0 ) Eloss_outgoing += Alshieldthick * rho_Al * dEdx_Al;

  TH1D *h_W2 = new TH1D("h_W2",";W2 (GeV^2);",250,0,2);
  TH1D *h_W = new TH1D("h_W",";W (GeV);",250,0,2);

  TTree *Tout = new TTree("Tout","Tree containing variables after final cuts placed on W2 and HCal");
  double T_ebeam, T_etheta, T_ephi, T_precon, T_pelastic, T_thetabend, T_dpel, T_W2, T_W;
  double T_pincident;
  double T_BBtr_vx, T_BBtr_vy, T_BBtr_vz, T_BBtr_px, T_BBtr_py, T_BBtr_pz;
  double T_BBtr_n, T_BBtr_p;
  double T_HCALx, T_HCALy, T_HCALe;
  double T_pp_expect, T_ptheta_expect, T_pphi_expect;
  double T_Q2;
  double T_BBps_x, T_BBps_y, T_BBps_e, T_BBsh_x, T_BBsh_y, T_BBsh_e;
  double T_BB_d, T_BB_th, T_HCal_d, T_HCal_th;

  Tout->Branch( "Ebeam", &T_ebeam, "Ebeam/D" );
  Tout->Branch( "Q2", &T_Q2, "Q2/D"); 
  Tout->Branch( "etheta", &T_etheta, "etheta/D");
  Tout->Branch( "ephi", &T_ephi, "ephi/D");
  Tout->Branch( "ep_recon", &T_precon, "ep_recon/D");
  Tout->Branch( "ep_elastic", &T_pelastic, "ep_elastic/D");
  Tout->Branch( "ep_incident", &T_pincident, "ep_incident/D");
  Tout->Branch( "thetabend", &T_thetabend, "thetabend/D");
  Tout->Branch( "dpel", &T_dpel, "dpel/D");
  Tout->Branch( "W2", &T_W2, "W2/D");
  Tout->Branch( "W", &T_W, "W/D");
  Tout->Branch( "BBtr_vx", &T_BBtr_vx, "BBtr_vx/D");
  Tout->Branch( "BBtr_vy", &T_BBtr_vy, "BBtr_vy/D");
  Tout->Branch( "BBtr_vz", &T_BBtr_vz, "BBtr_vz/D");
  Tout->Branch( "BBtr_px", &T_BBtr_px, "BBtr_px/D");
  Tout->Branch( "BBtr_py", &T_BBtr_py, "BBtr_py/D");
  Tout->Branch( "BBtr_pz", &T_BBtr_pz, "BBtr_pz/D");
  Tout->Branch( "BBtr_p", &T_BBtr_p, "BBtr_p/D");
  Tout->Branch( "BBtr_n", &T_BBtr_n, "BBtr_n/D");
  Tout->Branch( "BB_d", &T_BB_d, "BB_d/D");
  Tout->Branch( "BB_th", &T_BB_th, "BB_th/D");
  Tout->Branch( "HCal_d", &T_HCal_d, "HCal_d/D");
  Tout->Branch( "HCal_th", &T_HCal_th, "HCal_th/D");
  Tout->Branch( "HCALx", &T_HCALx, "HCALx/D");
  Tout->Branch( "HCALy", &T_HCALy, "HCALy/D");
  Tout->Branch( "HCALe", &T_HCALe, "HCALe/D");
  Tout->Branch( "pp_expect", &T_pp_expect, "pp_expect/D");
  Tout->Branch( "ptheta_expect", &T_ptheta_expect, "ptheta_expect/D");
  Tout->Branch( "pphi_expect", &T_pphi_expect, "pphi_expect/D");
  Tout->Branch( "BBps_x", &T_BBps_x, "BBps_x/D");
  Tout->Branch( "BBps_y", &T_BBps_y, "BBps_y/D");
  Tout->Branch( "BBps_e", &T_BBps_e, "BBps_e/D");
  Tout->Branch( "BBsh_x", &T_BBsh_x, "BBsh_x/D");
  Tout->Branch( "BBsh_y", &T_BBsh_y, "BBsh_y/D");
  Tout->Branch( "BBsh_e", &T_BBsh_e, "BBsh_e/D");

  // Set long int to keep track of total entries
  Long64_t Nevents = elist->GetN();

  cout << endl << "Initialization complete." << endl << endl;
  cout << "Opened up tree with nentries: " << C->GetEntries() << ", nentries passing globalcut: " << Nevents << "." << endl << endl;

  //Loop over events
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
      cout << "] " << int(progress*100.) << "%, elastic events: " << elasYield << ", time remaining: " << int(timeremains/60.) << "m \r";
      cout.flush();
      
      C->GetEntry( elist->GetEntry( nevent ) ); 

      //Cut on BBCal and HCal trigger coincidence
      double bbcal_time=0., hcal_time=0.;
      for(int ihit=0; ihit<TDCTndata; ihit++){
	if(TDCT_id[ihit]==5) bbcal_time=TDCT_tdc[ihit];
	if(TDCT_id[ihit]==0) hcal_time=TDCT_tdc[ihit];
      }
      double diff = hcal_time - bbcal_time; 

      //Cut here on BB trigger/HCal trigger coincidence time
      if( fabs(diff-510.)>10. ) continue;

      //Correct the beam energy with energy loss in target using vertex position
      double Eloss = (BBtr_vz[0]+l_tgt/2.0) * rho_tgt * dEdx_tgt + uwallthick_LH2 * rho_Al * dEdx_Al; //approximately 3 MeV
      double E_corr = E_e - Eloss;

      double p_corr = BBtr_p[0] - Eloss_outgoing; //Neglecting the mass of e'

      //Proceed only if at least one track exists in BB arm - lowest chi2 track always first element
      double etheta = acos( BBtr_pz[0]/BBtr_p[0] ); //Will use the uncorrected track momentum to reconstruct e' theta
      double ephi = atan2( BBtr_py[0], BBtr_px[0] );
      
      TVector3 vertex( 0, 0, BBtr_vz[0] ); // z location of vertex in hall coordinates
      TLorentzVector Pbeam( 0, 0, E_corr, E_corr ); //Mass of e negligable
      TLorentzVector kprime( BBtr_px[0], BBtr_py[0], BBtr_pz[0], BBtr_p[0] );
      TLorentzVector Ptarg( 0, 0, 0, M_p );
      
      TLorentzVector q = Pbeam - kprime;
      TLorentzVector PgammaN = Ptarg + q; //(-px, -py, ebeam - pz, Mp + ebeam - p)
      
      double pel = E_corr/( 1.+E_corr/M_p*( 1.-cos(etheta) ) );
      double nu = E_corr - BBtr_p[0];
      double pp = sqrt( pow(nu,2)+2.*M_p*nu );
      double phinucleon = ephi + PI; //assume coplanarity
      double thetanucleon = acos( (E_corr - BBtr_pz[0])/pp ); //use elastic constraint on nucleon kinematics
      
      TVector3 pNhat( sin(thetanucleon)*cos(phinucleon), sin(thetanucleon)*sin(phinucleon), cos(thetanucleon) );
      
      //Define HCal coordinate system
      TVector3 HCAL_zaxis( sin(-HCal_th), 0, cos(-HCal_th) );
      TVector3 HCAL_xaxis( 0, -1, 0 );
      TVector3 HCAL_yaxis = HCAL_zaxis.Cross(HCAL_xaxis).Unit();
      
      TVector3 HCAL_origin = HCal_d * HCAL_zaxis + hcalheight * HCAL_xaxis;
      
      //Define intersection points for hadron vector
      double sintersect = ( HCAL_origin - vertex ).Dot( HCAL_zaxis ) / (pNhat.Dot( HCAL_zaxis ) );
      
      TVector3 HCAL_intersect = vertex + sintersect * pNhat;
      
      //Define the expected position of hadron on HCal from BB track
      double yexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_yaxis );
      double xexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_xaxis );
      
      double E_ep = sqrt( pow(M_e,2) + pow(BBtr_p[0],2) ); // Obtain the scattered electron energy
      
      double p_ep = BBtr_p[0];
      
      double Q2 = 2*E_corr*E_ep*( 1-(BBtr_pz[0]/p_ep) ); // Obtain Q2 from beam energy, outgoing electron energy, and momenta
      
      //Get invarient mass transfer W from the four-momentum of the scattered nucleon
      double W = PgammaN.M();
      
      //Use the electron kinematics to predict the proton momentum assuming elastic scattering on free proton at rest (will need to correct for fermi motion):
      double E_pp = nu+M_p; // Get energy of the proton
      double Enucleon = sqrt(pow(pp,2)+pow(M_p,2)); // Check on E_pp, same
      //hE_pp->Fill( E_pp ); // Fill histogram
      
      double KE_p = nu; // For elastics
      //hKE_p->Fill( KE_p );

      double dx = HCALx - xexpect_HCAL;
      double dy = HCALy - yexpect_HCAL;
      
      
      //Check "elastic" events on center HCal for id with spot checks
      bool HCAL_on = false, is_p = false, is_n = false;
      if( HCALy>-0.75 && HCALy<0.75 && HCALx>-2.015 && HCALx<1.285 ){
	HCAL_on = true;
      }
      if( pow( (dx - dx0_p)/dx_sig_p,2) + pow( (dy - dy0_p)/dy_sig_p,2) <= pow(2.5,2) ){
	is_p = true;
      }
      if( pow( (dx - dx0_n)/dx_sig_n,2) + pow( (dy - dy0_n)/dy_sig_n,2) <= pow(2.5,2) ){
	is_n = true;
      }
      
      if( HCAL_on==true && is_n==true ){
	if( (HCALx-dx_max)>-2.015 ){
	  if( tar == "LD2" ){

	    T_BB_d = BB_d;
	    T_BB_th = BB_th;
	    T_HCal_d = HCal_d;
	    T_HCal_th = HCal_th;

	    T_BBtr_vz = BBtr_vz[0];
	    T_BBtr_vy = BBtr_vy[0];
	    T_BBtr_vx = BBtr_vx[0];
	    T_BBtr_pz = BBtr_pz[0];
	    T_BBtr_py = BBtr_py[0];
	    T_BBtr_px = BBtr_px[0];
	    T_BBtr_p = BBtr_p[0];

	    T_BBtr_n = BBtr_n;

	    T_BBps_x = BBps_x;
	    T_BBps_y = BBps_y;
	    T_BBps_e = BBps_e;
	    T_BBsh_x = BBsh_x;
	    T_BBsh_y = BBsh_y;
	    T_BBsh_e = BBsh_e;

	    T_HCALx = HCALx;
	    T_HCALy = HCALy;
	    T_HCALe = HCALe;

	    T_etheta = etheta;
	    T_ephi = ephi;

	    T_ebeam = E_corr;

	    double pelastic = E_corr/(1.+(E_corr/M_p)*(1.0-cos(etheta))); 
      
	    double precon = BBtr_p[0] + Eloss_outgoing; //reconstructed momentum, corrected for mean energy loss exiting the target (later we'll also want to include Al shielding on scattering chamber)

	    T_precon = precon;
	    T_pelastic = pelastic; 
	    T_pincident = pelastic - Eloss_outgoing;

	    double nu_recon = E_corr - precon;
	    double Q2recon = 2.0*E_corr*precon*(1.0-cos(etheta));
	    double W2recon = pow(M_p,2) + 2.0*M_p*nu_recon - Q2recon;

	    T_Q2 = Q2recon;
	    T_W2 = W2recon;
	    T_W = W;
	    h_W2->Fill(W2recon);
	    h_W->Fill(W);

	  }

	  elasYield++;
	  Tout->Fill();
	}
      }else if( HCAL_on==true && is_p==true ){
	if( (HCALx+dx_max)<1.285 ){

	  T_BB_d = BB_d;
	  T_BB_th = BB_th;
	  T_HCal_d = HCal_d;
	  T_HCal_th = HCal_th;

	  T_BBtr_vz = BBtr_vz[0];
	  T_BBtr_vy = BBtr_vy[0];
	  T_BBtr_vx = BBtr_vx[0];
	  T_BBtr_pz = BBtr_pz[0];
	  T_BBtr_py = BBtr_py[0];
	  T_BBtr_px = BBtr_px[0];
	  T_BBtr_p = BBtr_p[0];

	  T_BBtr_n = BBtr_n;

	  T_BBps_x = BBps_x;
	  T_BBps_y = BBps_y;
	  T_BBps_e = BBps_e;
	  T_BBsh_x = BBsh_x;
	  T_BBsh_y = BBsh_y;
	  T_BBsh_e = BBsh_e;

	  T_HCALx = HCALx;
	  T_HCALy = HCALy;
	  T_HCALe = HCALe;

	  T_etheta = etheta;
	  T_ephi = ephi;

	  T_ebeam = E_corr;

	  double pelastic = E_corr/(1.+(E_corr/M_p)*(1.0-cos(etheta))); 
      
	  double precon = BBtr_p[0] + Eloss_outgoing; //reconstructed momentum, corrected for mean energy loss exiting the target (later we'll also want to include Al shielding on scattering chamber)

	  T_precon = precon;
	  T_pelastic = pelastic; 
	  T_pincident = pelastic - Eloss_outgoing;

	  double nu_recon = E_corr - precon;
	  double Q2recon = 2.0*E_corr*precon*(1.0-cos(etheta));
	  double W2recon = pow(M_p,2) + 2.0*M_p*nu_recon - Q2recon;

	  T_Q2 = Q2recon;
	  T_W2 = W2recon;
	  T_W = W;
	  h_W2->Fill(W2recon);
	  h_W->Fill(W);

	  elasYield++;
	  Tout->Fill();
	}
      }
    }
  }
  
  cout << endl;
  
  ofstream logfile;
  logfile.open( logpath );
  logfile << "#Log of runs corresponding to root file with same date" << endl;
  for( int i=0; i<log.size(); i++ ){
    logfile << log[i] << endl;
  }

  logfile.close();

  cout << endl << endl << "Elastic yield for analyzed runs: " << elasYield << ". Total events analyzed: " << Nevents << "." << endl << endl;
  cout << "New root file generated with elastic cuts placed at " << outputfilename << endl;

  fout->Write();

  elist->Delete();
  st->Stop();

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

}
