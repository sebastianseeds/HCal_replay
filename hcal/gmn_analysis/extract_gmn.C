//SSeeds 4.1.22 - Post-production - Analysis code to extract first-order gmn from cross section ratio via LD2 data

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
#include "blind.h"

//Static Detector Parameters
const int maxTracks = 1000; // Reasonable limit on tracks to be stored per event
const int maxTdcChan = 10; // Set to accomodate original 5 TDCTrig channels with buffer
const double hcalheight = -0.2897; // Height of HCal above beamline

//Constants
const double PI = TMath::Pi();
const double M_e = 0.00051;
const double M_p = 0.938272;
const double M_n = 0.939565;

//Static Target Parameters
const double l_tgt = 0.15; // Length of the target (m)
const double rho_tgt = 0.0723; // Density of target (g/cc)
const double rho_Al = 2.7; // Density of aluminum windows (g/cc)

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

// Selective gaussian to subract inelastic background 30% SBS field
double BG_fit_30(double *x, double *par){
  double amp = par[0];
  double offset = par[1];
  double sigma = par[2];
  if (x[0]>-1.2 && x[0]<0.48) { 
    TF1::RejectPoint();
    return 0;
  }
  return amp*exp(-0.5*pow((x[0]-offset)/sigma,2.));
}

// Selective gaussian to subract inelastic background 30% SBS field
double BG_fit_50(double *x, double *par){
  double amp = par[0];
  double offset = par[1];
  double sigma = par[2];
  if (x[0]>-1.77 && x[0]<0.51){
    TF1::RejectPoint();
    return 0;
  }
  return amp*exp(-0.5*pow((x[0]-offset)/sigma,2.));
}

void extract_gmn( const char *configfilename="setup_extract_gmn.cfg", const char *outputfilename="outfiles/gmn_out.root" ){
  
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
    int ntokens = tokens->GetEntries();
    if( ntokens>1 ){
      TString skey = ( (TObjString*)(*tokens)[0] )->GetString();
      if( skey == "kine" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	kine = sval.Atoi();
	cout << "Loading kinematic setting: " << kine << endl;
      }
      if( skey == "E_e" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	E_e = sval.Atof();
	cout << "Loading beam energy: " << E_e << endl;
      }
      if( skey == "BB_d" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	BB_d = sval.Atof();
	cout << "Loading BigBite Spectrometer distance: " << BB_d << endl;
      }
      if( skey == "BB_th" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	BB_th = sval.Atof() * TMath::DegToRad();
	cout << "Loading BigBite Spectrometer angle (rad): " << BB_th << endl;
      }
      if( skey == "HCal_d" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	HCal_d = sval.Atof();
	cout << "Loading HCal distance: " << HCal_d << endl;
      }
      if( skey == "HCal_th" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	HCal_th = sval.Atof() * TMath::DegToRad();	
	cout << "Loading HCal angle (rad): " << HCal_th << endl;
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
      if( skey == "dx0_n" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	dx0_n = sval.Atof();
	cout << "Loading x position of neutron spot: " << dx0_n << endl;
      }
      if( skey == "dy0_n" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	dy0_n = sval.Atof();
	cout << "Loading y position of neutron spot: " << dy0_n << endl;
      }
      if( skey == "dx_sig_n" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	dx_sig_n = sval.Atof();
	cout << "Loading x sigma of neutron spot: " << dx_sig_n << endl;
      }
      if( skey == "dy_sig_n" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	dy_sig_n = sval.Atof();
	cout << "Loading y sigma of neutron spot: " << dy_sig_n << endl;
      }
      if( skey == "dx0_p" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	dx0_p = sval.Atof();
	cout << "Loading x position of proton spot: " << dx0_p << endl;
      }
      if( skey == "dy0_p" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	dy0_p = sval.Atof();
	cout << "Loading y position of proton spot: " << dy0_p << endl;
      }
      if( skey == "dx_sig_p" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	dx_sig_p = sval.Atof();
	cout << "Loading x sigma of proton spot: " << dx_sig_p << endl;
      }
      if( skey == "dy_sig_p" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	dy_sig_p = sval.Atof();
	cout << "Loading y sigma of proton spot: " << dy_sig_p << endl;
      }
      if( skey == "dx_max" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	dx_max = sval.Atof();
	cout << "Loading y max p-n displacement: " << dx_max << endl;
      }
      if( skey == "magSet" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	magSet = sval.Atof();
	cout << "Loading SBS magnet field strength: " << magSet << endl;
      }
    }
    delete tokens;
  }
  
  cout << "Setup parameters loaded." << endl;

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

  cout << "Tree variables linked." << endl;

  // Define a clock to check macro processing time
  TStopwatch *sw = new TStopwatch();
  sw->Start();
  
  // Declare outfile
  if( magSet==30 ) outputfilename = "outfiles/gmn_out_SBS30.root";
  if( magSet==50 ) outputfilename = "outfiles/gmn_out_SBS50.root";
  TFile *fout = new TFile( outputfilename, "RECREATE" );


  // Initialize misc. variables
  int elasYield = 0; // Keep track of total elastics analyzed
  double pBeam = E_e/(1.+E_e/M_p*(1.-cos(BB_th)));  

  // Initialize histograms	
  TH1D *hDiff = new TH1D( "hDiff","HCal time - BBCal time (ns)", 1300, -500, 800 );
  TH2D *hPAngleCorr_th = new TH2D( "hPAngCorr_th","BB theta vs HCal theta", 200, 0.55, 0.75, 300, 0.35, 0.65 );
  TH2D *hPAngleCorr_ph = new TH2D( "hPAngCorr_ph","BB phi vs HCal phi", 500, -0.4, 0.1, 500, 2.7, 3.2 );
  TH1D *hW = new TH1D( "W", "W", 250, 0.3, 1.8 );
  hW->GetXaxis()->SetTitle( "GeV" );
  TH1D *hW_cut = new TH1D( "W_cut", "W_cut", 250, 0.3, 1.8 );
  hW_cut->GetXaxis()->SetTitle( "GeV" );
  TH1D *hW_fcut = new TH1D( "W_fcut", "W_fcut", 250, 0.3, 1.8 );
  hW_fcut->GetXaxis()->SetTitle( "GeV" );
  TH1D *hQ2 = new TH1D( "Q2", "Q2", 250, 0.5, 3.0 );
  hQ2->GetXaxis()->SetTitle( "GeV" );
  TH1D *hE_ep = new TH1D( "Scattered Electron Energy","E_ep", 500, 0.0, E_e*1.5 ); 
  hE_ep->GetXaxis()->SetTitle( "GeV" );
  TH1D *hE_pp = new TH1D( "Scattered Proton Energy", "E_pp", 500, 0.0, E_e*1.5 );
  hE_pp->GetXaxis()->SetTitle( "GeV" );
  TH1D *hKE_p = new TH1D( "Scattered Proton Kinetic Energy", "KE_pp", 500, 0.0, E_e*1.5 );
  hKE_p->GetXaxis()->SetTitle( "GeV" );
  TH1D *hdx_HCAL = new TH1D("hdx_HCAL","; x_{HCAL} - x_{exp} (m);", 250, -2.5,2.5);
  TH1D *hdx_HCAL_cut = new TH1D("hdx_HCAL_cut","; x_{HCAL} - x_{exp} (m);", 250, -2.5,2.5);
  TH1D *hdx_HCAL_wcut = new TH1D("hdx_HCAL_wcut","; x_{HCAL} - x_{exp} (m);", 250, -2.5,2.5);
  TH1D *hdx_HCAL_fcut = new TH1D("hdx_HCAL_fcut","; x_{HCAL} - x_{exp} (m);", 250, -2.5,2.5);
  TH1D *hdy_HCAL = new TH1D("hdy_HCAL","; y_{HCAL} - y_{exp} (m);", 250, -1.25,1.25);
  TH1D *hdy_HCAL_cut = new TH1D("hdy_HCAL_cut","; y_{HCAL} - y_{exp} (m);", 250, -1.25,1.25);  
  TH2D *hdxdy_HCAL = new TH2D("hdxdy_HCAL",";y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)", 250, -1.25, 1.25, 250, -2.5, 2.5 );
  TH2D *hdxdy_HCAL_wcut = new TH2D("hdxdy_HCAL_wcut",";y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)", 250, -1.25, 1.25, 250, -2.5, 2.5 );
  TH2D *hdxdy_HCAL_cut = new TH2D("hdxdy_HCAL_cut",";y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)", 250, -1.25, 1.25, 250, -2.5, 2.5 );
  TH2D *hdxdy_HCAL_ncut = new TH2D("hdxdy_HCAL_ncut",";y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)", 250, -1.25, 1.25, 250, -2.5, 2.5 );
  TH2D *hdxdy_HCAL_pcut = new TH2D("hdxdy_HCAL_pcut",";y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)", 250, -1.25, 1.25, 250, -2.5, 2.5 );
  TH2D *hdxdy_HCAL_fcut = new TH2D("hdxdy_HCAL_fcut",";y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)", 250, -1.25, 1.25, 250, -2.5, 2.5 );
  TH2D *hxy_HCAL = new TH2D("hxy_HCAL",";y_{HCAL} (m); x_{HCAL} (m)",12,-0.9,0.9,24,-2.165,1.435);
  TH2D *hxy_HCAL_cut = new TH2D("hxy_HCAL_cut",";y_{HCAL} (m); x_{HCAL} (m)",12,-0.9,0.9,24,-2.165,1.435);
  //TH1D *hvz = new TH1D("hvz",";vertex z (m);", 250,-0.125,0.125);
  TH1D *hvz_cut = new TH1D("hvz_cut",";vertex z (m);", 250,-0.125,0.125);

  cout << "Variables and histograms defined." << endl;

  // Set long int to keep track of total entries
  Long64_t Nevents = elist->GetN();

  cout << endl << "All parameters loaded and initialization complete." << endl << endl;
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
      hDiff->Fill( diff );

      //Cut here on BB trigger/HCal trigger coincidence time
      if( fabs(diff-510.)>10. ) continue;

      //Sanity check
      if( (int)BBtr_n!=1 ) cout << endl << endl << "WARNING: Total tracks not as expected from global cut. Check globalcut for errors." << endl << endl;

      //Proceed only if at least one track exists in BB arm - lowest chi2 track always first element
      double etheta = acos( BBtr_pz[0]/BBtr_p[0] );
      double ephi = atan2( BBtr_py[0], BBtr_px[0] );
      
      TVector3 vertex( 0, 0, BBtr_vz[0] ); // z location of vertex in hall coordinates
      TLorentzVector Pbeam( 0, 0, E_e, E_e ); //Mass of e negligable
      TLorentzVector kprime( BBtr_px[0], BBtr_py[0], BBtr_pz[0], BBtr_p[0] );
      TLorentzVector Ptarg( 0, 0, 0, M_p );
      
      TLorentzVector q = Pbeam - kprime;
      TLorentzVector PgammaN = Ptarg + q; //(-px, -py, ebeam - pz, Mp + ebeam - p)
      
      double pel = E_e/( 1.+E_e/M_p*( 1.-cos(etheta) ) );
      double nu = E_e - BBtr_p[0];
      double pp = sqrt( pow(nu,2)+2.*M_p*nu );
      double phinucleon = ephi + PI; //assume coplanarity
      //double thetanucleon = acos( (E_e - BBtr_p[0]*cos(etheta))/pp ); //use elastic constraint on nucleon kinematics
      double thetanucleon = acos( (E_e - BBtr_pz[0])/pp ); //use elastic constraint on nucleon kinematics
      
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
      hE_ep->Fill( E_ep );
      
      double p_ep = BBtr_p[0];
      
      double Q2 = 2*E_e*E_ep*( 1-(BBtr_pz[0]/p_ep) ); // Obtain Q2 from beam energy, outgoing electron energy, and momenta
      hQ2->Fill( Q2 ); // Fill histogram
      
      //Get invarient mass transfer W from the four-momentum of the scattered nucleon
      double W = PgammaN.M();
      hW->Fill( W );
      
      //Use the electron kinematics to predict the proton momentum assuming elastic scattering on free proton at rest (will need to correct for fermi motion):
      double E_pp = nu+M_p; // Get energy of the proton
      double Enucleon = sqrt(pow(pp,2)+pow(M_p,2)); // Check on E_pp, same
      hE_pp->Fill( E_pp ); // Fill histogram
      
      double KE_p = nu; // For elastics
      hKE_p->Fill( KE_p );

      double dx = HCALx - xexpect_HCAL;
      double dy = HCALy - yexpect_HCAL;
      
      //Resolve the hadron spots without cuts
      hdx_HCAL->Fill( dx );
      hdy_HCAL->Fill( dy );
      hdxdy_HCAL->Fill( dy, dx );
      hxy_HCAL->Fill( HCALy, HCALx );
      
      // Coincidence timing cut and vertex cut to resolve W well
      if( fabs(diff-510)<40 && fabs(BBtr_vz[0])<0.08 ) hW_cut->Fill( W );
      
      // Preliminary HCal projections with single cut on W
      if( fabs(W-W_mean)<W_sig ){
	hdxdy_HCAL_wcut->Fill( dy, dx );
	hdx_HCAL_wcut->Fill( dx );
      }
      
      //Populate position histograms with cuts
      hdxdy_HCAL_cut->Fill( dy, dx );
      hdx_HCAL_cut->Fill( dx );
      
      //Populate BB/HCal correlation histograms from elastics
      hPAngleCorr_ph->Fill( ephi, phinucleon );
      hPAngleCorr_th->Fill( etheta, thetanucleon );
      
      //cout << ephi << " " << phinucleon << endl;
      //cout << etheta << " " << thetanucleon << endl;

      //Fill vertex position histogram for cut on tracks
      hvz_cut->Fill( BBtr_vz[0] );
      
      //Check "elastic" events on center HCal for id with spot checks
      bool HCAL_on = false, is_p = false, is_n = false;
      if( HCALy>-0.75 && HCALy<0.75 && HCALx>-2.015 && HCALx<1.285 ){
	HCAL_on = true;
      }
      if( pow( (HCALx-xexpect_HCAL - dx0_p)/dx_sig_p,2) + pow( (HCALy-yexpect_HCAL - dy0_p)/dy_sig_p,2) <= pow(2.5,2) ){
	is_p = true;
      }
      if( pow( (HCALx-xexpect_HCAL - dx0_n)/dx_sig_n,2) + pow( (HCALy-yexpect_HCAL - dy0_n)/dy_sig_n,2) <= pow(2.5,2) ){
	is_n = true;
      }
      
      //if( HCAL_on==true && (is_p==true || is_n==true) ) hdxdy_HCAL_cut->Fill( dy, dx );
      if( HCAL_on==true && is_n==true ) hdxdy_HCAL_ncut->Fill( dy, dx );
      if( HCAL_on==true && is_p==true ) hdxdy_HCAL_pcut->Fill( dy, dx );
      
      if( HCAL_on==true && is_n==true ){
	if( (HCALx-dx_max)>-2.015 ){
	  hdxdy_HCAL_fcut->Fill( dy, dx );
	  hdx_HCAL_fcut->Fill( dx );
	  hW_fcut->Fill( W );
	  elasYield++;
	}
      }else if( HCAL_on==true && is_p==true ){
	if( (HCALx+dx_max)<1.285 ){
	  hdxdy_HCAL_fcut->Fill( dy, dx );
	  hdx_HCAL_fcut->Fill( dx );
	  hW_fcut->Fill( W );
	  elasYield++;

	}
      }
    }
  }
  
  cout << endl;

  //Declare canvas
  TCanvas *c1 = new TCanvas("c1","c1",1600,1200);
  c1->Divide(2,1);

  //Write p-n distribution with spot cuts visualized
  c1->cd(1);
  hdxdy_HCAL_wcut->GetYaxis()->SetTitleOffset(1.3);
  hdxdy_HCAL_wcut->Draw("colz");

  TEllipse el_p;
  el_p.SetFillStyle(4005);
  el_p.SetLineColor(2);
  el_p.SetLineWidth(3);
  el_p.DrawEllipse(dy0_p,dx0_p,1.5*dy_sig_p,1.5*dx_sig_p,0,360,0);

  TEllipse el_n;
  el_n.SetFillStyle(4005);
  el_n.SetLineColor(4);
  el_n.SetLineWidth(3);
  el_n.DrawEllipse(dy0_n,dx0_n,1.5*dy_sig_n,1.5*dx_sig_n,0,360,0);

  //Write dx fit with visualization
  c1->cd(2);
  gStyle->SetOptFit(1111);

  //Fit dxdy histograms and extract xsig/ysig
  double p_low, p_hi, p_par[3];
  double n_low, n_hi, n_par[3];
  double bg_par[3];

  //defining fit parameters and range
  if( magSet==30 ){
    p_low = -1.13;
    p_hi = -0.2;
    n_low = -0.2;
    n_hi = 0.41;
    p_par[0] = 1967;
    p_par[1] = -0.58;
    p_par[2] = 0.16;
    n_par[0] = 662;
    n_par[1] = -0.047;
    n_par[2] = 0.184;
    bg_par[0] = 480;
    bg_par[1] = -0.3; 
    bg_par[2] = 0.5;
  }
  if( magSet==50 ){
    p_low=-1.4;
    p_hi=-0.78;
    n_low=-0.22;
    n_hi=0.26;
    p_par[0]=475;
    p_par[1]=-1.09;
    p_par[2]=0.184;
    n_par[0]=182;
    n_par[1]=0.009;
    n_par[2]=0.169;
    bg_par[0]=160;
    bg_par[1]=-0.42;
    bg_par[2]=0.7;
  }
  
  double par[9];
  TF1* fitBG;

  if( magSet==30 ){
    fitBG = new TF1("fitBG", BG_fit_30, -2.5, 2.5, 3 );
  }else if( magSet==50){
    fitBG = new TF1("fitBG", BG_fit_50, -2.5, 2.5, 3 );
  }else{
    cout << "ERROR: Poorly defined SBS field strength (not defined in code)." << endl;
    return 1;
  }
  fitBG->SetParameters(bg_par[0],bg_par[1],bg_par[2]);
  fitBG->SetLineColor(16);
  hdx_HCAL_fcut->Fit(fitBG,"QNR");
  fitBG->GetParameters(&par[0]);

  TF1* fitP = new TF1("fitP", "gaus", p_low, p_hi);
  fitP->SetParameters(p_par[0],p_par[1],p_par[2]);
  fitP->SetLineColor(2);
  hdx_HCAL_fcut->Fit(fitP,"QNR+");
  fitP->GetParameters(&par[3]);

  TF1* fitN = new TF1("fitN", "gaus", n_low, n_hi);
  fitN->SetParameters(n_par[0],n_par[1],n_par[2]);
  fitN->SetLineColor(4);
  hdx_HCAL_fcut->Fit(fitN,"QNR+");
  fitN->GetParameters(&par[6]);

  double parAll[9];
  TF1* fitAll = new TF1("fitAll","gaus(0)+gaus(3)+gaus(6)",-2.5,2.5);
  fitAll->SetLineColor(3); 
  fitAll->SetLineWidth(2);
  fitAll->SetParameters(par);
  hdx_HCAL_fcut->Fit(fitAll,"QR+");
  fitAll->GetParameters(parAll);

  //Extract GMn
  double UB_n_tot, UB_p_tot, n_tot, p_tot;
  TF1 *isolBG = new TF1("isolBG","gaus",-2.5,2.5);
  isolBG->SetLineWidth(1);
  isolBG->SetLineColor(15);
  isolBG->SetParameters(&par[0]);
  isolBG->Draw("same");
  TF1 *isolP = new TF1("isolP","gaus",-2.5,2.5);
  isolP->SetLineWidth(1);
  isolP->SetLineColor(kRed);
  isolP->SetParameters(&par[3]);
  isolP->Draw("same");
  TF1 *isolN = new TF1("isolP","gaus",-2.5,2.5);
  isolN->SetLineWidth(1);
  isolN->SetLineColor(kBlue);
  isolN->SetParameters(&par[6]);
  isolN->Draw("same");

  //Blind the yields
  UB_p_tot = isolP->Integral(-2.5,2.5)/hdx_HCAL_fcut->GetBinWidth(1);
  UB_n_tot = isolN->Integral(-2.5,2.5)/hdx_HCAL_fcut->GetBinWidth(1);
  p_tot = blind( UB_p_tot );
  n_tot = blind( UB_n_tot );

  //Add a legend
  auto legend = new TLegend(0.1,0.7,0.48,0.9);
  legend->SetTextSize(0.03);
  legend->AddEntry(isolBG,"Background Fit","l");
  legend->AddEntry(isolP,Form("Proton Fit, Yield: %f", p_tot),"l");
  legend->AddEntry(isolN,Form("Neutron Fit, Yield: %f", n_tot),"l");
  legend->AddEntry(fitAll,"Total Fit","l");
  legend->Draw();

  //Write to console
  cout << endl << endl;

  cout << "GMn ~ sqrt(n_tot/p_tot) = " << sqrt(n_tot/p_tot) << "." << endl;

  //Write out diagnostic histos, print to console, cleanup
  TString plotsfilename = outputfilename;
  plotsfilename.ReplaceAll(".root",".pdf");
  c1->Print(plotsfilename.Data(),"pdf");
  plotsfilename.ReplaceAll(".pdf",".png");
  c1->Print(plotsfilename.Data(),"png");
  fout->Write();
 
  cout << endl << endl << "Elastic yield for analyzed runs: " << elasYield << ". Total events analyzed: " << Nevents << "." << endl << endl;
  cout << "GMn analysis complete and histograms written to file." << endl;

  elist->Delete();
  st->Stop();

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

}

//Sample setup file, comment with #
////////////////////////////////////////////////////////
// #LD2, SBS4
// /lustre19/expphy/volatile/halla/sbs/pdbforce/gmn-replays/rootfiles/sbs4-sbs30p-ld2/e*
// endlist
// bb.tr.n==1&&bb.ps.e>0.2&&abs(bb.tr.vz[0])<.08&&bb.gem.track.nhits>3&&abs(bb.etot_over_p-0.92)<0.2&&sbs.hcal.e>0.025&&bb.ps.e+bb.sh.e>1.7
// endcut
// kine 4
// E_e 3.728
// BB_d 1.7988
// BB_th 36.
// HCal_d 11.
// HCal_th 31.9
// W_mean 0.73
// W_sig 0.1
// dx0 0.9
// dy0 0.62
// dx_sig 0.09
// dy_sig 0.15
// dx_max 0.65
/////////////////////////////////////////////////////////
