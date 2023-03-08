//SSeeds 6.23.22 - GMn Post-production - This script places tight elastic cuts on events and will extract TOF corrections to TDC values based on reconstructed hadron path lengths from the target to the face of HCal. Must take care to attend to all setup file parameters for proper configuration of the the output and one should not mix LH2 with LD2 or different magnet settings when stringing several files together.
//USAGE, configure setup file, then from command line with working root build: root -l 'TOF.C()')

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
const double HCal_Xi = -2.090; // Distance from beam center to top of HCal in analysis coordinates (m)
const double HCal_Xf = 1.360; // Distance from beam center to bottom of HCal in analysis coordinates (m)
const double HCal_Yi = -0.825; // Distance from beam center to opposite-beam side of HCal in analysis coordinates (m)
const double HCal_Yf = 0.825; // Distance from beam center to beam side of HCal in analysis coordinates (m)
const double HCal_div = 0.15; // Transverse width in x and y per cell
const Int_t xN = 24; //kNrows, total number of dispersive bins detection uni
const Int_t yN = 12; //kNcols, total number of transverse bins detection uni

//Constants
const double PI = TMath::Pi();
const double M_e = 0.00051; //Mass electron (GeV)
const double M_p = 0.938272; //Mass proton (GeV)
const double M_n = 0.939565; //Mass neutron (GeV)
const double v_gamma = 299792458; //Speed of light (m/s)


//Static Target Parameters
const double l_tgt = 0.15; // Length of the target (m)
const double rho_tgt = 0.0723; // Density of target (g/cc)
const double rho_Al = 2.7; // Density of aluminum windows (g/cc)
const double celldiameter = 1.6*2.54; //cm, right now this is a guess
const double Ztgt = 1.0; //Single proton for LH2, LD2
const double Atgt = 1.0; //Single proton for LH2, LD2
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

void TOF_alt( const char *configfilename="sTOF.cfg", const char *outputfilename="pelas_out.root" ){
  
  // Define a clock to check overall time
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
  int useAlshield = 1; // Use 1/8" aluminum shield on scattering chamber exit?
  string tar = "LH2"; // Target used during run set, LH2 liquid hydrogen, LD2 liquid deuterium
  vector<TString> log; // Log of runs analyzed

  // Declare variables for keeping track of proton and neutron yield. Declared outside of while loop in order to keep track of totals over all events.
  int Np = 0;
  int Nn = 0;

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
  // Applying simple global cuts for elastics from setup file
  TCut globalcut = "";
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endcut") ){
    if( !currentline.BeginsWith("#") ){
      globalcut += currentline;
    }    
  }

  // Loading all parameters in setup file
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

  // Populate new event list with globalcuts applied
  TEventList *elist = new TEventList("elist","Elastic Event List");
  C->Draw(">>elist",globalcut);

  cout << "Event list populated with loose cut placed on elastics." << endl;
  
  // Declare general tree parameters to read from TChain
  double BBtr_p[maxTracks], BBtr_px[maxTracks], BBtr_py[maxTracks], BBtr_pz[maxTracks];
  double BBtr_vx[maxTracks], BBtr_vy[maxTracks], BBtr_vz[maxTracks], BBtr_chi2[maxTracks];
  double BBfp_x[maxTracks], BBfp_y[maxTracks], BBfp_th[maxTracks], BBfp_ph[maxTracks];
  double BBtgt_x[maxTracks], BBtgt_y[maxTracks], BBtgt_th[maxTracks], BBtgt_ph[maxTracks];
  double BBtr_n, BBps_x, BBps_y, BBps_e, BBsh_x, BBsh_y, BBsh_e;

  double TDCT_id[maxTdcChan], TDCT_tdc[maxTdcChan]; 
  int TDCTndata;

  double HCALx, HCALy, HCALe;

  // Declare root tree variables and set values to memory locations in root file
  // Switch everything off (0) first to save memory
  C->SetBranchStatus( "*", 0 );
  // Switch only the branches we want on (1)
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

  // Set local variables with the memory addresses from TChain
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
  
  // Declare outfiles
  outputfilename = Form( "outfiles/TOF_SBS%d_tar%s_%s.root", magSet, tar.c_str(), date.c_str() );
  TFile *fout = new TFile( outputfilename, "RECREATE" );
  string logpath = Form( "outfiles/TOFLog_%s.log", date.c_str() );

  // Initialize misc. variables
  int elasYield = 0; // Keep track of total elastics analyzed
  double pBeam = E_e/(1.+E_e/M_p*(1.-cos(BB_th))); //Beam momentum

  //Mean energy loss of the beam prior to the scattering:
  double Eloss_outgoing = celldiameter/2.0/sin(BB_th) * rho_tgt * dEdx_tgt; //Energy lost in target. Approximately 1 MeV, could correct further with raster position
  if( useAlshield != 0 ) Eloss_outgoing += Alshieldthick * rho_Al * dEdx_Al; //Add energy loss from aluminum shield for runs where relevant

  TH1D *h_W = new TH1D("h_W",";W (GeV);",250,0,2); //Invarient mass 1D histo
  TH1D *h_Wcut = new TH1D("h_Wcut",";W (GeV);",250,0,2); //Invarient mass 1D histo

  //TH1D *Pelastic = new TH1D("Pelastic","Pelastic", 250,0,5);
  //TH1D *Precon = new TH1D("Precon","Precon",250,0,5);
  //TH1D *dx_dz = new TH1D("dx_dz","dx_dz",250,-5,5);
  //TH1D *dy_dz = new TH1D("dy_dz","dy_dz",250,-5,5);

  //TH1D *h_Q2 = new TH1D("h_Q2",";Q2 (GeV^2)",500,0,5); //Inverse momentum transfer before W2 and timing cuts are applied
  //TH1D *h_Q2cut = new TH1D("h_Q2cut",";Q2 (GeV^2)",500,0,5); //Inverse momentum transfer after cuts to compare

  TH2D *timep_vs_x = new TH2D("t_p_vs_x","t_px",xN,HCal_Xi,HCal_Xf,250,0,20);
  TH2D *timen_vs_x = new TH2D("t_n_vs_x","t_nx",xN,HCal_Xi,HCal_Xf,250,0,20); 
  TH2D *timep_vs_y = new TH2D("t_p_vs_y","t_py",yN,HCal_Yi,HCal_Yf,250,0,20);
  TH2D *timen_vs_y = new TH2D("t_n_vs_y","t_ny",yN,HCal_Yi,HCal_Yf,250,0,20);
  TH2D *timep_vs_row = new TH2D("t_p_vs_row","t_prow",500,-3,3,250,0,20);
  TH2D *timen_vs_row = new TH2D("t_n_vs_row","t_nrow",500,-3,3,250,0,20); 
  TH2D *timep_vs_col = new TH2D("t_p_vs_col","t_pcol",500,-3,3,250,0,20);
  TH2D *timen_vs_col = new TH2D("t_n_vs_col","t_ncol",500,-3,3,250,0,20);
  TH2D *ny_vs_x = new TH2D("ny_vs_x","y_x",500,-1,1,250,-3,1);
  TH2D *py_vs_x = new TH2D("py_vs_x","y_x",500,-1,1,250,-3,1);

  //Add hadron spot histogram to evaluate cuts on LD2
  TH2D *hdxdy_HCAL = new TH2D("hdxdy_HCAL",";y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)", 250, -2.0, 2.0, 500, -4, 4 );
  TH1D *hdx_HCAL = new TH1D("hdx_HCAL",";x_{HCAL}-x_{expect} (m)",500,-4,4);
  TH1D *hdy_HCAL = new TH1D("hdy_HCAL",";y_{HCAL}-y_{expect} (m)",250,-2,2);
  TH1D *h_TOFn = new TH1D("h_TOFn","Time of Flight: Neutron; ns",600,35,45);
  TH1D *h_TOFp = new TH1D("h_TOFp","Time of Flight: Proton; ns",600,35,45);
  TH1D *h_nD = new TH1D("h_nD","Neutron Path Length",1000,10.5,11.5);
  TH1D *h_pD = new TH1D("h_nP","Proton Path Length (Naive)",1000,10.5,11.5);

  TH1D *h_hcalX = new TH1D("h_hcalX","HCal Hits X",xN,HCal_Xi,HCal_Xf);
  TH1D *h_hcalY = new TH1D("h_hcalY","HCal Hits Y",yN,HCal_Yi,HCal_Yf);
  
  TH1D *timep = new TH1D("t_p","t_p",300,0,6);
  TH1D *timen = new TH1D("t_n","t_n",300,0,6);

  TH2D *t_nvID = new TH2D("tn_vs_ID","tn_vs_ID",288,0,288,200,37,40);
  TH2D *t_pvID = new TH2D("tp_vs_ID","tp_vs_ID",288,0,288,200,37,40);

  // Set long int to keep track of total entries after globalcut
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
      
      //Create a progress bar in terminal
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
      //Elements in dynamic arrays must be referenced from id
      for(int ihit=0; ihit<TDCTndata; ihit++){ 
	if(TDCT_id[ihit]==5) bbcal_time=TDCT_tdc[ihit];
	if(TDCT_id[ihit]==0) hcal_time=TDCT_tdc[ihit];
      }
      double diff = hcal_time - bbcal_time; 

      //Cut here on BB trigger/HCal trigger coincidence time. 510ns from cable lengths good for GMn generally, but probably shouldn't be hardcoded.
      if( fabs(diff-510.)>10. ) continue;

      //Correct the beam energy with energy loss in target using vertex position
      double Eloss = (BBtr_vz[0]+l_tgt/2.0) * rho_tgt * dEdx_tgt + uwallthick_LH2 * rho_Al * dEdx_Al; //approximately 3 MeV
      double E_corr = E_e - Eloss;

      double p_corr = BBtr_p[0] - Eloss_outgoing; //Neglecting the mass of e'

      //Proceed only if at least one track exists in BB arm - lowest chi2 (most probable) track always first element
      double etheta = acos( BBtr_pz[0]/BBtr_p[0] ); //Will use the uncorrected track momentum to reconstruct e' theta
      double ephi = atan2( BBtr_py[0], BBtr_px[0] );
      
      TVector3 vertex( 0, 0, BBtr_vz[0] ); // z location of vertex in hall coordinates
      TLorentzVector Pbeam( 0, 0, E_corr, E_corr ); //Mass of e negligable
      TLorentzVector kprime( BBtr_px[0], BBtr_py[0], BBtr_pz[0], BBtr_p[0] ); //e' track momentum four-vector in lab frame
      TLorentzVector Ptarg( 0, 0, 0, M_p ); //Target proton momentum four-vector in lab frame
      
      TLorentzVector q = Pbeam - kprime; // Momentum transfer. q^2 = -Q^2
      TLorentzVector PgammaN = Ptarg + q; //(-px, -py, ebeam - pz, Mp + ebeam - p)
      
      double nu = E_corr - BBtr_p[0]; //momentum transferred to proton
      double pp = sqrt( pow(nu,2)+2.*M_p*nu ); //proton momentum
      double pn = sqrt( pow(nu,2)+2.*M_n*nu ); //neutron momentum
      double phinucleon = ephi + PI; //nucleon phi angle from exit beamline assuming coplanarity between e' and p
      double thetanucleon = acos( (E_corr - BBtr_pz[0])/pp ); //nucleon theta angle from exit beamline using elastic constraint on nucleon kinematics
      
      TVector3 pNhat( sin(thetanucleon)*cos(phinucleon), sin(thetanucleon)*sin(phinucleon), cos(thetanucleon) ); //Scattered nucleon four-vector
      
      //Define HCal local coordinate system
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
      
      double p_ep = BBtr_p[0]; // Obtain scattered electron momentum
      
      double Q2 = 2*E_corr*E_ep*( 1-(BBtr_pz[0]/p_ep) ); // Obtain Q2 from beam energy, outgoing electron energy, and momentum
      //h_Q2->Fill(Q2); //Fill histogram before cut
      //HCalx->Fill(HCALx);
      
      //Get invarient mass transfer W from the four-momentum of the scattered nucleon
      double W = PgammaN.M();
      //Fill histogram to check elastic cut
      h_W->Fill( W );

      //Use the electron kinematics to predict the proton momentum assuming elastic scattering on free proton at rest (will need to correct for fermi motion):
      double E_pp = nu+M_p; // Get energy of the proton
      double Enucleon = sqrt(pow(pp,2)+pow(M_p,2)); // Check on E_pp, same
      double KE_p = nu; // Kinetic energy of scattered nucleon, for elastics

      double dx = HCALx - xexpect_HCAL; // (Energy weighted center of cluster (x component)) - (straight line nucleon projection obtained from e' track HCal location (x component))
      double dy = HCALy - yexpect_HCAL; // (Energy weighted center of cluster (y component)) - (straight line nucleon projection obtained from e' track HCal location (y component))
            
      //Calculate the hadron spot/s 
      hdxdy_HCAL->Fill( dy, dx );
      hdx_HCAL->Fill( dx );
      hdy_HCAL->Fill( dy );

      //Main cut on elastic invarient mass
      if( fabs(W-W_mean)>W_sig ) continue; 
      h_Wcut->Fill( W );

      //Check "elastic" events on center HCal for id with spot checks
      bool HCAL_on = false, is_p = false, is_n = false; //is_p = proton, is_n = neutron
      if( HCALy>-0.75 && HCALy<0.75 && HCALx>-2.015 && HCALx<1.285 ){ //Bounds set by acceptance of HCal in HCal coordinates, should probably not be hardcoded
	HCAL_on = true; //true when expected nucleon on hcal
      }
      if( pow( (dx - dx0_p)/dx_sig_p,2) + pow( (dy - dy0_p)/dy_sig_p,2) <= pow(2.5,2) ){
	is_p = true;
      }
      if( pow( (dx - dx0_n)/dx_sig_n,2) + pow( (dy - dy0_n)/dy_sig_n,2) <= pow(2.5,2) ){
	is_n = true;
      }
      
      if( HCAL_on==true && is_n==true ){
	if( (HCALx-dx_max)>-2.015 ){ //Fiducial check to verify that a scattered nucleon reconstructed from Bigbite would have landed on HCal
	  if( tar == "LD2" ){ //Consider neutrons only if we have a dueterium target
	   
	    double pelastic = E_corr/(1.+(E_corr/M_n)*(1.0-cos(etheta))); //Corrected momentum of elastic nucleon
      
	    //double precon = BBtr_p[0] + Eloss_outgoing; //reconstructed momentum, corrected for mean energy loss exiting the target (later we'll also want to include Al shielding on scattering chamber)

	    //double nu_recon = E_corr - precon;
	    //double Q2recon = 2.0*E_corr*precon*(1.0-cos(etheta));
	    //double W2recon = pow(M_p,2) + 2.0*M_p*nu_recon - Q2recon;
	    double nD = sqrt( pow(HCal_d,2) + pow(HCALx,2) + pow(HCALy,2) );
	    h_nD->Fill( nD );
	    double nGamma = pn/M_n;
	    double nBeta = sqrt(1-pow(1/nGamma,2));
	    double TOF_n = (nD/(nBeta*v_gamma))*pow(10,9);

	    h_TOFn->Fill( TOF_n );

	    double dxdz2 = pow((HCALx/HCal_d),2);
	    double dydz2 = pow((HCALy/HCal_d),2);
	    double part1 = (HCal_d)/(pow((1+(dxdz2)+(dydz2)),1/2));
	    double part2 = (pow(pow(pelastic,2)+pow(M_n,2),1/2))/pelastic;
	    double t_n = part1*part2;

	    //t_nvID->Fill(channelID,TOF_n);

	    timen->Fill(TOF_n);
	    timen_vs_x->Fill(HCALx,TOF_n);
	    timen_vs_y->Fill(HCALy,TOF_n);
	    ny_vs_x->Fill(HCALy,HCALx);
	    //cout << "Neutron time-of-flight =" << t_n << "." << endl;
	    //h_W2->Fill(W2recon);
	    //h_W->Fill(W);
	    //Pelastic->Fill(pelastic);
	    //Precon->Fill(precon);
	    //dx_dz->Fill(dxdz2);
	    //dy_dz->Fill(dydz2);
	    //h_Q2cut->Fill(Q2recon);

	    Nn++;

	    elasYield++;
	  }

	}
      }else if( HCAL_on==true && is_p==true ){
	if( (HCALx+dx_max)<1.285 ){ //Fiducial check to verify that a scattered nucleon reconstructed from Bigbite would have landed on HCal
	 
	  double pelastic = E_corr/(1.+(E_corr/M_p)*(1.0-cos(etheta))); //Corrected momentum of elastic nucleon
      
	  //double precon = BBtr_p[0] + Eloss_outgoing; //reconstructed momentum, corrected for mean energy loss exiting the target (later we'll also want to include Al shielding on scattering chamber)

	  //double nu_recon = E_corr - precon;
	  //double Q2recon = 2.0*E_corr*precon*(1.0-cos(etheta));
	  //double W2recon = pow(M_p,2) + 2.0*M_p*nu_recon - Q2recon;
	  double pD = sqrt( pow(HCal_d,2) + pow(HCALx,2) + pow(HCALy,2) );
	  h_pD->Fill( pD );
	  double pGamma = pp/M_p;
	  double pBeta = sqrt(1-pow(1/pGamma,2));
	  double TOF_p = (pD/(pBeta*v_gamma))*pow(10,9);

	  h_TOFp->Fill( TOF_p );

	  double dxdz2 = pow((HCALx/HCal_d),2);
	  double dydz2 = pow((HCALy/HCal_d),2);
	  double part1 = (HCal_d)/(pow((1+(dxdz2)+(dydz2)),1/2));
	  double part2 = (pow(pow(pelastic,2)+pow(M_p,2),1/2))/pelastic;
	  double t_p = part1*part2;
	  //cout << "Proton time-of-flight =" << t_p << "." << endl;

	  //t_pvID->Fill(channelID,TOF_p);

	  timep->Fill(TOF_p);
	  timep_vs_x->Fill(HCALx,TOF_p);
	  timep_vs_y->Fill(HCALy,TOF_p);
	  py_vs_x->Fill(HCALy,HCALx);
	  //h_W2->Fill(W2recon);
	  //h_W->Fill(W);
	  //Pelastic->Fill(pelastic);
	  //Precon ->Fill(precon);
	  //dx_dz->Fill(dxdz2);
	  //dy_dz->Fill(dydz2);


	  //h_Q2cut->Fill(Q2);
	 
	  Np++;

	  elasYield++;
	}
      }

      ////!!!!One could add additional code here populate histograms with TOF information from events that pass cuts
      ////
      ////
      ////asdfw
    }
  }
  
  cout << endl;
  
  /*
  //Construct graphs for uniformity check
  Double_t posErr[xN] = {0.};
  TF1 *f1;
  //For X (dispersive)
  Double_t X[xN];
  Double_t Xval[xN];
  Double_t Xerr[xN];
  TH1D *Xslice[xN];
  for( int x=0; x<xN; x++ ){
    X[x] = HCal_Xi + x*HCal_div;
    Xslice[x] = timep_vs_x->ProjectionY(Form("Xslice_%d",x+1),x+1,x+1);
    Xslice[x]->Fit("gaus","Q","",0.01,20);
    f1=Xslice[x]->GetFunction("gaus");
    if(Xslice[x]->GetEntries()>1){
      Xval[x] = f1->GetParameter(1);
      Xerr[x] = f1->GetParameter(2);
    }else{
      Xval[x] = -1;
      Xerr[x] = 0;
    }
  }
  TGraphErrors *cTOF_X = new TGraphErrors( xN, X, Xval, posErr, Xerr );
  cTOF_X->GetXaxis()->SetLimits(HCal_Xi-0.05,HCal_Xf+0.05);  
  cTOF_X->GetYaxis()->SetLimits(0.0,1.0);
  cTOF_X->SetTitle("Time of Flight Proton - Dispersive X");
  cTOF_X->GetXaxis()->SetTitle("X (m)");
  cTOF_X->GetYaxis()->SetTitle("E_{HCAL}/KE_{exp}");
  cTOF_X->SetMarkerStyle(20); // idx 20 Circles, idx 21 Boxes
  cTOF_X->Write("cTOF_X");

  //For Y (transverse)
  Double_t Y[yN];
  Double_t Yval[yN];
  Double_t Yerr[yN];
  TH1D *Yslice[yN];
  for( int x=0; x<yN; x++ ){
    Y[x] = HCal_Yi + x*HCal_div;
    Yslice[x] = timep_vs_y->ProjectionY(Form("Yslice_%d",x+1),x+1,x+1);
    Yslice[x]->Fit("gaus","Q","",0.01,20.);
    f1=Yslice[x]->GetFunction("gaus");
    if(Yslice[x]->GetEntries()>1){
      Yval[x] = f1->GetParameter(1);
      Yerr[x] = f1->GetParameter(2);
    }else{
      Yval[x] = -1;
      Yerr[x] = 0;
    }
  }
  TGraphErrors *cTOF_Y = new TGraphErrors( yN, Y, Yval, posErr, Yerr );
  cTOF_Y->GetXaxis()->SetLimits(HCal_Yi-0.05,HCal_Yf+0.05);  
  cTOF_Y->GetYaxis()->SetLimits(0.0,1.0);
  cTOF_Y->SetTitle("Time of Flight Proton - Transverse Y");
  cTOF_Y->GetXaxis()->SetTitle("Y (m)");
  cTOF_Y->GetYaxis()->SetTitle("E_{HCAL}/KE_{exp}");
  cTOF_Y->SetMarkerStyle(20); // idx 20 Circles, idx 21 Boxes
  cTOF_Y->Write("cTOF_Y");
  */
  ofstream logfile;
  logfile.open( logpath );
  logfile << "#Log of runs corresponding to root file with same date" << endl;
  for( int i=0; i<log.size(); i++ ){
    logfile << log[i] << endl;
  }

  logfile.close();

  cout << endl << endl << "Elastic yield for analyzed runs: " << elasYield << ". Total events analyzed: " << Nevents << "." << endl << endl;

  cout << endl << "Elastic proton yield: " << Np << "." << endl;

  cout << endl << "Elastic neutron yield: " << Nn << "." << endl;

  cout << endl << endl << "New root file generated with elastic cuts placed at " << outputfilename << endl;

  fout->Write();

  elist->Delete();
  st->Stop();

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

}
