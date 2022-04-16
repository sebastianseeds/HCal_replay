//SSeeds 3.4.22 - Post-production - Calibration code which employs best current cuts on elastic events to obtain ADC gain calibration parameters (pC/GeV).

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

const int kNcell = 288; // Total number of HCal modules
const int kNrows = 24; // Total number of HCal rows
const int kNcols = 12; // Total number of HCal columns
const int maxTracks = 1000; // Reasonable limit on tracks to be stored per event
const int maxTdcChan = 10; // Set to accomodate original 5 TDCTrig channels with buffer

double hcalheight = 0.365; //m The height of the center of HCAL above beam
double sampFrac = 0.0795; //HCal sampling frac (0.06588 GeV/0.8286 GeV) = 0.0795 = 7.95% -> (MC E_dep per proton) / (fit to data KE_p)

const double PI = TMath::Pi();
const double M_e = 0.00051;
const double M_p = 0.938272;

double oldGain[kNcell] = {0.0};
double oldRatio[kNcell] = {0.0};

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

void HCalECal( const char *configfilename = "setup_HCalECal.cfg", int run = -1 ){
  
  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  // Get the date
  string date = getDate();
  
  // Declare Chain for many root files
  TChain *C = new TChain("T");

  // Paths for input/output files
  string inConstPath = "/w/halla-scshelf2102/sbs/seeds/HCal_replay/hcal/setFiles/const_v2.txt"; //Settings from beam data for all runs prior to HV change
  string constPath = "/w/halla-scshelf2102/sbs/seeds/HCal_replay/hcal/setFiles/newRatio.txt";
  string ratioPath = "/w/halla-scshelf2102/sbs/seeds/HCal_replay/hcal/setFiles/ratio.txt";
  string oneBlockPath = "/w/halla-scshelf2102/sbs/seeds/HCal_replay/hcal/setFiles/oneBlock.txt";

  // Declare general physics parameters to be modified by input config file
  double kine = 8; // Keep track of kinematic calibrating from
  double E_e = 1.92; // Energy of beam (incoming electrons from accelerator)
  int minEventPerCell = 500; // Minimum number of scattered p in cell required to calibrate
  int maxEventPerCell = 4000; // Maximum number of scattered p events to contribute
  double highDelta = 0.1; // Minimum M(i,j)/b(i) factor allowed 
  double HCal_d = 14.5; // Distance to HCal from scattering chamber for comm1
  double HCal_th = 35.0; // Angle that the center of HCal is at  
  double W_mean = 0.93; // Mean of W at current kinematic
  double W_sig = 0.039; // Width of W at current kinematic
  double dx0 = 0.9; // Position of proton spot, x-x_expected
  double dy0 = 0.62; // Position of proton spot, y-y_expected
  double dx_sig = 0.09; // Max spread of proton spot, x-x_expected
  double dy_sig = 0.15; // Max spread of proton spot, y-y_expected

  int elasYield = 0; // Keep track of total elastics analyzed

  // Reading config file
  ifstream configfile(configfilename);
  TString currentline;
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
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
	kine = sval.Atof();
	cout << "Loading kinematic setting: " << kine << endl;
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

  cout << "Setup parameters loaded." << endl;

  TEventList *elist = new TEventList("elist","Elastic Event List");
  C->Draw(">>elist",globalcut);

  cout << "Event list populated with cut placed on elastics." << endl;

  // Declare general detector parameters
  double BBtr_p[maxTracks], BBtr_px[maxTracks], BBtr_py[maxTracks], BBtr_pz[maxTracks];
  double BBtr_vz[maxTracks], BBtr_chi2[maxTracks];
  double BBtr_n, BBps_x, BBps_y, BBps_e, BBsh_x, BBsh_y, BBsh_e;

  double TDCT_id[maxTdcChan], TDCT_tdc[maxTdcChan]; 
  int TDCTndata;

  double HCALx, HCALy, HCALe;
  double crow, ccol, nblk;
  double cblkid[kNcell], cblke[kNcell];


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

  C->SetBranchAddress( "sbs.hcal.x", &HCALx );
  C->SetBranchAddress( "sbs.hcal.y", &HCALy );
  C->SetBranchAddress( "sbs.hcal.e", &HCALe );
  C->SetBranchAddress( "sbs.hcal.rowblk", &crow );
  C->SetBranchAddress( "sbs.hcal.colblk", &ccol );
  C->SetBranchAddress( "sbs.hcal.nblk", &nblk ); // Total number of blocks in highest E clus
  C->SetBranchAddress( "sbs.hcal.clus_blk.id", cblkid ); // kNcell-1 index for each block
  C->SetBranchAddress( "sbs.hcal.clus_blk.e", cblke ); // Array of block energies
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

  // Create stopwatch to track processing time
  TStopwatch *sw = new TStopwatch();
  sw->Start();
  
  // Declare outfile
  TFile *fout = new TFile( "eCalEOut.root", "RECREATE" );
  
  // Initialize vectors and arrays
  double gOldConst[kNcell];
  double gRatio[kNcell] = {0.0}; //Old ratios from first iteration
  double oneBlock[kNcell] = {0.0};
  double GCoeff[kNcell] = {0.0};
  double GCoeff_oneblock[kNcell] = {0.0};
  double GCoeff_divide[kNcell] = {0.0};
  double A[kNcell] = {0.0}; // Array to keep track of ADC values per cell. Memory is reset on each event.

  // Initialize histograms	
  TH1D *hDeltaE = new TH1D( "hDeltaE","1.0-Eclus/p_rec", 100, -1.5, 1.5 );
  TH1D *hHCALe = new TH1D( "hHCALe","E HCal Cluster E", 400, 0., 4 );
  TH1D *hDiff = new TH1D( "hDiff","HCal time - BBCal time (ns)", 1300, -500, 800 );
  TH1D *hClusE = new TH1D( "hClusE","Best Cluster Energy", 100, 0.0, 2.0);
  TH2D *hPAngleCorr = new TH2D( "hPAngCorr","Track p vs Track ang", 100, 30, 60, 100, 0.4, 1.2 );
  TH2D *hPAngleCorr_2 = new TH2D( "hPAngCorr_2","Track p vs Track ang v2", 100, 30, 60, 100, 0.4, 1.2 );
  TH1D *hW = new TH1D( "W", "W", 250, 0.3, 1.5 );
  hW->GetXaxis()->SetTitle( "GeV" );
  TH1D *hNBlk = new TH1D( "hNBlk", "Number of Blocks in Primary Cluster", 25, 0, 25 );
  hW->GetXaxis()->SetTitle( "Number" );
  TH1D *hW_cuts = new TH1D( "W_cuts", "W_cuts", 250, 0.3, 1.5 );
  hW->GetXaxis()->SetTitle( "GeV" );
  TH1D *hQ2 = new TH1D( "Q2", "Q2", 250, 0.5, 3.0 );
  hQ2->GetXaxis()->SetTitle( "GeV" );
  TH1D *hE_ep = new TH1D( "Scattered Electron Energy","E_ep", 500, 0.0, E_e*1.5 ); 
  hE_ep->GetXaxis()->SetTitle( "GeV" );
  TH1D *hE_pp = new TH1D( "Scattered Proton Energy", "E_pp", 500, 0.0, E_e*1.5 );
  hE_pp->GetXaxis()->SetTitle( "GeV" );
  TH1D *hKE_p = new TH1D( "Scattered Proton Kinetic Energy", "KE_pp", 500, 0.0, E_e*1.5 );
  hKE_p->GetXaxis()->SetTitle( "GeV" );
  TH2D *hADC = new TH2D( "hADC", "HCal Int_ADC Spectra: W Cut", 288, 0, 288, 100., 0., 1. );
  TH2D *hADC_amp = new TH2D( "hADC_amp", "HCal ADC_amp Spectra: W Cut", 288, 0, 288, 100., 0., 10. );
  TH2D *hdxdy_HCAL = new TH2D("hdxdy_HCAL",";y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)", 250, -5.0, 5.0, 250, -10, 10 );
  TH2D *hXY_HCAL = new TH2D("hXY_HCAL",";y_{HCAL} (m); x_{HCAL} (m)", 250, -5.0, 5.0, 250, -10, 10 );
  TH1D *hvz = new TH1D("hvz",";vertex z (m);", 250,-0.125,0.125);
  TH1D *hvz_cut = new TH1D("hvz_cut",";vertex z (m);", 250,-0.125,0.125);
  TH2D *coefficients = new TH2D("coefficients",";Channel ;GeV/mV", 288,0,288,250,0,0.025);
  TH2D *coefficients_1b = new TH2D("coefficients_1b",";Channel ;GeV/mV", 288,0,288,250,0,0.025);
  TH2D *hEDiff_vs_X = new TH2D("hEDiff_vs_X",";X-pos (cm);(E_exp-E_dep)/E_exp (GeV)",360,-2,2,100,-3,3);  
  TH2D *hEDiff_vs_Y = new TH2D("hEDiff_vs_Y",";Y-pos (cm);(E_exp-E_dep)/E_exp (GeV)",180,-1,1,100,-3,3);
  TH2D *htDiff_vs_ADCint[kNcell+1];
  TH1D *htDiff[kNcell+1];
  TH1D *haDiff[kNcell+1];
  TH1D *hE_pp_exp_cell[kNcell+1];
  TH1D *hE_pp_cell[kNcell+1];
  TH1D *hE_pp_corr = new TH1D( "Deposited Elastic Proton Energy", "E_pp_corr", 500, 0.0, E_e*1.5 );

  for( int i=0; i<kNcell+1; i++ ){
    hE_pp_exp_cell[i] = new TH1D(Form("hE_pp_exp_cell_bl%d",i),Form(";E_{bl%d} (GeV)",i), 500, 0.0, E_e*1.5);
    hE_pp_cell[i] = new TH1D(Form("hE_pp_cell_bl%d",i),Form(";E_{bl%d} (GeV)",i), 500, 0.0, E_e*1.5);
  }

  // Set long int to keep track of total entries
  Long64_t Nevents = elist->GetN();

  // Read in previous ADC gain constants
  ifstream inConstFile( inConstPath );
  if( !inConstFile ){
    cerr << "No input constant file present -> setFiles/const_v2.txt expected." << endl;
    return 0;
  }

  cout << "Loading previous calibration constants.." << endl;
  int n1=0;
  double d1=0;
  string line;
  
  while( getline( inConstFile, line ) ){
    if( line.at( 0 )=='#' ) {
      continue;
    }
    istringstream iss( line );
    while( iss >> d1 ){
      gOldConst[n1] = d1;
      cout << "Previous ADC gain coeff for block " << n1 << ": " << d1 << endl;
      n1++;
    }
  }

  cout << endl << "All parameters loaded." << endl << endl;
  cout << "Opened up tree with nentries: " << Nevents << ".." << endl << endl;

  //Declare matrices for chi-square min calibration scheme and keep track of calibrated events
  TMatrixD Ma(kNcell,kNcell);
  TVectorD ba(kNcell);
  TVectorD bb(kNcell);
  int NEV[kNcell] = {0};
  int NEV_oneblock[kNcell] = {0};
  
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
      
      memset(A, 0, kNcell*sizeof(double));
	
      double etheta = acos( BBtr_pz[0]/BBtr_p[0] );
      double ephi = atan2( BBtr_py[0], BBtr_px[0] );

      TVector3 vertex(0,0,BBtr_vz[0]); // z location of vertex in hall coordinates
      TLorentzVector Pbeam(0,0,E_e,E_e); //Mass of e negligable
      TLorentzVector kprime(BBtr_px[0],BBtr_py[0],BBtr_pz[0],BBtr_p[0]);
      TLorentzVector Ptarg(0,0,0,M_p);

      TLorentzVector q = Pbeam - kprime;
      TLorentzVector PgammaN = Ptarg + q; //(-px, -py, ebeam - pz, Mp + ebeam - p)

      double pel = E_e/(1.+E_e/M_p*(1.-cos(etheta)));
      double nu = E_e - BBtr_p[0];
      double pp = sqrt(pow(nu,2)+2.*M_p*nu);
      double phinucleon = ephi + TMath::Pi(); //assume coplanarity
      double thetanucleon = acos( (E_e - BBtr_pz[0])/pp ); //use elastic constraint on nucleon kinematics
	
      TVector3 pNhat( sin(thetanucleon)*cos(phinucleon),sin(thetanucleon)*sin(phinucleon),cos(thetanucleon));

      //Define HCal coordinate system
      TVector3 HCAL_zaxis(sin(-HCal_th),0,cos(-HCal_th));
      TVector3 HCAL_xaxis(0,-1,0);
      TVector3 HCAL_yaxis = HCAL_zaxis.Cross(HCAL_xaxis).Unit();
	
      TVector3 HCAL_origin = HCal_d * HCAL_zaxis + hcalheight * HCAL_xaxis;

      //Plot 2D histo of position
      hXY_HCAL->Fill( HCALy, HCALx );

      //Define intersection points for hadron vector
      double sintersect = ( HCAL_origin - vertex ).Dot( HCAL_zaxis ) / (pNhat.Dot( HCAL_zaxis ) );

      TVector3 HCAL_intersect = vertex + sintersect * pNhat;

      double yexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_yaxis );
      double xexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_xaxis );

      //Calculate the proton spot - use for cut later on
      hdxdy_HCAL->Fill( HCALy - yexpect_HCAL, HCALx - xexpect_HCAL );

      //cout << HCALy - yexpect_HCAL << " " << HCALx - xexpect_HCAL << endl;

      double E_ep = sqrt( pow(M_e,2) + pow(BBtr_p[0],2) ); // Obtain the scattered electron energy
      hE_ep->Fill( E_ep ); // Fill histogram
	
      //double p_ep = opticsCorr*BBtr_p[0]; // Obtain the magnitude of scattered electron momentum with correction factor
      double p_ep = BBtr_p[0];
      double Q2 = 2*E_e*E_ep*( 1-(BBtr_pz[0]/p_ep) ); // Obtain Q2 from beam energy, outgoing electron energy, and momenta
      hQ2->Fill( Q2 ); // Fill histogram
	
      double W = PgammaN.M();
      hW->Fill( W );
	
      //Use the electron kinematics to predict the proton momentum assuming elastic scattering on free proton at rest:
      double E_pp = nu+M_p; // Get energy of the proton
      double Enucleon = sqrt(pow(pp,2)+pow(M_p,2)); // Check on E_pp, same
      hE_pp->Fill( E_pp ); // Fill histogram

      double KE_p = nu; //For elastics
      hKE_p->Fill( KE_p );   
	
      //Cut on BBCal and HCal trigger coincidence
      double bbcal_time=0., hcal_time=0.;
      for(int ihit=0; ihit<TDCTndata; ihit++){
	if(TDCT_id[ihit]==5) bbcal_time=TDCT_tdc[ihit];
	if(TDCT_id[ihit]==0) hcal_time=TDCT_tdc[ihit];
      }
      double diff = hcal_time - bbcal_time; 
      hDiff->Fill( diff ); // Fill histogram
	
      //Check that all branch addresses are set properly
      //cout << "BBtr_p[0]: " << BBtr_p[0] << endl;
      //cout << "BBtr_px[0]: " << BBtr_px[0] << endl;
      //cout << "BBtr_py[0]: " << BBtr_py[0] << endl;
      //cout << "BBtr_pz[0]: " << BBtr_pz[0] << endl;
      //cout << "BBtr_vz[0]: " << BBtr_vz[0] << endl;
      //cout << "BBtr_chi2[0]: " << BBtr_chi2[0] << endl;
      //cout << "BBtr_n: " << BBtr_n << endl;
      //cout << "BBps_x: " << BBps_x << endl;
      //cout << "BBps_y: " << BBps_y << endl;
      //cout << "BBps_e: " << BBps_e << endl;
      //cout << "BBsh_x: " << BBsh_x << endl;
      //cout << "BBsh_y: " << BBsh_y << endl;
      //cout << "BBsh_e: " << BBsh_e << endl;
      //cout << "TDCT_id[0]: " << TDCT_id[0] << endl;
      //cout << "TDCT_tdc[0]: " << TDCT_tdc[0] << endl;
      //cout << "TDCTndata: " << TDCTndata << endl;
      //cout << "HCALx: " << HCALx << endl;
      //cout << "HCALy: " << HCALy << endl;
      //cout << "HCALe: " << HCALe << endl;
      //cout << "crow: " << crow << endl;
      //cout << "ccol: " << ccol << endl;
      //cout << "nblk: " << nblk << endl;
      //cout << "cblkid[0]: " << cblkid[0] << endl;
      //cout << "cblke[0]: " << cblke[0] << endl;
      //cout << "HCALa[0]: " << HCALa[0] << endl;
      //cout << "HCALtdc[0]: " << HCALtdc[0] << endl;
      //cout << "HCALx: " << HCALx << endl;
      //cout << "xexpect_HCAL: " << xexpect_HCAL << endl;
      //cout << "HCALy: " << HCALy << endl;
      //cout << "yexpect_HCAL: " << yexpect_HCAL << endl;

      //Coincidence timing cut and vertex cut to resolve W well
      if( fabs(diff-510)<40 && fabs(BBtr_vz[0])<0.06 ) hW_cuts->Fill( W );

      //Additional cuts to select elastics
      if( fabs(W-W_mean)<W_sig && //Observed mean W cut
	  fabs(diff-510)<40 //Observed coincidence trigger HCal/BB
	  ){

	//Check energy deposited in HCal after elastics cuts
	hHCALe->Fill( HCALe );
	
	//Fill vertex position histogram for cut on tracks
	hvz_cut->Fill( BBtr_vz[0] );
	
	//Reject events where the primary block in the primary cluster is on the edge of the acceptance
	if( crow==0 ||
	    crow==23 ||
	    ccol==0 ||
	    ccol==11 )
	  continue;
	
	//Events that pass the above cuts constitute elastics
	elasYield++;
	
	double clusE = 0.0;
	
	hNBlk->Fill( nblk );
	
	// Get energies with simplest scheme from clusters only
	for( int blk = 0; blk<(int)nblk; blk++ ){
	  int blkid = int(cblkid[blk])-1; //-1 necessary since sbs.hcal.clus_blk.id ranges from 1 to 288 (different than just about everything else)
	  
	  clusE += cblke[blk];
	  A[blkid] += cblke[blk];
	  
	  hADC->Fill( blkid, cblke[blk] );
	  
	  // Simple estimation of the coefficients assuming 100% energy deposition in one block. Will check against the more sophisticated version with chi^2 reduction.
	  if(nblk==1)
	    NEV_oneblock[blkid]++;
	  
	  NEV[blkid]++;
	}
	
	//Fill some histograms
	double E_exp = E_pp*sampFrac;
	hDeltaE->Fill( 1.0-(clusE/KE_p) );
	hClusE->Fill( clusE );
	hE_pp_corr->Fill( E_exp );
	hPAngleCorr->Fill( ephi, p_ep );
	hPAngleCorr_2->Fill( etheta, p_ep );
	
	//Main Quality Check histograms
	double E_dev = (E_exp-clusE)/E_exp; //Energy deviation as percent of expected
	hEDiff_vs_X->Fill( HCALx, E_dev );
	hEDiff_vs_Y->Fill( HCALy, E_dev );
	
	//Build the matrix as simply as possible
	for(int icol = 0; icol<kNcell; icol++){
	  ba(icol)+= A[icol];
	  if(nblk==1){
	    bb(icol)+= A[icol];
	    oneBlock[icol]+= A[icol]*A[icol]/(KE_p*sampFrac);
	  }
	  for(int irow = 0; irow<kNcell; irow++){
	    Ma(icol,irow) += A[icol]*A[irow]/(KE_p*sampFrac);
	  
	  } 
	}
      }
    }
  }

  cout << endl << "Checking data, inverting matrix, and solving for coefficients.." << endl << endl;

  //Reject the bad cells and normalize the oneblock check
  int badcell[kNcell];
  int badcell_oneblock[kNcell];
  int cellBad = 0;
  double y[kNcell] = {0.0}; // For easy TGraphErrors build
  
  for(int i=0; i<kNcell; i++){
    badcell[i] = 0;
    y[i] = i;
    
    //Do not change ADC gain coeff if insufficient events or energy dep in cell
    if( NEV[i] < minEventPerCell || Ma(i,i) < 0.1*ba(i) ){ 

      cellBad = 1;

      double elemRatio = Ma(i,i)/ba(i);

      ba(i) = 1.0;  // Set RHS vector for cell i to 1.0 
      Ma(i,i) = 1.0; // Set diagonal element of Matrix M for cell i to 1.0 
      
      for(int j=0; j<kNcell; j++){
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
  }

  if( cellBad==0 ) cout << "No bad cells detected.." << endl << endl;

  //Perform same rejection on single block analysis
  for(int i=0; i<kNcell; i++){
    badcell_oneblock[i] = 0;
    if( NEV_oneblock[i] < minEventPerCell || oneBlock[i] < 0.1*bb[i] ){
      bb(i) = 1.0;
      oneBlock[i] = 1.0;
      badcell_oneblock[i] = 1;
    }
  }
  
  //Invert the matrix, solve for ratios
  TMatrixD M_inv = Ma.Invert();
  TVectorD Coeff = M_inv*ba; // Stays unmodified for reference
  double oneBlockCoeff[kNcell];

  for( int i=0; i<kNcell; i++ ){
    if(badcell_oneblock[i]==0){
      oneBlockCoeff[i] = gOldConst[i]*(bb[i]/oneBlock[i]);

    }else{
      oneBlockCoeff[i] = gOldConst[i];
    }
  }

  for(int i=0; i<kNcell; i++){
    if(badcell[i]==0){
      GCoeff[i]=gOldConst[i]*Coeff[i]; // The new gain coefficient is the old coefficient multiplied by the gain factor that we just solved for. We will take these coefficents as the input for another set of data to estimate the error per channel.
      GCoeff_divide[i]=GCoeff[i]/oneBlockCoeff[i];

    }else{
      GCoeff[i]=gOldConst[i]; // If the cell is bad, use the old coefficient
      GCoeff_divide[i]=-1.0;
    }
  }

  cout << "Inversion complete. Building histograms and writing out coefficients.." << endl << endl;

  for(int i=0; i<kNcell; i++){
    coefficients->Fill(i,GCoeff[i]);
    coefficients_1b->Fill(i,oneBlockCoeff[i]);
  }

  double yErr[kNcell] = {0.0};
  TGraphErrors *ccgraph = new TGraphErrors( kNcell, y, GCoeff, yErr, yErr ); 
  ccgraph->GetXaxis()->SetLimits(0.0,288.0);  
  ccgraph->GetYaxis()->SetLimits(0.0,0.25);
  ccgraph->SetTitle("Calibration Coefficients");
  ccgraph->GetXaxis()->SetTitle("Channel");
  ccgraph->GetXaxis()->SetTitle("Unitless");
  ccgraph->SetMarkerStyle(20); // idx 20 Circles, idx 21 Boxes
  ccgraph->Write("constants");

  TGraphErrors *ccgraph_oneBlock = new TGraphErrors( kNcell, y, oneBlockCoeff, yErr, yErr ); 
  ccgraph_oneBlock->GetXaxis()->SetLimits(0.0,288.0);  
  ccgraph_oneBlock->GetYaxis()->SetLimits(0.0,0.25);
  ccgraph_oneBlock->SetTitle("Calibration Coefficients One Block");
  ccgraph_oneBlock->GetXaxis()->SetTitle("Channel");
  ccgraph_oneBlock->GetXaxis()->SetTitle("Unitless");
  ccgraph_oneBlock->SetMarkerStyle(21);
  ccgraph_oneBlock->Write("constants_oneblock");

  TGraphErrors *ccgraph_divide = new TGraphErrors( kNcell, y, GCoeff_divide, yErr, yErr ); 
  ccgraph_divide->GetXaxis()->SetLimits(0.0,288.0);  
  ccgraph_divide->SetTitle("Calibration Coefficients / OneBlock Coeff");
  ccgraph_divide->GetXaxis()->SetTitle("Channel");
  ccgraph_divide->GetXaxis()->SetTitle("Unitless");
  ccgraph_divide->SetMarkerStyle(21);
  ccgraph_divide->Write("constants_divide");

  //Write out diagnostic histos and print to console
  fout->Write();

  cout << "Gain Coefficients: " << endl << endl;

  int cell = 0;

  for( int r = 0; r<kNrows; r++){
    for( int c = 0; c<kNcols; c++){
      if( GCoeff[cell]==1 ) GCoeff[cell] = 0.00175;
      cout << GCoeff[cell] << "  ";
      cell++;
    }
    cout << endl;
  }

  cell = 0;

  for( int i=0; i<kNcell; i++){
    GCoeff[i] = GCoeff[i] - 0.00175; //Cosmic value in mV/GeV for all channels
  }

  TGraphErrors *ccgraph_Cdiff = new TGraphErrors( kNcell, y, GCoeff, yErr, yErr ); 
  ccgraph_Cdiff->GetXaxis()->SetLimits(0.0,288.0);  
  ccgraph_Cdiff->GetYaxis()->SetLimits(0.0,0.025);
  ccgraph_Cdiff->SetTitle("Cosmic Delta Calibration Coefficients");
  ccgraph_Cdiff->GetXaxis()->SetTitle("Channel");
  ccgraph_Cdiff->GetXaxis()->SetTitle("Unitless");
  ccgraph_Cdiff->SetMarkerStyle(20); // idx 20 Circles, idx 21 Boxes
  ccgraph_Cdiff->Write("constants_diff");

  cout << endl << "Gain Coefficients diff from cosmic value: " << endl;

  for( int r = 0; r<kNrows; r++){
    for( int c = 0; c<kNcols; c++){
      cout << GCoeff[cell] << "  ";
      cell++;
    }
    cout << endl;
  }

  cell = 0;

  cout << endl << "Number of events available for calibration: " << endl << endl;

  for( int r = 0; r<kNrows; r++){
    for( int c = 0; c<kNcols; c++){
      cout << NEV[cell] << "  ";
      cell++;
    }
    cout << endl;
  }
  cell = 0;

  cout << endl << "One Block:" << endl;

  for( int r = 0; r<kNrows; r++){
    for( int c = 0; c<kNcols; c++){
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

  for( int i=0; i<kNcell; i++ ){   
    ratio << Coeff[i] << endl;
  }

  ratio.close();

  GainCoeff.open( constPath );
  GainCoeff << "#HCal gain coefficients from SBS-" << kine << " obtained " << date.c_str() << endl;
  GainCoeff << "#HCal_gainCoeff = " << endl;

  cell = 0;
  for( int r=0; r<kNrows; r++ ){
    for( int c=0; c<kNcols; c++ ){
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
  for( int r=0; r<kNrows; r++ ){
    for( int c=0; c<kNcols; c++ ){
      GainCoeff_oneblock << oneBlockCoeff[cell] << "  ";
      cell++;
    }
    GainCoeff_oneblock << endl;
  }
  cell = 0;

  GainCoeff_oneblock.close();

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
