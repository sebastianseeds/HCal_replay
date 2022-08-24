//SSeeds 8.12.22 - Post-production - Script to compare SBS-offline outputs and to verify fix to double peaking implemented by MKJ

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

//HCal position corrections - will add to DB file eventually
double hcalheight = 0.365; //m The height of the center of HCAL above beam

const double PI = TMath::Pi();
const double M_e = 0.00051;
const double M_p = 0.938272;

const double TDCCalib = 0.112;
//const double TDCCalib = 1.0;

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

void DP_comp( const char *configfilename="sPeak8.cfg", int run = -1 ){
  
  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );
  
  // Get the date
  string date = getDate();
  
  // Declare Chains for many root files of two sorts
  TChain *C = new TChain("T");
  C->Add("/lustre19/expphy/volatile/halla/sbs/seeds/rootfiles/hcal_general_13485_500001*");
  
  // Declare general physics parameters to be modified by input config file
  double test = 1; // Keep track of DB read for analysis. 0: ../SBS_REPLAY/SBS_Replay/DB   1: ../seeds/SBS-replay/DB
  double diag = 0; // Keep track of diagnostic canvas pdf option printed at end
  double kine = 8; // Keep track of kinematic calibrating from
  double spotcut = 0; // Option to cut on the expected location of elastic protons
  double waveform = 0; // Option to read in full HCal waveforms (MODE1)
  double cosmic = 0; // Option to configure for HCal cosmic data (HCal trig)
  double tFitMin = 30; // Minimum number of entries per channel to calibrate ADC/TDC time
  double t_trig = 510; // Mean tdc trig value (HCAL - BB) 
  double E_e = 1.92; // Energy of beam (incoming electrons from accelerator)
  int TDCmissC = 0; // Keep track of TDC misses
  int TDCmissD = 0; // Keep track of TDC misses
  double HCal_d = 14.5; // Distance to HCal from scattering chamber for comm1
  double HCal_th = 35.0; // Angle that the center of HCal is at  
  double W_mean = 0.93; // Mean of W at current kinematic
  double W_sig = 0.039; // Width of W at current kinematic
  double dx0 = 0.9; // Position of proton spot, x-x_expected
  double dy0 = 0.62; // Position of proton spot, y-y_expected
  double dx_sig = 0.09; // Max spread of proton spot, x-x_expected
  double dy_sig = 0.15; // Max spread of proton spot, y-y_expected
  int elasYieldC = 0; // Keep track of total elastics analyzed
  int elasYieldD = 0; // Keep track of total elastics analyzed
  
  // Declare arrays to hold old offset parameters
  double oldADCtoffsets[kNcell];
  double oldTDCoffsets[kNcell];
  
  cout << endl;
  
  // Reading config file
  ifstream configfile(configfilename);
  TString currentline;
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      if(!currentline) cout << "WARNING: No file exists at " << currentline << "." << endl;
      //C->Add(currentline);
      //cout << "Added file " << currentline << " to chain.." << endl;
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
      if( skey == "test" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	kine = sval.Atof();
	cout << "Loading test setting: " << test << endl;
      }
      if( skey == "diag" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	diag = sval.Atof();
	cout << "Loading diagnostic setting: " << diag << endl;
      }
      if( skey == "spotcut" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	spotcut = sval.Atof();
	cout << "Loading cut on hadron spot setting: " << spotcut << endl;
      }
      if( skey == "waveform" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	waveform = sval.Atof();
	cout << "Loading option to use full waveforms: " << waveform << endl;
      }
      if( skey == "cosmic" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	cosmic = sval.Atof();
	cout << "Loading HCal cosmic option: " << cosmic << endl;
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

  cout << endl << endl << "Setup parameters loaded." << endl;
  
  TEventList *elistC = new TEventList("elistC","Elastic Event List C");
  C->Draw(">>elistC",globalcut);
  
  cout << "Event lists populated with cut placed on elastics." << endl;
  
  // Declare general detector parameters
  double BBtr_p[maxTracks], BBtr_pz[maxTracks], BBtr_px[maxTracks], BBtr_py[maxTracks];
  double BBtr_vz[maxTracks], BBtr_chi2[maxTracks];
  double BBtr_x[maxTracks], BBtr_y[maxTracks];
  double BBtr_n, BBps_x, BBps_y, BBps_e, BBsh_x, BBsh_y, BBsh_e;
  
  double TDCT_id[maxTdcChan], TDCT_tdc[maxTdcChan]; 
  int TDCTndata;
  UInt_t runI, runN;
  UInt_t TBits;
  ULong64_t runT;
  
  double HCALx, HCALy, HCALe, HCALbit, HCALttb;
  double HCALtdc[kNcell], HCALtm[kNcell], HCALTtdc[kNcell], HCALa[kNcell], HCALap[kNcell];
  double crow, ccol, nblk;
  double cblkid[kNcell], cblke[kNcell];
  
  double tHODO;
  double nClusHODO;
  
  // Declare root tree variables and set values to memory locations in root file
  C->SetBranchStatus( "*", 0 );
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
  C->SetBranchStatus( "bb.tr.chi2", 1 );
  C->SetBranchStatus( "bb.tr.n", 1 );
  C->SetBranchStatus( "bb.tr.x", 1 );
  C->SetBranchStatus( "bb.tr.y", 1 );
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
  
  if( waveform==1 || cosmic==1 ){
    C->SetBranchStatus( "sbs.hcal.a_p", 1 );
    C->SetBranchStatus( "sbs.hcal.tdc_mult", 1 );
    C->SetBranchStatus( "sbs.hcal.ledbit", 1 );
    C->SetBranchStatus( "sbs.hcal.tdctimeblk", 1 );
    C->SetBranchStatus( "sbs.hcal.tdc", 1 );
  }
  
  C->SetBranchAddress( "sbs.hcal.x", &HCALx );
  C->SetBranchAddress( "sbs.hcal.y", &HCALy );
  C->SetBranchAddress( "sbs.hcal.e", &HCALe );
  C->SetBranchAddress( "sbs.hcal.rowblk", &crow );
  C->SetBranchAddress( "sbs.hcal.colblk", &ccol );
  C->SetBranchAddress( "sbs.hcal.nblk", &nblk ); // Total number of blocks in highest E clus
  C->SetBranchAddress( "sbs.hcal.clus_blk.id", cblkid ); // kNcell-1 index for each block
  C->SetBranchAddress( "sbs.hcal.clus_blk.e", cblke ); // Array of block energies
  C->SetBranchAddress( "sbs.hcal.clus_blk.tdctime", HCALtdc );
  C->SetBranchAddress( "sbs.hcal.clus_blk.atime", HCALa );
  C->SetBranchAddress( "bb.tr.n", &BBtr_n );
  C->SetBranchAddress( "bb.tr.chi2", BBtr_chi2 );
  C->SetBranchAddress( "bb.tr.x", BBtr_x );
  C->SetBranchAddress( "bb.tr.y", BBtr_y );
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
  
  if( waveform==1 || cosmic==1){
    C->SetBranchAddress( "sbs.hcal.a_p", HCALap );
    C->SetBranchAddress( "sbs.hcal.tdc_mult", HCALtm );
    C->SetBranchAddress( "sbs.hcal.ledbit", &HCALbit );
    C->SetBranchAddress( "sbs.hcal.tdctimeblk", &HCALttb );
    C->SetBranchAddress( "sbs.hcal.tdc", HCALTtdc);
  }

  C->SetBranchAddress( "bb.hodotdc.clus.tmean", &tHODO );
  C->SetBranchAddress( "bb.hodotdc.nclus", &nClusHODO );
  
  C->SetBranchAddress( "fEvtHdr.fRun", &runI );
  C->SetBranchAddress( "fEvtHdr.fEvtTime", &runT );
  C->SetBranchAddress( "fEvtHdr.fEvtNum", &runN );
  C->SetBranchAddress( "fEvtHdr.fTrigBits", &TBits );
  
  cout << "Tree variables linked." << endl;
  
  // Define a clock to check macro processing time
  TStopwatch *sw = new TStopwatch();
  sw->Start();
  
  // Declare outfile
  TString outputfilename = Form("DP_compOut_%s.root", date.c_str());
  TFile *fout = new TFile( outputfilename, "RECREATE" );
  

  // Initialize Histos
  TH2D *hHCALtdc_vs_HODOtdc_C = new TH2D("hHCALtdc_vs_HODOtdc_C",";TDC_{HODO} (ns);TDC_{HCAL} (ns)",150,-75,75,150,-150,0);
  TH2D *hHCALtdc_vs_HODOtdc_D = new TH2D("hHCALtdc_vs_HODOtdc_D",";TDC_{HODO} (ns);TDC_{HCAL} (ns)",150,-75,75,150,-150,0);
  TH2D *htHCAL_vs_HCALID_C = new TH2D("htHCAL_vs_HCALID_C",";Channel;TDC_{HCAL} (ns)",288,0,288,300,-150,0);
  TH2D *htHCAL_vs_HCALID_D = new TH2D("htHCAL_vs_HCALID_D",";Channel;TDC_{HCAL} (ns)",288,0,288,300,-150,0);
  TH2D *htcorr_HCAL_HODO_C = new TH2D("htcorr_HCAL_HODO_C",";TDC_{HCAL} (ns);TDC_{HODO} (ns)", 150, -150, 0, 150, -15, 15);  
  TH2D *htcorr_HCAL_HODO_D = new TH2D("htcorr_HCAL_HODO_D",";TDC_{HCAL} (ns);TDC_{HODO} (ns)", 150, -150, 0, 150, -15, 15);  
  TH1D *hclusblkcut_dt_C = new TH1D( "hclusblkcut_dt_C", "Time Difference Between Primary Hit and First Adjacent Block C", 4000, -1000, 1000);
  TH1D *hclusblkcut_dt_D = new TH1D( "hclusblkcut_dt_D", "Time Difference Between Primary Hit and First Adjacent Block D", 4000, -1000, 1000);
  TH1D *hclusblk_dt_C = new TH1D( "hclusblk_dt_C", "Time Difference Between Blocks in Cluster C (ns)", 1000, -50, 50);
  TH1D *hclusblk_dt_D = new TH1D( "hclusblk_dt_D", "Time Difference Between Blocks in Cluster D (ns)", 1000, -50, 50);

  cout << "Variables and histograms defined." << endl;
  
  // Set long int to keep track of total entries
  Long64_t NeventsC = elistC->GetN();

  UInt_t run_number = 0;
  
  cout << endl << "All parameters loaded and initialization complete." << endl << endl;
  cout << "Opened up tree with nentries on C: " << C->GetEntries() << ", nentries passing globalcut: " << NeventsC << "." << endl << endl;
  
  //First loop over events to obtain mean TOF and timewalk params by channel
  cout << "First loop over data subset commencing.." << endl;
  Double_t progress = 0.;
  Double_t timekeeper = 0., timeremains = 0.;
  while(progress<1.0){
    Int_t barwidth = 70;
    for(Long64_t nevent = 0; nevent<NeventsC; nevent++){
      
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
	timeremains = timekeeper*( double(NeventsC)/double(nevent) - 1. ); 
      sw->Reset();
      sw->Continue();
      
      progress = (double)((nevent+1.)/NeventsC);
      cout << "] " << int(progress*100.) << "%, elastic events: " << elasYieldC << ", time remaining: " << int(timeremains/60.) << "m \r";
      cout.flush();
      
      C->GetEntry( elistC->GetEntry( nevent ) ); 

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
      double pp = sqrt(pow(nu,2)+2.*M_p*nu); //momentum of the proton
      double phinucleon = ephi + PI; //assume coplanarity
      double thetanucleon = acos( (E_e - BBtr_p[0]*cos(etheta))/pp ); //use elastic constraint on nucleon kinematics
      
      TVector3 pNhat( sin(thetanucleon)*cos(phinucleon),sin(thetanucleon)*sin(phinucleon),cos(thetanucleon));
      
      //Define HCal coordinate system
      TVector3 HCAL_zaxis(sin(-HCal_th),0,cos(-HCal_th));
      TVector3 HCAL_xaxis(0,-1,0);
      TVector3 HCAL_yaxis = HCAL_zaxis.Cross(HCAL_xaxis).Unit();
      
      TVector3 HCAL_origin = HCal_d * HCAL_zaxis + hcalheight * HCAL_xaxis;
            
      //Define intersection points for hadron vector
      double sintersect = ( HCAL_origin - vertex ).Dot( HCAL_zaxis ) / ( pNhat.Dot( HCAL_zaxis ) );
      
      TVector3 HCAL_intersect = vertex + sintersect * pNhat;
      
      double yexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_yaxis );
      double xexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_xaxis );
      
      double E_ep = sqrt( pow(M_e,2) + pow(BBtr_p[0],2) ); // Obtain the scattered electron energy
      
      double p_ep = BBtr_p[0];
      double Q2 = 2*E_e*E_ep*( 1-(BBtr_pz[0]/p_ep) ); // Obtain Q2 from beam energy, outgoing electron energy, and momenta
      
      double W = PgammaN.M();
      
      //Use the electron kinematics to predict the proton momentum assuming elastic scattering on free proton at rest:
      double E_pp = nu+M_p; // Get total energy of the proton
      
      double KE_p = nu; //For elastics
      
      //Cut on BBCal and HCal trigger coincidence and plot all other trigger timing. All ref to tdcs need index adjustment (ihit+1).
      double bbcal_time=0., hcal_time=0., RF_time=0.; 
      double bbcalLO_time=0., bbcalHIveto_time=0., edtm_time=0.;
      
      for(int ihit=0; ihit<TDCTndata; ihit++){
	//cout << "TDCT_id[ihit]:" << TDCT_id[ihit] << ", TDCT_tdc:" << TDCT_tdc[ihit+1] << endl;
	if(TDCT_id[ihit]==5) bbcal_time=TDCT_tdc[ihit+1];
	if(TDCT_id[ihit]==4) RF_time=TDCT_tdc[ihit+1];
	if(TDCT_id[ihit]==3) edtm_time=TDCT_tdc[ihit+1];
	if(TDCT_id[ihit]==2) bbcalLO_time=TDCT_tdc[ihit+1];
	if(TDCT_id[ihit]==1) bbcalHIveto_time=TDCT_tdc[ihit+1];
	if(TDCT_id[ihit]==0) hcal_time=TDCT_tdc[ihit+1];
      }
      
      double diff = hcal_time - bbcal_time; 
      double RFmod = std::fmod(RF_time,4); //RF Signature measures the bunch crossing of beam electrons at 4ns intervals

      double xDiff = HCALx - xexpect_HCAL;
      double yDiff = HCALy - yexpect_HCAL;
      
      //Check how often the tdc failed to register for an otherwise good event
      if( HCALe>0.02 && HCALtdc[0]<-400 && nClusHODO<10 ) TDCmissC++;
      
      bool nogo = spotcut==1 && (fabs(xDiff-dx0)>dx_sig || fabs(yDiff-dy0)>dy_sig);
      
      //Primary cut on elastics
      if( fabs( W-W_mean )<W_sig ){
	
	//Get time difference between ADCt/TDC and HODO TDC mean time
	double HHdiff = HCALtdc[0]-tHODO; //time difference between HCal TDC and Hodo TDC
	double HBBdiff = HCALtdc[0]-bbcal_time; //time difference between HCal TDC and BBCal trigger
	double HBBdiff_rf = HCALtdc[0]-bbcal_time-RFmod; //time difference between HCal TDC and BBCal trigger
	double HRFdiff = HCALtdc[0]-RFmod; //time difference between HCal TDC and RF time modulo 4ns
	double HHAdiff = HCALa[0]-tHODO; //time difference between HCal ADC time and Hodo TDC
	double e = cblke[0]; //energy deposited in primary block
	int idx = (int)cblkid[0]-1; //index of primary block
	int col = (int)ccol; //column of primary block
	int row = (int)crow;
	if( idx<0 || idx>=288 ) cout << "ERROR: indexing out of bounds at idx = " << idx << endl;
	
	//Fill some unmodified time difference histos
	hHCALtdc_vs_HODOtdc_C->Fill( tHODO, HCALtdc[0] );
	htHCAL_vs_HCALID_C->Fill( cblkid[0], HCALtdc[0] );
	htcorr_HCAL_HODO_C->Fill( HCALtdc[0], tHODO );
	
	double dt = HCALtdc[0] - HCALtdc[1];

	// Check if next highest energy block after primary is immediately adjacent to primary block
	if( (cblkid[1]==cblkid[0]+1)||
	    (cblkid[1]==cblkid[0]-1)||
	    (cblkid[1]==cblkid[0]+12)||
	    (cblkid[1]==cblkid[0]-12) ){
	  hclusblkcut_dt_C->Fill( dt );
	}

	for( int i=0; i<nblk; i++){
	  if( i>0 ) hclusblk_dt_C->Fill( HCALtdc[0] - HCALtdc[i] );
	}

	//Events that pass the above cuts constitute elastics
	elasYieldC++;
      }
    }
  }

  cout << endl << "Loop over first data set complete. Proceeding with second data set." << endl << endl;

  //C->Delete();
  
  //TChain *D = new TChain("T");
  C->Add("/lustre19/expphy/volatile/halla/sbs/seeds/rootfiles/hcal_general_13485_500005*");

  //TEventList *elistD = new TEventList("elistD","Elastic Event List D");
  C->Draw(">>elistC",globalcut);

  Long64_t NeventsD = elistC->GetN();

  cout << "Opened up tree, now with nentries: " << C->GetEntries() << ", nentries passing globalcut: " << NeventsD << "." << endl << endl;

  //Second loop over events to obtain mean TOF and timewalk params by channel
  cout << "Second loop over all data commencing.." << endl;
  run_number = 0;
  progress = 0.;
  timekeeper = 0.; 
  timeremains = 0.;
  while(progress<1.0){
    Int_t barwidth = 70;
    for(Long64_t nevent = NeventsC; nevent<NeventsD; nevent++){
      
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
	timeremains = timekeeper*( double(NeventsD)/double(nevent) - 1. ); 
      sw->Reset();
      sw->Continue();
      
      progress = (double)((nevent+1.)/NeventsD);
      cout << "] " << int(progress*100.) << "%, elastic events: " << elasYieldD << ", time remaining: " << int(timeremains/60.) << "m \r";
      cout.flush();
      
      C->GetEntry( elistC->GetEntry( nevent ) ); 

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
      double pp = sqrt(pow(nu,2)+2.*M_p*nu); //momentum of the proton
      double phinucleon = ephi + PI; //assume coplanarity
      double thetanucleon = acos( (E_e - BBtr_p[0]*cos(etheta))/pp ); //use elastic constraint on nucleon kinematics
      
      TVector3 pNhat( sin(thetanucleon)*cos(phinucleon),sin(thetanucleon)*sin(phinucleon),cos(thetanucleon));
      
      //Define HCal coordinate system
      TVector3 HCAL_zaxis(sin(-HCal_th),0,cos(-HCal_th));
      TVector3 HCAL_xaxis(0,-1,0);
      TVector3 HCAL_yaxis = HCAL_zaxis.Cross(HCAL_xaxis).Unit();
      
      TVector3 HCAL_origin = HCal_d * HCAL_zaxis + hcalheight * HCAL_xaxis;
            
      //Define intersection points for hadron vector
      double sintersect = ( HCAL_origin - vertex ).Dot( HCAL_zaxis ) / ( pNhat.Dot( HCAL_zaxis ) );
      
      TVector3 HCAL_intersect = vertex + sintersect * pNhat;
      
      double yexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_yaxis );
      double xexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_xaxis );
      
      double E_ep = sqrt( pow(M_e,2) + pow(BBtr_p[0],2) ); // Obtain the scattered electron energy
      
      double p_ep = BBtr_p[0];
      double Q2 = 2*E_e*E_ep*( 1-(BBtr_pz[0]/p_ep) ); // Obtain Q2 from beam energy, outgoing electron energy, and momenta
      
      double W = PgammaN.M();
      
      //Use the electron kinematics to predict the proton momentum assuming elastic scattering on free proton at rest:
      double E_pp = nu+M_p; // Get total energy of the proton
      
      double KE_p = nu; //For elastics
      
      //Cut on BBCal and HCal trigger coincidence and plot all other trigger timing. All ref to tdcs need index adjustment (ihit+1).
      double bbcal_time=0., hcal_time=0., RF_time=0.; 
      double bbcalLO_time=0., bbcalHIveto_time=0., edtm_time=0.;
      
      for(int ihit=0; ihit<TDCTndata; ihit++){
	//cout << "TDCT_id[ihit]:" << TDCT_id[ihit] << ", TDCT_tdc:" << TDCT_tdc[ihit+1] << endl;
	if(TDCT_id[ihit]==5) bbcal_time=TDCT_tdc[ihit+1];
	if(TDCT_id[ihit]==4) RF_time=TDCT_tdc[ihit+1];
	if(TDCT_id[ihit]==3) edtm_time=TDCT_tdc[ihit+1];
	if(TDCT_id[ihit]==2) bbcalLO_time=TDCT_tdc[ihit+1];
	if(TDCT_id[ihit]==1) bbcalHIveto_time=TDCT_tdc[ihit+1];
	if(TDCT_id[ihit]==0) hcal_time=TDCT_tdc[ihit+1];
      }
      
      double diff = hcal_time - bbcal_time; 
      double RFmod = std::fmod(RF_time,4); //RF Signature measures the bunch crossing of beam electrons at 4ns intervals

      double xDiff = HCALx - xexpect_HCAL;
      double yDiff = HCALy - yexpect_HCAL;
      
      //Check how often the tdc failed to register for an otherwise good event
      if( HCALe>0.02 && HCALtdc[0]<-400 && nClusHODO<10 ) TDCmissD++;
            
      //Primary cut on elastics
      if( fabs( W-W_mean )<W_sig ){
	
	//Get time difference between ADCt/TDC and HODO TDC mean time
	double HHdiff = HCALtdc[0]-tHODO; //time difference between HCal TDC and Hodo TDC
	double HBBdiff = HCALtdc[0]-bbcal_time; //time difference between HCal TDC and BBCal trigger
	double HBBdiff_rf = HCALtdc[0]-bbcal_time-RFmod; //time difference between HCal TDC and BBCal trigger
	double HRFdiff = HCALtdc[0]-RFmod; //time difference between HCal TDC and RF time modulo 4ns
	double HHAdiff = HCALa[0]-tHODO; //time difference between HCal ADC time and Hodo TDC
	double e = cblke[0]; //energy deposited in primary block
	int idx = (int)cblkid[0]-1; //index of primary block
	int col = (int)ccol; //column of primary block
	int row = (int)crow;
	if( idx<0 || idx>=288 ) cout << "ERROR: indexing out of bounds at idx = " << idx << endl;
	
	//Fill some unmodified time difference histos
	hHCALtdc_vs_HODOtdc_D->Fill( tHODO, HCALtdc[0] );
	htHCAL_vs_HCALID_D->Fill( cblkid[0], HCALtdc[0] );
	htcorr_HCAL_HODO_D->Fill( HCALtdc[0], tHODO );
	
	double dt = HCALtdc[0] - HCALtdc[1];

	// Check if next highest energy block after primary is immediately adjacent to primary block
	if( (cblkid[1]==cblkid[0]+1)||
	    (cblkid[1]==cblkid[0]-1)||
	    (cblkid[1]==cblkid[0]+12)||
	    (cblkid[1]==cblkid[0]-12) ){
	  hclusblkcut_dt_D->Fill( dt );
	}

	for( int i=0; i<nblk; i++){
	  if( i>0 ) hclusblk_dt_D->Fill( HCALtdc[0] - HCALtdc[i] );
	}

	//Events that pass the above cuts constitute elastics
	elasYieldD++;
      }
    }
  }


  cout << "Analysis complete. Time elapsed: " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;
  
  cout << endl;
  
  //Write out diagnostic histos and print to console
  fout->Write();
  
  cout << endl << endl << "Elastic yield for analyzed runs: " << elasYieldC << ". Total events analyzed: " << NeventsC << ". Total TDC misses: " << TDCmissC << "." << endl << endl;
  cout << endl << endl << "Elastic yield for analyzed runs: " << elasYieldD << ". Total events analyzed: " << NeventsD << ". Total TDC misses: " << TDCmissD << "." << endl << endl;
  
  cout << "Timing offset analysis complete and parameters written to file." << endl;
  
  elistC->Delete();

  st->Stop();
  
  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;
  
}
