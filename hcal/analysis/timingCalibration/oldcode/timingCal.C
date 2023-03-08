//SSeeds 3.22.22 - Post-production - Calibration code which employs best current cuts on elastic events to obtain timing offsets for HCal ADC time and TDC

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
//The following are the positions of the "first" row and column from HCAL database (top right block as viewed from upstream)
double xoff_hcal = 0.92835;
double yoff_hcal = 0.47305;
double blockspace_hcal = 0.15254;

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

void timingCal( const char *configfilename="setup_timingCal.cfg", int run = -1 ){
  
  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  // Get the date
  string date = getDate();
  
  // Declare Chain for many root files
  TChain *C = new TChain("T");

  // Path for output file
  string paramsPath = "/w/halla-scshelf2102/sbs/seeds/HCal_replay/hcal/Hodo_HCal_Coorelations/timingCalParams.txt";

  // Declare general physics parameters to be modified by input config file
  double test = 1; // Keep track of DB read for analysis. 0: ../SBS_REPLAY/SBS_Replay/DB   1: ../seeds/SBS-replay/DB
  double diag = 0; // Keep track of diagnostic canvas pdf option printed at end
  double kine = 8; // Keep track of kinematic calibrating from
  double tFitMin = 30; // Minimum number of entries per channel to calibrate ADC/TDC time
  double E_e = 1.92; // Energy of beam (incoming electrons from accelerator)
  int TDCmiss = 0; // Keep track of TDC misses
  double HCal_d = 14.5; // Distance to HCal from scattering chamber for comm1
  double HCal_th = 35.0; // Angle that the center of HCal is at  
  double W_mean = 0.93; // Mean of W at current kinematic
  double W_sig = 0.039; // Width of W at current kinematic
  double dx0 = 0.9; // Position of proton spot, x-x_expected
  double dy0 = 0.62; // Position of proton spot, y-y_expected
  double dx_sig = 0.09; // Max spread of proton spot, x-x_expected
  double dy_sig = 0.15; // Max spread of proton spot, y-y_expected
  int elasYield = 0; // Keep track of total elastics analyzed

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
      C->Add(currentline);
      cout << "File loaded at " << currentline << endl;
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
  
  // Path for previous ADCt and TDC constants
  string inConstPath;
  if( test==0 ){
    inConstPath = "/w/halla-scshelf2102/sbs/SBS_REPLAY/SBS-replay/DB/db_sbs.hcal.dat";
  }else if( test==1 ){
    inConstPath = "/w/halla-scshelf2102/sbs/seeds/SBS-replay/DB/db_sbs.hcal.dat";
  }else{
    cout << "ERROR: No database file where specified" << endl;
    return 0;
  }

  // Reading ADC and TDC timing offsets from database
  cout << "Loading previous offsets.." << endl;
  ifstream inConstFile( inConstPath );
  if( !inConstFile ){
    cerr << endl << "ERROR: No input constant file present -> path to db_sbs.hcal.dat expected." << endl;
    return 0;
  }

  int n1=0;
  double d1=0;
  string line;
  bool skip_line = true;
  bool skip_one_line = true;
  bool pass_first_cond = false;
  bool pass_all_cond = false;

  while( getline( inConstFile, line ) ){

    if( pass_first_cond && n1==( kNcell ) ) break;

    TString Tline = (TString)line;

    if( Tline.BeginsWith("sbs.hcal.adc.timeoffset") && skip_line==true ) skip_line = false;
    if( Tline.BeginsWith("sbs.hcal.tdc.offset") && skip_line==true ) skip_line = false;

    if( skip_line==false && skip_one_line==true ){
      skip_one_line = false;
      continue;
    }

    if( n1==( kNcell ) && pass_first_cond==false ){
      skip_line = true;
      skip_one_line = true;
      n1=0;
      pass_first_cond = true;
    }

    if( skip_line==false && skip_one_line==false ){
      istringstream iss( line );
      while( iss >> d1 ){
	if( pass_first_cond==false ){
	  oldADCtoffsets[n1] = d1;
	}else{
	  oldTDCoffsets[n1] = d1;
	}
	n1++;
      }
    }
  }

  cout << endl << endl << "Old TDC offsets: " << endl;

  for( int r=0; r<kNrows; r++){
    for( int c=0; c<kNcols; c++){
      int i = r*kNcols+c;
      cout << oldTDCoffsets[i] << " ";
    }
    cout << endl;
  }

  cout << endl << endl << "Old ADC time offsets: " << endl;

  for( int r=0; r<kNrows; r++){
    for( int c=0; c<kNcols; c++){
      int i = r*kNcols+c;
      cout << oldADCtoffsets[i] << " ";
    }
    cout << endl;
  }

  cout << endl << endl << "Setup parameters loaded." << endl;

  TEventList *elist = new TEventList("elist","Elastic Event List");
  C->Draw(">>elist",globalcut);

  cout << "Event list populated with cut placed on elastics." << endl;

  // Declare general detector parameters
  double BBtr_p[maxTracks], BBtr_px[maxTracks], BBtr_py[maxTracks], BBtr_pz[maxTracks];
  double BBtr_vz[maxTracks], BBtr_chi2[maxTracks];
  double BBtr_n, BBps_x, BBps_y, BBps_e, BBsh_x, BBsh_y, BBsh_e;

  double TDCT_id[maxTdcChan], TDCT_tdc[maxTdcChan]; 
  int TDCTndata;
  UInt_t runI, runN;
  ULong64_t runT;

  double HCALx, HCALy, HCALe;
  double HCALtdc[kNcell], HCALa[kNcell];
  double crow, ccol, nblk;
  double cblkid[kNcell], cblke[kNcell];

  double tHODO;
  double nClusHODO;

  double kineW2;
  
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
  C->SetBranchStatus( "e.kine.W2", 1 );

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

  C->SetBranchAddress( "bb.hodotdc.clus.tmean", &tHODO );
  C->SetBranchAddress( "bb.hodotdc.nclus", &nClusHODO );

  C->SetBranchAddress( "fEvtHdr.fRun", &runI );
  C->SetBranchAddress( "fEvtHdr.fEvtTime", &runT );
  C->SetBranchAddress( "fEvtHdr.fEvtNum", &runN );
  C->SetBranchAddress( "e.kine.W2", &kineW2 );

  cout << "Tree variables linked." << endl;

  // Define a clock to check macro processing time
  TStopwatch *sw = new TStopwatch();
  sw->Start();
  
  // Declare outfile
  TString outputfilename = Form("timingCal_out_%d_corr.root", (int)kine);
  TFile *fout = new TFile( outputfilename, "RECREATE" );
  
  // Initialize vectors and arrays
  double TDCoffsets[kNcell] = {0.0};
  double combTDCoffsets[kNcell] = {0.0};
  double ADCtoffsets[kNcell] = {0.0};
  double combADCtoffsets[kNcell] = {0.0};
  double TDCsig[kNcell] = {0.0};
  double ADCtsig[kNcell] = {0.0};
  double TvsE[kNcell] = {0.0};
  double AtvsE[kNcell] = {0.0};

  // Initialize histograms	
  TH1D *hTOFcorr = new TH1D( "hTOFcorr", "Time of Flight Corrections; ns", 200, 0, 1 );
  TH1D *hCblkID = new TH1D( "hCblkID", "Block ID", 300, -10, 290 );
  TH1D *hElas_vp = new TH1D( "hElas_vp", "Elastic Proton Velocity", 400, 0, 4);
  TH1D *hCCol = new TH1D( "hCCol", "Cluster Seed Column", 15, -1, 14 );
  TH1D *hDeltaE = new TH1D( "hDeltaE","1.0-Eclus/p_rec", 100, -1.5, 1.5 );
  TH1D *hHCALe = new TH1D( "hHCALe","E HCal Cluster E", 400, 0., 4 );
  TH1D *hHODOnclus = new TH1D( "hHODOnclus","Number of Hodoscope Clusters", 50, 0., 50. );
  TH1D *hDiff = new TH1D( "hDiff","HCal time - BBCal time (ns)", 1300, -500, 800 );
  TH1D *hClusE = new TH1D( "hClusE","Best Cluster Energy", 100, 0.0, 2.0);
  TH2D *hPAngleCorr = new TH2D( "hPAngCorr","Track p vs Track ang", 100, 30, 60, 100, 0.4, 1.2 );
  TH2D *hPAngleCorr_2 = new TH2D( "hPAngCorr_2","Track p vs Track ang v2", 100, 30, 60, 100, 0.4, 1.2 );
  TH1D *hW = new TH1D( "W", "W; GeV", 250, 0.3, 1.5 );
  //hW->GetXaxis()->SetTitle( "GeV" );
  TH1D *hNBlk = new TH1D( "hNBlk", "Number of Blocks in Primary Cluster; Number", 25, 0, 25 );
  //hW->GetXaxis()->SetTitle( "Number" );
  TH1D *hW_cuts = new TH1D( "W_cuts", "W_cuts; GeV", 250, 0.3, 1.5 );
  //hW->GetXaxis()->SetTitle( "GeV" );
  TH1D *hQ2 = new TH1D( "Q2", "Q2; GeV", 250, 0.5, 3.0 );
  //hQ2->GetXaxis()->SetTitle( "GeV" );
  TH1D *hE_ep = new TH1D( "Scattered Electron Energy; GeV","E_ep", 500, 0.0, E_e*1.5 ); 
  //hE_ep->GetXaxis()->SetTitle( "GeV" );
  TH1D *hE_pp = new TH1D( "Scattered Proton Energy; GeV", "E_pp", 500, 0.0, E_e*1.5 );
  //hE_pp->GetXaxis()->SetTitle( "GeV" );
  TH1D *hKE_p = new TH1D( "Scattered Proton Kinetic Energy; GeV", "KE_pp", 500, 0.0, E_e*1.5 );
  //hKE_p->GetXaxis()->SetTitle( "GeV" );
  TH2D *hADC = new TH2D( "hADC", "HCal Int_ADC Spectra: W Cut", 288, 0, 288, 100., 0., 1. );
  TH2D *hTOFvID = new TH2D( "hTOFvID", "Time of Flight vs HCal ID", 288, 0, 288, 200., 0., 20. );
  TH2D *hADC_amp = new TH2D( "hADC_amp", "HCal ADC_amp Spectra: W Cut", 288, 0, 288, 100., 0., 10. );
  TH2D *hdxdy_HCAL = new TH2D("hdxdy_HCAL",";y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)", 250, -5.0, 5.0, 250, -10, 10 );
  TH2D *hXY_HCAL_ps = new TH2D("hXY_HCAL_ps",";y_{HCAL} (m); x_{HCAL} (m)", 250, -5.0, 5.0, 250, -10, 10 );
  TH2D *hXY_HCAL = new TH2D("hXY_HCAL",";y_{HCAL} (m); x_{HCAL} (m)", 250, -5.0, 5.0, 250, -10, 10 );
  TH1D *hvz = new TH1D("hvz",";vertex z (m);", 250,-0.125,0.125);
  TH1D *hvz_cut = new TH1D("hvz_cut",";vertex z (m);", 250,-0.125,0.125);
 
  // Timing histograms
  TH2D *htcorr_HCAL_HODO = new TH2D("htcorr_HCAL_HODO",";TDC_{HCAL} (ns);TDC_{HODO} (ns)", 150, -150, 0, 150, -15, 15);
  TH1D *htDiff_HODO_HCAL = new TH1D("htDiff_HODO_HCAL","",150,-150,0);
  TH1D *haDiff_HODO_HCAL = new TH1D("haDiff_HODO_HCAL","",300,0,150);
  TH1D *haDiff_HODO_HCAL_JLAB = new TH1D("haDiff_HODO_HCAL_JLAB","",300,0,150);
  TH1D *haDiff_HODO_HCAL_CMU = new TH1D("haDiff_HODO_HCAL_CMU","",300,0,150);
  TH1D *htDiff_HODO_HCAL_JLAB = new TH1D("htDiff_HODO_HCAL_JLAB","",300,-150,0);
  TH1D *htDiff_HODO_HCAL_CMU = new TH1D("htDiff_HODO_HCAL_CMU","",300,-150,0);
  TH2D *htDiff_vs_ADCint_JLAB = new TH2D("htDiff_vs_ADCint_JLAB",";E_{JLAB} (GeV);TDC_{HCAL}-TDC_{HODO} (ns)",270,0.0,0.9,600,-250,50);
  TH2D *haDiff_vs_ADCint_JLAB = new TH2D("haDiff_vs_ADCint_JLAB",";E_{JLAB} (GeV);ADCt_{HCAL}-TDC_{HODO} (ns)",270,0.0,0.9,600,0,150);
  TH2D *htDiff_vs_ADCint_CMU = new TH2D("htDiff_vs_ADCint_CMU",";E_{CMU} (GeV);TDC_{HCAL}-TDC_{HODO} (ns)",270,0.0,0.9,600,-250,50);
  TH2D *haDiff_vs_ADCint_CMU = new TH2D("haDiff_vs_ADCint_CMU",";E_{CMU} (GeV);ADCt_{HCAL}-TDC_{HODO} (ns)",270,0.0,0.9,600,0,150);
  TH1D *hTDCTimewalk = new TH1D("hTDCTimewalk","",500,0,50);
  TH1D *hADCtTimewalk = new TH1D("hADCtTimewalk","",500,0,50);
  TH2D *htDiff_vs_HCALID = new TH2D("htDiff_vs_HCALID",";Channel;TDC_{HCAL}-TDC_{HODO} (ns)",288,0,288,300,-150,0);
  TH2D *haDiff_vs_HCALID = new TH2D("haDiff_vs_HCALID",";Channel;ADCtime_{HCAL}-TDC_{HODO} (ns)",288,0,288,300,0,150);
  TH2D *hTDCoffsets_vs_HCALID = new TH2D("hTDCoffsets_vs_HCALID",";Channel;TDC Offset (ns)",288,0,288,200,-50,50);
  TH2D *holdTDCoffsets_vs_HCALID = new TH2D("holdTDCoffsets_vs_HCALID",";Channel;TDC Offset (ns)",288,0,288,200,-50,50);
  TH2D *hcombTDCoffsets_vs_HCALID = new TH2D("hcombTDCoffsets_vs_HCALID",";Channel;TDC Offset (ns)",288,0,288,200,-50,50);
  TH2D *hADCtoffsets_vs_HCALID = new TH2D("hADCtoffsets_vs_HCALID",";Channel;ADC Time Offset (ns)",288,0,288,200,-50,50);
  TH2D *hTDCsig_vs_HCALID = new TH2D("hTDCsig_vs_HCALID",";Channel;TDC Std Dev (ns)",288,0,288,100,0,10);
  TH2D *hADCtsig_vs_HCALID = new TH2D("hADCtsig_vs_HCALID",";Channel;ADC Time Std Dev (ns)",288,0,288,500,0,50);
  TH2D *htDiff_vs_HCALID_corr = new TH2D("htDiff_vs_HCALID_corr",";Channel;TDC_{HCAL}-TDC_{HODO} (ns)",288,0,288,300,-150,0);
  TH2D *htDiff_vs_ADCint[kNcell+1];
  TH2D *haDiff_vs_ADCint[kNcell+1];
  TH1D *htDiff[kNcell+1];
  TH1D *haDiff[kNcell+1];
  TH1D *haDiff_vs_col[kNcols];
  TH1D *htDiff_vs_col[kNcols];
  TH2D *htDiff_vs_ADCint_cols[kNcols];
  TH2D *haDiff_vs_ADCint_cols[kNcols];

  cout << endl;

  for( int i=0; i<kNcols; i++){
    haDiff_vs_col[i] = new TH1D(Form("haDiff_col%d",i),";ADCt_{HCAL}-TDC_{HODO} (ns)",300,0,150);
    htDiff_vs_col[i] = new TH1D(Form("htDiff_col%d",i),";ADCt_{HCAL}-TDC_{HODO} (ns)",300,-150,0);
    haDiff_vs_ADCint_cols[i] = new TH2D(Form("htDiff_vs_ADCint_col%d",i),Form(";E_{col%d} (GeV);TDC_{HCAL}-TDC_{HODO} (ns)",i),270,0.0,0.9,600,-250,50);
    htDiff_vs_ADCint_cols[i] = new TH2D(Form("haDiff_vs_ADCint_col%d",i),Form(";E_{col%d} (GeV);ADCt_{HCAL}-TDC_{HODO} (ns)",i),270,0.0,0.9,600,0,150);
  }

  for( int i=0; i<kNcell; i++ ){
    htDiff_vs_ADCint[i] = new TH2D(Form("htDiff_vs_ADCint_bl%d",i),Form(";E_{bl%d} (GeV);TDC_{HCAL}-TDC_{HODO} (ns)",i),270,0.0,0.9,600,-250,50);
    haDiff_vs_ADCint[i] = new TH2D(Form("haDiff_vs_ADCint_bl%d",i),Form(";E_{bl%d} (GeV);ADCt_{HCAL}-TDC_{HODO} (ns)",i),270,0.0,0.9,600,0,150);
    haDiff_vs_ADCint[i]->Fill( 0.1, 55 ); //To test fits

    htDiff[i] = new TH1D(Form("htDiff_bl%d",i),";TDC_{HCAL}-TDC_{HODO} (ns)",300,-150,0);
    haDiff[i] = new TH1D(Form("haDiff_bl%d",i),";ADCt_{HCAL}-TDC_{HODO} (ns)",1400,-400,300);
    haDiff[i]->Fill( 55 ); //To test fits
  }

  cout << "Variables and histograms defined." << endl;

  // Set long int to keep track of total entries
  Long64_t Nevents = elist->GetN();
  UInt_t run_number = 0;

  cout << endl << "All parameters loaded and initialization complete." << endl << endl;
  cout << "Opened up tree with nentries: " << C->GetEntries() << ", nentries passing globalcut: " << Nevents << "." << endl << endl;

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
		
      if( run_number!=runI ){
	run_number=runI;
	//cout << "Now analyzing run " << run_number << "." << endl;
      }

      //Sort all tracks by lowest Chi^2 (highest confidence track)
      int track_tot = (int)BBtr_n;
      int track = 0;
      double min = 1000.0; // Set arbitrarily high chi^2 to minimize with loop over tracks
      for( int elem=0; elem<track_tot; elem++ ){            
	if( BBtr_chi2[elem] < min )
	  track = elem;      
      }

      double etheta = acos( BBtr_pz[track]/BBtr_p[track] );
      double ephi = atan2( BBtr_py[track], BBtr_px[track] );

      TVector3 vertex(0,0,BBtr_vz[track]); // z location of vertex in hall coordinates
      TLorentzVector Pbeam(0,0,E_e,E_e); //Mass of e negligable
      TLorentzVector kprime(BBtr_px[track],BBtr_py[track],BBtr_pz[track],BBtr_p[track]);
      TLorentzVector Ptarg(0,0,0,M_p);

      TLorentzVector q = Pbeam - kprime;
      TLorentzVector PgammaN = Ptarg + q; //(-px, -py, ebeam - pz, Mp + ebeam - p)

      double pel = E_e/(1.+E_e/M_p*(1.-cos(etheta)));
      double nu = E_e - BBtr_p[track];
      double pp = sqrt(pow(nu,2)+2.*M_p*nu); //momentum of the proton
      double phinucleon = ephi + TMath::Pi(); //assume coplanarity
      double thetanucleon = acos( (E_e - BBtr_p[track]*cos(etheta))/pp ); //use elastic constraint on nucleon kinematics
	
      TVector3 pNhat( sin(thetanucleon)*cos(phinucleon),sin(thetanucleon)*sin(phinucleon),cos(thetanucleon));

      //Define HCal coordinate system
      TVector3 HCAL_zaxis(sin(-HCal_th),0,cos(-HCal_th));
      TVector3 HCAL_xaxis(0,-1,0);
      TVector3 HCAL_yaxis = HCAL_zaxis.Cross(HCAL_xaxis).Unit();
	
      TVector3 HCAL_origin = HCal_d * HCAL_zaxis + hcalheight * HCAL_xaxis;

      //Sanity check on position after offsets applied
      hXY_HCAL->Fill( HCALy, HCALx );

      //Define intersection points for hadron vector
      double sintersect = ( HCAL_origin - vertex ).Dot( HCAL_zaxis ) / ( pNhat.Dot( HCAL_zaxis ) );

      TVector3 HCAL_intersect = vertex + sintersect * pNhat;

      double yexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_yaxis );
      double xexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_xaxis );

      //Calculate the proton spot - use for cut later on
      hdxdy_HCAL->Fill( HCALy - yexpect_HCAL, HCALx - xexpect_HCAL );

      double E_ep = sqrt( pow(M_e,2) + pow(BBtr_p[track],2) ); // Obtain the scattered electron energy
      hE_ep->Fill( E_ep ); // Fill histogram
	
      double p_ep = BBtr_p[track];
      double Q2 = 2*E_e*E_ep*( 1-(BBtr_pz[track]/p_ep) ); // Obtain Q2 from beam energy, outgoing electron energy, and momenta
      hQ2->Fill( Q2 ); // Fill histogram
	
      double W = PgammaN.M();
      double W_2 = -pow(M_p,2)-2*pow(M_p,2)*(E_e-E_ep)+Q2;

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
	
      //cout << "BBtr_p[track]: " << BBtr_p[track] << endl;
      //cout << "BBtr_px[track]: " << BBtr_px[track] << endl;
      //cout << "BBtr_py[track]: " << BBtr_py[track] << endl;
      //cout << "BBtr_pz[track]: " << BBtr_pz[track] << endl;
      //cout << "BBtr_vz[track]: " << BBtr_vz[track] << endl;
      //cout << "BBtr_chi2[track]: " << BBtr_chi2[track] << endl;
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

      //cout << "runI: " << (int)runI << endl;
      //cout << "runN: " << (int)runN << endl;
      //cout << "runT: " << (int)runT << endl;

      // Coincidence timing cut and vertex cut to resolve W well
      if( fabs(diff+370)<40 && fabs(BBtr_vz[track])<0.08 ) hW_cuts->Fill( W );
	
      if( pow( (HCALx-xexpect_HCAL - dx0)/dx_sig,2) + pow( (HCALy-yexpect_HCAL - dy0)/dy_sig,2) <= pow(2.5,2) ){
	//cout << "Position cut passed.." << endl;
      }
	
      hHCALe->Fill( HCALe );
      hHODOnclus->Fill( nClusHODO );
      if( HCALe>0.02 && HCALtdc[0]<-400 && nClusHODO<10 ) TDCmiss++;
	
      //Get timing offsets for tdc and adc from primary block in cluster
      if( fabs(W-W_mean)<W_sig&&fabs(diff+370)<40 ){

	cout << "W^2:" << pow(W,2) << " W_2^2:" << pow(W_2,2) << " kineW2:" << kineW2 << endl;

	// Get the TOF difference from beam plane normal to HCal
	double vp_el = pp/M_p;
	//cout << vp_el << endl;
	hElas_vp->Fill(vp_el);
	double dd = sqrt( pow(HCal_d,2) + pow(HCALx,2) + pow(HCALy,2) ) - HCal_d;
	double dt = dd/vp_el;
	hTOFcorr->Fill( dt );
	hTOFvID->Fill( cblkid[0], dt );

	double HHdiff = HCALtdc[0]-tHODO; //time difference between HCal TDC and Hodo TDC
	double HHAdiff = HCALa[0]-tHODO; //time difference between HCal ADC time and Hodo TDC
	double e = cblke[0]; //energy deposited in primary block
	int idx = (int)cblkid[0]-1; //index of primary block
	int col = (int)ccol; //column of primary block
	  
	htcorr_HCAL_HODO->Fill( HCALtdc[0], tHODO );
	htDiff_HODO_HCAL->Fill( HHdiff );
	haDiff_HODO_HCAL->Fill( HHAdiff );
	if(ccol>3&&ccol<8){
	  htDiff_HODO_HCAL_JLAB->Fill( HHdiff );
	  haDiff_HODO_HCAL_JLAB->Fill( HHAdiff );
	  htDiff_vs_ADCint_JLAB->Fill( e, HHdiff );
	  haDiff_vs_ADCint_JLAB->Fill( e, HHAdiff );
	}else{
	  htDiff_HODO_HCAL_CMU->Fill( HHdiff );
	  haDiff_HODO_HCAL_CMU->Fill( HHAdiff );
	  htDiff_vs_ADCint_CMU->Fill( e, HHdiff );
	  haDiff_vs_ADCint_CMU->Fill( e, HHAdiff );
	}
	
	//if( cblkid[0]==135 ) cout << HHAdiff << endl;
	hCblkID->Fill( cblkid[0] );
	hCCol->Fill( ccol );
	htDiff_vs_HCALID->Fill( cblkid[0], HHdiff );
	haDiff_vs_HCALID->Fill( cblkid[0], HHAdiff );

	htDiff[ idx ]->Fill( HHdiff );
	haDiff[ idx ]->Fill( HHAdiff );
	if( idx<0 || idx>288 ) cout << "ERROR: TDC out of bounds at tdcID = " << idx << endl;
	if( HHdiff<-65. && HHdiff>-80. ){
	  htDiff_vs_ADCint[idx]->Fill( e, HHdiff );
	  htDiff_vs_ADCint_cols[col]->Fill( e, HHdiff );
	}
	if( HHdiff<-65. && HHdiff>-80. ) htDiff_vs_HCALID_corr->Fill( cblkid[0], HHdiff+22.0*e ); //Timewalk corrected from linear fit to time diff vs clus e
	haDiff_vs_ADCint[idx]->Fill( e, HHAdiff );
	haDiff_vs_ADCint_cols[col]->Fill( e, HHAdiff );
	  
	haDiff_vs_col[col]->Fill( HHAdiff );
	htDiff_vs_col[col]->Fill( HHdiff );
	  
      }
	
      //Fill vertex position histogram for cut on tracks
      hvz_cut->Fill( BBtr_vz[track] );
	
      //Events that pass the above cuts constitute elastics
      elasYield++;
	
      double clusE = 0.0;
	
      hNBlk->Fill( nblk );
	
      // Get energies with simplest scheme from clusters only
      for( int blk = 0; blk<(int)nblk; blk++ ){
	int blkid = int(cblkid[blk])-1; //-1 necessary since sbs.hcal.clus_blk.id ranges from 1 to 288 (different than just about everything else)
	  
	clusE += cblke[blk];
	  
	hADC->Fill( blkid, cblke[blk] );
      }
	
      hDeltaE->Fill( 1.0-(clusE/KE_p) );
      hClusE->Fill( clusE );
      hPAngleCorr->Fill( ephi, p_ep );
      hPAngleCorr_2->Fill( etheta, p_ep );
    }
  }
  
  cout << endl;

  //Fit timing data and extract offsets
  TF1 *f01;
  TF1 *f02;
  TF1 *f03;
  TF1 *f04;
  TF1 *f05;
  TF1 *f06;
  TF1 *f07;
  TF1 *f08;


  double JLAB_mean_TDC = 0.;
  double CMU_mean_TDC = 0.;
  double JLAB_mean_ADCt = 0.;
  double CMU_mean_ADCt = 0.;
  double JLAB_timewalkMean_TDC = 0.;
  double CMU_timewalkMean_TDC = 0.;
  double JLAB_timewalkMean_ADCt = 0.;
  double CMU_timewalkMean_ADCt = 0.;

  if( haDiff_HODO_HCAL_JLAB->GetEntries()>0 ){
    haDiff_HODO_HCAL_JLAB->Fit("gaus","Q","",20,80);
    f01=haDiff_HODO_HCAL_JLAB->GetFunction("gaus");
    JLAB_mean_ADCt = f01->GetParameter(1);
    cout << "JLAB ADCt mean =" << JLAB_mean_ADCt << "." << endl;
  }else{
    cout << "WARNING: Insufficient entries in haDiff_HODO_HCAL_JLAB to fit. Entries = " << haDiff_HODO_HCAL_JLAB->GetEntries() << "." << endl;
  }

  if( haDiff_HODO_HCAL_CMU->GetEntries()>0 ){
    haDiff_HODO_HCAL_CMU->Fit("gaus","Q","",20,80);
    f02=haDiff_HODO_HCAL_CMU->GetFunction("gaus");
    CMU_mean_ADCt = f02->GetParameter(1);
    cout << "CMU ADCt mean =" << CMU_mean_ADCt << "." << endl;
  }else{
    cout << "WARNING: Insufficient entries in haDiff_HODO_HCAL_CMU to fit. Entries = " << haDiff_HODO_HCAL_CMU->GetEntries() << "." << endl;
  }

  if( htDiff_HODO_HCAL_JLAB->GetEntries()>0 ){
    htDiff_HODO_HCAL_JLAB->Fit("gaus","Q","",-150,-65);
    f03=htDiff_HODO_HCAL_JLAB->GetFunction("gaus");
    JLAB_mean_TDC = f03->GetParameter(1);
    cout << "JLAB TDC mean =" << JLAB_mean_TDC << "." << endl;
  }else{
    cout << "WARNING: Insufficient entries in htDiff_HODO_HCAL_JLAB to fit. Entries = " << htDiff_HODO_HCAL_JLAB->GetEntries() << "." << endl;
  }

  if( htDiff_HODO_HCAL_CMU->GetEntries()>0 ){
    htDiff_HODO_HCAL_CMU->Fit("gaus","Q","",-150,-65);
    f04=htDiff_HODO_HCAL_CMU->GetFunction("gaus");
    CMU_mean_TDC = f04->GetParameter(1);
    cout << "CMU TDC mean =" << CMU_mean_TDC << "." << endl;
  }else{
    cout << "WARNING: Insufficient entries in htDiff_HODO_HCAL_CMU to fit. Entries = " << htDiff_HODO_HCAL_CMU->GetEntries() << "." << endl;
  }

  if( htDiff_vs_ADCint_JLAB->GetEntries()>0 ){
    htDiff_vs_ADCint_JLAB->Fit("pol1","Q","",0.0,0.5);
    f05=htDiff_vs_ADCint_JLAB->GetFunction("pol1");
    JLAB_timewalkMean_TDC = f05->GetParameter(1);
    cout << "JLAB TDC Timewalk mean =" << JLAB_timewalkMean_TDC << "." << endl;
  }else{
    cout << "WARNING: Insufficient entries in htDiff_vs_ADCint_JLAB to fit. Entries = " << htDiff_vs_ADCint_JLAB->GetEntries() << "." << endl;
  }

  if( htDiff_vs_ADCint_CMU->GetEntries()>0 ){
    htDiff_vs_ADCint_CMU->Fit("pol1","Q","",0.0,0.5);
    f06=htDiff_vs_ADCint_CMU->GetFunction("pol1");
    CMU_timewalkMean_TDC = f06->GetParameter(1);
    cout << "CMU TDC Timewalk mean =" << CMU_timewalkMean_TDC << "." << endl;
  }else{
    cout << "WARNING: Insufficient entries in htDiff_vs_ADCint_CMU to fit. Entries = " << htDiff_vs_ADCint_CMU->GetEntries() << "." << endl;
  }

  if( haDiff_vs_ADCint_JLAB->GetEntries()>0 ){
    haDiff_vs_ADCint_JLAB->Fit("pol1","Q","",0.0,0.5);
    f07=haDiff_vs_ADCint_JLAB->GetFunction("pol1");
    JLAB_timewalkMean_ADCt = f07->GetParameter(1);
    cout << "JLAB ADCt Timewalk mean =" << JLAB_timewalkMean_ADCt << "." << endl;
  }else{
    cout << "WARNING: Insufficient entries in haDiff_vs_ADCint_JLAB to fit. Entries = " << haDiff_vs_ADCint_JLAB->GetEntries() << "." << endl;
  }

  if( haDiff_vs_ADCint_CMU->GetEntries()>0 ){
    haDiff_vs_ADCint_CMU->Fit("pol1","Q","",0.0,0.5);
    f08=haDiff_vs_ADCint_CMU->GetFunction("pol1");
    CMU_timewalkMean_ADCt = f08->GetParameter(1);
    cout << "CMU ADCt Timewalk mean =" << CMU_timewalkMean_ADCt << "." << endl;
  }else{
    cout << "WARNING: Insufficient entries in haDiff_vs_ADCint_CMU to fit. Entries = " << haDiff_vs_ADCint_CMU->GetEntries() << "." << endl;
  }

  cout << "tFitMin: " << tFitMin << endl;

  for(int i=0; i<kNcell; i++){

    int r = (i)/kNcols;
    int c = (i)%kNcols;

    TF1 *f1;
    TF1 *f2;
    TF1 *f3;
    TF1 *f4;

    //cout << "chan " << i << " htDiff " << htDiff[i]->GetEntries() << endl;

    //Fit HCal-TDC/Hodo-TDC difference histogram by cell and obtain offset parameters
    if( htDiff[i]->GetEntries()>tFitMin && i!=244 && i!=12 && i!=264 ){
      //Tuning 042522
      
      /*
      if( i==5 || i==22 || i==35 || i==47 || i==72 || i==95 || i==120 || i==144 || i==155 || i==167 || i==179 || i==191 || i==251 || i==263 || i==275 || i==287 ){
	htDiff[i]->Fit("gaus","Q","",-70,-50);
      }else if( i==11 || i==203 || i==215 || i==227 || i==228 || i==239 || i==240 ){
	htDiff[i]->Fit("gaus","Q","",-75,-55);
      }else if( i==252 ){
	htDiff[i]->Fit("gaus","Q","",-55,-35);
      }else if( i==276 ){
	htDiff[i]->Fit("gaus","Q","",-45,-30);
      }else{
	htDiff[i]->Fit("gaus","Q","",-80,-60);
      }
      */
      //Tuning with zero offsets
      /*
      if( i==4 || i==5 || i==7 || i==16 || i==17 || i==18 || i==19 || i==28 || i==29 || i==30 || i==31 || i==67 || i==77 || i==78 || i==79 || i==91 || i==102 || i==103 || i==114 || i==115 || i==138 || i==139 || i==147 || i==196 || i==199 || i==212 || i==226 || i==227 || i==250 || i==251 ){
	htDiff[i]->Fit("gaus","Q","",-85,-65);
      }else if( i==6 || i==40 || i==41 || i==42 || i==43 || i==48 || i==49 || i==50 || i==51 || i==56 || i==57 || i==58 || i==59 || i==66 || i==70 || i==71 || i==74 || i==75 || i==88 || i==89 || i==90 || i==92 || i==93 || i==100 || i==101 || i==112 || i==113 || i==124 || i==125 || i==126 || i==127 || i==136 || i==137 || i==144 || i==145 || i==146 || i==152 || i==153 || i==154 || i==155 || i==160 || i==161 || i==162 || i==163 || i==172 || i==173 || i==174 || i==175 || i==180 || i==181 || i==184 || i==185 || i==197 || i==198 || i==206 || i==208 || i==209 || i==222 || i==223 || i==235 || i==246 || i==247 || i==256 || i==257 || i==258 || i==259 || i==268 || i==269 || i==270 || i==271 || i==277 || i==280 || i==281 ){
	htDiff[i]->Fit("gaus","Q","",-80,-60);
      }else if( i==12 ){
	htDiff[i]->Fit("gaus","Q","",-45,-25);
      }else if( i==36 || i==120 || i==156 || i==166 || i==167 || i==214 || i==276 ){
	htDiff[i]->Fit("gaus","Q","",-60,-40);
      }else if( i==37 ){
	htDiff[i]->Fit("gaus","Q","",-70,-40);
      }else if( i==52 || i==53 || i==54 || i==55 || i==64 || i==65 || i==76 || i==186 || i==187 || i==207 || i==213 || i==238 ){
	htDiff[i]->Fit("gaus","Q","",-90,-70);
      }else if( i==148 || i==149 || i==150 || i==151 || i==210 || i==211 || i==220 || i==221 || i==232 || i==233 || i==245 || i==282 || i==283 ){
	htDiff[i]->Fit("gaus","Q","",-95,-80);
      }else if( i==1 || i==9 ){
	htDiff[i]->Fit("gaus","Q","",-66,-56);
      }else if( i==23 ){
	htDiff[i]->Fit("gaus","Q","",-52,-49);
      }else if( i==38 ){
	htDiff[i]->Fit("gaus","Q","",-64,-55);
      }else if( i==60 ){
	htDiff[i]->Fit("gaus","Q","",-120,0);
      }else if( i==73 ){
	htDiff[i]->Fit("gaus","Q","",-68,-62);
      }else if( i==130 ){
	htDiff[i]->Fit("gaus","Q","",-70,-54);
      }else if( i==180 ){
	htDiff[i]->Fit("gaus","Q","",-70,-62);
      }else if( i==204 || i==205 ){
	htDiff[i]->Fit("gaus","Q","",-81,-69);
      }else if( i==228 ){
	htDiff[i]->Fit("gaus","Q","",-60,-49);
      }else if( i==239 ){
	htDiff[i]->Fit("gaus","Q","",-130,-40);
      }else{
	htDiff[i]->Fit("gaus","Q","",-70,-50);
      }
      */

      //cout << "Fitting block: " << i << endl;
      
      //if( i==149 ) continue;

      //Tuning 050322
      /*
      if( i==0 || i==11 ){
	htDiff[i]->Fit("gaus","Q","",-60,-50);
      }else if( i==10 || i==36 || i==60 ){
	htDiff[i]->Fit("gaus","Q","",-60,-40);
      }else if( i==12 ){
	htDiff[i]->Fit("gaus","Q","",-150,0);
      }else if( i==18 || i==30 || i==31 || i==39 || i==61 || i==62 || i==63 || i==66 || i==68 || i==69 || i==74 || i==75 || i==76 || i==77 || i==78 || i==79 || i==80 || i==81 || i==86 || i==87 || i==88 || i==89 || i==90 || i==91 || i==92 || i==93 || i==94 || i==96 || i==97 || i==98 || i==99 || i==100 || i==101 || i==104 || i==105 || i==106 || i==107 || i==108 || i==109 || i==110 || i==111 || i==112 || i==113 || i==117 || i==118 || i==119 || i==121 || i==122 || i==123 || i==124 || i==128 || i==129 || i==134 || i==135 || i==136 || i==137 || i==138 || i==140 || i==141 || i==147 || i==148 || i==149 || i==150 || i==151 || i==152 || i==153 || i==162 || i==163 || i==165 || i==172 || i==173 || i==174 || i==175 || i==176 || i==177 || i==184 || i==185 || i==186 || i==187 || i==188 || i==189 || i==194 || i==195 || i==196 || i==197 || i==198 || i==199 || i==200 || i==201 || i==206 || i==207 || i==208 || i==209 || i==210 || i==211 || i==219 || i==220 || i==222 || i==223 || i==224 || i==225 || i==226 || i==229 || i==230 || i==231 || i==232 || i==233 || i==236 || i==237 || i==238 || i==245 || i==246 || i==247 || i==257 || i==258 || i==259 || i==268 || i==269 || i==270 || i==271 || i==281 || i==282 || i==283 ){
	htDiff[i]->Fit("gaus","Q","",-80,-60);
      }else if( i==24 || i==276 ){
	htDiff[i]->Fit("gaus","Q","",-45,-25);
      }else if( i==102 || i==103 || i==114 || i==115 || i==116 || i==125 || i==126 || i==127 || i==139 || i==212 || i==213 || i==221 || i==234 || i==235 ){
	htDiff[i]->Fit("gaus","Q","",-90,-70);
      }else if( i==19 || i==67 ){
	htDiff[i]->Fit("gaus","Q","",-80,-70);
      }else if( i==35 || i==37 || i==47 || i==155 || i==180 || i==192 ){
	htDiff[i]->Fit("gaus","Q","",-65,-50);
      }else if( i==130 ){
	htDiff[i]->Fit("gaus","Q","",-80,-65);
      }else if( i==240 ){
	htDiff[i]->Fit("gaus","Q","",-65,-40);
      }else if( i==263 ){
	htDiff[i]->Fit("gaus","Q","",-65,-45);
      }else if( i==239 ){
	htDiff[i]->Fit("gaus","Q","",-64,-59);
      }else if( i==95 ){
	htDiff[i]->Fit("gaus","Q","",-61,-50);
      }else{
	htDiff[i]->Fit("gaus","Q","",-80,-50);
      }
      */
      //cout << i << " " << htDiff[i]->GetEntries() << endl;
      if((i>47&&i<60)||i==64||i==65||i==70||i==71||(i>143&&i<156)||i==180||i==181||i==186||i==187){
	htDiff[i]->Fit("gaus","Q","",-85,-64);
      }else if(i==0){
	htDiff[i]->Fit("gaus","Q","",-120,-118);
      }else if(i==24){
	htDiff[i]->Fit("gaus","Q","",-118,-112);
      }else if(i==25){
	htDiff[i]->Fit("gaus","Q","",-130,-115);
      }else if(i==36){
	htDiff[i]->Fit("gaus","Q","",-125,-123);
      }else if(i==96){
	htDiff[i]->Fit("gaus","Q","",-122,-110);
      }else if(i==107){
	htDiff[i]->Fit("gaus","Q","",-125,-105);
      }else if(i==108){
	htDiff[i]->Fit("gaus","Q","",-124,-114);
      }else if(i==131){
	htDiff[i]->Fit("gaus","Q","",-126,-112);
      }else if(i==192){
	htDiff[i]->Fit("gaus","Q","",-125,-105);
      }else if(i==193){
	htDiff[i]->Fit("gaus","Q","",-126,-112);
      }else if(i==203||i==204||i==218||i==228||i==240||i==265||i==267){
	htDiff[i]->Fit("gaus","Q","",-124,-110);
      }else{
	htDiff[i]->Fit("gaus","Q","",-140,-110);
      }
      f1=htDiff[i]->GetFunction("gaus");
      htDiff[i]->SetTitle(Form("TDCGoodFitMean:%f",f1->GetParameter(1)));
      if( f1->GetParameter(1)<-140||f1->GetParameter(1)>-60 ){
	cout << "WARNING: Bad fit in cell " << i << ", defaulting to old calibration value: " << oldTDCoffsets[i] << "." << endl;
	if(c>3&&c<8){
	  TDCoffsets[i] = 0.;
	  htDiff[i]->SetName(Form("htDiff_bl%d_r%d_c%d_JLABm%f",i,r,c,JLAB_mean_TDC));	
	}else{
	  TDCoffsets[i] = 0.;
	  htDiff[i]->SetName(Form("htDiff_bl%d_r%d_c%d_CMUm%f",i,r,c,CMU_mean_TDC));	
	}
      }
      TDCoffsets[i] = -75. - f1->GetParameter(1);
      if(i==0){
	TDCoffsets[i] = -75. - -119.;
      }
      if(i==24){
	TDCoffsets[i] = -75. - -116.;
      }
      if(i==36){
	TDCoffsets[i] = -75. - -124.;
      }
      if(i==48){
	TDCoffsets[i] = -75. - -72.;
      }
      if(i==72){
	TDCoffsets[i] = -75. - -121.;
      }
      if(i==97){
	TDCoffsets[i] = -75. - -116.;
      }
      if(i==155){
	TDCoffsets[i] = -75. - -78.;
      }
      if(i==156){
	TDCoffsets[i] = -75. - -113.;
      }
      if(i==240){
	TDCoffsets[i] = -75. - -116.;
      }
      if(i==242){
	TDCoffsets[i] = -75. - -120.;
      }
      if(i==265){
	TDCoffsets[i] = -75. - -115.;
      }
      if(i==275){
	TDCoffsets[i] = -75. - -114.;
      }
      TDCsig[i] = f1->GetParameter(2);
      double cs = f1->GetChisquare();
      htDiff[i]->SetName(Form("htDiff_bl%d_r%d_c%d_mean%f",i,r,c,f1->GetParameter(1)));	
    }else if( htDiff[i]->GetEntries()>1 ){
      //TDCoffsets[i] = 0.; //Set to zero and default to cosmic-calibrated value
      cout << "WARNING: Low statistics in cell " << i << ", defaulting to old calibration value: " << oldTDCoffsets[i] << "." << endl;
      htDiff[i]->SetName(Form("htDiff_bl%d_r%d_c%d",i,r,c));
    }else{
      cout << "WARNING: Low statistics in cell " << i << ", defaulting to old calibration value: " << oldTDCoffsets[i] << "." << endl;
      if(c>3&&c<8){
	//TDCoffsets[i] = 0.;
	htDiff[i]->SetName(Form("htDiff_bl%d_r%d_c%d_JLABm%f",i,r,c,JLAB_mean_TDC));	
      }else{
	//TDCoffsets[i] = 0.;
	htDiff[i]->SetName(Form("htDiff_bl%d_r%d_c%d_CMUm%f",i,r,c,CMU_mean_TDC));	
      }
    }
    htDiff[i]->SetTitle(Form("TDCDiff_Offset%f",TDCoffsets[i]));

    int j=i+1; //Ham-handed fix to indexing problem and manual fits

    //Fit HCal-ADC-time/Hodo-TDC difference histogram by cell and obtain offset parameters
    if( haDiff[i]->GetEntries()>tFitMin ){
      //Tuning 042522
      /*
      if( j==1 || j==2 || j==3 || j==4 || j==11 || j==12 || j==13 || j==14 || j==23 || j==86 || j==133){
	haDiff[i]->Fit("gaus","Q","",65,80);
      }else if( j==20 || j==26 || j==36 || j==144 || j==157 || j==217 ){
	haDiff[i]->Fit("gaus","Q","",40,55);
      }else if( j==24 ){
	haDiff[i]->Fit("gaus","Q","",66,71);
      }else if( j==37 || j==48 || j==61 || j==62 || j==72 || j==84 || j==97 || j==120 || j==121 || j==132 || j==156 || j==168 || j==180 || j==181 || j==192 || j==193 || j==205 || j==229 || j==266 || j==276 || j==288 ){
	haDiff[i]->Fit("gaus","Q","",45,60);
      }else if( j==85 || j==169 || j==277 ){
	haDiff[i]->Fit("gaus","Q","",60,80);
      }else if( j==96 ){
	haDiff[i]->Fit("gaus","Q","",-25,0);
      }else if( j==109 ){
	haDiff[i]->Fit("gaus","Q","",-80,-65);
      }else if( j==112 ){
	haDiff[i]->Fit("gaus","Q","",-195,-175);
      }else if( j==117 ){
	haDiff[i]->Fit("gaus","Q","",75,90);
      }else if( j==134 ){
	haDiff[i]->Fit("gaus","Q","",-90,-60);
      }else if( j==135 ){ 
	haDiff[i]->Fit("gaus","Q","",-245,-230);
      }else if( j==253 || j==264 ){
	haDiff[i]->Fit("gaus","Q","",0,25);
      }else if( j==278 ){
	haDiff[i]->Fit("gaus","Q","",-125,-100);
      }else if( j==280 ){
	haDiff[i]->Fit("gaus","Q","",-400,-280);
      }else{
	haDiff[i]->Fit("gaus","Q","",0,80);
      }
      */
      haDiff[i]->Fit("gaus","Q","",0,150);
      f2=haDiff[i]->GetFunction("gaus");
      haDiff[i]->SetTitle(Form("ADCtGoodFitMean:%f",f2->GetParameter(1)));
      ADCtoffsets[i] = 50. - f2->GetParameter(1);
      ADCtsig[i] = f2->GetParameter(2);
      //cout << f2->GetParameter(1) << endl;
      double cs = f2->GetChisquare();
      haDiff[i]->SetName(Form("haDiff_bl%d_r%d_c%d_cs%f",i,r,c,cs));
    }else if( htDiff[i]->GetEntries()>10 ){
      ADCtoffsets[i] = 50. - haDiff[i]->GetMean();
      cout << "WARNING: Low statistics in cell " << i << ", defaulting to 50. - arithmetic mean for ADCt offset:" << 50. - haDiff[i]->GetMean() << "." << endl;
      haDiff[i]->SetName(Form("haDiff_bl%d_r%d_c%d",i,r,c));
    }else{
      cout << "WARNING: Low statistics in cell " << i << ", defaulting to 50. - PMT type mean for ADCt offset:" << 50. - haDiff[i]->GetMean() << "." << endl;
      if(c>3&&c<8){
	ADCtoffsets[i] = 50. - JLAB_mean_ADCt;
	haDiff[i]->SetName(Form("haDiff_bl%d_r%d_c%d_JLABm%f",i,r,c,JLAB_mean_ADCt));	
      }else{
        ADCtoffsets[i] = 50. - CMU_mean_ADCt;
	haDiff[i]->SetName(Form("haDiff_bl%d_r%d_c%d_CMUm%f",i,r,c,CMU_mean_ADCt));	
      }
    }
    haDiff[i]->SetTitle(Form("ADCtDiff_Offset%f",ADCtoffsets[i]));

    //Fit HCal-TDC/Hodo-TDC difference vs energy deposited in HCal to obtain TDC timewalk
    if( htDiff_vs_ADCint[i]->GetEntries()>tFitMin ){
      htDiff_vs_ADCint[i]->Fit("pol1","Q","",0.0,0.2);
      f3=htDiff_vs_ADCint[i]->GetFunction("pol1");
      TvsE[i]=f3->GetParameter(1);
      hTDCTimewalk->Fill(f3->GetParameter(1));
    }else{
      cout << "Low/absent statistics requires using mean TDCvsADCint slope by PMT type: ";
      if(ccol>3&&ccol<8){
	TvsE[i] = JLAB_timewalkMean_TDC;
      }else{
	TvsE[i] = CMU_timewalkMean_TDC;
      }
      cout << TvsE[i] << "." << endl;
    }

    //Fit HCal-ADC-time/Hodo-TDC difference vs energy deposited in HCal to obtain ADCt timewalk
    if( haDiff_vs_ADCint[i]->GetEntries()>tFitMin ){
      haDiff_vs_ADCint[i]->Fit("pol1","Q","",0.0,0.2);
      f4=haDiff_vs_ADCint[i]->GetFunction("pol1");
      AtvsE[i]=f4->GetParameter(1);
      hADCtTimewalk->Fill(f4->GetParameter(1));
    }else{
      cout << "Low/absent statistics requires using mean ADCtvsADCint slope by PMT type: ";
      if(ccol>3&&ccol<8){
	AtvsE[i] = JLAB_timewalkMean_ADCt;
      }else{
	AtvsE[i] = CMU_timewalkMean_ADCt;
      }
      cout << AtvsE[i] << "." << endl;
    }
  }

  //Make offsets input-ready and improve readability
  for(int i=0; i<kNcell; i++){
    hTDCoffsets_vs_HCALID->Fill(i,TDCoffsets[i]);
    holdTDCoffsets_vs_HCALID->Fill(i,oldTDCoffsets[i]*TDCCalib);  
    combTDCoffsets[i] = -TDCoffsets[i]/TDCCalib+oldTDCoffsets[i]; //TODO: Must double check this with new data
    hcombTDCoffsets_vs_HCALID->Fill(i,combTDCoffsets[i]);     
    hTDCsig_vs_HCALID->Fill(i,TDCsig[i]);
    
    hADCtoffsets_vs_HCALID->Fill(i,ADCtoffsets[i]);
    combADCtoffsets[i] = ADCtoffsets[i]+oldADCtoffsets[i];
    hADCtsig_vs_HCALID->Fill(i,ADCtsig[i]);
  }
 
  //Make canvas to hold all fits for comparison to HCal geometry
  TCanvas *TDC_top = new TCanvas("TDC_top","TDC_top",1600,1200);
  TCanvas *TDC_bot = new TCanvas("TDC_bot","TDC_bot",1600,1200);
  TCanvas *ADCt_top = new TCanvas("ADCt_top","ADCt_top",1600,1200);
  TCanvas *ADCt_bot = new TCanvas("ADCt_bot","ADCt_bot",1600,1200);

  TDC_top->Divide(12,12);
  TDC_bot->Divide(12,12);
  ADCt_top->Divide(12,12);
  ADCt_bot->Divide(12,12);

  gStyle->SetOptStat(0);
  for(int i=0; i<kNcell; i++){
    TDC_top->cd(i+1);
    if( i>=144 ){
      TDC_bot->cd(i-143);
      gStyle->SetOptStat(0);
    }
    if(htDiff[i]->GetEntries()<tFitMin){
      htDiff[i]->SetAxisColor(2);
    }else{
      htDiff[i]->SetAxisColor(1);
    }
    htDiff[i]->Draw();
  }

  gStyle->SetOptStat(0);
  for(int i=0; i<kNcell; i++){
    ADCt_top->cd(i+1);
    if( i>=144 ){
      ADCt_bot->cd(i-143);
      gStyle->SetOptStat(0);
    }
    if(haDiff[i]->GetEntries()<tFitMin){
      haDiff[i]->SetAxisColor(2);
    }else{
      haDiff[i]->SetAxisColor(1);
    }
    haDiff[i]->Draw();
  }

  //Write out diagnostic histos and print to console
  fout->Write();

  //Declare outfile
  ofstream params;
  params.open( paramsPath );
  params << "#Timewalk parameters from SBS-" << kine << " obtained " << date.c_str() << endl;
  params << "#" << endl << "#" << endl << "#HCal_TDC_timewalk = " << endl;
  cout << endl << "TDC Timewalk Corrections: " << endl << endl;

  int cell = 0;
  for( int r=0; r<kNrows; r++ ){
    for( int c=0; c<kNcols; c++ ){
      params << TvsE[cell] << "  ";
      cout << TvsE[cell] << "  ";
      cell++;
    }
    params << endl;
    cout << endl;
  }

  params << "#" << endl << "#" << endl << "#HCal_ADCt_timewalk = " << endl;
  cout << endl << endl << "ADCt Timewalk Corrections: " << endl << endl;

  cell = 0;
  for( int r=0; r<kNrows; r++ ){
    for( int c=0; c<kNcols; c++ ){
      params << AtvsE[cell] << "  ";
      cout << AtvsE[cell] << "  ";
      cell++;
    }
    params << endl;
    cout << endl;
  }

  //params << "#" << endl << "#" << endl << "#HCal_TDC_offset = " << endl;
  cout << endl << endl << "Base TDC Offsets: " << endl << endl;
  
  cell = 0;
  vector<int> TDC_noDataCell;
  for( int r=0; r<kNrows; r++ ){
    for( int c=0; c<kNcols; c++ ){
      if( htDiff[cell]->GetEntries()<1 ) TDC_noDataCell.push_back(cell);
      
      //params << TDCoffsets[cell] << " ";
      cout << TDCoffsets[cell] << " ";
      
      cell++;
    }
    //params << endl;
    cout << endl;
  }

  params << "#" << endl << "#" << endl << "#HCal_TDC_offset = " << endl;
  cout << endl << endl << "TDC Offsets: " << endl << endl;
  
  cell = 0;
  //vector<int> TDC_noDataCell;
  for( int r=0; r<kNrows; r++ ){
    for( int c=0; c<kNcols; c++ ){
      //if( htDiff[cell]->GetEntries()<1 ) TDC_noDataCell.push_back(cell);
      
      params << combTDCoffsets[cell] << " ";
      cout << combTDCoffsets[cell] << " ";
      
      cell++;
    }
    params << endl;
    cout << endl;
  }

  if( TDC_noDataCell.size()>0 ){
    cout << endl << endl << "WARNING: Missing TDC Data in Cells (idx from 0): ";
    for( int i=0; i<TDC_noDataCell.size(); i++ ){
      cout << TDC_noDataCell[i] << " ";
    }
    cout << endl;
  }

  params << "#" << endl << "#" << endl << "#HCal_ADCtime_offset = " << endl;
  cout << endl << endl << "ADC Time Offsets: " << endl << endl;

  cell = 0;
  vector<int> ADCt_noDataCell;
  for( int r=0; r<kNrows; r++ ){
    for( int c=0; c<kNcols; c++ ){
      if( haDiff[cell]->GetEntries()<1 ) ADCt_noDataCell.push_back(cell); 

      params << combADCtoffsets[cell] << " ";
      cout << combADCtoffsets[cell] << " ";
	
      cell++;
    }
    params << endl;
    cout << endl;
  }

  if( ADCt_noDataCell.size()>0 ){
    cout << endl << endl << "WARNING: Missing ADCt Data in Cells (idx from 0): ";
    for( int i=0; i<ADCt_noDataCell.size(); i++ ){
      cout << ADCt_noDataCell[i] << " ";
    }
    cout << endl;
  }

  params.close();

  if( diag==1 ){
    TString plotsfilename = outputfilename;
    plotsfilename.ReplaceAll(".root","_1.pdf");
    TDC_top->Print(plotsfilename.Data(),"pdf");
    plotsfilename.ReplaceAll(".root","_2.pdf");
    TDC_bot->Print(plotsfilename.Data(),"pdf");
    plotsfilename.ReplaceAll(".root","_3.pdf");
    ADCt_top->Print(plotsfilename.Data(),"pdf");
    plotsfilename.ReplaceAll(".root","_4.pdf");
    ADCt_bot->Print(plotsfilename.Data(),"pdf");
  }

  cout << endl << endl << "Elastic yield for analyzed runs: " << elasYield << ". Total events analyzed: " << Nevents << ". Total TDC misses: " << TDCmiss << "." << endl << endl;

  cout << "Timing offset analysis complete and parameters written to file." << endl;

  elist->Delete();
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
//TDC offsets 050222
/*
51.6052 53.1381 75.711 62.3824 -64.8101 -109.37 -83.1522 -84.9016 20.6541 15.1315 74.9802 46.4545 
78.6722 78.2396 63.5186 78.6833 -88.5921 -82.0225 -83.3094 -49.5275 85.6185 15.701 2.15878 56.5909 
64.0445 52.6813 87.4801 78.7829 -94.516 -81.0376 -86.6389 -70.9349 64.2656 36.1346 8.05743 -0.205335 
95.8128 98.0925 122.216 65.0213 -31.9919 -45.628 -41.1135 -47.1724 63.1208 84.5394 58.4076 22.4098 
-71.5665 -9.2373 -8.41619 -42.0759 -108.828 -129.119 -137.452 -140.078 -33.0511 -36.4883 -40.984 -52.639 
-12.3159 113.306 98.7955 107.096 -105.371 -109.768 -14.4186 7.78215 144.632 134.126 4.81561 -26.0724 
-3.63963 77.3066 64.1368 85.0212 -32.1225 -45.6403 -39.5834 -69.7827 103.482 116.736 80.5107 27.1408 
96.7496 88.5692 99.0824 76.6522 -28.3414 -38.6806 -48.1009 -37.0115 83.5445 102.798 66.9183 20.0184 
111.863 102.581 110.541 99.2195 -15.5739 -5.38606 5.35004 36.44 122.533 112.388 111.026 104.949 
56.5101 95.6989 72.8023 88.6887 6.61235 13.7108 21.1924 32.6727 133.175 94.1777 108.589 111.913 
94.3001 108.904 78.226 78.4414 20.1627 13.6701 40.2308 55.2586 132.028 116.528 112.243 99.1502 
75.4167 75.8384 95.8968 94.4299 -10.852 -30.9235 -13.2526 14.2203 136.567 122.191 94.9797 59.2745 
-80.4373 -12.4915 -81.7404 -51.1083 -133.647 -117.525 -126.039 -128.933 -9.862 -28.549 -30.9201 -55.9466 
128.323 125.449 83.2571 51.8873 8.25293 -11.8253 23.2971 31.9416 109.228 75.6646 109.545 64.6167 
97.4597 121.844 74.3539 65.1699 -13.5127 2.91005 14.676 33.2594 109.414 54.919 87.1137 73.8852 
1.23497 -10.884 68.6788 86.5581 -7.34149 -9.09077 -111.255 -111.246 108.57 115.432 85.3385 66.5992 
95.3282 92.8283 131.052 116.433 13.6149 40.6905 -8.12992 -43.928 121.406 122.156 124.906 49.996 
-109.833 -72.157 -45.4263 -44.8109 50.9042 34.6486 -144.167 -155.273 -9.62435 -9.40513 153.829 99.2835 
55.4669 93.79 100.023 76.4594 -138.89 -98.5467 29.1635 7.32819 113.922 93.196 -56.2945 -73.2276 
78.7678 141.75 112.669 136.596 -142.095 -138.187 48.2404 37.4863 154.622 154.294 -55.3294 -84.8173 
79.6345 105.374 99.1642 115.317 -141.702 -104.695 22.2861 11.8522 92.3068 67.9966 -88.9921 -81.4167 
21.8488 85.6216 88.1268 61.6097 -30.4453 0.119774 24.8667 2.4181 93.3816 52.5699 52.9011 25.551 
103.489 135.58 122.72 57.343 -20.2896 19.593 -0.512462 17.5389 110.286 55.6388 42.3875 52.3349 
-56.248 -51.1446 98.0804 87.691 -32.4454 0.858946 -156.045 -166.215 90.6417 43.3894 29.6404 44.151
*/
