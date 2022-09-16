//SSeeds 9.1.22 - Post-production - Pass1 - Calibration code which employs basic cuts on elastic events to obtain timing offsets for HCal ADC time and TDC. Leaving out proton spot cuts at this stage to obtain better statistics (will improve with these cuts after first pass mass replay).

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

// Dying exponential fit for timewalk corrections
double TW_fit(double *x, double *par){
  double amp = par[0];
  double str = par[1];
  double offset = par[2];
  return amp*exp(str*x[0])+offset;
}

void tcal( const char *configfilename="setup_timingCal.cfg", int run = -1 ){
  
  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  // Get the date
  string date = getDate();
  
  // Declare Chain for many root files
  TChain *C = new TChain("T");

  // Path for output file
  string paramsPath = "/w/halla-scshelf2102/sbs/seeds/HCal_replay/hcal/Hodo_HCal_Coorelations/timingCalParams_v2.txt";

  // Declare general physics parameters to be modified by input config file
  double test = 1; // Keep track of DB read for analysis. 0: ../SBS_REPLAY/SBS_Replay/DB   1: ../seeds/SBS-replay/DB
  double diag = 0; // Keep track of diagnostic canvas pdf option printed at end
  double kine = 8; // Keep track of kinematic calibrating from
  double spotcut = 0;
  double tFitMin = 30; // Minimum number of entries per channel to calibrate ADC/TDC time
  double t_trig = 510; // Mean tdc trig value (HCAL - BB) 
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
	test = sval.Atof();
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
  cout << endl << endl << "Loading previous offsets from file: " << inConstPath << ".." << endl;
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

  cout << "Tree variables linked." << endl;

  // Define a clock to check macro processing time
  TStopwatch *sw = new TStopwatch();
  sw->Start();
  
  // Declare outfile
  TString outputfilename = Form("tcalout%d_%s.root", (int)kine, date.c_str());
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
  TH1D *hTOF = new TH1D( "hTOF", "Time of Flight Corrections; ns", 200, 0, 10 );
  TH1D *hTOF2 = new TH1D( "hTOF2", "Time of Flight Corrections (nu method); ns", 200, 0, 10 );
  TH1D *hTOFvp = new TH1D( "hTOFvp", "Proton Velocity", 400, 0, 4 );
  TH1D *hTOFvp2 = new TH1D( "hTOFvp2", "Proton Velocity (nu method)", 400, 0, 4 );
  TH1D *hPP = new TH1D( "hPP", "Proton Momentum; GeV", 200, 0, 10 );
  TH1D *hPdp = new TH1D( "hPdp", "Proton Path Distance; m", 200, 0, 20 );
  TH1D *hdx = new TH1D( "hdx", "HCal X - X Expected", 400, -2, 2 );
  TH1D *hdy = new TH1D( "hdy", "HCal Y - Y Expected", 400, -2, 2 );
  TH1D *hCblkID = new TH1D( "hCblkID", "Block ID", 300, -10, 290 );
  TH1D *hCCol = new TH1D( "hCCol", "Cluster Seed Column", 15, -1, 14 );
  TH1D *hHCALe = new TH1D( "hHCALe","E HCal Cluster E", 400, 0., 4 );
  TH1D *hHODOnclus = new TH1D( "hHODOnclus","Number of Hodoscope Clusters", 50, 0., 50. );
  TH1D *hTBBt = new TH1D( "hTBBt","BBCal Trigger (L1A) (ns)", 500, 340, 390 );
  TH1D *hTRF = new TH1D( "hTRF","RF Signature (ns)", 1700, -10, 160 );
  TH1D *hTRFmod = new TH1D( "hTRFmod","RF Signature Modulo 4ns (ns)", 200, -10, 10 );
  TH1D *hTEDTM = new TH1D( "hTEDTM","Electronic Dead Time Monitor Signal (ns)", 2000, -1, 2 );
  TH1D *hTBBlo = new TH1D( "hTBBlo","BBCal Lo Trigger (ns)", 199, 1, 200 );
  TH1D *hTBBhiv = new TH1D( "hTBBhiv","BBCal Lo (ns)", 2010, -10, 2000 );
  TH1D *hTHCAL = new TH1D( "hTHCAL","HCal time (ns)", 1000, 800, 900 ); //trigger
  TH1D *hTHODO = new TH1D( "hTHODO","Hodo mean TDC time (ns)", 2000, -1000, 1000 );
  TH1D *hTHCALvRF = new TH1D( "hTHCALvRF","HCal time - RF Signature Modulo 4ns (ns)", 1000, 800, 900 );
  TH1D *hDiff = new TH1D( "hDiff","HCal time - BBCal time (ns)", 1300, -500, 800 );
  TH1D *hW = new TH1D( "W", "W; GeV", 250, 0.3, 1.5 );
  TH1D *hNBlk = new TH1D( "hNBlk", "Number of Blocks in Primary Cluster; Number", 25, 0, 25 );
  TH1D *hW_cuts = new TH1D( "W_cuts", "W_cuts; GeV", 250, 0.3, 1.5 );
  TH1D *hQ2 = new TH1D( "Q2", "Q2; GeV", 250, 0.5, 3.0 );
  TH1D *hE_ep = new TH1D( "Scattered Electron Energy","E_ep; GeV", 500, 0.0, E_e*1.5 ); 
  TH1D *hE_pp = new TH1D( "Scattered Proton Energy", "E_pp; GeV", 500, 0.0, E_e*1.5 );
  TH1D *hKE_p = new TH1D( "Scattered Proton Kinetic Energy", "KE_pp; GeV", 500, 0.0, E_e*1.5 );
  TH2D *hdxdy_HCAL = new TH2D("hdxdy_HCAL",";y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)", 250, -5.0, 5.0, 250, -10, 10 );
  TH2D *hXY_HCAL = new TH2D("hXY_HCAL",";y_{HCAL} (m); x_{HCAL} (m)", 250, -5.0, 5.0, 250, -10, 10 );
  TH1D *hvz_cut = new TH1D("hvz_cut",";vertex z (m);", 250,-0.125,0.125);
 
  // Timing histograms
  TH2D *htcorr_HCAL_HODO = new TH2D("HCAL/HODO TDC/TDC Correlation",";TDC_{HCAL} (ns);TDC_{HODO} (ns)", 150, -150, 0, 150, -15, 15);  
  TH2D *hacorr_HCAL_HODO = new TH2D("HCAL/HODO ADCt/TDC Correlation",";ADCt_{HCAL} (ns);TDC_{HODO} (ns)", 300, 0, 150, 150, -15, 15);
  TH2D *htcorr_HCAL_HODO_corr = new TH2D("HCAL/HODO TDC/TDC Correlation (corrected)",";TDC_{HCAL} (ns);TDC_{HODO} (ns)", 150, -150, 0, 150, -15, 15);  
  TH2D *hacorr_HCAL_HODO_corr = new TH2D("HCAL/HODO ADCt/TDC Correlation (corrected)",";ADCt_{HCAL} (ns);TDC_{HODO} (ns)", 300, 0, 150, 150, -15, 15);
  TH1D *htHCAL = new TH1D("HCal_TDC","HCal TDC (All Channels);ns",300,-150,0);
  TH1D *htDiff_HODO_HCAL = new TH1D("HCalHODO_TDC","HCal TDC (All Channels);ns",300,-150,0);
  TH1D *haDiff_HODO_HCAL = new TH1D("HCal_ADCt-HODO_TDC","HCal ADCt (All Channels);ns",300,0,150);
  TH1D *htDiff_HODO_HCAL_corr = new TH1D("HCal_TDC-HODO_TDC Corrected","ns",300,-150,0);
  TH1D *haDiff_HODO_HCAL_corr = new TH1D("HCal_ADCt-HODO_TDC Corrected","ns",300,0,150);
  TH1D *haDiff_HODO_HCAL_JLAB = new TH1D("HCal_ADCt-HODO_TDC JLAB","ns",300,0,150);
  TH1D *haDiff_HODO_HCAL_CMU = new TH1D("HCal_ADCt-HODO_TDC CMU","ns",300,0,150);
  TH1D *htDiff_HODO_HCAL_JLAB = new TH1D("HCal_TDC-HODO_TDC JLAB","ns",300,-150,0);
  TH1D *htDiff_HODO_HCAL_CMU = new TH1D("HCal_TDC-HODO_TDC CMU","ns",300,-150,0);
  TH1D *haDiff_HODO_HCAL_JLAB_corr = new TH1D("HCal_ADCt-HODO_TDC JLAB Corrected","ns",300,0,150);
  TH1D *haDiff_HODO_HCAL_CMU_corr = new TH1D("HCal_ADCt-HODO_TDC CMU Corrected","ns",300,0,150);
  TH1D *htDiff_HODO_HCAL_JLAB_corr = new TH1D("HCal_TDC-HODO_TDC JLAB Corrected","ns",300,-150,0);
  TH1D *htDiff_HODO_HCAL_CMU_corr = new TH1D("HCal_TDC-HODO_TDC CMU Corrected","ns",300,-150,0);
  TH2D *htDiff_vs_ADCint_JLAB = new TH2D("htDiff_vs_ADCint_JLAB",";E_{JLAB} (GeV);TDC_{HCAL}-TDC_{HODO} (ns)",270,0.0,0.9,600,-250,50);
  TH2D *haDiff_vs_ADCint_JLAB = new TH2D("haDiff_vs_ADCint_JLAB",";E_{JLAB} (GeV);ADCt_{HCAL}-TDC_{HODO} (ns)",270,0.0,0.9,600,0,150);
  TH2D *htDiff_vs_ADCint_CMU = new TH2D("htDiff_vs_ADCint_CMU",";E_{CMU} (GeV);TDC_{HCAL}-TDC_{HODO} (ns)",270,0.0,0.9,600,-250,50);
  TH2D *haDiff_vs_ADCint_CMU = new TH2D("haDiff_vs_ADCint_CMU",";E_{CMU} (GeV);ADCt_{HCAL}-TDC_{HODO} (ns)",270,0.0,0.9,600,0,150);
  TH1D *hTDCTimewalk = new TH1D("hTDCTimewalk","ns",500,0,50);
  TH1D *hADCtTimewalk = new TH1D("hADCtTimewalk","ns",500,0,50);
  TH2D *hHCALtdc_vs_HODOtdc = new TH2D("hHCALtdc_vs_HODOtdc",";TDC_{HODO} (ns);TDC_{HCAL} (ns)",150,-75,75,150,-150,0);
  TH2D *htDiff_vs_HCALID = new TH2D("htDiff_vs_HCALID",";Channel;TDC_{HCAL}-TDC_{HODO} (ns)",288,0,288,300,-150,0);
  TH2D *htHCAL_vs_HCALID = new TH2D("htHCAL_vs_HCALID",";Channel;TDC_{HCAL} (ns)",288,0,288,300,-150,0);
  TH2D *htHODO_vs_HCALID = new TH2D("htHODO_vs_HCALID",";Channel;TDC_{HODO} (ns)",288,0,288,300,-75,75);
  TH2D *htDiff_bb_vs_HCALID = new TH2D("htDiff_bb_vs_HCALID",";Channel;TDC_{HCAL}-BB_{trigtime} (ns)",288,0,288,400,-500,-300);
  TH2D *htDiff_bbrf_vs_HCALID = new TH2D("htDiff_bbrf_vs_HCALID",";Channel;TDC_{HCAL}-BB_{trigtime}-RFtime_mod4 (ns)",288,0,288,400,-500,-300);
  TH2D *htDiffBB_vs_HCALID = new TH2D("htDiffBB_vs_HCALID",";Channel;TDC_{HCAL}-BBCalt (ns)",288,0,288,1000,-500,500);
  TH2D *htDiffRF_vs_HCALID = new TH2D("htDiffRF_vs_HCALID",";Channel;(TDC_{HCAL}-RF)%4 (ns)",288,0,288,300,-150,0);
  TH2D *haDiff_vs_HCALID = new TH2D("haDiff_vs_HCALID",";Channel;ADCtime_{HCAL}-TDC_{HODO} (ns)",288,0,288,300,0,150);
  TH2D *haHCAL_vs_HCALID = new TH2D("haHCAL_vs_HCALID",";Channel;ADCtime_{HCAL} (ns)",288,0,288,300,0,150);
  TH2D *hTDCoffsets_vs_HCALID = new TH2D("hTDCoffsets_vs_HCALID",";Channel;TDC Offset (ns)",288,0,288,200,-50,50);
  TH2D *holdTDCoffsets_vs_HCALID = new TH2D("holdTDCoffsets_vs_HCALID",";Channel;TDC Offset (ns)",288,0,288,200,-50,50);
  TH2D *hcombTDCoffsets_vs_HCALID = new TH2D("hcombTDCoffsets_vs_HCALID",";Channel;TDC Offset (ns)",288,0,288,200,-50,50);
  TH2D *hADCtoffsets_vs_HCALID = new TH2D("hADCtoffsets_vs_HCALID",";Channel;ADC Time Offset (ns)",288,0,288,200,-50,50);
  TH2D *hTDCsig_vs_HCALID = new TH2D("hTDCsig_vs_HCALID",";Channel;TDC Std Dev (ns)",288,0,288,100,0,10);
  TH2D *hADCtsig_vs_HCALID = new TH2D("hADCtsig_vs_HCALID",";Channel;ADC Time Std Dev (ns)",288,0,288,500,0,50);
  TH2D *htDiff_vs_HCALID_corr = new TH2D("htDiff_vs_HCALID_corr",";Channel;TDC_{HCAL}-TDC_{HODO} (ns)",288,0,288,300,-150,0);
  TH2D *haDiff_vs_HCALID_corr = new TH2D("haDiff_vs_HCALID_corr",";Channel;ADCt_{HCAL}-TDC_{HODO} (ns)",288,0,288,300,0.0,150.0);
  TH2D *htDiff_vs_ADCint[kNcell];
  TH2D *haDiff_vs_ADCint[kNcell];
  TH1D *htDiff[kNcell];
  TH1D *htDiff_bb[kNcell];
  TH1D *htDiff_bbrf[kNcell];
  TH1D *haDiff[kNcell];
  TH1D *htDiff_corr[kNcell];
  TH1D *haDiff_corr[kNcell];
  TH1D *haDiff_vs_col[kNcols];
  TH1D *htDiff_vs_col[kNcols];
  TH1D *haDiff_vs_col_corr[kNcols];
  TH1D *htDiff_vs_col_corr[kNcols];
  TH2D *htDiff_vs_ADCint_cols[kNcols];
  TH2D *haDiff_vs_ADCint_cols[kNcols];

  cout << endl;

  for( int i=0; i<kNcols; i++){
    haDiff_vs_col[i] = new TH1D(Form("haDiff_col%d",i),";ADCt_{HCAL}-TDC_{HODO} (ns)",300,0,150);
    htDiff_vs_col[i] = new TH1D(Form("htDiff_col%d",i),";TDC_{HCAL}-TDC_{HODO} (ns)",300,-150,0);
    haDiff_vs_col_corr[i] = new TH1D(Form("haDiff_col%d_corr",i),";ADCt_{HCAL}-TDC_{HODO} (ns)",300,0,150);
    htDiff_vs_col_corr[i] = new TH1D(Form("htDiff_col%d_corr",i),";TDC_{HCAL}-TDC_{HODO} (ns)",300,-150,0);
    haDiff_vs_ADCint_cols[i] = new TH2D(Form("htDiff_vs_ADCint_col%d",i),Form(";E_{col%d} (GeV);TDC_{HCAL}-TDC_{HODO} (ns)",i),270,0.0,0.9,600,-250,50);
    htDiff_vs_ADCint_cols[i] = new TH2D(Form("haDiff_vs_ADCint_col%d",i),Form(";E_{col%d} (GeV);ADCt_{HCAL}-TDC_{HODO} (ns)",i),270,0.0,0.9,600,0,150);
  }

  for( int i=0; i<kNcell; i++ ){
    htDiff_vs_ADCint[i] = new TH2D(Form("htDiff_vs_ADCint_bl%d",i),Form(";E_{bl%d} (GeV);TDC_{HCAL}-TDC_{HODO} (ns)",i),270,0.0,0.9,600,-250,50);
    htDiff_vs_ADCint[i]->Fill( 0.1, -150 );
    haDiff_vs_ADCint[i] = new TH2D(Form("haDiff_vs_ADCint_bl%d",i),Form(";E_{bl%d} (GeV);ADCt_{HCAL}-TDC_{HODO} (ns)",i),270,0.0,0.9,600,0,150);
    haDiff_vs_ADCint[i]->Fill( 0.1, 55 ); //To test fits

    htDiff[i] = new TH1D(Form("htDiff_bl%d",i),";TDC_{HCAL}-TDC_{HODO} (ns)",300,-150,0);
    htDiff_bb[i] = new TH1D(Form("htDiff_bb_bl%d",i),";TDC_{HCAL}-Trig_{BB} (ns)",400,-500,-300);
    htDiff_bbrf[i] = new TH1D(Form("htDiff_bbrf_bl%d",i),";TDC_{HCAL}-Trig_{BB,rf corr} (ns)",400,-500,-300);
    haDiff[i] = new TH1D(Form("haDiff_bl%d",i),";ADCt_{HCAL}-TDC_{HODO} (ns)",1400,-400,300);
    htDiff_corr[i] = new TH1D(Form("htDiff_corr_bl%d",i),";TDC_{HCAL}-TDC_{HODO} Corrected (ns)",300,-150,0);
    haDiff_corr[i] = new TH1D(Form("haDiff_corr_bl%d",i),";ADCt_{HCAL}-TDC_{HODO} Corrected (ns)",1400,-400,300);
    haDiff[i]->Fill( 55 ); //To test fits
    htDiff_corr[i]->Fill( -75. );
  }

  cout << "Variables and histograms defined." << endl;

  // Set long int to keep track of total entries
  Long64_t Nevents = elist->GetN();
  UInt_t run_number = 0;

  cout << endl << "All parameters loaded and initialization complete." << endl << endl;
  cout << "Opened up tree with nentries: " << C->GetEntries() << ", nentries passing globalcut: " << Nevents << "." << endl << endl;

  //First loop over events to obtain mean TOF and timewalk params by channel
  cout << "First loop over all data commencing.." << endl;
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

      //Plot 2D histo of position
      hXY_HCAL->Fill( HCALy, HCALx );

      //Define intersection points for hadron vector
      double sintersect = ( HCAL_origin - vertex ).Dot( HCAL_zaxis ) / ( pNhat.Dot( HCAL_zaxis ) );

      TVector3 HCAL_intersect = vertex + sintersect * pNhat;

      double yexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_yaxis );
      double xexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_xaxis );

      double E_ep = sqrt( pow(M_e,2) + pow(BBtr_p[0],2) ); // Obtain the scattered electron energy
      hE_ep->Fill( E_ep ); // Fill histogram
	
      double p_ep = BBtr_p[0];
      double Q2 = 2*E_e*E_ep*( 1-(BBtr_pz[0]/p_ep) ); // Obtain Q2 from beam energy, outgoing electron energy, and momenta
      hQ2->Fill( Q2 ); // Fill histogram
	
      double W = PgammaN.M();
      hW->Fill( W );
	
      //Use the electron kinematics to predict the proton momentum assuming elastic scattering on free proton at rest:
      double E_pp = nu+M_p; // Get total energy of the proton
      //double Enucleon = sqrt(pow(pp,2)+pow(M_p,2)); // Check on E_pp, same
      hE_pp->Fill( E_pp ); // Fill histogram
	
      double KE_p = nu; //For elastics
      hKE_p->Fill( KE_p );   
	
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
      hTBBt->Fill( bbcal_time );
      hTRF->Fill( RF_time );
      hTRF->Fill( RFmod );
      hTEDTM->Fill( edtm_time );
      hTBBlo->Fill( bbcalLO_time );
      hTBBhiv->Fill( bbcalHIveto_time );
      hTHCAL->Fill( hcal_time );
      hTHCALvRF->Fill( hcal_time - RFmod );

      hDiff->Fill( diff );

      //Fill some HCal histos
      hHCALe->Fill( HCALe );
      hHODOnclus->Fill( nClusHODO );

      double xDiff = HCALx - xexpect_HCAL;
      double yDiff = HCALy - yexpect_HCAL;

      //Check how often the tdc failed to register for an otherwise good event
      if( HCALe>0.02 && HCALtdc[0]<-400 && nClusHODO<10 ) TDCmiss++;
	
      bool nogo = spotcut==1 && (fabs(xDiff-dx0)>dx_sig || fabs(yDiff-dy0)>dy_sig);

      //Get timing offsets for tdc and adc from primary block in cluster
      if( fabs( W-W_mean )<W_sig&&fabs( diff+t_trig )<40&&nogo==false ){

	//Fill some basic histos
	hW_cuts->Fill( W );
	hCblkID->Fill( cblkid[0] );
	hCCol->Fill( ccol );
	hvz_cut->Fill( BBtr_vz[0] );
	//Draw some location plots for better cuts on elastics later on
	hdxdy_HCAL->Fill( yDiff, xDiff );
	hdy->Fill( yDiff );
	hdx->Fill( xDiff );
	//cout << "RF_time:" << RF_time << ", RFmod:" << RFmod << ", hcal_time:" << hcal_time << endl;

	//Get total time of flight per elastic proton (very naive - need to use beta)
	//Recall: p=(m*v)/sqrt(1-pow( beta, 2 )), beta = v/c, c=1
	double vp = pp/M_p; //Assuming elastic
	hPP->Fill( pp );
	double vp2 = sqrt( 1/pow( M_p/pp, 2 )); //Via beta for comparison (and measure of inelastic contamination)?
	hTOFvp->Fill( vp );
	hTOFvp2->Fill( vp2 );

	double dp = sqrt( pow(HCal_d,2) + pow(HCALx,2) + pow(HCALy,2) ); //Distance from vertex to center of cluster
	hPdp->Fill( dp );
	double tp = dp/vp;
	double tp2 = dp/vp2;
	hTOF->Fill( tp );
	hTOF2->Fill( tp2 );;

	//Get time difference between ADCt/TDC and HODO TDC mean time
	double HHdiff = HCALtdc[0]-tHODO; //time difference between HCal TDC and Hodo TDC
	double HBBdiff = HCALtdc[0]-bbcal_time; //time difference between HCal TDC and BBCal trigger
	double HBBdiff_rf = HCALtdc[0]-bbcal_time-RFmod; //time difference between HCal TDC and BBCal trigger
	double HRFdiff = HCALtdc[0]-RFmod; //time difference between HCal TDC and RF time modulo 4ns
	double HHAdiff = HCALa[0]-tHODO; //time difference between HCal ADC time and Hodo TDC
	double e = cblke[0]; //energy deposited in primary block
	int idx = (int)cblkid[0]-1; //index of primary block
	int col = (int)ccol; //column of primary block
	if( idx<0 || idx>=288 ) cout << "ERROR: indexing out of bounds at idx = " << idx << endl;

	//Fill some unmodified time difference histos
	hTHODO->Fill( tHODO );
	hHCALtdc_vs_HODOtdc->Fill( tHODO, HCALtdc[0] );
	htHCAL->Fill( HCALtdc[0] );
	htHCAL_vs_HCALID->Fill( cblkid[0], HCALtdc[0] );
	htHODO_vs_HCALID->Fill( cblkid[0], tHODO );
	htcorr_HCAL_HODO->Fill( HCALtdc[0], tHODO );
	htDiff_HODO_HCAL->Fill( HHdiff );
	haDiff_HODO_HCAL->Fill( HHAdiff );
	haHCAL_vs_HCALID->Fill( cblkid[0], HCALa[0] );
	htDiffBB_vs_HCALID->Fill( cblkid[0], HBBdiff );
	htDiffRF_vs_HCALID->Fill( cblkid[0], HRFdiff );
	htDiff_vs_HCALID->Fill( cblkid[0], HHdiff );
	htDiff_bb_vs_HCALID->Fill( cblkid[0], HBBdiff );
	htDiff_bbrf_vs_HCALID->Fill( cblkid[0], HBBdiff_rf );

	haDiff_vs_HCALID->Fill( cblkid[0], HHAdiff );
	htDiff[ idx ]->Fill( HHdiff );
	htDiff_bb[ idx ]->Fill( HBBdiff );
	//cout << HBBdiff << endl;
	htDiff_bbrf[ idx ]->Fill( HBBdiff_rf );
	haDiff[ idx ]->Fill( HHAdiff );
	//To compare columns 
	haDiff_vs_ADCint_cols[col]->Fill( e, HHAdiff ); 
	haDiff_vs_col[col]->Fill( HHAdiff );
	htDiff_vs_col[col]->Fill( HHdiff );
	//To compare JLAB and CMU and obtain averages for timewalk corrections with low stats
	if(ccol>3&&ccol<8){
	  htDiff_HODO_HCAL_JLAB->Fill( HHdiff );
	  haDiff_HODO_HCAL_JLAB->Fill( HHAdiff );
	  if(HHdiff>-82&&HHdiff<-69) htDiff_vs_ADCint_JLAB->Fill( e, HHdiff );
	  if(HHAdiff>42&&HHAdiff<60) haDiff_vs_ADCint_JLAB->Fill( e, HHAdiff );
	}else{
	  htDiff_HODO_HCAL_CMU->Fill( HHdiff );
	  haDiff_HODO_HCAL_CMU->Fill( HHAdiff );
	  if(HHdiff>-82&&HHdiff<-69) htDiff_vs_ADCint_CMU->Fill( e, HHdiff );
	  if(HHAdiff>42&&HHAdiff<60) haDiff_vs_ADCint_CMU->Fill( e, HHAdiff );
	}
	//Fill timing difference correlations with energy deposited for timewalk corrections
	htDiff_vs_ADCint[idx]->Fill( e, HHdiff );
	haDiff_vs_ADCint[idx]->Fill( e, HHAdiff );
      }
	
      //Events that pass the above cuts constitute elastics
      elasYield++;
    }
  }
  cout << "First loop complete. Time elapsed: " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;


  //Need to fit tdiff vs E histos and extract slope
  cout << endl << endl << "Extracting timewalk slopes and TOF mean.. " << endl;
  
  //By JLAB/CMU PMT construction
  TF1 *fitTW_tJLAB = new TF1( "fitTW_tJLAB", TW_fit, 0, 1 );
  htDiff_vs_ADCint_JLAB->Fit( "expo", "QNR+" );
  double TDC_JLAB_amp = fitTW_tJLAB->GetParameter(0);
  double TDC_JLAB_str = fitTW_tJLAB->GetParameter(1);
  TF1 *fitTW_aJLAB = new TF1( "fitTW_aJLAB", TW_fit, 0, 1 );
  haDiff_vs_ADCint_JLAB->Fit( "expo", "QNR+" );
  double ADC_JLAB_amp = fitTW_aJLAB->GetParameter(0);
  double ADC_JLAB_str = fitTW_aJLAB->GetParameter(1);
  TF1 *fitTW_tCMU = new TF1( "fitTW_tCMU", TW_fit, 0, 1 );
  htDiff_vs_ADCint_CMU->Fit( "expo", "QNR+" );
  double TDC_CMU_amp = fitTW_tCMU->GetParameter(0);
  double TDC_CMU_str = fitTW_tCMU->GetParameter(1);
  TF1 *fitTW_aCMU = new TF1( "fitTW_aCMU", TW_fit, 0, 1 );
  haDiff_vs_ADCint_CMU->Fit( "expo", "QNR+" );
  double ADC_CMU_amp = fitTW_aCMU->GetParameter(0);
  double ADC_CMU_str = fitTW_aCMU->GetParameter(1);

  //By Cell
  for(int i=0; i<kNcell; i++){

    int r = (i-1)/kNcols;
    int c = (i-1)%kNcols;

    TF1 *tw1;
    TF1 *tw2;

    //Fit HCal-TDC/Hodo-TDC difference vs energy deposited in HCal to obtain TDC timewalk
    if( htDiff_vs_ADCint[i]->GetEntries()>0 ){
      //Will need to add by-hand tuning here

      htDiff_vs_ADCint[i]->Fit("pol1","Q","",0.0,0.2);
      tw1=htDiff_vs_ADCint[i]->GetFunction("pol1");
      TvsE[i]=tw1->GetParameter(1);
      hTDCTimewalk->Fill(tw1->GetParameter(1));
    }else{
      TvsE[i]=0;
      cout << "No TDC data in channel " << i << "." << endl;
    }

    //Fit HCal-ADC-time/Hodo-TDC difference vs energy deposited in HCal to obtain ADCt timewalk
    if( haDiff_vs_ADCint[i]->GetEntries()>0 ){
      //Will need to add by-hand tuning here

      haDiff_vs_ADCint[i]->Fit("pol1","Q","",0.0,0.2);
      tw2=haDiff_vs_ADCint[i]->GetFunction("pol1");
      AtvsE[i]=tw2->GetParameter(1);
      hADCtTimewalk->Fill(tw2->GetParameter(1));
    }else{
      AtvsE[i]=0;
      cout << "No ADCt data in channel " << i << "." << endl;
    }
  }

  //Need to fit TOF and extract mean
  TF1 *tof1;
  hTOF2->Fit("gaus","Q","",0.0,10.0);
  tof1=hTOF2->GetFunction("gaus");
  double TOFmean = tof1->GetParameter(1);
  cout << endl << endl << "Slopes extracted. TOF mean = " << TOFmean << "." << endl;


  //Second loop over events with timewalk and TOF corrections applied
  cout << "Second loop over all data commencing.." << endl;
  progress = 0.;
  timekeeper = 0., timeremains = 0.;
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

      //Reconstructing the events with timewalk/TOF corrections by event
      double etheta = acos( BBtr_pz[0]/BBtr_p[0] );
      double ephi = atan2( BBtr_py[0], BBtr_px[0] );
      TVector3 vertex(0,0,BBtr_vz[0]);
      TLorentzVector Pbeam(0,0,E_e,E_e);
      TLorentzVector kprime(BBtr_px[0],BBtr_py[0],BBtr_pz[0],BBtr_p[0]);
      TLorentzVector Ptarg(0,0,0,M_p);
      TLorentzVector q = Pbeam - kprime;
      TLorentzVector PgammaN = Ptarg + q;
      double pel = E_e/(1.+E_e/M_p*(1.-cos(etheta)));
      double nu = E_e - BBtr_p[0];
      double pp = sqrt(pow(nu,2)+2.*M_p*nu);
      double phinucleon = ephi + PI;
      double thetanucleon = acos( (E_e - BBtr_p[0]*cos(etheta))/pp );
      double E_ep = sqrt( pow(M_e,2) + pow(BBtr_p[0],2) );
      double p_ep = BBtr_p[0];
      double Q2 = 2*E_e*E_ep*( 1-(BBtr_pz[0]/p_ep) );	
      double W = PgammaN.M();
      double E_pp = nu+M_p;
      double Enucleon = sqrt(pow(pp,2)+pow(M_p,2));
      double KE_p = nu;
      TVector3 pNhat( sin(thetanucleon)*cos(phinucleon),sin(thetanucleon)*sin(phinucleon),cos(thetanucleon));
      TVector3 HCAL_zaxis(sin(-HCal_th),0,cos(-HCal_th));
      TVector3 HCAL_xaxis(0,-1,0);
      TVector3 HCAL_yaxis = HCAL_zaxis.Cross(HCAL_xaxis).Unit();
      TVector3 HCAL_origin = HCal_d * HCAL_zaxis + hcalheight * HCAL_xaxis;
      double sintersect = ( HCAL_origin - vertex ).Dot( HCAL_zaxis ) / ( pNhat.Dot( HCAL_zaxis ) );      TVector3 HCAL_intersect = vertex + sintersect * pNhat;
      double yexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_yaxis );
      double xexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_xaxis );

      //Could improve by parsing the elist on first loop...
      double bbcal_time=0., hcal_time=0.;
      for(int ihit=0; ihit<TDCTndata; ihit++){
	if(TDCT_id[ihit]==5) bbcal_time=TDCT_tdc[ihit];
	if(TDCT_id[ihit]==0) hcal_time=TDCT_tdc[ihit];
      }
      double diff = hcal_time - bbcal_time; 

      double xDiff = HCALx - xexpect_HCAL;
      double yDiff = HCALy - yexpect_HCAL;
      bool nogo = spotcut==1 && (fabs(xDiff-dx0)>dx_sig || fabs(yDiff-dy0)>dy_sig);

      //Get timing offsets for tdc and adc from primary block in cluster
      if( fabs(W-W_mean)<W_sig&&fabs(diff+t_trig)<40&&nogo==false ){
	double vp2 = sqrt( 1/pow( M_p/pp, 2 ));
	double dp = sqrt( pow(HCal_d,2) + pow(HCALx,2) + pow(HCALy,2) );
	double tp2 = dp/vp2;
	double tp_diff = tp2 - TOFmean; //Obtain difference from TOF mean to obtain relative correction

	double HHdiff = HCALtdc[0]-tHODO;
	double HHAdiff = HCALa[0]-tHODO;
	double e = cblke[0];
	int idx = (int)cblkid[0]-1;
	int col = (int)ccol;

	//Calculate timewalk corrections via exponential fit
	double TW_JLAB_TDCcorr = TDC_JLAB_amp*exp( TDC_JLAB_str*e );
	//cout << "TW_JLAB_TDCcorr:" << TW_JLAB_TDCcorr << endl;
	double TW_JLAB_ADCcorr = ADC_JLAB_amp*exp( ADC_JLAB_str*e );
	//cout << "TW_JLAB_ADCcorr:" << TW_JLAB_ADCcorr << endl;
	double TW_CMU_TDCcorr = TDC_CMU_amp*exp( TDC_CMU_str*e );
	//cout << "TW_CMU_TDCcorr:" << TW_CMU_TDCcorr << endl;
	double TW_CMU_ADCcorr = ADC_CMU_amp*exp( ADC_CMU_str*e );
	//cout << "TW_CMU_ADCcorr:" << TW_CMU_ADCcorr << endl;

	//Fill histograms with the corrected time differences (base + timewalk + TOF)
	/*
	htDiff_corr[ idx ]->Fill( HHdiff + TvsE[ idx ]*e + tp_diff );
	haDiff_corr[ idx ]->Fill( HHAdiff + AtvsE[ idx ]*e + tp_diff );
	htDiff_vs_HCALID_corr->Fill( idx, HHdiff + TvsE[ idx ]*e + tp_diff  );
	haDiff_vs_HCALID_corr->Fill( idx, HHAdiff + AtvsE[ idx ]*e + tp_diff  );

	htcorr_HCAL_HODO_corr->Fill( HCALtdc[0] + TvsE[ idx ]*e + tp_diff, tHODO );
	hacorr_HCAL_HODO_corr->Fill( HCALa[0] + AtvsE[ idx ]*e + tp_diff, tHODO );
	htDiff_HODO_HCAL_corr->Fill( HHdiff + TvsE[ idx ]*e + tp_diff );
	haDiff_HODO_HCAL_corr->Fill( HHAdiff + AtvsE[ idx ]*e + tp_diff );
	if(ccol>3&&ccol<8){
	  htDiff_HODO_HCAL_JLAB_corr->Fill( HHdiff + TvsE[ idx ]*e + tp_diff );
	  haDiff_HODO_HCAL_JLAB_corr->Fill( HHAdiff + AtvsE[ idx ]*e + tp_diff, tHODO );
	}else{
	  htDiff_HODO_HCAL_CMU_corr->Fill( HHdiff + TvsE[ idx ]*e + tp_diff );
	  haDiff_HODO_HCAL_CMU_corr->Fill( HHAdiff + AtvsE[ idx ]*e + tp_diff, tHODO );
	}
	*/
	//For now, with low stats, correct with timewalk by PMT type
	htDiff_corr[ idx ]->Fill( HHdiff + TvsE[ idx ]*e + tp_diff );
	haDiff_corr[ idx ]->Fill( HHAdiff + AtvsE[ idx ]*e + tp_diff );
	htDiff_vs_HCALID_corr->Fill( idx, HHdiff + TvsE[ idx ]*e + tp_diff  );
	haDiff_vs_HCALID_corr->Fill( idx, HHAdiff + AtvsE[ idx ]*e + tp_diff  );

	htcorr_HCAL_HODO_corr->Fill( HCALtdc[0] + TvsE[ idx ]*e + tp_diff, tHODO );
	hacorr_HCAL_HODO_corr->Fill( HCALa[0] + AtvsE[ idx ]*e + tp_diff, tHODO );
	htDiff_HODO_HCAL_corr->Fill( HHdiff + TvsE[ idx ]*e + tp_diff );
	haDiff_HODO_HCAL_corr->Fill( HHAdiff + AtvsE[ idx ]*e + tp_diff );
	if(ccol>3&&ccol<8){
	  htDiff_HODO_HCAL_JLAB_corr->Fill( HHdiff + TvsE[ idx ]*e + tp_diff );
	  haDiff_HODO_HCAL_JLAB_corr->Fill( HHAdiff + AtvsE[ idx ]*e + tp_diff, tHODO );
	}else{
	  htDiff_HODO_HCAL_CMU_corr->Fill( HHdiff + TvsE[ idx ]*e + tp_diff );
	  haDiff_HODO_HCAL_CMU_corr->Fill( HHAdiff + AtvsE[ idx ]*e + tp_diff, tHODO );
	}

	haDiff_vs_col_corr[col]->Fill( HHAdiff );
	htDiff_vs_col_corr[col]->Fill( HHdiff );
      }
    }
  }
  
  cout << endl;

  cout << "tFitMin: " << tFitMin << endl;

  for(int i=0; i<kNcell; i++){
    
    int r = (i-1)/kNcols;
    int c = (i-1)%kNcols;
    
    TF1 *f1;
    TF1 *f2;
    
    //Fit HCal-TDC/Hodo-TDC difference histogram by cell and obtain offset parameters
    if( htDiff_corr[i]->GetEntries()>0 ){
      
      htDiff_corr[i]->Fit("gaus","Q","",-150,0); //May need to tune fits here
      f1=htDiff_corr[i]->GetFunction("gaus");
      htDiff_corr[i]->SetTitle(Form("TDCGoodFitMean:%f",f1->GetParameter(1)));
      TDCoffsets[i] = -75. - f1->GetParameter(1);
      TDCsig[i] = f1->GetParameter(2);
      //double cs = f1->GetChisquare();
      htDiff_corr[i]->SetName(Form("htDiff_corr_bl%d",i));		
    }else{   
      TDCoffsets[i] = 0.;
    }
    
    //Fit HCal-ADC-time/Hodo-TDC difference histogram by cell and obtain offset parameters
    if( haDiff[i]->GetEntries()>0 ){
      
      haDiff[i]->Fit("gaus","Q","",0,150); //May need to tune fits here
      f2=haDiff[i]->GetFunction("gaus");
      haDiff[i]->SetTitle(Form("ADCtGoodFitMean:%f",f2->GetParameter(1)));
      ADCtoffsets[i] = 50. - f2->GetParameter(1);
      ADCtsig[i] = f2->GetParameter(2);
      //double cs = f2->GetChisquare();
      haDiff[i]->SetName(Form("haDiff_bl%d",i));
    }else{
      ADCtoffsets[i] = 0.;
    }
  }
  

  //Make offsets input-ready and improve readability
  for(int i=0; i<kNcell; i++){
    hTDCoffsets_vs_HCALID->Fill(i,TDCoffsets[i]);
    holdTDCoffsets_vs_HCALID->Fill(i,oldTDCoffsets[i]*TDCCalib);  
    combTDCoffsets[i] = TDCoffsets[i]/TDCCalib+oldTDCoffsets[i];
    hcombTDCoffsets_vs_HCALID->Fill(i,combTDCoffsets[i]);     
    hTDCsig_vs_HCALID->Fill(i,TDCsig[i]);
    
    hADCtoffsets_vs_HCALID->Fill(i,ADCtoffsets[i]);
    combADCtoffsets[i] = ADCtoffsets[i]+oldADCtoffsets[i];
    hADCtsig_vs_HCALID->Fill(i,ADCtsig[i]);
  }
 
  //Make canvas to hold all fits for comparison to HCal geometry
  if( diag==1 ){
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
  }

  //Construct graphs for TDC timing resolution
  const Int_t xN = 288; //total channels
  Double_t posErr[xN] = {0.};
  TF1 *f1;
  //For TDC corrected
  Double_t X[xN];
  Double_t Xval[xN];
  Double_t Xerr[xN];
  TH1D *Xslice[xN];
  Double_t sigmaX[2] = {0,0};
  Double_t minSigmaX[2] = {0,100};
  for( int x=0; x<xN; x++ ){
    X[x] = x;
    Xslice[x] = htDiff_vs_HCALID->ProjectionY(Form("Xslice_%d",x+1),x+1,x+1);
    Xslice[x]->Fill(-75); //Add one point to fix fits
    Xslice[x]->Fit("gaus","Q","",-80,-65);
    f1=Xslice[x]->GetFunction("gaus");
    if(Xslice[x]->GetEntries()>100){
      Xval[x] = f1->GetParameter(1);
      Xerr[x] = f1->GetParameter(2);
      sigmaX[0]++;
      sigmaX[1]+=Xerr[x];
      if( Xerr[x]<minSigmaX[1] ){
	minSigmaX[0] = x;
	minSigmaX[1] = Xerr[x];
      }
    }else{
      Xval[x] = -10;
      Xerr[x] = 0.;
    }
  }
  cout << endl << "Avg TDC sigma: " << sigmaX[1]/sigmaX[0] << endl << endl;
  cout << "Best TDC: " << minSigmaX[0] << " sig: " << minSigmaX[1] << endl;
  TGraphErrors *ctdcres_Ch = new TGraphErrors( xN, X, Xval, posErr, Xerr );
  ctdcres_Ch->GetXaxis()->SetLimits(0,xN);  
  ctdcres_Ch->GetYaxis()->SetLimits(-100,-40);
  ctdcres_Ch->SetTitle("Timing Resolution TDC");
  ctdcres_Ch->GetXaxis()->SetTitle("Channel");
  ctdcres_Ch->GetYaxis()->SetTitle("TDC_{HCal}-TDC_{HODO} (ns)");
  ctdcres_Ch->SetMarkerStyle(20); // idx 20 Circles, idx 21 Boxes
  ctdcres_Ch->Write("ctdcres_Ch");

  //For TDC uncorrected
  Double_t Xu[xN];
  Double_t Xvalu[xN];
  Double_t Xerru[xN];
  TH1D *Xsliceu[xN];
  Double_t sigmaXu[2] = {0,0};
  Double_t minSigmaXu[2] = {0,100};
  for( int x=0; x<xN; x++ ){
    Xu[x] = x-0.5;
    Xsliceu[x] = htHCAL_vs_HCALID->ProjectionY(Form("Xsliceu_%d",x+1),x+1,x+1);
    Xsliceu[x]->Fill(-75); //Add one point to fix fits
    Xsliceu[x]->Fit("gaus","Q","",-80,-65);
    f1=Xsliceu[x]->GetFunction("gaus");
    if(Xsliceu[x]->GetEntries()>100){
      Xvalu[x] = f1->GetParameter(1);
      Xerru[x] = f1->GetParameter(2);
      sigmaXu[0]++;
      sigmaXu[1]+=Xerr[x];
      if( Xerru[x]<minSigmaXu[1] ){
	minSigmaXu[0] = x;
	minSigmaXu[1] = Xerr[x];
      }
    }else{
      Xvalu[x] = -10;
      Xerru[x] = 0.;
    }
  }
  cout << endl << "Avg TDC uncorrected sigma: " << sigmaXu[1]/sigmaXu[0] << endl << endl;
  cout << "Best TDC: " << minSigmaXu[0] << " sig: " << minSigmaXu[1] << endl;
  TGraphErrors *ctdcresu_Ch = new TGraphErrors( xN, Xu, Xvalu, posErr, Xerru );
  ctdcresu_Ch->GetXaxis()->SetLimits(0,xN);  
  ctdcresu_Ch->GetYaxis()->SetLimits(-100,-40);
  ctdcresu_Ch->SetTitle("Timing Resolution TDC Uncorrected");
  ctdcresu_Ch->GetXaxis()->SetTitle("Channel");
  ctdcresu_Ch->GetYaxis()->SetTitle("TDC_{HCal} (ns)");
  ctdcresu_Ch->SetMarkerColor(2);
  ctdcresu_Ch->SetLineColor(2);
  ctdcresu_Ch->SetMarkerStyle(24); // idx 24 open circles
  ctdcresu_Ch->Write("ctdcresu_Ch");

  //For ADC time
  Double_t Y[xN];
  Double_t Yval[xN];
  Double_t Yerr[xN];
  TH1D *Yslice[xN];
  Double_t sigmaY[2] = {0,0};
  Double_t minSigmaY[2] = {0,100};
  for( int x=0; x<xN; x++ ){
    Y[x] = x;
    Yslice[x] = haDiff_vs_HCALID->ProjectionY(Form("Yslice_%d",x+1),x+1,x+1);
    Yslice[x]->Fill(55); //Add point to fix fits
    Yslice[x]->Fit("gaus","Q","",40,60);
    f1=Yslice[x]->GetFunction("gaus");
    if(Yslice[x]->GetEntries()>100){
      Yval[x] = f1->GetParameter(1);
      Yerr[x] = f1->GetParameter(2);
      sigmaY[0]++;
      sigmaY[1]+=Yerr[x];
      if( Yerr[x]<minSigmaY[1] ){
	minSigmaY[0] = x;
	minSigmaY[1] = Yerr[x];
      }
    }else{
      Yval[x] = 80;
      Yerr[x] = 0.;
    }
  }
  cout << endl << "Avg ADCt sigma: " << sigmaY[1]/sigmaY[0] << endl << endl;
  cout << "Best ADC: " << minSigmaY[0] << " sig: " << minSigmaY[1] << endl;
  TGraphErrors *cadcres_Ch = new TGraphErrors( xN, Y, Yval, posErr, Yerr );
  cadcres_Ch->GetXaxis()->SetLimits(0,xN);  
  cadcres_Ch->GetYaxis()->SetLimits(20,80);
  cadcres_Ch->SetTitle("Timing Resolution ADCt");
  cadcres_Ch->GetXaxis()->SetTitle("Channel");
  cadcres_Ch->GetYaxis()->SetTitle("ADCt_{HCal}-TDC_{HODO} (ns)");
  cadcres_Ch->SetMarkerStyle(20); // idx 20 Circles, idx 21 Boxes
  cadcres_Ch->Write("cadcres_Ch");


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
  /*
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
  */
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
