//SSeeds 6.30.22 - Post-production - Copy of timing calibration code which employs basic cuts on elastic events to obtain timing offsets for HCal ADC time and TDC. Leaving out proton spot cuts at this stage to obtain better statistics. Developed to characterize the double peaking structure in HCal data

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

void doublePeak( const char *configfilename="sPeak4.cfg", int run = -1 ){
  
  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );
  
  // Get the date
  string date = getDate();
  
  // Declare Chain for many root files
  TChain *C = new TChain("T");
  
  // Path for output file
  string paramsPath = "/w/halla-scshelf2102/sbs/seeds/HCal_replay/hcal/Hodo_HCal_Coorelations/sPeak.txt";
  
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
      cout << "Added file " << currentline << " to chain.." << endl;
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
  TString outputfilename = Form("doublePeakOut_%s.root", date.c_str());
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
  // Peak cut histos
  TH1D *hTBits = new TH1D( "hTBits", "TrigBits", 40, 0, 40);  
  TH1D *honeblk_dt = new TH1D( "honeblk_dt", "TDC All Channels, One Block Clusters", 300, -150, 0);
  TH1D *hclusblk_dt = new TH1D( "hclusblk_dt", "Time Difference Between Blocks in Cluster (ns)", 1000, -50, 50);
  TH1D *hblk_dt = new TH1D( "hblk_dt", "Time Difference Between Primary Hit and Next Adjacent Block", 1000, -50, 50);
  TH1D *hclusblkcut_dt = new TH1D( "hclusblkcut_dt", "Time Difference Between Primary Hit and First Adjacent Block", 4000, -1000, 1000);
  TH1D *hclusblkcutv_dt = new TH1D( "hclusblkcutv_dt", "Time Difference Between Primary Hit and First Vert. Adjacent Block", 4000, -1000, 1000);
  TH1D *hclusblkcuth_dt = new TH1D( "hclusblkcuth_dt", "Time Difference Between Primary Hit and First Horiz. Adjacent Block", 4000, -1000, 1000);
  TH2D *hclusdt_vs_x = new TH2D("hclusdt_vs_x","; x_{HCAL} (m);tdc_{block_1}-tdc{block_2} (ns)", 500, -2.25, 1.5, 1000, -50, 50 );
  TH2D *hclusdt_vs_y = new TH2D("hclusdt_vs_y","; y_{HCAL} (m);tdc_{block_1}-tdc{block_2} (ns)", 250, -1, 1, 1000, -50, 50 ); 
  TH2D *hclusdt_vs_row = new TH2D("hclusdt_vs_row","; ROW_{HCAL};tdc_{block_1}-tdc{block_2} (ns)", 24, 0, 24, 1000, -50, 50 );
  TH2D *hclusdt_vs_col = new TH2D("hclusdt_vs_col","; COL_{HCAL} (m);tdc_{block_1}-tdc{block_2} (ns)", 12, 0, 12, 1000, -50, 50 );  
  TH2D *hclusdt_vs_HCALID = new TH2D("hclusdt_vs_HCALID",";channel;tdc_{block_1}-tdc_{block_2} (ns)", 288, 0, 288, 1000, -50, 50 );

  TH1D *h1p_E = new TH1D( "h1p_E", "HCal E, First Peak; GeV", 1000, 0, 0.5);
  TH1D *h2p_E = new TH1D( "h2p_E", "HCal E, Second Peak; GeV", 1000, 0, 0.5);
  TH1D *h1pb = new TH1D( "h1pb", "LED Bit; Bit", 20, 0, 20);
  TH1D *h2pb = new TH1D( "h2pb", "LED Bit; Bit", 20, 0, 20);
  TH1D *h1ptm = new TH1D( "h1ptm", "TDC Multiplicity; Number Hits", 20, 0, 20);
  TH1D *h2ptm = new TH1D( "h2ptm", "TDC Multiplicity; Number Hits", 20, 0, 20);
  TH1D *h1pttb = new TH1D( "h1pttb", "TDC Time Highest E Blk (SBSOFFLINE); ns", 1000, -500, 500);
  TH1D *h2pttb = new TH1D( "h2pttb", "TDC Time Highest E Blk (SBSOFFLINE); ns", 1000, -500, 500);
  TH1D *h1p_shE = new TH1D( "h1p_shE", "Shower E, First Peak; GeV", 1000, 0, 5);
  TH1D *h2p_shE = new TH1D( "h2p_shE", "Shower E, Second Peak; GeV", 1000, 0, 5);
  TH1D *h1p_psE = new TH1D( "h1p_psE", "Preshower E, First Peak; GeV", 1000, 0, 5);
  TH1D *h2p_psE = new TH1D( "h2p_psE", "Preshower E, Second Peak; GeV", 1000, 0, 5);
  TH1D *h1p_nblk = new TH1D( "h1p_nblk", "HCal Cluster Size, First Peak; N", 10, 0, 10);
  TH1D *h2p_nblk = new TH1D( "h2p_nblk", "HCal Cluster Size, Second Peak; N", 10, 0, 10);  
  TH1D *h1p_Eexp = new TH1D( "h1p_Eexp", "HCal KE Expected, First Peak; GeV", 500, 1, 2);
  TH1D *h2p_Eexp = new TH1D( "h2p_Eexp", "HCal KE Expected, Second Peak; GeV", 500, 1, 2);
  TH2D *h1p_Hpos = new TH2D("h1p_Hpos",";y_{HCAL} (m); x_{HCAL} (m)", 250, -1, 1, 250, -3, 2 );  
  TH2D *h2p_Hpos = new TH2D("h2p_Hpos",";y_{HCAL} (m); x_{HCAL} (m)", 250, -1, 1, 250, -3, 2 );
  TH1D *h1p_Hx = new TH1D("h1p_Hx",";x_{HCAL} (m)", 250, -3, 2 );  
  TH1D *h2p_Hx = new TH1D("h2p_Hx",";x_{HCAL} (m)", 250, -3, 2 );
  TH1D *h1p_Hy = new TH1D("h1p_Hy",";y_{HCAL} (m)", 250, -1, 1 );  
  TH1D *h2p_Hy = new TH1D("h2p_Hy",";y_{HCAL} (m)", 250, -1, 1 );
  TH1D *h1p_Bx = new TH1D("h1p_Bx",";x_{BBtr} (m)", 100, -1, 1 );  
  TH1D *h2p_Bx = new TH1D("h2p_Bx",";x_{BBtr} (m)", 100, -1, 1 );
  TH1D *h1p_By = new TH1D("h1p_By",";y_{BBtr} (m)", 100, -0.5, 0.5 );  
  TH1D *h2p_By = new TH1D("h2p_By",";y_{BBtr} (m)", 100, -0.5, 0.5 );
  TH2D *h1p_BBpos = new TH2D("h1p_BBpos",";y_{expect} (m); x_{expect} (m)", 250, -2, 1, 250, -3, 2 );  
  TH2D *h2p_BBpos = new TH2D("h2p_BBpos",";y_{expect} (m); x_{expect} (m)", 250, -2, 1, 250, -3, 2 );
  TH1D *h1p_BBx = new TH1D("h1p_BBx",";x_{HCAL} (m)", 250, -3, 2 );  
  TH1D *h2p_BBx = new TH1D("h2p_BBx",";x_{HCAL} (m)", 250, -3, 2 );
  TH1D *h1p_BBy = new TH1D("h1p_BBy",";y_{HCAL} (m)", 250, -2, 1 );  
  TH1D *h2p_BBy = new TH1D("h2p_BBy",";y_{HCAL} (m)", 250, -2, 1 );
  TH2D *h1p_tDiff = new TH2D("h1p_tDiff",";Channel;TDC_{HCAL}-TDC_{HODO} (ns)",288,0,288,300,-150,0);
  TH2D *h2p_tDiff = new TH2D("h2p_tDiff",";Channel;TDC_{HCAL}-TDC_{HODO} (ns)",288,0,288,300,-150,0);
  TH2D *h1p_aDiff = new TH2D("h1p_aDiff",";Channel;ADCtime_{HCAL}-TDC_{HODO} (ns)",288,0,288,300,20,80);  
  TH2D *h2p_aDiff = new TH2D("h2p_aDiff",";Channel;ADCtime_{HCAL}-TDC_{HODO} (ns)",288,0,288,300,20,80);
  TH2D *h1pap_ID = new TH2D("h1pap_ID",";Channel; Ped Sub intADC",288,0,288,2500,-500,2000);  
  TH2D *h2pap_ID = new TH2D("h2pap_ID",";Channel; Ped Sub intADC",288,0,288,2500,-500,2000);  
  TH2D *h1p_t = new TH2D("h1p_t",";Channel;TDC_{HCAL}-TDC_{HODO} (ns)",288,0,288,300,-150,0);
  TH2D *h2p_t = new TH2D("h2p_t",";Channel;TDC_{HCAL}-TDC_{HODO} (ns)",288,0,288,300,-150,0);
  TH2D *h1p_a = new TH2D("h1p_a",";Channel;ADCtime_{HCAL}-TDC_{HODO} (ns)",288,0,288,300,20,80);  
  TH2D *h2p_a = new TH2D("h2p_a",";Channel;ADCtime_{HCAL}-TDC_{HODO} (ns)",288,0,288,300,20,80);
  TH2D *h1pT_vs_nblk = new TH2D("h1pT_vs_nblk",";Block N;TDC_{HCAL} (ns)",10,0,10,300,-150,0);
  TH2D *h2pT_vs_nblk = new TH2D("h2pT_vs_nblk",";Block N;TDC_{HCAL} (ns)",10,0,10,300,-150,0);
  TH2D *hT_vs_nblk = new TH2D("hT_vs_nblk",";Block N;TDC_{HCAL} (ns)",10,0,10,300,-150,0);
  
  
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
  TH1D *hTBBt_elas = new TH1D( "hTBBt_elas","BBCal Trigger Elastic Cut (L1A) (ns)", 500, 340, 390 );
  TH1D *hTRF = new TH1D( "hTRF","RF Signature (ns)", 1700, -10, 160 );
  TH1D *hTRFmod = new TH1D( "hTRFmod","RF Signature Modulo 4ns (ns)", 200, -10, 10 );
  TH1D *hTEDTM = new TH1D( "hTEDTM","Electronic Dead Time Monitor Signal (ns)", 2000, -1, 2 );
  TH1D *hTBBlo = new TH1D( "hTBBlo","BBCal Lo Trigger (ns)", 199, 1, 200 );
  TH1D *hTBBhiv = new TH1D( "hTBBhiv","BBCal Lo (ns)", 2010, -10, 2000 );
  TH1D *hTHCAL = new TH1D( "hTHCAL","HCal time (ns)", 1000, 800, 900 ); //trigger
  TH1D *hTHODO = new TH1D( "hTHODO","Hodo mean TDC time (ns)", 2000, -1000, 1000 );
  TH2D *hTHODO_vs_BBtrig = new TH2D( "hTHODO_vs_BBtrig",";BB Trigger Time (ns); Hodo mean TDC time (ns)" , 500, 340, 390, 200, -20, 20 );
  TH1D *hTHCALvRF = new TH1D( "hTHCALvRF","HCal time - RF Signature Modulo 4ns (ns)", 1000, 800, 900 );
  TH1D *hDiff = new TH1D( "hDiff","HCal time - BBCal time (ns)", 1300, -500, 800 );
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
  TH2D *hTOFvID = new TH2D( "hTOFvID", "Time of Flight vs HCal ID", 288, 0, 288, 200., 0., 20. );
  TH2D *hTOF2vID = new TH2D( "hTOFvID2", "Time of Flight vs HCal ID (nu method)", 288, 0, 288, 200., 0., 20. );
  TH2D *hdxdy_HCAL = new TH2D("hdxdy_HCAL",";y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)", 250, -5.0, 5.0, 250, -10, 10 );
  TH2D *hXY_HCAL = new TH2D("hXY_HCAL",";y_{HCAL} (m); x_{HCAL} (m)", 250, -5.0, 5.0, 250, -10, 10 );
  TH1D *hvz_cut = new TH1D("hvz_cut",";vertex z (m);", 250,-0.125,0.125);
  
  // Timing histograms
  TH2D *htcorr_HCAL_HODO = new TH2D("htcorr_HCAL_HODO",";TDC_{HCAL} (ns);TDC_{HODO} (ns)", 150, -150, 0, 150, -15, 15);  
  TH2D *hacorr_HCAL_HODO = new TH2D("hacorr_HCAL_HODO",";ADCt_{HCAL} (ns);TDC_{HODO} (ns)", 300, 0, 150, 150, -15, 15);
  TH2D *htcorr_HCAL_HODO_corr = new TH2D("HCAL/HODO TDC/TDC Correlation (corrected)",";TDC_{HCAL} (ns);TDC_{HODO} (ns)", 150, -150, 0, 150, -15, 15);  
  TH2D *hacorr_HCAL_HODO_corr = new TH2D("HCAL/HODO ADCt/TDC Correlation (corrected)",";ADCt_{HCAL} (ns);TDC_{HODO} (ns)", 300, 0, 150, 150, -15, 15);
  TH1D *htHCAL = new TH1D("HCal_TDC","HCal TDC (All Channels);ns",300,-150,0);
  TH1D *htHCAL_second = new TH1D("HCal_TDC_second","HCal TDC, Second Block (All Channels);ns",300,-150,0);
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
  TH2D *htDiffBB_vs_HCALID = new TH2D("htDiffBB_vs_HCALID",";Channel;TDC_{HCAL}-BBCalt (ns)",288,0,288,1000,-1000,0);
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
      
      hTBits->Fill(TBits);

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
      
      //cout << diff-t_trig << endl;

      //Get timing offsets for tdc and adc from primary block in cluster
      //if( cosmic==1 || (fabs( W-W_mean )<W_sig&&fabs( diff-t_trig )<40 && nogo==false) ){
      if( fabs( W-W_mean )<W_sig ){
	
	//Fill some basic histos
	hTBBt_elas->Fill( bbcal_time );
	hW_cuts->Fill( W );
	hCblkID->Fill( cblkid[0] );
	hCCol->Fill( ccol );
	hvz_cut->Fill( BBtr_vz[0] );
	//Draw some location plots for better cuts on elastics later on
	hdxdy_HCAL->Fill( yDiff, xDiff );
	hdy->Fill( yDiff );
	hdx->Fill( xDiff );
	
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
	hTOF2->Fill( tp2 );
	hTOFvID->Fill( cblkid[0], tp );
	hTOF2vID->Fill( cblkid[0], tp );
	
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
	hTHODO->Fill( tHODO );
	hTHODO_vs_BBtrig->Fill( bbcal_time, tHODO );
	hHCALtdc_vs_HODOtdc->Fill( tHODO, HCALtdc[0] );
	htHCAL->Fill( HCALtdc[0] );
	htHCAL_second->Fill( HCALtdc[1] );
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
	
	if( waveform==1 ){
	  double adj = cblkid[0] + 1;
	  if( adj<=288 ) hblk_dt->Fill( HCALtdc[0] - HCALTtdc[ (int)adj ] );
	}
	double dt = HCALtdc[0] - HCALtdc[1];

	// Check if next highest energy block after primary is immediately adjacent to primary block
	if( (cblkid[1]==cblkid[0]+1)||
	    (cblkid[1]==cblkid[0]-1)||
	    (cblkid[1]==cblkid[0]+12)||
	    (cblkid[1]==cblkid[0]-12) ){
	  //double dt = HCALtdc[0] - HCALtdc[1];
	  //cout << dt << endl;
	  hclusblkcut_dt->Fill( dt );
	  hclusdt_vs_x->Fill( HCALx, dt );
	  hclusdt_vs_y->Fill( HCALy, dt );
	  hclusdt_vs_row->Fill( row, dt );
	  hclusdt_vs_col->Fill( col, dt );
	  hclusdt_vs_HCALID->Fill( cblkid[0], dt );
	}
	if( (cblkid[1]==cblkid[0]+1)||
	    (cblkid[1]==cblkid[0]-1) ){
	  hclusblkcuth_dt->Fill( dt );
	}
	if( (cblkid[1]==cblkid[0]+12)||
	    (cblkid[1]==cblkid[0]-12) ){
	  hclusblkcutv_dt->Fill( dt );
	}
	if( nblk==1 ) honeblk_dt->Fill( HCALtdc[0] );

	for( int i=0; i<nblk; i++){
	  hT_vs_nblk->Fill( i, HCALtdc[i] );
	  if( i>0 ) hclusblk_dt->Fill( HCALtdc[0] - HCALtdc[i] );
	}
	if( HHdiff < -66 && HHdiff > -82){
	  h1p_E->Fill( HCALe );
	  h1p_psE->Fill( BBps_e );
	  h1p_shE->Fill( BBsh_e );
	  h1p_Eexp->Fill( KE_p );
	  h1p_Hpos->Fill( HCALy, HCALx);
	  h1p_BBpos->Fill( yexpect_HCAL, xexpect_HCAL);
	  h1p_tDiff->Fill( cblkid[0], HHdiff );
	  h1p_t->Fill( cblkid[0], HCALtdc[0]);
	  h1p_aDiff->Fill( cblkid[0], HHAdiff );
	  h1p_a->Fill( cblkid[0], HCALa[0]);
	  h1p_nblk->Fill( nblk );
	  
	  if( waveform==1 ){
	    h1pap_ID->Fill( cblkid[0], HCALap[ (int)cblkid[0] ] );
	    h1pb->Fill( HCALbit );
	    h1ptm->Fill( HCALtm[0] );
	    h1pttb->Fill( HCALttb );
	  }
	  
	  if( nblk>0 ){
	    h1p_Hx->Fill( HCALx );  
	    h1p_Hy->Fill( HCALy );
	    h1p_Bx->Fill( BBtr_x[0] );  
	    h1p_By->Fill( BBtr_y[0] ); 
	    h1p_BBx->Fill( xexpect_HCAL );  
	    h1p_BBy->Fill( yexpect_HCAL ); 
	  }
	  for( int i=0; i<nblk; i++){
	    h1pT_vs_nblk->Fill( i, HCALtdc[i] );
	  }
	}
	
	if( HHdiff >= -66 && HHdiff < -54){
	  h2p_E->Fill( HCALe );
	  h2p_psE->Fill( BBps_e );
	  h2p_shE->Fill( BBsh_e );
	  h2p_Eexp->Fill( KE_p );
	  h2p_Hpos->Fill( HCALy, HCALx);
	  h2p_BBpos->Fill( yexpect_HCAL, xexpect_HCAL);
	  h2p_tDiff->Fill( cblkid[0], HHdiff );
	  h2p_t->Fill( cblkid[0], HCALtdc[0]);
	  h2p_aDiff->Fill( cblkid[0], HHAdiff );
	  h2p_a->Fill( cblkid[0], HCALa[0]);
	  h2p_nblk->Fill( nblk );
	  
	  if( waveform==1 ){
	    h2pap_ID->Fill( cblkid[0], HCALap[ (int)cblkid[0] ] );
	    h2pb->Fill( HCALbit );
	    h2ptm->Fill( HCALtm[0] );
	    h2pttb->Fill( HCALttb );
	  }
	  
	  if( nblk>0 ){
	    h2p_Hx->Fill( HCALx ); 
	    h2p_Hy->Fill( HCALy );
	    h2p_Bx->Fill( BBtr_x[0] ); 
	    h2p_By->Fill( BBtr_y[0] ); 
	    h2p_BBx->Fill( xexpect_HCAL );  
	    h2p_BBy->Fill( yexpect_HCAL ); 
	  }
	  for( int i=0; i<nblk; i++){
	    h2pT_vs_nblk->Fill( i, HCALtdc[i] );
	  }
	}
	
	//Events that pass the above cuts constitute elastics
	elasYield++;
	
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
      
    }
  }

  cout << "Analysis complete. Time elapsed: " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;
  
  cout << endl;
  
  //Write out diagnostic histos and print to console
  fout->Write();
  
  cout << endl << endl << "Elastic yield for analyzed runs: " << elasYield << ". Total events analyzed: " << Nevents << ". Total TDC misses: " << TDCmiss << "." << endl << endl;
  
  cout << "Timing offset analysis complete and parameters written to file." << endl;
  
  elist->Delete();
  st->Stop();
  
  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;
  
}
