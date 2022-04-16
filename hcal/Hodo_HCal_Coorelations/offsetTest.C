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
const double kdBlock = 0.152; // Width and height of each module including distance between
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

void timingCal( const char *configfilename, int run = -1 ){
  
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
  
  cout << "Setup parameters loaded." << endl;

  // Declare general detector parameters
  double BBtr_p[maxTracks], BBtr_px[maxTracks], BBtr_py[maxTracks], BBtr_pz[maxTracks];
  double BBtr_vz[maxTracks], BBtr_chi2[maxTracks];
  double BBtr_n, BBps_x, BBps_y, BBps_e, BBsh_x, BBsh_y, BBsh_e;

  double TDCT_id[maxTdcChan], TDCT_tdc[maxTdcChan]; 
  int TDCTndata;

  double HCALx, HCALy, HCALe;
  double HCALtdc[kNcell], HCALa[kNcell];
  double crow, ccol, nblk;
  double cblkid[kNcell], cblke[kNcell];

  double tHODO;
  double nClusHODO;

  double oldTDCOffsets[kNcell] = 
    {51.6052,53.1381,75.711,62.3824,-64.8101,-49.5039,-16.1844,-33.8902,67.763,65.2803,74.9802,46.4545,78.6722,78.2396,63.5186,78.6833,-36.4883,-30.7707,-24.0491,-16.4224,85.6185,68.1837,71.3112,56.5909,64.0445,52.6813,87.4801,78.7829,-40.4086,-33.1198,-37.1459,-42.5131,97.7638,74.9225,60.3911,54.1738,95.8128,98.0925,122.216,111.835,-4.14176,-9.11937,0,-0.0246738,109.739,122.598,92.6682,74.059,-12.3159,-9.2373,-8.41619,-0.39065,-72.9534,-99.4989,-105.177,-108.026,-3.92294,5.20389,-1.74325,-10.1903,110.999,113.306,126.067,120.815,-100.918,-84.1958,19.4866,16.1717,137.576,127.848,7.37405,8.3169,42.849,77.3066,88.1684,93.4548,-36.7794,-24.315,-21.1369,-33.6198,126.699,114.999,79.33,63.7887,96.7496,88.5692,99.0824,92.4448,-6.87795,-13.2878,-19.3181,-10.1457,87.3303,111.9,66.7838,59.101,111.863,125.724,113.91,110.842,7.28227,4.90495,21.2512,12.3495,112.563,121.519,113.138,105.109,84.9052,95.6989,107.758,100.946,18.9589,20.2386,20.9004,20.8449,122.96,91.8478,110.502,109.098,132.579,108.904,107.074,106.125,36.3715,23.224,40.6818,38.4267,114.892,120.195,107.725,88.8613,105.714,112.342,121.262,109.64,3.5718,-5.55847,14.629,16.5227,127.141,121.331,94.1712,92.4037,-28.3475,-12.4915,-16.176,-2.35038,-118.309,-102.759,-112.91,-105.626,11.6323,3.00557,-2.78524,-14.6378,128.323,125.449,120.715,104.829,51.3949,14.4151,54.3329,41.447,112.008,114.84,142.813,108.737,97.4597,121.844,115.635,103.866,23.6189,30.5603,44.9171,33.1805,112.155,89.7849,104.902,117.556,1.23497,-10.884,128.598,137.733,28.8473,7.63449,-81.1601,-81.3845,127.848,144.716,112.554,110.205,95.3282,92.8283,131.052,139.466,21.4253,46.3608,5.66095,-14.1183,150.785,131.072,134.173,106.736,-73.7123,-72.157,0.825379,-28.8035,59.844,53.0085,-122.239,-130.565,9.17359,-13.9012,150.787,135.087,102.953,93.79,137.682,119.49,-107.541,-90.8528,29.1073,41.3303,137.553,126.15,-12.7449,-39.9297,117.381,141.75,151.81,150.345,-130.425,-122.595,63.285,42.3564,153.701,166.062,-39.0284,-43.4599,79.6345,105.374,99.1642,115.317,-108.623,-104.695,34.7333,23.3209,114.342,112.427,-39.0258,-46.3406,81.9736,85.6216,122.161,110.249,16.1559,37.152,63.5,23.8763,111.813,100.931,105.722,74.6418,103.489,135.58,122.72,113.029,34.5295,48.8976,31.2923,36.1094,126.375,109.23,101.896,103.769,-56.248,-51.1446,98.0804,87.691,22.5598,39.2716,-117.507,-145.363,114.921,94.8606,85.7814,109.389};
  

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

  cout << "Tree variables linked." << endl;

  // Define a clock to check macro processing time
  TStopwatch *sw = new TStopwatch();
  sw->Start();
  
  // Declare outfile
  TFile *fout = new TFile( "timingCal_out.root", "RECREATE" );
  
  // Initialize vectors and arrays
  double TDCoffsets[kNcell] = {0.0};
  double ADCToffsets[kNcell] = {0.0};
  double TDCsig[kNcell] = {0.0};
  double ADCTsig[kNcell] = {0.0};
  double TvsE[kNcell] = {0.0};
  double ATvsE[kNcell] = {0.0};

  // Initialize histograms	
  TH1D *hDeltaE = new TH1D( "hDeltaE","1.0-Eclus/p_rec", 100, -1.5, 1.5 );
  TH1D *hHCALe = new TH1D( "hHCALe","E HCal Cluster E", 400, 0., 4 );
  TH1D *hHODOnclus = new TH1D( "hHODOnclus","Number of Hodoscope Clusters", 50, 0., 50. );
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
  TH1D *hTimewalk = new TH1D("hTimewalk","",500,0,50);
  TH2D *htDiff_vs_HCALID = new TH2D("htDiff_vs_HCALID",";Channel;TDC_{HCAL}-TDC_{HODO} (ns)",288,0,288,300,-150,0);
  TH2D *haDiff_vs_HCALID = new TH2D("haDiff_vs_HCALID",";Channel;ADCtime_{HCAL}-TDC_{HODO} (ns)",288,0,288,300,0,150);
  TH2D *hTDCoffsets_vs_HCALID = new TH2D("hTDCoffsets_vs_HCALID",";Channel;TDC Offset (ns)",288,0,288,200,-50,50);
  TH2D *hADCToffsets_vs_HCALID = new TH2D("hADCToffsets_vs_HCALID",";Channel;ADC Time Offset (ns)",288,0,288,200,-50,50);
  TH2D *hTDCsig_vs_HCALID = new TH2D("hTDCsig_vs_HCALID",";Channel;TDC Std Dev (ns)",288,0,288,100,0,10);
  TH2D *hADCTsig_vs_HCALID = new TH2D("hADCTsig_vs_HCALID",";Channel;ADC Time Std Dev (ns)",288,0,288,500,0,50);
  TH2D *htDiff_vs_HCALID_corr = new TH2D("htDiff_vs_HCALID_corr",";Channel;TDC_{HCAL}-TDC_{HODO} (ns)",288,0,288,300,-150,0);
  TH2D *htDiff_vs_ADCint[kNcell+1];
  TH1D *htDiff[kNcell+1];
  TH1D *haDiff[kNcell+1];
  TH1D *haDiff_vs_col[kNcols];

  cout << endl;

  for( int i=0; i<kNcols; i++){
    haDiff_vs_col[i] = new TH1D(Form("haDiff_col%d",i),";ADCt_{HCAL}-TDC_{HODO} (ns)",300,0,150);
  }

  for( int i=0; i<kNcell+1; i++ ){
    htDiff_vs_ADCint[i] = new TH2D(Form("htDiff_vs_ADCint_bl%d",i),Form(";E_{bl%d} (GeV);TDC_{HCAL}-TDC_{HODO} (ns)",i),270,0.0,0.9,600,-250,50);
    htDiff[i] = new TH1D(Form("htDiff_bl%d",i),";TDC_{HCAL}-TDC_{HODO} (ns)",300,-150,0);
    haDiff[i] = new TH1D(Form("haDiff_bl%d",i),";ADCt_{HCAL}-TDC_{HODO} (ns)",300,0,150);

    //cout << oldTDCOffsets[i] << endl;
  }

  cout << "Variables and histograms defined." << endl;

  // Set long int to keep track of total entries
  Long64_t Nevents = C->GetEntries();

  cout << endl << "All parameters loaded and initialization complete." << endl << endl;
  cout << "Opened up tree with nentries: " << Nevents << ".." << endl << endl;
  
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
      
      C->GetEntry( nevent ); 
      
      //Proceed only if at least one track exists in BB arm
      if( BBtr_n > 0 ){
		
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
	double pp = sqrt(pow(nu,2)+2.*M_p*nu);
	double phinucleon = ephi + TMath::Pi(); //assume coplanarity
	double thetanucleon = acos( (E_e - BBtr_p[track]*cos(etheta))/pp ); //use elastic constraint on nucleon kinematics
	
	TVector3 pNhat( sin(thetanucleon)*cos(phinucleon),sin(thetanucleon)*sin(phinucleon),cos(thetanucleon));

	//Define HCal coordinate system
	TVector3 HCAL_zaxis(-sin(HCal_th),0,cos(HCal_th));
	TVector3 HCAL_xaxis(0,1,0);
	TVector3 HCAL_yaxis = HCAL_zaxis.Cross(HCAL_xaxis).Unit();
	
	TVector3 HCAL_origin = HCal_d * HCAL_zaxis + hcalheight * HCAL_xaxis;
	
	TVector3 TopRightBlockPos_DB(xoff_hcal,yoff_hcal,0);
	
	TVector3 TopRightBlockPos_Hall( hcalheight + (kNrows/2-0.5)*blockspace_hcal,
					(kNcols/2-0.5)*blockspace_hcal, 0 );

	//Sanity check on position
	hXY_HCAL_ps->Fill( HCALy, HCALx );

	//Assume that HCAL origin is at the vertical and horizontal midpoint of HCAL and locate cluster position
	HCALx += TopRightBlockPos_Hall.X() - TopRightBlockPos_DB.X();
	HCALy += TopRightBlockPos_Hall.Y() - TopRightBlockPos_DB.Y();

	//Sanity check on position after offsets applied
	hXY_HCAL->Fill( HCALy, HCALx );

	//Define intersection points for hadron vector
	double sintersect = ( HCAL_origin - vertex ).Dot( HCAL_zaxis ) / (pNhat.Dot( HCAL_zaxis ) );

	TVector3 HCAL_intersect = vertex + sintersect * pNhat;

	double yexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_yaxis );
	double xexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_xaxis );

	//Calculate the proton spot - use for cut later on
	hdxdy_HCAL->Fill( HCALy - yexpect_HCAL, HCALx - xexpect_HCAL );

	//cout << HCALy - yexpect_HCAL << " " << HCALx - xexpect_HCAL << endl;

	double E_ep = sqrt( pow(M_e,2) + pow(BBtr_p[track],2) ); // Obtain the scattered electron energy
	hE_ep->Fill( E_ep ); // Fill histogram
	
	//double p_ep = opticsCorr*BBtr_p[track]; // Obtain the magnitude of scattered electron momentum with correction factor
	double p_ep = BBtr_p[track];
	double Q2 = 2*E_e*E_ep*( 1-(BBtr_pz[track]/p_ep) ); // Obtain Q2 from beam energy, outgoing electron energy, and momenta
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

	// Coincidence timing cut and vertex cut to resolve W well
	if( fabs(diff-510)<40 && fabs(BBtr_vz[track])<0.06 ) hW_cuts->Fill( W );

	//Main cut on elastics
	if( fabs(W-W_mean)<W_sig && //Observed mean W cut
	    fabs(diff-510)<40 && //Observed coincidence trigger HCal/BB
	    BBps_e >= 0.15 && //Pion rejection in preshower
	    fabs( (BBps_e+BBsh_e)/BBtr_p[track] - 1.0 ) <= 0.3 && //Electron cut BB (E/p)
	    fabs(BBtr_vz[track]) < 0.06 //Position of track vertex cut
	    ){
	  
	  if( pow( (HCALx-xexpect_HCAL - dx0)/dx_sig,2) + pow( (HCALy-yexpect_HCAL - dy0)/dy_sig,2) <= pow(2.5,2) ){
	    //cout << "Position cut passed.." << endl;
	  }

	  hHCALe->Fill( HCALe );
	  hHODOnclus->Fill( nClusHODO );
	  if( HCALe>0.02 && nClusHODO==1 && HCALtdc[0]<-400) TDCmiss++;

	  //Get timing offsets for tdc and adc from primary block in cluster
	  if( HCALe>0.02 && HCALtdc[0]>-400 && nClusHODO<10 ){
	    htcorr_HCAL_HODO->Fill( HCALtdc[0], tHODO );
	    double HHdiff = HCALtdc[0]-tHODO;
	    double HHAdiff = HCALa[0]-tHODO;
	    htDiff_HODO_HCAL->Fill( HHdiff );
	    haDiff_HODO_HCAL->Fill( HHAdiff );
	    if(ccol>3&&ccol<8){
	      haDiff_HODO_HCAL_JLAB->Fill( HHAdiff );
	    }else{
	      haDiff_HODO_HCAL_CMU->Fill( HHAdiff );
	    }
	    
	    htDiff_vs_HCALID->Fill( cblkid[0], HHdiff );
	    haDiff_vs_HCALID->Fill( cblkid[0], HHAdiff );
	    int idx = (int)cblkid[0];
	    htDiff[ idx ]->Fill( HHdiff );
	    haDiff[ idx ]->Fill( HHAdiff );
	    if( idx<0 || idx>288 ) cout << "ERROR: TDC out of bounds at tdcID = " << idx << endl;
	    if( HHdiff<-65. && HHdiff>-80. ) htDiff_vs_ADCint[idx]->Fill( cblke[0], HHdiff );
	    if( HHdiff<-65. && HHdiff>-80. ) htDiff_vs_HCALID_corr->Fill( cblkid[0], HHdiff+22.0*cblke[0] ); //Timewalk corrected from linear fit to time diff vs clus e
	    
	    int col = (int)ccol;
	    haDiff_vs_col[col]->Fill( HHAdiff );

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
    }
  }

  //Fit timing data and extract offsets
  TF1 *f01;
  TF1 *f02;
  double JLAB_mean = 0.;
  double CMU_mean = 0.;

  haDiff_HODO_HCAL_JLAB->Fit("gaus","Q","",20,80);
  f01=haDiff_HODO_HCAL_JLAB->GetFunction("gaus");
  JLAB_mean = f01->GetParameter(1);
  cout << "JLAB mean =" << JLAB_mean << "." << endl;

  haDiff_HODO_HCAL_CMU->Fit("gaus","Q","",20,80);
  f02=haDiff_HODO_HCAL_CMU->GetFunction("gaus");
  CMU_mean = f02->GetParameter(1);
  cout << "CMU mean =" << CMU_mean << "." << endl;

  for(int i=0; i<kNcell; i++){
    TF1 *f1;
    TF1 *f2;
    TF1 *f3;
    if( htDiff[i]->GetEntries()>tFitMin ){
      htDiff[i]->Fit("gaus","Q","",-80,-65);
      f1=htDiff[i]->GetFunction("gaus");
      TDCoffsets[i] = -75. - f1->GetParameter(1);
      TDCsig[i] = f1->GetParameter(2);
      //cout << f1->GetParameter(1) << endl;
    }
    if( haDiff[i]->GetEntries()>tFitMin ){
      haDiff[i]->Fit("gaus","Q","",20,80);
      f2=haDiff[i]->GetFunction("gaus");
      ADCToffsets[i] = 50. - f2->GetParameter(1);
      ADCTsig[i] = f2->GetParameter(2);
      //cout << f2->GetParameter(1) << endl;
    }
    if( htDiff_vs_ADCint[i]->GetEntries()>tFitMin ){
      htDiff_vs_ADCint[i]->Fit("pol1","Q","",0.0,0.2);
      f3=htDiff_vs_ADCint[i]->GetFunction("pol1");
      TvsE[i]=f3->GetParameter(1);
      hTimewalk->Fill(f3->GetParameter(1));
    }
  }

  for(int i=0; i<kNcell; i++){
    if( htDiff[i]->GetEntries()>tFitMin ){
      hTDCoffsets_vs_HCALID->Fill(i,TDCoffsets[i]);
      hTDCsig_vs_HCALID->Fill(i,TDCsig[i]);
    }
    if( haDiff[i]->GetEntries()>tFitMin ){
      hADCToffsets_vs_HCALID->Fill(i,ADCToffsets[i]);
      hADCTsig_vs_HCALID->Fill(i,ADCTsig[i]);
    }
  }
 
  //Write out diagnostic histos and print to console
  fout->Write();

  cout << endl << "TDC Offsets: " << endl << endl;

  int cell = 0;
  for( int r = 0; r<kNrows; r++){
    for( int c = 0; c<kNcols; c++){
      if( htDiff[cell]->GetEntries()<tFitMin ){
	cout << "--" << "  ";
      }else{
	cout << TDCoffsets[cell] + oldTDCOffsets[cell]/TDCCalib << "  ";
      }
      cell++;
    }
    cout << endl;
  }

  cout << endl << "ADC Time Offsets: " << endl << endl;

  cell = 0;

  for( int r = 0; r<kNrows; r++){
    for( int c = 0; c<kNcols; c++){
      if( haDiff[cell]->GetEntries()<tFitMin ){
	if(c>3&&c<8){
	  cout << 50. - JLAB_mean << "  ";
	}else{
	  cout << 50. - CMU_mean << "  ";
	}
      }else{
	cout << ADCToffsets[cell] << "  ";
      }
      cell++;
    }
    cout << endl;
  }

  cout << endl << "Timewalk Corrections: " << endl << endl;

  cell = 0;

  for( int r = 0; r<kNrows; r++){
    for( int c = 0; c<kNcols; c++){
      if( htDiff[cell]->GetEntries()<tFitMin ){
	cout << hTimewalk->GetMean() << "  ";
      }else{
	cout << TvsE[cell] << "  ";
      }
      cell++;
    }
    cout << endl << endl;
  }

  //Declare outfile
  ofstream params;
  params.open( paramsPath );
  params << "#Timewalk parameters from SBS-" << kine << " obtained " << date.c_str() << endl;
  params << "#HCal_timewalk = " << endl;

  cell = 0;
  for( int r=0; r<kNrows; r++ ){
    for( int c=0; c<kNcols; c++ ){
      if( htDiff[cell]->GetEntries()<tFitMin ){
	params << hTimewalk->GetMean() << "  ";
      }else{
	params << TvsE[cell] << "  ";
      }
      cell++;
    }
    params << endl;
  }

  params << "#" << endl << "#" << endl << "#HCal_TDC_offset = " << endl;
  
  cell = 0;
  for( int r=0; r<kNrows; r++ ){
    for( int c=0; c<kNcols; c++ ){
      if( htDiff[cell]->GetEntries()<tFitMin ){
	params << 0.0 << "  ";
      }else{
	params << htDiff[cell] << "  ";
      }
      cell++;
    }
    params << endl;
  }

  params << "#" << endl << "#" << endl << "#HCal_ADCTime_offset = " << endl;

  cell = 0;
  for( int r=0; r<kNrows; r++ ){
    for( int c=0; c<kNcols; c++ ){
      if( htDiff[cell]->GetEntries()<tFitMin ){
	if(ccol>3&&ccol<8){
	  cout << 50. - JLAB_mean << "  ";
	}else{
	  cout << 50. - CMU_mean << "  ";
	}
      }else{
	params << haDiff[cell] << "  ";
      }
      cell++;
    }
    params << endl;
  }

  params.close();

  cout << "Elastic yield for analyzed runs: " << elasYield << ". Total events analyzed: " << Nevents << ". Total TDC misses: " << TDCmiss << "." << endl << endl;

  cout << "Timing offset analysis complete and parameters written to file." << endl;

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
