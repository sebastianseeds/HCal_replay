//SSeeds 1.27.23 Script written to investigate the hodoscope/hcal tdc times to mine out the RF signature from the accelerator

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

const Int_t kNcell = 288; // Total number of HCal modules
const Int_t kNrows = 24; // Total number of HCal rows
const Int_t kNcols = 12; // Total number of HCal columns
const Int_t maxTdcChan = 10; // Set to accomodate original 5 TDCTrig channels with buffer
const Int_t maxTracks = 1000; // Reasonable limit on tracks to be stored per event
const Double_t hcalheight = 0.365; //m The height of the center of HCAL above beam
const Double_t scint_gamma_v = 12.61; //estimated speed of light in plastic scintillator cm/ns 
const Double_t PI = TMath::Pi();
const Double_t M_e = 0.00051;
const Double_t M_p = 0.938272;
const Double_t hIDvtrXoff = 46.99; //bb.hodotdc.clus.id:bb.tr.x pol1 fit p0
const Double_t hIDvtrXm = -60.24; //bb.hodotdc.clus.id:bb.tr.x pol1 fit p1
const Double_t hIDvtrXsig = 0.02; //bb.tr.x, bb.hodotdc.clus.id==hIDvtrXoff, 2sig
const Double_t htrXshift = 0.02; //empyrical check on bb.tr.x fit dists

void rf_study( Int_t kine = 4 ){

  // Declare Chain for many root files
  TChain *C = new TChain("T");

  // Declare general physics parameters to be modified by input config file
  //Double_t kine = 4; // Keep track of kinematic calibrating from
  Double_t tFitMin = 30; // Minimum number of entries per channel to calibrate ADC/TDC time
  Double_t t_trig = 510; // Mean tdc trig value (HCAL - BB) 
  Double_t E_e = 1.92; // Energy of beam (incoming electrons from accelerator)
  Double_t HCal_d = 14.5; // Distance to HCal from scattering chamber for comm1
  Double_t HCal_th = 35.0; // Angle that the center of HCal is at  
  Double_t W2_mean = 0.93; // Mean of W at current kinematic
  Double_t W2_sig = 0.039; // Width of W at current kinematic
  // Additional global variables
  Long64_t elasYield = 0; // Keep track of total elastics analyzed
  Long64_t clusBarMiss = 0; // Keep track of total elastics analyzed
  Long64_t hodoClusYield = 0; // Keep track of total elastics analyzed

  Int_t TDCmiss = 0; // Keep track of TDC misses

  cout << endl;

  string configfilename = Form("config/srf_study_sbs%d.cfg",kine);

  // Reading config file
  ifstream configfile(configfilename);
  TString currentline;
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      if(!currentline) cout << "WARNING: No file exists at " << currentline << "." << endl;
      C->Add(currentline);
      cout << "Loaded file at: " << currentline << endl;
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
    Int_t ntokens = tokens->GetEntries();
    if( ntokens>1 ){
      TString skey = ( (TObjString*)(*tokens)[0] )->GetString();
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
      if( skey == "W2_mean" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        W2_mean = sval.Atof();
	cout << "Loading W2 mean cut: " << W2_mean << endl;
      }
      if( skey == "W2_sig" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        W2_sig = sval.Atof();
	cout << "Loading W2 sigma cut: " << W2_sig << endl;
      }
    }
    delete tokens;
  }
  
  TEventList *elist = new TEventList("elist","Elastic Event List");
  C->Draw(">>elist",globalcut);

  cout << endl << "Event list populated with cut placed on elastics to constrain TOF." << endl;

  // Declare general detector and physics parameters
  Double_t TDCT_id[maxTdcChan], TDCT_tdc[maxTdcChan]; 
  Int_t TDCTndata;
  Int_t HODOndata;
  Int_t HODOndata_Rle;
  Int_t HODOndata_Rtot;
  Int_t HODOndata_id;
  Int_t HODOndata_cid;
  Int_t HODOndata_tmean;
  Int_t HODOndata_ymean;
  UInt_t runI=1; 
  UInt_t runN=1;
  UInt_t TBits;
  ULong64_t runT;
  Double_t kineW2;

  // BB params
  Double_t BBtr_p[maxTracks], BBtr_px[maxTracks], BBtr_py[maxTracks], BBtr_pz[maxTracks];
  Double_t BBtr_x[maxTracks], BBtr_y[maxTracks];
  Double_t BBtr_vz[maxTracks], BBtr_chi2[maxTracks];
  Double_t BBtr_n, BBps_x, BBps_y, BBps_e, BBsh_x, BBsh_y, BBsh_e;

  // HCal params
  Double_t HCALx, HCALy, HCALe;
  Double_t HCALblktdc[kNcell], HCALblka[kNcell];
  Double_t HCALtdc[kNcell], HCALatime[kNcell], HCALa_p[kNcell];
  Double_t crow, ccol, nblk;
  Double_t cblkid[kNcell], cblke[kNcell];
  Double_t ce[kNcell], catime[kNcell], cid[kNcell], cnblk[kNcell], ctdc[kNcell];
  Double_t nclus;

  // Hodo params
  Double_t HODOtmean[kNcell];
  Double_t HODOnclus;
  Double_t HODOtdc[kNcell];
  Double_t HODOtdcid[kNcell];
  Double_t HODOtdc_tot[kNcell];
  Double_t HODOtdc_mult[kNcell];
  Double_t HODOymean[kNcell];
  Double_t HODOcid[kNcell];
  Double_t HODObarid[kNcell];
  Double_t HODORle[kNcell];
  Double_t HODORtot[kNcell];
  
  // Declare root tree variables and set values to memory locations in root file
  // Turn off all branches
  C->SetBranchStatus( "*", 0 );
  // Turn on specific branches for analysis
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
  C->SetBranchStatus( "sbs.hcal.clus.e", 1 );
  C->SetBranchStatus( "sbs.hcal.clus.atime", 1 );
  C->SetBranchStatus( "sbs.hcal.clus.id", 1 );
  C->SetBranchStatus( "sbs.hcal.clus.nblk", 1 );
  C->SetBranchStatus( "sbs.hcal.clus.tdctime", 1 );
  C->SetBranchStatus( "sbs.hcal.nclus", 1 );  
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
  C->SetBranchStatus( "bb.hodotdc.clus.id", 1 );
  C->SetBranchStatus( "bb.hodotdc.nclus", 1 );
  C->SetBranchStatus( "bb.hodotdc.tdc", 1 );
  C->SetBranchStatus( "bb.hodotdc.tdc_tot", 1 );
  C->SetBranchStatus( "bb.hodotdc.tdcelemID", 1 );
  C->SetBranchStatus( "bb.hodotdc.clus.ymean", 1 );
  C->SetBranchStatus( "bb.hodotdc.bar.tdc.R.leW", 1 );
  C->SetBranchStatus( "bb.hodotdc.bar.tdc.R.tot", 1 );
  C->SetBranchStatus( "bb.hodotdc.bar.tdc.id", 1 );

  C->SetBranchStatus( "Ndata.bb.hodotdc.clus.id", 1 );
  C->SetBranchStatus( "Ndata.bb.hodotdc.clus.ymean", 1 );
  C->SetBranchStatus( "Ndata.bb.hodotdc.clus.tmean", 1 );
  C->SetBranchStatus( "Ndata.bb.hodotdc.bar.tdc.R.leW", 1 );
  C->SetBranchStatus( "Ndata.bb.hodotdc.tdc", 1 );
  C->SetBranchStatus( "e.kine.W2", 1 );

  // Map the branches to variables
  C->SetBranchAddress( "sbs.hcal.x", &HCALx );
  C->SetBranchAddress( "sbs.hcal.y", &HCALy );
  C->SetBranchAddress( "sbs.hcal.e", &HCALe );
  C->SetBranchAddress( "sbs.hcal.rowblk", &crow );
  C->SetBranchAddress( "sbs.hcal.colblk", &ccol );
  C->SetBranchAddress( "sbs.hcal.nblk", &nblk );
  C->SetBranchAddress( "sbs.hcal.clus_blk.id", cblkid );
  C->SetBranchAddress( "sbs.hcal.clus_blk.e", cblke );
  C->SetBranchAddress( "sbs.hcal.clus_blk.tdctime", HCALblktdc );
  C->SetBranchAddress( "sbs.hcal.clus_blk.atime", HCALblka );
  C->SetBranchAddress( "sbs.hcal.clus.e", ce );
  C->SetBranchAddress( "sbs.hcal.clus.atime", catime );
  C->SetBranchAddress( "sbs.hcal.clus.id", cid );
  C->SetBranchAddress( "sbs.hcal.clus.nblk", cnblk );
  C->SetBranchAddress( "sbs.hcal.clus.tdctime", ctdc );
  C->SetBranchAddress( "sbs.hcal.nclus", &nclus );
  C->SetBranchAddress( "bb.tr.chi2", BBtr_chi2 );
  C->SetBranchAddress( "bb.tr.n", &BBtr_n );
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
  C->SetBranchAddress( "bb.hodotdc.clus.tmean", HODOtmean );
  C->SetBranchAddress( "Ndata.bb.hodotdc.clus.tmean", &HODOndata_tmean );
  C->SetBranchAddress( "bb.hodotdc.clus.id", HODOcid );
  C->SetBranchAddress( "Ndata.bb.hodotdc.clus.id", &HODOndata_cid );
  C->SetBranchAddress( "bb.hodotdc.clus.ymean", HODOymean );
  C->SetBranchAddress( "Ndata.bb.hodotdc.clus.ymean", &HODOndata_ymean );
  C->SetBranchAddress( "bb.hodotdc.bar.tdc.R.le", HODORle );
  C->SetBranchAddress( "bb.hodotdc.bar.tdc.R.tot", HODORtot );
  C->SetBranchAddress( "bb.hodotdc.bar.tdc.id", HODObarid ); //Ndata same as tdc.R.leW
  C->SetBranchAddress( "Ndata.bb.hodotdc.bar.tdc.R.le", &HODOndata_Rle );
  C->SetBranchAddress( "Ndata.bb.hodotdc.bar.tdc.R.tot", &HODOndata_Rtot );
  C->SetBranchAddress( "Ndata.bb.hodotdc.bar.tdc.id", &HODOndata_id );
  C->SetBranchAddress( "bb.hodotdc.nclus", &HODOnclus );
  C->SetBranchAddress( "bb.hodotdc.tdc", HODOtdc );
  C->SetBranchAddress( "bb.hodotdc.tdc_tot", HODOtdc_tot );
  C->SetBranchAddress( "bb.hodotdc.tdc_mult", HODOtdc_mult );
  C->SetBranchAddress( "bb.hodotdc.tdcelemID", HODOtdcid );
  C->SetBranchAddress( "Ndata.bb.hodotdc.tdc", &HODOndata );
  C->SetBranchAddress( "e.kine.W2", &kineW2 );

  cout << "Tree variables linked." << endl;
  
  // Declare outfile
  TFile *fout = new TFile( Form("output/rfStudyOUT_sbs%d.root",kine), "RECREATE" );

  // // Physics and trigger histograms
  // TH1D *hW2 = new TH1D( "hW2", "W2; GeV", 250, 0, 6 );
  // TH1D *hW2_cut = new TH1D( "hW2_cut", "W2 Cut; GeV", 250, 0, 6 );
  // TH1D *hL1A = new TH1D( "hL1A", "BBCal HI Time (L1A); ns", 200, 250, 450 );
  // TH1D *hL1A_hodocorr = new TH1D( "hL1A_hodocorr", "BBCal HI Time Hodoscope Corrected (L1A); ns", 200, 250, 450 );
  // TH1D *hRFsig = new TH1D( "hRFsig", "L1A Corrected - RF time; ns", 100, 340, 390 );
  // TH1D *hTRF = new TH1D( "hTRF","RF Signature (ns)", 1700, -10, 160 );
  // TH1D *hTRFmod = new TH1D( "hTRFmod","RF Signature Modulo 4ns (ns)", 200, -10, 10 );
  // TH1D *hTHCALvRF = new TH1D( "hTHCALvRF","HCal time - RF Signature Modulo 4ns (ns)", 1650, -450, 1200 );
  // TH1D *hHODOnclus = new TH1D( "hHODOnclus","Number of Hodoscope Clusters", 50, 0., 50. );
  // TH2D *hHodoVsTrackP = new TH2D( "hHodoVsTrackP", "Hodo Mean Time vs Elastic Track Mom", 60, 1.8, 2.4, 60, -6, 12 );
  // TH1D *hHODOtdc = new TH1D( "hHODOtdc", "Hodoscope ALL TDC", 2000, -500, 1500);
  // TH1D *hHODOtmean = new TH1D( "hHODOtmean", "Hodoscope MEAN TDC", 400, -20, 20);
  // TH1D *hHODOtmean_RFmod = new TH1D( "hHODOtmean_RFmod", "Hodoscope MEAN TDC - RFmod4", 400, -20, 20);
  // TH1D *hHODORle = new TH1D( "hHODORle", "Hodoscope PClus Right Leading Edge", 400, -20, 20);
  // TH1D *hHODORle_trYslice = new TH1D( "hHODORle_trYslice", "Hodoscope PClus Right Leading Edge, yslice, bar var", 400, -20, 20);
  // TH1D *hHODORle_trYslice_RFcorr = new TH1D( "hHODORle_trYslice_RFcorr", "Hodoscope PClus Right Leading Edge, yslice, rfcorr, bar var", 400, -20, 20);  
  // TH1D *hHODORle_trYtotslice = new TH1D( "hHODORle_trYtotslice", "Hodoscope PClus Right Leading Edge, yslice, totslice, bar var", 400, -20, 20);
  // TH1D *hHODORle_gcorr = new TH1D( "hHODORle_gcorr", "Hodoscope PClus Right Leading Edge Light V Corrected", 600, -20, 40);
  // TH1D *htrY = new TH1D( "htrY", "Track Y; m", 100, -0.5, 0.5 );
  // TH1D *hhodoY = new TH1D( "hhodoY", "Hodoscope Cluster Y mean; m", 100, -0.5, 0.5 );
  // TH1D *htot = new TH1D( "htot", "Hodo Right PMT time over threshold; ns", 120, 0, 40 );
  // TH1D *htot_B2 = new TH1D( "htot_B2", "Hodo Right PMT time over threshold, cluster block 2; ns", 120, 0, 40 );
  // TH1D *htot_B3 = new TH1D( "htot_B3", "Hodo Right PMT time over threshold, cluster block 3; ns", 120, 0, 40 );
  // TH1D *hB1B2diff_trYslice = new TH1D( "hB1B2diff_trYslice", "Hodoscope PClus Right Leading Edge, B1-B2, yslice, bar var", 400, -20, 20);
  // TH1D *hB1B2diff_trYslice_RFcorr = new TH1D( "hB1B2diff_trYslice_RFcorr", "Hodoscope PClus Right Leading Edge, B1-B2, yslice, bar var, RF corr", 400, -20, 20);
  // TH1D *hB1B3diff_trYslice = new TH1D( "hB1B3diff_trYslice", "Hodoscope PClus Right Leading Edge, B1-B3, yslice, bar var", 400, -20, 20);
  // TH1D *hB1B3diff_trYslice_RFcorr = new TH1D( "hB1B3diff_trYslice_RFcorr", "Hodoscope PClus Right Leading Edge, B1-B3, yslice, bar var, RF corr", 400, -20, 20);
  // TH1D *hB1B2diff_trYtotslice = new TH1D( "hB1B2diff_trYtotslice", "Hodoscope PClus Right Leading Edge, B1-B2, yslice, totslice, bar var", 400, -20, 20);
  // TH1D *hB1B2diff_trYtotslice_RFcorr = new TH1D( "hB1B2diff_trYtotslice_RFcorr", "Hodoscope PClus Right Leading Edge, B1-B2, yslice, totslice, bar var, RF corr", 400, -20, 20);  
  // TH1D *hB1B3diff_trYtotslice = new TH1D( "hB1B3diff_trYtotslice", "Hodoscope PClus Right Leading Edge, B1-B3, yslice, totslice, bar var", 400, -20, 20); 
  // TH1D *hB1B3diff_trYtotslice_RFcorr = new TH1D( "hB1B3diff_trYtotslice_RFcorr", "Hodoscope PClus Right Leading Edge, B1-B3, yslice, totslice, bar var, RF corr", 400, -20, 20);  TH2D *htrXvHID = new TH2D("htrXvHID", "Track X vs Hodo bar ID", 120, -0.6, 0.6, 90, 0, 90);
  // TH1D *hFIT1 = new TH1D("hFIT1", "Fit map from track to HODO cid", 90, 0, 90);
  // TH1D *hFIT2 = new TH1D("hFIT2", "Fit map from HODO cid to track", 120, -0.6, 0.6);
  
  // TH1D *hB1B2diff_hcal = new TH1D( "hB1B2diff_hcal", "HCal PClus TDC, B1-B2", 400, -20, 20);
  // TH1D *hB1B2diff_hcal_RFmod = new TH1D( "hB1B2diff_hcal_RFmod", "HCal PClus TDC, B1-B2, RFmod", 400, -20, 20);

  TH1D *hrawmaxtot = new TH1D( "hrawmaxtot", "Hodoscope raw signal max time over threshold", 400, 0, 100 );
  TH1D *hrawB1B2 = new TH1D( "hrawB1B2", "Hodoscope raw signal high tot block - high tot block - 1", 400, 0, 100 ); 
  TH1D *hrawB1B3 = new TH1D( "hrawB1B3", "Hodoscope raw signal high tot block - high tot block - 2", 400, 0, 100 ); 
  TH1D *hrawB1B4 = new TH1D( "hrawB1B4", "Hodoscope raw signal high tot block - high tot block - 3", 400, 0, 100 ); 
  TH1D *hrawB1B5 = new TH1D( "hrawB1B5", "Hodoscope raw signal high tot block - high tot block - 4", 400, 0, 100 );  


  cout << "Variables and histograms defined." << endl;

  // Set long int to keep track of total entries
  Long64_t Nevents = elist->GetN();
  UInt_t run_number = 0;

  cout << "Opened up TChain with nentries: " << C->GetEntries() << ". After globalcut: " << Nevents << "." << endl << endl;

  //Loop over events
  cout << "Main loop over all data commencing.." << endl;
  Double_t progress = 0.;
  Double_t timekeeper = 0., timeremains = 0.;

  for(Long64_t nevent = 0; nevent<Nevents; nevent++){
    
    if( nevent == Nevents ){
      cout << "Hit limit on event " << nevent << endl;
      break;
    }
    
    cout << nevent << "/" << Nevents << " \r";
    cout.flush();
      
    C->GetEntry( elist->GetEntry( nevent ) ); 
      
    //Basic checks, only global/W2 cuts.
    
    bool trn1 = BBtr_n==1;
    if( !trn1 ) continue;

    bool elastic = fabs( kineW2-W2_mean )<W2_sig;
    if( !elastic ) continue;

    //Add slice on bb.tr.y
    bool inY = fabs(BBtr_y[0])<0.02;
    if( !inY ) continue;

    Double_t RF_time = -1000;
    bool rfmult = false;
    Int_t rfnhits = 0;
    for(Int_t ihit=0; ihit<TDCTndata; ihit++){
      if(rfmult) cout << "Multiple hits in RF time channel, up to: " << rfnhits << endl;
      if(TDCT_id[ihit]==4){
	RF_time=TDCT_tdc[ihit];
	rfmult=true;
	rfnhits++;
      }
    }

    //Obtain block with highest ToT
    bool multskip = false;
    Double_t maxtot = 0.;
    Int_t maxtotid = -1;
    Double_t maxtottdc = 0.;
    Double_t maxtottdcB5 = 0.;
    Double_t maxtottdcB4 = 0.;
    Double_t maxtottdcB3 = 0.;
    Double_t maxtottdcB2 = 0.;
    for( Int_t h=1; h<HODOndata; h++ ){
      Int_t blkid = HODOtdcid[h];
      Int_t blkmult = HODOtdc_mult[h];
      Double_t blktdc = HODOtdc[h];
      Double_t blktot = HODOtdc_tot[h];

      if( blkmult!=1 ) continue;
      if( blktot>maxtot ){
	maxtot = blktot;
	maxtotid = blkid;
	maxtottdc = blktdc;
      }

    }
    
    //Fill tdc difference histos with high ToT as seed
    for( Int_t h=1; h<HODOndata; h++ ){
      Int_t blkid = HODOtdcid[h];
      Int_t blkmult = HODOtdc_mult[h];
      Double_t blktdc = HODOtdc[h];
      Double_t blktot = HODOtdc_tot[h];

      if( (blkid-4)==maxtotid ) maxtottdcB5=blktdc;
      if( (blkid-3)==maxtotid ) maxtottdcB4=blktdc;
      if( (blkid-2)==maxtotid ) maxtottdcB3=blktdc;
      if( (blkid-1)==maxtotid ) maxtottdcB2=blktdc;

    }

    hrawmaxtot->Fill(maxtot);

    if( maxtottdcB2>0 ) hrawB1B2->Fill(maxtottdc-maxtottdcB2);
    if( maxtottdcB3>0 ) hrawB1B3->Fill(maxtottdc-maxtottdcB3);
    if( maxtottdcB4>0 ) hrawB1B4->Fill(maxtottdc-maxtottdcB4);
    if( maxtottdcB5>0 ) hrawB1B5->Fill(maxtottdc-maxtottdcB5);


    // //Continue if no clusters exist in Hodoscope
    // if( HODOndata_tmean==0 ) continue;
    // hodoClusYield++;

    //BEGIN coincidence timing cut
    //Cut on BBCal and HCal trigger coincidence and plot all other trigger timing. All ref to tdcs need index adjustment (ihit+1).
    
    // Double_t bbcal_time=0., hcal_time=0., RF_time=0., coin_time=0.; 
    // Double_t bbcalLO_time=0., bbcalHI_time=0., bbcalHIveto_time=0., edtm_time=0.;
    //GEn
    /*
    for(Int_t ihit=0; ihit<TDCTndata; ihit++){
      if(TDCT_id[ihit]==5) bbcal_time=TDCT_tdc[ihit];
      if(TDCT_id[ihit]==4) RF_time=TDCT_tdc[ihit];
      if(TDCT_id[ihit]==3) edtm_time=TDCT_tdc[ihit];
      if(TDCT_id[ihit]==2) bbcalHI_time=TDCT_tdc[ihit];
      if(TDCT_id[ihit]==1) coin_time=TDCT_tdc[ihit];
      if(TDCT_id[ihit]==0) hcal_time=TDCT_tdc[ihit];
    }
    */
    // //GMn
    // for(Int_t ihit=0; ihit<TDCTndata; ihit++){
    //   if(TDCT_id[ihit]==5) bbcal_time=TDCT_tdc[ihit];
    //   if(TDCT_id[ihit]==4) RF_time=TDCT_tdc[ihit];
    //   if(TDCT_id[ihit]==3) edtm_time=TDCT_tdc[ihit];
    //   if(TDCT_id[ihit]==2) bbcalLO_time=TDCT_tdc[ihit];
    //   if(TDCT_id[ihit]==1) bbcalHIveto_time=TDCT_tdc[ihit];
    //   if(TDCT_id[ihit]==0) hcal_time=TDCT_tdc[ihit];
    // }
    // Double_t diff = hcal_time - bbcal_time; 
    // //Double_t RFmod = std::fmod(RF_time,2); //RF Signature measures the bunch crossing of beam electrons at 4 or 2 ns intervals
    // //Double_t RFmod = HODOtmean[0] -(RF_time-bbcal_time);
    // Double_t RFmod = std::fmod(RF_time,4); //Confirmed 4ns bunch crossing time GMn
    // hTRF->Fill(RF_time);
    // hTRFmod->Fill(RFmod);
    // hTHCALvRF->Fill( hcal_time - RFmod );

    // Double_t L1A = bbcal_time;
    // Double_t L1A_hodocorr = bbcal_time - HODOtmean[0];

    // Double_t RFsig = L1A_hodocorr - RF_time;
  
    //cout << L1A_hodocorr << endl;

    //if( fabs( diff-t_trig )>30 ) continue;
    //cout << "hcal_time:" << hcal_time << ", BBCal time:" << bbcal_time << endl;

    //END coincidence timing cut

    // //BEGIN elastic projections/cuts
    // Double_t etheta = acos( BBtr_pz[0]/BBtr_p[0] );
    // Double_t ephi = atan2( BBtr_py[0], BBtr_px[0] );
    // TVector3 vertex(0,0,BBtr_vz[0]); // z location of vertex in hall coordinates
    // TLorentzVector Pbeam(0,0,E_e,E_e); //Mass of e negligable
    // TLorentzVector kprime(BBtr_px[0],BBtr_py[0],BBtr_pz[0],BBtr_p[0]);
    // TLorentzVector Ptarg(0,0,0,M_p);
    // TLorentzVector q = Pbeam - kprime;
    // TLorentzVector PgammaN = Ptarg + q; //(-px, -py, ebeam - pz, Mp + ebeam - p)
    // Double_t pel = E_e/(1.+E_e/M_p*(1.-cos(etheta)));
    // Double_t nu = E_e - BBtr_p[0];
    // Double_t pp = sqrt(pow(nu,2)+2.*M_p*nu); //momentum of the proton
    // Double_t phinucleon = ephi + PI; //assume coplanarity
    // Double_t thetanucleon = acos( (E_e - BBtr_p[0]*cos(etheta))/pp ); //use elastic constraint on nucleon kinematics
    // TVector3 pNhat( sin(thetanucleon)*cos(phinucleon),sin(thetanucleon)*sin(phinucleon),cos(thetanucleon));

    // //Define HCal coordinate system
    // TVector3 HCAL_zaxis(sin(-HCal_th),0,cos(-HCal_th));
    // TVector3 HCAL_xaxis(0,-1,0);
    // TVector3 HCAL_yaxis = HCAL_zaxis.Cross(HCAL_xaxis).Unit();	
    // TVector3 HCAL_origin = HCal_d * HCAL_zaxis + hcalheight * HCAL_xaxis;

    // //Define intersection points for hadron vector
    // Double_t sintersect = ( HCAL_origin - vertex ).Dot( HCAL_zaxis ) / ( pNhat.Dot( HCAL_zaxis ) );
    // TVector3 HCAL_intersect = vertex + sintersect * pNhat;
    // Double_t yexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_yaxis );
    // Double_t xexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_xaxis );
    // Double_t E_ep = sqrt( pow(M_e,2) + pow(BBtr_p[0],2) ); // Obtain the scattered electron energy
    // Double_t p_ep = BBtr_p[0];
    // Double_t Q2 = 2*E_e*E_ep*( 1-(BBtr_pz[0]/p_ep) ); // Obtain Q2 from beam energy, outgoing electron energy, and momenta
    // //Use the electron kinematics to predict the proton momentum assuming elastic scattering on free proton at rest:
    // Double_t E_pp = nu+M_p; // Get total energy of the proton
    // //Double_t Enucleon = sqrt(pow(pp,2)+pow(M_p,2)); // Check on E_pp, same
    // Double_t KE_p = nu; //For elastics
 
    // Double_t xDiff = HCALx - xexpect_HCAL;
    // Double_t yDiff = HCALy - yexpect_HCAL;

    //hW2->Fill(kineW2);
    //hHODOnclus->Fill(HODOnclus);

    // //Check how often the tdc failed to register for an otherwise good event
    // if( HCALe>0.02 && HCALblktdc[0]<-400 && HODOnclus<10 ) TDCmiss++;

    // //track x hodo cluster primary bar id cut bool
    // Double_t hodofittrX = (BBtr_x[0]-htrXshift)*hIDvtrXm + hIDvtrXoff;
    // hFIT1->Fill(hodofittrX);
    // Double_t trXfithodo = (HODOcid[0]-hIDvtrXoff)/hIDvtrXm + htrXshift;
    // hFIT2->Fill(trXfithodo);
    // //bool inXtrack = (HODOcid[0] > hodofittrX - hIDvtrXsig) && (HODOcid[0] < hodofittrX + hIDvtrXsig); // HODOcid == BBtr_x[0]*hIDvtrXm + hIDvtrXoff with hIDvtrXsig
    // bool inXtrack = (BBtr_x[0] > trXfithodo - 1.) && (BBtr_x[0] < trXfithodo + 1.); // HODOcid == BBtr_x[0]*hIDvtrXm + hIDvtrXoff with hIDvtrXsig
    
    //ELASTIC CUT
    //bool elastic = fabs( kineW2-W2_mean )<W2_sig;

    //bool elastic = fabs( kineW2-W2_mean )<W2_sig && fabs( diff-t_trig )<30;
    //bool elastic = fabs( diff-t_trig )<30 && HCALblka[0]<70 && HCALblka[0]>50;
    //if( elastic==false ) continue;
    // elasYield++;
    // hW2_cut->Fill(kineW2);

    // hL1A->Fill(L1A);
    // hL1A_hodocorr->Fill(L1A_hodocorr);
    // if(RF_time>0.) hRFsig->Fill(RFsig);

    // for( Int_t i=0; i<HODOndata; i++ ){
    //   //cout << HODOtdc[i] << endl;
    //   hHODOtdc->Fill(HODOtdc[i]);
    // }

    // htrXvHID->Fill( BBtr_x[0], HODOcid[0]);

    // if( inXtrack==false ) continue;

    // Double_t hodoRcorr = 0.;
    // Double_t Rle_time = -1000.;
    // Double_t Rle_trYtime = -1000.;
    // Double_t Rle_trYtime_B2 = -1000.;
    // Double_t Rle_trYtime_B3 = -1000.;
    // Double_t Rle_tot = -1000.;
    // Double_t Rle_tot_B2 = -1000.;
    // Double_t Rle_tot_B3 = -1000.;
    // Double_t gcorr = 0.;
    // //Select the right element of the bar variables using the cluster id
    // if( HODOndata_cid>0 && fabs(HODOymean[0])<3.0){
    //   for( Int_t i=1; i<HODOndata_Rle; i++ ){
    // 	//cout << HODOcid[0] << " " << HODObarid[i] << ", ";
    // 	if( HODOcid[0]==HODObarid[i] ){ //Apparent mismatch between cluster and bar indexes much to my chagrin
    // 	  Rle_time = HODORle[i];
    // 	  //cout << "HODOcid[0]:" << HODOcid[0] << " " << "HODObarid[i]:" << HODObarid[i] << endl;
    // 	  //cout << "HODOcid[1]:" << HODOcid[1] << " " << "HODObarid[i]:" << HODObarid[i] << endl;
    // 	}
    //   }
    //   //cout << endl;
    //   if( Rle_time==-1000. ){
    // 	//cout << "DERP: Bar/Cluster Mismatch" << endl;
    // 	clusBarMiss++;
    // 	continue;
    //   }

    //   gcorr = (10.*HODOymean[0]-3.0)/scint_gamma_v; //convert to cm, then get corrected time with light velocity in scintillator and non-dispersive position in hodoscope paddle

    //   //cout << gcorr << endl;

    //   hodoRcorr = Rle_time-gcorr;
    //   hHODORle_gcorr->Fill(hodoRcorr);
    //   hHODORle->Fill(Rle_time);

    // }
    
    // hHODOtmean->Fill(HODOtmean[0]);
    // hHODOtmean_RFmod->Fill(HODOtmean[0]-RFmod);

    // //cout << BBtr_y[0] << " " << HODOndata_cid << " " << endl;
    // htrY->Fill(BBtr_y[0]);
    // hhodoY->Fill(HODOymean[0]);

    // //cout << endl << HODORtot[0] << endl;

    // //select information on bar which matches the id of the primary cluster element
    // // if( HODOndata_cid>0 && fabs(BBtr_y[0])<0.02){
    // //   for( Int_t i=1; i<HODOndata_Rle; i++ ){
    // // 	if( HODOcid[0]==HODObarid[i] ){
    // // 	  Rle_trYtime = HODORle[i];
    // // 	  Rle_tot = HODORtot[i];
    // // 	}
    // //   }
    // // }

    // if( HODOndata_cid>0 && fabs(BBtr_y[0])<0.02){
    //   for( Int_t i=1; i<HODOndata_Rle; i++ ){
    // 	if( HODObarid[i]==HODOcid[0] ){ //if the bar id matches the cluster first element id
    // 	  Rle_trYtime = HODORle[i];
    // 	  Rle_tot = HODORtot[i];

    // 	  //cout << Rle_trYtime << endl;
    // 	  //cout << HODObarid[i] << " " << HODObarid[i+1] << endl;

    // 	  if( HODObarid[i]-HODObarid[i+1]==-1 ) Rle_trYtime_B2 = HODORle[i+1];
    // 	  if( HODObarid[i]-HODObarid[i+1]==-1 ) Rle_tot_B2 = HODORtot[i+1];

    // 	  if( HODObarid[i]-HODObarid[i+2]==-2 ) Rle_trYtime_B3 = HODORle[i+1];
    // 	  if( HODObarid[i]-HODObarid[i+2]==-2 ) Rle_tot_B3 = HODORtot[i+1];

    // 	  //cout << "B1 index: " << HODObarid[i] << " value:" << Rle_trYtime << endl;
    // 	  //cout << "B2 index: " << HODObarid[i+1] << " value:" << Rle_trYtime_B2 << endl;
    // 	}   
    // 	// cout << "cid index: " << 0 << ", value: " << HODOcid[0] << endl;
    // 	// cout << "cid index: " << 1 << ", value: " << HODOcid[1] << endl;
    // 	// cout << "bar index: " << i << ", value: " << HODObarid[i] << endl;

    // 	// if( HODOndata_cid>1&&HODObarid[i]==HODOcid[1] ){ //if the bar id matches the cluster second element id
    // 	//   Rle_trYtime_B2 = HODORle[i];
    // 	//   Rle_tot_B2 = HODORtot[i];
    // 	//   //cout << "B2: " << Rle_trYtime_B2 << endl;
    // 	// }
    //   }
    // }
    // htot->Fill(Rle_tot);
    // htot_B2->Fill(Rle_tot_B2);
    // htot_B3->Fill(Rle_tot_B3);


    // // for( Int_t i=0; i<HODOndata_cid; i++ ){
    // //   cout << "cid index: " << i << ", value: " << HODOcid[i] << endl;
    // // }
    // // for( Int_t i=0; i<HODOndata_Rle; i++ ){
    // //   cout << "bar Rtdc index: " << i << ", value: " << HODORle[i] << endl;
    // // }
    // // for( Int_t i=0; i<HODOndata_Rle; i++ ){
    // //   cout << "bar Rtdc index: " << i << ", value: " << HODObarid[i] << endl;
    // // }
    // // for( Int_t i=0; i<HODOndata_Rle; i++ ){
    // //   cout << "bar Rtot index: " << i << ", value: " << HODORtot[i] << endl;
    // // }

    // //Check index matching on barid and cluster id. Suspect mismatch.
    // // if ( HODOndata_Rle!=HODOndata_Rtot ) cout << "NDATA MISMATCH 1" << endl;
    // // if ( HODOndata_Rle!=HODOndata_id ) cout << "NDATA MISMATCH 2" << endl;
    // // bool nomatch = true;
    // // Int_t cid = -1;
    // // Int_t barid = -1;
    // // for( Int_t i=0; i<HODOndata_cid; i++ ){
    // //   for( Int_t j=0; j<HODOndata_Rle; j++ ){
    // // 	if( HODOcid[i]==HODObarid[j] ){
    // // 	  nomatch = false;
    // // 	  cid = HODOcid[i];
    // // 	  barid = HODObarid[j];
    // // 	}
    // //   }
    // // }

    // // if( nomatch==true ) cout << "NO MATCH on cluster/bar indexes. cid:" << cid << " barid:" << barid << endl;

    // //Fill histos with bar variable for first element with non-dispersive cut and with added time over thresh cut
    // if( Rle_trYtime!=-1000. ){
    //   hHODORle_trYslice->Fill(Rle_trYtime);
    //   hHODORle_trYslice_RFcorr->Fill(Rle_trYtime-RFmod);
    //   if( abs(Rle_tot-11.)<1. ) hHODORle_trYtotslice->Fill(Rle_trYtime);
    // }
    
    // //cout << "B1B2 diff: " << Rle_trYtime-Rle_trYtime_B2 << " = " << "Rle_trYtime:" << Rle_trYtime << " - " << "Rle_trYtime_B2:" << Rle_trYtime_B2 << endl;

    // //Fill similar histo with the first - second cluster element difference from bar variables
    // if( Rle_trYtime!=-1000.&&Rle_trYtime_B2!=-1000. ){
    //   hB1B2diff_trYslice->Fill(Rle_trYtime-Rle_trYtime_B2);
    //   hB1B2diff_trYslice_RFcorr->Fill(Rle_trYtime-Rle_trYtime_B2-RFmod);
    //   //cout << Rle_tot << " " << Rle_tot_B2 << endl;
    //   if( abs(Rle_tot-25.)<3.&&abs(Rle_tot_B2-22.)<3. ){
    // 	hB1B2diff_trYtotslice->Fill(Rle_trYtime-Rle_trYtime_B2);
    // 	hB1B2diff_trYtotslice_RFcorr->Fill(Rle_trYtime-Rle_trYtime_B2-RFmod);
    // 	//cout << "GOT THERE" << endl;
    //   }
    // }

    // //Fill similar histo with the first - third cluster element difference from bar variables
    // //if( Rle_trYtime!=-1000.&&Rle_trYtime_B2!=-1000.&&Rle_trYtime_B3!=-1000. ){
    // if( Rle_trYtime!=-1000.&&Rle_trYtime_B3!=-1000. ){
    //   hB1B3diff_trYslice->Fill(Rle_trYtime-Rle_trYtime_B3);
    //   hB1B3diff_trYslice_RFcorr->Fill(Rle_trYtime-Rle_trYtime_B3-RFmod);
    //   //cout << Rle_tot << " " << Rle_tot_B2 << endl;
    //   //if( abs(Rle_tot-25.)<3.&&abs(Rle_tot_B2-22.)<3.&&abs(Rle_tot_B3-21.)<3. ){
    //   if( abs(Rle_tot-25.)<3.&&abs(Rle_tot_B3-21.)<3. ){
    // 	hB1B3diff_trYtotslice->Fill(Rle_trYtime-Rle_trYtime_B3);
    // 	hB1B3diff_trYtotslice_RFcorr->Fill(Rle_trYtime-Rle_trYtime_B3-RFmod);
    // 	//cout << "GOT THERE" << endl;
    //   }
    // }

    // hHodoVsTrackP->Fill( BBtr_p[0], HODOtmean[0] );
    // hB1B2diff_hcal->Fill(HCALblktdc[0]-HCALblktdc[1]);    

    // //Just for the lololololols
    // //if( RFmod==0 ) continue;

    // hB1B2diff_hcal->Fill(HCALblktdc[0]-HCALblktdc[1]-RFmod);
    

  }
  
  cout << endl << endl;

  cout << "Total elastics on sample: " << elasYield << "/" << Nevents << endl << endl;
  cout << "Total hodo cluster/bar ID mismatches after cuts: " << clusBarMiss << "/" << elasYield << endl;

  fout->Write();
}
