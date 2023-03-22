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
  Int_t HODOndata_nblk;
  Int_t HODOndata_tmean;
  Int_t HODOndata_ymean;
  Int_t HODOndata_cbid;
  Int_t HODOndata_cbtmean;
  Int_t HODOndata_cbtot;
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
  Double_t HCALblktdc[kNcell], HCALblka[kNcell], HCALtdc[kNcell], HCALamp_p[kNcell], HCALtdcid[kNcell];
  Double_t HCALatime[kNcell], HCALa_p[kNcell];
  Double_t crow, ccol, nblk;
  Double_t cblkid[kNcell], cblke[kNcell];
  Double_t ce[kNcell], catime[kNcell], cid[kNcell], cnblk[kNcell], ctdc[kNcell];
  Double_t nclus;

  // Hodo params
  Double_t HODOtmean[kNcell];
  Double_t HODOnclus;
  Double_t HODOnblk[kNcell];
  Double_t HODOtdc[kNcell];
  Double_t HODOtdcid[kNcell];
  Double_t HODOtdc_tot[kNcell];
  Double_t HODOtdc_mult[kNcell];
  Double_t HODOymean[kNcell];
  Double_t HODOcid[kNcell];
  Double_t HODObarid[kNcell]; // bar id
  Double_t HODORle[kNcell];
  Double_t HODORtot[kNcell];
  Double_t HODOcbid[kNcell]; //cluster bar id
  Double_t HODOcbtmean[kNcell]; //cluster bar tdc mean time
  Double_t HODOcbtot[kNcell]; //cluster bar tdc mean time
  
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
  C->SetBranchStatus( "sbs.hcal.tdc", 1 );  
  C->SetBranchStatus( "sbs.hcal.a_amp_p", 1 );  
  C->SetBranchStatus( "sbs.hcal.tdcelemID", 1 );  
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
  C->SetBranchStatus( "bb.hodotdc.clus.size", 1 );
  C->SetBranchStatus( "bb.hodotdc.nclus", 1 );
  C->SetBranchStatus( "bb.hodotdc.tdc", 1 );
  C->SetBranchStatus( "bb.hodotdc.tdc_tot", 1 );
  C->SetBranchStatus( "bb.hodotdc.tdcelemID", 1 );
  C->SetBranchStatus( "bb.hodotdc.clus.ymean", 1 );
  C->SetBranchStatus( "bb.hodotdc.clus.bar.tdc.id", 1 );
  C->SetBranchStatus( "bb.hodotdc.clus.bar.tdc.meantime", 1 );
  C->SetBranchStatus( "bb.hodotdc.clus.bar.tdc.meantot", 1 );
  C->SetBranchStatus( "bb.hodotdc.bar.tdc.R.leW", 1 );
  C->SetBranchStatus( "bb.hodotdc.bar.tdc.R.tot", 1 );
  C->SetBranchStatus( "bb.hodotdc.bar.tdc.id", 1 );

  C->SetBranchStatus( "Ndata.bb.hodotdc.clus.id", 1 );
  C->SetBranchStatus( "Ndata.bb.hodotdc.clus.bar.tdc.id", 1 );
  C->SetBranchStatus( "Ndata.bb.hodotdc.clus.bar.tdc.meantime", 1 );
  C->SetBranchStatus( "Ndata.bb.hodotdc.clus.bar.tdc.meantot", 1 );
  C->SetBranchStatus( "Ndata.bb.hodotdc.clus.ymean", 1 );
  C->SetBranchStatus( "Ndata.bb.hodotdc.clus.size", 1 );
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
  C->SetBranchAddress( "sbs.hcal.tdc", HCALtdc );
  C->SetBranchAddress( "sbs.hcal.a_amp_p", HCALamp_p );
  C->SetBranchAddress( "sbs.hcal.tdcelemID", HCALtdcid );
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
  C->SetBranchAddress( "bb.hodotdc.clus.bar.tdc.id", HODOcbid );
  C->SetBranchAddress( "Ndata.bb.hodotdc.clus.bar.tdc.id", &HODOndata_cbid );
  C->SetBranchAddress( "bb.hodotdc.clus.bar.tdc.meantime", HODOcbtmean );
  C->SetBranchAddress( "Ndata.bb.hodotdc.clus.bar.tdc.meantime", &HODOndata_cbtmean );
  C->SetBranchAddress( "bb.hodotdc.clus.bar.tdc.meantot", HODOcbtot );
  C->SetBranchAddress( "Ndata.bb.hodotdc.clus.bar.tdc.meantot", &HODOndata_cbtot );
  C->SetBranchAddress( "bb.hodotdc.clus.id", HODOcid );
  C->SetBranchAddress( "Ndata.bb.hodotdc.clus.id", &HODOndata_cid );
  C->SetBranchAddress( "bb.hodotdc.clus.size", HODOnblk );
  C->SetBranchAddress( "Ndata.bb.hodotdc.clus.size", &HODOndata_nblk );
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

  //HODO histos
  TH1D *hrawmaxtot = new TH1D( "hrawmaxtot", "Hodoscope raw signal max time over threshold", 400, 0, 100 );
  TH1D *hrawB1B2 = new TH1D( "hrawB1B2", "Hodoscope raw signal high tot block - high tot block - 1", 400, -20, 20 ); 
  TH1D *hrawB1B3 = new TH1D( "hrawB1B3", "Hodoscope raw signal high tot block - high tot block - 2", 400, -20, 20 ); 
  TH1D *hrawB1B4 = new TH1D( "hrawB1B4", "Hodoscope raw signal high tot block - high tot block - 3", 400, -20, 20 ); 
  TH1D *hrawB1B5 = new TH1D( "hrawB1B5", "Hodoscope raw signal high tot block - high tot block - 4", 400, -20, 20 );  
  TH2D *hclusdiffID = new TH2D( "hclusdiffID", "Hodo Cluster Seed - Block Diff, (Primary Block Reftime, nblk>1), vs ID", 100, -10, 90, 400, -10, 10 );
  TH2D *hB1B2diffID = new TH2D( "hB1B2diffID", "Hodo Cluster Bar Seed - Next Highest ToT Bar Diff vs ID (uncorrected)", 100, -10, 90, 400, -10, 10 );
  TH1D *hpvabardiff = new TH1D( "hpvabardiff", "Hodo Cluster Bar Seed - Adjacent Bar Diff (uncorrected)", 400, -10, 10 );
  TH2D *hpvabardiffID = new TH2D( "hpvabardiffID", "Hodo Cluster Bar Seed - Adjacent Bar Diff vs ID (uncorrected)", 100, -10, 90, 400, -10, 10 );
  TH1D *hclusdiff = new TH1D( "hclusdiff", "Hodo Cluster Seed - Block Diff, (Primary Block Reftime, nblk>1)", 400, -10, 10 );
  TH2D *hclusmeanID = new TH2D( "hclusmeanID", "Hodo Cluster Mean Time, (Primary Block Reftime/ID, nblk>1), vs ID", 100, -10, 90, 400, -10, 10 );
  TH1D *hclusmean = new TH1D( "hclusmean", "Hodo Cluster Mean Time, (Primary Block Reftime/ID, nblk>1)", 400, -2, 2 );
  TH2D *hclusbardiffID = new TH2D( "hclusbardiffID", "Hodo Cluster Bar Seed - Block Diff, (Primary Block Reftime, nblk>1, no corrections), vs ID", 100, -10, 90, 400, -10, 10 );
  TH2D *hclusbarmeanID = new TH2D( "hclusbarmeanID", "Hodo Cluster Bar Mean Time, (Primary Block Reftime/ID, nblk>1, no corrections), vs ID", 100, -10, 90, 400, -10, 10 );
  TH1D *hclusbarmean = new TH1D( "hclusbarmean", "Hodo Cluster Bar Mean Time, (Primary Block Reftime/ID, nblk>1, no corrections)", 400, -10, 10 );
  TH2D *hclusbarB1vB2 = new TH2D( "hclusbarB1vB2", "Hodo Cluster Bar 1 vs Bar 2, (no corrections)", 400, -20, 20, 400, -20, 20 );
  TH2D *hbarB40vB41 = new TH2D( "hbarB40vB41", "Hodo Bar 40 vs Bar 41, (no corrections); bar40 (ns); bar41 (ns)", 400, -20, 20, 400, -20, 20 );



  //TDCTRIG histos
  TH1D *hrf = new TH1D( "hrf", "RF signal", 800, -100, 300 );

  //HCAL histos
  //TH1D *hhclusdiff = new TH1D( "hhclusdiff", "HCal Cluster Seed - Block Diff (Primary Block Reftime, nblk>1)", 400, -100, 100 );
  TH1D *hhclusmean = new TH1D( "hhclusmean", "HCal Cluster Mean Time (Primary Block Reftime, nblk>1)", 400, -20, 20 );
  TH2D *hhclusmeanID = new TH2D( "hhclusmeanID", "HCal Cluster Mean Time, (Primary Block Reftime/ID, nblk>1), vs ID", 288, 0, 288, 400, -20, 20 );
  TH1D *hhclusdiff = new TH1D( "hhclusdiff", "HCal Cluster Seed - Block Diff, (Primary Block Reftime, nblk>1)", 400, -20, 20 );  
  TH2D *hhclusdiffID = new TH2D( "hhclusdiffID", "HCal Cluster Seed - Block Diff, (Primary Block Reftime, nblk>1), vs ID", 288, 0, 288, 400, -20, 20 );
  TH2D *hhB1B2diffID = new TH2D( "hhB1B2diffID", "HCal Cluster Seed - Next Highest E Block Diff vs ID", 288, 0, 288, 400, -20, 20 );
  TH2D *hhpvablkdiffID = new TH2D( "hhpvablkdiffID", "HCal Cluster Seed - Adjacent Block Diff vs ID", 288, 0, 288, 400, -20, 20 );
  TH1D *hhpvadiff = new TH1D( "hhpvadiff", "HCal Block - Adjacent Block Diff", 400, -40, 40 );
  TH2D *hhpvadiffID = new TH2D( "hhpvadiffID", "HCal Block - Adjacent Block Diff vs ID", 288, 0, 288, 400, -40, 40 );
  TH2D *hhpvadiffcol = new TH2D( "hhpvadiffcol", "HCal Block - Adjacent Block Diff vs Column", kNcols, 0, kNcols, 400, -40, 40 );
  TH2D *hhpvadiffrow = new TH2D( "hhpvadiffrow", "HCal Block - Adjacent Block Diff vs Row", kNrows, 0, kNrows, 400, -40, 40 );
  TH2D *hhclusB1vB2 = new TH2D( "hhclusB1vB2", "HCal Cluster Block1 vs Block2", 1000, -300, 50, 1000, -300, 50 );
  TH2D *hhB127vB128 = new TH2D( "hhB127vB128", "HCal block 127 vs block 128; block 127 (ns); block 128 (ns)", 1000, -300, 50, 1000, -300, 50 );
  TH2D *hhampB127vB128 = new TH2D( "hhampB127vB128", "HCal block 127 vs block 128; block 127 (mV); block 128 (ns)", 500, 0, 100, 500, 0, 100 );

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
    
    //Cut on single bigbite tracks
    bool trn1 = BBtr_n==1;
    //if( !trn1 ) continue;

    //Cut on elastic peak in W2
    bool elastic = fabs( kineW2-W2_mean )<W2_sig;
    //if( !elastic ) continue;

    //Cut on hcal energy
    bool hcale = HCALe>0.025;
    if( !hcale ) continue;

    //Add slice on bb.tr.y
    bool inY = fabs(BBtr_y[0])<0.02;
    //if( !inY ) continue;

    Double_t RF_time = -1000;
    bool rfmult = false;
    Int_t rfnhits = 0;
    for(Int_t ihit=0; ihit<TDCTndata; ihit++){
      if(TDCT_id[ihit]==4){
	if(rfmult) cout << "Multiple hits in RF time channel, up to: " << rfnhits << endl;
	RF_time=TDCT_tdc[ihit];
	rfmult=true;
	rfnhits++;
      }
    }

    //HODO

    //Using existing cluster information
    Int_t hpblkid = HODOcbid[0];
    Double_t hblkseed = 0.;
    Double_t htavg = 0;
    Int_t hpnblk = HODOndata_cbid;
    //if( hpnblk>0 ) cout << "cluster block id, tmean, mtot - event " << nevent << " - cluster bars " << hpnblk << endl;
    for( Int_t b=0; b<hpnblk; b++ ){
      //if( hpnblk<2 ) continue; //for diagnostics
      Double_t hcblktime = HODOcbtmean[b];
      Double_t hcblkid = HODOcbid[b];
      Double_t hcblktot = HODOcbtot[b];
      Double_t htcdiff = hcblktime-hblkseed;

      //cout << hcblkid << " " << hcblktime << " " << hcblktot << endl;

      if( b==0 ){
	hblkseed = hcblktime;
      }else{
	htavg += htcdiff;
	//if( b==1 ) hB1B2diffID->Fill( hpblkid, htcdiff );
	//if( hcblkid == hpblkid+1 ) hpvabardiffID->Fill( hpblkid, htcdiff );
	hclusdiffID->Fill( hpblkid, htcdiff );
	hclusdiff->Fill( htcdiff );
      }
    }

    htavg /= (hpnblk-1);
    if( htavg!=0 ){
      hclusmeanID->Fill( hpblkid, htavg );
      hclusmean->Fill(htavg);
    }

    //cout << endl;

    //cout << "bar tdc id, time, tot" << endl;
    Double_t hbarseed = 0.;
    Double_t htbavg = 0;
    Double_t hbar40 = -100.;
    Double_t hbar41 = -100.;
    for( Int_t b=0; b<HODOndata_id; b++ ){ //loop over all bar.tdc data (larger set)
      Double_t hbartime = HODORle[b];
      Double_t hbarid = HODObarid[b];
      Double_t hbartot = HODORtot[b];

      //cout << "hbarid " << hbarid << endl;

      if(hbarid==40) hbar40 = hbartime;
      if(hbarid==41) hbar41 = hbartime;

      //cout << hbarid << " " << hbartime << " " << hbartot << endl;
      for( Int_t i=0; i<hpnblk; i++ ){ //loop over clus.bar.tdc data (smaller set)
	Double_t hcblkid = HODOcbid[i];
	Double_t hpblkid = HODOcbid[0];
	Double_t htcbdiff = hbartime-hbarseed;
	
	//cout << "!!hcblkid " << hcblkid << endl;

	if( hpnblk<2 ) continue;

	if(hbarid == hcblkid){

	  //cout << "got here " << hbarid << " " << hcblkid << endl;

	  if( i==0 ){ //the hodo cluster elements are ranked by tot, the bar elements are not
	    hbarseed = hbartime;
	    //cout << hbarseed << endl;
	  }else{
	    if( i==1 ){
	      hclusbarB1vB2->Fill(hbarseed,hbartime);
	      hB1B2diffID->Fill(hpblkid,htcbdiff);
	    }
	    if( hcblkid == hpblkid+1 ){
	      hpvabardiff->Fill( htcbdiff );
	      hpvabardiffID->Fill( hpblkid, htcbdiff );
	    }
	    htbavg += htcbdiff;
	    hclusbardiffID->Fill( hcblkid, htcbdiff );
	    //cout << htcbdiff << endl;
	  }	

	  //cout << "hbarseed " << hbarseed << " hbartime " << hbartime << endl;
  
	}
      }
    }

    if(hbar40!=-100.&&hbar41!=-100.) hbarB40vB41->Fill(hbar40,hbar41);

    htbavg /= (hpnblk-1);
    if( htbavg!=0 ){
      hclusbarmeanID->Fill( hpblkid, htbavg );
      hclusbarmean->Fill(htbavg);
    }


    //cout << HCALtdc[128] << endl;


    //cout << endl;

    //cout << endl;

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

    //HCAL
    Int_t pblkid = cblkid[0];
    Double_t blkseed = 0.;
    Double_t tavg = 0;

    //get block/adjacent block time without any clustering. sbs.hcal.tdc is a static array with 288 elements
    for( Int_t b=0; b<kNcell; b++ ){
      Int_t row = b/kNcols;
      Int_t col = b%kNcols;
      if( HCALtdc[b]< -900 || HCALblktdc[b] > 1000 ) continue; //check if data exists on block
      if( HCALtdc[b+1]< -900 || HCALblktdc[b+1] > 1000 ) continue; //check if data exists on adj. block
      Double_t tdiff = HCALtdc[b]-HCALtdc[b+1];
      if( tdiff==0. ) continue;
      hhpvadiff->Fill( tdiff );
      hhpvadiffID->Fill( b, tdiff );
      hhpvadiffcol->Fill( col, tdiff );
      hhpvadiffrow->Fill( row, tdiff );
      
    }


    if( nblk<2 || HCALblktdc[0]< -900 || HCALblktdc[0] > 1000 ) continue; //for diagnostics at first

    if( HCALtdc[128]<500&&HCALtdc[128]>-500&&HCALtdc[129]<500&&HCALtdc[129]>-500 ){
      hhB127vB128->Fill(HCALtdc[128],HCALtdc[129]);
      hhampB127vB128->Fill(HCALamp_p[128],HCALamp_p[129]);
    }

    for( Int_t b=0; b<nblk; b++ ){
      Double_t blktime = HCALblktdc[b];
      Double_t blkid = cblkid[b];

      if( b==0 ){
	blkseed = blktime;
      }else{
	Double_t tcdiff = blktime-blkseed;

	if( b==1 ){
	  hhclusB1vB2->Fill(blkseed,blktime);
	  hhB1B2diffID->Fill(pblkid,tcdiff);
	}
	if( blkid == pblkid+1 ) hhpvablkdiffID->Fill( pblkid, tcdiff );
	tavg += tcdiff;
	hhclusdiffID->Fill( blkid, tcdiff );
	hhclusdiff->Fill( tcdiff );
      }
    }

    tavg /= (nblk-1);
    if( tavg!=0 ){
      hhclusmeanID->Fill( pblkid, tavg );
      hhclusmean->Fill( tavg );
    }

  }
  
  cout << endl << endl;

  cout << "Total elastics on sample: " << elasYield << "/" << Nevents << endl << endl;
  cout << "Total hodo cluster/bar ID mismatches after cuts: " << clusBarMiss << "/" << elasYield << endl;

  fout->Write();
}
