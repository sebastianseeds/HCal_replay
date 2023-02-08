//SSeeds 2.6.23 - Simple script to obtain tdc signal efficiencies over all kinematics

#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TString.h"
#include <iostream>
#include <fstream>
#include "TGraph.h"
#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <TChain.h>
#include <TSystem.h>
#include <TStopwatch.h>
#include "TMath.h"
#include "TTreeFormula.h"
//BB
const Int_t maxTracks = 1000; // Reasonable limit on tracks to be stored per event
//HCal - Note that actual measurement of vertical length is 381.6cm, indicating that the MC figures are correct
const Int_t maxHCalChan = 288; // Total HCal channels
const Int_t maxHCalRows = 24; // Total HCal rows
const Int_t maxHCalCols = 12; // Total HCal cols
const Double_t HCALHeight = 0.365; // Height of HCal above beamline in m
const Double_t HCalblk_l_h_MC = 0.15494; // Horizontal length of all HCAL blocks in m from MC database
const Double_t HCalblk_l_v_MC = 0.15875; // Vertical length of all HCAL blocks in m from MC database
const Double_t posHCalXi_MC = -2.355005; // Distance from beam center to top of HCal in m from MC database
const Double_t posHCalXf_MC = 1.454995; // Distance from beam center to bottom of HCal in m from MC database
const Double_t posHCalYi_MC = -0.92964; // Distance from beam center to opposite-beam side of HCal in m from MC database
const Double_t posHCalYf_MC = 0.92964; // Distance from beam center to beam side of HCal in m from MC database
const Int_t nkine = 6; // Total number of kinematic points
const Int_t maxruns = 7; // Maximum number of available runs zero field LH2 per kinematic
const Int_t maxblock = 12; // Maximum number of blocks observed in clusters

//Physics/Math
const Double_t PI = TMath::Pi();
const Double_t M_e = 0.00051;
const Double_t M_p = 0.938272;
const Double_t M_n = 0.939565;
//Static Target/Scattering Chamber Parameters
const Double_t l_tgt = 0.15; // Length of the target (m)
const Double_t rho_tgt = 0.0723; // Density of target (g/cc)
const Double_t rho_Al = 2.7; // Density of aluminum windows (g/cc)
const Double_t celldiameter = 1.6*2.54; //cm, right now this is a guess
const Double_t dEdx_tgt=0.00574; //According to NIST ESTAR, the collisional stopping power of hydrogen is about 5.74 MeV*cm2/g at 2 GeV energy
const Double_t dEdx_Al = 0.0021; //According to NIST ESTAR, the collisional stopping power of Aluminum is about 2.1 MeV*cm2/g between 1-4 GeV
const Double_t uwallthick_LH2 = 0.0145; //cm

//Use opt=0 for all
void hcaltdc_eff( Int_t opt=0 ){

  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  TChain *C = new TChain("T");

  // Declare simple variables and strings to keep track of file extensions
  TString outputfilename = "hcaltdc_eff_allkin.root";

  // Declare general physics parameters over all kinematics
  Int_t kine[nkine] = {4,7,11,14,8,9}; // GMn nominal kinematic
  Double_t E_e[nkine] = {3.7278,7.906,9.91,5.965,4.013,6.373}; // Energy of beam (incoming electrons from accelerator)
  Double_t BB_d[nkine] = {1.7988,1.84896,1.55146,1.84787,1.97473,1.550}; // Distance to bigbite spectrometer from target chamber (m)
  Double_t BB_th[nkine] = {36.,40.,42.,46.5,26.5,49.}; // Angle BB spectrometer makes with exit beamline
  Double_t HCal_d[nkine] = {11.,14.,14.5,14.,11.,11.}; // Distance to HCal from scattering chamber for comm1
  Double_t HCal_th[nkine] = {31.9,16.1,13.3,17.3,29.4,22.}; // Angle HCal center makes with exit beamline  
  Double_t W2_mean[nkine] = {0.95,0.88,0.925,0.94,0.91,0.91}; // Mean of W at current kinematic
  Double_t W2_sig[nkine] = {0.25,0.5,0.325,0.2,0.21,0.21}; // Width of W at current kinematic
  Int_t zfruns[nkine] = {3,0,2,4,7,0}; // Number of zero field runs available for analysis

  // Declare global cuts for all kinematics
  TCut globalcut[nkine] = {
    "bb.tr.n==1&&bb.ps.e>0.2&&abs(bb.tr.vz[0])<0.08&&bb.gem.track.nhits>3&&abs(bb.etot_over_p-0.92)<0.2&&sbs.hcal.e>0.10&&bb.ps.e+bb.sh.e>1.7",
    "bb.tr.n>0&&abs(bb.tr.vz[0])<0.08&&bb.gem.track.nhits>2&&bb.tr.p[0]>2.0&&sbs.hcal.e>0.10&&bb.ps.e>0.2",
    "bb.tr.n>0&&abs(bb.tr.vz[0])<0.08&&bb.gem.track.nhits>2&&bb.tr.p[0]>2.0&&sbs.hcal.e>0.10&&bb.ps.e>0.2",
    "bb.tr.n>0&&abs(bb.tr.vz[0])<0.08&&bb.gem.track.nhits>2&&bb.tr.p[0]>1.6&&sbs.hcal.e>0.10&&bb.ps.e>0.2",
    "bb.tr.n==1&&bb.ps.e>0.2&&abs(bb.tr.vz[0])<0.08&&bb.gem.track.nhits>3&&sbs.hcal.e>0.10&&abs(bb.tr.tg_th[0])<0.15&&abs(bb.tr.tg_ph[0])<0.3",
    "bb.tr.n==1&&bb.ps.e>0.2&&abs(bb.tr.vz[0])<0.08&&bb.gem.track.nhits>3&&sbs.hcal.e>0.10&&abs(bb.tr.tg_th[0])<0.15&&abs(bb.tr.tg_ph[0])<0.3"
  };

  //Read in data from zero field runs only to ensure expected detections from straight line projections to hcal
  Int_t ninfile[nkine][maxruns] = {0};
  //SBS4
  ninfile[0][0] = 11573;
  ninfile[0][1] = 11587;
  ninfile[0][2] = 11588;
  //SBS11
  ninfile[2][0] = 12730;
  ninfile[2][1] = 12731;
  //SBS14
  ninfile[3][0] = 13375;
  ninfile[3][1] = 13376;
  ninfile[3][2] = 13377;
  ninfile[3][3] = 13378;
  //SBS8
  ninfile[4][0] = 13459;
  ninfile[4][1] = 13460;
  ninfile[4][2] = 13461;
  ninfile[4][3] = 13463;
  ninfile[4][4] = 13464;
  ninfile[4][5] = 13465;
  ninfile[4][6] = 13466;

  TTreeFormula *GlobalCut[nkine];
  TEventList *elist[nkine];
  //Long64_t Nevents;

  // Base histograms
  TH1D *hW2_nocut[nkine];
  TH1D *hdy_BBcut[nkine];
  TH1D *hdx_BBcut[nkine];
  TH1D *httime_BBcut[nkine];
  TH1D *hblkID[nkine];
  
  TH1D *httime_blkN_BBcut[nkine][maxblock];

  // Run params
  UInt_t runI; 

  // BB params
  Double_t BBtr_p[maxTracks], BBtr_px[maxTracks], BBtr_py[maxTracks], BBtr_pz[maxTracks];
  Double_t BBtr_vz[maxTracks];
  Double_t BB_eOp, BBgem_nplanes, BBps_e, BBsh_e; 
  Double_t BBtgt_th[maxTracks], BBtgt_ph[maxTracks];
  
  // HCal params 
  Double_t cid[maxHCalChan], crow[maxHCalRows], ccol[maxHCalCols];
  Double_t ce[maxHCalChan], catime[maxHCalChan], ctdc[maxHCalChan];
  Double_t nblk;
  Double_t hcalx, hcaly, hcale, kineW2;

  for( Int_t k=0; k<nkine; k++ ){

    cout << "Adding events from kinematic " << k << " on sbs " << kine[k] << endl;

    for( Int_t r=0; r<maxruns; r++){
      Int_t runnum = ninfile[k][r];
      if( runnum==0 ) continue;
      //cout << runnum << endl;
      if( k<2 ) {
	C->Add(Form("/work/halla/sbs/sbs-gmn/pass0/SBS%d/LH2/rootfiles/e1209019_fullreplay_%d*",kine[k],runnum));
	cout << "Added run: " << Form("/work/halla/sbs/sbs-gmn/pass0/SBS%d/LH2/rootfiles/e1209019_fullreplay_%d*",kine[k],runnum) << endl;
      }else{
	C->Add(Form("/work/halla/sbs/sbs-gmn/pass1/SBS%d/LH2/rootfiles/e1209019_fullreplay_%d*",kine[k],runnum));
	cout << "Added run: " << Form("/work/halla/sbs/sbs-gmn/pass1/SBS%d/LH2/rootfiles/e1209019_fullreplay_%d*",kine[k],runnum);
      }
    }

  }

  //elist[k] = new TEventList(Form("elist_sbs%d",kine[k]),Form("Elastic Event List, SBS%d",kine[k]));

  //cout << "For kinematic " << kine[k] << ", making globalcut: " << globalcut[k] << endl;
  //C[k]->Draw(Form(">>elist_sbs%d",kine[k]),globalcut[k]);
  //C[k]->Draw(Form(">>elist[%d]",k),globalcut[k]);
  //Nevents[k] = elist[k]->GetN();
  //cout << "On kinematic " << kine[k] << " " << Nevents[k] << " passed global cut." << endl;
    
  Long64_t Nevents = C->GetEntries();
  cout << "Got " << Nevents << " events over all kinematics." << endl;

  //Switch on branches and link them
  C->SetBranchStatus( "*", 0 );
  C->SetBranchStatus( "sbs.hcal.x", 1 );
  C->SetBranchStatus( "sbs.hcal.y", 1 );
  C->SetBranchStatus( "sbs.hcal.e", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.row", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.col", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.tdctime", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.atime", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.id", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.e", 1 );
  C->SetBranchStatus( "sbs.hcal.nblk", 1 );
  C->SetBranchStatus( "bb.tr.p", 1 );
  C->SetBranchStatus( "bb.tr.px", 1 );
  C->SetBranchStatus( "bb.tr.py", 1 );
  C->SetBranchStatus( "bb.tr.pz", 1 );
  C->SetBranchStatus( "bb.tr.vz", 1 );
  C->SetBranchStatus( "bb.tr.n", 1 );
  C->SetBranchStatus( "bb.ps.e", 1 );
  C->SetBranchStatus( "bb.sh.e", 1 );
  C->SetBranchStatus( "e.kine.W2", 1 );
  C->SetBranchStatus( "bb.etot_over_p", 1 );
  C->SetBranchStatus( "bb.gem.track.nhits", 1 );
  C->SetBranchStatus( "bb.tr.tg_th", 1 );
  C->SetBranchStatus( "bb.tr.tg_ph", 1 );
  C->SetBranchStatus( "Event_Branch", 1 );
  C->SetBranchStatus( "fEvtHdr.fRun", 1 );

  //Link the branches to vars
  C->SetBranchAddress( "sbs.hcal.x", &hcalx );
  C->SetBranchAddress( "sbs.hcal.y", &hcaly );
  C->SetBranchAddress( "sbs.hcal.e", &hcale );
  C->SetBranchAddress( "sbs.hcal.clus_blk.row", crow );
  C->SetBranchAddress( "sbs.hcal.clus_blk.col", ccol );
  C->SetBranchAddress( "sbs.hcal.clus_blk.tdctime", ctdc );
  C->SetBranchAddress( "sbs.hcal.clus_blk.atime", catime );
  C->SetBranchAddress( "sbs.hcal.clus_blk.id", cid );
  C->SetBranchAddress( "sbs.hcal.clus_blk.e", ce );
  C->SetBranchAddress( "sbs.hcal.nblk", &nblk );
  C->SetBranchAddress( "bb.tr.p", BBtr_p );
  C->SetBranchAddress( "bb.tr.px", BBtr_px );
  C->SetBranchAddress( "bb.tr.py", BBtr_py );
  C->SetBranchAddress( "bb.tr.pz", BBtr_pz );
  C->SetBranchAddress( "bb.tr.vz", BBtr_vz );
  C->SetBranchAddress( "bb.ps.e", &BBps_e );
  C->SetBranchAddress( "bb.sh.e", &BBsh_e );
  C->SetBranchAddress( "e.kine.W2", &kineW2 );
  C->SetBranchAddress( "bb.etot_over_p", &BB_eOp );
  C->SetBranchAddress( "bb.gem.track.nhits", &BBgem_nplanes );
  C->SetBranchAddress( "bb.tr.tg_th", BBtgt_th );
  C->SetBranchAddress( "bb.tr.tg_ph", BBtgt_ph );
  C->SetMakeClass(1);
  C->SetBranchAddress( "fEvtHdr.fRun", &runI );

  // Set up global cut and histograms for each kinematic
  for( Int_t k=0; k<nkine; k++ ){

    //Set up the global cut to be evaluated on each loop
    GlobalCut[k]= new TTreeFormula( Form("GlobalCut_sbs%d",kine[k]), globalcut[k], C );

    //Set up histograms      
    hW2_nocut[k] = new TH1D( Form("hW2_nocut_sbs%d",kine[k]), Form("W^2 sbs%d; (GeV)",kine[k]), 400, 0.0, 2.0 );
    hdy_BBcut[k] = new TH1D( Form("hdy_BBcut_sbs%d",kine[k]), Form("dy expect cut sbs%d;y_{HCAL}-y_{expect} (m)",kine[k]), 500, posHCalYi_MC-2*HCalblk_l_h_MC, posHCalYf_MC+2*HCalblk_l_h_MC );
    hdx_BBcut[k] = new TH1D( Form("hdx_BBcut_sbs%d",kine[k]), Form("dx expect cut sbs%d;x_{HCAL}-x_{expect} (m)",kine[k]), 500, posHCalXi_MC-2*HCalblk_l_v_MC, posHCalXf_MC+2*HCalblk_l_v_MC );  
    httime_BBcut[k] = new TH1D( Form("httime_BBcut_sbs%d",kine[k]), Form("Cluster TDC sbs%d, primary blk, expect cut; ns",kine[k]), 1200, -2000, 10000 );
    hblkID[k] = new TH1D( Form("hblkID_sbs%d",kine[k]), Form("Primary Block ID, sbs%d",k), 300, -10, 290 );
    
    for( Int_t b=1; b<=maxblock; b++ ){
      httime_blkN_BBcut[k][b-1] = new TH1D( Form( "httime_sbs%d_blkN_%d",kine[k],b ), Form( "Cluster TDC sbs%d, blk %d, BB cut, Active Cut, Expect Cut; ns",kine[k],b ), 1200, -2000, 10000 );  
    
    }
  }
  
  //Outfiles
  TFile *fout = new TFile(outputfilename,"RECREATE");

  //Accounting and diagnostic variables
  Long64_t nexpect[nkine] = {0};
  Long64_t nphit[nkine] = {0};
  Long64_t nrhit[nkine] = {0};
  
  cout << endl << "Preamble complete. Proceeding to analysis." << endl;
    
  // Set variables to keep track of event sequence for ttreeformula
  Int_t treenum = 0, currenttreenum = 0;
  Int_t cs = 0; //Counter to keep track of current kinematic step

  cout << endl << "Processing events.." << endl << endl;
  
  for( Long64_t nevent = 1; nevent <Nevents; nevent++){

    // Should be able to do better with a map
    if( runI<11589 ) cs = 0; //sbs4
    if( runI>12729&&runI<12732 ) cs = 2; //sbs11
    if( runI>13374&&runI<13379 ) cs = 3; //sbs14
    if( runI>13458 ) cs = 4; //sbs8

    if ( nevent%1000==0 ) cout << "Kinematic: " << kine[cs] << ", entry: " << nevent << "/" << Nevents << "/r";
    cout.flush();

    //C->GetEntry( elist[cs]->GetEntry(nevent) );
    C->GetEntry( nevent );

    currenttreenum = C->GetTreeNumber();
    if( nevent == 1 || currenttreenum != treenum ){
      treenum = currenttreenum; 
      GlobalCut[cs]->UpdateFormulaLeaves();
    }

    bool failedglobal = GlobalCut[cs]->EvalInstance(0) == 0;
    if( failedglobal ) continue;


    Double_t pBeam = E_e[cs]/(1.+E_e[cs]/M_p*(1.-cos(BB_th[cs]))); // Momentum of beam calculated from theta
    Double_t Eloss_outgoing = celldiameter/2.0/sin(BB_th[cs]) * rho_tgt * dEdx_tgt; //Approximately 1 MeV, could correct further with raster position
    //Correct the beam energy with energy loss in target using vertex position
    Double_t Eloss = (BBtr_vz[0]+l_tgt/2.0) * rho_tgt * dEdx_tgt + uwallthick_LH2 * rho_Al * dEdx_Al; //approximately 3 MeV
    Double_t E_corr = E_e[cs] - Eloss;

    //Physics Calculations
    Double_t p_corr = BBtr_p[0] - Eloss_outgoing; //Neglecting the mass of e'
    Double_t etheta = acos( BBtr_pz[0]/BBtr_p[0] ); //Will use the uncorrected track momentum to reconstruct e' theta
    Double_t ephi = atan2( BBtr_py[0], BBtr_px[0] );   
    TVector3 vertex( 0, 0, BBtr_vz[0] ); // z location of vertex in hall coordinates
    TLorentzVector Pbeam( 0, 0, E_corr, E_corr ); //Mass of e negligable
    TLorentzVector kprime( BBtr_px[0], BBtr_py[0], BBtr_pz[0], BBtr_p[0] );
    TLorentzVector Ptarg( 0, 0, 0, M_p );   
    TLorentzVector q = Pbeam - kprime;
    TLorentzVector PgammaN = Ptarg + q; //(-px, -py, ebeam - pz, Mp + ebeam - p)     
    Double_t pel = E_corr/( 1.+E_corr/M_p*( 1.-cos(etheta) ) );
    Double_t nu = E_corr - BBtr_p[0];
    Double_t pp = sqrt( pow(nu,2)+2.*M_p*nu );
    Double_t phinucleon = ephi + PI; //assume coplanarity
    Double_t thetanucleon = acos( (E_corr - BBtr_pz[0])/pp ); //use elastic constraint on nucleon kinematics     
    TVector3 pNhat( sin(thetanucleon)*cos(phinucleon), sin(thetanucleon)*sin(phinucleon), cos(thetanucleon) );
      
    //Define HCal coordinate system
    TVector3 HCAL_zaxis( sin(-HCal_th[cs]), 0, cos(-HCal_th[cs]) );
    TVector3 HCAL_xaxis( 0, -1, 0 );
    TVector3 HCAL_yaxis = HCAL_zaxis.Cross(HCAL_xaxis).Unit();      
    TVector3 HCAL_origin = HCal_d[cs] * HCAL_zaxis + HCALHeight * HCAL_xaxis;     
    TVector3 HCALpos = HCAL_origin + hcalx * HCAL_xaxis + hcaly * HCAL_yaxis;
    Double_t sintersect = ( HCAL_origin - vertex ).Dot( HCAL_zaxis ) / (pNhat.Dot( HCAL_zaxis ) ); 
    TVector3 HCAL_intersect = vertex + sintersect * pNhat;   

    //Calculate quantities of interest
    Double_t yexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_yaxis );
    Double_t xexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_xaxis ); 
    Double_t dx = hcalx - xexpect_HCAL;
    Double_t dy = hcaly - yexpect_HCAL;
    Double_t E_ep = sqrt( pow(M_e,2) + pow(BBtr_p[0],2) ); // Obtain the scattered electron energy 
    Double_t p_ep = BBtr_p[0];     
    Double_t Q2 = 2*E_corr*E_ep*( 1-(BBtr_pz[0]/p_ep) ); // Obtain Q2 from beam energy, outgoing electron energy, and momenta  
	
    Double_t W2 = kineW2;

    //Record best tdc time per event
    Int_t pblkid = (int)(cid[0]-1); //Must correct for index switch on blocks

    //Check blk index
    hblkID[cs]->Fill(pblkid);
	
    //Fill histograms without any cuts
    hW2_nocut[cs]->Fill( W2 );
    //httime_nocut[cs]->Fill( ctdc[0] );

    //perform w2 cut to tighten elastic selection
    bool elastic = abs( W2 - W2_mean[cs] ) < W2_sig[cs];
	
    bool exphcalactive = false;
    if(  yexpect_HCAL > (posHCalYf_MC-HCalblk_l_h_MC) &&
	 yexpect_HCAL < (posHCalYi_MC+HCalblk_l_h_MC) &&
	 xexpect_HCAL > (posHCalXf_MC-HCalblk_l_v_MC) &&
	 xexpect_HCAL < (posHCalXi_MC+HCalblk_l_v_MC)){
	  
      exphcalactive=true;
    }
    bool hcalactive = false;
    if( crow[0] > 0 &&
	crow[0] < 23 &&
	ccol[0] > 0 &&
	ccol[0] < 11 ){
      hcalactive = true;
    }

    //Cut all events where no good tdc is expected in hcal and a cluster exists in HCal from fADC (BBcut)
    if( elastic == false || 
	exphcalactive == false || 
	hcalactive == false 
	) 
      continue;

    nexpect[cs]++;

    //Fill histograms with BigBite and HCal active area cuts
    hdx_BBcut[cs]->Fill( dx );
    hdy_BBcut[cs]->Fill( dy );

    // Record primary block tdc detections. This includes reconstructed hits for efficiency calc
    if( ctdc[0]>-500&&ctdc[0]<500 ){
      httime_blkN_BBcut[cs][0]->Fill( ctdc[0] );
      nphit[cs]++;
      nrhit[cs]++;
    }else{ // See if time can be reconstructed from additional blocks in cluster
      bool ghit = false;
      for( Int_t b=0; b<nblk; b++ ){ //Blocks are sorted in descending energy
	httime_blkN_BBcut[cs][b]->Fill( ctdc[b] );

	if( ctdc[b]>-500&&ctdc[b]<500&&ghit==false ){
	  nrhit[cs]++;
	  ghit=true;
	  continue;
	}
      }
    }
  }

    //}
  
  //Int_t n = nkine;
  Double_t TDC_peff[nkine];
  Double_t TDC_reff[nkine];

  for( Int_t k=0; k<nkine; k++ ){
    if( nexpect[k]==0 ) nexpect[k]=1; //Return zero if no events exist on the kinematic

    TDC_peff[k] = (Double_t)nphit[k] / (Double_t)nexpect[k];
    TDC_reff[k] = (Double_t)nrhit[k] / (Double_t)nexpect[k];
  }

  //TGraph *cpeff = new TGraph( nkine, kine, TDC_peff );
  //TGraph *creff = new TGraph( nkine, kine, TDC_reff );


  TCanvas *c1 = new TCanvas("c1","HCal TDC Detection Efficiency by Kinematic",1600,1000);
  
  c1->cd();
  
  //cpeff->SetMarkerStyle(20);
  //cpeff->Draw("AP");
  
  //creff->SetMarkerStyle(21);
  //creff->Draw();

  for( Int_t k=0; k<nkine; k++ ){
    
    cout << "TDC primary detection efficiency for kin " << kine[k] << ": " << TDC_peff[k] << endl;

  }

  cout << endl << endl;

  for( Int_t k=0; k<nkine; k++ ){
    
    cout << "TDC reconstructed detection efficiency for kin " << kine[k] << ": " << TDC_reff[k] << endl;

  }

  fout->Write();

  st->Stop();

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

}
