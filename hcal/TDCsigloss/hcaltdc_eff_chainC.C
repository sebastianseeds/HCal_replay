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
#include "TGraphErrors.h"
#include <ctime>
#include <fstream>
#include <string>
#include <TChain.h>
#include <TSystem.h>
#include <TStopwatch.h>
#include <stdlib.h>
#include <stdio.h>

#include "TMath.h"
#include "TTreeFormula.h"
//BB
const Int_t maxTracks = 1000; // Reasonable limit on tracks to be stored per event
//HCal - Note that actual measurement of vertical length is 381.6cm, indicating that the MC figures are correct
const Int_t maxHCalChan = 288; // Total HCal channels
const Int_t maxHCalRows = 24; // Total HCal rows
const Int_t maxHCalCols = 12; // Total HCal cols
//const Double_t HCALHeight = 0.365; // Height of HCal above beamline in m
const Double_t HCALHeight = -0.2897;
const Double_t HCalblk_l_h_MC = 0.15494; // Horizontal length of all HCAL blocks in m from MC database
const Double_t HCalblk_l_v_MC = 0.15875; // Vertical length of all HCAL blocks in m from MC database
const Double_t posHCalXi_MC = -2.355005; // Distance from beam center to top of HCal in m from MC database
const Double_t posHCalXf_MC = 1.454995; // Distance from beam center to bottom of HCal in m from MC database
const Double_t posHCalYi_MC = -0.92964; // Distance from beam center to opposite-beam side of HCal in m from MC database
const Double_t posHCalYf_MC = 0.92964; // Distance from beam center to beam side of HCal in m from MC database
const Int_t nkine = 6; // Total number of kinematic points
const Int_t maxruns = 20; // Maximum number of available runs zero field LH2 per kinematic
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
void hcaltdc_eff_chainC( Int_t opt=0 ){

  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  TChain *C[nkine];

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
  Double_t magoff[nkine] = {0.01621,-0.6924,-0.6307,-0.7245,0.1118,-0.8206}; // P peak fit offset to HCal X expected to account for SBS magnet
  Double_t cenKE_n[nkine] = {1.61,5.24,7.24,3.96,2.39,2.39}; // Central kinetic energy of the scattered nucleon
  Double_t atdiffc[nkine] = {120.,125.,110.,105.,110.,105.}; // ADC time - TDC central value
  Double_t atdiffs[nkine] = {15.,15.,30.,30.,35.,25.}; // ADC time - TDC max spread on each side

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
  //SBS4 - all zero field
  ninfile[0][0] = 11573;
  ninfile[0][1] = 11587;
  ninfile[0][2] = 11588;
  //SBS7 - all 50% field
  ninfile[1][0] = 11994;
  ninfile[1][1] = 12000;
  ninfile[1][2] = 12008;
  ninfile[1][3] = 12022;
  ninfile[1][4] = 12035;
  ninfile[1][5] = 12049;
  ninfile[1][6] = 12053;
  //SBS11 - all zero field
  //ninfile[2][0] = 12730;
  //ninfile[2][1] = 12731;
  //SBS11 - all 100% field
  ninfile[2][0] = 12313;
  ninfile[2][1] = 12320;
  ninfile[2][2] = 12345;
  ninfile[2][3] = 12355;
  ninfile[2][4] = 12358;
  ninfile[2][5] = 12363;
  ninfile[2][6] = 12471;
  ninfile[2][7] = 12959;
  ninfile[2][8] = 12525;
  ninfile[2][9] = 12931;
  ninfile[2][10] = 12910;
  ninfile[2][11] = 13027;
  ninfile[2][12] = 13028;
  ninfile[2][13] = 13042;
  ninfile[2][14] = 13053;
  //SBS14 - all zero field
  //ninfile[3][0] = 13375;
  //ninfile[3][1] = 13376;
  //ninfile[3][2] = 13377;
  //ninfile[3][3] = 13378;
  //SBS14 - all 70% field
  ninfile[3][0] = 13312;
  ninfile[3][1] = 13313;
  ninfile[3][2] = 13321;
  ninfile[3][3] = 13345;
  ninfile[3][4] = 13351;
  ninfile[3][5] = 13379;
  ninfile[3][6] = 13397;
  //SBS8 - all zero field
  ninfile[4][0] = 13459;
  ninfile[4][1] = 13460;
  ninfile[4][2] = 13461;
  //ninfile[4][3] = 13463;
  //ninfile[4][4] = 13464;
  //ninfile[4][5] = 13465;
  //ninfile[4][6] = 13466;
  //SBS9 - all 70% field
  ninfile[5][0] = 13676;
  ninfile[5][1] = 13683;
  ninfile[5][2] = 13696;

  TTreeFormula *GlobalCut[nkine];
  TEventList *elist[nkine];
  Long64_t Nevents[nkine];
  Long64_t NTevents[nkine];

  // Base histograms
  TH1D *hW2_nocut[nkine];
  TH1D *hdy_BBcut[nkine];
  TH1D *hdx_BBcut[nkine];
  TH1D *hblkID[nkine];
  TH1D *httime_blkN_BBcut[nkine][maxblock];
  TH1D *hatime_blkN_BBcut[nkine][maxblock];
  TH2D *httime_blkID[nkine];
  TH2D *hatdiff_blkID[nkine];

  // BB params
  Double_t BBtr_p[nkine][maxTracks], BBtr_px[nkine][maxTracks], BBtr_py[nkine][maxTracks], BBtr_pz[nkine][maxTracks];
  Double_t BBtr_vz[nkine][maxTracks];
  
  // HCal params 
  Double_t cid[nkine][maxHCalChan], crow[nkine][maxHCalRows], ccol[nkine][maxHCalCols];
  Double_t ce[nkine][maxHCalChan], catime[nkine][maxHCalChan], ctdc[nkine][maxHCalChan];
  Double_t nblk[nkine];
  Double_t hcalx[nkine], hcaly[nkine], hcale[nkine], kineW2[nkine];

  for( Int_t k=0; k<nkine; k++ ){

    cout << "Setting up chain " << k << " on sbs " << kine[k] << " with runs " << endl;

    C[k] = new TChain("T");
    for( Int_t r=0; r<maxruns; r++){
      Int_t runnum = ninfile[k][r];
      if( runnum==0 ) continue;
      //cout << runnum << endl;
      if( k<2 ){
	C[k]->Add(Form("/work/halla/sbs/sbs-gmn/pass0/SBS%d/LH2/rootfiles/e1209019_fullreplay_%d*",kine[k],runnum));
      }else if(k==2){
	C[k]->Add(Form("/volatile/halla/sbs/sbs-gmn/GMN_REPLAYS/pass1/SBS%d/LH2/rootfiles/e1209019_fullreplay_%d*",kine[k],runnum));
      }else{
	C[k]->Add(Form("/work/halla/sbs/sbs-gmn/pass1/SBS%d/LH2/rootfiles/e1209019_fullreplay_%d*",kine[k],runnum));
      }
      cout << runnum << " ";
    }
    cout << endl;

      //cout << "finished." << endl;

    // for( Int_t k=0; k<nkine; k++){
    //   NTevents[k] = C[k]->GetEntries();
    // }

    NTevents[k] = C[k]->GetEntries();

    elist[k] = new TEventList(Form("elist_sbs%d",kine[k]),Form("Elastic Event List, SBS%d",kine[k]));

    cout << "For kinematic " << kine[k] << ", making globalcut: " << endl << globalcut[k] << endl;
    C[k]->Draw(Form(">>elist_sbs%d",kine[k]),globalcut[k]);
    //C[k]->Draw(Form(">>elist[%d]",k),globalcut[k]);
    Nevents[k] = elist[k]->GetN();
    cout << "On kinematic " << kine[k] << ", " << Nevents[k] << " out of " << NTevents[k] << " passed global cut." << endl << endl;
    
    //Nevents[k] = C[k]->GetEntries();
    //cout << "On kinematic " << kine[k] << " got " << Nevents[k] << " events." << endl;

    //Switch on branches and link them
    C[k]->SetBranchStatus( "*", 0 );
    C[k]->SetBranchStatus( "sbs.hcal.x", 1 );
    C[k]->SetBranchStatus( "sbs.hcal.y", 1 );
    C[k]->SetBranchStatus( "sbs.hcal.e", 1 );
    C[k]->SetBranchStatus( "sbs.hcal.clus_blk.row", 1 );
    C[k]->SetBranchStatus( "sbs.hcal.clus_blk.col", 1 );
    C[k]->SetBranchStatus( "sbs.hcal.clus_blk.tdctime", 1 );
    C[k]->SetBranchStatus( "sbs.hcal.clus_blk.atime", 1 );
    C[k]->SetBranchStatus( "sbs.hcal.clus_blk.id", 1 );
    C[k]->SetBranchStatus( "sbs.hcal.clus_blk.e", 1 );
    C[k]->SetBranchStatus( "sbs.hcal.nblk", 1 );
    C[k]->SetBranchStatus( "bb.tr.p", 1 );
    C[k]->SetBranchStatus( "bb.tr.px", 1 );
    C[k]->SetBranchStatus( "bb.tr.py", 1 );
    C[k]->SetBranchStatus( "bb.tr.pz", 1 );
    C[k]->SetBranchStatus( "bb.tr.vz", 1 );
    C[k]->SetBranchStatus( "bb.tr.n", 1 );
    C[k]->SetBranchStatus( "e.kine.W2", 1 );

    //Link the branches to vars
    C[k]->SetBranchAddress( "sbs.hcal.x", &hcalx[k] );
    C[k]->SetBranchAddress( "sbs.hcal.y", &hcaly[k] );
    C[k]->SetBranchAddress( "sbs.hcal.e", &hcale[k] );
    C[k]->SetBranchAddress( "sbs.hcal.clus_blk.row", crow[k] );
    C[k]->SetBranchAddress( "sbs.hcal.clus_blk.col", ccol[k] );
    C[k]->SetBranchAddress( "sbs.hcal.clus_blk.tdctime", ctdc[k] );
    C[k]->SetBranchAddress( "sbs.hcal.clus_blk.atime", catime[k] );
    C[k]->SetBranchAddress( "sbs.hcal.clus_blk.id", cid[k] );
    C[k]->SetBranchAddress( "sbs.hcal.clus_blk.e", ce[k] );
    C[k]->SetBranchAddress( "sbs.hcal.nblk", &nblk[k] );
    C[k]->SetBranchAddress( "bb.tr.p", BBtr_p[k] );
    C[k]->SetBranchAddress( "bb.tr.px", BBtr_px[k] );
    C[k]->SetBranchAddress( "bb.tr.py", BBtr_py[k] );
    C[k]->SetBranchAddress( "bb.tr.pz", BBtr_pz[k] );
    C[k]->SetBranchAddress( "bb.tr.vz", BBtr_vz[k] );
    C[k]->SetBranchAddress( "e.kine.W2", &kineW2[k] );

    //Set up the global cut to be evaluated on each loop
    //GlobalCut[k]= new TTreeFormula( Form("GlobalCut_sbs%d",kine[k]), globalcut[k], C[k] );

    //Set up histograms      
    hW2_nocut[k] = new TH1D( Form("hW2_nocut_sbs%d",kine[k]), Form("W^2 sbs%d; (GeV)",kine[k]), 400, 0.0, 2.0 );
    hdy_BBcut[k] = new TH1D( Form("hdy_BBcut_sbs%d",kine[k]), Form("dy expect cut sbs%d;y_{HCAL}-y_{expect} (m)",kine[k]), 500, posHCalYi_MC-2*HCalblk_l_h_MC, posHCalYf_MC+2*HCalblk_l_h_MC );
    hdx_BBcut[k] = new TH1D( Form("hdx_BBcut_sbs%d",kine[k]), Form("dx expect cut sbs%d;x_{HCAL}-x_{expect} (m)",kine[k]), 500, posHCalXi_MC-2*HCalblk_l_v_MC, posHCalXf_MC+2*HCalblk_l_v_MC );  
    hblkID[k] = new TH1D( Form("hblkID_sbs%d",kine[k]), Form("Primary Block ID, sbs%d",k), 300, -10, 290 );
    httime_blkID[k] = new TH2D( Form("httime_blkID_sbs%d",kine[k]), Form("Primary blk tdc vs ID, expect cut, sbs%d; blk; ns",kine[k]), 300, -10, 280, 400, -300, 100);
    hatdiff_blkID[k] = new TH2D( Form("hatdiff_blkID_sbs%d",kine[k]), Form("(atime - tdc)_{pblk} vs ID, expect cut, sbs%d; blk; ns",kine[k]), 300, -10, 280, 600, -300, 300);

    for( Int_t b=1; b<=maxblock; b++ ){
      httime_blkN_BBcut[k][b-1] = new TH1D( Form( "httime_sbs%d_blkN_%d",kine[k],b ), Form( "Cluster TDC sbs%d, blk %d, BB cut, Expect Cut; ns",kine[k],b ), 400, -300, 100 );
      hatime_blkN_BBcut[k][b-1] = new TH1D( Form( "hatime_sbs%d_blkN_%d",kine[k],b ), Form( "Cluster ADC time sbs%d, blk %d, BB cut, Expect Cut; ns",kine[k],b ), 200, -20, 180 );  
    }
  }
  
  //Outfiles
  TFile *fout = new TFile(outputfilename,"RECREATE");

  //Accounting and diagnostic variables
  Long64_t nexpect[nkine] = {0};
  Long64_t nphit[nkine] = {0};
  Long64_t nrhit[nkine] = {0};
  Long64_t nrchit[nkine] = {0};
  
  cout << endl << "All chains set up. Proceeding to analysis." << endl;

  for( Int_t k=0; k<nkine; k++ ){
    
    //cout << C[k]->GetEntries() << endl;

    // Set variables to keep track of event sequence for ttreeformula
    //Int_t treenum = 0, currenttreenum = 0;

    Double_t pBeam = E_e[k]/(1.+E_e[k]/M_p*(1.-cos(BB_th[k]*TMath::DegToRad()))); // Momentum of beam calculated from theta
    Double_t Eloss_outgoing = celldiameter/2.0/sin(BB_th[k]*TMath::DegToRad()) * rho_tgt * dEdx_tgt; //Approximately 1 MeV, could correct further with raster position


    cout << endl << "Processing events on kinematic " << kine[k] << ".." << endl << endl;

    for( Long64_t nevent = 1; nevent <Nevents[k]; nevent++){

      if ( nevent%10000==0 ) cout << "Kinematic: " << kine[k] << ", entry: " << nevent << "/" << Nevents[k] << ", tdc miss ratio: " << nphit[k] << "/" << nexpect[k] << " \r";
      cout.flush();

      C[k]->GetEntry( elist[k]->GetEntry(nevent) );

      //Correct the beam energy with energy loss in target using vertex position
      Double_t Eloss = (BBtr_vz[k][0]+l_tgt/2.0) * rho_tgt * dEdx_tgt + uwallthick_LH2 * rho_Al * dEdx_Al; //approximately 3 MeV
      Double_t E_corr = E_e[k] - Eloss;

      //Physics Calculations
      Double_t p_corr = BBtr_p[k][0] - Eloss_outgoing; //Neglecting the mass of e'
      Double_t etheta = acos( BBtr_pz[k][0]/BBtr_p[k][0] ); //Will use the uncorrected track momentum to reconstruct e' theta
      Double_t ephi = atan2( BBtr_py[k][0], BBtr_px[k][0] );   
      TVector3 vertex( 0, 0, BBtr_vz[k][0] ); // z location of vertex in hall coordinates
      TLorentzVector Pbeam( 0, 0, E_corr, E_corr ); //Mass of e negligable
      TLorentzVector kprime( BBtr_px[k][0], BBtr_py[k][0], BBtr_pz[k][0], BBtr_p[k][0] );
      TLorentzVector Ptarg( 0, 0, 0, M_p );   
      TLorentzVector q = Pbeam - kprime;
      TLorentzVector PgammaN = Ptarg + q; //(-px, -py, ebeam - pz, Mp + ebeam - p)     
      Double_t pel = E_corr/( 1.+E_corr/M_p*( 1.-cos(etheta) ) );
      Double_t nu = E_corr - BBtr_p[k][0];
      Double_t pp = sqrt( pow(nu,2)+2.*M_p*nu );
      Double_t phinucleon = ephi + PI; //assume coplanarity
      Double_t thetanucleon = acos( (E_corr - BBtr_pz[k][0])/pp ); //use elastic constraint on nucleon kinematics     
      TVector3 pNhat( sin(thetanucleon)*cos(phinucleon), sin(thetanucleon)*sin(phinucleon), cos(thetanucleon) );
      
      //Define HCal coordinate system
      TVector3 HCAL_zaxis( sin(-HCal_th[k]*TMath::DegToRad()), 0, cos(-HCal_th[k]*TMath::DegToRad()) );
      TVector3 HCAL_xaxis( 0, -1, 0 );
      TVector3 HCAL_yaxis = HCAL_zaxis.Cross(HCAL_xaxis).Unit();      
      TVector3 HCAL_origin = HCal_d[k] * HCAL_zaxis + HCALHeight * HCAL_xaxis;     
      TVector3 HCALpos = HCAL_origin + hcalx[k] * HCAL_xaxis + hcaly[k] * HCAL_yaxis;
      Double_t sintersect = ( HCAL_origin - vertex ).Dot( HCAL_zaxis ) / (pNhat.Dot( HCAL_zaxis ) ); 
      TVector3 HCAL_intersect = vertex + sintersect * pNhat;   

      //Calculate quantities of interest
      Double_t yexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_yaxis );
      Double_t xexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_xaxis ); 
      Double_t dx = hcalx[k] - xexpect_HCAL;
      Double_t dy = hcaly[k] - yexpect_HCAL;
      Double_t E_ep = sqrt( pow(M_e,2) + pow(BBtr_p[k][0],2) ); // Obtain the scattered electron energy 
      Double_t p_ep = BBtr_p[k][0];     
      Double_t Q2 = 2*E_corr*E_ep*( 1-(BBtr_pz[k][0]/p_ep) ); // Obtain Q2 from beam energy, outgoing electron energy, and momenta  
	
      Double_t W2 = kineW2[k];

      //Record best tdc time per event
      Int_t pblkid = (int)(cid[k][0]-1); //Must correct for index switch on blocks

      //Check blk index
      hblkID[k]->Fill( pblkid );
	
      //Fill histograms without any cuts
      hW2_nocut[k]->Fill( W2 );

      //perform w2 cut to tighten elastic selection
      bool elastic = abs( W2 - W2_mean[k] ) < W2_sig[k];
	
      //check if the BB track would predict a hit on HCal. Account for SBS field in x.
      bool exphcalactive = false;
      if(  yexpect_HCAL < (posHCalYf_MC-HCalblk_l_h_MC) &&
	   yexpect_HCAL > (posHCalYi_MC+HCalblk_l_h_MC) &&
	   xexpect_HCAL-magoff[k] < (posHCalXf_MC-HCalblk_l_v_MC) &&
	   xexpect_HCAL-magoff[k] > (posHCalXi_MC+HCalblk_l_v_MC)){
	  
	exphcalactive=true;
      }

      //check if HCal ADC saw a hit and resultant cluster exists in active area (avoid periphery)
      bool hcalactive = false;
      if( crow[k][0] > 0 &&
	  crow[k][0] < 23 &&
	  ccol[k][0] > 0 &&
	  ccol[k][0] < 11 ){
	hcalactive = true;
      }

      //Cut all events where no good tdc is expected in hcal and a cluster exists in HCal from fADC (BBcut)
      if( elastic == false || 
	  exphcalactive == false || 
	  hcalactive == false 
	  ) 
	continue;

      nexpect[k]++;

      //Fill histograms with BigBite and HCal active area cuts
      //Get difference in time between adc time and tdc time for reference and later t cuts on blk
      Double_t atdiff = catime[k][0]-ctdc[k][0];
      hdx_BBcut[k]->Fill( dx );
      hdy_BBcut[k]->Fill( dy );
      httime_blkID[k]->Fill( pblkid, ctdc[k][0] );
      hatdiff_blkID[k]->Fill( pblkid, atdiff );

      // Record primary block tdc detections. This includes reconstructed hits for efficiency calc
      if( ctdc[k][0]>-500&&ctdc[k][0]<500 ){
	httime_blkN_BBcut[k][0]->Fill( ctdc[k][0] );
	hatime_blkN_BBcut[k][0]->Fill( catime[k][0] );
	nphit[k]++;
	nrhit[k]++;
	nrchit[k]++;
      }else{ // See if time can be reconstructed from additional blocks in cluster
	bool ghit = false;
	bool gchit = false;
	for( Int_t b=0; b<nblk[k]; b++ ){ //Blocks are sorted in descending energy
	  Double_t atdiffblk = catime[k][b]-ctdc[k][b];
	  httime_blkN_BBcut[k][b]->Fill( ctdc[k][b] );
	  hatime_blkN_BBcut[k][b]->Fill( catime[k][b] );

	  if( ctdc[k][b]>-500&&
	      ctdc[k][b]<500&&
	      ghit==false ){
	    nrhit[k]++;
	    ghit=true;
	    //continue;
	  }

	  if( ctdc[k][b]>-500&&
	      ctdc[k][b]<500&&
	      gchit==false&&
	      abs( atdiffblk - atdiffc[k] )<atdiffs[k] ){
	    nrchit[k]++;
	    gchit=true;
	    //continue;
	  }
	}
      }
    }

    //fout->Write();
    //cout << nphit[k] << endl;
  }
  
  //Write histos to outfile
  for( Int_t k=0; k<nkine; k++ ){
    hW2_nocut[k]->Write();
    hdy_BBcut[k]->Write();
    hdx_BBcut[k]->Write();
    hblkID[k]->Write();
    httime_blkID[k]->Write();
    hatdiff_blkID[k]->Write();
    for( Int_t b=0; b<maxblock; b++ ){
      if( httime_blkN_BBcut[k][b]->GetEntries()>0 ) httime_blkN_BBcut[k][b]->Write();
      if( hatime_blkN_BBcut[k][b]->GetEntries()>0 ) hatime_blkN_BBcut[k][b]->Write();
    }

  }

  Double_t TDC_peff[nkine]={0.};
  Double_t TDC_reff[nkine]={0.};
  Double_t TDC_rceff[nkine]={0.};
  Double_t xerr[nkine]={0.};
  Double_t kineD[nkine]={0.};
  Double_t pblkerr[nkine]={0.};
  Double_t reconerr[nkine]={0.};
  Double_t reconcuterr[nkine]={0.};
  
  for( Int_t k=0; k<nkine; k++ ){
    kineD[k]=(Double_t)kine[k];

    if( nexpect[k]==0 ) nexpect[k]=1; //Return zero if no events exist on the kinematic

    TDC_peff[k] = (Double_t)nphit[k] / (Double_t)nexpect[k];
    TDC_reff[k] = (Double_t)nrhit[k] / (Double_t)nexpect[k];
    TDC_rceff[k] = (Double_t)nrchit[k] / (Double_t)nexpect[k];

    //Error from binomial statistics - sqrt(n eff(1-eff))
    pblkerr[k] = sqrt(TDC_peff[k]*(1-TDC_peff[k])/(Double_t)nphit[k]);
    reconerr[k] = sqrt(TDC_reff[k]*(1-TDC_reff[k])/(Double_t)nrhit[k]);
    reconcuterr[k] = sqrt(TDC_rceff[k]*(1-TDC_rceff[k])/(Double_t)nrhit[k]);
  }

  //cout << TDC_peff[0] << " " << TDC_reff[0] << endl;

  TGraphErrors *cpeff = new TGraphErrors( nkine, kineD, TDC_peff, xerr, pblkerr );
  TGraphErrors *creff = new TGraphErrors( nkine, kineD, TDC_reff, xerr, reconerr );
  TGraphErrors *crceff = new TGraphErrors( nkine, kineD, TDC_rceff, xerr, reconcuterr );

  TGraphErrors *cpeff_nKE = new TGraphErrors( nkine, cenKE_n, TDC_peff, xerr, pblkerr );
  TGraphErrors *creff_nKE = new TGraphErrors( nkine, cenKE_n, TDC_reff, xerr, reconerr );
  TGraphErrors *crceff_nKE = new TGraphErrors( nkine, cenKE_n, TDC_rceff, xerr, reconcuterr );

  TGraphErrors *cpeff_bE = new TGraphErrors( nkine, E_e, TDC_peff, xerr, pblkerr );
  TGraphErrors *creff_bE = new TGraphErrors( nkine, E_e, TDC_reff, xerr, reconerr );
  TGraphErrors *crceff_bE = new TGraphErrors( nkine, E_e, TDC_rceff, xerr, reconcuterr );
  
  cpeff->SetMarkerStyle(20);
  cpeff->SetTitle("TDC efficiency, primary cluster, primary block; kinematic; detected/expected");
  cpeff->Draw("AC*");
  cpeff->Write("tdceff_pblk");

  creff->SetMarkerStyle(21);
  creff->SetTitle("TDC efficiency, primary cluster, first block w/data; kinematic; detected/expected");
  creff->Draw("CP");
  creff->Write("tdceff_recon");

  crceff->SetMarkerStyle(22);
  crceff->SetTitle("TDC efficiency, primary cluster, first block w/data, adct-tdc cut; kinematic; detected/expected");
  crceff->Draw("CP");
  crceff->Write("tdceff_reconcut");

  cpeff_nKE->SetMarkerStyle(20);
  cpeff_nKE->SetTitle("TDC efficiency, primary block; nucleon KE (GeV); detected/expected");
  cpeff_nKE->Draw("AC*");
  cpeff_nKE->Write("tdceff_pblk_nKE");

  creff_nKE->SetMarkerStyle(21);
  creff_nKE->SetTitle("TDC efficiency, primary cluster, first block w/data; nucleon KE (GeV); detected/expected");
  creff_nKE->Draw("CP");
  creff_nKE->Write("tdceff_recon_nKE");

  crceff_nKE->SetMarkerStyle(22);
  crceff_nKE->SetTitle("TDC efficiency, primary cluster, first block w/data, adct-tdc cut; nucleon KE (GeV); detected/expected");
  crceff_nKE->Draw("CP");
  crceff_nKE->Write("tdceff_reconcut_nKE");

  cpeff_bE->SetMarkerStyle(20);
  cpeff_bE->SetTitle("TDC efficiency, primary block; beam E (GeV); detected/expected");
  cpeff_bE->Draw("AC*");
  cpeff_bE->Write("tdceff_pblk_bE");

  creff_bE->SetMarkerStyle(21);
  creff_bE->SetTitle("TDC efficiency, primary cluster, first block w/data; beam E (GeV); detected/expected");
  creff_bE->Draw("CP");
  creff_bE->Write("tdceff_recon_bE");

  crceff_bE->SetMarkerStyle(22);
  crceff_bE->SetTitle("TDC efficiency, primary cluster, first block w/data, adct-tdc cut; beam E (GeV); detected/expected");
  crceff_bE->Draw("CP");
  crceff_bE->Write("tdceff_reconcut_bE");

  //c1->Write();

  for( Int_t k=0; k<nkine; k++ ){
    
    cout << "TDC primary detection efficiency for kin " << kine[k] << ": " << TDC_peff[k] << " +/- " << pblkerr[k] << endl;

  }

  cout << endl << endl;

  for( Int_t k=0; k<nkine; k++ ){
    
    cout << "TDC reconstructed detection efficiency for kin " << kine[k] << ": " << TDC_reff[k] << " +/- " << reconerr[k] << endl;

  }

  cout << endl << endl;

  for( Int_t k=0; k<nkine; k++ ){
    
    cout << "TDC reconstructed with time difference cut detection efficiency for kin " << kine[k] << ": " << TDC_rceff[k] << " +/- " << reconcuterr[k] << endl;

  }

  fout->Write();

  st->Stop();

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

}
