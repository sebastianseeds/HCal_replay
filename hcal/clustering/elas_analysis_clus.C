//SSeeds 8.23.22 - Post-production - Script written to evaluate performance of clustering and evaluate how to cut effectively on clusters for elastic events

#include <ctime>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iostream>
#include <TSystem.h>
#include <TStopwatch.h>
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"

TChain *T = 0;

const Int_t kNcell = 288; // Total number of HCal modules
const Int_t kNrows = 24; // Total number of HCal rows
const Int_t kNcols = 12; // Total number of HCal columns
const Int_t maxTracks = 1000; // Reasonable limit on tracks to be stored per event
const Int_t maxTdcChan = 10; // Set to accomodate original 5 TDCTrig channels with buffer

const Double_t PI = TMath::Pi();
const Double_t M_e = 0.00051;
const Double_t M_p = 0.938272;
const Double_t M_n = 0.939565;
const Double_t c_light = 299792458.0;

Double_t hcalheight = 0.365; //m The height of the center of HCAL above beam

string getDate();

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

void elas_analysis_clus( Int_t run = -1, Double_t tmax = 20 ){
  
  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  //Declare branch variables
  Double_t cid[kNcell] = {0};
  Double_t cx[kNcell] = {0};
  Double_t cy[kNcell] = {0};
  Double_t crow[kNcell] = {0};
  Double_t ccol[kNcell] = {0};
  Double_t cbrow[kNcell] = {0};
  Double_t cbcol[kNcell] = {0};
  Double_t ce[kNcell] = {0};
  Double_t ceblk[kNcell] = {0};
  Double_t cnblk[kNcell] = {0};
  Double_t cblke[kNcell] = {0};
  Double_t cblkatime[kNcell] = {0};
  Double_t catime[kNcell] = {0};
  Double_t ctdctime[kNcell] = {0};
  Double_t nclus; 
  Double_t cblkid[kNcell] = {0};
  Double_t BBtr_p[maxTracks], BBtr_px[maxTracks], BBtr_py[maxTracks], BBtr_pz[maxTracks];
  Double_t BBtr_vz[maxTracks], BBtr_chi2[maxTracks];
  Double_t BBtr_n, BBps_x, BBps_y, BBps_e, BBsh_x, BBsh_y, BBsh_e;
  Double_t HCALx, HCALy, HCALe;

  Double_t TDCT_id[maxTdcChan], TDCT_tdc[maxTdcChan]; 
  Int_t TDCTndata;

  if( !T ) {
    T = new TChain( "T" );

    T->Add(Form("/volatile/halla/sbs/seeds/rootfiles/hcal_general_%d*",run));

    T->SetBranchStatus( "*", 0 ); //Turn off all branches for faster execution
    //Select only necessary branches for processing
    T->SetBranchStatus( "sbs.hcal.x", 1 );
    T->SetBranchStatus( "sbs.hcal.y", 1 );
    T->SetBranchStatus( "sbs.hcal.e", 1 );
    T->SetBranchStatus( "sbs.hcal.clus.id", 1 );
    T->SetBranchStatus( "sbs.hcal.clus.row", 1 );
    T->SetBranchStatus( "sbs.hcal.clus.col", 1 );
    T->SetBranchStatus( "sbs.hcal.clus.atime", 1 );
    T->SetBranchStatus( "sbs.hcal.clus.tdctime", 1 );
    T->SetBranchStatus( "sbs.hcal.clus_blk.atime", 1 );
    T->SetBranchStatus( "sbs.hcal.clus_blk.row", 1 );
    T->SetBranchStatus( "sbs.hcal.clus_blk.col", 1 );
    T->SetBranchStatus( "sbs.hcal.clus.e", 1 );
    T->SetBranchStatus( "sbs.hcal.clus.x", 1 );
    T->SetBranchStatus( "sbs.hcal.clus.y", 1 );
    T->SetBranchStatus( "sbs.hcal.clus.eblk", 1 );
    T->SetBranchStatus( "sbs.hcal.clus.nblk", 1 );
    T->SetBranchStatus( "sbs.hcal.nclus", 1 );
    T->SetBranchStatus( "sbs.hcal.clus_blk.id", 1 );
    T->SetBranchStatus( "sbs.hcal.clus_blk.e", 1 );
    T->SetBranchStatus( "bb.tr.chi2", 1 );
    T->SetBranchStatus( "bb.tr.n", 1 );
    T->SetBranchStatus( "bb.tr.px", 1 );
    T->SetBranchStatus( "bb.tr.py", 1 );
    T->SetBranchStatus( "bb.tr.pz", 1 );    
    T->SetBranchStatus( "bb.tr.vz", 1 );
    T->SetBranchStatus( "bb.tr.p", 1 );
    T->SetBranchStatus( "bb.ps.e", 1 );
    T->SetBranchStatus( "bb.ps.x", 1 );
    T->SetBranchStatus( "bb.ps.y", 1 );
    T->SetBranchStatus( "bb.sh.e", 1 );
    T->SetBranchStatus( "bb.sh.x", 1 );
    T->SetBranchStatus( "bb.sh.y", 1 );
    T->SetBranchStatus( "bb.tdctrig.tdc", 1 );
    T->SetBranchStatus( "bb.tdctrig.tdcelemID", 1 );
    T->SetBranchStatus( "Ndata.bb.tdctrig.tdcelemID", 1 );
    
    //Set branches to variables
    T->SetBranchAddress( "sbs.hcal.x", &HCALx );
    T->SetBranchAddress( "sbs.hcal.y", &HCALy );
    T->SetBranchAddress( "sbs.hcal.e", &HCALe );
    T->SetBranchAddress( "sbs.hcal.clus.id", cid );
    T->SetBranchAddress( "sbs.hcal.clus.x", cx );
    T->SetBranchAddress( "sbs.hcal.clus.y", cy );
    T->SetBranchAddress( "sbs.hcal.clus.row", crow );
    T->SetBranchAddress( "sbs.hcal.clus.col", ccol );
    T->SetBranchAddress( "sbs.hcal.clus.atime", catime ); //a_time for main block in each cluster
    T->SetBranchAddress( "sbs.hcal.clus.tdctime", ctdctime ); //a_time for main block in each cluster
    T->SetBranchAddress( "sbs.hcal.clus_blk.atime", cblkatime ); //a_time for each block in main cluster
    T->SetBranchAddress( "sbs.hcal.clus_blk.row", cbrow );
    T->SetBranchAddress( "sbs.hcal.clus_blk.col", cbcol );
    T->SetBranchAddress( "sbs.hcal.clus.e", ce );
    T->SetBranchAddress( "sbs.hcal.clus.eblk", ceblk ); // Highest energy block
    T->SetBranchAddress( "sbs.hcal.clus.nblk", cnblk ); //number of blocks in cluster per event
    T->SetBranchAddress( "sbs.hcal.nclus", &nclus ); //number of clusters
    T->SetBranchAddress( "sbs.hcal.clus_blk.id", cblkid );
    T->SetBranchAddress( "sbs.hcal.clus_blk.e", cblke ); // Array of block energies
    T->SetBranchAddress( "bb.tr.n", &BBtr_n );
    T->SetBranchAddress( "bb.tr.chi2", BBtr_chi2 );
    T->SetBranchAddress( "bb.tr.px", BBtr_px );
    T->SetBranchAddress( "bb.tr.py", BBtr_py );
    T->SetBranchAddress( "bb.tr.pz", BBtr_pz );
    T->SetBranchAddress( "bb.tr.vz", BBtr_vz );
    T->SetBranchAddress( "bb.tr.p", BBtr_p );
    T->SetBranchAddress( "bb.ps.e", &BBps_e );
    T->SetBranchAddress( "bb.ps.x", &BBps_x );
    T->SetBranchAddress( "bb.ps.y", &BBps_y );
    T->SetBranchAddress( "bb.sh.e", &BBsh_e );
    T->SetBranchAddress( "bb.sh.x", &BBsh_x );
    T->SetBranchAddress( "bb.sh.y", &BBsh_y );
    T->SetBranchAddress( "bb.tdctrig.tdcelemID", TDCT_id );
    T->SetBranchAddress( "bb.tdctrig.tdc", TDCT_tdc );
    T->SetBranchAddress( "Ndata.bb.tdctrig.tdcelemID", &TDCTndata );  

    cout << endl << "Opened up tree from run " << run << " with nentries = " << T->GetEntries() << "." << endl;
  }
 
  cout << endl << "Maximum time between blocks in same cluster from db_sbs.hcal.dat = " << tmax << "." << endl;

  // Define a clock to check macro processing time
  TStopwatch *sw = new TStopwatch();
  sw->Start();

  // Get the date
  string date = getDate();
  
  // Declare outfile
  TFile *fout = new TFile( "outfiles/elas_analysis4_clus_out.root", "RECREATE" );
  /*
  // Declare general vars (SBS4)
  Double_t E_e = 3.728; //GeV
  Double_t HCal_d = 11; //m
  Double_t HCal_th = 31.9 * TMath::DegToRad(); //rad
  Double_t W_mean = 0.97; //GeV
  Double_t W_sig = 0.2; //GeV
  */
  
  // Declare general vars (SBS8)
  Double_t E_e = 5.965; //GeV
  Double_t HCal_d = 11; //m
  Double_t HCal_th = 29.4 * TMath::DegToRad(); //rad
  Double_t W_mean = 0.875; //GeV
  Double_t W_sig = 0.2; //GeV
  Double_t Ycut = 0.0; //m
  Double_t Ysig = 0.4; //m
  Double_t Xcut = 0.73; //m
  Double_t Xsig = 0.35; //m

  // Other general vars
  Int_t Ntot = 0;
  Int_t Nelas = 0;

  // Initialize global histograms
  TH1D *hTDC_clu_tcut = new TH1D("hTDC_clu_tcut", "TDC Time, Highest E Cluster after TDC Cut; ns", 400, -300, 100 );
  TH1D *hADCt_clu_tcut = new TH1D("hADCt_clu_tcut", "ADCt Time, Highest E Cluster after ADCt Cut; ns", 200, -20, 180 );
  TH1D *hTDC_clu_tcutelse = new TH1D("hTDC_clu_tcutelse", "TDC Time, All Other Clusters after TDC Cut; ns", 400, -300, 100 );
  TH1D *hADCt_clu_tcutelse = new TH1D("hADCt_clu_tcutelse", "ADCt Time, All Other Clusters after ADCt Cut; ns", 200, -20, 180 );

  TH1D *hClusE_nocut = new TH1D( "ClusE_nocut", "Primary Cluster E (no cuts); GeV", 100, 0., 1. );
  TH1D *hNblk_nocut = new TH1D( "hNblk_nocut", "Blocks in Primary Cluster (no cuts)", 20, 0, 20 );
  TH1D *hNclus_nocut = new TH1D( "hNclus_nocut", "Total Clusters (no cuts)", 20, 0, 20 );

  TH1D *hTDC_clu0_nocut = new TH1D("hTDC_clu0_nocut", "TDC Time, Primary Cluster (no cuts); ns", 400, -300, 100 );
  TH1D *hTDC_clu1_nocut = new TH1D("hTDC_clu1_nocut", "TDC Time, Cluster 2 (no cuts); ns", 400, -300, 100 );
  TH1D *hTDC_clu2_nocut = new TH1D("hTDC_clu2_nocut", "TDC Time, Cluster 3 (no cuts); ns", 400, -300, 100 );
  TH1D *hTDC_clu3_nocut = new TH1D("hTDC_clu3_nocut", "TDC Time, Cluster 4 (no cuts); ns", 400, -300, 100 );
  TH1D *hTDC_clu4_nocut = new TH1D("hTDC_clu4_nocut", "TDC Time, Cluster 5 (no cuts); ns", 400, -300, 100 );
  TH1D *hADCt_clu0_nocut = new TH1D("hADCt_clu0_nocut", "ADCt Time, Primary Cluster (no cuts); ns", 200, -20, 180 );
  TH1D *hADCt_clu1_nocut = new TH1D("hADCt_clu1_nocut", "ADCt Time, Cluster 2 (no cuts); ns", 200, -20, 180 );
  TH1D *hADCt_clu2_nocut = new TH1D("hADCt_clu2_nocut", "ADCt Time, Cluster 3 (no cuts); ns", 200, -20, 180 );
  TH1D *hADCt_clu3_nocut = new TH1D("hADCt_clu3_nocut", "ADCt Time, Cluster 4 (no cuts); ns", 200, -20, 180 );
  TH1D *hADCt_clu4_nocut = new TH1D("hADCt_clu4_nocut", "ADCt Time, Cluster 5 (no cuts); ns", 200, -20, 180 );

  TH1D *hTDC_clu0 = new TH1D("hTDC_clu0", "TDC Time, Primary Cluster; ns", 400, -300, 100 );
  TH1D *hTDC_clu1 = new TH1D("hTDC_clu1", "TDC Time, Cluster 2; ns", 400, -300, 100 );
  TH1D *hTDC_clu2 = new TH1D("hTDC_clu2", "TDC Time, Cluster 3; ns", 400, -300, 100 );
  TH1D *hTDC_clu3 = new TH1D("hTDC_clu3", "TDC Time, Cluster 4; ns", 400, -300, 100 );
  TH1D *hTDC_clu4 = new TH1D("hTDC_clu4", "TDC Time, Cluster 5; ns", 400, -300, 100 );
  TH1D *hADCt_clu0 = new TH1D("hADCt_clu0", "ADCt Time, Primary Cluster; ns", 200, -20, 180 );
  TH1D *hADCt_clu1 = new TH1D("hADCt_clu1", "ADCt Time, Cluster 2; ns", 200, -20, 180 );
  TH1D *hADCt_clu2 = new TH1D("hADCt_clu2", "ADCt Time, Cluster 3; ns", 200, -20, 180 );
  TH1D *hADCt_clu3 = new TH1D("hADCt_clu3", "ADCt Time, Cluster 4; ns", 200, -20, 180 );
  TH1D *hADCt_clu4 = new TH1D("hADCt_clu4", "ADCt Time, Cluster 5; ns", 200, -20, 180 );

  TH1D *hMaxClusT = new TH1D( "maxClusT", Form("maxClusT tmax%f",tmax), 200, -20, 180 );
  TH2D *hMaxClusTvC = new TH2D( "maxClusTvC", Form("maxClusTvC tmax%f",tmax), 288, 0, 288, 200, -20, 180 );
  TH1D *hClusE = new TH1D( "ClusE", "Primary Cluster E; GeV", 100, 0., 1. );
  TH1D *hNblk = new TH1D( "hNblk", "Blocks in Primary Cluster", 20, 0, 20 );
  TH1D *hNclus = new TH1D( "hNclus", "Total Clusters", 20, 0, 20 );
  TH1D *hClusBlkE = new TH1D( "ClusBlkE", "ClusBlkE", 100, 0., 1. );
  TH1D *hClusPrimeBlkE = new TH1D( "ClusPrimeBlkE", "ClusPrimeBlkE", 100, 0., 1. );
  TH2D *hXY_pclus = new TH2D("hXY_pclus",";y_{HCAL} (m); x_{HCAL} (m)", 250, -5.0, 5.0, 250, -10, 10 );
  TH2D *hXY = new TH2D("hXY",";y_{expect} (m); x_{expect} (m)", 250, -5.0, 5.0, 250, -10, 10 );
  TH2D *hdxdy_HCAL = new TH2D("hdxdy_HCAL",";y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)", 250, -2.0, 2.0, 500, -4, 4 );
  TH1D *hdx_HCAL = new TH1D("hdx_HCAL",";x_{HCAL}-x_{expect} (m)",500,-4,4);
  TH1D *hdy_HCAL = new TH1D("hdy_HCAL",";y_{HCAL}-y_{expect} (m)",250,-2,2);
  TH1D *hE_ep = new TH1D( "Scattered Electron Energy","E_ep;GeV", 500, 0.0, E_e*1.5 ); 
  TH1D *hKE_p = new TH1D( "Scattered Proton Kinetic Energy", "KE_pp;GeV", 500, 0.0, E_e*1.5 );
  TH1D *hTrigDiff = new TH1D( "hTrigDiff","HCal time - BBCal time (ns)", 1300, -500, 800 );
  TH1D *hW = new TH1D( "W", "W;GeV", 250, 0.3, 1.5 );
  TH1D *hQ2 = new TH1D( "Q2", "Q2;GeV2", 250, 0.5, 3.0 );
  TH1D *hE_pp = new TH1D( "Scattered Proton Energy", "E_pp;GeV", 500, 0.0, E_e*1.5 );

  // Initialize cell histograms
  /*  
  TH1D *hMaxClusTpC[kNcell];
  TH1D *hClusBlkEpC[kNcell];

  for( Int_t c=0; c<kNcell+1; c++ ){  
    hMaxClusTpC[c] = new TH1D( Form( "ADC Time Max Diff Cell%d", c ), Form( "ADC Time Max Diff Cell%d tmax%f", c, tmax ), 200, -100, 100 );
    hClusBlkEpC[c] = new TH1D( Form( "Cluster Block E Cell%d", c ), Form( "Cluster Block E Cell%d", c ), 100, 0., 1. );
  }
  */
  Long64_t Nevents = T->GetEntries();
  cout << "Total events in tree: " << Nevents << ".." << endl << endl;
  
  //Loop over events
  cout << "Main loop over all data commencing.." << endl << endl;
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
	timeremains = timekeeper*( Double_t(Nevents)/Double_t(nevent) - 1. ); 
      sw->Reset();
      sw->Continue();
      
      progress = (Double_t)((nevent+1.)/Nevents);
      cout << "] " << Int_t(progress*100.) << "%, elastics: " << Nelas << "/" << Ntot << ", time remaining: " << Int_t(timeremains/60.) << "m \r";
      cout.flush();
      
      T->GetEntry( nevent ); 
      
      // Fill histograms without any cuts
      hClusE_nocut->Fill(ce[0]);
      hNblk_nocut->Fill(cnblk[0]);
      hNclus->Fill(nclus);

      for( Int_t c=0; c<nclus; c++ ){
	if( c==0 ) {
	  hTDC_clu0_nocut->Fill( ctdctime[0] );
	  hADCt_clu0_nocut->Fill( catime[0] );
	}
	if( c==1 ) {
	  hTDC_clu1_nocut->Fill( ctdctime[1] );
	  hADCt_clu1_nocut->Fill( catime[1] );
	}
	if( c==2 ) {
	  hTDC_clu2_nocut->Fill( ctdctime[2] );
	  hADCt_clu2_nocut->Fill( catime[2] );
	}
	if( c==3 ) {
	  hTDC_clu3_nocut->Fill( ctdctime[3] );
	  hADCt_clu3_nocut->Fill( catime[3] );
	}
	if( c==4 ) {
	  hTDC_clu4_nocut->Fill( ctdctime[4] );
	  hADCt_clu4_nocut->Fill( catime[4] );
	}
      }

      //ELASTIC CUT//////
	
      Double_t etheta = acos( BBtr_pz[0]/BBtr_p[0] );
      Double_t ephi = atan2( BBtr_py[0], BBtr_px[0] );

      TVector3 vertex(0,0,BBtr_vz[0]); // z location of vertex in hall coordinates
      TLorentzVector Pbeam(0,0,E_e,E_e); //Mass of e negligable
      TLorentzVector kprime(BBtr_px[0],BBtr_py[0],BBtr_pz[0],BBtr_p[0]);
      TLorentzVector Ptarg(0,0,0,M_p);

      TLorentzVector q = Pbeam - kprime;
      TLorentzVector PgammaN = Ptarg + q; //(-px, -py, ebeam - pz, Mp + ebeam - p)

      Double_t pel = E_e/(1.+E_e/M_p*(1.-cos(etheta)));
      Double_t nu = E_e - BBtr_p[0];
      Double_t pp = sqrt(pow(nu,2)+2.*M_p*nu);
      Double_t phinucleon = ephi + TMath::Pi(); //assume coplanarity
      Double_t thetanucleon = acos( (E_e - BBtr_pz[0])/pp ); //use elastic constraint on nucleon kinematics
	
      TVector3 pNhat( sin(thetanucleon)*cos(phinucleon),sin(thetanucleon)*sin(phinucleon),cos(thetanucleon));

      //Define HCal coordinate system
      TVector3 HCAL_zaxis(sin(-HCal_th),0,cos(-HCal_th));
      TVector3 HCAL_xaxis(0,-1,0);
      TVector3 HCAL_yaxis = HCAL_zaxis.Cross(HCAL_xaxis).Unit();
	
      TVector3 HCAL_origin = HCal_d * HCAL_zaxis + hcalheight * HCAL_xaxis;

      //Plot 2D histo of position
      hXY_pclus->Fill( cy[0], cx[0] );

      //Define intersection points for hadron vector
      Double_t sintersect = ( HCAL_origin - vertex ).Dot( HCAL_zaxis ) / ( pNhat.Dot( HCAL_zaxis ) );

      TVector3 HCAL_intersect = vertex + sintersect * pNhat;

      Double_t yexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_yaxis );
      Double_t xexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_xaxis );

      //Calculate the proton spot - use for cut later on
      //hdxdy_HCAL->Fill( HCALy - yexpect_HCAL, HCALx - xexpect_HCAL );

      hXY->Fill( yexpect_HCAL, xexpect_HCAL );

      Double_t E_ep = sqrt( pow(M_e,2) + pow(BBtr_p[0],2) ); // Obtain the scattered electron energy
      hE_ep->Fill( E_ep ); // Fill histogram
	
      Double_t p_ep = BBtr_p[0];
      Double_t Q2 = 2*E_e*E_ep*( 1-(BBtr_pz[0]/p_ep) ); // Obtain Q2 from beam energy, outgoing electron energy, and momenta
      hQ2->Fill( Q2 ); // Fill histogram
	
      Double_t W = PgammaN.M();
      hW->Fill( W );
	
      //Use the electron kinematics to predict the proton momentum assuming elastic scattering on free proton at rest:
      Double_t E_pp = nu+M_p; // Get energy of the proton
      Double_t Enucleon = sqrt(pow(pp,2)+pow(M_p,2)); // Check on E_pp, same
      hE_pp->Fill( E_pp ); // Fill histogram

      Double_t KE_p = nu; //For elastics
      hKE_p->Fill( KE_p );   

      //Cut on BBCal and HCal trigger coincidence
      Double_t bbcal_time=0., hcal_time=0.;
      for(Int_t ihit=0; ihit<TDCTndata; ihit++){
	if(TDCT_id[ihit]==5) bbcal_time=TDCT_tdc[ihit];
	if(TDCT_id[ihit]==0) hcal_time=TDCT_tdc[ihit];
      }
      Double_t diff = hcal_time - bbcal_time; 
      hTrigDiff->Fill( diff ); // Fill histogram
	
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
      //if( fabs(diff-510)<40 && fabs(BBtr_vz[0])<0.06 ) hW_cuts->Fill( W );

      bool elaspass = false;
      Ntot++;

      //Additional cuts to select elastics
      if( fabs(W-W_mean)<W_sig && //Observed mean W cut
	  fabs(diff-510)<10 && //Observed coincidence trigger HCal/BB
	  BBtr_n==1 && //Single track in BigBite
	  BBps_e>0.2 && //Preshower cut to remove pions
	  abs(BBtr_vz[0])<0.08 && //Track points back to vertex inside target
	  (BBps_e+BBsh_e)>1.7 //Total energy deposited in BigBite e' range
	  ) elaspass = true;

      if(elaspass == false) continue;

      bool spotCut = false;

      // Evaluate p and n on HCal
      Double_t dx = HCALx - xexpect_HCAL; // (Energy weighted center of cluster (x component)) - (straight line nucleon projection obtained from e' track HCal location (x component))
      Double_t dy = HCALy - yexpect_HCAL; // (Energy weighted center of cluster (y component)) - (straight line nucleon projection obtained from e' track HCal location (y component))
            
      //Calculate the hadron spot/s 
      hdxdy_HCAL->Fill( dy, dx );
      hdx_HCAL->Fill( dx );
      hdy_HCAL->Fill( dy );

      if( abs(dx-Xcut)<Xsig && abs(dy-Ycut)<Ysig ) spotCut = true;

      //if( spotCut==false ) continue;

      Nelas++;

      //END ELASTIC CUT//

      // Fill primary cell histograms
      hClusE->Fill(ce[0]); //First element is primary cell in cluster
      hNblk->Fill(cnblk[0]);
      hNclus->Fill(nclus);
      hClusPrimeBlkE->Fill(cblke[0]);
      
      Double_t maxDiff = 0.0;
      Int_t maxDiffID = -1;
      
      for( Int_t b=0; b<cnblk[0]; b++ ){
	
	// Data validation
	/*
	cout << "event:" << nevent << ", block: " << b << " cid: " << cid[b] << endl;
	cout << "event:" << nevent << ", block: " << b << " crow: " << crow[b] << endl;
	cout << "event:" << nevent << ", block: " << b << " ccol: " << ccol[b] << endl;
	cout << "event:" << nevent << ", block: " << b << " catime: " << catime[b] << endl;
	cout << "event:" << nevent << ", block: " << b << " cblkatime: " << cblkatime[b] << endl;
	cout << "event:" << nevent << ", block: " << b << " cbrow: " << cbrow[b] << endl;
	cout << "event:" << nevent << ", block: " << b << " cbcol: " << cbcol[b] << endl;
	cout << "event:" << nevent << ", block: " << b << " ce: " << ce[b] << endl;
	cout << "event:" << nevent << ", block: " << b << " ceblk: " << ceblk[b] << endl;
	cout << "event:" << nevent << ", block: " << b << " cnblk: " << cnblk[b] << endl;
	cout << "event:" << nevent << ", block: " << b << " nclus: " << nclus[b] << endl;
	cout << "event:" << nevent << ", block: " << b << " cblkid: " << cblkid[b] << endl;
	cout << "event:" << nevent << ", block: " << b << " cblke: " << cblke[b] << endl;
	*/

	Int_t id = cblkid[b];
	//if(id>200) cout << id << endl;
	Double_t diffADCt = fabs(cblkatime[b]-cblkatime[0]);
	
	if( diffADCt>maxDiff ){
	  maxDiff = diffADCt;
	  maxDiffID = id;
	}

	hClusBlkE->Fill(cblke[b]);
	hMaxClusT->Fill(cblkatime[b]);
	//hClusBlkEpC[id]->Fill(cblke[b]);
      }
      
      if( maxDiff!=0.0 ){
	//hMaxClusTpC[maxDiffID]->Fill(maxDiff);
	hMaxClusTvC->Fill(maxDiffID,maxDiff);
      }

      bool TDCclusAdd = false;
      bool ADCtclusAdd = false;

      for( Int_t c=0; c<nclus; c++ ){
	if( c==0 ) {
	  hTDC_clu0->Fill( ctdctime[0] );
	  hADCt_clu0->Fill( catime[0] );
	}
	if( c==1 ) {
	  hTDC_clu1->Fill( ctdctime[1] );
	  hADCt_clu1->Fill( catime[1] );
	}
	if( c==2 ) {
	  hTDC_clu2->Fill( ctdctime[2] );
	  hADCt_clu2->Fill( catime[2] );
	}
	if( c==3 ) {
	  hTDC_clu3->Fill( ctdctime[3] );
	  hADCt_clu3->Fill( catime[3] );
	}
	if( c==4 ) {
	  hTDC_clu4->Fill( ctdctime[4] );
	  hADCt_clu4->Fill( catime[4] );
	}

	bool clusterpass = false;
	double cdx = cx[c] - xexpect_HCAL;
	double cdy = cy[c] - yexpect_HCAL;
	if( abs(cdx-Xcut)<Xsig && abs(cdy-Ycut)<Ysig ) clusterpass = true;

	if( abs(ctdctime[c]+92.41)<6.4 && TDCclusAdd == false && clusterpass==true ){ //Tuned from fits to sharp peak in all clusters (-92.41 =/- 3.2 (1 sig)
	  TDCclusAdd = true;
	  hTDC_clu_tcut->Fill( ctdctime[c] );
	}else{
	  hTDC_clu_tcutelse->Fill( ctdctime[c] );
	}

	if( abs(catime[c]-46.12)<7.8 && ADCtclusAdd == false && clusterpass==true ){ //Tuned from fits to sharp peak in all clusters (46.12 +/- 3.9 (1 sig)
	  ADCtclusAdd = true;
	  hADCt_clu_tcut->Fill( catime[c] );
	}else{
	  hADCt_clu_tcutelse->Fill( catime[c] );
	}

      }
    }
  }
  
  fout->Write();

  cout << endl << endl << "Elastic rate: " << (Double_t)Nelas/Ntot << endl;

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;
  
}

