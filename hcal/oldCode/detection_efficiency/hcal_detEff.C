//SSeeds 11.03.21 - Production - Code designed with help from A.Puckett to project scattered proton during zero-field running to face of HCal and determine detection efficiency. Currently measures relative to BB elastic track detection and sets a limit on good hits in HCal at 2*sigma the proton peak at zero field.
//Updates 10.6.22 - sseeds

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
#include "hcal.h"
#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <TChain.h>
#include <TSystem.h>
#include <TStopwatch.h>
#include "TMath.h"
#include "hcal.h"

const Int_t kNcell = 288; // Total number of HCal modules
const Int_t kNrows = 24; // Total number of HCal rows
const Int_t kNcols = 12; // Total number of HCal columns
const Double_t kdBlock = 0.15; // Width and height of each module from SBS4 DB
const Double_t hcaloffset = 0.365; 
const Int_t maxTdcChan = 10;
const Int_t maxTracks = 100;

//const Double_t Yl = -0.825; //From data pass0 SBS4 DB
//const Double_t Yh = 0.825; //From data pass0 SBS4 DB
//const Double_t Xl = -2.09; //From data pass0 SBS4 DB
//const Double_t Xh = 1.36; //From data pass0 SBS4 DB

const Double_t Yl = -0.75; //From data pass0 SBS4 DB
const Double_t Yh = 0.75; //From data pass0 SBS4 DB
const Double_t Xl = -2.015; //From data pass0 SBS4 DB
const Double_t Xh = 1.285; //From data pass0 SBS4 DB

const Double_t PI = TMath::Pi();
const Double_t Mp = 0.938272;
const Double_t Mn = 0.939565;

Double_t gausFit_alt(Double_t *x, Double_t *par) {
  Double_t amp = par[0];
  Double_t loc = par[1];
  Double_t sigma = par[2];
  Double_t expo_0 = par[3];
  Double_t expo_1 = par[4];
  Double_t ADC = x[0];
  return amp * TMath::Exp(-0.5 * pow(((ADC - loc) / sigma), 2)) * TMath::Exp( expo_0 + expo_1*ADC );
}

Double_t expoDFit(Double_t *x, Double_t *par){
  Double_t amp = par[0];
  Double_t str = par[1];
  Double_t offset = par[2];
  return amp*exp(-str*x[0])+offset;
}

Double_t expoFit(Double_t *x, Double_t *par) {
  
  Double_t expo_0 = par[0];
  Double_t expo_1 = par[1];
  Double_t ADC = x[0];
  
  if ( ADC > 0.75 && ADC < 1.05 ) {
    TF1::RejectPoint();
    return 0;
  }
  
  return TMath::Exp( expo_0*ADC - expo_1 );
}

/*
Double_t gausFit_alt(Double_t *X,Double_t *p) //Single parameter fits - GOING FOR THE GOLD!
{
  Double_t fitval = p[0]*pow(X[0]/p[1],p[2]); //X normalized to lowest HV value
  return fitval;
}
*/

void hcal_detEff( const char *configfilename="shcal_detEff.cfg", Int_t run = -1 ){
  
  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );
  
  TChain *C = new TChain("R");
  
  //Declare configuration variables
  Double_t E_e = 1.92; // Energy of the beam
  Double_t emin = 0.05; // Minimum energy of a block to be considered in clusters
  Double_t HCal_d = 14.5; // Distance to HCal from scattering chamber for comm1
  Double_t HCal_th = 35.0; // Angle that the center of HCal is at  
  Double_t trigdiff = 510.0; // Observed BBCal/HCal trigger difference
  Double_t opticsCorr = 1.05; // Correction to magnitude of p_e to account for misaligned optics
  Double_t W_mean = 0.93; // Mean of W at current kinematic
  Double_t W_sig = 0.039; // Width of W at current kinematic
  Double_t BB_th = 36.0; // Angle that the center of BB is at
  Double_t efficiency_rel; // Relative detection efficiency of HCAL (elastic events detected by HCal) / (elastic events as defined by BB tracking)
  Int_t hits_elBB = 0; // Count of total elastic hits detected in BB
  Int_t hits_gHCAL = 0; // Count of all elastic events that HCal detected
  
  // Reading config file
  cout << "Opening the following files.." << endl;
  ifstream configfile(configfilename);
  TString currentline;
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      C->Add(currentline);
      cout << currentline << endl;
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
      if( skey == "E_e" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	E_e = sval.Atof();
	cout << "Beam energy: " << E_e << endl;
      }
      if( skey == "HCal_d" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	HCal_d = sval.Atof();
	cout << "HCal distance: " << HCal_d << endl;

      }
      if( skey == "HCal_th" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	HCal_th = sval.Atof();
	cout << "HCal theta: " << HCal_th << endl;
	HCal_th *= TMath::DegToRad();
  
      }
      if( skey == "trigdiff" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	trigdiff = sval.Atof();
	cout << "BBCal/HCal trigger difference: " << trigdiff << endl;

      }
      if( skey == "opticsCorr" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	opticsCorr = sval.Atof();
	cout << "Optics correction factor: " << opticsCorr << endl;

      }
      if( skey == "W_mean" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        W_mean = sval.Atof();
	cout << "W mean: " << W_mean << endl;

      }
      if( skey == "W_sig" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        W_sig = sval.Atof();
	cout << "W sigma: " << W_sig << endl;
      }
      if( skey == "BB_th" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        BB_th = sval.Atof();
	cout << "Big Bite theta: " << BB_th << endl;
	BB_th *= TMath::DegToRad();

      }
      if( skey == "emin" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        emin = sval.Atof();
	cout << "Minimum Block Energy " << emin << endl;

      }
    }
    delete tokens;
  }
  cout << endl;

  // Declare general detector and physics parameters
  Double_t TDCT_id[maxTdcChan], TDCT_tdc[maxTdcChan]; 
  Int_t TDCTndata;

  // BB params
  Double_t BBtr_p[maxTracks], BBtr_px[maxTracks], BBtr_py[maxTracks], BBtr_pz[maxTracks];
  Double_t BBtr_vx[maxTracks], BBtr_vy[maxTracks], BBtr_vz[maxTracks], BBtr_chi2[maxTracks];
  Double_t BBtr_n, BBps_x, BBps_y, BBps_e, BBsh_x, BBsh_y, BBsh_e;

  // HCal params
  Double_t HCALblktdc[kNcell], HCALblka[kNcell];
  Double_t HCALtdc[kNcell], HCALatime[kNcell], HCALa_p[kNcell];
  Double_t crow, ccol, nblk;
  Double_t cblkid[kNcell], cblke[kNcell];
  Double_t ce[kNcell], cx[kNcell], cy[kNcell], catime[kNcell], cid[kNcell], cnblk[kNcell], ctdc[kNcell];
  Double_t nclus;
  Double_t hcalx, hcaly, hcale;

  C->SetBranchStatus( "*", 0 );
  C->SetBranchStatus( "sbs.hcal.x", 1 );
  C->SetBranchStatus( "sbs.hcal.y", 1 );
  C->SetBranchStatus( "sbs.hcal.e", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.row", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.col", 1 );
  C->SetBranchStatus( "sbs.hcal.nclus", 1 );
  C->SetBranchStatus( "bb.tr.n", 1 );
  C->SetBranchStatus( "bb.tr.p", 1 );
  C->SetBranchStatus( "bb.tr.px", 1 );
  C->SetBranchStatus( "bb.tr.py", 1 );
  C->SetBranchStatus( "bb.tr.pz", 1 );
  C->SetBranchStatus( "bb.tr.vx", 1 );
  C->SetBranchStatus( "bb.tr.vy", 1 );
  C->SetBranchStatus( "bb.tr.vz", 1 );
  C->SetBranchStatus( "bb.ps.e", 1 );
  C->SetBranchStatus( "bb.ps.x", 1 );
  C->SetBranchStatus( "bb.ps.y", 1 );
  C->SetBranchStatus( "bb.sh.e", 1 );
  C->SetBranchStatus( "bb.sh.x", 1 );
  C->SetBranchStatus( "bb.sh.y", 1 );
  C->SetBranchStatus( "bb.tdctrig.tdc", 1 );
  C->SetBranchStatus( "bb.tdctrig.tdcelemID", 1 );
  C->SetBranchStatus( "Ndata.bb.tdctrig.tdcelemID", 1 );
  C->SetBranchStatus( "bb.sh.nclus", 1 );

  C->SetBranchAddress( "bb.tr.n", &BBtr_n );
  C->SetBranchAddress( "bb.tr.px", BBtr_px );
  C->SetBranchAddress( "bb.tr.py", BBtr_py );
  C->SetBranchAddress( "bb.tr.pz", BBtr_pz );
  C->SetBranchAddress( "bb.tr.vx", BBtr_vx );
  C->SetBranchAddress( "bb.tr.vy", BBtr_vy );
  C->SetBranchAddress( "bb.tr.vz", BBtr_vz );
  C->SetBranchAddress( "bb.tr.p", BBtr_p );
  C->SetBranchAddress( "bb.sh.x", &BBsh_x );
  C->SetBranchAddress( "bb.sh.y", &BBsh_y );
  C->SetBranchAddress( "bb.sh.e", &BBsh_e );
  C->SetBranchAddress( "bb.ps.x", &BBps_x );
  C->SetBranchAddress( "bb.ps.y", &BBps_y );
  C->SetBranchAddress( "bb.ps.e", &BBps_e );
  C->SetBranchAddress( "sbs.hcal.rowblk", &crow );
  C->SetBranchAddress( "sbs.hcal.colblk", &ccol );
  C->SetBranchAddress( "sbs.hcal.x", &hcalx );
  C->SetBranchAddress( "sbs.hcal.y", &hcaly );
  C->SetBranchAddress( "sbs.hcal.e", &hcale );
  C->SetBranchAddress( "bb.tdctrig.tdcelemID", TDCT_id );
  C->SetBranchAddress( "bb.tdctrig.tdc", TDCT_tdc );
  C->SetBranchAddress( "Ndata.bb.tdctrig.tdcelemID", &TDCTndata );
  
  Int_t nentries = C->GetEntries();
  cout << endl;

  cout << "Opened tree with " << nentries << " entries." << endl;
  
  //Keep a wide cut on W for allowed elastics
  Double_t Wmin_elastic = W_mean - 2.*W_sig;
  Double_t Wmax_elastic = W_mean + 2.*W_sig;
  
  //BigBite tracks
  Double_t ntrack;
  
  //Outfile
  TFile *fout = new TFile("eff.root","RECREATE");

  //Histograms
  TH1D *hdpel = new TH1D("hdpel",";p/p_{elastic}(#theta)-1;", 250, -1.0, 0.5);
  TH1D *hW = new TH1D( "hW",";W (GeV);", 400, 0.0, 4.0 );
  TH1D *hW_pCut = new TH1D( "hW_pCut",";W (GeV);", 400, 0.0, 4.0 );
  TH1D *hW_cut_emin = new TH1D( "hW_cut_emin",";W (GeV);", 400, 0., 4.0 );
  TH1D *htDiff = new TH1D( "htDiff",";TDC_{HCal}-TDC_{BBCal} (ns);", 1300, -500, 800 );

  //TH1D *hW_cut = new TH1D("hW_cut",";W (GeV);", 400,Wmin_elastic,Wmax_elastic);
  TH1D *hW_cut = new TH1D( "hW_cut","W with Cuts; W (GeV);", 1150., 0., 1.15 ); //Tuned s.t. bins/1000=range

  TH1D *hdx_HCAL = new TH1D("hdx_HCAL",";x_{HCAL}-x_{expect} (m);", 500, -Xl-2*kdBlock, Xh+2*kdBlock);
  TH1D *hX_HCAL = new TH1D("hX_HCAL",";x_{HCAL} (m);", 500, -Xl-2*kdBlock, Xh+2*kdBlock);
  TH1D *hdy_HCAL = new TH1D("hdy_HCAL",";y_{HCAL}-y_{expect} (m);", 500, -Yl-2*kdBlock, Yh+2*kdBlock);
  TH1D *hY_HCAL = new TH1D("hY_HCAL",";y_{HCAL} (m);", 500, -1.5, 1.5);
  TH1D *hdx_HCAL_emin = new TH1D("hdx_HCAL_emin",";x_{HCAL}-x_{expect} (m);", 500, -Xl-2*kdBlock, Xh+2*kdBlock);
  TH1D *hdy_HCAL_emin = new TH1D("hdy_HCAL_emin",";y_{HCAL}-y_{expect} (m);", 500, -Yl-2*kdBlock, Yh+2*kdBlock);
  TH1D *hdr_HCAL = new TH1D("hdr_HCAL",";r_{HCAL}-r_{expect} (m);", 500, 0.0, Xh-Xl);
  TH2D *hdxdy_HCAL = new TH2D("hdxdy_HCAL",";y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)", 250, -1.25, 1.25, 250, -2.5, 2.5 );  
  TH2D *hdxdy_hcaldiag = new TH2D("hdxdy_hcaldiag",";y_{HCAL} (m); x_{HCAL} (m)", 250, -Yl-2*kdBlock, Yh+2*kdBlock, 250, -Xl-2*kdBlock, Xh+2*kdBlock );  
  TH2D *hdxdy_projdiag = new TH2D("hdxdy_projdiag",";y_{expect} (m); x_{expect} (m)", 250,  -Yl-2*kdBlock, Yh+2*kdBlock, 250, -Xl-2*kdBlock, Xh+2*kdBlock );  
  TH2D *hXY_HCAL = new TH2D("hXY_HCAL",";y_{HCAL} (m); x_{HCAL} (m)", 250, -2, 2, 250, -4, 4 );
  TH2D *hxcorr_HCAL = new TH2D("hxcorr_HCAL",";x_{expect} (m);x_{HCAL} (m)", 250, -2.5, 2.5, 250, -2.5, 2.5 );
  TH2D *hycorr_HCAL = new TH2D("hycorr_HCAL",";y_{expect} (m);y_{HCAL} (m)", 250, -1.25, 1.25, 250, -1.25, 1.25);

  TH1D *hdx_HCAL_dyCut = new TH1D("hdx_HCAL_dyCut",";x_{HCAL}-x_{expect} (m);", 500, -2.5, 2.5);
  TH1D *hdx_HCAL_dyCut_emin = new TH1D("hdx_HCAL_dyCut_emin",";x_{HCAL}-x_{expect} (m);", 500, -2.5, 2.5);

  TH1D *hvz = new TH1D("hvz","",250,-0.15,0.15);

  TH2D *hdy_HCAL_vs_z = new TH2D("hdy_HCAL_vs_z","",250,-0.15,0.15,250,-1.25,1.25);
  TH2D *hdy_HCAL_vs_ptheta = new TH2D("hdy_HCAL_vs_ptheta","",250,HCal_th-0.3,HCal_th+0.3,250,-1.25,1.25);
  TH2D *hdxdy_HCAL_cut = new TH2D("hdxdy_HCAL_cut",";y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)", 250, -1.25, 1.25, 250, -2.5, 2.5 );


  long nevent = 0;
  Double_t hcalXmax = 0;
  Double_t hcalXmin = 0;
  Double_t hcalYmax = 0;
  Double_t hcalYmin = 0;
  Double_t projXmax = 0;
  Double_t projXmin = 0;
  Double_t projYmax = 0;
  Double_t projYmin = 0;

  

  cout << "Processing events.." << endl;

  Double_t progress = 0.;
  
  while( progress<1.0 ){
    
    Int_t barwidth = 70;
    Int_t step = 1;
    
    
    
    while( C->GetEntry( nevent++ ) ){
      //if( nevent % 10000 == 0 ) cout << nevent << " / " << nentries << endl;
      
      //Only consider event if tracks exist
      ntrack = BBtr_n;
      if( ntrack > 0 ){
	
	//Define vectors from track in BB
	Double_t etheta = acos( BBtr_pz[0]/BBtr_p[0] );
	Double_t ephi = atan2( BBtr_py[0], BBtr_px[0] );
	
	//Reaction origin
	TVector3 vertex(0,0,BBtr_vz[0]);
	
	TLorentzVector Pbeam(0,0,E_e,E_e);
	TLorentzVector kprime(BBtr_px[0],BBtr_py[0],BBtr_pz[0],BBtr_p[0]);
	TLorentzVector Ptarg(0,0,0,Mp);
	
	TLorentzVector q = Pbeam - kprime;
	
	TLorentzVector PgammaN = Ptarg + q; //(-px, -py, ebeam - pz, Mp + ebeam - p)
	
	Double_t pel = E_e/(1.+E_e/Mp*(1.-cos(etheta)));
	
	hdpel->Fill( BBtr_p[0]/pel - 1.0 );
	
	Double_t W = PgammaN.M();

	hW->Fill( W );
	
	hvz->Fill( vertex.Z() );
	
	//Now project to HCAL and compare to best HCAL cluster:
	//Assume neutron (straight-line), quasi-elastic kinematics:, and assume BB/HCal trigger coincidence <20ns
	
	//Cut on BBCal and HCal trigger coincidence
	Double_t bbcal_time=0., hcal_time=0.;
	for(Int_t ihit=0; ihit<TDCTndata; ihit++){
	  if(TDCT_id[ihit]==5) bbcal_time=TDCT_tdc[ihit];
	  if(TDCT_id[ihit]==0) hcal_time=TDCT_tdc[ihit];
	}
	Double_t diff = hcal_time - bbcal_time; 
	
	htDiff->Fill(diff);

	//Count the number of elastic events with acceptable tracks in BB (reduction of pions with PS)
	if( PgammaN.M() >= Wmin_elastic && 
	    PgammaN.M() <= Wmax_elastic && 
	    fabs( vertex.Z() )<=0.08 ) {
	  hits_elBB++;
	  //hW_cut->Fill(PgammaN.M());
	}
	
	//if( PgammaN.M()<1.15 &&
	//fabs( vertex.Z() )<=0.08 &&
	//  BBps_hcale>0.1 && 
	//  BBsh_hcale>0.1 ) hW_cut->Fill(PgammaN.M());
	
	//Require that track corresponds to elastic event, energy deposited in preshower, and coincidence trigger between HCal and BB
	//if( PgammaN.M() >= Wmin_elastic && PgammaN.M() <= Wmax_elastic && BBps_hcale>0.1 && BBsh_hcale>0.1 && fabs( vertex.Z() )<=0.08 && fabs(diff-trigdiff)<30 ){
	
	//Now determine projected expectation for proton from BB e' track
	Double_t nu = E_e - BBtr_p[0];
	Double_t pp = sqrt(pow(nu,2)+2.*Mp*nu);
	Double_t phinucleon = ephi + TMath::Pi(); //Coplanarity
	Double_t thetanucleon = acos( (E_e - BBtr_p[0]*cos(etheta))/pp ); //Elastic constraInt_t on nucleon kinematics
	
	//Nucleon momentum unit vector - no field
	TVector3 pNhat( sin(thetanucleon)*cos(phinucleon),sin(thetanucleon)*sin(phinucleon),cos(thetanucleon));
	
	//Get the coordinate system for HCal from the kinematic HCal theta
	TVector3 HCAL_zaxis(-sin(HCal_th),0,cos(HCal_th));
	TVector3 HCAL_xaxis(0,1,0);
	TVector3 HCAL_yaxis = HCAL_zaxis.Cross(HCAL_xaxis).Unit();
	
	//Project the origin of the detector coordinate system from the scattering chamber
	TVector3 HCAL_origin = HCal_d * HCAL_zaxis;
	
	//
	Double_t sInt_tersect = ( HCAL_origin - vertex ).Dot( HCAL_zaxis ) / ( pNhat.Dot( HCAL_zaxis ) );
	TVector3 HCAL_Int_tersect = vertex + sInt_tersect * pNhat;
	
	
	Double_t yexpect_HCAL = (HCAL_Int_tersect - HCAL_origin).Dot( HCAL_yaxis );
	Double_t xexpect_HCAL = (HCAL_Int_tersect - HCAL_origin).Dot( HCAL_xaxis );
	
	//Double_t dx = (cx[0]+1.550) - xexpect_HCAL; //Adding offset to correct for error in database file. Will adjust later.
	//Double_t dy = (cy[0]+0.323) - yexpect_HCAL; //Adding offset
	
	//Double_t dx = (cx[0]) - xexpect_HCAL;
	//Double_t dy = (cy[0]) - yexpect_HCAL;

	Double_t dx = hcalx - xexpect_HCAL;
	Double_t dy = hcaly - yexpect_HCAL;

	Double_t dr = sqrt( pow( dx, 2 )+pow( dy, 2 ));
	
	//cout << "dx=" << dx << " dy=" << dy << " dr=" << dr << endl;
	//if( chcale > emin ) hW_cut_emin->Fill(PgammaN.M());

	// Can add cut on proton spot here to populate W distribution
	// Add fiducial cut here on both BB elastics and HCal elastics
	// Fill denominator if conditions are met

	bool elastic = abs( W - W_mean ) < W_sig;

	//DENOMINATOR - Fill, fit, and extract total number of elastics from peak (all Earm)
	if(  //PgammaN.M()<1.15 &&
	     //PgammaN.M()>0. &&
	     elastic &&
	     fabs( vertex.Z() )<=0.08 &&
	     //BBps_hcale>0.1 && 
	     //BBsh_hcale>0.1 &&
	     yexpect_HCAL < Yh &&
	     yexpect_HCAL > Yl &&
	     xexpect_HCAL < Xh &&
	     xexpect_HCAL > Xl){
	  hW_cut->Fill(PgammaN.M());
	  
	  
	  
	  //if( chcale > emin ) hW_cut_emin->Fill(PgammaN.M());
	}
	/*
	if(  PgammaN.M()<1.15 &&
	     PgammaN.M()>0. &&
	     fabs( vertex.Z() )<=0.08 &&
	     BBps_hcale>0.1 && 
	     BBsh_hcale>0.1 &&
	     yexpect_HCAL < 0.885 &&
	     yexpect_HCAL > -0.885 &&
	     xexpect_HCAL < 1.785 + hcalheight &&
	     xexpect_HCAL > -1.785 + hcalheight ){
	  hW_cut->Fill(PgammaN.M());
	  
	  //if( chcale > emin ) hW_cut_emin->Fill(PgammaN.M());
	}
	*/

	//if(crow[0]==0||crow[0]==24||ccol[0]==0||ccol[0]==12) cout << crow[0] << " " << ccol[0] << endl;	    
	//cx[0] > -2.5 &&
	//cx[0] < 2.5 &&
	//cy[0] > -1.2 &&
	//cy[0] < 1.2 &&
	hdxdy_hcaldiag->Fill( dy, dx );
	
	hdxdy_projdiag->Fill( yexpect_HCAL, xexpect_HCAL );

	//cout << "BBCal-HCal tDiff: " << diff << endl;

	// Fill Numerator if conditions are met
	if( PgammaN.M() >= Wmin_elastic && 
	    PgammaN.M() <= Wmax_elastic && 
	    //BBps_hcale > 0.1 &&
	    //BBsh_hcale > 0.1 &&
	    fabs( vertex.Z() ) <= 0.08 &&
	    fabs( diff-trigdiff )<30 ){
	  
	  hdx_HCAL->Fill( dx );
	  hdy_HCAL->Fill( dy );
	  if( hcale > emin ) hdx_HCAL_emin->Fill( dx ); //Only fill these if the minimum energy of the cluster correlates well with more than a single block
	  //cout << hcale << endl;

	  if( hcale > emin ) hdy_HCAL_emin->Fill( dy );
	  hdr_HCAL->Fill( dr );
	  
	  if(dy>-0.5&&dy<0.5) {
	    if( hcale > emin ) hdx_HCAL_dyCut_emin->Fill( dx );
	    hdx_HCAL_dyCut->Fill( dx );
	  }
	  hdxdy_HCAL->Fill( dy, dx );

	  hX_HCAL->Fill( cx[0] );
	  hY_HCAL->Fill( cy[0] );
	  hXY_HCAL->Fill( cx[0], cy[0] );

	  if( crow != 0 && crow != 23 && ccol != 0 && ccol != 11 ){
	    if( hcalXmax < cx[0] ) hcalXmax = cx[0];
	    if( hcalXmin > cx[0] ) hcalXmin = cx[0];
	    if( hcalYmax < cy[0] ) hcalYmax = cy[0];
	    if( hcalYmin > cy[0] ) hcalYmin = cy[0];

	    if( projXmax < xexpect_HCAL ) projXmax = xexpect_HCAL;
	    if( projXmin > xexpect_HCAL ) projXmin = xexpect_HCAL;
	    if( projYmax < yexpect_HCAL ) projYmax = yexpect_HCAL;
	    if( projYmin > yexpect_HCAL ) projYmin = yexpect_HCAL;
	    
	  }

	  //cout << "dx dy: " << dx << " " << dy << endl;
	  
	  hxcorr_HCAL->Fill( xexpect_HCAL, cx[0] );
	  hycorr_HCAL->Fill( yexpect_HCAL, cy[0] );
	  
	  hdy_HCAL_vs_z->Fill( vertex.Z(), cy[0] - yexpect_HCAL );
	  hdy_HCAL_vs_ptheta->Fill( thetanucleon, cy[0] - yexpect_HCAL );
	  
	}
	
	//Plot another distribution of W cutting on the region which should contain elastics on HCal
	if(dy>-0.5&&dy<0.5 && dx>-0.5&&dx<0.5) {
	  hW_pCut->Fill(PgammaN.M());
	  if( hcale > emin ) hW_cut_emin->Fill(PgammaN.M());
	}
      }
      
      cout << "[";
      Int_t pos = barwidth*progress;
      for( Int_t i=0; i<barwidth; ++i ){
	if( i<pos ) cout << "_";
	else if( i==pos ){ 
	  
	  if( step%4==0 ){
	    cout << "(>^o^)>";
	  }
	  if( step%4==1 ){
	    cout << "<(^o^)>";
	  }
	  if( step%4==2 ){
	    cout << "<(^o^<)";
	  }
	  if( step%4==3 ){
	    cout << "<( ; )>";
	  }
	  
	}
	else cout << " ";
      }
      progress = (Double_t)( ( nevent+1. )/nentries );
      
      cout << "]" << Int_t( progress*100 ) << "%\r";
      cout.flush();
      if( nevent%10000==0 ) step++;
      
      
      
      
    }
    
    
    
  }
  
  TF1 *fg = new TF1("fg",gausFit_alt,0.0,1.15,5);
  fg->SetRange(0.0,1.15);
  fg->SetParameter(0,1.);
  fg->SetParameter(1,0.9);
  fg->FixParameter(2,0.5);
  fg->SetParameter(3,1);
  fg->SetParameter(4,2.5);
  
  TF1 *fe = new TF1("fe",expoFit,0.0,1.15,2);
  //fe->SetRange(0.0,1.15);
  fe->SetParameter(0,7);
  fe->SetParameter(1,1);
  //fe->FixParameter(0,-0.613473);
  //fe->FixParameter(1,5.27608);
  //fe->FixParameter(0,-0.585486);
  //fe->FixParameter(0,7);
  //fe->FixParameter(1,1);
  //fe->FixParameter(1,7.15219);

  //Can do sideband analysis with RejectPoInt_t() on TF1

  TF1 *f1;
  TF1 *f2;
  TF1 *f3;
  TF1 *f4;
  //TF1 *f5;
  Double_t sig_x;
  Double_t sig_y;
  Double_t sig_r;
  Double_t sig_x_dyCut;
  Double_t BB_elas;
  
  hdx_HCAL->Fit("gaus","Q");
  f1=hdx_HCAL->GetFunction("gaus");
  sig_x = f1->GetParameter(2);
  hdx_HCAL->SetTitle(Form("Fit Sigma = %f",sig_x));

  cout << "Sigma x = " << sig_x << endl;

  hdy_HCAL->Fit("gaus","Q");
  f2=hdy_HCAL->GetFunction("gaus");
  sig_y = f2->GetParameter(2);
  hdy_HCAL->SetTitle(Form("Fit Sigma = %f",sig_y));

  cout << "Sigma y = " << sig_y << endl;

  hdr_HCAL->Fit("gaus","Q");
  f3=hdr_HCAL->GetFunction("gaus");
  sig_r = f3->GetParameter(2);

  cout << "Sigma r = " << sig_r << endl;

  hdx_HCAL_dyCut->Fit("gaus","Q");
  f4=hdx_HCAL_dyCut->GetFunction("gaus");
  sig_x_dyCut = f4->GetParameter(2);

  cout << "Sigma x dyCut = " << sig_x_dyCut << endl;

  hW_cut->Fit(fe,"Q");

  cout << "Extracted BG fit params (y=exp(Ax-B)): A:" << fe->GetParameter(0) << ", B:" << fe->GetParameter(1) << endl;
  
  //hW_cut->Fit("expo","Q");
  //f5=hW_cut->GetFunction("expo");
  //f5->FixParameter(1,2.5);
  //hW_cut->Fit(f5,"Q");

  //TF1 *f5 = new TF1( "f5", "Background Fit", exp([0]*x-[1]), 0, 1.15 ); 
  TF1 *f5 = new TF1("f5",expoFit,0.0,1.15,2);
  f5->FixParameter(0,fe->GetParameter(0));
  f5->FixParameter(1,fe->GetParameter(1));

  BB_elas = hW_cut->GetEntries() - (f5->Integral( 0., 1.15 )*1000.); //Include conversion factor, bins/GeV

  cout << "BB elastics from Integral = " << sig_x_dyCut << endl;

  cout << "Total elastic events detected in BB = " << hits_elBB << endl;

  cout << "Or, alternatively with fits = " << BB_elas << endl;
  
  fout->Write();


  //Simple Integral method
  hits_gHCAL = 100.*f4->Integral(-2.5,2.5); //500 bins to 5.0 range ratio

  cout << endl << "Total elastic events detected in HCal = " << hits_gHCAL << endl;

  Double_t efficiency = Double_t(hits_gHCAL)/Double_t(BB_elas);

  cout << endl << "Detection efficiency for run = " << efficiency*100. << " percent." << endl << endl;

  cout << "Analysis complete. File written to eff.root " << endl;


  cout << "hcalXmax:" << hcalXmax << " hcalXmin:" << hcalXmin << " hcalYmax:" << hcalYmax << " hcalYmin:" << hcalYmin << endl;
  cout << "projXmax:" << projXmax << " projXmin:" << projXmin << " projYmax:" << projYmax << " projYmin:" << projYmin << endl;

  st->Stop();

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

}
