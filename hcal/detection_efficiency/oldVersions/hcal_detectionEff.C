//SSeeds 11.03.21 - Production - Code designed with help from A.Puckett to project scattered proton during zero-field running to face of HCal and determine detection efficiency. Currently measures relative to BB elastic track detection and sets a limit on good hits in HCal at 2*sigma the proton peak at zero field.

#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
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

TChain *T = 0;

const int kNcell = 288; // Total number of HCal modules
const int kNrows = 24; // Total number of HCal rows
const int kNcols = 12; // Total number of HCal columns
const double kdBlock = 0.152; // Width and height of each module including distance between

const double PI = TMath::Pi();
const double Mp = 0.938272;
const double Mn = 0.939565;

void hcal_detectionEff( const char *configfilename, int run = -1 ){

  //const char *outfilename, double ebeam=3.7278, double bbtheta=36.0, double sbstheta=31.9, double hcaldist=11.0 ){

  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  //TChain *C = new TChain("T");
  //C->Add(rootfilename);

  if( !T ) {
    T = new TChain( "T" );
    T->SetBranchStatus( "*", 0 );
    T->SetBranchStatus( "sbs.hcal.x", 1 );
    T->SetBranchStatus( "sbs.hcal.y", 1 );
    T->SetBranchStatus( "sbs.hcal.e", 1 );
    T->SetBranchStatus( "sbs.hcal.nclus", 1 );
    T->SetBranchStatus( "sbs.hcal.clus.e", 1 );
    T->SetBranchStatus( "sbs.hcal.clus.nblk", 1 );
    T->SetBranchStatus( "bb.tr.n", 1 );
    T->SetBranchStatus( "bb.tr.p", 1 );
    T->SetBranchStatus( "bb.tr.px", 1 );
    T->SetBranchStatus( "bb.tr.py", 1 );
    T->SetBranchStatus( "bb.tr.pz", 1 );
    T->SetBranchStatus( "bb.tr.vx", 1 );
    T->SetBranchStatus( "bb.tr.vy", 1 );
    T->SetBranchStatus( "bb.tr.vz", 1 );
    T->SetBranchStatus( "bb.ps.e", 1 );
    T->SetBranchStatus( "bb.ps.x", 1 );
    T->SetBranchStatus( "bb.ps.y", 1 );
    T->SetBranchStatus( "bb.sh.e", 1 );
    T->SetBranchStatus( "bb.sh.x", 1 );
    T->SetBranchStatus( "bb.sh.y", 1 );
    T->SetBranchStatus( "bb.tdctrig.tdc", 1 );
    T->SetBranchStatus( "bb.tdctrig.tdcelemID", 1 );
    T->SetBranchStatus( "Ndata.bb.tdctrig.tdcelemID", 1 );
    T->SetBranchStatus( "bb.sh.nclus", 1 );

    T->SetBranchAddress( "bb.tr.n", &hcalt::BBtr_n );
    T->SetBranchAddress( "bb.tr.px", hcalt::BBtr_px );
    T->SetBranchAddress( "bb.tr.py", hcalt::BBtr_py );
    T->SetBranchAddress( "bb.tr.pz", hcalt::BBtr_pz );
    T->SetBranchAddress( "bb.tr.vx", hcalt::BBtr_vx );
    T->SetBranchAddress( "bb.tr.vy", hcalt::BBtr_vy );
    T->SetBranchAddress( "bb.tr.vz", hcalt::BBtr_vz );
    T->SetBranchAddress( "bb.tr.p", hcalt::BBtr_p );
    T->SetBranchAddress( "bb.sh.x", hcalt::BBsh_x );
    T->SetBranchAddress( "bb.sh.y", hcalt::BBsh_y );
    T->SetBranchAddress( "bb.sh.e", hcalt::BBsh_e );
    T->SetBranchAddress( "bb.ps.x", hcalt::BBps_x );
    T->SetBranchAddress( "bb.ps.y", hcalt::BBps_y );
    T->SetBranchAddress( "bb.ps.e", hcalt::BBps_e );
    T->SetBranchAddress( "sbs.hcal.x", hcalt::cx );
    T->SetBranchAddress( "sbs.hcal.y", hcalt::cy );
    T->SetBranchAddress( "sbs.hcal.e", hcalt::ce );
    T->SetBranchAddress( "sbs.hcal.nclus", hcalt::nclus );
    T->SetBranchAddress( "sbs.hcal.clus.e", hcalt::ce );
    T->SetBranchAddress( "sbs.hcal.clus.nblk", hcalt::cnblk );
    T->SetBranchAddress( "bb.tdctrig.tdcelemID", hcalt::TDCT_id );
    T->SetBranchAddress( "bb.tdctrig.tdc", hcalt::TDCT_tdc );
    T->SetBranchAddress( "Ndata.bb.tdctrig.tdcelemID", &hcalt::TDCTndata );
  }
  
  double E_e = 1.92; // Energy of the beam
  double emin = 0.5; // Minimum energy of a block to be considered in a cluster
  double HCal_d = 14.5; // Distance to HCal from scattering chamber for comm1
  double HCal_th = 35.0; // Angle that the center of HCal is at  
  double opticsCorr = 1.05; // Correction to magnitude of p_e to account for misaligned optics
  double W_mean = 0.93; // Mean of W at current kinematic
  double W_sig = 0.039; // Width of W at current kinematic
  double BB_th = 36.0; // Angle that the center of BB is at
  double efficiency_rel; // Relative detection efficiency of HCAL (elastic events detected by HCal) / (elastic events as defined by BB tracking)
  int hits_elBB = 0; // Count of total elastic hits detected in BB
  int hits_gHCAL = 0; // Count of all elastic events that HCal detected

  // Reading config file
  cout << "Opening the following files.." << endl;
  ifstream configfile(configfilename);
  TString currentline;
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      T->Add(currentline);
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
    int ntokens = tokens->GetEntries();
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
      }
      if( skey == "BB_th" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        BB_th = sval.Atof();
	cout << "Big Bite theta: " << BB_th << endl;
      }
      if( skey == "emin" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        emin = sval.Atof();
	cout << "Cluster energy minimum " << emin << endl;

      }
    }
    delete tokens;
  }
  cout << endl;

  int nentries = T->GetEntries();
  cout << endl;

  cout << "Opened tree with " << nentries << " entries." << endl;

  BB_th *= TMath::DegToRad();
  HCal_th *= TMath::DegToRad();
  double hcalheight = 0.45; //m (we are guessing)
  
  //Keep a wide cut on W for allowed elastics
  double Wmin_elastic = W_mean - 3.*W_sig;
  double Wmax_elastic = W_mean + 3.*W_sig;
  
  //BigBite tracks
  double ntrack;

  //Outfile
  TFile *fout = new TFile("outfiles/efficiency.root","RECREATE");

  //Histograms
  TH1D *hdpel = new TH1D("hdpel",";p/p_{elastic}(#theta)-1;", 250, -1.0, 0.5);
  TH1D *hW = new TH1D("hW",";W (GeV);", 400,0.0,4.0);
  
  TH1D *hdx_HCAL = new TH1D("hdx_HCAL",";x_{HCAL}-x_{expect} (m);", 500, -2.5, 2.5);
  TH1D *hdy_HCAL = new TH1D("hdy_HCAL",";y_{HCAL}-y_{expect} (m);", 500, -1.25, 1.25);
  TH1D *hdr_HCAL = new TH1D("hdr_HCAL",";r_{HCAL}-r_{expect} (m);", 500, -2.5, 2.5);
  TH2D *hdxdy_HCAL = new TH2D("hdxdy_HCAL",";y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)", 250, -1.25, 1.25, 250, -2.5, 2.5 );
  TH2D *hxcorr_HCAL = new TH2D("hxcorr_HCAL",";x_{expect} (m);x_{HCAL} (m)", 250, -2.5, 2.5, 250, -2.5, 2.5 );
  TH2D *hycorr_HCAL = new TH2D("hycorr_HCAL",";y_{expect} (m);y_{HCAL} (m)", 250, -1.25, 1.25, 250, -1.25, 1.25);

  TH1D *hvz = new TH1D("hvz","",250,-0.15,0.15);

  TH2D *hdy_HCAL_vs_z = new TH2D("hdy_HCAL_vs_z","",250,-0.15,0.15,250,-1.25,1.25);
  TH2D *hdy_HCAL_vs_ptheta = new TH2D("hdy_HCAL_vs_ptheta","",250,HCal_th-0.3,HCal_th+0.3,250,-1.25,1.25);
  TH2D *hdxdy_HCAL_cut = new TH2D("hdxdy_HCAL_cut",";y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)", 250, -1.25, 1.25, 250, -2.5, 2.5 );


  long nevent = 0;

  cout << "Processing events.." << endl;

  while( T->GetEntry( nevent++ ) ){
    if( nevent % 10000 == 0 ) cout << nevent << " / " << nentries << endl;

    //Only consider event if tracks exist
    ntrack = hcalt::BBtr_n;
    if( ntrack > 0 ){
      
      //Define vectors from track in BB
      double etheta = acos( hcalt::BBtr_pz[0]/hcalt::BBtr_p[0] );
      double ephi = atan2( hcalt::BBtr_py[0], hcalt::BBtr_px[0] );

      //Reaction origin
      TVector3 vertex(0,0,hcalt::BBtr_vz[0]);

      TLorentzVector Pbeam(0,0,E_e,E_e);
      TLorentzVector kprime(hcalt::BBtr_px[0],hcalt::BBtr_py[0],hcalt::BBtr_pz[0],hcalt::BBtr_p[0]);
      TLorentzVector Ptarg(0,0,0,Mp);

      TLorentzVector q = Pbeam - kprime;

      TLorentzVector PgammaN = Ptarg + q; //(-px, -py, ebeam - pz, Mp + ebeam - p)
      
      double pel = E_e/(1.+E_e/Mp*(1.-cos(etheta)));

      hdpel->Fill( hcalt::BBtr_p[0]/pel - 1.0 );

      hW->Fill( PgammaN.M() );

      hvz->Fill( vertex.Z() );

      //Now project to HCAL and compare to best HCAL cluster:
      //Assume neutron (straight-line), quasi-elastic kinematics:, and assume BB/HCal trigger coincidence <20ns
      
      //Cut on BBCal and HCal trigger coincidence
      double bbcal_time=0., hcal_time=0.;
      for(int ihit=0; ihit<hcalt::TDCTndata; ihit++){
	if(hcalt::TDCT_id[ihit]==5) bbcal_time=hcalt::TDCT_tdc[ihit];
	if(hcalt::TDCT_id[ihit]==0) hcal_time=hcalt::TDCT_tdc[ihit];
      }
      double diff = hcal_time - bbcal_time; 
      
      //Count the number of elastic events with acceptable tracks in BB (reduction of pions with PS)
      if( PgammaN.M() >= Wmin_elastic && PgammaN.M() <= Wmax_elastic && hcalt::BBps_e[0]>0.1 && fabs( vertex.Z() )<=0.1 ) 
	hits_elBB++;

      //Require that track corresponds to elastic event, energy deposited in preshower, and coincidence trigger between HCal and BB
      if( PgammaN.M() >= Wmin_elastic && PgammaN.M() <= Wmax_elastic && hcalt::BBps_e[0]>0.1 && fabs( vertex.Z() )<=0.1 && fabs(diff-510.)<20 ){

	//Now determine projected expectation for proton from BB e' track
	double nu = E_e - hcalt::BBtr_p[0];
	double pp = sqrt(pow(nu,2)+2.*Mp*nu);
	double phinucleon = ephi + TMath::Pi(); //Coplanarity
	double thetanucleon = acos( (E_e - hcalt::BBtr_p[0]*cos(etheta))/pp ); //Elastic constraint on nucleon kinematics

	//Nucleon momentum unit vector - no field
	TVector3 pNhat( sin(thetanucleon)*cos(phinucleon),sin(thetanucleon)*sin(phinucleon),cos(thetanucleon));

	//Get the coordinate system for HCal from the kinematic HCal theta
	TVector3 HCAL_zaxis(-sin(HCal_th),0,cos(HCal_th));
	TVector3 HCAL_xaxis(0,1,0);
	TVector3 HCAL_yaxis = HCAL_zaxis.Cross(HCAL_xaxis).Unit();
	
	//Project the origin of the detector coordinate system from the scattering chamber
	TVector3 HCAL_origin = HCal_d * HCAL_zaxis;

	//
	double sintersect = ( HCAL_origin - vertex ).Dot( HCAL_zaxis ) / ( pNhat.Dot( HCAL_zaxis ) );
	TVector3 HCAL_intersect = vertex + sintersect * pNhat;


	double yexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_yaxis );
	double xexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_xaxis );
	
	double dx = hcalt::cx[0] - xexpect_HCAL;
	double dy = hcalt::cy[0] - yexpect_HCAL;

	double dr = sqrt( pow( dx, 2 )+pow( dy, 2 ));

	//cout << "dx=" << dx << " dy=" << dy << " dr=" << dr << endl;
	
	hdx_HCAL->Fill( dx );
	hdy_HCAL->Fill( dy );
	hdr_HCAL->Fill( dr );

	hdxdy_HCAL->Fill( dy, dx );

	hxcorr_HCAL->Fill( xexpect_HCAL, hcalt::cx[0] );
	hycorr_HCAL->Fill( yexpect_HCAL, hcalt::cy[0] );

	hdy_HCAL_vs_z->Fill( vertex.Z(), hcalt::cy[0] - yexpect_HCAL );
	hdy_HCAL_vs_ptheta->Fill( thetanucleon, hcalt::cy[0] - yexpect_HCAL );

      }

    }
  }
  
  TF1 *f1;
  TF1 *f2;
  TF1 *f3;
  double sig_x;
  double sig_y;
  double sig_r;
  double r_max;

  hdx_HCAL->Fit("gaus","Q");
  f1=hdx_HCAL->GetFunction("gaus");
  sig_x = f1->GetParameter(2);

  cout << "Sigma x = " << sig_x << endl;

  hdy_HCAL->Fit("gaus","Q");
  f2=hdy_HCAL->GetFunction("gaus");
  sig_y = f2->GetParameter(2);

  cout << "Sigma y = " << sig_y << endl;

  hdr_HCAL->Fit("gaus","Q");
  f3=hdr_HCAL->GetFunction("gaus");
  sig_r = f3->GetParameter(2);

  cout << "Sigma r = " << sig_r << endl;

  double expect_factor = 3.0; //As of now, arbitrary

  //r_max = sqrt( pow( expect_factor*sig_x, 2 )+pow( expect_factor*sig_y, 2 ) );
  r_max = expect_factor*sig_r;

  cout << "Max distance from expected location to be included as good hit = " << r_max << endl;
  cout << "Total elastic events detected in BB = " << hits_elBB << endl;

  fout->Write();

  //Reset nevent counter
  nevent = 0;

  cout << "Reprocessing events to obtain efficiency.." << endl;

  while( T->GetEntry( nevent++ ) ){
    if( nevent % 10000 == 0 ) cout << nevent << " / " << nentries << endl;

    //Only consider event if tracks exist
    ntrack = hcalt::BBtr_n;
    if( ntrack > 0 ){
      
      //Define vectors from track in BB
      double etheta = acos( hcalt::BBtr_pz[0]/hcalt::BBtr_p[0] );
      double ephi = atan2( hcalt::BBtr_py[0], hcalt::BBtr_px[0] );

      TVector3 vertex(0,0,hcalt::BBtr_vz[0]);

      TLorentzVector Pbeam(0,0,E_e,E_e);
      TLorentzVector kprime(hcalt::BBtr_px[0],hcalt::BBtr_py[0],hcalt::BBtr_pz[0],hcalt::BBtr_p[0]);
      TLorentzVector Ptarg(0,0,0,Mp);

      TLorentzVector q = Pbeam - kprime;

      TLorentzVector PgammaN = Ptarg + q; //(-px, -py, ebeam - pz, Mp + ebeam - p)
      
      double pel = E_e/(1.+E_e/Mp*(1.-cos(etheta)));
      
      double bbcal_time=0., hcal_time=0.;
      for(int ihit=0; ihit<hcalt::TDCTndata; ihit++){
	if(hcalt::TDCT_id[ihit]==5) bbcal_time=hcalt::TDCT_tdc[ihit];
	if(hcalt::TDCT_id[ihit]==0) hcal_time=hcalt::TDCT_tdc[ihit];
      }
      double diff = hcal_time - bbcal_time; 

      // Need to keep all cuts to eliminate incidentals from numerator
      if( PgammaN.M() >= Wmin_elastic && PgammaN.M() <= Wmax_elastic && hcalt::BBps_e[0]>0.1 && fabs( vertex.Z() )<=0.1 && fabs(diff-510.)<20 ){

	double nu = E_e - hcalt::BBtr_p[0];
	double pp = sqrt(pow(nu,2)+2.*Mp*nu);
	double phinucleon = ephi + TMath::Pi();
	double thetanucleon = acos( (E_e - hcalt::BBtr_p[0]*cos(etheta))/pp );

	TVector3 pNhat( sin(thetanucleon)*cos(phinucleon),sin(thetanucleon)*sin(phinucleon),cos(thetanucleon));

	TVector3 HCAL_zaxis(-sin(HCal_th),0,cos(HCal_th));
	TVector3 HCAL_xaxis(0,1,0);
	TVector3 HCAL_yaxis = HCAL_zaxis.Cross(HCAL_xaxis).Unit();
	
	TVector3 HCAL_origin = HCal_d * HCAL_zaxis;

	double sintersect = ( HCAL_origin - vertex ).Dot( HCAL_zaxis ) / (pNhat.Dot( HCAL_zaxis ) );

	TVector3 HCAL_intersect = vertex + sintersect * pNhat;

	double yexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_yaxis );
	double xexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_xaxis );
	
	double dx = hcalt::cx[0] - xexpect_HCAL;
	double dy = hcalt::cy[0] - yexpect_HCAL;

	double dr = sqrt( pow( dx, 2 )+pow( dy, 2 ));

	//cout << "dx=" << dx << " dy=" << dy << " dr=" << dr << endl;

	if( dr<r_max ){
	  hdxdy_HCAL_cut->Fill( hcalt::cy[0] - yexpect_HCAL, hcalt::cx[0] - xexpect_HCAL );
	  hits_gHCAL++;
	}
      }
    }
  }

  cout << endl << "Total elastic events detected in HCal = " << hits_gHCAL << endl;

  double efficiency = double(hits_gHCAL)/double(hits_elBB);

  cout << endl << "Detection efficiency for run = " << efficiency*100. << " percent." << endl << endl;

  cout << "Analysis complete. File written to outfiles/efficiency.root " << endl;

  st->Stop();

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

}
