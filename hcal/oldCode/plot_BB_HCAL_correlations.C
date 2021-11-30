#include "TChain.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include <iostream>
#include <fstream>

const double Mp = 0.938272;
const double Mn = 0.939565;

void plot_BB_HCAL_correlations( const char *rootfilename, const char *outfilename, double ebeam=3.7278, 
				double bbtheta=36.0, double sbstheta=31.9, double hcaldist=11.0 ){
  //ifstream infile(configfilename);

  TChain *C = new TChain("T");

  C->Add(rootfilename);

  // TString currentline;

  // while( currentline.ReadLine(infile) && !currentline.BeginsWith("endlist") ){
  //   if( !currentline.BeginsWith("#") ){
  //     C->Add(currentline.Data());
  //   }
  // }

  // double ebeam = 3.7278;
  bbtheta *= TMath::DegToRad();
  sbstheta *= TMath::DegToRad();
  //double hcaldist = 11.0; //m
  double hcalheight = 0.45; //m (we are guessing)

  double Wmin_elastic = 0.8;
  double Wmax_elastic = 1.1;
  
  //HCAL variables:
  double xHCAL, yHCAL, EHCAL;
  
  //BigBite track variables:
  double ntrack;
  int MAXNTRACKS=1000;
  //double trackX[MAXNTRACKS], trackY[MAXNTRACKS], trackXp[MAXNTRACKS],
  double px[MAXNTRACKS], py[MAXNTRACKS], pz[MAXNTRACKS], p[MAXNTRACKS];
  
  double vx[MAXNTRACKS], vy[MAXNTRACKS], vz[MAXNTRACKS];

  //BigBite shower/preshower variables:
  double Eps_BB, Esh_BB;
  double xps_BB, yps_BB, xsh_BB, ysh_BB;

  C->SetBranchStatus("*",0);
  C->SetBranchStatus("sbs.hcal.x",1);
  C->SetBranchStatus("sbs.hcal.y",1);
  C->SetBranchStatus("sbs.hcal.e",1);
  C->SetBranchStatus("bb.tr.n",1);
  C->SetBranchStatus("bb.tr.px",1);
  C->SetBranchStatus("bb.tr.py",1);
  C->SetBranchStatus("bb.tr.pz",1);
  C->SetBranchStatus("bb.tr.vx",1);
  C->SetBranchStatus("bb.tr.vy",1);
  C->SetBranchStatus("bb.tr.vz",1);
  C->SetBranchStatus("bb.ps.e",1);
  C->SetBranchStatus("bb.ps.x",1);
  C->SetBranchStatus("bb.ps.y",1);
  C->SetBranchStatus("bb.sh.e",1);
  C->SetBranchStatus("bb.sh.x",1);
  C->SetBranchStatus("bb.sh.y",1);

  C->SetBranchAddress("sbs.hcal.x",&xHCAL);
  C->SetBranchAddress("sbs.hcal.y",&yHCAL);
  C->SetBranchAddress("sbs.hcal.e",&EHCAL);

  C->SetBranchAddress("bb.tr.n",&ntrack);
  C->SetBranchAddress("bb.tr.px",px);
  C->SetBranchAddress("bb.tr.py",py);
  C->SetBranchAddress("bb.tr.pz",pz);
  C->SetBranchAddress("bb.tr.p",p);
  C->SetBranchAddress("bb.tr.vx",vx);
  C->SetBranchAddress("bb.tr.vy",vy);
  C->SetBranchAddress("bb.tr.vz",vz);
  C->SetBranchAddress("bb.ps.e",&Eps_BB);
  C->SetBranchAddress("bb.sh.e",&Esh_BB);
  C->SetBranchAddress("bb.ps.x",&xps_BB);
  C->SetBranchAddress("bb.ps.y",&yps_BB);
  C->SetBranchAddress("bb.sh.x",&xsh_BB);
  C->SetBranchAddress("bb.sh.y",&ysh_BB);

  TFile *fout = new TFile(outfilename,"RECREATE");


  TH1D *hdpel = new TH1D("hdpel",";p/p_{elastic}(#theta)-1;", 250, -1.0, 0.5);
  TH1D *hW = new TH1D("hW",";W (GeV);", 400,0.0,4.0);
  
  TH1D *hdx_HCAL = new TH1D("hdx_HCAL",";x_{HCAL}-x_{expect} (m);", 500, -2.5, 2.5);
  TH1D *hdy_HCAL = new TH1D("hdy_HCAL",";y_{HCAL}-y_{expect} (m);", 500, -1.25, 1.25);
  TH2D *hdxdy_HCAL = new TH2D("hdxdy_HCAL",";y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)", 250, -1.25, 1.25, 250, -2.5, 2.5 );
  TH2D *hxcorr_HCAL = new TH2D("hxcorr_HCAL",";x_{expect} (m);x_{HCAL} (m)", 250, -2.5, 2.5, 250, -2.5, 2.5 );
  TH2D *hycorr_HCAL = new TH2D("hycorr_HCAL",";y_{expect} (m);y_{HCAL} (m)", 250, -1.25, 1.25, 250, -1.25, 1.25);

  TH1D *hvz = new TH1D("hvz","",250,-0.15,0.15);

  TH2D *hdy_HCAL_vs_z = new TH2D("hdy_HCAL_vs_z","",250,-0.15,0.15,250,-1.25,1.25);
  TH2D *hdy_HCAL_vs_ptheta = new TH2D("hdy_HCAL_vs_ptheta","",250,sbstheta-0.3,sbstheta+0.3,250,-1.25,1.25);


  long nevent = 0;

  while( C->GetEntry( nevent++ ) ){
    if( nevent % 1000 == 0 ) cout << nevent << endl;

    if( ntrack > 0 ){
      double etheta = acos( pz[0]/p[0] );
      double ephi = atan2( py[0], px[0] );

      TVector3 vertex(0,0,vz[0]);

      TLorentzVector Pbeam(0,0,ebeam,ebeam);
      TLorentzVector kprime(px[0],py[0],pz[0],p[0]);
      TLorentzVector Ptarg(0,0,0,Mp);

      TLorentzVector q = Pbeam - kprime;

      TLorentzVector PgammaN = Ptarg + q; //(-px, -py, ebeam - pz, Mp + ebeam - p)
      
      double pel = ebeam/(1.+ebeam/Mp*(1.-cos(etheta)));

      hdpel->Fill( p[0]/pel - 1.0 );

      hW->Fill( PgammaN.M() );

      hvz->Fill( vertex.Z() );

      //Now project to HCAL and compare to best HCAL cluster:
      //Assume neutron (straight-line): 
      //Also assume quasi-elastic kinematics:

      if( PgammaN.M() >= Wmin_elastic && PgammaN.M() <= Wmax_elastic && Eps_BB>0.1 && 
	  fabs( vertex.Z() )<=0.1 ){
	//TVector3 pnucleon_expect = PgammaN.Vect();
	//TVector3 pNhat = pnucleon_expect.Unit();

	//pnucleon_expect gives suspect results. Not clear why. Let's try direct calculation:
	double nu = ebeam - p[0];
	double pp = sqrt(pow(nu,2)+2.*Mp*nu);
	double phinucleon = ephi + TMath::Pi(); //assume coplanarity
	double thetanucleon = acos( (ebeam - p[0]*cos(etheta))/pp ); //use elastic constraint on nucleon kinematics

	TVector3 pNhat( sin(thetanucleon)*cos(phinucleon),sin(thetanucleon)*sin(phinucleon),cos(thetanucleon));

	

	TVector3 HCAL_zaxis(-sin(sbstheta),0,cos(sbstheta));
	TVector3 HCAL_xaxis(0,1,0);
	TVector3 HCAL_yaxis = HCAL_zaxis.Cross(HCAL_xaxis).Unit();
	
	TVector3 HCAL_origin = hcaldist * HCAL_zaxis;

	double sintersect = ( HCAL_origin - vertex ).Dot( HCAL_zaxis ) / (pNhat.Dot( HCAL_zaxis ) );

	TVector3 HCAL_intersect = vertex + sintersect * pNhat;

	double yexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_yaxis );
	double xexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_xaxis );

	hdx_HCAL->Fill( xHCAL - xexpect_HCAL );
	hdy_HCAL->Fill( yHCAL - yexpect_HCAL );

	hdxdy_HCAL->Fill( yHCAL - yexpect_HCAL, xHCAL - xexpect_HCAL );

	hxcorr_HCAL->Fill( xexpect_HCAL, xHCAL );
	hycorr_HCAL->Fill( yexpect_HCAL, yHCAL );

	hdy_HCAL_vs_z->Fill( vertex.Z(), yHCAL - yexpect_HCAL );
	hdy_HCAL_vs_ptheta->Fill( thetanucleon, yHCAL - yexpect_HCAL );

      }

    }
  }
  
  fout->Write();


}
