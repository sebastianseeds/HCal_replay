//SSeeds 10.6.22 Attempt to reduce the check to simplest possible algorithm and test on SBS4 data

#include <ctime>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <TH2F.h>
#include <TChain.h>
#include <TCanvas.h>
#include <iostream>
#include <TSystem.h>
#include <TStopwatch.h>
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TString.h"
#include "TCut.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TVector3.h"
#include "TGraphErrors.h"

const Int_t kNcell = 288; // Total number of HCal modules
const Int_t kNrows = 24; // Total number of HCal rows
const Int_t kNcols = 12; // Total number of HCal columns
const Int_t kNtrack = 100; // Reasonable number of tracks in BigBite
const Double_t Xi = -2.20; // Distance from beam center to top of HCal in m
const Double_t Xf = 1.47; // Distance from beam center to bottom of HCal in m
const Double_t Yi = -0.853; // Distance from beam center to opposite-beam side of HCal in m
const Double_t Yf = 0.853; // Distance from beam center to beam side of HCal in m
const Double_t M_p = 0.938272;
const Double_t PI = TMath::Pi();

void hcal_elasDispTest( const char *configfilename="shcal_elasDispTest.cfg", int run = -1 ){

  // Start the chain for root files passed with config file
  TChain *C = new TChain("T");
  
  Double_t E_e = 1.92; // Energy of beam (incoming electrons from accelerator)
  Int_t SBS_field = 0; // Strength in percent of 2100A of SBS magnet
  Double_t W_mean = M_p; // With perfect optics, this should be true. Will read in until then by fitting W distribution on each run
  Double_t W_sigma = 0.03; // Reasonable value by default, but will read in from W distribution

  // Reading config file
  ifstream configfile(configfilename);
  TString currentline;
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      C->Add(currentline);
      cout << "Loaded file: " << currentline <<  endl;
    }    
  }
  TCut globalcut = "";
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endcut") ){
    if( !currentline.BeginsWith("#") ){
      globalcut += currentline;
      cout << "Applying the following cut to all data: " << globalcut << endl;
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
      }
      if( skey == "SBS_field" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	SBS_field = sval.Atoi();
      }
      if( skey == "W_mean" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	W_mean = sval.Atof();
      }
      if( skey == "W_sigma" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	W_sigma = sval.Atof();
      }
      delete tokens;
    }
  }

  //Declare vars
  Double_t atime[kNcell], row[kNcell], col[kNcell], tdctime[kNcell];
  Double_t nblk, nclus, SHnclus, PSnclus, hcalx, hcaly, hcale;
  Double_t Tn;
  Double_t Tpx[kNtrack], Tpy[kNtrack], Tpz[kNtrack], Tp[kNtrack];

  //Setup leaves
  C->SetBranchStatus( "*", 0 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.atime", 1 );  
  C->SetBranchStatus( "sbs.hcal.clus_blk.row", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.col", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.tdctime", 1 );
  C->SetBranchStatus( "sbs.hcal.x", 1 );
  C->SetBranchStatus( "sbs.hcal.y", 1 );
  C->SetBranchStatus( "sbs.hcal.e", 1 );
  C->SetBranchStatus( "sbs.hcal.nblk", 1 );
  C->SetBranchStatus( "sbs.hcal.nclus", 1 );
  C->SetBranchStatus( "bb.sh.nclus", 1 );
  C->SetBranchStatus( "bb.ps.nclus", 1 );
  C->SetBranchStatus( "bb.tr.n", 1 );
  C->SetBranchStatus( "bb.tr.px", 1 );
  C->SetBranchStatus( "bb.tr.py", 1 );
  C->SetBranchStatus( "bb.tr.pz", 1 );
  C->SetBranchStatus( "bb.tr.p", 1 );

  C->SetBranchAddress( "sbs.hcal.clus_blk.atime", atime );
  C->SetBranchAddress( "sbs.hcal.clus_blk.row", row );
  C->SetBranchAddress( "sbs.hcal.clus_blk.col", col);
  C->SetBranchAddress( "sbs.hcal.clus_blk.tdctime", tdctime );
  C->SetBranchAddress( "sbs.hcal.x", &hcalx );
  C->SetBranchAddress( "sbs.hcal.y", &hcaly );
  C->SetBranchAddress( "sbs.hcal.e", &hcale );
  C->SetBranchAddress( "sbs.hcal.nblk", &nblk );
  C->SetBranchAddress( "sbs.hcal.nclus", &nclus );
  C->SetBranchAddress( "bb.sh.nclus", &SHnclus );
  C->SetBranchAddress( "bb.ps.nclus", &PSnclus );
  C->SetBranchAddress( "bb.tr.n", &Tn );
  C->SetBranchAddress( "bb.tr.px", Tpx );
  C->SetBranchAddress( "bb.tr.py", Tpy );
  C->SetBranchAddress( "bb.tr.pz", Tpz );
  C->SetBranchAddress( "bb.tr.p", Tp );

  // Declare outfile
  TFile *fout = new TFile( Form("pSpot_mag%d.root", SBS_field), "RECREATE" );

  TEventList *elist = new TEventList("elist","Elastic Event List");
  C->Draw(">>elist",globalcut);

  // Initialize histogram
  TH1D *h_atime = new TH1D( "atime", "HCal ADC Time, All Channels; ns", 160, 0, 160 );
  TH1D *h_W = new TH1D( "W", "W; GeV", 250, 0.7, E_e*1.5 );
  TH2D *h_phXcorr = new TH2D( "phXcorr", "Azymuthal Correlation BB/HCal; phi_{BB} (deg); X_{HCal}", 100, 0, 0, 100, 0, 0);
  TH2D *h_thYcorr = new TH2D( "thYcorr", "Polar Correlation BB/HCal; theta_{BB} (deg); Y_{HCal}", 100, 0, 0, 100, 0, 0);
  TH2D *elas_rowcol = new TH2D( "elas_rowcol", "HCal Block Position Elastics, HCal; Col, Row", kNcols, -(kNcols+0.5), -0.5, kNrows, -(kNrows+0.5), -0.5 );
  TH2D *elas_pos = new TH2D( "elas_pos", "HCal Position Elastics, HCal; Y(m), X(m)", 100, -Yf, -Yi, 200, -Xf, -Xi );

  Long64_t Nevents = elist->GetN();
  Int_t yield = 0;
  bool minfill = false;
  
  cout << endl << "Proceeding to loop over all events in chain.." << endl;

  for(Long64_t nevent = 0; nevent<Nevents; nevent++){

    if( nevent%10 == 0 ) cout << "Loop: " << nevent << "/" << Nevents << ". Yield: " << yield << ". \r";
    cout.flush();

    C->GetEntry( elist->GetEntry( nevent ) ); 

    if( minfill==false ){
      for( Int_t r=0; r<kNrows; r++){
	for ( Int_t c=0; c<kNcols; c++){
	  elas_rowcol->Fill( -(c+1), -(r+1) );
	}
      }
      minfill=true;
    }

    double E_ep = Tp[0]; // Obtain the scattered electron energy, neglect mass e
    double p_ep = Tp[0]; // Obtain the magnitude of scattered electron momentum
    double Q2 = 2*E_e*E_ep*( 1-(Tpz[0]/p_ep) );
    double nu = E_e-E_ep; // Obtain energy transfer
    double W2 = pow( M_p,2 )+2*M_p*nu-Q2; // Obtain W2 from Q2 and nu
    double W = 0;

    if( W2>0 ){
      W = sqrt( W2 ); // Obtain W for real events by throwing out negative W2 solns
      h_W->Fill( W );
    }

    // Liberal cut on W
    if( fabs( W-W_mean ) > W_sigma ) continue;
    yield++;

    Double_t eth = acos( Tpz[0]/Tp[0] );
    Double_t eph = atan2( Tpy[0], Tpx[0] );
      
    h_phXcorr->Fill( hcalx, eph );
    h_thYcorr->Fill( hcaly, eth );

    h_atime->Fill( atime[0] );

    for( Int_t b=0; b<nblk; b++){
      
      if( atime[b]>40. && atime[b]<65. ){
	elas_rowcol->Fill( -((Double_t)col[b]+1.), -((Double_t)row[b]+1.));
      }
    }

    if( atime[0]>40. && atime[0]<65. && nblk>2 ) elas_pos->Fill( -hcaly, -hcalx );
      
  }
    
  fout->Write();
  fout->Close();
  
  cout << "Histograms populated and written to file: protonDisplay_short.root" << endl;
  
}

