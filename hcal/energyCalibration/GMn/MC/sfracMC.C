//SSeeds 3.2.23 GMn: Script to evaluate HCal sampling fraction by kinematic from MC data with primary cluster constraint

#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <TChain.h>
#include <TSystem.h>
#include <TStopwatch.h>
#include <unistd.h>
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
#include "TString.h"
#include "TF1.h"
#include "../../include/MC.C"

//#include <bits/stdc++.h>

//Detector constants
const Int_t ifac = 3; // Inclusion factor, number of sigma to keep around cut peaks
const Int_t kNcell = 288; // Total number of HCal modules
const Int_t kNrows = 24; // Total number of HCal rows
const Int_t kNcols = 12; // Total number of HCal columns
const Int_t maxTracks = 1000; // Reasonable limit on tracks to be stored per event
const Int_t uni_N = 400; // Total number of bins used to measure detection uniformity (hSampFrac histos)
const Int_t xN = 48; //2*kNrows, total number of dispersive bins detection uni
const Int_t yN = 24; //2*kNcols, total number of transverse bins detection uni
//Double_t hcalheight = 0.365; //m The height of the center of HCAL above beam
const Double_t hcalheight = -0.2897;
const UInt_t second = 1000000; 
const Int_t nkine = 6; // Total number of kinematic points
const Double_t Rmax = 0.30; //Max x,y,r to be included in cluster from cluster center

//Math
const Double_t PI = TMath::Pi();
const Double_t M_e = 0.00051;
const Double_t M_p = 0.938272;

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

//Keep this around just in case it's needed later
// Double_t Gfit(Double_t *x, Double_t *par){
//   Double_t amp = par[0];
//   Double_t offset = par[1];
//   Double_t sigma = par[2];
//   return amp*exp(-0.5*pow((x[0]-offset)/sigma,2.));
// }

//Main.
void sfracMC( Int_t kine=-1 ){

  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );
  
  // Get the date
  string date = getDate();

  TChain *C = new TChain("T");

  C->Add(Form("/lustre19/expphy/volatile/halla/sbs/seeds/simulation/gmn_sbs%d*",kine));

  // Declare Chain for many root files
  // First obtain index for nfset constants
  Int_t kIdx;
  Double_t beamE; //GeV
  if( kine == 4 ){
    kIdx=0;
    beamE=3.728;
  }
  if( kine == 7 ){
    kIdx=1;
    beamE=7.907;
  } 
  if( kine == 11 ){
    kIdx=2;
    beamE=9.859;
  }
  if( kine == 14 ){
    kIdx=3;
    beamE=5.965;
  }
  if( kine == 8 ){
    kIdx=4;
    beamE=5.965;
  }
  if( kine == 9 ){
    kIdx=5;
    beamE=4.015;
  }

  cout << endl << "Setup parameters loaded and chains linked." << endl;

  // Declare outfile
  TFile *fout = new TFile( Form("output/sfracMC_sbs%d.root", kine ), "RECREATE" );;

  MC *G = new MC(C);

  Int_t Nevents = C->GetEntries();

  //Declare diagnostic histograms (as sparse as possible)
  TH1D *hSampFrac_det = new TH1D( "hSampFrac_det","HCal Det Sum E Dep / (Beam E - P Track E)", 400, 0., 1. );
  TH1D *hSampFrac_clus = new TH1D( "hSampFrac_clus","HCal Cluster E Dep / (Beam E - P Track E)", 400, 0., 1. );
  TH1D *hep = new TH1D( "hep", "Event Electron Momentum", 500, 0., 5. );
  TH2D *hcEvID = new TH2D( "hcEvID", "HCal Cluster Element E vs Cell", 288, 0, 288, 400, 0., 1.);
  TH1D *csize = new TH1D( "csize", "P Cluster size, hit thresh 0.01 GeV", 30, 0, 30 );
  TH1D *bhits = new TH1D( "bhits", "Number of blocks over 1MeV thresh", 30, 0, 30 );

  //Loop over events
  cout << "Main loop over all data commencing.." << endl;

  for( Long64_t nevent = 1; nevent <Nevents; nevent++){
      
    if ( nevent%1000==0 ) cout << "MC LD2 kinematic " << kine << ", entry: " << nevent << "/" << Nevents << " \r";
      cout.flush();

      G->GetEntry( nevent ); 
      ///////////////////////////////////////////////////
      //General cuts
      Int_t ntracks = G->PTrack_ntracks;
      //if( ntracks!=1 ) continue; //Cut on number of tracks
      //cout << ntracks << endl;

      Int_t PID = (*(G->PTrack_PID))[0]; //electron, PID==11
      //if( PID!=11 ) continue; //Cut on electrons only
      //cout << PID << endl;

      Double_t xptrack = (*(G->PTrack_momx))[0]; //primary track momentum
      if( xptrack<0 ) continue;
      //cout << xptrack << endl;

      Int_t nhits = G->Harm_HCalScint_hit_nhits;
      if( nhits<1 ) continue; //Cut on hits in hcal
      //cout << nhits << endl;

      Double_t PSE = G->Earm_BBPSTF1_det_esum;
      if( PSE<0.2 ) continue; //Cut on BBCal preshower energy
      //cout << PSE << endl;

      ///////////////////////////////////////////////////

      Double_t hcalE = G->Harm_HCalScint_det_esum;
      Double_t ep = G->ev_ep;

      Double_t detSF = hcalE/(beamE-ep);
      hSampFrac_det->Fill(detSF);

      Int_t blockhits = 0;
      Double_t clustermax = 0.;
      Int_t highblock = -1;
      Int_t highblockR = -1;
      Int_t highblockC = -1;
      for( Int_t h=0; h<nhits; h++ ){
	Double_t blockE = (*(G->Harm_HCalScint_hit_sumedep))[h];
	if( blockE < 0.01 ) continue; //Min E to be considered for a cluster in SBS-offline
	Int_t cell = (*(G->Harm_HCalScint_hit_cell))[h];
	Int_t row = (*(G->Harm_HCalScint_hit_row))[h];
	Int_t col = (*(G->Harm_HCalScint_hit_col))[h];

	blockhits++;
	if( blockE>clustermax ){
	  clustermax=blockE;
	  highblock=cell;
	  highblockR=row;
	  highblockC=col;
	}
      }
      if( highblock==-1 ) continue; //Cut event if no cluster seed of sufficient E is found

      if( highblockR==0 || 
	  highblockR==23 || 
	  highblockC==0 || 
	  highblockC==11 ) continue; //All events with primary cluster element on edge blocks cut
      
      //Make check on block position to build primary cluster
      Double_t clusterE = 0.;
      Int_t clustersize = 0;
      for( Int_t h=0; h<nhits; h++ ){
	Double_t blockE = (*(G->Harm_HCalScint_hit_sumedep))[h];
	if( blockE < 0.01 ) continue; //Min E to be considered for a cluster in SBS-offline
	Int_t cell = (*(G->Harm_HCalScint_hit_cell))[h];
	Int_t row = (*(G->Harm_HCalScint_hit_row))[h];
	Int_t col = (*(G->Harm_HCalScint_hit_col))[h];

	if( abs(highblockR-row)>3 ) continue; //Cluster row cut
	if( abs(highblockC-col)>3 ) continue; //Cluster col cut
	clustersize++;
	clusterE+=blockE;
	hcEvID->Fill(cell,blockE);
      }

      if( clusterE<0.01 ) continue; //Mirror min cluster energy required to be considered in calibration

      hSampFrac_clus->Fill(clusterE/(beamE-ep));
      hep->Fill(ep);
      csize->Fill(clustersize);
      bhits->Fill(blockhits);

  }

  fout->Write();

  st->Stop();

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

}

//int MID = (*(G->SDTrack_MID))[(*(G->Harm_HCalScint_hit_sdtridx))[ihit]];
