//SSeeds 10.7.22 Check on energy of each participating member of the primary cluster over cuts on trigger sums

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
#include "TLorentzVector.h"

const Int_t kNcell = 288; // Total number of HCal modules
const Int_t kNrows = 24; // Total number of HCal rows
const Int_t kNcols = 12; // Total number of HCal columns
const Int_t kNtrack = 100; // Reasonable number of tracks in BigBite
const Int_t kNsum = 10; // Total number of HCal 2x2 sector sums
const Int_t kNsect = 18; // Total number of HCal 4x4 block sectors
const Int_t kNsrows = 6; // Total number of sector rows
const Int_t kNscols = 3; // Total number of sector columns
const Double_t Xi = -2.20; // Distance from beam center to top of HCal in m
const Double_t Xf = 1.47; // Distance from beam center to bottom of HCal in m
const Double_t Yi = -0.853; // Distance from beam center to opposite-beam side of HCal in m
const Double_t Yf = 0.853; // Distance from beam center to beam side of HCal in m
const Double_t M_p = 0.938272;
const Double_t PI = TMath::Pi();

void sumE( const char *configfilename="ssumE.cfg", int run = -1 ){

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
  Double_t atime[kNcell], row[kNcell], col[kNcell], id[kNcell], tdctime[kNcell], en[kNcell], Tr_ap[kNcell], TrID[kNcell];
  Double_t nblk, nclus, SHnclus, PSnclus, hcalx, hcaly, hcale;
  Double_t Tn;
  Double_t Tpx[kNtrack], Tpy[kNtrack], Tpz[kNtrack], Tp[kNtrack];
  Int_t Trdata;
  Double_t eW2;

  //Setup leaves
  C->SetBranchStatus( "*", 0 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.atime", 1 );  
  C->SetBranchStatus( "sbs.hcal.clus_blk.row", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.col", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.id", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.tdctime", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.e", 1 ); //check
  C->SetBranchStatus( "sbs.hcal.x", 1 );
  C->SetBranchStatus( "sbs.hcal.y", 1 );
  C->SetBranchStatus( "sbs.hcal.e", 1 );
  C->SetBranchStatus( "sbs.hcal.nblk", 1 ); //check
  C->SetBranchStatus( "sbs.hcal.nclus", 1 );
  C->SetBranchStatus( "sbs.trig.a_p", 1 ); //check
  C->SetBranchStatus( "sbs.trig.adcelemID", 1 ); //check
  C->SetBranchStatus( "Ndata.sbs.trig.adcelemID", 1 ); //check
  C->SetBranchStatus( "bb.sh.nclus", 1 );
  C->SetBranchStatus( "bb.ps.nclus", 1 );
  C->SetBranchStatus( "bb.tr.n", 1 );
  C->SetBranchStatus( "bb.ps.e", 1 );
  C->SetBranchStatus( "bb.sh.e", 1 );
  C->SetBranchStatus( "bb.tr.px", 1 );
  C->SetBranchStatus( "bb.tr.py", 1 );
  C->SetBranchStatus( "bb.tr.pz", 1 );
  C->SetBranchStatus( "bb.tr.vz", 1 );
  C->SetBranchStatus( "bb.tr.p", 1 );
  C->SetBranchStatus( "e.kine.W2", 1);

  C->SetBranchAddress( "sbs.hcal.clus_blk.atime", atime );
  C->SetBranchAddress( "sbs.hcal.clus_blk.row", row );
  C->SetBranchAddress( "sbs.hcal.clus_blk.col", col);
  C->SetBranchAddress( "sbs.hcal.clus_blk.id", id);
  C->SetBranchAddress( "sbs.hcal.clus_blk.tdctime", tdctime );
  C->SetBranchAddress( "sbs.hcal.clus_blk.e", en ); //check
  C->SetBranchAddress( "sbs.hcal.x", &hcalx );
  C->SetBranchAddress( "sbs.hcal.y", &hcaly );
  C->SetBranchAddress( "sbs.hcal.e", &hcale ); 
  C->SetBranchAddress( "sbs.hcal.nblk", &nblk ); //check
  C->SetBranchAddress( "sbs.hcal.nclus", &nclus );
  C->SetBranchAddress( "sbs.trig.a_p", Tr_ap ); //check
  C->SetBranchAddress( "sbs.trig.adcelemID", TrID ); //check
  C->SetBranchAddress( "Ndata.sbs.trig.adcelemID", &Trdata ); //check
  C->SetBranchAddress( "bb.sh.nclus", &SHnclus );
  C->SetBranchAddress( "bb.ps.nclus", &PSnclus );
  C->SetBranchAddress( "bb.tr.n", &Tn );
  C->SetBranchAddress( "bb.tr.px", Tpx );
  C->SetBranchAddress( "bb.tr.py", Tpy );
  C->SetBranchAddress( "bb.tr.pz", Tpz );
  C->SetBranchAddress( "bb.tr.p", Tp );
  C->SetBranchAddress( "e.kine.W2", &eW2 );

  // Declare outfile
  TFile *fout = new TFile( "sumEout.root", "RECREATE" );

  TEventList *elist = new TEventList("elist","Elastic Event List");
  C->Draw(">>elist",globalcut);

  // Initialize histograms
  TH1D *h_W = new TH1D( "h_W", "W", 50, 0, 5 );

  TH1D *h_Etot_s[kNsum];
  TH1D *h_NCEtot_s[kNsum];
  TH1D *h_Eblk_s[kNsum];
  TH1D *h_NCEblk_s[kNsum];
  TH2D *h_eVsId_s[kNsum];
  TH2D *h_NCeVsId_s[kNsum];

  for( int s=0; s<kNsum; s++){
    h_Etot_s[s] = new TH1D( Form("h_Etot_s%d",s), Form("HCal Clus Tot E - Sum %d; GeV",s), 200, 0.001, 0.5 );
    h_NCEtot_s[s] = new TH1D( Form("NCh_Etot_s%d",s), Form("NoCut HCal Clus Tot E - Sum %d; GeV",s), 200, 0.001, 0.5 );
    h_Eblk_s[s] = new TH1D( Form("h_Eblk_s%d",s), Form("HCal Clus Blk E - Sum %d; GeV",s), 200, 0., 0.5 );
    h_NCEblk_s[s] = new TH1D( Form("h_NCEblk_s%d",s), Form("NoCut HCal Clus Blk E - Sum %d; GeV",s), 200, 0., 0.5 );
    h_eVsId_s[s] = new TH2D( Form("h_eVsId_s%d",s), Form("HCal Cluster E vs Block - Sum %d; GeV",s), 288, 0.5, 288.5, 200, 0, 0.5 );
    h_NCeVsId_s[s] = new TH2D( Form("NCh_eVsId_s%d",s), Form("NoCut HCal Cluster E vs Block - Sum %d; GeV",s), 288, 0.5, 288.5, 200, 0, 0.5 );
  }
  /*
  TH1D *h_Etot_s0 = new TH1D( "h_Etot_s0", "HCal Clus Tot E - Sum 0; GeV", 200, 0.001, 0.5 );
  TH1D *h_Etot_s1 = new TH1D( "h_Etot_s1", "HCal Clus Tot E - Sum 1; GeV", 200, 0.001, 0.5 );
  TH1D *h_Etot_s2 = new TH1D( "h_Etot_s2", "HCal Clus Tot E - Sum 2; GeV", 200, 0.001, 0.5 );
  TH1D *h_Etot_s3 = new TH1D( "h_Etot_s3", "HCal Clus Tot E - Sum 3; GeV", 200, 0.001, 0.5 );
  TH1D *h_Etot_s4 = new TH1D( "h_Etot_s4", "HCal Clus Tot E - Sum 4; GeV", 200, 0.001, 0.5 );
  TH1D *h_Etot_s5 = new TH1D( "h_Etot_s5", "HCal Clus Tot E - Sum 5; GeV", 200, 0.001, 0.5 );
  TH1D *h_Etot_s6 = new TH1D( "h_Etot_s6", "HCal Clus Tot E - Sum 6; GeV", 200, 0.001, 0.5 );
  TH1D *h_Etot_s7 = new TH1D( "h_Etot_s7", "HCal Clus Tot E - Sum 7; GeV", 200, 0.001, 0.5 );
  TH1D *h_Etot_s8 = new TH1D( "h_Etot_s8", "HCal Clus Tot E - Sum 8; GeV", 200, 0.001, 0.5 );
  TH1D *h_Etot_s9 = new TH1D( "h_Etot_s9", "HCal Clus Tot E - Sum 9; GeV", 200, 0.001, 0.5 );
  TH1D *h_Etot_s10 = new TH1D( "h_Etot_s10", "HCal Clus Tot E - Sum 10; GeV", 200, 0.001, 0.5 );

  TH1D *h_NCEtot_s0 = new TH1D( "h_NCEtot_s0", "NoCut HCal Clus Tot E - Sum 0; GeV", 200, 0.001, 0.5 );
  TH1D *h_NCEtot_s1 = new TH1D( "h_NCEtot_s1", "NoCut HCal Clus Tot E - Sum 1; GeV", 200, 0.001, 0.5 );
  TH1D *h_NCEtot_s2 = new TH1D( "h_NCEtot_s2", "NoCut HCal Clus Tot E - Sum 2; GeV", 200, 0.001, 0.5 );
  TH1D *h_NCEtot_s3 = new TH1D( "h_NCEtot_s3", "NoCut HCal Clus Tot E - Sum 3; GeV", 200, 0.001, 0.5 );
  TH1D *h_NCEtot_s4 = new TH1D( "h_NCEtot_s4", "NoCut HCal Clus Tot E - Sum 4; GeV", 200, 0.001, 0.5 );
  TH1D *h_NCEtot_s5 = new TH1D( "h_NCEtot_s5", "NoCut HCal Clus Tot E - Sum 5; GeV", 200, 0.001, 0.5 );
  TH1D *h_NCEtot_s6 = new TH1D( "h_NCEtot_s6", "NoCut HCal Clus Tot E - Sum 6; GeV", 200, 0.001, 0.5 );
  TH1D *h_NCEtot_s7 = new TH1D( "h_NCEtot_s7", "NoCut HCal Clus Tot E - Sum 7; GeV", 200, 0.001, 0.5 );
  TH1D *h_NCEtot_s8 = new TH1D( "h_NCEtot_s8", "NoCut HCal Clus Tot E - Sum 8; GeV", 200, 0.001, 0.5 );
  TH1D *h_NCEtot_s9 = new TH1D( "h_NCEtot_s9", "NoCut HCal Clus Tot E - Sum 9; GeV", 200, 0.001, 0.5 );
  TH1D *h_NCEtot_s10 = new TH1D( "h_NCEtot_s10", "NoCut HCal Clus Tot E - Sum 10; GeV", 200, 0.001, 0.5 );

  TH1D *h_Eblk_s0 = new TH1D( "h_Eblk_s0", "HCal Clus Blk E - Sum 0; GeV", 200, 0., 0.5 );
  TH1D *h_Eblk_s1 = new TH1D( "h_Eblk_s1", "HCal Clus Blk E - Sum 1; GeV", 200, 0., 0.5 );
  TH1D *h_Eblk_s2 = new TH1D( "h_Eblk_s2", "HCal Clus Blk E - Sum 2; GeV", 200, 0., 0.5 );
  TH1D *h_Eblk_s3 = new TH1D( "h_Eblk_s3", "HCal Clus Blk E - Sum 3; GeV", 200, 0., 0.5 );
  TH1D *h_Eblk_s4 = new TH1D( "h_Eblk_s4", "HCal Clus Blk E - Sum 4; GeV", 200, 0., 0.5 );
  TH1D *h_Eblk_s5 = new TH1D( "h_Eblk_s5", "HCal Clus Blk E - Sum 5; GeV", 200, 0., 0.5 );
  TH1D *h_Eblk_s6 = new TH1D( "h_Eblk_s6", "HCal Clus Blk E - Sum 6; GeV", 200, 0., 0.5 );
  TH1D *h_Eblk_s7 = new TH1D( "h_Eblk_s7", "HCal Clus Blk E - Sum 7; GeV", 200, 0., 0.5 );
  TH1D *h_Eblk_s8 = new TH1D( "h_Eblk_s8", "HCal Clus Blk E - Sum 8; GeV", 200, 0., 0.5 );
  TH1D *h_Eblk_s9 = new TH1D( "h_Eblk_s9", "HCal Clus Blk E - Sum 9; GeV", 200, 0., 0.5 );
  TH1D *h_Eblk_s10 = new TH1D( "h_Eblk_s10", "HCal Clus Blk E - Sum 10; GeV", 200, 0., 0.5 );

  TH1D *h_NCEblk_s0 = new TH1D( "h_NCEblk_s0", "NoCut HCal Clus Blk E - Sum 0; GeV", 200, 0., 0.5 );
  TH1D *h_NCEblk_s1 = new TH1D( "h_NCEblk_s1", "NoCut HCal Clus Blk E - Sum 1; GeV", 200, 0., 0.5 );
  TH1D *h_NCEblk_s2 = new TH1D( "h_NCEblk_s2", "NoCut HCal Clus Blk E - Sum 2; GeV", 200, 0., 0.5 );
  TH1D *h_NCEblk_s3 = new TH1D( "h_NCEblk_s3", "NoCut HCal Clus Blk E - Sum 3; GeV", 200, 0., 0.5 );
  TH1D *h_NCEblk_s4 = new TH1D( "h_NCEblk_s4", "NoCut HCal Clus Blk E - Sum 4; GeV", 200, 0., 0.5 );
  TH1D *h_NCEblk_s5 = new TH1D( "h_NCEblk_s5", "NoCut HCal Clus Blk E - Sum 5; GeV", 200, 0., 0.5 );
  TH1D *h_NCEblk_s6 = new TH1D( "h_NCEblk_s6", "NoCut HCal Clus Blk E - Sum 6; GeV", 200, 0., 0.5 );
  TH1D *h_NCEblk_s7 = new TH1D( "h_NCEblk_s7", "NoCut HCal Clus Blk E - Sum 7; GeV", 200, 0., 0.5 );
  TH1D *h_NCEblk_s8 = new TH1D( "h_NCEblk_s8", "NoCut HCal Clus Blk E - Sum 8; GeV", 200, 0., 0.5 );
  TH1D *h_NCEblk_s9 = new TH1D( "h_NCEblk_s9", "NoCut HCal Clus Blk E - Sum 9; GeV", 200, 0., 0.5 );
  TH1D *h_NCEblk_s10 = new TH1D( "h_NCEblk_s10", "NoCut HCal Clus Blk E - Sum 10; GeV", 200, 0., 0.5 );

  TH2D *h_eVsId_s0 = new TH2D( "h_eVsId_s0", "HCal Cluster E vs Block - Sum 0; GeV", 288, 0.5, 288.5, 200, 0, 0.5 );
  TH2D *h_eVsId_s1 = new TH2D( "h_eVsId_s1", "HCal Cluster E vs Block - Sum 1; GeV", 288, 0.5, 288.5, 200, 0, 0.5 );
  TH2D *h_eVsId_s2 = new TH2D( "h_eVsId_s2", "HCal Cluster E vs Block - Sum 2; GeV", 288, 0.5, 288.5, 200, 0, 0.5 );
  TH2D *h_eVsId_s3 = new TH2D( "h_eVsId_s3", "HCal Cluster E vs Block - Sum 3; GeV", 288, 0.5, 288.5, 200, 0, 0.5 );
  TH2D *h_eVsId_s4 = new TH2D( "h_eVsId_s4", "HCal Cluster E vs Block - Sum 4; GeV", 288, 0.5, 288.5, 200, 0, 0.5 );
  TH2D *h_eVsId_s5 = new TH2D( "h_eVsId_s5", "HCal Cluster E vs Block - Sum 5; GeV", 288, 0.5, 288.5, 200, 0, 0.5 );
  TH2D *h_eVsId_s6 = new TH2D( "h_eVsId_s6", "HCal Cluster E vs Block - Sum 6; GeV", 288, 0.5, 288.5, 200, 0, 0.5 );
  TH2D *h_eVsId_s7 = new TH2D( "h_eVsId_s7", "HCal Cluster E vs Block - Sum 7; GeV", 288, 0.5, 288.5, 200, 0, 0.5 );
  TH2D *h_eVsId_s8 = new TH2D( "h_eVsId_s8", "HCal Cluster E vs Block - Sum 8; GeV", 288, 0.5, 288.5, 200, 0, 0.5 );
  TH2D *h_eVsId_s9 = new TH2D( "h_eVsId_s9", "HCal Cluster E vs Block - Sum 9; GeV", 288, 0.5, 288.5, 200, 0, 0.5 );
  TH2D *h_eVsId_s10 = new TH2D( "h_eVsId_s10", "HCal Cluster E vs Block - Sum 10; GeV", 288, 0.5, 288.5, 200, 0, 0.5 );

  TH2D *h_NCeVsId_s0 = new TH2D( "h_NCeVsId_s0", "NoCut HCal Cluster E vs Block - Sum 0; GeV", 288, 0.5, 288.5, 200, 0, 0.5 );
  TH2D *h_NCeVsId_s1 = new TH2D( "h_NCeVsId_s1", "NoCut HCal Cluster E vs Block - Sum 1; GeV", 288, 0.5, 288.5, 200, 0, 0.5 );
  TH2D *h_NCeVsId_s2 = new TH2D( "h_NCeVsId_s2", "NoCut HCal Cluster E vs Block - Sum 2; GeV", 288, 0.5, 288.5, 200, 0, 0.5 );
  TH2D *h_NCeVsId_s3 = new TH2D( "h_NCeVsId_s3", "NoCut HCal Cluster E vs Block - Sum 3; GeV", 288, 0.5, 288.5, 200, 0, 0.5 );
  TH2D *h_NCeVsId_s4 = new TH2D( "h_NCeVsId_s4", "NoCut HCal Cluster E vs Block - Sum 4; GeV", 288, 0.5, 288.5, 200, 0, 0.5 );
  TH2D *h_NCeVsId_s5 = new TH2D( "h_NCeVsId_s5", "NoCut HCal Cluster E vs Block - Sum 5; GeV", 288, 0.5, 288.5, 200, 0, 0.5 );
  TH2D *h_NCeVsId_s6 = new TH2D( "h_NCeVsId_s6", "NoCut HCal Cluster E vs Block - Sum 6; GeV", 288, 0.5, 288.5, 200, 0, 0.5 );
  TH2D *h_NCeVsId_s7 = new TH2D( "h_NCeVsId_s7", "NoCut HCal Cluster E vs Block - Sum 7; GeV", 288, 0.5, 288.5, 200, 0, 0.5 );
  TH2D *h_NCeVsId_s8 = new TH2D( "h_NCeVsId_s8", "NoCut HCal Cluster E vs Block - Sum 8; GeV", 288, 0.5, 288.5, 200, 0, 0.5 );
  TH2D *h_NCeVsId_s9 = new TH2D( "h_NCeVsId_s9", "NoCut HCal Cluster E vs Block - Sum 9; GeV", 288, 0.5, 288.5, 200, 0, 0.5 );
  TH2D *h_NCeVsId_s10 = new TH2D( "h_NCeVsId_s10", "NoCut HCal Cluster E vs Block - Sum 10; GeV", 288, 0.5, 288.5, 200, 0, 0.5 );
  */
  Long64_t Nevents = elist->GetN();
  Int_t yield = 0;
  bool minfill = false;
  
  cout << endl << "Proceeding to loop over all events in chain.." << endl;

  for(Long64_t nevent = 0; nevent<Nevents; nevent++){

    //if( nevent%10 == 0 ) cout << "Loop: " << nevent << "/" << Nevents << ". Yield: " << yield << ". \r";
    cout.flush();

    C->GetEntry( elist->GetEntry( nevent ) ); 

    bool elas = false;
    double W = eW2;
    h_W->Fill( W );

    // Liberal cut on W
    if( fabs( W-W_mean ) < W_sigma ){
      elas=true;
      yield++;
    }

    for( Int_t b=0; b<Trdata; b++){
      
      int trigID = (int)TrID[b];

      Double_t totE = 0.;
      Double_t NCtotE = 0.;

      if( Tr_ap[b]<=0. ){
	cout << "no trig" << endl;
	continue;
      }

      //cout << endl << "clusEN=";

      for( Int_t i=0; i<nblk; i++ ){
	
	Double_t clusID = id[i];
	Double_t clusEN = en[i];
	
	//cout << clusEN << "+";

	bool coin = atime[b]>40. && atime[b]<65.;
		
	Int_t sector = -1;
		
	/////////////////////
	///HCal Sector Map///
	//Each sector 4x4bl//
	//bl (0,0) topright//
	//Front side f/beam//
	/////////////////////
	/*  3  |  2  |  1  */
	/* ----2-----1---- */
	/*  6  |  5  |  4  */
	/* ----4-----3---- */
	/*  9  |  8  |  7  */
	/* ----6-----5---- */
	/* 12  | 11  | 10  */
	/* ----8-----7---- */
	/* 15  | 14  | 13  */
	/* ---10-----9---- */
	/* 18  | 17  | 16  */
	/////////////////////
	///Each # between ///
	///  is a trigID  ///
	/////////////////////

	//int sRows = 6;
	//int sCols = 3;
	int sMap[kNsrows][kNscols] = {{3,2,1},{6,5,4},{9,8,7},{12,11,10},{15,14,13},{18,17,16}};

	//Reverse the columns to match the map (blocks start at top right)
	int tcol = (int)col[i]/4;
	if( tcol==0 ){
	  tcol=2;
	}else if( tcol==2 ){
	  tcol=0;
	}

	sector = sMap[(int)(row[i]/4)][tcol];

	/*
	if( col[i]>=0 && col[i]<4 ){
	  if( row[i]>=0 && row[i]<4 ) sector=1;
	  if( row[i]>3 && row[i]<8 ) sector=4;
	  if( row[i]>7 && row[i]<12 ) sector=7;
	  if( row[i]>11 && row[i]<16 ) sector=10;
	  if( row[i]>15 && row[i]<20 ) sector=13;
	  if( row[i]>19 && row[i]<24 ) sector=16;
	}
	if( col[i]>3 && col[i]<8 ){
	  if( row[i]>=0 && row[i]<4 ) sector=2;
	  if( row[i]>3 && row[i]<8 ) sector=5;
	  if( row[i]>7 && row[i]<12 ) sector=8;
	  if( row[i]>11 && row[i]<16 ) sector=11;
	  if( row[i]>15 && row[i]<20 ) sector=14;
	  if( row[i]>19 && row[i]<24 ) sector=17;
	}
	if( col[i]>7 && col[i]<12 ){
	  if( row[i]>=0 && row[i]<4 ) sector=3;
	  if( row[i]>3 && row[i]<8 ) sector=6;
	  if( row[i]>7 && row[i]<12 ) sector=9;
	  if( row[i]>11 && row[i]<16 ) sector=12;
	  if( row[i]>15 && row[i]<20 ) sector=15;
	  if( row[i]>19 && row[i]<24 ) sector=18;
	}
	*/
	if( sector<0 || sector>kNsect ){
	  cerr << "ERROR: Cluster member out of acceptance. Exiting..." << endl;
	  break;
	}

	/*
	if( trigID==0 ) {
	  h_NCeVsId_s[0]->Fill( clusID, clusEN );
	  h_NCEblk_s[0]->Fill( clusEN );
	  NCtotE += clusEN;
	}
	if( trigID==1 ) {
	  h_NCeVsId_s[1]->Fill( clusID, clusEN );
	  if( sector==1 || sector==2 || sector==4 || sector==5 ) { 
	    h_NCEblk_s[1]->Fill( clusEN );
	    NCtotE += clusEN;
	  }
	}
	if( trigID==2 ) {
	  h_NCeVsId_s[2]->Fill( clusID, clusEN );
	  if( sector==2 || sector==3 || sector==5 || sector==6 ) {
	    h_NCEblk_s[2]->Fill( clusEN );
	    NCtotE += clusEN;
	  }
	}
	if( trigID==3 ) {
	  h_NCeVsId_s[3]->Fill( clusID, clusEN );
	  if( sector==4 || sector==5 || sector==7 || sector==8 ) {
	    h_NCEblk_s[3]->Fill( clusEN );
	    NCtotE += clusEN;
	  }
	}
	if( trigID==4 ) {
	  h_NCeVsId_s[4]->Fill( clusID, clusEN );
	  if( sector==5 || sector==6 || sector==8 || sector==9 ) {
	    h_NCEblk_s[4]->Fill( clusEN );
	    NCtotE += clusEN;
	  }
	}
	if( trigID==5 ) {
	  h_NCeVsId_s[5]->Fill( clusID, clusEN );
	  if( sector==7 || sector==8 || sector==10 || sector==11 ) {
	    h_NCEblk_s[5]->Fill( clusEN );
	    NCtotE += clusEN;
	  }
	}
	if( trigID==6 ) {
	  h_NCeVsId_s[6]->Fill( clusID, clusEN );
	  if( sector==8 || sector==9 || sector==11 || sector==12 ) {
	    h_NCEblk_s[6]->Fill( clusEN );
	    NCtotE += clusEN;
	  }
	}
	if( trigID==7 ) {
	  h_NCeVsId_s[7]->Fill( clusID, clusEN );
	  if( sector==10 || sector==11 || sector==13 || sector==14 ) {
	    h_NCEblk_s[7]->Fill( clusEN );
	    NCtotE += clusEN;
	  }
	}
	if( trigID==8 ) {
	  h_NCeVsId_s[8]->Fill( clusID, clusEN );
	  if( sector==11 || sector==12 || sector==14 || sector==15 ) {
	    h_NCEblk_s[8]->Fill( clusEN );
	    NCtotE += clusEN;
	  }
	}
	if( trigID==9 ) {
	  h_NCeVsId_s[9]->Fill( clusID, clusEN );
	  if( sector==13 || sector==14 || sector==16 || sector==17 ) {
	    h_NCEblk_s[9]->Fill( clusEN );
	    NCtotE += clusEN;
	  }
	}
	if( trigID==10 ) {
	  h_NCeVsId_s[10]->Fill( clusID, clusEN );
	  if( sector==14 || sector==15 || sector==17 || sector==18 ) {
	    h_NCEblk_s[10]->Fill( clusEN );
	    NCtotE += clusEN;
	  }
	};
	*/	

	if( trigID==0 ) {
	  h_eVsId_s[0]->Fill( clusID, clusEN );
	  h_Eblk_s[0]->Fill( clusEN );
	  NCtotE += clusEN;
	  if( elas==true && coin==true ){
	    h_eVsId_s[0]->Fill( clusID, clusEN );
	    h_Eblk_s[0]->Fill( clusEN );
	    totE += clusEN;
	  }
	}

	for( int s=1; s<kNsum; s++){
	  if( trigID==s ){ 
	    int cStart = (trigID+1)%2;
	    if( cStart == 0 ) cStart = 2;
	    int rStart = (trigID-1)/2;

	    if( sector==sMap[rStart][cStart] || 
		sector==sMap[rStart][cStart-1] || 
		sector==sMap[rStart+1][cStart] || 
		sector==sMap[rStart+1][cStart-1] ){
	      h_NCEblk_s[1]->Fill( clusEN );
	      NCtotE += clusEN;
	      if( elas==true && coin==true ){
		h_Eblk_s[1]->Fill( clusEN );
		totE += clusEN;
	      }
	    }
	    h_NCeVsId_s[s]->Fill( clusID, clusEN );
	    if( elas==true && coin==true )
	      h_eVsId_s[s]->Fill( clusID, clusEN );
	  } 
	}

	/*
	if( elas==false || coin==false ) continue;
	
	if( trigID==0 ) {
	  h_eVsId_s[0]->Fill( clusID, clusEN );
	  h_Eblk_s[0]->Fill( clusEN );
	  totE += clusEN;
	}

	for( int s=1; s<kNsum; s++){
	  if( trigID==s ){ 
	    h_eVsId_s[s]->Fill( clusID, clusEN );
	    
	    int cStart = (trigID+1)%2;
	    if( cStart == 0 ) cStart = 2;
	    int rStart = (trigID-1)/2;
	    if( sector==sMap[rStart][cStart] || 
		sector==sMap[rStart][cStart-1] || 
		sector==sMap[rStart+1][cStart] || 
		sector==sMap[rStart+1][cStart-1] ){
	      h_Eblk_s[1]->Fill( clusEN );
	      totE += clusEN;
	    }
	  } 
	}
	*/
	/*
	if( trigID==0 ) {
	  h_eVsId_s[0]->Fill( clusID, clusEN );
	  h_Eblk_s[0]->Fill( clusEN );
	  totE += clusEN;
	}
	if( trigID==1 ) {
	  h_eVsId_s[1]->Fill( clusID, clusEN );
	  if( sector==1 || sector==2 || sector==4 || sector==5 ) {
	    h_Eblk_s[1]->Fill( clusEN );
	    totE += clusEN;
	  }
	}
	if( trigID==2 ) {
	  h_eVsId_s[2]->Fill( clusID, clusEN );
	  if( sector==2 || sector==3 || sector==5 || sector==6 ) {
	    h_Eblk_s[2]->Fill( clusEN );
	    totE += clusEN;
	  }
	}
	if( trigID==3 ) {
	  h_eVsId_s[3]->Fill( clusID, clusEN );
	  if( sector==4 || sector==5 || sector==7 || sector==8 ) {
	    h_Eblk_s[3]->Fill( clusEN );
	    totE += clusEN;
	  }
	}
	if( trigID==4 ) {
	  h_eVsId_s[4]->Fill( clusID, clusEN );
	  if( sector==5 || sector==6 || sector==8 || sector==9 ) {
	    h_Eblk_s[4]->Fill( clusEN );
	    totE += clusEN;
	  }
	}
	if( trigID==5 ) {
	  h_eVsId_s[5]->Fill( clusID, clusEN );
	  if( sector==7 || sector==8 || sector==10 || sector==11 ) {
	    h_Eblk_s[5]->Fill( clusEN );
	    totE += clusEN;
	  }
	}
	if( trigID==6 ) {
	  h_eVsId_s[6]->Fill( clusID, clusEN );
	  if( sector==8 || sector==9 || sector==11 || sector==12 ) {
	    h_Eblk_s[6]->Fill( clusEN );
	    totE += clusEN;
	  }
	}
	if( trigID==7 ) {
	  h_eVsId_s[7]->Fill( clusID, clusEN );
	  if( sector==10 || sector==11 || sector==13 || sector==14 ) {
	    h_Eblk_s[7]->Fill( clusEN );
	    totE += clusEN;
	  }
	}
	if( trigID==8 ) {
	  h_eVsId_s[8]->Fill( clusID, clusEN );
	  if( sector==11 || sector==12 || sector==14 || sector==15 ) {
	    h_Eblk_s[8]->Fill( clusEN );
	    totE += clusEN;
	  }
	}
	if( trigID==9 ) {
	  h_eVsId_s[9]->Fill( clusID, clusEN );
	  if( sector==13 || sector==14 || sector==16 || sector==17 ) {
	    h_Eblk_s[9]->Fill( clusEN );
	    totE += clusEN;
	  }
	}
	if( trigID==10 ) {
	  h_eVsId_s[10]->Fill( clusID, clusEN );
	  if( sector==14 || sector==15 || sector==17 || sector==18 ) {
	    h_Eblk_s[10]->Fill( clusEN );
	    totE += clusEN;
	  }
	}
	*/

      }	  
      
      //cout << endl << "For TrigID:" << trigID << ", Etot:" << totE << ", nocut:" << NCtotE << endl;

      for( int s=0; s<kNsum; s++){
	if( trigID==s ){
	  h_Etot_s[s]->Fill(totE);
	  h_NCEtot_s[s]->Fill(NCtotE);
	}
      }

      /*
      if( trigID==0 ) {
	h_Etot_s[0]->Fill(totE);
	h_NCEtot_s[0]->Fill(NCtotE);
      }
      if( trigID==1 ) {
	h_Etot_s[1]->Fill(totE);
	h_NCEtot_s[1]->Fill(NCtotE);
      }
      if( trigID==2 ) {
	h_Etot_s[2]->Fill(totE);
	h_NCEtot_s[2]->Fill(NCtotE);
      }
      if( trigID==3 ) {
	h_Etot_s[3]->Fill(totE);
	h_NCEtot_s[3]->Fill(NCtotE);
      }
      if( trigID==4 ) {
	h_Etot_s[4]->Fill(totE);
	h_NCEtot_s[4]->Fill(NCtotE);
      }
      if( trigID==5 ) {
	h_Etot_s[5]->Fill(totE);
	h_NCEtot_s[5]->Fill(NCtotE);
      }
      if( trigID==6 ) {
	h_Etot_s[6]->Fill(totE);
	h_NCEtot_s[6]->Fill(NCtotE);
      }
      if( trigID==7 ) {
	h_Etot_s[7]->Fill(totE);
	h_NCEtot_s[7]->Fill(NCtotE);
      }
      if( trigID==8 ) {
	h_Etot_s[8]->Fill(totE);
	h_NCEtot_s[8]->Fill(NCtotE);
      }
      if( trigID==9 ) {
	h_Etot_s[9]->Fill(totE);
	h_NCEtot_s[9]->Fill(NCtotE);
      }
      if( trigID==10 ) {
	h_Etot_s[10]->Fill(totE);
	h_NCEtot_s[10]->Fill(NCtotE);
      }
      */
    } 
  }
   
  //Make a few canvases to hold all important plots
  /*
  TCanvas *sum_Etot = new TCanvas("sum_Etot","Energy Total by Sum",1600,1200);
  sum_Etot->Divide(2,5);
  gStyle->SetOptStat(0);
  for( int i=0; i<10; i++ ){
    



  }

  */
  /*
  TCanvas *TDC_top = new TCanvas("TDC_top","TDC_top",1600,1200);
  TCanvas *TDC_bot = new TCanvas("TDC_bot","TDC_bot",1600,1200);
  TCanvas *ADCt_top = new TCanvas("ADCt_top","ADCt_top",1600,1200);
  TCanvas *ADCt_bot = new TCanvas("ADCt_bot","ADCt_bot",1600,1200);

  TDC_top->Divide(12,12);
  TDC_bot->Divide(12,12);
  ADCt_top->Divide(12,12);
  ADCt_bot->Divide(12,12);

  gStyle->SetOptStat(0);
  for(Int_t i=0; i<kNcell; i++){
    TDC_top->cd(i+1);
    if( i>=144 ){
      TDC_bot->cd(i-143);
      gStyle->SetOptStat(0);
    }
    if(htdcDiff_corr[i]->GetEntries()<tFitMin){
      htdcDiff_corr[i]->SetAxisColor(2);
    }else{
      htdcDiff_corr[i]->SetAxisColor(1);
    }
    htdcDiff_corr[i]->Draw();
  }

  gStyle->SetOptStat(0);
  for(Int_t i=0; i<kNcell; i++){
    ADCt_top->cd(i+1);
    if( i>=144 ){
      ADCt_bot->cd(i-143);
      gStyle->SetOptStat(0);
    }
    if(hatime[i]->GetEntries()<tFitMin){
      hatime[i]->SetAxisColor(2);
    }else{
      hatime[i]->SetAxisColor(1);
    }
    hatime[i]->Draw();
  }
  */

  fout->Write();
  fout->Close();
  
  cout << "Histograms populated and written to file: sumEout.root" << endl;
  
}
