//SSeeds 3.24.22 - Post-production - Script written to evaluate performance of cluster timing cut added to SBS-offline <SBSCalorimeter.cxx> using HCal as a basis detector

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

const int kNcell = 288; // Total number of HCal modules
const int kNrows = 24; // Total number of HCal rows
const int kNcols = 12; // Total number of HCal columns

const double PI = TMath::Pi();
const double M_e = 0.00051;
const double M_p = 0.938272;
const double M_n = 0.939565;
const double c_light = 299792458.0;

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

void clus_time_check( int run = -1, double tmax = 20 ){
  
  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  //Declare branch variables
  Double_t cid[kNcell] = {0};
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
  Double_t nclus[kNcell] = {0}; 
  Double_t cblkid[kNcell] = {0};

  if( !T ) {
    T = new TChain( "T" );

    T->Add(Form("/volatile/halla/sbs/seeds/rootfiles/hcal_general_%d*",run));

    T->SetBranchStatus( "*", 0 ); //Turn off all branches for faster execution
    //Select only necessary branches for processing
    T->SetBranchStatus( "sbs.hcal.clus.id", 1 );
    T->SetBranchStatus( "sbs.hcal.clus.row", 1 );
    T->SetBranchStatus( "sbs.hcal.clus.col", 1 );
    T->SetBranchStatus( "sbs.hcal.clus.atime", 1 );
    T->SetBranchStatus( "sbs.hcal.clus_blk.atime", 1 );
    T->SetBranchStatus( "sbs.hcal.clus_blk.row", 1 );
    T->SetBranchStatus( "sbs.hcal.clus_blk.col", 1 );
    T->SetBranchStatus( "sbs.hcal.clus.e", 1 );
    T->SetBranchStatus( "sbs.hcal.clus.eblk", 1 );
    T->SetBranchStatus( "sbs.hcal.clus.nblk", 1 );
    T->SetBranchStatus( "sbs.hcal.nclus", 1 );
    T->SetBranchStatus( "sbs.hcal.clus_blk.id", 1 );
    T->SetBranchStatus( "sbs.hcal.clus_blk.e", 1 );

    //Set branches to variables
    T->SetBranchAddress("sbs.hcal.clus.id",cid);
    T->SetBranchAddress("sbs.hcal.clus.row",crow);
    T->SetBranchAddress("sbs.hcal.clus.col",ccol);
    T->SetBranchAddress("sbs.hcal.clus.atime",catime); //a_time for main block in each cluster
    T->SetBranchAddress("sbs.hcal.clus_blk.atime",cblkatime); //a_time for each block in main cluster
    T->SetBranchAddress("sbs.hcal.clus_blk.row",cbrow);
    T->SetBranchAddress("sbs.hcal.clus_blk.col",cbcol);
    T->SetBranchAddress("sbs.hcal.clus.e",ce);
    T->SetBranchAddress("sbs.hcal.clus.eblk",ceblk); // Highest energy block
    T->SetBranchAddress("sbs.hcal.clus.nblk",cnblk); //number of blocks in cluster per event
    T->SetBranchAddress("sbs.hcal.nclus",nclus); //number of clusters
    T->SetBranchAddress("sbs.hcal.clus_blk.id",cblkid);
    T->SetBranchAddress("sbs.hcal.clus_blk.e",cblke); // Array of block energies

    cout << endl << "Opened up tree from run " << run << " with nentries = " << T->GetEntries() << "." << endl;
  }
 
  cout << endl << "Maximum time between blocks in same cluster from db_sbs.hcal.dat = " << tmax << "." << endl;

  // Define a clock to check macro processing time
  TStopwatch *sw = new TStopwatch();
  sw->Start();

  // Get the date
  string date = getDate();
  
  // Declare outfile
  TFile *fout = new TFile( Form("outfiles/clus_time_check_out_run%d_%s.root", run, date.c_str()), "RECREATE" );

  // Initialize global histograms
  TH1D *hMaxClusT = new TH1D( "maxClusT", Form("maxClusT tmax%f",tmax), 310, -10, 300 );
  TH2D *hMaxClusTvC = new TH2D( "maxClusTvC", Form("maxClusTvC tmax%f",tmax), 288, 0, 288, 310, -10, 300);
  TH1D *hClusE = new TH1D( "ClusE", "ClusE", 100, 0., 1. );
  TH1D *hClusBlkE = new TH1D( "ClusBlkE", "ClusBlkE", 100, 0., 1. );
  TH1D *hClusPrimeBlkE = new TH1D( "ClusPrimeBlkE", "ClusPrimeBlkE", 100, 0., 1. );

  // Initialize cell histograms  
  TH1D *hMaxClusTpC[kNcell];
  TH1D *hClusBlkEpC[kNcell];

  for( int c=0; c<kNcell+1; c++ ){  
    hMaxClusTpC[c] = new TH1D( Form( "ADC Time Max Diff Cell%d", c ), Form( "ADC Time Max Diff Cell%d tmax%f", c, tmax ), 200, -100, 100 );
      hClusBlkEpC[c] = new TH1D( Form( "Cluster Block E Cell%d", c ), Form( "Cluster Block E Cell%d", c ), 100, 0., 1. );
  }
  
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
	timeremains = timekeeper*( double(Nevents)/double(nevent) - 1. ); 
      sw->Reset();
      sw->Continue();
      
      progress = (double)((nevent+1.)/Nevents);
      cout << "] " << int(progress*100.) << ", time remaining: " << int(timeremains/60.) << "m \r";
      cout.flush();
      
      T->GetEntry( nevent ); 
      
      // Fill primary cell histograms
      hClusE->Fill(ce[0]); //First element is primary cell in cluster
      hClusPrimeBlkE->Fill(cblke[0]);
      
      double maxDiff = 0.0;
      int maxDiffID = -1;
      
      for( int b=0; b<cnblk[0]; b++ ){
	
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

	int id = cblkid[b];
	//if(id>200) cout << id << endl;
	double diff = fabs(cblkatime[b]-cblkatime[0]);
	
	if( diff>maxDiff ){
	  maxDiff = diff;
	  maxDiffID = id;
	}

	hClusBlkE->Fill(cblke[b]);
	hMaxClusT->Fill(cblkatime[b]);
	hClusBlkEpC[id]->Fill(cblke[b]);
      }
      
      if( maxDiff!=0.0 ){
	hMaxClusTpC[maxDiffID]->Fill(maxDiff);
	hMaxClusTvC->Fill(maxDiffID,maxDiff);
      }
    }
  }
  
  fout->Write();

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;
  
}

