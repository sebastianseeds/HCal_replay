//SSeeds 12.16.21 - Production - Script to search over general replays and plot all BBCal trigger times relative to L1A 
// relevant database file can be found $SBS_REPLAY/DB/db_bb.tdctrig.dat

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
#include "hcal.h"

TChain *T = 0;
/*
const int kNcell = 288; // Total number of HCal modules
const int kNrows = 24; // Total number of HCal rows
const int kNcols = 12; // Total number of HCal columns
const int minSample = 0.0;
const int maxSample = 40.0;
const int totSample = (maxSample-minSample); //Should be the total fADC window size with each samp = 4ns
*/
//void hcal_signal_centering();

//gErrorIgnoreLevel = kFatal;

  //The following runs are from SBS11
const int runStart = 12236;
//const int runEnd = 12246;
const int runEnd = 12935;
//const int runTotal = 10;
//int runTotal = runEnd - runStart;

//TH1D *hTrigBBCal[runTotal];

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
/*
// Create generic histogram function
TH1D* MakeHisto(int row, int col, int bins){
  TH1D *h = new TH1D(Form("h%02d%02d",row,col),Form("%d-%d",row+1,col+1),bins,minSample,maxSample);
  h->SetStats(0);
  h->SetLineWidth(2);
  return h;
}
*/
// Main
void bbTrigShift(){
  
  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );
		  
  // Get the date
  string date = getDate();
  
  // Declare outfile
  TFile *fout = new TFile( Form("outfiles/bbTrigShift_%s.root",date.c_str()), "RECREATE" );
   
  // Declare histograms and arrays
  TH2D *hbbTrigShift = new TH2D( "hbbTrigShift","BBCal Trigger vs Run Number", 1000,12000,13000, 1000, 0, 1000 );
  TH1D *hTrigBBCal = new TH1D( "hTrigBBCal","BBCal Trig (ns)", 2000, -1000, 1000 );
  
  //TH1D *hTrigHCal = new TH1D( "hTrigHCal","1190 HCal Trig (ns)", 2000, -1000, 1000 );

  //TH1D *ATimeDist = new TH1D( "ADC time over all channels", "ADCTime", 4*totSample, 4*minSample, 4*maxSample );
  //TF1 *f1;
  //int gEvt[kNrows][kNcols] = {0};
  int runTotal = runEnd - runStart;

  /*
  for( int i=0; i<runTotal; i++){
    hTrigBBCal[i] = new TH1D( Form("hTrigBBCal_%d",i+runStart),Form("BBCal Trig %d (ns)",i+runStart), 2000, -1000, 1000 );
  }
  */
  double progress = 0.;

  while( progress<1.0 ){
    
    int barwidth = 70;
    int step = 1;
    //Loop over all runs investigated
    for( int i=runStart; i<runEnd; i++ ){
      
      hTrigBBCal->Reset("ICES");

      //if( !T ) {
      T = new TChain( "T" );
      //cout << "Loading root file from run " << i << ".." << endl;

      //if(!Form("/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_%d_stream0_seg0_0.root",i)) cout << "FILE DNE" << endl;

      //if(Form("/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_%d_stream0_seg0_0.root",i)) {

      T->Add( Form("/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_%d_stream0_seg0_0.root",i) );
      T->Add( Form("/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_%d_stream0_seg0_0_firstevent0_nevent10000.root",i) );
      T->Add( Form("/adaqfs/home/a-onl/sbs/Rootfiles/hcal_trim_%d_500000.root",i) );
      T->Add( Form("/adaqfs/home/a-onl/sbs/Rootfiles/hcal_trim_%d_50000.root",i) );
      T->Add( Form("/adaqfs/home/a-onl/sbs/Rootfiles/hcal_trim_%d_10000.root",i) );
      T->Add( Form("/adaqfs/home/a-onl/sbs/Rootfiles/hcal_trim_%d_100000.root",i) );
      //T->Add( Form("/adaqfs/home/a-onl/sbs/Rootfiles/bbshower_gmn_%d_stream0_seg0_0.root",i) );

	//}
      //if(!T) cout << "NO FILE EXISTS 1" << endl;
      if(T->GetEntries()==0) {
	cout << "File from run " << i << " does not exist." << endl;
	//continue;
      }
      //T->Add( Form("/adaqfs/home/a-onl/sbs/Rootfiles/hcal_gmn_%d_4100000_70.root",run) );
      //T->Add( Form("/adaqfs/home/a-onl/sbs/Rootfiles/hcal_gmn_%d_4100000_71.root",run) );
      //T->Add( Form("/adaqfs/home/a-onl/sbs/Rootfiles/hcal_gmn_%d_4100000_72.root",run) );
      //T->Add( Form("/adaqfs/home/a-onl/sbs/Rootfiles/hcal_gmn_%d_4100000_73.root",run) );
      //T->Add( Form("/adaqfs/home/a-onl/sbs/Rootfiles/hcal_gmn_%d_4100000_74.root",run) );

      //T->Add( Form("/adaqfs/home/a-onl/sbs/Rootfiles/hcal_trim_%d*.root",run) );
      T->SetBranchStatus( "*", 0 );
      //T->SetBranchStatus( "sbs.hcal.*", 1 );
      T->SetBranchStatus( "bb.tdctrig.tdc", 1 );
      T->SetBranchStatus( "bb.tdctrig.tdcelemID", 1 );
      T->SetBranchStatus( "Ndata.bb.tdctrig.tdcelemID", 1 );

      //T->SetBranchAddress( "sbs.hcal.tdcrow", hcalt::trow );
      //T->SetBranchAddress( "sbs.hcal.tdccol", hcalt::tcol );    
      //T->SetBranchAddress( "sbs.hcal.adcrow", hcalt::row );
      //T->SetBranchAddress( "sbs.hcal.adccol", hcalt::col );
      //T->SetBranchAddress( "sbs.hcal.a_time", hcalt::a_time );
      T->SetBranchAddress( "bb.tdctrig.tdcelemID", hcalt::TDCT_id );
      T->SetBranchAddress( "bb.tdctrig.tdc", hcalt::TDCT_tdc );
      T->SetBranchAddress( "Ndata.bb.tdctrig.tdcelemID", &hcalt::TDCTndata );
      //T->SetBranchStatus( "Ndata.sbs.hcal.adcrow", 1 );
      //T->SetBranchAddress( "Ndata.sbs.hcal.adcrow", &hcalt::ndata );
      //}
      
      Long64_t Nevents = T->GetEntries();
      cout << "Total events in tree: " << Nevents << ".." << endl;
      
      // Need to get electron energy and W for cuts on elastic events
      long mevent = 0;
      
      //cout << "Building timing and trigger histograms.." << endl;
      
      //double progress = 0.;
      /*
	while( progress<1.0 ){
	
	int barwidth = 70;
	int step = 1;
      */
      while( T->GetEntry( mevent++ ) ){ 
	
	//Cut on BBCal and HCal trigger coincidence
	double bbcal_time=-1000., hcal_time=-1000.;
	for(int ihit=0; ihit<hcalt::TDCTndata; ihit++){
	  if(hcalt::TDCT_id[ihit]==5) {
	    bbcal_time=hcalt::TDCT_tdc[ihit];
	    hTrigBBCal->Fill(bbcal_time);
	    //hTrigBBCal[i]->Fill(bbcal_time);
	    //cout << hTrigBBCal[i]->GetEntries() << endl;
	  }
	  /*
	    if(hcalt::TDCT_id[ihit]==0) {
	    hcal_time=hcalt::TDCT_tdc[ihit];
	    hTrigHCal->Fill(hcal_time);
	    }
	  */
	}
	/*
	  if( bbcal_time > -900. && hcal_time > -900. ) {
	  double diff = hcal_time - bbcal_time; 
	  hTrigDiff->Fill(diff);
	  //ATimeDist->Fill(hcalt::a_time[m]);
	  //cout << diff << endl;
	  }
	  
	*/
	/*
	  cout << "[";
	  int pos = barwidth*progress;
	  for( int i=0; i<barwidth; ++i ){
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
	    progress = (double)( ( mevent+1. )/Nevents );
	    
	    cout << "]" << int( progress*100 ) << "%\r";
	    cout.flush();
	    if( mevent%1000==0 ) step++;
	*/
      }
    
    
      //hTrigBBCal->Write(Form("hTrigBBCal_%d",i));
      if( hTrigBBCal->GetEntries()>0) hbbTrigShift->Fill(i,hTrigBBCal->GetMean());
      //hbbTrigShift->Fill(i,hTrigBBCal[i]->GetMean());
    
      T->Clear();
    
      //if( T->GetEntries()>0 ) cout << "Tree still has entries.." << endl;

      cout << "[";
      int pos = barwidth*progress;
      for( int i=0; i<barwidth; ++i ){
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

      progress = (double)( ( i+1.-runStart )/runTotal );
      cout << "]" << int( progress*100 ) << "%\r";
      cout.flush();
      step++;
    }
  }

  cout << endl << endl;
  fout->cd();
  //hTrigHCal->Write();
  //hTrigBBCal->Write();
  /*
  for( int i=0; i<runTotal; i++){
    hTrigBBCal[i]->Write();
  }
  */
  //hTrigDiff->Write();
  hbbTrigShift->Write();
  

  fout->Close();

  cout << "Wrote results to file: outfiles/bbTrigShift_" << date.c_str() << ".root" << endl;
  
  st->Stop();
  
  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;
  
}

