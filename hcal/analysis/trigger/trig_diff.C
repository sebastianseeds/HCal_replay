//SSeeds 11.12.21 - Production - Script to find coincidence signals between BB and HCal and calculate latency parameter at 
// ssh daq@enpcamsonne
// ssh intelsbshcal1
// open cfg/fadc250/hcal_ts.cnf
// change FADC250_W_OFFSET for both FADC250_CRATE intelsbshcal1.jlab.org and FADC250_CRATE intelsbshcal2.jlab.org

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

const int kNcell = 288; // Total number of HCal modules
const int kNrows = 24; // Total number of HCal rows
const int kNcols = 12; // Total number of HCal columns
const int minSample = 0.0;
const int maxSample = 40.0;
const int totSample = (maxSample-minSample); //Should be the total fADC window size with each samp = 4ns

void hcal_signal_centering();

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

// Create generic histogram function
TH1D* MakeHisto(int row, int col, int bins){
  TH1D *h = new TH1D(Form("h%02d%02d",row,col),Form("%d-%d",row+1,col+1),bins,minSample,maxSample);
  h->SetStats(0);
  h->SetLineWidth(2);
  return h;
}

// Main
void trig_diff( int run = -1 ){
  
  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );
		  
  // Get the date
  string date = getDate();
  
  // Declare outfile
  TFile *fout = new TFile( Form("outfiles/trig_diff_%d.root",run), "RECREATE" );
   
  // Declare histograms and arrays
  TH1D *hTrigDiff = new TH1D( "hTrigDiff","1190 HCal Trig - BBCal Trig (ns)", 2000, -1000, 1000 );
  TH1D *hTrigBBCal = new TH1D( "hTrigBBCal","1190 BBCal Trig (ns)", 2000, -1000, 1000 );
  TH1D *hTrigHCal = new TH1D( "hTrigHCal","1190 HCal Trig (ns)", 2000, -1000, 1000 );

  //TH1D *ATimeDist = new TH1D( "ADC time over all channels", "ADCTime", 4*totSample, 4*minSample, 4*maxSample );
  TF1 *f1;
  int gEvt[kNrows][kNcols] = {0};


  if( !T ) {
    T = new TChain( "T" );
    T->Add( Form("/adaqfs/home/a-onl/sbs/Rootfiles/hcal_trim_%d*.root",run) );
    //T->Add( Form("/adaqfs/home/a-onl/sbs/Rootfiles/hcal_gmn_%d_4100000_70.root",run) );
    //T->Add( Form("/adaqfs/home/a-onl/sbs/Rootfiles/hcal_gmn_%d_4100000_71.root",run) );
    //T->Add( Form("/adaqfs/home/a-onl/sbs/Rootfiles/hcal_gmn_%d_4100000_72.root",run) );
    //T->Add( Form("/adaqfs/home/a-onl/sbs/Rootfiles/hcal_gmn_%d_4100000_73.root",run) );
    //T->Add( Form("/adaqfs/home/a-onl/sbs/Rootfiles/hcal_gmn_%d_4100000_74.root",run) );
    //T->Add( Form("/adaqfs/home/a-onl/sbs/Rootfiles/hcal_trim_%d*.root",run) );
    T->SetBranchStatus( "*", 0 );
    T->SetBranchStatus( "sbs.hcal.*", 1 );
    T->SetBranchStatus( "bb.tdctrig.tdc", 1 );
    T->SetBranchStatus( "bb.tdctrig.tdcelemID", 1 );
    T->SetBranchStatus( "Ndata.bb.tdctrig.tdcelemID", 1 );

    T->SetBranchAddress( "sbs.hcal.tdcrow", hcalt::trow );
    T->SetBranchAddress( "sbs.hcal.tdccol", hcalt::tcol );    
    T->SetBranchAddress( "sbs.hcal.adcrow", hcalt::row );
    T->SetBranchAddress( "sbs.hcal.adccol", hcalt::col );
    T->SetBranchAddress( "sbs.hcal.a_time", hcalt::a_time );
    T->SetBranchAddress( "bb.tdctrig.tdcelemID", hcalt::TDCT_id );
    T->SetBranchAddress( "bb.tdctrig.tdc", hcalt::TDCT_tdc );
    T->SetBranchAddress( "Ndata.bb.tdctrig.tdcelemID", &hcalt::TDCTndata );
    T->SetBranchStatus( "Ndata.sbs.hcal.adcrow", 1 );
    T->SetBranchAddress( "Ndata.sbs.hcal.adcrow", &hcalt::ndata );
  }
  
  Long64_t Nevents = T->GetEntries();
  cout << "Total events in tree: " << Nevents << ".." << endl;
  
  // Need to get electron energy and W for cuts on elastic events
  long mevent = 0;
  
  cout << "Building timing and trigger histograms.." << endl;
  
  double progress = 0.;
  
  while( progress<1.0 ){
    
    int barwidth = 70;
    int step = 1;
    while( T->GetEntry( mevent++ ) ){ 
      
      //Cut on BBCal and HCal trigger coincidence
      double bbcal_time=-1000., hcal_time=-1000.;
      for(int ihit=0; ihit<hcalt::TDCTndata; ihit++){
	if(hcalt::TDCT_id[ihit]==5) {
	  bbcal_time=hcalt::TDCT_tdc[ihit];
	  hTrigBBCal->Fill(bbcal_time);
	}
	if(hcalt::TDCT_id[ihit]==0) {
	  hcal_time=hcalt::TDCT_tdc[ihit];
	  hTrigHCal->Fill(hcal_time);
	}
      }
      if( bbcal_time > -900. && hcal_time > -900. ) {
	double diff = hcal_time - bbcal_time; 
	hTrigDiff->Fill(diff);
	//ATimeDist->Fill(hcalt::a_time[m]);
	//cout << diff << endl;
      }
	  
      
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

    }
  }
  
  cout << endl << endl;
  fout->cd();
  hTrigHCal->Write();
  hTrigBBCal->Write();
  hTrigDiff->Write();
  fout->Close();

  cout << "Wrote results to file: outfiles/signalCentering_" << run << ".root" << endl;
  
  st->Stop();
  
  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;
  
}

