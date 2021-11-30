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
void hcal_signal_centering( int run = -1 ){
  
  //Debug with run 11587

  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );
		  
  // Get the date
  string date = getDate();
  
  // Declare outfile
  //TFile *fout = new TFile( Form("signalCentering_%d.root", run, "RECREATE" ) );
  TFile *fout = new TFile( Form("outfiles/SignalCentering_%d.root",run), "RECREATE" );
   
  // Declare histograms and arrays
  TH1D *gHistos[kNrows][kNcols];
  TH1D *gATimeSpec[kNrows][kNcols];
  TH2D *Heatmap_rowcol = new TH2D( "Heatmap_rowcol", "", kNcols, 0.5, kNcols+0.5, kNrows, 0.5, kNrows+0.5 );
  TF1 *f1;
  bool gSaturated[kNrows][kNcols];
  int gEvt[kNrows][kNcols] = {0};

  //string inFile = Form("/adaqfs/home/a-onl/sbs/Rootfiles/hcal_gmn_%d*.root",run);

  if( !T ) {
    T = new TChain( "T" );
    T->Add(Form("/adaqfs/home/a-onl/sbs/Rootfiles/hcal_gmn_%d*.root",run));
    //T->Add( inFile );
    //T->Add( "/adaqfs/home/a-onl/sbs/Rootfiles/hcal_gmn_11587_-1.root" );
    //cout << "Loading branches from file: " << inFile << ".." << endl;
    T->SetBranchStatus( "*", 0 );
    T->SetBranchStatus( "sbs.hcal.*", 1 );
    T->SetBranchAddress( "sbs.hcal.a_p", hcalt::a_p );
    T->SetBranchAddress( "sbs.hcal.a_amp", hcalt::a_amp );
    T->SetBranchAddress( "sbs.hcal.a_amp_p", hcalt::a_amp_p );
    T->SetBranchAddress( "sbs.hcal.a_c", hcalt::a_c );
    T->SetBranchAddress( "sbs.hcal.adcrow", hcalt::row );
    T->SetBranchAddress( "sbs.hcal.adccol", hcalt::col );
    T->SetBranchAddress( "sbs.hcal.a_time", hcalt::a_time );
    T->SetBranchStatus( "bb.tdctrig.tdc", 1 );
    T->SetBranchStatus( "bb.tdctrig.tdcelemID", 1 );
    T->SetBranchStatus( "Ndata.bb.tdctrig.tdcelemID", 1 );

    T->SetBranchAddress("sbs.hcal.nsamps",hcalt::nsamps);
    T->SetBranchAddress("sbs.hcal.samps",hcalt::samps);
    T->SetBranchAddress("sbs.hcal.samps_idx",hcalt::samps_idx);

    T->SetBranchAddress( "bb.tdctrig.tdcelemID", hcalt::TDCT_id );
    T->SetBranchAddress( "bb.tdctrig.tdc", hcalt::TDCT_tdc );
    T->SetBranchAddress( "Ndata.bb.tdctrig.tdcelemID", &hcalt::TDCTndata );
    T->SetBranchStatus( "Ndata.sbs.hcal.adcrow", 1 );
    T->SetBranchAddress( "Ndata.sbs.hcal.adcrow", &hcalt::ndata );
    for(int r = 0; r < kNrows; r++) {
      for(int c = 0; c < kNcols; c++) {
        gHistos[r][c] = MakeHisto(r,c,totSample);
      }
    }
  }

  for( int r=0; r<kNrows; r++ ){
    for( int c=0; c<kNcols; c++ ){  
      int m = r*kNcols+c;
      gATimeSpec[r][c] = new TH1D( Form( "ADC Time Spect R%d C%d PMT%d", r, c, m ), Form( "ADC Time Spect R%d C%d PMT%d", r, c, m ), 4*totSample, 4*minSample, 4*maxSample ); //Each bin 4 ns wide
      gATimeSpec[r][c]->GetXaxis()->SetTitle( "ns" );
    }
  }
  
  Long64_t Nevents = T->GetEntries();
  cout << "Total events in tree: " << Nevents << ".." << endl;
  
  // Need to get electron energy and W for cuts on elastic events
  long mevent = 0;
  
  cout << "Cutting on BB/HCal timing coincidence and building ADC time histograms.." << endl;
  
  double progress = 0.;
  
  while( progress<1.0 ){
    
    int barwidth = 70;
    
    while( T->GetEntry( mevent++ ) ){ 
      
      //Cut on BBCal and HCal trigger coincidence
      double bbcal_time=0., hcal_time=0.;
      for(int ihit=0; ihit<hcalt::TDCTndata; ihit++){
	if(hcalt::TDCT_id[ihit]==5) bbcal_time=hcalt::TDCT_tdc[ihit];
	if(hcalt::TDCT_id[ihit]==0) hcal_time=hcalt::TDCT_tdc[ihit];
      }
      double diff = hcal_time - bbcal_time; 
      
      // Cut on BB/HCal trigger coincidence
      if( fabs(diff-510.)<95 ){
	
	int r,c,idx,n;
	double peak[kNrows][kNcols];
	
	// Clear old histograms
	for(r = 0; r < kNrows; r++) {
	  for(c = 0; c < kNcols; c++) {
	    gHistos[r][c]->Reset("ICES M");
	    peak[r][c] = 0.0;
	  }
	}
	
	for( int m = 0; m < hcalt::ndata; m++ ) { //Process event with m data 
	  r = hcalt::row[m];
	  c = hcalt::col[m];
	  
	  if( r < 0 || c < 0 ) {
	    cerr << "Error: row or col negative." << endl;
	    continue;
	  }
	  
	  if( r >= kNrows || c >= kNcols ) {
	    cerr << "Error: row or col outside of range." << endl;
	    continue;
	  }
	  
	  idx = hcalt::samps_idx[m];
	  n = hcalt::nsamps[m];
	  bool processed = false;
	  bool saturated = false;
	  
	  //cout << r << " " << c << " " << idx << " " << n << endl;


	  //Fill signal histogram from samps and mark saturated array if applicable
	  for(int s = minSample; s < maxSample && s < n; s++) {
	    processed = true;
	    gHistos[r][c]->SetBinContent(s+1-minSample,hcalt::samps[idx+s]);
	    if(peak[r][c]<hcalt::samps[idx+s])
	      peak[r][c]=hcalt::samps[idx+s];
	    if(peak[r][c]>2.2) {
	      saturated = true;
	    }
	  }
	  
	  //Write out saturated waveforms (one per channel max)
	  if( gSaturated[r][c]==false && saturated==true ){
	    fout->cd();
	    gHistos[r][c]->SetTitle(Form("Sample Saturated R%d C%d",r,c));
	    gHistos[r][c]->GetYaxis()->SetTitle("mV");
	    gHistos[r][c]->GetYaxis()->CenterTitle();
	    gHistos[r][c]->GetXaxis()->SetTitle("ADC Sample");
	    gHistos[r][c]->GetXaxis()->CenterTitle();
	    gHistos[r][c]->Write(Form("saturHisto_r%d_c%d",r,c));
	    gSaturated[r][c]=true;
	  }
	  
	  //cout << peak[r][c] << endl;
	  //cout << hcalt::a_p[m] << endl;

	  //T->Add(Form("/adaqfs/home/a-onl/sbs/Rootfiles/hcal_gmn_%d*.root",run));
	  //cout << gEvt[r][c] << endl;

	  //Write out 3 good sample histograms per channel
	  if( gEvt[r][c]<3 && peak[r][c]>0.03 && peak[r][c]<3000){
	    
	    //cout << "Good sample" << endl;
	   
	    fout->cd();
	    gHistos[r][c]->SetTitle(Form("Sample R%d C%d",r,c));
	    gHistos[r][c]->GetYaxis()->SetTitle("mV");
	    gHistos[r][c]->GetYaxis()->CenterTitle();
	    gHistos[r][c]->GetXaxis()->SetTitle("ADC Sample");
	    gHistos[r][c]->GetXaxis()->CenterTitle();
	    gHistos[r][c]->Write(Form("Histo%d_r%d_c%d",gEvt[r][c],r,c));
	    gEvt[r][c]++;
	    
	    gATimeSpec[r][c]->Fill( hcalt::a_time[m] );
	    
	    //cout << "gATimeSpec[" << r << "][" << c << "]=" << gATimeSpec[r][c]->GetEntries() << endl;
	  }
	  
	  //if( hcalt::a_time[m]>0 ) gATimeSpec[r][c]->Fill( hcalt::a_time[m] );
	  //cout << hcalt::a_time[m] << endl;
	  
	}		
      }
      
      cout << "[";
      int pos = barwidth*progress;
      for( int i=0; i<barwidth; ++i ){
	if( i<pos ) cout << "=";
	else if( i==pos ) cout << ">";
	else cout << " ";
      }
      progress = (double)( ( mevent+1. )/Nevents );
      
      cout << "]" << int( progress*100 ) << "%\r";
      cout.flush();
      
    }
  }
  
  cout << endl << endl;
  
  TH1D *ATimeDist = new TH1D( "Mean ADC time over all channels", "meanADCTime", 4*totSample, 4*minSample, 4*maxSample );
  
  for( int r=0; r<kNrows; r++){
    for( int c=0; c<kNcols; c++){
      
      if(gATimeSpec[r][c]->GetEntries()==0) 
	cout << "Empty a_time spectra at " << r << " " << c << "." << endl; 

      if(gATimeSpec[r][c]->GetEntries()>0){
	gATimeSpec[r][c]->Fit("gaus","Q");
	f1=gATimeSpec[r][c]->GetFunction("gaus");
	fout->cd();
	gATimeSpec[r][c]->Write(Form("ADC Time Spect R%d C%d FitMean%f FitSigma%f",r,c,f1->GetParameter(1),f1->GetParameter(2)));
	
	ATimeDist->Fill(f1->GetParameter(1));

      }
    }
  }
  
  fout->cd();
  ATimeDist->Write();

  fout->Close();

  cout << "Wrote results to file: outfiles/signalCentering_" << run << ".root" << endl;
  
  st->Stop();
  
  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;
  
}

