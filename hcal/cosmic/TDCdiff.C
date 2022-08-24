//SSeeds 8.24.22 - Cosmic - Assumes that the ADC and TDC offsets for replayed data are set to zero. This script gets the absolute fADC time for each channel and the absolute difference between TDC time for each channel and the reference time.

#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <TChain.h>
#include <TSystem.h>
#include <TStopwatch.h>
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

const Int_t kNcell = 288; // Total number of HCal modules
const Int_t kNrows = 24; // Total number of HCal rows
const Int_t kNcols = 12; // Total number of HCal columns
const Int_t maxTdcChan = 10; // Set to accomodate original 5 TDCTrig channels with buffer
const Int_t maxTracks = 1000; // Reasonable limit on tracks to be stored per event

const Double_t PI = TMath::Pi();

const Double_t TDCCalib = 0.112;

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

// Inputs: run (run number to be analyzed), fresh (toggle to check previous ADC/TDC offsets or not, 1=no)
void TDCdiff( Int_t run = -1, Int_t fresh=1 ){
  
  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  // Get the date
  string date = getDate();
  
  // Declare Chain for many root files
  TChain *C = new TChain("T");

  // Declare arrays to hold old offset parameters
  Double_t oldADCtoffsets[kNcell]={0.};
  Double_t oldTDCoffsets[kNcell]={0.};

  cout << endl;

  // Loads data file. WARNING: Must be configured to point to data on new environments
  TString datafile = Form("/lustre19/expphy/volatile/halla/sbs/seeds/rootfiles/hcal_general_%i*",run);
  C->Add(datafile);
  
  // Path for previous ADCt and TDC constants
  string inConstPath;
  inConstPath = "/w/halla-scshelf2102/sbs/seeds/SBS-replay/DB/db_sbs.hcal.dat";

  // Reading ADC and TDC timing offsets from database if fresh=0
  if( fresh==0 ){
    cout << "Loading previous offsets.." << endl;
    ifstream inConstFile( inConstPath );
    if( !inConstFile ){
      cerr << endl << "ERROR: No input constant file present -> path to db_sbs.hcal.dat expected." << endl;
      return 0;
    }
    
    Int_t n1=0;
    Double_t d1=0;
    string line;
    bool skip_line = true;
    bool skip_one_line = true;
    bool pass_first_cond = false;
    bool pass_all_cond = false;
    
    while( getline( inConstFile, line ) ){
      
      if( pass_first_cond && n1==( kNcell ) ) break;
      
      TString Tline = (TString)line;
      
      if( Tline.BeginsWith("sbs.hcal.adc.timeoffset") && skip_line==true ) skip_line = false;
      if( Tline.BeginsWith("sbs.hcal.tdc.offset") && skip_line==true ) skip_line = false;
      
      if( skip_line==false && skip_one_line==true ){
	skip_one_line = false;
	continue;
      }
      
      if( n1==( kNcell ) && pass_first_cond==false ){
	skip_line = true;
	skip_one_line = true;
	n1=0;
	pass_first_cond = true;
      }
      
      if( skip_line==false && skip_one_line==false ){
	istringstream iss( line );
	while( iss >> d1 ){
	  if( pass_first_cond==false ){
	    oldADCtoffsets[n1] = d1;
	  }else{
	    oldTDCoffsets[n1] = d1;
	  }
	  n1++;
	}
      }
    }
    
    cout << endl << endl << "Old TDC offsets: " << endl;
    
    for( Int_t r=0; r<kNrows; r++){
      for( Int_t c=0; c<kNcols; c++){
	Int_t i = r*kNcols+c;
	cout << oldTDCoffsets[i] << " ";
      }
      cout << endl;
    }
    
    cout << endl << endl << "Old ADC time offsets: " << endl;
    
    for( Int_t r=0; r<kNrows; r++){
      for( Int_t c=0; c<kNcols; c++){
	Int_t i = r*kNcols+c;
	cout << oldADCtoffsets[i] << " ";
      }
      cout << endl;
    }
    
    cout << endl << endl << "Setup parameters loaded." << endl;
  }

  // Declare general detector parameters
  Double_t TDCT_id[maxTdcChan], TDCT_tdc[maxTdcChan]; 
  Int_t TDCTndata;
  UInt_t runI = -1; 
  UInt_t runN = -1;
  ULong64_t runT;
  Int_t tFitMin = 5;

  Double_t HCALx, HCALy, HCALe;
  Double_t HCALblktdc[kNcell], HCALblka[kNcell];
  Double_t HCALtdc[kNcell], HCALatime[kNcell];

  Double_t crow, ccol, nblk;
  Double_t cblkid[kNcell], cblke[kNcell];

  Double_t tHODO;
  Double_t nClusHODO;
  
  // Declare root tree variables and set values to memory locations in root file
  C->SetBranchStatus( "*", 0 );

  C->SetBranchStatus( "sbs.hcal.x", 1 );
  C->SetBranchStatus( "sbs.hcal.y", 1 );
  C->SetBranchStatus( "sbs.hcal.e", 1 );
  C->SetBranchStatus( "sbs.hcal.tdc", 1 );
  C->SetBranchStatus( "sbs.hcal.a_time", 1 );
  C->SetBranchStatus( "sbs.hcal.rowblk", 1 );
  C->SetBranchStatus( "sbs.hcal.colblk", 1 );
  C->SetBranchStatus( "sbs.hcal.nblk", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.tdctime", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.atime", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.id", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.e", 1 );
  C->SetBranchStatus( "bb.tdctrig.tdc", 1 );
  C->SetBranchStatus( "bb.tdctrig.tdcelemID", 1 );
  C->SetBranchStatus( "Ndata.bb.tdctrig.tdcelemID", 1 );
  C->SetBranchStatus( "bb.hodotdc.clus.tmean", 1 );
  C->SetBranchStatus( "bb.hodotdc.nclus", 1 );
  C->SetBranchStatus( "fEvtHdr.fRun", 1 );
  C->SetBranchStatus( "fEvtHdr.fEvtTime", 1 );
  C->SetBranchStatus( "fEvtHdr.fEvtNum", 1 );

  C->SetBranchAddress( "sbs.hcal.x", &HCALx );
  C->SetBranchAddress( "sbs.hcal.y", &HCALy );
  C->SetBranchAddress( "sbs.hcal.e", &HCALe );
  C->SetBranchAddress( "sbs.hcal.tdc", HCALtdc );
  C->SetBranchAddress( "sbs.hcal.a_time", HCALatime );
  C->SetBranchAddress( "sbs.hcal.rowblk", &crow );
  C->SetBranchAddress( "sbs.hcal.colblk", &ccol );
  C->SetBranchAddress( "sbs.hcal.nblk", &nblk );
  C->SetBranchAddress( "sbs.hcal.clus_blk.id", cblkid );
  C->SetBranchAddress( "sbs.hcal.clus_blk.e", cblke );
  C->SetBranchAddress( "sbs.hcal.clus_blk.tdctime", HCALblktdc );
  C->SetBranchAddress( "sbs.hcal.clus_blk.atime", HCALblka );
  C->SetBranchAddress( "bb.tdctrig.tdcelemID", TDCT_id );
  C->SetBranchAddress( "bb.tdctrig.tdc", TDCT_tdc );
  C->SetBranchAddress( "Ndata.bb.tdctrig.tdcelemID", &TDCTndata );  
  C->SetBranchAddress( "bb.hodotdc.clus.tmean", &tHODO );
  C->SetBranchAddress( "bb.hodotdc.nclus", &nClusHODO );
  C->SetBranchAddress( "fEvtHdr.fRun", &runI );
  C->SetBranchAddress( "fEvtHdr.fEvtTime", &runT );
  C->SetBranchAddress( "fEvtHdr.fEvtNum", &runN );

  cout << "Tree variables linked." << endl;

  // Define a clock to check macro processing time
  TStopwatch *sw = new TStopwatch();
  sw->Start();
  
  // Declare outfile
  TString outputfilename = "TDCdiff_out.root";
  TFile *fout = new TFile( outputfilename, "RECREATE" );
  
  // Initialize vectors and arrays
  Double_t TDCoffsets[kNcell] = {0.0};
  Double_t ADCtoffsets[kNcell] = {0.0};
  Double_t TDCsig[kNcell] = {0.0};
  Double_t ADCtsig[kNcell] = {0.0};
  Double_t TvsE[kNcell] = {0.0};
  Double_t AtvsE[kNcell] = {0.0};

  // Initialize histograms
  TH1D *htdc[kNcell];
  TH1D *hatime[kNcell];
  for( int i=0; i<kNcell; i++ ){
    htdc[i] = new TH1D(Form("htdc_bl%d",i),";TDC_{HCAL} (ns)",800,-400,400);
    hatime[i] = new TH1D(Form("hatime_bl%d",i),";ADCt_{HCAL} (ns)",180,0,180);
  }

  cout << "Variables and histograms defined." << endl;

  // Set long int to keep track of total entries
  Long64_t Nevents = C->GetEntries();
  UInt_t run_number = 0;

  cout << endl << "All parameters loaded and initialization complete." << endl << endl;
  cout << "Opened up TChain with nentries: " << C->GetEntries() << "." << endl << endl;

  //Loop over events
  cout << "Main loop over all data commencing.." << endl;
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
	timeremains = timekeeper*( Double_t(Nevents)/Double_t(nevent) - 1. ); 
      sw->Reset();
      sw->Continue();
      
      progress = (Double_t)((nevent+1.)/Nevents);
      cout << "] " << int(progress*100.) << "%, time remaining: " << int(timeremains/60.) << "m \r";
      cout.flush();
      
      C->GetEntry( nevent ); 
		
      if( run_number!=runI ){
	run_number=runI;
	cout << "Now analyzing run " << run_number << "." << endl;
      }
	
      //Cut on BBCal and HCal trigger coincidence
      Double_t bbcal_time=0., hcal_time=0.;
      for(int ihit=0; ihit<TDCTndata; ihit++){
	if(TDCT_id[ihit]==5) bbcal_time=TDCT_tdc[ihit];
	if(TDCT_id[ihit]==0) hcal_time=TDCT_tdc[ihit];
      }
      Double_t diff = hcal_time - bbcal_time; 

      //Fill histograms. Note that this assumes that atime and tdc are fixed size arrays!
      for( int i=0; i<kNcell; i++ ){

	if( HCALatime[i]>0 ) hatime[i]->Fill( HCALatime[i] );
	if( HCALtdc[i]<10000 ) htdc[i]->Fill( HCALtdc[i] );
      }
    }
  }
  
  cout << endl;

  //Fit timing data and extract offsets
  for(int i=0; i<kNcell; i++){

    int r = (i)/kNcols;
    int c = (i)%kNcols;

    TF1 *f1;
    TF1 *f2;

    if( htdc[i]->GetEntries()>tFitMin ){
      htdc[i]->Fit("gaus","Q","",-400,400);
      f1=htdc[i]->GetFunction("gaus");
      htdc[i]->SetTitle(Form("TDCFitMean:%f",f1->GetParameter(1)));
      TDCoffsets[i] = f1->GetParameter(1);
      TDCsig[i] = f1->GetParameter(2);
      //Double_t cs = f1->GetChisquare();
      htdc[i]->SetName(Form("htdc_bl%d_r%d_c%d_mean%f",i,r,c,f1->GetParameter(1)));	
    }

    if( hatime[i]->GetEntries()>tFitMin ){
      hatime[i]->Fit("gaus","Q","",0,150);
      f2=hatime[i]->GetFunction("gaus");
      hatime[i]->SetTitle(Form("ADCtFitMean:%f",f2->GetParameter(1)));
      ADCtoffsets[i] = f2->GetParameter(1);
      ADCtsig[i] = f2->GetParameter(2);
      //Double_t cs = f2->GetChisquare();
      hatime[i]->SetName(Form("hatime_bl%d_r%d_c%d_mean%f",i,r,c,f2->GetParameter(1)));
    }
  }

  //Make canvas to hold all fits for comparison to HCal geometry
  TCanvas *TDC_top = new TCanvas("TDC_top","TDC_top",1600,1200);
  TCanvas *TDC_bot = new TCanvas("TDC_bot","TDC_bot",1600,1200);
  TCanvas *ADCt_top = new TCanvas("ADCt_top","ADCt_top",1600,1200);
  TCanvas *ADCt_bot = new TCanvas("ADCt_bot","ADCt_bot",1600,1200);

  TDC_top->Divide(12,12);
  TDC_bot->Divide(12,12);
  ADCt_top->Divide(12,12);
  ADCt_bot->Divide(12,12);

  gStyle->SetOptStat(0);
  for(int i=0; i<kNcell; i++){
    TDC_top->cd(i+1);
    if( i>=144 ){
      TDC_bot->cd(i-143);
      gStyle->SetOptStat(0);
    }
    if(htdc[i]->GetEntries()<tFitMin){
      htdc[i]->SetAxisColor(2);
    }else{
      htdc[i]->SetAxisColor(1);
    }
    htdc[i]->Draw();
  }

  gStyle->SetOptStat(0);
  for(int i=0; i<kNcell; i++){
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

  //Write out diagnostic histos and print to console
  fout->Write();

  cout << endl << endl << "TDC Offsets: " << endl << endl;
  
  Int_t cell = 0;
  for( int r=0; r<kNrows; r++ ){
    for( int c=0; c<kNcols; c++ ){
      cout << TDCoffsets[cell] << " ";
      cell++;
    }
    cout << endl;
  }

  cout << endl << endl << "ADC Time Offsets: " << endl << endl;

  cout << endl << endl << "Total events analyzed: " << Nevents << "." << endl << endl;

  cout << "Timing offset analysis complete." << endl;

  st->Stop();

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

}
