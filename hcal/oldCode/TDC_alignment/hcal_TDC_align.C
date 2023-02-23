//SSeeds 10.12.21 - Production - Code to fit TDC spectra, fit, extract mean, and give alignment parameters intended for db_sbs.hcal.dat

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
const double kdBlock = 0.152; // Width and height of each module including distance between

const double PI = TMath::Pi();
const double M_e = 0.00051;
const double M_p = 0.938272;
const double M_n = 0.939565;
const double c_light = 299792458.0;

void hcal_TDC_align();
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

void hcal_TDC_align( const char *configfilename, int run = -1 ){
  
  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );
  
  // Start the chain for root files passed with config file
  //TChain *C = new TChain("T");
  
  double E_e = 1.92; // Energy of beam (incoming electrons from accelerator)
  double current = 1.0; // Current of beam
  double W_mean = M_p; // With perfect optics, this should be true. Will read in until then by fitting W distribution on each run
  double W_sigma = 0.03; // Reasonable value by default, but will read in from W distribution
  int cosmic = 0; // Configure alignment script for cosmic signals
  double tdc_calib = 0.112;
  vector<double> zero_offset;

  if( !T ) {
    T = new TChain( "T" );
    T->SetBranchStatus( "*", 0 );
    T->SetBranchStatus( "sbs.hcal.*", 1 );
    T->SetBranchAddress( "sbs.hcal.a_p", hcalt::a_p );
    T->SetBranchAddress( "sbs.hcal.tdc", hcalt::tdc );
    T->SetBranchAddress( "sbs.hcal.a_amp", hcalt::a_amp );
    T->SetBranchAddress( "sbs.hcal.a_c", hcalt::a_c );
    T->SetBranchAddress( "sbs.hcal.adcrow", hcalt::row );
    T->SetBranchAddress( "sbs.hcal.adccol", hcalt::col );
    T->SetBranchAddress( "sbs.hcal.a_time", hcalt::a_time );
    //T->SetBranchAddress( "sbs.hcal.e", hcalt::e ); //En
    //T->SetBranchAddress( "sbs.hcal.eblk", hcalt::e ); //En
    T->SetBranchStatus( "bb.tdctrig.tdcelemID", 1 );
    T->SetBranchStatus( "Ndata.bb.tdctrig.tdcelemID", 1 );
    T->SetBranchStatus( "bb.sh.nclus", 1 );

    T->SetBranchAddress( "bb.tr.chi2", hcalt::BBtr_chi2 );
    T->SetBranchAddress( "bb.tr.n", hcalt::BBtr_n );
    T->SetBranchAddress( "bb.tr.px", hcalt::BBtr_px );
    T->SetBranchAddress( "bb.tr.py", hcalt::BBtr_py );
    T->SetBranchAddress( "bb.tr.pz", hcalt::BBtr_pz );
    T->SetBranchAddress( "bb.tr.p", hcalt::BBtr_p );
    T->SetBranchAddress("sbs.hcal.clus.id",hcalt::cid);
    T->SetBranchAddress("sbs.hcal.clus.row",hcalt::crow);
    T->SetBranchAddress("sbs.hcal.clus.col",hcalt::ccol);
    T->SetBranchAddress("sbs.hcal.clus.atime",hcalt::catime); //a_time for main block in each cluster
    T->SetBranchAddress("sbs.hcal.clus_blk.atime",hcalt::cblkatime); //a_time for each block in main cluster

    T->SetBranchAddress("sbs.hcal.clus_blk.row",hcalt::cbrow);
    T->SetBranchAddress("sbs.hcal.clus_blk.col",hcalt::cbcol);
    T->SetBranchAddress("sbs.hcal.clus_blk.atime",hcalt::cblkatime);
    T->SetBranchAddress("sbs.hcal.clus.e",hcalt::ce);
    T->SetBranchAddress("sbs.hcal.clus.eblk",hcalt::ceblk); // Highest energy block
    T->SetBranchAddress("sbs.hcal.clus.nblk",hcalt::cnblk);
    T->SetBranchAddress("sbs.hcal.nclus",hcalt::nclus);
    T->SetBranchAddress("sbs.hcal.nblk",hcalt::nblk);
    T->SetBranchAddress("sbs.hcal.clus_blk.id",hcalt::cblkid);
    T->SetBranchAddress("sbs.hcal.clus_blk.e",hcalt::cblke); // Array of block energies
    //T->SetBranchAddress("sbs.hcal.clus.atimeblk",hcalt::atimeblk); //a_time for main block in main cluster 
    //T->SetBranchAddress("sbs.hcal.atime",hcalt::atimeblk); // atime for block with highest energy in highest energy cluster
    T->SetBranchAddress( "bb.tdctrig.tdcelemID", hcalt::TDCT_id );
    T->SetBranchAddress( "bb.tdctrig.tdc", hcalt::TDCT_tdc );
    T->SetBranchAddress( "Ndata.bb.tdctrig.tdcelemID", &hcalt::TDCTndata );
    T->SetBranchStatus( "Ndata.sbs.hcal.adcrow", 1 );
    T->SetBranchAddress( "Ndata.sbs.hcal.adcrow", &hcalt::ndata );
    //cout << "Opened up tree with nentries=" << T->GetEntries() << endl;
  }
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
      if( skey == "current" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	current = sval.Atof();
	cout << "Beam current: " << current << endl;
      }
      if( skey == "W_mean" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	W_mean = sval.Atof();
	cout << "W mean: " << W_mean << endl;
      }
      if( skey == "W_sigma" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	W_sigma = sval.Atof();
	cout << "W std dev: " << W_sigma << endl;
      }
      if( skey == "cosmic" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	cosmic = sval.Atoi();
	cout << "Cosmic data? (0 = false; 1 = true): " << cosmic << endl;
      }
      if( skey == "tdc_calib" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	tdc_calib = sval.Atof();
	cout << "Operating TDC calibration constant: " << tdc_calib << endl;
      }
      delete tokens;
    }
  }

  // Get the date
  string date = getDate();
  
  // Declare outfile

  TFile *fout = new TFile( Form("outfiles/TDCalign%i_%s.root", cosmic, date.c_str()), "RECREATE" );

  // Initialize histograms
  TH1D *h_W = new TH1D( "W", "W", 250, 0.0, E_e*1.5 );
  h_W->GetXaxis()->SetTitle( "GeV" );
  TH1D *h_Q2 = new TH1D( "Q2", "Q2", 250, 0.5, E_e*1.5 );
  h_Q2->GetXaxis()->SetTitle( "GeV" );
  TH1D *h_E_ep = new TH1D( "Scattered Electron Energy","E_ep", 500, 0.0, E_e*1.5 ); // Set limits based on energy conservation
  h_E_ep->GetXaxis()->SetTitle( "GeV" );
  TH1D *h_E_pp = new TH1D( "Scattered Proton Energy", "E_pp", 500, 0.0, E_e*1.5 ); // Set limits based on energy conservation
  h_E_pp->GetXaxis()->SetTitle( "GeV" );
  TH1D *h_KE_p = new TH1D( "Scattered Proton Kinetic Energy", "KE_pp", 500, 0.0, E_e*1.5 ); // Set limits based on energy conservation
  h_KE_p->GetXaxis()->SetTitle( "GeV" );
  h_W->GetXaxis()->SetTitle( "GeV" );

  cout << "Global cut is: " << globalcut << ".." << endl;
  
  // Build ADC histograms - Set integrated ADC spectrum histogram limits 
  // Build histograms
  TH1F *gTDCSpec[kNrows][kNcols];
  TF1 *f1;

  for( int r=0; r<kNrows; r++ ){
    for( int c=0; c<kNcols; c++ ){  
      int m = r*kNcols+c;
      gTDCSpec[r][c] = new TH1F( Form( "ADC Time Spect R%d C%d PMT%d", r, c, m ), Form( "ADC Time Spect R%d C%d PMT%d", r, c, m ), 200, -300, 100 );
      gTDCSpec[r][c]->GetXaxis()->SetTitle( "ns" );
    }
  }
  
  Long64_t Nevents = T->GetEntries();
  cout << "Total events in tree: " << Nevents << ".." << endl;
  
  // Need to get electron energy and W for cuts on elastic events
  long mevent = 0;
  
  // e' momentum correction factor
  double eCorr = 1.04;
  
  if( cosmic==0 ) cout << "Loop cutting on elastics and plotting HCal TDC values commencing.." << endl;
  if( cosmic==0 ){
    while( T->GetEntry( mevent++ ) ){ 
    
      if( mevent%1000 == 0 ) cout << "Loop: " << mevent << " / " << Nevents << endl;
    
      // Sort all tracks by lowest Chi^2 (highest confidence track)
      int track_tot = int(hcalt::BBtr_n[0]);
      int track = 0;
      double min = 1000.0; //Set arbitrarily high chi^2 to minimize with loop over tracks

      for( int elem=0; elem<track_tot; elem++ ){            
	if( hcalt::BBtr_chi2[elem] < min )
	  track = elem;      
      }
    
      TVector3 lVep( hcalt::BBtr_px[track], hcalt::BBtr_py[track], hcalt::BBtr_pz[track] ); // Construct the scattered electron momentum vector from the tree (GEM reconstructed)
    
      double E_ep = sqrt( pow(M_e,2) + pow(hcalt::BBtr_p[track],2) ); // Obtain the scattered electron energy
      h_E_ep->Fill( E_ep ); // Fill histogram

      double p_ep = eCorr*hcalt::BBtr_p[track]; // Obtain the magnitude of scattered electron momentum with correction factor

      double theta_e = acos( ( hcalt::BBtr_pz[track])/p_ep ) * 180.0 / PI; // Obtain the electron's scattering angle in degrees
    
      double phi_e = TMath::ASin( hcalt::BBtr_py[track]/( p_ep*TMath::Sin( theta_e ) ) );
    
      double Q2 = 2*E_e*E_ep*( 1-(hcalt::BBtr_pz[track]/p_ep) ); // Obtain Q2 from beam energy, outgoing electron energy, and momenta
      h_Q2->Fill( Q2 ); // Fill histogram

      double nu = E_e-E_ep; // Obtain energy transfer

      double W2 = pow( M_p,2 )+2*M_p*nu-Q2; // Obtain W2 from Q2 and nu

      double W = 0;
      if( W2>0 ){
	W = sqrt( W2 ); // Obtain W for real events by throwing out negative W2 solns
	h_W->Fill( W );
      }
    
      // Use the electron kinematics to predict the proton momentum assuming elastic scattering on free proton at rest:

      double E_pp = nu+M_p; // Get energy of the proton
      h_E_pp->Fill( E_pp ); // Fill histogram

      double p_p = sqrt( pow( E_pp,2 )-W2 ); // Magnitude of the scattered proton

      //double KE_p = pow( p_p,2 )/(2*M_p);
      double KE_p = nu; //For elastics
      h_KE_p->Fill( KE_p );

      // Each component in the form p_p_z = p_p*cos(phi)
      double phi_p = TMath::ACos( ( E_e-hcalt::BBtr_pz[track])/p_p ) * 180.0 / PI;

      //Cut on BBCal and HCal trigger coincidence
      double bbcal_time=0., hcal_time=0.;
      for(int ihit=0; ihit<hcalt::TDCTndata; ihit++){
	if(hcalt::TDCT_id[ihit]==5) bbcal_time=hcalt::TDCT_tdc[ihit];
	if(hcalt::TDCT_id[ihit]==0) hcal_time=hcalt::TDCT_tdc[ihit];
      }
      double diff = hcal_time - bbcal_time; 
      
      // Liberal cut on W (W=Mp for elastics)
      if( W < (W_mean+W_sigma) && W > (W_mean-W_sigma) && fabs(diff-510.)<20 ){

	for( int m = 0; m < hcalt::ndata; m++ ) { //Process event with m data 
	  int r = hcalt::row[m];
	  int c = hcalt::col[m];
	  if( r < 0 || c < 0 ) {
	    cerr << "Error: row or col negative." << endl;
	    continue;
	  }
	
	  if( r >= kNrows || c >= kNcols ) {
	    cerr << "Error: row or col outside of range." << endl;
	    continue;
	  }

	  // if(globalcut){
	  if( hcalt::ce[0]>0 && hcalt::tdc[m]<1000 ){
	    gTDCSpec[r][c]->Fill( hcalt::tdc[m] );
	    
	    // cout << "Good event " << r << " " << c << " " << hcalt::tdc[m] << endl;
	    
	  }
	}		
      }
    }
  }

  if( cosmic==0 ){
    fout->cd();
    
    h_W->Write( "W" );
    h_W->Draw( "AP" );
    
    h_Q2->Write( "Q2" );
    h_Q2->Draw( "AP" );
    
    h_E_ep->Write( "E_ep" );
    h_E_ep->Draw( "AP" );
    
    h_KE_p->Write( "KE_p" );
    h_KE_p->Draw( "AP" );
    
    h_E_pp->Write( "E_pp" );
    h_E_pp->Draw( "AP" );
  }

  if( cosmic==1 ){
    cout << "Loop on cosmic data commencing.." << endl;

    while( T->GetEntry( mevent++ ) ){ 

      if( mevent%1000 == 0 ) cout << "Loop: " << mevent << "/" << Nevents << endl;

      for( int m = 0; m < hcalt::ndata; m++ ) { //Process event with m data 
	int r = hcalt::row[m];
	int c = hcalt::col[m];
	if( r < 0 || c < 0 ) {
	  cerr << "Error: row or col negative." << endl;
	  continue;
	}

	if( r >= kNrows || c >= kNcols ) {
	  cerr << "Error: row or col outside of range." << endl;
	  continue;
	}
	
	//Require that the tdc is not off (1e38 registered when tdc off) and that some event occurs in HCal
	if( hcalt::tdc[m]<1000&&hcalt::ce[0]>0 ){ 
	  
	  gTDCSpec[r][c]->Fill( hcalt::tdc[m] );
	  
	  //cout << hcalt::tdc[m] << endl;

	}
      }
    }
  }
  
  vector<double> tdc_res;

  fout->cd();
  for( int r=0; r<kNrows; r++){
    for( int c=0; c<kNcols; c++){
      
      if(gTDCSpec[r][c]->GetEntries()>0){
	gTDCSpec[r][c]->Fit("gaus","Q");
	f1=gTDCSpec[r][c]->GetFunction("gaus");
	fout->cd();
	gTDCSpec[r][c]->Write(Form("TDC Spect R%d C%d FitMean%f FitSigma%f",r,c,f1->GetParameter(1),f1->GetParameter(2)));
	
	zero_offset.push_back(f1->GetParameter(1));
	tdc_res.push_back(f1->GetParameter(2));
	//gTDCSpec[r][c]->Write();
      }
    }
  }
  
  fout->Close();

  ofstream TDCoffsets;
  TDCoffsets.open( "outfiles/TDCoffsets.txt" );
  int cell = 0;
  for( int r=0; r<kNrows; r++ ){
    for( int c=0; c<kNcols; c++ ){
      TDCoffsets << (zero_offset[cell]-zero_offset[42])/tdc_calib << "  "; //Align to PMT 42 (for no reason whatsoever) and correct for raw TDC value with tdc calibration constant
      cell++;
    }
    TDCoffsets << endl;
  }

  TDCoffsets.close();

  cout << "Wrote results to file: outfiles/TDCalign" << cosmic << "_" << date.c_str() << ".root" << endl;
  
  st->Stop();
  
  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;
  
}

