//SSeeds 10.12.21 - Production - One-pass code combining elements from Andrew Puckett's BigCal calibration script and PDatta/EFuchey BB calibration code used to calibrate HCal energy deposition by PMT module. Tested with Eric Fuchey's digitized g4sbs data and used during detector commissioning in GMn run group.

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

void hcal_enCal( const char *configfilename, int run = -1 ){
  
  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  // Start the chain for root files passed with config file

  if( !T ) {
    T = new TChain( "T" );
    T->SetBranchStatus( "*", 0 );
    T->SetBranchStatus( "sbs.hcal.*", 1 );
    T->SetBranchAddress( "bb.tr.chi2", hcalt::BBtr_chi2 );
    T->SetBranchAddress( "bb.tr.n", hcalt::BBtr_n );
    T->SetBranchAddress( "bb.tr.px", hcalt::BBtr_px );
    T->SetBranchAddress( "bb.tr.py", hcalt::BBtr_py );
    T->SetBranchAddress( "bb.tr.pz", hcalt::BBtr_pz );
    T->SetBranchAddress( "bb.tr.p", hcalt::BBtr_p );
    T->SetBranchAddress("sbs.hcal.clus.id",hcalt::cid);
    T->SetBranchAddress("sbs.hcal.clus.row",hcalt::crow);
    T->SetBranchAddress("sbs.hcal.clus.col",hcalt::ccol);
    T->SetBranchAddress("sbs.hcal.clus.e",hcalt::ce);
    T->SetBranchAddress("sbs.hcal.clus.eblk",hcalt::ceblk); // Highest energy block
    T->SetBranchAddress("sbs.hcal.clus.nblk",hcalt::cnblk);
    T->SetBranchAddress("sbs.hcal.nclus",hcalt::nclus); 
    T->SetBranchAddress("sbs.hcal.nblk",hcalt::nblk); // Total number of blocks in cell
    T->SetBranchAddress("sbs.hcal.clus_blk.id",hcalt::cblkid); // 1-kNcell index for each block
    T->SetBranchAddress("sbs.hcal.clus_blk.e",hcalt::cblke); // Array of block energies

    cout << "Opened up tree with nentries=" << T->GetEntries() << endl;
  }

  double E_e = 1.92; // Energy of beam (incoming electrons from accelerator)
  int minEventPerCell = 100; // Minimum number of scattered p in cell required to calibrate
  int maxEventPerCell = 4000; // Maximum number of scattered p events to contribute
  double highDelta = 0.1; // Minimum M(i,j)/b(i) factor allowed 
  double HCal_d = 14.5; // Distance to HCal from scattering chamber for comm1
  double HCal_th = 35.0; // Angle that the center of HCal is at  
  string HVfilePath = "path/to/HV/settings";
  string alphasFilePath = "path/to/alpha/parameters";
  string constPath = "path/to/output/constants";
  string HVoutPath = "path/to/output/HV/settings";
  string inConstPath = "path/to/input/constants";
  string ratioPath = "path/to/output/ratios";

  // Reading config file
  ifstream configfile(configfilename);
  TString currentline;
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      T->Add(currentline);
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
      }
      if( skey == "minEventPerCell" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	minEventPerCell = sval.Atof();
      }
      if( skey == "maxEventPerCell" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	maxEventPerCell = sval.Atof();
      }
      if( skey == "highDelta" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	highDelta = sval.Atof();
      }
      if( skey == "HCal_d" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	HCal_d = sval.Atof();
      }
      if( skey == "HCal_th" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	HCal_th = sval.Atof();
      }
      if( skey == "inConstPath" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        inConstPath = sval;
      }
      if( skey == "HVfilePath" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	HVfilePath = sval;
      }
      if( skey == "alphasFilePath" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	alphasFilePath = sval;
      }
      if( skey == "HVoutPath" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	HVoutPath = sval;
      }
      if( skey == "ratioPath" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	ratioPath = sval;
      }
      if( skey == "constPath" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	constPath = sval;
      }
    }
    delete tokens;
  }
  
  // Declare outfile
  TFile *fout = new TFile( "eCalOut.root", "RECREATE" );
  
  // Initialize vectors and arrays
  double gHV[kNrows][kNcols];
  double gAlphas[kNrows][kNcols];
  double gTargetHV[kNrows][kNcols];
  double gOldConst[kNcell];
  double ga_c[kNcell] = {0.0};
  double gRatio[kNcell] = {0.0}; //Old ratios from first iteration
  double oneBlock[kNcell] = {0.0};
  
  // Initialize histograms
  TH1D *h_W = new TH1D( "W", "W", 250, 0.7, 1.5 );
  h_W->GetXaxis()->SetTitle( "GeV" );
  TH1D *h_Q2 = new TH1D( "Q2", "Q2", 250, 0.5, 3.0 );
  h_Q2->GetXaxis()->SetTitle( "GeV" );
  TH1D *h_E_ep = new TH1D( "Scattered Electron Energy","E_ep", 500, 0.0, E_e*1.5 ); // Set limits based on energy conserv
  h_E_ep->GetXaxis()->SetTitle( "GeV" );
  TH1D *h_E_pp = new TH1D( "Scattered Proton Energy", "E_pp", 500, 0.0, E_e*1.5 ); // Set limits based on energy conserv
  h_E_pp->GetXaxis()->SetTitle( "GeV" );
  TH1D *h_KE_p = new TH1D( "Scattered Proton Kinetic Energy", "KE_pp", 500, 0.0, E_e*1.5 ); // Set limits based on energy conserv
  h_KE_p->GetXaxis()->SetTitle( "GeV" );
  TH1F *h_HCalCol_p = new TH1F( "HCalCol_p", "Projection of Scattered p on HCal Columns", kNcols, 0, kNcols );
  h_HCalCol_p->GetXaxis()->SetTitle( "column" );

  Long64_t Nevents = T->GetEntries();
  cout << "Opened up tree with nentries: " << Nevents << ".." << endl;

  // Read in HV settings for PMTs
  ifstream HVFile( HVfilePath );
  
  if( !HVFile ){
    cerr << "No HV settings from run present -> setFiles/HV_run" << run << ".txt expected." << endl;
    return 0;
  }
  
  cout << "Getting HV settings for each pmt.." << endl;
  int n1=0;
  double d1;
  int rval, cval;
  string line;
  
  while( getline( HVFile, line ) ){
    if( line.at(0) == '#' ) {
      continue;
    }
    stringstream ss( line );
    ss >> d1;
    rval = floor( n1/kNcols );
    cval = n1 % kNcols;
    gHV[rval][cval] = -d1; 

    cout << "HV" << n1 << " = " << -d1 << "." << endl;

    n1++;
  }
  
  // Read in alpha parameters for PMTs
  ifstream alphaFile( alphasFilePath );
  if( !alphaFile ){
    cerr << "No PMT alphas file present -> setFiles/alphas.txt expected." << endl;
    return 0;
  }

  cout << "Getting alpha parameters for each pmt.." << endl;
  n1=0;
  d1=0;
  string line2;
  
  while( getline( alphaFile, line2 ) ){
    if( line2.at( 0 )=='#' ) {
      continue;
    }
    stringstream ss( line2 );
    ss >> d1;
    rval = floor( n1/kNcols );
    cval = n1 % kNcols;
    gAlphas[rval][cval] = d1;

    cout << "alpha" << n1 << " = " << d1 << "." << endl;

    n1++;
  }
 
  // Read in previous constants for PMTs
  ifstream inConstFile( inConstPath );
  if( !inConstFile ){
    cerr << "No input constant file present -> setFiles/const.txt expected." << endl;
    return 0;
  }

  cout << "Getting previous calibration constants.." << endl;
  n1=0;
  d1=0;
  string line4;
  
  while( getline( inConstFile, line4 ) ){
    if( line4.at( 0 )=='#' ) {
      continue;
    }
    stringstream ss( line4 );
    ss >> d1;
    gOldConst[n1] = d1;

    cout << "const" << n1 << " = " << d1 << "." << endl;

    n1++;
  }

  // Read in previous constants for PMTs
  ifstream ratioFile( ratioPath );
  if( !ratioFile ){
    cerr << "No input constant file present -> setFiles/ratio.txt expected." << endl;
    return 0;
  }

  cout << "Getting previous calibration constants.." << endl;
  n1=0;
  d1=0;
  string line5;
  
  while( getline( ratioFile, line4 ) ){
    if( line4.at( 0 )=='#' ) {
      continue;
    }
    stringstream ss( line4 );
    ss >> d1;
    gRatio[n1] = d1;

    cout << "ratio" << n1 << " = " << d1 << "." << endl;

    n1++;
  }

  cout << "All parameters loaded." << endl;

  // Need to get electron energy and W for cuts on elastic events
  long mevent = 0;

  // e' momentum correction factor
  double eCorr = 1.04;
  
  cout << "Main loop over all data commencing.." << endl;
    
  TMatrixD Ma(kNcell,kNcell);
  TVectorD ba(kNcell);
  int NEV[kNcell] = {0};
  int NEV_oneblock[kNcell] = {0};

  while( T->GetEntry( mevent++ ) ){ 

    if( mevent%1000 == 0 ) cout << "Main loop: " << mevent << "/" << Nevents << endl;
    
    // Sort all tracks by lowest Chi^2 (highest confidence track)
    int track_tot = (int)hcalt::BBtr_n[0];
    int track = 0;
    double min = 1000.0; //Set arbitrarily high chi^2 to minimize with loop over tracks

    for( int elem=0; elem<track_tot; elem++ ){            
      if( hcalt::BBtr_chi2[elem] < min )
	track = elem;      
    }

    double E_ep = sqrt( pow(M_e,2) + pow(hcalt::BBtr_p[track],2) ); // Obtain the scattered electron energy
    h_E_ep->Fill( E_ep ); // Fill histogram

    double p_ep = eCorr*hcalt::BBtr_p[track]; // Obtain the magnitude of scattered electron momentum with correction factor
    double Q2 = 2*E_e*E_ep*( 1-(hcalt::BBtr_pz[track]/p_ep) ); // Obtain Q2 from beam energy, outgoing electron energy, and momenta
    h_Q2->Fill( Q2 ); // Fill histogram

    double nu = E_e-E_ep; // Obtain energy transfer
    double W2 = pow( M_p,2 )+2*M_p*nu-Q2; // Obtain W2 from Q2 and nu
    double W = 0; // Should be Mp for elastic events. Will use to estimate correction factor while GEM tracking ongoing
    if( W2>0 ){
      W = sqrt( W2 ); // Obtain W for real events by throwing out negative W2 solns
      h_W->Fill( W );
    }
    
    // Use the electron kinematics to predict the proton momentum assuming elastic scattering on free proton at rest:
    double E_pp = nu+M_p; // Get energy of the proton
    h_E_pp->Fill( E_pp ); // Fill histogram

    double p_p = sqrt( pow( E_pp,2 )-W2 ); // Magnitude of the scattered proton momentum

    //double KE_p = pow( p_p,2 )/(2*M_p);
    double KE_p = nu; //For elastics
    h_KE_p->Fill( KE_p );

    //Each component in the form p_p_z = p_p*cos(phi)
    double phi_p = TMath::ACos( ( E_e-hcalt::BBtr_pz[track])/p_p ) * 180.0 / PI;

    //Get the projection of each event onto the face of HCal - consider only X (column) components to avoid calculations with 48D48 field map. Should still be useful
    double relPosX = HCal_d*TMath::Tan( HCal_th-phi_p );
    //int relPosCol = std::floor( relPosX*HCal_th+(kNcols/2) ); // Return the column where the p is expected to land in HCal
    h_HCalCol_p->Fill( relPosX );

    double Wsigma = 0.02;
    if( W<1.1&&W>0.3 ){ //Primary cut on elastics, kept wide to 
    //if( W<(0.97+10*Wsigma)&&W>(0.97-10*Wsigma) ){ //Tighter cut on W for final constants )

      int nblk = hcalt::cnblk[0];
      memset(ga_c, 0, kNcell*sizeof(double)); // Reset the memory of a_c on each event  

      // Get energies with simplest scheme from clusters only
      for( int blk = 0; blk<nblk; blk++ ){
	int blkid = int(hcalt::cblkid[blk])-1; //-1 necessary since sbs.hcal.clus_blk.id ranges from 1 to 288 (different than just about everything else)
	
	bool edge_max = false;

	if(nblk==1) {
	  
	  oneBlock[blkid] += hcalt::cblke[blk]/KE_p; // Simple estimation of the coefficients assuming 100% energy deposition in one block. Will check against the more sophisticated version with chi^2 reduction.

	  //cout << "oneblock ratio" << oneBlock[blkid] << " blkid" << blkid << endl;

	  NEV_oneblock[blkid]++;
	}

	// Cut on edge blocks since much of the energy is lost out of the edges of the modules
	int row = blkid/kNcols;
	int col = blkid%kNcols;

	//cout << "row" << row << " col" << col << endl;

	if( row == 0 || row == 23 || col == 0 || col == 11 ){
	  edge_max = true;

	  //cout << "Edge on r" << row << " col" << col << endl;
	}

	if(edge_max==false){
	  //if( blkid==0 ) cout << "0/0 event on number " << mevent << endl;
	  //ga_c[blkid] += hcalt::cblke[blk]*gRatio[blkid]; //Multiply by the ratio (unitless) of old constants to simulate replay of data
	  ga_c[blkid] += hcalt::cblke[blk];
	  NEV[blkid]++;
	  
	  // Build the matrix as simply as possible
	  for(int icol = 0; icol<kNcell; icol++){
	    ba(icol)+= ga_c[icol];
	    for(int irow = 0; irow<kNcell; irow++){
	      Ma(icol,irow)+= ga_c[icol]*ga_c[irow]/KE_p;
	      
	    }
	  }
	} 
      }
    }
  }

  // Reject the bad cells
  int badcell[kNcell];
  double y[kNcell] = {0.0}; // Build for easy TGraphErrors build

  for(int i=0; i<kNcell; i++){
    badcell[i] = 0;
    y[i] = i;

    oneBlock[i] /= NEV_oneblock[i]; // 

    cout << "oneBlock[i]" << oneBlock[i] << "NEV_oneblock[i]" << NEV_oneblock[i] << "i" << i << endl;

    if( NEV[i] < minEventPerCell || Ma(i,i) < 0.1*ba(i) ){ 
      
      ba(i) = 1.0;  // set RHS vector for cell i to 1.0 
      Ma(i,i) = 1.0; //set diagonal element of Matrix M for cell i to 1.0 
      
      for(int j=0; j<kNcell; j++){
	if( j != i ){
	  
	  Ma(i,j) = 0.0; 
	  Ma(j,i) = 0.0;
	}
      }
      badcell[i] = 1;
      cout << "cell" << i << " bad" << endl;
      cout << "Number of events in bad cell =" << NEV[i] << endl;
      cout << "Matrix element/vector element ratio =" << Ma(i,i)/ba(i) << endl;
    }  
  }

  // Invert the matrix, solve for coefficients
  TMatrixD M_inv = Ma.Invert();
  TVectorD Coeff = M_inv*ba; // Stayes unmodified for reference
  double GCoeff[kNcell] = {0.0};
  double GCoeff_oneblock[kNcell] = {0.0};
  double HVCoeff[kNcell] = {0.0}; // Seperate for HV such that bad cells stay the same.

  for(int i=0; i<kNcell; i++){
    if(badcell[i]==0){
      GCoeff[i]=gOldConst[i]*Coeff[i]; // The new gain coefficient is the old coefficient multiplied by the gain factor that we just solved for. We will take these coefficents as the input for another set of data to estimate the error per channel.
      //GCoeff[i]=gOldConst[i]*Coeff[i]-gOldConst[i]*gRatio[i];
      GCoeff_oneblock[i]=gOldConst[i]*oneBlock[i];
      HVCoeff[i]=Coeff[i];
    }else{
      GCoeff[i]=gOldConst[i];
      HVCoeff[i]=1.0;
    }
  }

  double yErr[kNcell] = {0.0};
  TGraphErrors *ccgraph = new TGraphErrors( kNcell, y, GCoeff, yErr, yErr ); 
  ccgraph->GetXaxis()->SetLimits(0.0,288.0);
  TGraphErrors *ccgraph_oneBlock = new TGraphErrors( kNcell, y, GCoeff_oneblock, yErr, yErr ); 
  ccgraph->GetXaxis()->SetLimits(0.0,288.0);

  // Write out diagnostic histos
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

  h_HCalCol_p->Write( "HCalCol_p" );
  h_HCalCol_p->Draw( "AP" );

  ccgraph->SetTitle("Calibration Coefficients");
  ccgraph->GetXaxis()->SetTitle("Channel");
  ccgraph->GetXaxis()->SetTitle("Unitless");
  ccgraph->SetMarkerStyle(21); //Boxes
  ccgraph->Write("CONSTANTS");
  ccgraph->Draw("AP");

  ccgraph_oneBlock->SetTitle("Calibration Coefficients One Block");
  ccgraph_oneBlock->GetXaxis()->SetTitle("Channel");
  ccgraph_oneBlock->GetXaxis()->SetTitle("Unitless");
  ccgraph_oneBlock->SetMarkerStyle(21); //Boxes
  ccgraph_oneBlock->Write("CONSTANTS_ONEBLOCK");
  ccgraph_oneBlock->Draw("AP");

  // Declare outfiles
  ofstream ratio;
  ofstream ECal_HV;

  ofstream GainCoeff;
  ofstream GainCoeff_oneblock;

  ECal_HV.open( HVoutPath );
  ECal_HV << "Module" << '\t' << " HV" << endl;

  for( int i=0; i<kNcell; i++ ){   
    int r = i/kNcols;
    int c = i%kNcols;

    ECal_HV << r*kNcols+c+1 << '\t' << -gHV[r][c]*pow(HVCoeff[i],(1/gAlphas[r][c])) << endl;
  }

  ECal_HV.close();

  ratio.open( constPath );
  ratio << "HCal_ratio = " << endl;

  for( int i=0; i<kNcell; i++ ){   
    ratio << Coeff[i] << endl;
  }

  ratio.close();

  GainCoeff.open( "gainCoefficients.txt" );
  GainCoeff << "HCal_gainCoeff = " << endl;

  for( int i=0; i<kNcell; i++ ){   
    GainCoeff << GCoeff[i] << endl;
  }

  GainCoeff_oneblock.close();

  GainCoeff_oneblock.open( "gainCoefficientsONEBLOCK.txt" );
  GainCoeff_oneblock << "HCal_gainCoeff = " << endl;

  for( int i=0; i<kNcell; i++ ){   
    GainCoeff_oneblock << GCoeff_oneblock[i] << endl;
  }

  GainCoeff_oneblock.close();

  cout << "Calibration complete and constants written to file." << endl;

  st->Stop();

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

}


