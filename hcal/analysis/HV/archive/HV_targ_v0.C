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

const int kNcell = 288; // Total number of HCal modules
const int kNrows = 24; // Total number of HCal rows
const int kNcols = 12; // Total number of HCal columns


int HV_targ( int run = -1 ){
  
  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  string HVfilePath = "oldHV.txt";
  string alphasFilePath = "alphas.txt";
  //string inConstPath = "hcal_eng_cal_gainCoeff_1.txt";
  string inConstPath = "setFiles/const.txt";
  string HVoutPath = "HVOUT.txt";
  
  // Declare outfile
  TFile *fout = new TFile( "HVTarg.root", "RECREATE" );
  
  // Initialize vectors and arrays
  double gHV[kNrows][kNcols];
  double gAlphas[kNrows][kNcols];
  double gTargetHV[kNrows][kNcols];
  double gOldConst[kNcell];
  
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

  cout << "All parameters loaded." << endl;


  // Print to console
  int cell = 0;
  for( int r = 0; r<kNrows; r++){
    for( int c = 0; c<kNcols; c++){
      cout << -gHV[r][c]*pow(gOldConst[cell],(1/gAlphas[r][c])) << "  ";
      cell++;
    }
    cout << endl;
  }
  
  TH1D *hHVtarg = new TH1D( "hHVtarg","HV Targets", 288, 0, 288 );

  // Declare outfiles
  ofstream ECal_HV;

  ECal_HV.open( HVoutPath );
  ECal_HV << "Module" << '\t' << " HV" << endl;

  for( int i=0; i<kNcell; i++ ){   
    int r = i/kNcols;
    int c = i%kNcols;

    gTargetHV[r][c] = -gHV[r][c]*pow(gOldConst[i],(1/gAlphas[r][c]));

    hHVtarg->Fill(i,gTargetHV[r][c]);

    ECal_HV << r*kNcols+c+1 << '\t' << -gHV[r][c]*pow(gOldConst[i],(1/gAlphas[r][c])) << endl;
  }

  ECal_HV.close();

  fout->Write();

  cout << "HV targets written to file." << endl;

  st->Stop();

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

  return 0;

}


