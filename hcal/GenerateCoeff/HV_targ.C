//SSeeds 10.12.21 - Production - Generate HV settings from calibration constants and print graphs

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
  string inConstPath = "goldGC_102221.txt";
  string HVoutPath = "goldHVOUT.txt";
  string errorPath = "goldErrorEstim.txt";
  
  // Declare outfile
  TFile *fout = new TFile( "HVTarg.root", "RECREATE" );
  
  // Initialize vectors and arrays
  double gHV[kNcell];
  double gAlphas[kNcell];
  double gTargetHV[kNcell];
  double gHVdiff[kNcell];
  double gConst[kNcell];
  double oldHVconst = 0.00175; //Results from a constant existing for cosmic calibration of the detector, but not in the replayed data file
  double errors[kNcell];
  double xErr[kNcell] = {0.0};

  double y[kNcell] = {0.0}; // For easy TGraphErrors build
  for(int i=0; i<kNcell; i++)
    y[i] = i;
  //int badModules = {5,19,20,31,32,77,78,103,149,150,284,21,22,23,24,26,27,33,35,73,74,75,76,84,109,131,132,145,146,147,148,154,158,169,178,204};


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

    gHV[n1] = -d1; 

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
    gAlphas[n1] = d1;

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
    gConst[n1] = d1;

    cout << "const " << n1 << " = " << d1 << "." << endl;

    n1++;
  }

  // Read in previous constants for PMTs
  ifstream errorFile( errorPath );
  if( !errorFile ){
    cerr << "No error file present -> setFiles/error.txt expected." << endl;
    return 0;
  }

  cout << "Getting estimated errors.." << endl;
  n1=0;
  d1=0;
  string line5;
  
  while( getline( errorFile, line5 ) ){
    if( line5.at( 0 )=='#' ) {
      continue;
    }
    stringstream ss( line5 );
    ss >> d1;
    errors[n1] = d1;

    cout << "error " << n1 << " = " << d1 << "." << endl;

    n1++;
  }

  cout << "All parameters loaded." << endl;

  for( int i=0; i<kNcell; i++ ){   
    bool err = false;

    if(errors[i]>1.5||errors[i]<0.5)
      err = true;

    if(errors[i]==-1) {
      gConst[i]=oldHVconst; //If not enough data, set old values
    }else{
      errors[i] *= gConst[i];
    }
            
    gTargetHV[i] = -gHV[i]*pow(gConst[i]/oldHVconst,(1/gAlphas[i]));
    gHVdiff[i] = -gHV[i]+gHV[i]*pow(gConst[i]/oldHVconst,(1/gAlphas[i]));

    if(err==true){
      gTargetHV[i] = -gHV[i];
      gHVdiff[i] = 10000;
    }

    if(err==true) cout << "ERROR " << gTargetHV[i] << " " << gHVdiff[i] << endl;
  }

  //TGraphErrors *ccgraph = new TGraphErrors( kNcell, y, gConst, xErr, errors ); 
  TGraphErrors *ccgraph = new TGraphErrors( kNcell, y, gConst, xErr, xErr ); //Leave errors out until they are understood 
  ccgraph->GetXaxis()->SetLimits(0.0,288.0);  
  ccgraph->SetTitle("Calibration Coefficients");
  ccgraph->GetXaxis()->SetTitle("Channel");
  ccgraph->SetMarkerStyle(21); //Boxes
  ccgraph->Write("cconst");

  //TMultiGraph *hvBadMod = new TMultigraph( "hvBadMod","hvBadMod" );

  //TGraphErrors *badMod = new TGraphErrors();


  TGraphErrors *hvgraph = new TGraphErrors( kNcell, y, gTargetHV, xErr, xErr ); //Graph differences between predicted HV and previous for comparison
  hvgraph->GetXaxis()->SetLimits(0.0,288.0);  
  hvgraph->SetTitle("Target HV");
  hvgraph->GetXaxis()->SetTitle("Channel");
  hvgraph->SetMarkerStyle(21); //Boxes
  hvgraph->Write("hv");

  TGraphErrors *hvdiffgraph = new TGraphErrors( kNcell, y, gHVdiff, xErr, xErr ); //Graph differences between predicted HV and previous for comparison
  hvdiffgraph->GetXaxis()->SetLimits(0.0,288.0);  
  hvdiffgraph->SetTitle("Delta HV");
  hvdiffgraph->GetXaxis()->SetTitle("Channel");
  hvdiffgraph->SetMarkerStyle(21); //Boxes
  hvdiffgraph->Write("hvdiff");

  // Print to console
  for( int r = 0; r<kNrows; r++){
    for( int c = 0; c<kNcols; c++){
      int i = r*kNcols+c;
      cout << -gHV[i] << "  ";
    }
    cout << endl;
  }
  
  cout << endl << endl;

  // Print to console
  for( int r = 0; r<kNrows; r++){
    for( int c = 0; c<kNcols; c++){
      int i = r*kNcols+c;
      cout << -gHV[i]*pow(gConst[i]/oldHVconst,(1/gAlphas[i])) << "  ";
    }
    cout << endl;
  }
  
  cout << endl << endl;

  // Print to console
  for( int r = 0; r<kNrows; r++){
    for( int c = 0; c<kNcols; c++){
      int i = r*kNcols+c;
      //cout << -gHV[i]+gHV[i]*pow(gConst[i]/oldHVconst,(1/gAlphas[i])) << "  ";

      
    }
    cout << endl;
  }
  
  cout << endl << endl;
  for(int i=0; i<kNcell; i++){
    if( i==0 ) cout << -gHV[i] << endl;
    if( -gHV[i]+gHV[i]*pow(gConst[i]/oldHVconst,(1/gAlphas[i])) < -50 ){
      cout << "Module" << i << " change" << -gHV[i]*pow(gConst[i]/oldHVconst,(1/gAlphas[i]))+gHV[i] << " TargetHV " << -gHV[i]*pow(gConst[i]/oldHVconst,(1/gAlphas[i])) << " oldHV " << -gHV[i] << endl;
    }
  }
  


  fout->Write();
  
  // Declare outfiles
  ofstream ECal_HV;
  ECal_HV.open( HVoutPath );
  ECal_HV << "Module" << '\t' << " HV" << endl;
  for(int i=0; i<kNcell; i++){
    int r = i/kNcols;
    int c = i%kNcols;
    ECal_HV << r*kNcols+c+1 << '\t' << gTargetHV[i] << endl;
  }
  ECal_HV.close();
  
  cout << "HV targets written to file." << endl;

  st->Stop();

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

  return 0;

}


