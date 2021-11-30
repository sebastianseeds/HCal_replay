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


int Error( int run = -1 ){
  
  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  //string firstHalfPath = "goldFirstHalf.txt";
  //string secondHalfPath = "goldSecondHalf.txt";
  //string errorOutPath = "goldErrorEstim.txt";

  string firstHalfPath = "halfOne.txt";
  string secondHalfPath = "halfTwo.txt";
  string errorOutPath = "error.txt";


  // Initialize vectors and arrays
  double gfirstHalf[kNcell];
  double gsecondHalf[kNcell];

  
  // Read in constants from first half of data set
  ifstream firstFile( firstHalfPath );
  
  if( !firstFile ){
    cerr << "No firstHalf.txt exists.." << endl;
    return 0;
  }
  
  int n1=0;
  double d1;
  int rval, cval;
  string line;
  
  while( getline( firstFile, line ) ){
    if( line.at(0) == '#' ) {
      continue;
    }
    istringstream iss( line );
    
    while( iss >> d1 ){

      gfirstHalf[n1] = d1; 
    
      cout << "firstHalf " << n1 << " = " << d1 << "." << endl;
      
      n1++;
    }
  }
  
  // Read in alpha parameters for PMTs
  ifstream secondFile( secondHalfPath );
  if( !secondFile ){
    cerr << "No secondHalf.txt exists.." << endl;
    return 0;
  }
  
  n1=0;
  d1=0;
  string line2;
  
  while( getline( secondFile, line2 ) ){
    if( line2.at( 0 )=='#' ) {
      continue;
    }
    istringstream iss( line2 );
    
    while( iss >> d1 ){
      
      gsecondHalf[n1] = d1;
      
      cout << "secondHalf " << n1 << " = " << d1 << "." << endl;
      
      n1++;
    }
  }

  // Print to console

  cout << "First half minus second half.." << endl;
  for( int i = 0; i<kNcell; i++){
    
    if(gfirstHalf[i]==1||gsecondHalf[i]==1) {
      cout << -1 << endl;
    }else{
      cout << abs(gfirstHalf[i]-gsecondHalf[i]) << endl;
    }
  }
  cout << endl << endl;
    
  cout << "First half divided by second half.." << endl;
  // Print to console
  for( int i = 0; i<kNcell; i++){
    if(gfirstHalf[i]==1||gsecondHalf[i]==1) {
      cout << -1 << endl;
    }else{
      cout << abs(gfirstHalf[i]/gsecondHalf[i]) << endl;
    }
  }
  
  cout << endl << endl;

  // Declare outfiles
  ofstream HCal_Error;

  HCal_Error.open( errorOutPath );

  for( int i=0; i<kNcell; i++ ){   
    if(gfirstHalf[i]==1||gsecondHalf[i]==1) {
      HCal_Error << -1 << endl;
    }else{
      HCal_Error << abs(gfirstHalf[i]/gsecondHalf[i]) << endl;
    }    
  }



  HCal_Error.close();

  st->Stop();

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

  return 0;

}


