//sseeds 5.2.23 - Short script to correct db params with vertical and horizontal corrections adapted to work without sbs-ana framework
//NOTE: Requires properly configured $DB_DIR environment variable to function

#include <vector>
#include <iostream>
#include <iomanip>
#include "TCut.h"
#include "TLorentzVector.h"
#include "TTreeFormula.h"
#include "TChain.h"
#include "TLatex.h"
#include "TH1D.h"
#include "TStopwatch.h"

const Int_t kNcell = 288;                // Total number of HCal channels
const Int_t kNrows = 24;                 // Total number of HCal rows
const Int_t kNcols = 12;                 // Total number of HCal columns
const Double_t hcalblk_div_h = 0.15494;  //m, horizontal center-to-center dist.
const Double_t hcalblk_div_v = 0.15875;  //m, vertical center-to-center dist.
const Double_t hcal_hrange = 1.85928;    //m, total range in horizontal direction of HCal (end-to-end)
const Double_t hcal_vrange = 3.81;       //m, total range in vertical direction of HCal (end-to-end)
const Double_t second = 1000000.;

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

//v_corr should be positive for towards the sky, h_offset should be positive for away from beamline
//v_corr for vertical correction to current parameter (x, dispersive, neutron peak location for SBS)
//h_corr for horizontal correction to current parameter (y, non-dispersive) 
void HCal_corrDBoffsets( Double_t v_corr = -1000., Double_t h_corr = -1000.){ //main

  //Require that user pass offsets
  if( v_corr==-1000. || h_corr==-1000. ){
    std::cout << "Error: User must pass two args <vertical offset (+,sky)> <horizontal offset (+,away-from-beam)>" << std::endl;
    return;
  }

  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  // Get the date
  std::string date = getDate();

  // Create files to store corrected hcal positions
  ofstream hcal_db_corrected;
  std::string corrPath = "hcal_db_corrected.txt";

  hcal_db_corrected.open( corrPath );
  
  // Create arrays to hold old positions. Initialize to zero.
  Double_t oldXpos[kNcell] = {0.}; //dispersive
  Double_t oldYpos[kNcell] = {0.}; //transverse

  // Read in previous offsets
  TString DBpath = gSystem->Getenv("DB_DIR");

  // Keep seperate for testing
  TString XConstPath = DBpath + "/db_sbs.hcal.dat";
  TString YConstPath = DBpath + "/db_sbs.hcal.dat";

  cout << "Loading current database position offsets from: " << XConstPath << ".." << endl;
  ifstream XConstFile( XConstPath );
  if( !XConstFile ){
    cerr << endl << "ERROR: No input constant file present -> path to db_sbs.hcal.dat expected." << endl;
    return 0;
  }

  Int_t n1=0;
  Double_t d1=0;
  std::string line;
  bool skip_line = true;

  // Get X position offsets
  while( getline( XConstFile, line ) ){
      
    if( n1==( kNcell ) ) break;
      
    TString Tline = (TString)line;
      
    if( Tline.BeginsWith("sbs.hcal.xpos") && skip_line==true ){	
      skip_line = false;
      continue;
    }

    if( skip_line==false ){
      istringstream iss( line );
      while( iss >> d1 ){
	oldXpos[n1] = d1;
	n1++;
      }
    }
  }

  cout << "Loading current database position offsets from: " << YConstPath << ".." << endl;
  ifstream YConstFile( YConstPath );
  if( !YConstFile ){
    cerr << endl << "ERROR: No input constant file present -> path to db_sbs.hcal.dat expected." << endl;
    return 0;
  }

  n1=0;
  d1=0;
  skip_line = true;
  std::string line2;

  // Get Y position offsets
  while( getline( YConstFile, line2 ) ){
      
    if( n1==( kNcell ) ) break;
      
    TString Tline = (TString)line2;
      
    if( Tline.BeginsWith("sbs.hcal.ypos") && skip_line==true ){	
      skip_line = false;
      continue;
    }

    if( skip_line==false ){
      istringstream iss( line2 );
      while( iss >> d1 ){
	oldYpos[n1] = d1;
	n1++;
      }
    }
  }
    
  // Console outs to verify offsets
  std::cout << std::setprecision(7); //std::cout default precision is 5 sig figs, set to 7
  std::cout << endl << endl << "Current X offsets from database: " << std::endl;
    
  for( Int_t r=0; r<kNrows; r++){
    for( Int_t c=0; c<kNcols; c++){
      Int_t i = r*kNcols+c;
      std::cout << oldXpos[i] << " ";
    }
    std::cout << endl;
  }
  
  usleep( 3*second ); //Give some time for review of step
  
  std::cout << std::endl << std::endl << "Current Y offsets from database: " << std::endl;
    
  for( Int_t r=0; r<kNrows; r++){
    for( Int_t c=0; c<kNcols; c++){
      Int_t i = r*kNcols+c;
      std::cout << oldYpos[i] << " ";
    }
    std::cout << endl;
  }
  
  usleep( 3*second ); //Give some time for review of step

  // Declare arrays for storing corrected positions
  Double_t corrXpos[kNcell] = {0.};
  Double_t corrYpos[kNcell] = {0.};

  // Calculate new offsets and console out for verification
  std::cout << std::endl << std::endl << "Corrected X offsets from database: " << std::endl;
  hcal_db_corrected << std::setprecision(7); //std::cout default precision is 5 sig figs, set to 7
  hcal_db_corrected << "sbs.hcal.xpos =" << std::endl;

  for( Int_t r=0; r<kNrows; r++){
    for( Int_t c=0; c<kNcols; c++){
      Int_t i = r*kNcols+c;
      corrXpos[i] = oldXpos[i]-v_corr;
      std::cout << corrXpos[i] << " ";
      hcal_db_corrected << corrXpos[i] << " "; 
    }
    std::cout << std::endl;
    hcal_db_corrected << std::endl;
  }
  
  usleep( 3*second ); //Give some time for review of step
  
  std::cout << std::endl << std::endl << "Corrected Y offsets from database: " << std::endl;
  hcal_db_corrected << std::endl << std::endl << "sbs.hcal.ypos =" << std::endl;

  for( Int_t r=0; r<kNrows; r++){
    for( Int_t c=0; c<kNcols; c++){
      Int_t i = r*kNcols+c;
      corrYpos[i] = oldYpos[i]-h_corr;
      std::cout << corrYpos[i] << " ";
      hcal_db_corrected << corrYpos[i] << " "; 
    }
    std::cout << std::endl;
  }
  
  usleep( 3*second ); //Give some time for review of step

  
  hcal_db_corrected.close();

  std::cout << std::endl << "Done. CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << std::endl;

  std::cout << "Database parameters generated and stored in working directory: " << corrPath << std::endl;

}//end main
