//SSeeds 10.12.21 - Production - Code combining elements from Andrew Puckett's BigCal calibration script and PDatta/EFuchey BB calibration code used to calibrate HCal energy deposition by PMT module. Tested with Eric Fuchey's digitized g4sbs data and used during detector commissioning in GMn run group.

#include <ctime>
#include <iomanip>
#include <iostream>
#include <fstream>
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
#include "gmn_rec_tree.C"

const int kNcell = 288; // Total number of HCal modules
const int kNrows = 24; // Total number of HCal rows
const int kNcols = 12; // Total number of HCal columns
const double kdBlock = 0.152; // Width and height of each module including distance between

const double PI = TMath::Pi();
const double M_e = 0.00051;
const double M_p = 0.938272;
const double M_n = 0.939565;
const double c_light = 299792458.0;

void hcal_eCal_p( const char *configfilename, int run = -1 ){
  
  // Start the chain for root files passed with config file
  TChain *C = new TChain("T");
  
  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );
  
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
  
  // Reading config file
  ifstream configfile(configfilename);
  TString currentline;
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      C->Add(currentline);
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
      if( skey == "HVfilePath" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	HVfilePath = sval;
      }
      if( skey == "alphasFilePath" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	alphasFilePath = sval;
      }
      if( skey == "constPath" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	constPath = sval;
      }
      if( skey == "HVoutPath" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	HVoutPath = sval;
      }
    }
    delete tokens;
  }
  
  //Declare outfile
  TFile *fout = new TFile( "hcalECalResults.root", "RECREATE" );
  
  //initialize vectors and arrays
  vector<int> hitlist;
  vector<double> energy;
  double gHV[kNrows][kNcols];
  double gAlphas[kNrows][kNcols];
  double gTargetHV[kNrows][kNcols];
  double gCRatio[kNrows][kNcols];
  
  //initialize histograms
  TH1D *h_W = new TH1D( "W", "W", 200, 0, 2 );
  h_W->GetXaxis()->SetTitle( "GeV" );
  TH1D *h_Q2 = new TH1D( "Q2", "Q2", 200, 0, 2 );
  h_Q2->GetXaxis()->SetTitle( "GeV" );
  TH1D *h_E_ep = new TH1D( "Scattered Electron Energy","E_ep",500,-2,8);
  h_E_ep->GetXaxis()->SetTitle( "GeV" );
  TH1D *h_E_pp = new TH1D("Scattered Proton Energy","E_pp",500,-2,15);
  h_E_pp->GetXaxis()->SetTitle( "GeV" );
  TH1F *h_HCalCol_p = new TH1F( "HCalCol_p", "Projection of Scattered p on HCal Columns", kNcols, 0, kNcols );
  h_HCalCol_p->GetXaxis()->SetTitle( "column" );
  TH1F *CvCh = new TH1F( "CvCh", "Calibration Constant vs Channel", kNcols*kNrows, 0, kNcols*kNrows );
  TH1F *NEVvCh = new TH1F( "NEVvChannel", "Number Events vs Channel", kNcols*kNrows, 0, kNcols*kNrows );
  TH1F *badCvCh = new TH1F( "badCvChannel", "Bad Cell vs Channel", kNcols*kNrows, 0, kNcols*kNrows );
  TH1F *CRatiovCh = new TH1F( "CRatiovChannel", "Constant Ratio (new/old) vs Channel", kNcols*kNrows, 0, kNcols*kNrows );

  
  // Build ADC histograms - Set integrated ADC spectrum histogram limits
  int SpecBins = 520;
  double IntSpec_min = -20.0, IntSpec_max = 500.0; 
  // Build histograms
  TH1F *gIntSpec[kNrows][kNcols];
  for( int r=0; r<kNrows; r++ ){
    for( int c=0; c<kNcols; c++ ){  
      int m = r*12+c;
      gIntSpec[r][c] = new TH1F( Form( "Int ADC Spect R%d C%d PMT%d", r, c, m ), Form( "Int ADC Spect R%d C%d PMT%d", r, c, m ), SpecBins, IntSpec_min, IntSpec_max );
      gIntSpec[r][c]->GetXaxis()->SetTitle( "pC" );
    }
  }

  Long64_t Nevents = C->GetEntries();
  cout << "Opened up tree with nentries: " << Nevents << ".." << endl;

  // Read in previous high voltage values
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
    
    n1++;
  }
  
  ifstream alphaFile( alphasFilePath );
  
  if( !alphaFile ){
    cerr << "No PMT alphas file present -> setFiles/alphas.txt expected." << endl;
    return 0;
  }

  cout << "Getting alpha parameters for each pmt.." << endl;
  
  n1=0;
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
    
    n1++;
  }
  
  cout << "HV settings and alpha parameters loaded." << endl;

  // Need to get electron energy and W for cuts on elastic events
  long mevent = 0;
  
  cout << "Preliminary loop over all events in tree to get physics quantities commencing.." << endl;
    
  gmn_rec_tree *T = new gmn_rec_tree(C);

  while( T->GetEntry( mevent++ ) ){ 
    
    if( mevent%1000 == 0 ) cout << "Preliminary loop: " << mevent << "/" << Nevents << endl;
    
    // Sort all tracks by lowest Chi^2 (highest confidence track)
    int track_tot = (int)T->bb_tr_n;
    int track = 0;
    double min = 1000.0; //Set arbitrarily high chi^2 to minimize with loop over tracks

    for( int elem=0; elem<track_tot; elem++ ){            
      if( T->bb_tr_chi2[elem] < min )
	track = elem;      
    }
        
    TVector3 lVep( T->bb_tr_px[track], T->bb_tr_py[track], T->bb_tr_pz[track] ); // Construct the scattered electron momentum vector from the tree (GEM reconstructed)
    
    double E_ep = sqrt( pow(M_e,2) + lVep.Mag2() ); // Obtain the scattered electron energy
    h_E_ep->Fill( E_ep ); // Fill histogram

    double p_ep = lVep.Mag(); // Obtain the magnitude of scattered electron momentum

    double theta_e = acos( ( T->bb_tr_pz[track])/p_ep ) * 180.0 / PI; // Obtain the electron's scattering angle in degrees
    
    double phi_e = TMath::ASin( T->bb_tr_py[track]/( p_ep*TMath::Sin( theta_e ) ) );
    
    double Q2 = 2*E_e*E_ep*( 1-(lVep.Z()/p_ep) ); // Obtain Q2 from beam energy, outgoing electron energy, and momenta

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

    double p_p = sqrt( pow( E_pp,2 )+W2 );

    //Each component in the form p_p_z = p_p*cos(phi)
    double phi_p = TMath::ACos( ( E_e-T->bb_tr_pz[track])/p_p ) * 180.0 / PI;

    //Get the projection of each event onto the face of HCal - consider only X (column) components to avoid calculations with 48D48 field map. Should still be useful
    double relPosX = HCal_d*TMath::Tan( HCal_th-phi_p );
    int relPosCol = std::floor( relPosX*HCal_th+6 ); // Return the column where the p is expected to land in HCal
    h_HCalCol_p->Fill( relPosCol );

    if( W>0&&W<2 ){
      //if( W>( 0.945-0.02 )&&W<( 0.945+0.02 ) ){
      hitlist.push_back( mevent );
      energy.push_back( E_pp );
    }
    
  }

  fout->cd();

  h_W->Write( "W" );
  h_W->Draw( "AP" );

  h_Q2->Write( "Q2" );
  h_Q2->Draw( "AP" );

  h_E_ep->Write( "E_ep" );
  h_E_ep->Draw( "AP" );

  h_E_pp->Write( "E_pp" );
  h_E_pp->Draw( "AP" );

  h_Q2->Write( "Q2" );
  h_Q2->Draw( "AP" );

  h_HCalCol_p->Write( "HCalCol_p" );
  h_HCalCol_p->Draw( "AP" );

  cout << "Cut and electron energy histograms written to file." << endl;
  cout << "Number of events passing global cut = " << hitlist.size() << endl;
  
  // Initialize matrix M and vector b to set up the linear system:
  // kNcell here is the total number of channels to be calibrated:
  TMatrixD M( kNcell,kNcell );
  TVectorD b( kNcell );
  
  // Array to count the number of good events in each cell:
  int nevents[kNcell];
  
  for( int i=0; i<kNcell; i++ ){
    for( int j=0; j<kNcell; j++ ){
      M(i,j) = 0.0; // Initialize sums for matrix M to zero.
    } 
    b(i) = 0.0; // Initialize sums for vector b to zero
    nevents[i] = 0; // Initialize event counters to zero.
  }
  
  long nevent=0;
  
  //Initialize old constants (for updating) to zero
  double oldconstants_nevent[kNcell];
  double oldconstants_sum[kNcell];
  double oldconstants_sum2[kNcell];
  for( int i=0; i<kNcell; i++ ){
    oldconstants_sum[i] = 0.0;
    oldconstants_nevent[i] = 0.0;
    oldconstants_sum2[i] = 0.0;
  }
 
  //////////////////////////////
  //Main loop over the events://
  //////////////////////////////
  
  //Loop over only elastic events as passed by BB
  //Assume that a clean sample of elastic events has already been selected by global_cut:
  //Minimize chi^2 = sum_i=1,nevent (E_i - sum_j c_j A^i_j)^2/sig_E^2, where sig_E_i^2 ~ E_i; in this case, dchi^2/dc_k = 2 * sum_i=1,nevent 1/E_i * (E_i - sum_j c_j A_j) * A_k, set the difference on the RHS to zero to minimize.
  // This is a system of linear equations for c_k, with RHS = sum_i A_k and LHS = sum_i sum_j A_j c_j A_k/E_i 
  
  cout << "Beginning main loop over events.." << endl;

  while( T->GetEntry( hitlist[nevent++] ) ){ 
    
    if( hitlist[nevent]>C->GetEntries() ) continue;

    if( nevent%1000 == 0 ) cout << "Main loop: " << nevent << "/" << hitlist.size() << endl;
    
    int r,c;

    double adc_p[kNrows][kNcols] = {0.0};
    double e[kNrows][kNcols] = {0.0};
    double amp[kNrows][kNcols] = {0.0};

    int saturated = 0;

    for( int m = 0; m < T->Ndata_sbs_hcal_adcrow; m++ ) { //Process event with m data 
      r = T->sbs_hcal_adcrow[m];
      c = T->sbs_hcal_adccol[m];
      if( r < 0 || c < 0 ) {
	cerr << "Error: row or col negative." << endl;
	continue;
      }
      
      if( r >= kNrows || c >= kNcols ) continue;
      
      // Fill adc and energy arrays from tree
      adc_p[r][c] = T->sbs_hcal_a_p[m]; 
      amp[r][c] = T->sbs_hcal_a_amp[m];
      e[r][c] = T->sbs_hcal_a_c[m];
      
      // Mark saturated array when amplitude meets max RAU
      if( amp[r][c] > 4095 ) {
	saturated = 1;
      }
    }
    
    int best = T->sbs_hcal_clus_id[0]; // Index in the array of the "best" block in the best cluster found in this event. Currently only one cluster is recorded per event, so it is the best

    int ncellclust = T->sbs_hcal_clus_nblk[0]; // Total number of hits in the "best" cluster, default to 9 (3x3) for now
    
    if( best >= 0 && ncellclust >= 1 && saturated == 0 ){ // Require at least 2 hits in this cluster to use as a "good" event for calibration purposes:

      double E_p = energy[nevent]; // Energy of e_p from first loop for this event

      // Check whether this cluster has its maximum at the edge of the calorimeter
      int rowmax = best/kNcols;
      int colmax = best%kNcols;
      bool edge_max = false;

      if( rowmax <= 0 || rowmax >= kNrows || colmax <= 0 || colmax >= kNcols ){
	edge_max = true;
      }
     
      if( !edge_max ){  //Only consider clusters with maximum at least one block away from the edge for calibration:
	for(int ihit=0; ihit<ncellclust; ihit++ ){ 

	  int cell_i = T->sbs_hcal_clus_blk_id[ihit]-1;

	  if(cell_i<0) cout << "cell" << cell_i << " less than zero." << endl;

	  int row_i = cell_i/kNcols;
	  int col_i = cell_i%kNcols;

	  float Ahit_i = adc_p[row_i][col_i]; // Ahit_i is the ADC value for this hit (ped-subtracted)
	  float Ehit_i = e[row_i][col_i]; // Ehit_i is the reconstructed energy for this hit (using previous calibration constants where a_c = const*a_p)

	  gIntSpec[row_i][col_i]->Fill( Ahit_i ); // Opportune time to fill ADC histogram for this run (diagnostic)
	  
	  nevents[cell_i] += 1; // Increment recorded event counter for cell i
	  if(nevents[cell_i] > maxEventPerCell) continue; // Added to keep statistics similar across all channels

	  //Check old calibration constants:
	  oldconstants_sum[cell_i] += Ehit_i/Ahit_i;
	  oldconstants_nevent[cell_i] += 1.0;
	  oldconstants_sum2[cell_i] += pow(Ehit_i/Ahit_i,2);

	  for(int jhit=0; jhit<ncellclust; jhit++ ){ 
	  
	    int cell_j = T->sbs_hcal_clus_blk_id[jhit]-1;

	    int row_j = cell_j/kNcols;
	    int col_j = cell_j%kNcols;

	    float Ahit_j = adc_p[row_j][col_j]; //ADC value of hit j 
	    
	    //increment Matrix element M_{i,j} with Ai * Aj / E_p, where, recall A_i and A_j are ADC values of hit i and hit j, and E_p is the predicted energy of the elastically scattered proton

	    if( cell_i >= kNcell || cell_j >= kNcell ) {
	      cout << "Error: Cell number greater than expected geometry." << endl;
	      continue;
	    }
	    
	    M( cell_i,cell_j ) += Ahit_i * Ahit_j / E_p;
	    
	  }
	  b(cell_i) += Ahit_i; //increment vector element i with ADC value of hit i. 

	} 
      }
    }
  }

  cout << "Matrix populated.." << endl;

  // IF the number of events in a cell exceeds the minimum, then we can calibrate 

  // Make adjustments to the matrix M and vector b to exclude "bad" cells from calibration: 
  int badcells[kNcell];
  double xErr[kNcell]; // Set up error for TGraphErrors

  for(int i=0; i<kNcell; i++){
    badcells[i] = 0;

    if( nevents[i] < minEventPerCell || M(i,i) < 0.1*b(i) ){ //If we have fewer than the minimum (100) events to calibrate, OR the diagonal element of matrix M for this cell is less than 0.1 * the corresponding vector element, exclude this block from the calibration. To do so without affecting the solution of the system, we do the following: 
      if(M(i,i)<0.1*b(i)) cout << "M(i,i)<0.1*b(i)" << endl;
      b(i) = 1.0;  // set RHS vector for cell i to 1.0 
      M(i,i) = 1.0; //set diagonal element of Matrix M for cell i to 1.0 
      for(int j=0; j<kNcell; j++){ //Set all off-diagonal elements for Matrix M for cell i to 0:
	if( j != i ){
	  M(i,j) = 0.0; 
	  M(j,i) = 0.0;
	}
      }
      badcells[i] = 1;
    }
  }

  cout << "Solving for calibration constants..." << endl;

  // Invert the matrix and multiply M^{-1} by the "RHS" vector b: 
  TVectorD solution = M.Invert() * b;

  cout << "Matrix inverted.." << endl;
  cout << endl;
  cout << "..done." << endl;
  
  // Declare outfiles
  ofstream HCalECal_Const;
  ofstream HCalECal_HV;

  HCalECal_Const.open( constPath );
  HCalECal_HV.open( HVoutPath );
  HCalECal_HV << "Module" << '\t' << " HV" << endl;

  // Declare output calibration histograms
  TH1D *Constants_chan = new TH1D("Constants_chan","",kNcell,0.5,kNcell+0.5);
  TH1D *Constants = new TH1D("Constants","",200,0.0,2.0);
  TH2D *Constants_rowcol = new TH2D("Constants_rowcol","",kNcols,0.5,kNcols+0.5,kNrows,0.5,kNrows+0.5);

  TH1D *Old_chan = new TH1D( "Old_chan","",kNcell, 0.5, kNcell+0.5 );
  TH1D *Old_chan_rms = new TH1D( "Old_chan_rms","",kNcell,0.5,kNcell+0.5 );
  TH1D *Ratio = new TH1D( "Ratio","new/old",200,0.5, 1.5 );
  TH1D *Ratio_chan = new TH1D( "Ratio_chan","new/old by chan.",kNcell,0.5, kNcell+0.5 );
  TH2D *Ratio_rowcol = new TH2D( "Ratio_rowcol","new/old by row/col",kNcols,0.5,kNcols+0.5,kNrows,0.5,kNrows+0.5 );

  // Declare calibration variables
  float ncalibrated = 0.0;
  float average_constant = 0.0;
  double y[kNcell] = {0.0};

  // Loop over all the cells, fill diagnostic histograms, and write out new constants:

  double E_targ = 571; //Average peak for 14 MeV is 8 pC, so our conversion is 8/0.014=571 pC/GeV to get GeV. Similarly 14 MeV is 18 mV for maxADC (for reference). This is the target gain of the PMT

  for( int i=0; i<kNcell; i++ ){

    y[i] = i; // Set up a vector for easy TGraphErrors construction later
    
    xErr[i] = E_targ*sqrt(fabs(M(i,i))); // Overestimate the error

    NEVvCh->SetBinContent(i+1,nevents[i]); // Fill the number of events histogram

    int r = i/kNcols;
    int c = i%kNcols;
    
    if( nevents[i] >= minEventPerCell ){
      
      
      double oldconst_avg = oldconstants_sum[i]/oldconstants_nevent[i];
      double oldconst_rms = sqrt(fabs( oldconstants_sum2[i]/oldconstants_nevent[i] - pow(oldconst_avg,2) ) );
	
      Old_chan->Fill(i+1, oldconst_avg*E_targ );
      Old_chan_rms->Fill(i+1, oldconst_rms*E_targ);
	
      average_constant += solution(i);
      ncalibrated += 1.0;
	
      Constants_chan->Fill( i+1, solution(i)*E_targ );
      Constants_chan->SetBinError( i+1, E_targ*sqrt(fabs(M(i,i))) );	

      Constants->Fill( solution(i)*E_targ );
      
      gCRatio[r][c] = solution[i]/oldconst_avg; // Ratio of new to old constants
      Ratio->Fill( gCRatio[r][c] );
      Ratio_chan->Fill( i+1, gCRatio[r][c] );
      CRatiovCh->SetBinContent( kNcols*r+c+1, gCRatio[r][c] );

      int irow, icol;
	
      if( (i+1) <= kNcols*kNrows ){
	irow = i/kNcols;
	icol = i%kNcols;
      }
	
      Constants_rowcol->Fill( icol+1, irow+1, E_targ*solution(i) );
      Constants_rowcol->SetBinError( icol+1, irow+1, E_targ*sqrt( fabs( M(i,i) ) ) );
      Ratio_rowcol->Fill( icol+1, irow+1, solution(i)/oldconst_avg );
    }
  }

  average_constant /= ncalibrated;
  
  HCalECal_Const << "HCal_cfac = " << endl;

  for( int i=0; i<kNcell; i++ ){
    char cconst[20];
    
    if( badcells[i] == 0 ){ //if the cell was calibrated and the calibration result is between 0.1 and 10, write out the new constant to the file:
      sprintf( cconst, "%15.4g", solution(i) );

      CvCh->SetBinContent( i+1, E_targ*solution(i) );
      
    } else { //Otherwise, use the average calibration constant:
      if( i < kNrows*kNcols ){
	//sprintf( cconst, "%15.4g", average_constant );
	sprintf( cconst, "%15.4g", oldconstants_sum[i]/oldconstants_nevent[i]); //Use old constant instead of average
	//sprintf( cconst, "%15.4g", -10 ); //Make these failures very obvious
      } else {
	cout << "Error: looping outside of normal geometry for cconst" << endl;
      }
    }
    
    badCvCh->SetBinContent(i+1,badcells[i]);
    
    int r = i/kNcols;
    int c = i%kNcols;

    cout << "Constant r" << r << " c" << c << " = " << cconst << endl;

    HCalECal_Const << cconst;
    if( i+1 != kNrows*kNcols && i+1 != kNcell ) HCalECal_Const << ",";
    if( (i+1)%kNcols == 0 ) HCalECal_Const << endl;
    if( i == kNrows*kNcols ) HCalECal_Const << "Error: kNcell > kNrows * kNcols" << endl;
  }

  HCalECal_HV << endl << "HCal_gain_cor = " << endl;

  for(int i=0; i<kNcell; i++){

    int r = i/kNcols;
    int c = i%kNcols;

    gTargetHV[r][c] = gHV[r][c]*pow(gCRatio[r][c],(1/gAlphas[r][c]));

    HCalECal_HV << r*12+c+1 << '\t' << -gTargetHV[r][c] << endl;

    char cconst[20];
    sprintf( cconst, "%15.4f", gTargetHV[r][c] );

    HCalECal_HV << cconst;
    if( i+1 != kNrows*kNcols && i+1 != kNcell ) HCalECal_HV << ",";
    if( (i+1)%kNcols == 0 ) HCalECal_HV << endl;
    if( i == kNrows*kNcols ) HCalECal_HV << "Error: kNcell > kNrows * kNcols" << endl;
  }

  double yErr[kNcell] = {0.0};
  TGraphErrors *ccgraph = new TGraphErrors( kNcell, y, &solution[0], yErr, xErr );
  ccgraph->GetXaxis()->SetLimits(0.0,288.0);

  //Write out root file for diagnostics and close
  fout->cd();
  CvCh->Write( "CvChannel" );
  CvCh->Draw( "AP" );

  ccgraph->SetTitle("Calibration Constants");
  ccgraph->GetXaxis()->SetTitle("Channel");
  ccgraph->GetXaxis()->SetTitle("GeV/pC");
  ccgraph->SetMarkerStyle(21); //Boxes
  ccgraph->Write("CONSTANTS");
  ccgraph->Draw("AP");

  NEVvCh->Write( "NEVvChannel" );
  NEVvCh->Draw( "AP" );
  
  badCvCh->Write( "badCvChannel" );
  badCvCh->Draw( "AP" );
  
  CRatiovCh->Write( "CRatiovChannel" );
  CRatiovCh->Draw( "AP" );

  Constants_chan->Write( "Constants_chan" );
  Constants_chan->Draw( "AP" );

  Constants->Write( "Constants" );
  Constants->Draw( "AP" );

  Constants_rowcol->Write( "Constants_rowcol" );
  Constants_rowcol->Draw( "AP" );

  Old_chan->Write( "Old_chan" );
  Old_chan->Draw( "AP" );

  Old_chan_rms->Write( "Old_chan_rms" );
  Old_chan_rms->Draw( "AP" );
  
  Ratio->Write( "Ratio" );
  Ratio->Draw( "AP" );
  
  Ratio_chan->Write( "Ratio_chan" );
  Ratio_chan->Draw( "AP" );

  Ratio_rowcol->Write( "Ratio_rowcol" );
  Ratio_rowcol->Draw( "AP" );

  for( int r=0; r<kNrows; r++ ){
    for( int c=0; c<kNcols; c++ ){
      gIntSpec[r][c]->GetXaxis()->SetRange(-20,150); //For easier viewing
      gIntSpec[r][c]->Write(Form("intADC r%d c%d",r,c));
      gIntSpec[r][c]->Draw("AP");
    }
  }
  fout->Write();
  fout->Close();

  st->Stop();

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

}
  
