//Code to calibrate HCal energy deposition by PMT module - proof of concept with assumptions made by SSeeds pending tree structure update from BB tracking. SSeeds 9.1.21

#include "hcal.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TString.h"
#include "TEventList.h"
#include "TCut.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TClonesArray.h"
#include <TH2F.h>
#include <TChain.h>
#include <TCanvas.h>
#include <iostream>
#include <TSystem.h>
#include <TStopwatch.h>
#include <iomanip>

TChain *T = 0;

const int ncell = 288; 
const int kNrows = 24;
const int kNcols = 12;

void hcal_e_cal_POC( int run = -1 ){

  //User input .root file for sample analysis
  cout << "Enter run number for toy analysis." << endl;
  cin >> run;

  //Empirical limits
  double min_const=0.1, max_const=10.0;

  //Add generic detector specs
  int nrows = 24, ncols = 12;

  //Add output file
  TFile *fout = new TFile( Form("hcalECal_%d.root",run), "RECREATE" );
     
  //Hard code this with arbitrary cut for now, just to check the integrity of the algorithm
  TCut global_cut = ""; 

  //Create event list to loop over
  TEventList *elist = new TEventList("elist");

  // Read in data produced by analyzer in root format
  cout << "Reading replayed root file.." << endl;
  if( !T ) { 
    T = new TChain( "T" );
    T->Add( Form( "/volatile/halla/sbs/seeds/rootfiles/hcal_%d*.root", run ) );
    T->SetBranchStatus( "*", 0 );
    T->SetBranchStatus( "sbs.hcal.*", 1 );
    T->SetBranchAddress( "sbs.hcal.a", hcalt::a );
    T->SetBranchAddress( "sbs.hcal.a_amp", hcalt::a_amp );
    T->SetBranchAddress( "sbs.hcal.a_p", hcalt::a_p );
    T->SetBranchAddress( "sbs.hcal.a_amp_p", hcalt::a_amp_p );
    T->SetBranchAddress( "sbs.hcal.tdc", hcalt::tdc );
    T->SetBranchAddress( "sbs.hcal.ped", hcalt::ped );
    T->SetBranchAddress( "sbs.hcal.adcrow", hcalt::row );
    T->SetBranchAddress( "sbs.hcal.adccol", hcalt::col );
    T->SetBranchStatus( "Ndata.sbs.hcal.adcrow", 1 );
    T->SetBranchAddress( "Ndata.sbs.hcal.adcrow", &hcalt::ndata );

    // Add clustering branches
    T->SetBranchAddress("sbs.hcal.clus.id",hcalt::cid);
    T->SetBranchAddress("sbs.hcal.clus.row",hcalt::crow);
    T->SetBranchAddress("sbs.hcal.clus.col",hcalt::ccol);
    T->SetBranchAddress("sbs.hcal.clus.e",hcalt::ce);
    T->SetBranchAddress("sbs.hcal.clus.eblk",hcalt::ceblk);
    T->SetBranchAddress("sbs.hcal.clus.nblk",hcalt::cnblk);
    T->SetBranchAddress("sbs.hcal.nclus",hcalt::nclus);
    T->SetBranchAddress("sbs.hcal.nblk",hcalt::nblk);
    T->SetBranchAddress("sbs.hcal.clus_blk.id",hcalt::cblkid);
    T->SetBranchAddress("sbs.hcal.clus.eblk",hcalt::ceblk);
    
    cout << "Opened up tree with nentries=" << T->GetEntries() << endl;
  }

  //Parse the tree with events passing the cut
  T->Draw( ">>elist", global_cut );
  
  cout << "Number of events passing global cut = " << elist->GetN() << endl;

  //Initialize matrix M and vector b to set up the linear system:
  //ncell here is the total number of channels to be calibrated:
  TMatrixD M(ncell,ncell);
  TVectorD b(ncell);

  //Array to count the number of good events in each cell:
  int nevents[ncell];

  for(int i=0; i<ncell; i++){
    for(int j=0; j<ncell; j++){
      M(i,j) = 0.0; //Initialize sums for matrix M to zero.
    } 
    b(i) = 0.0; //Initialize sums for vector b to zero
    nevents[i] = 0; //Initialize event counters to zero.
  }

  int min_events_per_cell=100; //require at least 100 events to calculate a new calibration constant for the cell:

  long nevent=0;

  //Initialize old constants (for updating) to zero
  double oldconstants_nevent[ncell];
  double oldconstants_sum[ncell];
  double oldconstants_sum2[ncell];
  for(int i=0; i<ncell; i++){
    oldconstants_sum[i] = 0.0;
    oldconstants_nevent[i] = 0.0;
    oldconstants_sum2[i] = 0.0;
  }

  TH1D *hdpel_hms = new TH1D("hdpel_hms","",200,-0.04,0.04);
  TH1D *hdx = new TH1D("hdx","",200,-30.0,30.0);
  TH1D *hdy = new TH1D("hdy","",200,-60.0,60.0);
  TH1D *hemiss = new TH1D("hemiss","(eclust/e_hms-1)",200,-1.0,1.0);
  TH1D *htdiff = new TH1D("htdiff","tclust-t_hms",200,-30,30);
  
  TH1D *hncell = new TH1D("hncell","ncells",25,0.5,25.5);
  TH1D *hnx    = new TH1D("hnx","size in x",6,0.5,6.5);
  TH1D *hny    = new TH1D("hny","size in y",6,0.5,6.5);
  TH2D *hnxy   = new TH2D("hnxy","size y vs. size x",6,0.5,6.5,6,0.5,6.5);
  
  TH1D *hxmom = new TH1D("hxmom", "x moment", 200, -5.0,5.0 );
  TH1D *hymom = new TH1D("hymom", "y moment", 200, -5.0,5.0 );
  TH1D *hxdiff = new TH1D("hxdiff","x_{clust} - x_{cell}", 100, -5.0, 5.0 );
  TH1D *hydiff = new TH1D("hydiff","y_{clust} - y_{cell}", 100, -5.0, 5.0 );

  TH2D *hxdiff_ixcell_prot = new TH2D("hxdiff_ixcell_prot","",32,0.5,32.5,100,-5.,5.);
  TH2D *hydiff_iycell_prot = new TH2D("hydiff_iycell_prot","",32,0.5,32.5,100,-5.,5.);
  
  hxdiff->GetXaxis()->SetTitle("x_{clust}-x_{cell} (cm)");
  hydiff->GetXaxis()->SetTitle("y_{clust}-y_{cell} (cm)");

  //////////////////////////////
  //Main loop over the events://
  //////////////////////////////

  //Loop over only elastic events as passed by BB
  //Assume that a clean sample of elastic events has already been selected by global_cut:
  //In this case, we are trying to minimize chi^2 = sum_i=1,nevent (E_i - sum_j c_j A^i_j)^2/sig_E^2, where sig_E_i^2 ~ E_i; in this case, dchi^2/dc_k = 2 * sum_i=1,nevent 1/E_i * (E_i - sum_j c_j A_j) * A_k 
  // This is a system of linear equations for c_k, with RHS = sum_i A_k and LHS = sum_i sum_j A_j c_j A_k/E_i 

  while( T->GetEntry( elist->GetEntry(nevent++)) ){ 

    if( nevent%1000 == 0 ) cout << nevent << endl;

    int r,c;
    double adc_p[kNrows][kNcols] = {0.0};
    double e[kNrows][kNcols] = {0.0};
    double amp[kNrows][kNcols] = {0.0};

    int saturated = 0;

    // Process event with m data
    for( int m = 0; m < hcalt::ndata; m++ ) {
      // Define row and column
      r = hcalt::row[m];
      c = hcalt::col[m];
      if( r < 0 || c < 0 ) {
	cerr << "Error: row or col negative." << endl;
	continue;
      }
      
      if( r >= kNrows || c >= kNcols ) continue;
      
      // Fill adc and energy arrays
      adc_p[r][c] = hcalt::a_p[m];
      amp[r][c] = hcalt::a_amp[m];
      e[r][c] = hcalt::cblke[m];
      
      // Mark saturated array when amplitude meets max RAU
      if( amp[r][c] > 3900 ) {
	saturated = 1;
      }
    }
    
    int best = hcalt::cid[0] - 1; //index in the array of the "best" block in the best cluster found in this event. Currently only one cluster is recorded per event

    //cout << "best = " << best << endl;

    int ncellclust = hcalt::nblk[0]; //total number of hits in the "best" cluster, default to 9 (3x3) for now
    

    if( best >= 0 && ncellclust >= 1 && saturated == 0 ){ // require at least 2 hits in this cluster to use as a "good" event for calibration purposes:

      float E_e = 14; //Predicted energy of elastic electron from BB. Here using avg energy deposition from cosmics (in MeV) as a stand-in. Will need to look at energy ratio between electron and proton for current kinematic and project for more accurate results

      //Check whether this cluster has its maximum at the edge of the calorimeter - very unrealistic test assumptions made where "best" is max energy:
      int rowmax = best/kNcols;
      int colmax = best%kNcols;

      //cout << "rowmax " << rowmax << " colmax " << colmax << endl; 

      bool edge_max = false;

      //Corrected for HCal geometry
      if( rowmax <= 2 || rowmax >= 23 || colmax <= 2 || colmax >= 11 ){
	edge_max = true;
      }
     
      if( !edge_max ){  //Only consider clusters with maximum at least one block away from the edge for calibration:
	for(int ihit=0; ihit<ncellclust; ihit++ ){ //outer loop over all cells in the cluster: 
	  //int dr_i = ihit/3-1;
	  //int dc_i = ihit%3-1;

	  //int row_i = rowmax+dr_i;
	  //int col_i = colmax+dc_i;

	  int cell_i = hcalt::cblkid[ihit]-1;

	  int row_i = cell_i/kNcols;
	  int col_i = cell_i%kNcols;

	  //cout << "OUT blk row id " << row_i << " OUT blk col id " << col_i << endl;

	  float Ahit_i = adc_p[row_i][col_i]; //Ahit_i is the ADC value for this hit (ped-subtracted)
	  float Ehit_i = e[row_i][col_i]; //Ehit_i is the reconstructed energy for this hit (using previous calibration constants)

	  //updated for HCal
	  //int cell_i;
	  //cell_i = col_i + kNcols*row_i;
	  
	  //Check old calibration constants:
	  oldconstants_sum[cell_i] += Ehit_i/Ahit_i;
	  oldconstants_nevent[cell_i] += 1.0;
	  oldconstants_sum2[cell_i] += pow(Ehit_i/Ahit_i,2);
	  
	  for(int jhit=0; jhit<ncellclust; jhit++ ){ //inner loop over all cells in the cluster: 
	    
	    //int dr_j = jhit/3-1;
	    //int dc_j = jhit%3-1;

	    //int row_j = rowmax+dr_j;
	    //int col_j = colmax+dc_j;

	    int cell_j = hcalt::cblkid[jhit]-1;

	    //cout << "cell_i = " << cell_i << " cell_j = " << cell_j << endl; 

	    int row_j = cell_j/kNcols;
	    int col_j = cell_j%kNcols;

	    //cout << "IN blk row id " << row_j << " IN blk col id " << col_j << endl;

	    float Ahit_j = adc_p[row_j][col_j]; //ADC value of hit j 

	    //updated for HCal
	    //int cell_j;
	    //cell_j = col_j + kNcols*row_j;

	    //increment Matrix element M_{i,j} with Ai * Aj / E_e, where, recall A_i and A_j are ADC values of hit i and hit j, and E_e is the predicted energy
	    // of the elastically scattered electron:
	    M(cell_i,cell_j) += Ahit_i * Ahit_j / E_e;
	    
	    if( cell_i >= 288 || cell_j >= 288) {
	      cout << "Error: Cell number greater than expected geometry." << endl;
	      return 0;
	    }

	    //cout << "Matrix element i:j " << cell_i << ":" << cell_j << " = " << M(cell_i,cell_j) << endl;

	  }
	  b(cell_i) += Ahit_i; //increment vector element i with ADC value of hit i. 
	  nevents[cell_i] += 1; //increment event counter for cell i:
	} 
      }
    }
  }

  cout << "Matrix populated.." << endl;


  //IF the number of events in a cell exceeds the minimum, then we can calibrate 

  int smalld[ncell];

  for(int i=0; i<ncell; i++){
    smalld[i] = 0;
    if( nevents[i] == 0 || M(i,i) <= 0.1*b(i) ){ //this cell has less than the minimal number of events required; set diagonal elements to 1, off-diagonals to zero:
      b(i) = 1.0;
      M(i,i) = 1.0;
      for(int j=0; j<ncell; j++){
   	if( j != i ) {
   	  M(i,j) = 0.0;
   	  M(j,i) = 0.0;
   	}
      }
      
      if( nevents[i] > 0 ) smalld[i] = 1;
      
    }
  }
  
  //Make adjustments to the matrix M and vector b to exclude "bad" cells from calibration: 
  int badcells[ncell];
  
  for(int i=0; i<ncell; i++){
    badcells[i] = 0;
    if( nevents[i] < min_events_per_cell || M(i,i) < 0.1*b(i) ){ //If we have fewer than the minimum (100) events to calibrate, OR the diagonal element of matrix M for this cell is less than 0.1 * the corresponding vector element, exclude this block from the calibration:
      //To do so without affecting the solution of the system, we do the following: 
      b(i) = 1.0;  // set RHS vector for cell i to 1.0 
      M(i,i) = 1.0; //set diagonal element of Matrix M for cell i to 1.0 
      for(int j=0; j<ncell; j++){ //Set all off-diagonal elements for Matrix M for cell i to 0:
	if( j != i ){
	  M(i,j) = 0.0; 
	  M(j,i) = 0.0;
	}
      }
      badcells[i] = 1;
    }
  }

  cout << "Solving for calibration constants..." << endl;

  //Invert the matrix and multiply M^{-1} by the "RHS" vector b: 
  TVectorD solution = M.Invert() * b;

  cout << "done." << endl;
  
  elist->Delete();

  ofstream constants_file;

  constants_file.open( Form( "outFiles/hcalECalParams_%d.txt", run ) );

  //No idea where these constants come from - SS
  
  TH1D *hconstants_chan = new TH1D("hconstants_chan","",ncell,0.5,ncell+0.5);
  TH1D *hconstants = new TH1D("hconstants","",200,0.0,2.0);
  TH2D *hconstants_rowcol = new TH2D("hconstants_rowcol","",32,0.5,32.5,56,0.5,56.5);
  
  float ncalibrated_prot = 0.0;
  float average_constant_prot = 0.0;

  TH1D *hold_chan = new TH1D("hold_chan","",ncell, 0.5, ncell+0.5 );
  TH1D *hold_chan_rms = new TH1D("hold_chan_rms","",ncell,0.5,ncell+0.5);
  TH1D *hratio = new TH1D("hratio","new/old",200,0.5, 1.5);
  TH1D *hratio_chan = new TH1D("hratio_chan","new/old by chan.",ncell,0.5, ncell+0.5);
  TH2D *hratio_rowcol = new TH2D("hratio_rowcol","new/old by row/col",32,0.5,32.5,56,0.5,56.5);

  //loop over all the cells, fill diagnostic histograms, and write out new constants:

  for(int i=0; i<ncell; i++){
    //    if( nevents[i] >= min_events_per_cell ){

    //If this cell was calibrated:
    if( badcells[i] == 0 && solution(i)*1000.0 >= min_const && solution(i)*1000.0 <= max_const ){

      //We multiply by 1000 here because the energy units were in GeV and the PMT HVs were typically set to give
      // 1,000 ADC = 1 GeV, so the constants were typically of order 10^-3. Here we multiply them by 1,000 to give a number of O(1).  
      //Different for HCal. We have max_adc = 4095 and max_E for GMn = 700 MeV. May need to update as comparitive orders of magnitude are 10 not 10^3
      
      double oldconst_avg = oldconstants_sum[i]/oldconstants_nevent[i];
      double oldconst_rms = sqrt(fabs( oldconstants_sum2[i]/oldconstants_nevent[i] - pow(oldconst_avg,2) ) );

      hold_chan->Fill(i+1, oldconst_avg*1000.0 );
      hold_chan_rms->Fill(i+1, oldconst_rms*1000.0);

      //compute average calibration constants for both sections of the calorimeter and use the average result for the channels which could not be calibrated:
      //Geometry dependent - SS
      //if( i < 32*32 ){
      if( i < ncols*nrows ){
      average_constant_prot += solution(i);
	ncalibrated_prot += 1.0;
      } else {
        cout << "Error: looping outside of HCal geometry" << endl;
      }
      
      hconstants_chan->Fill( i+1, solution(i)*1000.0 );
      //hconstants_chan->SetBinError( i+1, 1000.0*sqrt(fabs(M(i,i))) );
      
      hconstants->Fill( solution(i)*1000.0 );

      //compute the ratio of new / old calibration constants:
      hratio->Fill( solution(i)/oldconst_avg );
      hratio_chan->Fill( i+1, solution(i)/oldconst_avg );

      int irow, icol;
      
      if( (i+1) <= ncols*nrows ){
	irow = i/ncols;
	icol = i%ncols;
      }
      
      hconstants_rowcol->Fill( icol, irow, 1000.0*solution(i) );
      //hconstants_rowcol->SetBinError( icol, irow, 1000.0*sqrt(fabs(M(i,i))));
      hratio_rowcol->Fill( icol, irow, solution(i)/oldconst_avg );
    }
  }
  
  average_constant_prot /= ncalibrated_prot;
  
  constants_file << "HCal_prot_cfac = " << endl;
  
  for( int i=0; i<ncell; i++ ){
    char cconst[20];
    
    if( badcells[i] == 0 && solution(i)*1000.0 >= min_const && solution(i)*1000.0 <= max_const ){ //if the cell was calibrated and the calibration result is between 0.1 and 10, write out the new constant to the file:
      sprintf( cconst, "%15.4g", solution(i) );
    } else { //Otherwise, use the average calibration constant:
      if( i < nrows*ncols ){
	sprintf( cconst, "%15.4g", average_constant_prot );
      } else {
	cout << "Error: looping outside of normal geometry for cconst" << endl;
      }
    }

    cout << "Constant: " << cconst << endl;

    constants_file << cconst;
    if( i+1 != nrows*ncols && i+1 != ncell ) constants_file << ",";
    if( (i+1)%8 == 0 ) constants_file << endl;
    if( (i+1) == nrows*ncols ) constants_file << "Error: ncell > nrows * ncols" << endl;
  }
  
  constants_file << endl << "HCal_prot_gain_cor = " << endl;

  for(int i=0; i<ncell; i++){
    char cconst[20];
    sprintf( cconst, "%15.4f", 1.0 );

    cout << "Gain factor: " << cconst << endl;

    constants_file << cconst;
    if( i+1 != nrows*ncols && i+1 != ncell ) constants_file << ",";
    if( (i+1)%8 == 0 ) constants_file << endl;
    if( i+1 == nrows*ncols ) constants_file << "Error: ncell > nrows * ncols" << endl;
  }

  fout->Write();
  fout->Close();
}
