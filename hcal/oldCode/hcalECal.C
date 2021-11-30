//Code to calibrate HCal energy deposition by PMT module modified from BigCal script courtesy A.Puckett. SSeeds 8.20.21

//#include "bigcal_ntuple.C"
//#include "HCal_ntuple.C"
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

const int ncell = 1744; //Why? For bigCal 32 rows by 32 cols = 1024. Later, the difference (1744 - 1024 = 720) is attributed to RCS (?) where the bigCal cells are protvino (also ?)

void hcalECal( int run, const char *setupfilename, const char *outputfilename ){

  //set some limits on the calibration constants:
  //where do these limits come from? - SS
  double min_const=0.1, max_const=10.0;

  //Add generic detector specs
  int nrows = 24, ncols = 12;
 
  TString outfilename = outputfilename;
  if( !outfilename.EndsWith(".root") ){
    outfilename.ReplaceAll(".","_");
    outfilename += ".root";
  }

  TFile *fout = new TFile( outfilename.Data(), "RECREATE" );

  //Need to read in run number
  //Need to read in list of events from T tree which correspond to elastic scattering (TEventList)
  //Need to read in global cut parameters per event

  /*
  //TChain *C = new TChain("h9030");

  TString currentline;
  ifstream setupfile(setupfilename);

  if( !setupfile ){ 
    return; 
  }
  
  while( currentline.ReadLine(setupfile) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ) 
      C->Add(currentline.Data());
  }
  
  TCut global_cut = ""; 

  while( currentline.ReadLine(setupfile) && !currentline.BeginsWith("endcut") ){
    if( !currentline.BeginsWith("#") )
      global_cut += currentline.Data();
  }

  TEventList *elist = new TEventList("elist");

  C->Draw( ">>elist", global_cut );

  //populate the list of events passing the cut

  HCal_ntuple *T = new HCal_ntuple(C);
  */

  if( !T ){
    T = new TChain( "T" );
    T->Add( Form( "/volatile/halla/sbs/seeds/rootFiles/hcal_%d_*.root", run));
    T->SetBranchStatus( "*", 0 );
    T->SetBranchStatus( "sbs.hcal.*", 1 );
    T->SetBranchAddress( "sbs.hcal.nsamps", hcalt::nsamps );
    T->SetBranchAddress( "sbs.hcal.a", hcalt::a );
    T->SetBranchAddress( "sbs.hcal.tdc", hcalt::tdc );
    T->SetBranchAddress( "sbs.hcal.samps", hcalt::samps );
    T->SetBranchAddress( "sbs.hcal.samps_idx", hcalt::samps_idx );
    T->SetBranchAddress( "sbs.hcal.adcrow", hcalt::row );
    T->SetBranchAddress( "sbs.hcal.adccol", hcalt::col );
    //Add all clustering variables here
    //Need (at minimum): nclust, ncellclust, ibest, iycell, ixcell, ablock
    //
    //
    T->SetBranchStatus( "Ndata.sbs.hcal.adcrow", 1 );
    T->SetBranchAddress( "Ndata.sbs.hcal.adcrow", &hcalt::ndata );
    //Add proton energy from tree here
    //
    //
    if( T->GetEntries == 0 ){
      cout << "Run data is empty. Returning.." << endl;
      return;
    }
    cout << "Opened tree from run " << run << " with nentries = " << T->GetEntries() << endl;

  }

  //Need to read in E_p from BB

  TChain *C = new TChain("h9030");

  TEventList *elist = new TEventList("elist");

  TCut elistCut = ""; 

  TString currentline;
  ifstream setupfile(setupfilename);

  if( !setupfile ){ 
    cout << "Cannot locate setup file. Returning.." << endl;
    return; 
  }

  while( currentline.ReadLine(setupfile) && !currentline.BeginsWith("endcut") ){
    if( !currentline.BeginsWith("#") )
      elistCut += currentline.Data();
  }


  T->Draw( ">>elist", elistCut );

  cout << "Number of events passing global cut = " << elist->GetN() << endl;

  //initialize matrix M and vector b to set up the linear system:
  //ncell here is the total number of channels to be calibrated:
  TMatrixD M(ncell,ncell);
  TVectorD b(ncell);

  //array to count the number of good events in each cell:
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
  TH2D *hxdiff_ixcell_rcs = new TH2D("hxdiff_ixcell_rcs","",30,0.5,30.5,100,-5.,5.);
  TH2D *hydiff_iycell_rcs = new TH2D("hydiff_iycell_rcs","",24,32.5,56.5,100,-5.,5.);
  
  
  hxdiff->GetXaxis()->SetTitle("x_{clust}-x_{cell} (cm)");
  hydiff->GetXaxis()->SetTitle("y_{clust}-y_{cell} (cm)");

  //MAIN LOOP over the events:

  //Loop over only elastic events as passed by BB

  while( T->GetEntry( elist->GetEntry(nevent++)) ){ //Assume that a clean sample of elastic events has already been selected by global_cut:
    //In this case, we are trying to minimize 
    // chi^2 = sum_i=1,nevent (E_i - sum_j c_j A^i_j)^2/sig_E^2, where sig_E_i^2 ~ E_i; in this case, 
    // dchi^2/dc_k = 2 * sum_i=1,nevent 1/E_i * (E_i - sum_j c_j A_j) * A_k 
    // This is a system of linear equations for c_k, with RHS = sum_i A_k and LHS = sum_i sum_j A_j c_j A_k/E_i 
    if( nevent%1000 == 0 ) cout << nevent << endl;
    int best = int(T->ibest) - 1; //index in the array of the "best" cluster found in this event
    int ncellclust = int(T->ncellclust[best]); //total number of hits in the "best" cluster
    
    if( best >= 0 && T->nclust > 0 && ncellclust >= 1 ){ // require at least 2 hits in this cluster to use as a "good" event for calibration purposes:
      //******************************************************************
      //* Starting here, this information is irrelevant for calibration: *
      //******************************************************************
      float dpel_h = T->dpel_hms; 
      float dx = T->xclust[best]-T->X_hms;
      float dy = T->yclust[best]-T->Y_hms;
      float Emiss = T->eclust[best]/T->E_hms - 1.0;
      float tdiff = T->ctime_clust[best] - T->T_hms;

      hdpel_hms->Fill(dpel_h);
      hdx->Fill(dx);
      hdy->Fill(dy);
      hemiss->Fill(Emiss);
      htdiff->Fill(tdiff);

      hxdiff->Fill( T->xclust[best] - T->xcell[best][0] );
      hydiff->Fill( T->yclust[best] - T->ycell[best][0] );

      if( T->iycell[best][0] <= 32 ){
      
	hxdiff_ixcell_prot->Fill( T->ixcell[best][0], T->xclust[best] - T->xcell[best][0] );
	hydiff_iycell_prot->Fill( T->iycell[best][0], T->yclust[best] - T->ycell[best][0] );
      } else {
	hxdiff_ixcell_rcs->Fill( T->ixcell[best][0], T->xclust[best] - T->xcell[best][0] );
	hydiff_iycell_rcs->Fill( T->iycell[best][0], T->yclust[best] - T->ycell[best][0] );
      }
      
	
      hncell->Fill( T->ncellclust[best]);
      hnx->Fill( T->ncellx[best] );
      hny->Fill( T->ncelly[best] );
      hnxy->Fill( T->ncellx[best], T->ncelly[best] );
      
      hxmom->Fill( T->xmoment[best] );
      hymom->Fill( T->ymoment[best] );

      //******************************************************************
      //* End of irrelevant lines                                        *
      //******************************************************************
            
      float E_e = T->E_hms; //E_hms is the predicted electron energy in BigCal from the measured elastically scattered proton kinematics (assuming elastic ep scattering)

      //Check whether this cluster has its maximum at the edge of the calorimeter:
      int rowmax = T->iycell[best][0];
      int colmax = T->ixcell[best][0];
      bool edge_max = false;

      //Needs updating to HCal geometry (12 x 24 = col_tot x row_tot = 288 channels)
      //if( rowmax <= 2 || rowmax >= 55 || colmax <= 2 || (rowmax<=32&&colmax>=31) || (rowmax > 32 && colmax>=29) ){
      //edge_max = true;
      //}

      //Corrected for HCal geometry
      if( rowmax <= 2 || rowmax >= 23 || colmax <= 2 || colmax >= 11 ){
	edge_max = true;
      }
      
      // cout << "ibest,rowmax,colmax,ncell=" << best << ", " 
      // 	   << rowmax << ", " << colmax << ", " << ncell << endl;
     
      if( !edge_max ){  //Only consider clusters with maximum at least one block away from the edge for calibration:
	for(int ihit=0; ihit<ncellclust; ihit++ ){ //outer loop over all cells in the cluster: 
	  float Ahit_i = (T->ablock)[best][ihit]; //Ahit_i is the ADC value for this hit (ped-subtracted)
	  float Ehit_i = (T->eblock)[best][ihit]; //Ehit_i is the reconstructed energy for this hit (using previous calibration constants)

	  int row_i = (T->iycell)[best][ihit]; 
	  int col_i = (T->ixcell)[best][ihit];

	  /*
	  //the following lines map row and column indices to unique cell number for the BigCal geometry:
	  int cell_i;
	  if( row_i <= 32 ){
	    cell_i = col_i + 32*(row_i - 1) - 1;
	  } else {
	    cell_i = 32*32 + col_i + 30*(row_i - 33) - 1;
	  }
	  */

	  //updated for HCal
	  int cell_i;
	  cell_i = col_i + 12*row_i;
	  
	  //Check old calibration constants:
	  oldconstants_sum[cell_i] += Ehit_i/Ahit_i;
	  oldconstants_nevent[cell_i] += 1.0;
	  oldconstants_sum2[cell_i] += pow(Ehit_i/Ahit_i,2);

	  // cout << "i, row_i, col_i, cell_i=" << ihit << ", " << row_i 
	  //      << ", " << col_i << ", " 
	  //      << cell_i << endl;

	  
	  for(int jhit=0; jhit<ncellclust; jhit++ ){ //inner loop over all cells in the cluster: 
	    float Ahit_j = (T->ablock)[best][jhit]; //ADC value of hit j 
	    
	    int row_j = (T->iycell)[best][jhit];
	    int col_j = (T->ixcell)[best][jhit];

	    /*
	    //again, this is mapping row and column indices to a unique cell number:
	    int cell_j;
	    if( row_j <= 32 ){
	      cell_j = col_j + 32*(row_j - 1) - 1;
	    } else {
	      cell_j = 32*32 + col_j + 30*(row_j - 33) - 1;
	    }
	    */
	    
	    //updated for HCal
	    int cell_j;
	    cell_j = col_j + 12*row_j;
	    
	    // cout << "j, row_j, col_j, cell_j=" << jhit << ", " 
	    // 	 << row_j << ", " << col_j << ", " 
	    // 	 << cell_j << endl;

	    //increment Matrix element M_{i,j} with Ai * Aj / E_e, where, recall A_i and A_j are ADC values of hit i and hit j, and E_e is the predicted energy
	    // of the elastically scattered electron:
	    M(cell_i,cell_j) += Ahit_i * Ahit_j / E_e;
	  }
	  b(cell_i) += Ahit_i; //increment vector element i with ADC value of hit i. 
	  nevents[cell_i] += 1; //increment event counter for cell i:
	} 
      }
    }
  }
  //IF the number of events in a cell exceeds the minimum, then we can calibrate 

  // int smalld[ncell];

  // for(int i=0; i<ncell; i++){
  //   smalld[i] = 0;
  //   if( nevents[i] == 0 || M(i,i) <= 0.1*b(i) ){ //this cell has less than the minimal number of events required; set diagonal elements to 1, off-diagonals to zero:
  //     b(i) = 1.0;
  //     M(i,i) = 1.0;
  //     for(int j=0; j<ncell; j++){
  // 	if( j != i ) {
  // 	  M(i,j) = 0.0;
  // 	  M(j,i) = 0.0;
  // 	}
  //     }
      
  //     if( nevents[i] > 0 ) smalld[i] = 1;
      
  //   }
  // }

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

  //solution.Print();
  
  elist->Delete();

  TString textfilename = outfilename;
  textfilename.ReplaceAll(".root",".param");

  ofstream constants_file(textfilename.Data());

  //No idea where these constants come from - SS
  
  TH1D *hconstants_chan = new TH1D("hconstants_chan","",ncell,0.5,ncell+0.5);
  TH1D *hconstants = new TH1D("hconstants","",200,0.0,2.0);
  TH2D *hconstants_rowcol = new TH2D("hconstants_rowcol","",32,0.5,32.5,56,0.5,56.5);
  
  float ncalibrated_prot = 0.0;
  float average_constant_prot = 0.0;
  float ncalibrated_rcs = 0.0;
  float average_constant_rcs = 0.0;

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
	average_constant_rcs += solution(i);
	ncalibrated_rcs += 1.0;
      }
      
      hconstants_chan->Fill( i+1, solution(i)*1000.0 );
      //hconstants_chan->SetBinError( i+1, 1000.0*sqrt(fabs(M(i,i))) );
      
      hconstants->Fill( solution(i)*1000.0 );

      //compute the ratio of new / old calibration constants:
      hratio->Fill( solution(i)/oldconst_avg );
      hratio_chan->Fill( i+1, solution(i)/oldconst_avg );

      int irow, icol;

      /*
      //Geometry dependent. Also, what does it mean for the cell to be off of the detector here?  - SS
      if( (i+1) <= 32*32 ){ //Protvino
        irow = i/32 + 1;
	icol = i%32 + 1;
      } else { //RCS
	irow = (i-32*32)/30 + 33;
	icol = (i-32*32)%30 + 1;
      }
      */
      
      if( (i+1) <= ncols*nrows ){
	irow = i/nrows;
	icol = i%ncols;
      }
      
      hconstants_rowcol->Fill( icol, irow, 1000.0*solution(i) );
      //hconstants_rowcol->SetBinError( icol, irow, 1000.0*sqrt(fabs(M(i,i))));
      hratio_rowcol->Fill( icol, irow, solution(i)/oldconst_avg );
    }
  }
  
  // average_constant_prot /= ncalibrated_prot;
  // average_constant_rcs /= ncalibrated_rcs;
  average_constant /= ncalibrated;
  
  constants_file << "bigcal_prot_cfac = " << endl;
  
  for( int i=0; i<ncell; i++ ){
    char cconst[20];
    //Geometry dependent. What is rcs? -SS
    /*
    if( badcells[i] == 0 && solution(i)*1000.0 >= min_const && solution(i)*1000.0 <= max_const ){ //if the cell was calibrated and the calibration result is between 0.1 and 10, write out the new constant to the file:
      sprintf( cconst, "%15.4g", solution(i) );
    } else { //Otherwise, use the average calibration constant:
      if( i < 32*32 ){
	sprintf( cconst, "%15.4g", average_constant_prot );
      } else {
	sprintf( cconst, "%15.4g", average_constant_rcs );
      }
    }
    */
    
    if( badcells[i] == 0 && solution(i)*1000.0 >= min_const && solution(i)*1000.0 <= max_const ){ //if the cell was calibrated and the calibration result is between 0.1 and 10, write out the new constant to the file:
      sprintf( cconst, "%15.4g", solution(i) );
    } else { //Otherwise, use the average calibration constant:
      if( i < nrows*ncols ){
	sprintf( cconst, "%15.4g", average_constant_prot );
      } else {
	sprintf( cconst, "%15.4g", average_constant_rcs );
      }
    }
    /*
    constants_file << cconst;
    if( i+1 != 32*32 && i+1 != ncell ) constants_file << ",";
    if( (i+1)%8 == 0 ) constants_file << endl;
    if( (i+1) == 32*32 ) constants_file << endl << "bigcal_rcs_cfac = " << endl;
    */
    constants_file << cconst;
    if( i+1 != nrows*ncols && i+1 != ncell ) constants_file << ",";
    if( (i+1)%8 == 0 ) constants_file << endl;
    if( (i+1) == nrows*ncols ) constants_file << endl << "HCal_rcs_cfac = " << endl;
  
  }
  

  //  constants_file << endl << "bigcal_prot_gain_cor = " << endl;
  constants_file << endl << "HCal_prot_gain_cor = " << endl;

  for(int i=0; i<ncell; i++){
    char cconst[20];
    sprintf( cconst, "%15.4f", 1.0 );
    /*
    constants_file << cconst;
    if( i+1 != 32*32 && i+1 != ncell ) constants_file << ",";
    if( (i+1)%8 == 0 ) constants_file << endl;
    if( i+1 == 32*32 ) constants_file << endl << "bigcal_rcs_gain_cor = " << endl;
    */
    constants_file << cconst;
    if( i+1 != nrows*ncols && i+1 != ncell ) constants_file << ",";
    if( (i+1)%8 == 0 ) constants_file << endl;
    if( i+1 == nrows*ncols ) constants_file << endl << "HCal_rcs_gain_cor = " << endl;
  }

  fout->Write();
  fout->Close();
}

//For shower/preshower energy calibration

void bigcal_shower_profiles( const char *setupfilename, const char *outputfilename ){

  double csxp=3.8152;
  double csyp=3.8092;
  double csxr=4.0199;
  double csyr=4.0162;

  TString outfilename = outputfilename;
  if( !outfilename.EndsWith(".root") ){
    outfilename.ReplaceAll(".","_");
    outfilename += ".root";
  }
  
  TFile *fout = new TFile( outfilename.Data(), "RECREATE" );

  TChain *C = new TChain("h9030");
  ifstream setupfile(setupfilename);
  if( !setupfile ) return;

  TString currentline;
  
  while( currentline.ReadLine(setupfile) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") )
      C->Add(currentline.Data());
  }

  TCut global_cut = "";
  while( currentline.ReadLine(setupfile) && !currentline.BeginsWith("endcut") ){
    if( !currentline.BeginsWith("#") )
      global_cut += currentline.Data();
  }

  TEventList *elist = new TEventList("elist");

  C->Draw( ">>elist", global_cut );

  bigcal_ntuple *T = new bigcal_ntuple(C);

  cout << "number of events passing global cut = " << elist->GetN() << endl;

  TClonesArray *xmom_histos = new TClonesArray("TH1D", 29);
  TClonesArray *ymom_histos = new TClonesArray("TH1D", 29);

  for(int i=1; i<=29; i++){
    //TString histname ="h";
    TString histname;
    histname.Form("h%d",i);
    new( (*xmom_histos)[i-1] ) TH1D(histname.Data(),"",100,-1.0,1.0);
    
    histname.Form("h%d",100+i);
    new( (*ymom_histos)[i-1] ) TH1D(histname.Data(),"",100,-1.0,1.0);

  }
  
  long nevent=0;

  while( T->GetEntry(elist->GetEntry(nevent++))){
    if( nevent%1000 == 0 ) cout << nevent << endl;
    int best = int(T->ibest)-1;
    int ncellclust = int(T->ncellclust[best]);

    int rowmax = (T->iycell)[best][0];
    int colmax = (T->ixcell)[best][0];

    bool edge_max = false;
    if( rowmax <= 2 || rowmax >= 55 || colmax <= 2 || (rowmax<=32&&colmax>=31) || (rowmax > 32 && colmax>=29) ){
      edge_max = true;
    }

    if( !edge_max ){
      int xsection = (colmax-1)/8 + 1;
      int ysection = (rowmax-1)/8 + 1;

      int section = xsection + 4*(ysection-1);
      
      double csx = csxp, csy = csyp;
      if( rowmax > 32 ){
	csx = csxr;
	csy = csyr;
      }

      ( (TH1D*) (*xmom_histos)[section-1] )->Fill( T->xmoment[best]/csx );
      ( (TH1D*) (*ymom_histos)[section-1] )->Fill( T->ymoment[best]/csy );

      ( (TH1D*) (*xmom_histos)[28] )->Fill(   T->xmoment[best]/csx );
      ( (TH1D*) (*ymom_histos)[28] )->Fill(   T->ymoment[best]/csy );

    }
    
  }

  xmom_histos->Write();
  ymom_histos->Write();
  
  fout->Close();

}
