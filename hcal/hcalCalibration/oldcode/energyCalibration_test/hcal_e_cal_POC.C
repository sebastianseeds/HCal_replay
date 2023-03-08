//SSeeds 10.7.21 - Code modified from Andrew Puckett's BigCal calibration script to calibrate HCal energy deposition by PMT module. Currently proof of concept, tested with Eric Fuchey's digitized g4sbs data and with assumptions pending data from beam.

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
#include <ctime>
#include <iomanip>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TGraphErrors.h"

TChain *T = 0;

const int ncell = 288; 
const int kNrows = 24;
const int kNcols = 12;

const double PI = TMath::Pi();
const double M_e = 0.00051;
const double M_p = 0.938272;
const double M_n = 0.939565;
const double c_light = 299792458.0;

void hcal_e_cal_POC( int run = -1 ){
  
  //User input .root file for sample analysis
  cout << "Enter run number for toy analysis." << endl;
  cin >> run;
  
  double E_e = 0;
  cout << "Enter incident electron (beam) energy in GeV" << endl;
  cin >> E_e;
  
  int debug = 0;
  cout << "Debug mode?" << endl;
  cin >> debug;
  
  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );
  
  //Declare outfile
  TFile *foutTEST = new TFile( Form("hcalECal_%dTEST.root",run), "RECREATE" );  
  
  //initialize vectors and arrays
  vector<int> hitlist;
  vector<double> energy;
  double gHV[kNrows][kNcols];
  double gAlphas[kNrows][kNcols];
  double targetHV[kNrows][kNcols];
  double cRatio[kNrows][kNcols];

  //initialize histograms
  TH1D *W2_p = new TH1D("W2 Proton","W2_p",500,-10,10);
  W2_p->GetXaxis()->SetTitle( "GeV" );
  TH1D *Energy_ep = new TH1D("Scattered Electron Energy","Energy_ep",500,-2,8);
  Energy_ep->GetXaxis()->SetTitle( "GeV" );
  TH1F *CvChannel = new TH1F( "CvChannel", "CConst vs Channel", kNcols*kNrows, 0, kNcols*kNrows );
  TH1F *NEVvChannel = new TH1F( "NEVvChannel", "Number Events vs Channel", kNcols*kNrows, 0, kNcols*kNrows );
  TH1D *Phi_p = new TH1D("Scattering Angle Phi: proton","Phi_p",100,0,45);
  Phi_p->GetXaxis()->SetTitle( "Phi" );
  TH1D *Theta_e = new TH1D("Scattering Angle Theta: electron","Theta_e",100,0,45);
  Theta_e->GetXaxis()->SetTitle( "Theta" );
  TH1D *Phi_e = new TH1D("Scattering Angle Phi: electron","Phi_e",100,-10,10);
  Phi_e->GetXaxis()->SetTitle( "phi" );

  //Empirical limits, GeV/pC, may need updating
  double min_const=0.1, max_const=10.0;
  
  // Read in data produced by analyzer in root format
  cout << "Reading replayed root file.." << endl;
  if( !T ) { 
    T = new TChain( "T" );
    //T->Add( Form( "/volatile/halla/sbs/seeds/rootfiles/hcal_%d*.root", run ) );
    T->Add( "/volatile/halla/sbs/seeds/rootfiles/replayed_simdigtest_2_20211004.root" );
    T->SetBranchStatus( "*", 0 );
    T->SetBranchStatus( "sbs.hcal.*", 1 );
    T->SetBranchAddress( "sbs.hcal.a", hcalt::a );
    T->SetBranchAddress( "sbs.hcal.a_amp", hcalt::a_amp );
    T->SetBranchAddress( "sbs.hcal.a_p", hcalt::a_p );
    T->SetBranchAddress( "sbs.hcal.a_c", hcalt::a_c );
    T->SetBranchAddress( "sbs.hcal.a_amp_p", hcalt::a_amp_p );
    T->SetBranchAddress( "sbs.hcal.tdc", hcalt::tdc );
    T->SetBranchAddress( "sbs.hcal.ped", hcalt::ped );
    T->SetBranchAddress( "sbs.hcal.adcrow", hcalt::row );
    T->SetBranchAddress( "sbs.hcal.adccol", hcalt::col );
    T->SetBranchStatus( "Ndata.sbs.hcal.adcrow", 1 );
    T->SetBranchAddress( "Ndata.sbs.hcal.adcrow", &hcalt::ndata );

    // Add clustering branches
    T->SetBranchAddress("sbs.hcal.clus.id",hcalt::cid); //ID of highest energy block in cluster
    T->SetBranchAddress("sbs.hcal.clus.row",hcalt::crow);
    T->SetBranchAddress("sbs.hcal.clus.col",hcalt::ccol);
    T->SetBranchAddress("sbs.hcal.clus.e",hcalt::ce);
    T->SetBranchAddress("sbs.hcal.clus.eblk",hcalt::ceblk); //Energy highest energy block in cluster
    T->SetBranchAddress("sbs.hcal.clus.nblk",hcalt::cnblk);
    T->SetBranchAddress("sbs.hcal.nclus",hcalt::nclus);
    T->SetBranchAddress("sbs.hcal.nblk",hcalt::nblk);
    T->SetBranchAddress("sbs.hcal.clus_blk.id",hcalt::cblkid);

    // Add BB cal vars
    T->SetBranchAddress("bb.tr.px",hcalt::BBtr_px);
    T->SetBranchAddress("bb.tr.py",hcalt::BBtr_py);
    T->SetBranchAddress("bb.tr.pz",hcalt::BBtr_pz);
    
    cout << "Opened up tree with nentries=" << T->GetEntries() << endl;
  }


  ifstream HVFile( Form( "setFiles/HV_run%d.txt", run ) );
  
  if( !HVFile ){
    cerr << "No HV settings from run present -> setFiles/HV_run" << run << ".txt expected." << endl;
    return 0;
  }
  
  cout << "Getting HV settings for each pmt for run " << run << "." << endl;
  
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
  
  ifstream alphaFile( "setFiles/alphas.txt" );
  
  if( !alphaFile ){
    cerr << "No PMT alphas file present -> setFiles/alphas.txt expected." << endl;
    return 0;
  }
  
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

  //Need to get electron energy and W for cuts
  long mevent=0;

  cout << "Preliminary loop over all events in tree to get physics quantities commencing.." << endl;

  // Add W2 histogram for analysis

  while( T->GetEntry( mevent++ ) ){ 
    
    if( mevent%1000 == 0 ) cout << "Preliminary loop: " << mevent << "/" << T->GetEntries() << endl;

    //Get electron energy
    //Initialize vectors
    TVector3 lVep(hcalt::BBtr_px[0],hcalt::BBtr_py[0],hcalt::BBtr_pz[0]); //3Vector of scattered electron in lab frame, zero element most probable track

    if(debug==1) cout << "px = " << hcalt::BBtr_px[0] << ", py = " << hcalt::BBtr_py[0] << ", pz = " << hcalt::BBtr_pz[0] << ", lVep.Mag() = " << lVep.Mag() << "." << endl;
    
    //TLorentzVector pLVep(0.0,0.0,0.0,M_e); //Lorentz vector of track scattered electron in electron frame
    //Boost back to lab frame
    //pLVep.Boost(-lVep); //Doesn't work!
    //Take scattered electron energy from energy component
    //double E_ep = pLVep.E(); //nan
    
    double E_ep = sqrt(pow(M_e,2) + lVep.Mag2());  //Energy is obtained brute force since I cannot figure out Boost()

    if(debug==1) cout << sqrt(pow(M_e,2) + lVep.Mag2()) << " = sqrt(" <<  pow(M_e,2) << " + " << lVep.Mag2() << "))" << endl;

    //if(debug==1) cout << pLVep.Z() << " " << pLVep.X() << " " << pLVep.Y() << endl;
    
    if(debug==1) cout << "E_ep = " << E_ep << endl;
    
    //Get magnitude of scattered electron momentum
    double p_ep = lVep.Mag();

    double theta = acos((hcalt::BBtr_pz[0])/p_ep) * 180.0 / PI;
    Theta_e->Fill(theta);
    
    double phi_e = asin(hcalt::BBtr_py[0]/(p_ep*sin(theta)));
    Phi_e->Fill(phi_e);

    if(debug==1) cout << "p_ep = " << p_ep << endl;
    
    //Get Q2 from beam energy, outgoing electron energy, and momenta
    double Q2_ep = 2*E_e*E_ep*(1-(lVep.Z()/p_ep));
    
    if(debug==1) cout << "Cos theta = " << lVep.Z()/p_ep << endl;
    if(debug==1) cout << "Q2 = " << Q2_ep << endl;

    //Get W2 from Q2
    double W2_ep = pow(M_p,2)+2*M_p*(E_e-E_ep)-Q2_ep;

    if(debug==1) cout << "W2 = " << W2_ep << endl;

    W2_p->Fill( W2_ep );
    Energy_ep->Fill( E_ep );
    
    //W2 should be equal to Mp2 plus some smearing estimated at 30 MeV for now. Need to plot and get an informed value with elastic events.
    //if(W2_ep > pow(M_p,2)-0.03 && W2_ep < pow(M_p,2)+0.03){ 
    if(W2_ep > 0 && W2_ep < 2){ //Totally arbitrary for testing inelastic data until we have elastic data
      hitlist.push_back( mevent );
      energy.push_back( E_ep );
    }

    //Get proton momentum in zero SBS field
    double nu = E_e-E_ep;

    //Use the electron kinematics to predict the proton momentum assuming elastic scattering on free proton at rest:
    double p_p = sqrt(Q2_ep*(1.+Q2_ep/(4.*pow(M_p,2))));
    
    TVector3 lVpp(-lVep.X(),-lVep.Y(),E_e-lVep.Z());
    TLorentzVector lLVpp(lVpp, M_p+nu);

    //Each component in the form p_p_z = p_p*cos(phi)
    double phi_p = acos((E_e-hcalt::BBtr_pz[0])/p_p) * 180.0 / PI;
    Phi_p->Fill(phi_p);

  }

  foutTEST->cd();
  W2_p->SetTitle( "W2" );
  W2_p->Write( "W2" );
  W2_p->Draw( "AP" );

  foutTEST->cd();
  Energy_ep->SetTitle( "Scattered E_e" );
  Energy_ep->Write( "E_e" );
  Energy_ep->Draw( "AP" );

  cout << "Cut and electron energy histograms written to file." << endl;
  cout << "Number of events passing global cut = " << hitlist.size() << endl;
  
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
  
  int min_events_per_cell=1; //require at least 100 events to calculate a new calibration constant for the cell:
  
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
  
  /*
  //Histograms from Andrew's code. Will keep for reference.
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
  */
  //////////////////////////////
  //Main loop over the events://
  //////////////////////////////
  
  //Loop over only elastic events as passed by BB
  //Assume that a clean sample of elastic events has already been selected by global_cut:
  //In this case, we are trying to minimize chi^2 = sum_i=1,nevent (E_i - sum_j c_j A^i_j)^2/sig_E^2, where sig_E_i^2 ~ E_i; in this case, dchi^2/dc_k = 2 * sum_i=1,nevent 1/E_i * (E_i - sum_j c_j A_j) * A_k, set the difference on the RHS to zero to minimize.
  // This is a system of linear equations for c_k, with RHS = sum_i A_k and LHS = sum_i sum_j A_j c_j A_k/E_i 
  
  cout << "Beginning main loop over events.." << endl;

  //int nrows = kNrows;
  //int ncols = kNcols;

  while( T->GetEntry( hitlist[nevent++] ) ){ 
    
    if(hitlist[nevent]>T->GetEntries()) continue;

    if(debug==1) cout << "hitlist event: " << hitlist[nevent] << endl;
    if(debug==1) cout << "electron energy: " << energy[nevent] << endl;

    if( nevent%1000 == 0 ) cout << "Main loop: " << nevent << "/" << hitlist.size() << endl;
    
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
      e[r][c] = hcalt::a_c[m];  //Assuming this is the previously calibrated ADC
      
      // Mark saturated array when amplitude meets max RAU
      if( amp[r][c] > 4095 ) {
	saturated = 1;
      }
    }
    
    int best = hcalt::cid[0]; //index in the array of the "best" block in the best cluster found in this event. Currently only one cluster is recorded per event

    if(debug==1) cout << "best = " << best << endl;

    int ncellclust = hcalt::nblk[0]; //total number of hits in the "best" cluster, default to 9 (3x3) for now
    
    if(debug==1) cout << "nblk = " << ncellclust << endl;

    if( best >= 0 && ncellclust >= 1 && saturated == 0 ){ // require at least 2 hits in this cluster to use as a "good" event for calibration purposes:

      double E_e = energy[nevent]; //Energy of e' from first loop

      //Check whether this cluster has its maximum at the edge of the calorimeter
      int rowmax = best/kNcols;
      int colmax = best%kNcols;

      //if(debug==1) cout << "rowmax " << rowmax << " colmax " << colmax << endl; 

      bool edge_max = false;

      //Corrected for HCal geometry
      if( rowmax <= 2 || rowmax >= 23 || colmax <= 2 || colmax >= 11 ){
	edge_max = true;
      }
     
      if( !edge_max ){  //Only consider clusters with maximum at least one block away from the edge for calibration:
	for(int ihit=0; ihit<ncellclust; ihit++ ){ //outer loop over all cells in the cluster: 

	  int cell_i = hcalt::cblkid[ihit]-1;

	  if(debug==1) cout << "cell = " << cell_i << endl;

	  int row_i = cell_i/kNcols;
	  int col_i = cell_i%kNcols;

	  if(debug==1) cout << "r_i:c_i = " << row_i << ":" << col_i << endl;

	  float Ahit_i = adc_p[row_i][col_i]; //Ahit_i is the ADC value for this hit (ped-subtracted)
	  float Ehit_i = e[row_i][col_i]; //Ehit_i is the reconstructed energy for this hit (using previous calibration constants)

	  if(debug==1) cout << "Ahit_i = " << Ahit_i << "; Ehit_i = " << Ehit_i << endl;

	  //Check old calibration constants:
	  oldconstants_sum[cell_i] += Ehit_i/Ahit_i;
	  oldconstants_nevent[cell_i] += 1.0;
	  oldconstants_sum2[cell_i] += pow(Ehit_i/Ahit_i,2);
	  
	  for(int jhit=0; jhit<ncellclust; jhit++ ){ //inner loop over all cells in the cluster: 

	    int cell_j = hcalt::cblkid[jhit]-1;

	    int row_j = cell_j/kNcols;
	    int col_j = cell_j%kNcols;

	    float Ahit_j = adc_p[row_j][col_j]; //ADC value of hit j 
	    
	    //increment Matrix element M_{i,j} with Ai * Aj / E_e, where, recall A_i and A_j are ADC values of hit i and hit j, and E_e is the predicted energy of the elastically scattered electron:

	    if(debug==1) cout << "r_j:c_j; cell_i:cell_j = " << row_j << ":" << col_j << "; " << cell_i << ":" << cell_j << endl;


	    if( cell_i >= 288 || cell_j >= 288) {
	      cout << "Error: Cell number greater than expected geometry." << endl;
	      return 0;
	    }
	    
	    M(cell_i,cell_j) += Ahit_i * Ahit_j / E_e;
	    
	  }
	  b(cell_i) += Ahit_i; //increment vector element i with ADC value of hit i. 
	  nevents[cell_i] += 1; //increment event counter for cell i

	  if(b(cell_i) < 0) cout << "b vector element < 0, i=" << cell_i << ": " << b(cell_i) << endl; 

	} 
      }
    }
  }

  cout << "Matrix populated.." << endl;

  ////////////////
  ////CRITICAL////
  ////////////////
  //Critical debug and possible loss of function - this section should be removed before production. Suspect that not enough data yet exists to populate the matrix effectively
  //if(debug==1){
  /*
  for( int i=0; i<ncell; i++ ){
    for( int j=0; j<ncell; j++){

      //if(M(i,j)<0) cout << "Matrix element below zero: " << M(i,j) << ".." << endl;
      
      if(M(i,j)==0) M(i,j) = 1.0;

      if(debug==1) cout << "Final matrix elements i:j -> " << i << ":" << j << " = " << M(i,j) << endl; 

    }
    if(debug==1) cout << "Final ADC sum vector elements i -> " << i << " = " << b(i) << endl; 
  }
  */
  //}

  //IF the number of events in a cell exceeds the minimum, then we can calibrate 

  //Make adjustments to the matrix M and vector b to exclude "bad" cells from calibration: 
  int badcells[ncell];
  double xErr[ncell];

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

  cout << "..done." << endl;
  
  ofstream constants_file;

  constants_file.open( Form( "/work/halla/sbs/seeds/outFiles/hcalECalParams_%d.txt", run ) );

  TH1D *hconstants_chan = new TH1D("hconstants_chan","",ncell,0.5,ncell+0.5);
  TH1D *hconstants = new TH1D("hconstants","",200,0.0,2.0);
  TH2D *hconstants_rowcol = new TH2D("hconstants_rowcol","",kNcols,0.5,kNcols+0.5,kNrows,0.5,kNrows+0.5);
  
  float ncalibrated_prot = 0.0;
  float average_constant_prot = 0.0;

  TH1D *hold_chan = new TH1D("hold_chan","",ncell, 0.5, ncell+0.5 );
  TH1D *hold_chan_rms = new TH1D("hold_chan_rms","",ncell,0.5,ncell+0.5);
  TH1D *hratio = new TH1D("hratio","new/old",200,0.5, 1.5);
  TH1D *hratio_chan = new TH1D("hratio_chan","new/old by chan.",ncell,0.5, ncell+0.5);
  TH2D *hratio_rowcol = new TH2D("hratio_rowcol","new/old by row/col",kNcols,0.5,kNcols+0.5,kNrows,0.5,kNrows+0.5);

  //loop over all the cells, fill diagnostic histograms, and write out new constants:

  double y[ncell] = {0.0};

  double E_targ = 1000.0; //Previously set to 1000.0, may need to be set with the old target pC in mind (406 pC)

  for(int i=0; i<ncell; i++){

    y[i] = i;

    xErr[i] = E_targ*sqrt(fabs(M(i,i)));

    NEVvChannel->SetBinContent(i,nevents[i]);

    int r = i/kNcols;
    int c = i%kNcols;

    if( nevents[i] >= min_events_per_cell ){
      
      //If this cell was calibrated:
      //if( badcells[i] == 0 && solution(i)*E_targ >= min_const && solution(i)*E_targ <= max_const ){
	
	//We multiply by 1000 here because the energy units were in GeV and the PMT HVs were typically set to give
	// 1,000 ADC = 1 GeV, so the constants were typically of order 10^-3. Here we multiply them by 1,000 to give a number of O(1).  
	//We have max_adc = 4095 and max_E for GMn = 700 MeV. May need to update as comparitive orders of magnitude are 10 not 10^3
	//Now we have 406 pC and max_E for GMn = 0.7 GeV
	
	double oldconst_avg = oldconstants_sum[i]/oldconstants_nevent[i];
	double oldconst_rms = sqrt(fabs( oldconstants_sum2[i]/oldconstants_nevent[i] - pow(oldconst_avg,2) ) );
	
	cRatio[r][c] = solution[i]/( oldconst_avg*E_targ );

	if( cRatio[r][c] <= 0 ) cout << "cRatio[" << r << "][" << c << "] = " << solution[i] << "/" << oldconst_avg << "*" << E_targ << endl;

	hold_chan->Fill(i+1, oldconst_avg*E_targ );
	hold_chan_rms->Fill(i+1, oldconst_rms*E_targ);
	
	//compute average calibration constants for both sections of the calorimeter and use the average result for the channels which could not be calibrated:
	
	if( i < kNcols*kNrows ){
	  average_constant_prot += solution(i);
	  ncalibrated_prot += 1.0;
	} else {
	  cout << "Error: looping outside of HCal geometry" << endl;
	}
	
	hconstants_chan->Fill( i+1, solution(i)*E_targ );
	hconstants_chan->SetBinError( i+1, E_targ*sqrt(fabs(M(i,i))) );
	
	hconstants->Fill( solution(i)*E_targ );
	
	//compute the ratio of new / old calibration constants:
	hratio->Fill( solution(i)/oldconst_avg );
	hratio_chan->Fill( i+1, solution(i)/oldconst_avg );
	
	int irow, icol;
	
	if( (i+1) <= kNcols*kNrows ){
	  irow = i/kNcols;
	  icol = i%kNcols;
	}
	
	hconstants_rowcol->Fill( icol, irow, E_targ*solution(i) );
	//hconstants_rowcol->SetBinError( icol, irow, E_targ*sqrt(fabs(M(i,i))));
	hratio_rowcol->Fill( icol, irow, solution(i)/oldconst_avg );
      }
    }
  //}

  average_constant_prot /= ncalibrated_prot;
  
  constants_file << "HCal_prot_cfac = " << endl;
  
  for( int i=0; i<ncell; i++ ){
    char cconst[20];
    
    if( badcells[i] == 0 && solution(i)*E_targ >= min_const && solution(i)*E_targ <= max_const ){ //if the cell was calibrated and the calibration result is between 0.1 and 10, write out the new constant to the file:
      sprintf( cconst, "%15.4g", solution(i) );
    } else { //Otherwise, use the average calibration constant:
      if( i < kNrows*kNcols ){
	sprintf( cconst, "%15.4g", average_constant_prot );
      } else {
	cout << "Error: looping outside of normal geometry for cconst" << endl;
      }
    }

    if(debug==1) cout << "Constant: " << cconst << endl;

    CvChannel->SetBinContent( i+1, solution(i) );

    constants_file << cconst;
    if( i+1 != kNrows*kNcols && i+1 != ncell ) constants_file << ",";
    if( (i+1)%kNcols == 0 ) constants_file << endl;
    if( i == kNrows*kNcols ) constants_file << "Error: ncell > kNrows * kNcols" << endl;
  }
  
  constants_file << endl << "HCal_prot_gain_cor = " << endl;

  for(int i=0; i<ncell; i++){

    int r = i/kNcols;
    int c = i%kNcols;

    targetHV[r][c] = gHV[r][c]*pow(cRatio[r][c],(1/gAlphas[r][c]));

    char cconst[20];
    sprintf( cconst, "%15.4f", targetHV[r][c] );

    if(debug==1) cout << "Target HV: " << cconst << endl;

    constants_file << cconst;
    if( i+1 != kNrows*kNcols && i+1 != ncell ) constants_file << ",";
    if( (i+1)%kNcols == 0 ) constants_file << endl;
    if( i == kNrows*kNcols ) constants_file << "Error: ncell > kNrows * kNcols" << endl;
  }

  double yErr[ncell] = {0.0};
  //TGraphErrors *ccgraph = new TGraphErrors( ncell, &solution[0], y, xErr, yErr );
  TGraphErrors *ccgraph = new TGraphErrors( ncell, y, &solution[0], yErr, xErr );
  ccgraph->GetXaxis()->SetLimits(0.0,288.0);

  //Write out root file for diagnostics and close
  foutTEST->cd();
  CvChannel->SetTitle( "CConst vs Channel" );
  CvChannel->Write( "CvChannel" );
  CvChannel->Draw( "AP" );

  Phi_p->SetTitle( "Proton Scattering Angle Phi" );
  Phi_p->Write( "Phi_p" );
  Phi_p->Draw( "AP" );

  Theta_e->SetTitle( "Electron Scattering Angle Theta" );
  Theta_e->Write( "Theta_e" );
  Theta_e->Draw( "AP" );

  Phi_e->SetTitle( "Electron Scattering Angle Phi" );
  Phi_e->Write( "Phi_e" );
  Phi_e->Draw( "AP" );

  ccgraph->SetTitle("Calibration Constants");
  ccgraph->GetXaxis()->SetTitle("Channel");
  ccgraph->GetXaxis()->SetTitle("GeV/pC");
  ccgraph->SetMarkerStyle(21); //Boxes
  ccgraph->Write("AP");
  ccgraph->Draw("AP");

  NEVvChannel->SetTitle( "CConst vs Channel" );
  NEVvChannel->Write( "NEVvChannel" );
  NEVvChannel->Draw( "AP" );

  hconstants_chan->Write( "hconstants_chan" );
  hconstants_chan->Draw( "AP" );

  hconstants->Write( "hconstants" );
  hconstants->Draw( "AP" );

  hconstants_rowcol->Write( "hconstants_rowcol" );
  hconstants_rowcol->Draw( "AP" );

  hold_chan->Write( "hold_chan" );
  hold_chan->Draw( "AP" );

  hold_chan_rms->Write( "hold_chan_rms" );
  hold_chan_rms->Draw( "AP" );
  
  hratio->Write( "hratio" );
  hratio->Draw( "AP" );
  
  hratio_chan->Write( "hratio_chan" );
  hratio_chan->Draw( "AP" );

  hratio_rowcol->Write( "hratio_rowcol" );
  hratio_rowcol->Draw( "AP" );
  foutTEST->Close();

  st->Stop();

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

}
