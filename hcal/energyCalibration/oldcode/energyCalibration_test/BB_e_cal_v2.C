//SSeeds 10.12.21 - Production - Code modified from Andrew Puckett's BigCal calibration script to calibrate HCal energy deposition by PMT module. Currently proof of concept, tested with Eric Fuchey's digitized g4sbs data and with assumptions pending data from beam.

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

const int ncell = 241; 
const int kNrowsSH = 27;
const int kNcolsSH = 7;
const int kNrowsPS = 26;
const int kNcolsPS = 2;

const double PI = TMath::Pi();
const double M_e = 0.00051;
const double M_p = 0.938272;
const double M_n = 0.939565;
const double c_light = 299792458.0;

void BB_e_cal_v2( int run = -1 ){
  
  //User input .root file for sample analysis
  cout << "Enter run number for toy analysis." << endl;
  cin >> run;
  
  double E_e = 0;
  cout << "Enter incident electron (beam) energy in GeV" << endl;
  cin >> E_e;
  
  int debug = 0;
  cout << "Debug mode?" << endl;
  cin >> debug;
  
  //double HCal_d = 10.0; //Distance to HCal from scattering chamber for comm1
  //double HCal_th = 39.0; //Angle that the center of HCal is at
  //double dBlock = 0.15; //Transverse dimensions of HCal blocks (x and y)
  
  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );
  
  //Declare outfile
  TFile *fout = new TFile( Form("BBECal_%d.root",run), "RECREATE" );
  
  //initialize vectors and arrays
  //vector<int> eventlist;
  //vector<double> energy_e;
  double gHV[kNrowsSH][kNcolsSH];
  double gAlphas[kNrowsSH][kNcolsSH];
  double targetHV[kNrowsSH][kNcolsSH];
  double cRatio[kNrowsSH][kNcolsSH];
  
  //initialize histograms
  TH1D *W2_p = new TH1D("W2 Proton","W2_p",500,-10,10);
  W2_p->GetXaxis()->SetTitle( "GeV" );
  TH1D *Energy_ep = new TH1D("Scattered Electron Energy","Energy_ep",500,-2,8);
  Energy_ep->GetXaxis()->SetTitle( "GeV" );
  TH1D *Energy_pp = new TH1D("Scattered Proton Energy","Energy_pp",500,-2,15);
  Energy_pp->GetXaxis()->SetTitle( "GeV" );
  TH1F *CvChannel = new TH1F( "CvChannel", "CConst vs Channel", kNcolsSH*kNrowsSH, 0, kNcolsSH*kNrowsSH );
  TH1F *NEVvChannel = new TH1F( "NEVvChannel", "Number Events vs Channel", kNcolsSH*kNrowsSH, 0, kNcolsSH*kNrowsSH );
  TH1F *badCvChannel = new TH1F( "badCvChannel", "Bad Cell vs Channel", kNcolsSH*kNrowsSH, 0, kNcolsSH*kNrowsSH );
  TH1F *CRatiovChannel = new TH1F( "CRatiovChannel", "Constant Ratio (new/old) vs Channel", kNcolsSH*kNrowsSH, 0, kNcolsSH*kNrowsSH );
  TH1D *Phi_p = new TH1D("Scattering Angle Phi: proton","Phi_p",100,0,45);
  Phi_p->GetXaxis()->SetTitle( "Phi" );
  TH1D *Theta_e = new TH1D("Scattering Angle Theta: electron","Theta_e",100,0,45);
  Theta_e->GetXaxis()->SetTitle( "Theta" );
  TH1D *Phi_e = new TH1D("Scattering Angle Phi: electron","Phi_e",100,-10,10);
  Phi_e->GetXaxis()->SetTitle( "phi" );
  //TH1D *hcal_col = new TH1D("Expected Col","Col",100,-50,50);
  //Phi_e->GetXaxis()->SetTitle( "col" );
  //TH1D *hcal_row = new TH1D("Expected Row","Row",100,-50,50);
  //Phi_e->GetXaxis()->SetTitle( "row" );
  TH1D *p_prot = new TH1D("Proton Momentum","p_prot",100,-2,16);
  p_prot->GetXaxis()->SetTitle( "GeV" );
  TH1D *BBbesttr = new TH1D("BB Track with Lowest Chi^2","BBbesttr",10,0,10);
  BBbesttr->GetXaxis()->SetTitle( "Track Index" );
  TH1D *BBtr_chi2 = new TH1D("Chi^2 among all BB tracks","BBtr_chi2",100,0,50);
  BBtr_chi2->GetXaxis()->SetTitle( "Chi^2" );
  TH1D *BBtr_bestchi2 = new TH1D("Chi^2 among best BB tracks","BBtr_bestchi2",100,0,50);
  BBtr_bestchi2->GetXaxis()->SetTitle( "Chi^2" );
  
  
  //ADC histograms
  //Set ADC spectrum histogram limits
  int SpecBins = 520;
  //INT
  double IntSpec_min = -20.0, IntSpec_max = 500.0; 
  //Build histograms
  TH1F *gIntSpec[kNrowsSH][kNcolsSH];
  for( int r=0; r<kNrowsSH; r++ ){
    for( int c=0; c<kNcolsSH; c++ ){  
      int m = r*12+c;
      gIntSpec[r][c] = new TH1F( Form( "Int ADC Spect R%d C%d PMT%d", r, c, m ), Form( "Int ADC Spect R%d C%d PMT%d", r, c, m ), SpecBins, IntSpec_min, IntSpec_max );
      gIntSpec[r][c]->GetXaxis()->SetTitle( "pC" );
      gIntSpec[r][c]->GetXaxis()->CenterTitle();
    }
  }
  
  
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
    T->SetBranchAddress("bb.tr.p",hcalt::BBtr_p); //Magnitude momentum BB
    T->SetBranchAddress("bb.tr.chi2",hcalt::BBtr_chi2); //Chi^2 per track in BB
    T->SetBranchAddress("bb.tr.n",hcalt::BBtr_n); //number of tracks in BB
    T->SetBranchStatus( "Ndata.bb.sh.adcrow", 1 );
    T->SetBranchAddress( "Ndata.bb.sh.adcrow", &hcalt::BBndata );
    T->SetBranchAddress( "bb.sh.adcrow", hcalt::BBrow );
    T->SetBranchAddress( "bb.sh.adccol", hcalt::BBcol );
    T->SetBranchAddress( "bb.sh.a_p", hcalt::BBa_p );
    T->SetBranchAddress( "bb.sh.a_c", hcalt::BBa_c );
    T->SetBranchAddress( "bb.sh.a_amp", hcalt::BBa_amp );

    T->SetBranchAddress("bb.sh.clus.nblk",hcalt::BBcnblk);
    T->SetBranchAddress("bb.sh.nclus",hcalt::BBnclus);
    T->SetBranchAddress("bb.sh.nblk",hcalt::BBnblk);
    T->SetBranchAddress("bb.sh.clus_blk.id",hcalt::BBcblkid);
    T->SetBranchAddress("bb.sh.clus.id",hcalt::BBcid); //ID of highest energy block in cluster

    cout << "Opened up tree with nentries=" << T->GetEntries() << endl;
  }
  
  double energy_e[T->GetEntries()];

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
    
    rval = floor( n1/kNcolsSH );
    cval = n1 % kNcolsSH;
    
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
    
    rval = floor( n1/kNcolsSH );
    cval = n1 % kNcolsSH;
    
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
    int track_tot = (int)hcalt::BBtr_n[0];
    int track = 0;
    
    for(int elem=0; elem<track_tot; elem++){
      double min = 1000.0;
      
      BBtr_chi2->Fill(hcalt::BBtr_chi2[elem]);
      
      if(hcalt::BBtr_chi2[elem] < min)
	track = elem;      
    }
    
    BBbesttr->Fill(track);
    BBtr_bestchi2->Fill(hcalt::BBtr_chi2[track]);
    
    TVector3 lVep(hcalt::BBtr_px[track],hcalt::BBtr_py[track],hcalt::BBtr_pz[track]);
    //TVector3 lVep(hcalt::BBtr_px[0],hcalt::BBtr_py[0],hcalt::BBtr_pz[0]); //3Vector of scattered electron in lab frame, zero element most probable track
    
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

    //Get energy transfer
    double nu = E_e-E_ep;

    //Get W2 from Q2
    double W2_ep = pow(M_p,2)+2*M_p*nu-Q2_ep;

    if(debug==1) cout << "W2 = " << W2_ep << endl;

    W2_p->Fill( W2_ep );
    Energy_ep->Fill( E_ep );

    //Use the electron kinematics to predict the proton momentum assuming elastic scattering on free proton at rest:
    //double p_p = sqrt(Q2_ep*(1.+Q2_ep/(4.*pow(M_p,2))));

    //Get energy of the proton
    double E_pp = nu+M_p;
    Energy_pp->Fill( E_pp );

    double p_p = sqrt(pow(E_pp,2)+W2_ep);

    //double p_p = sqrt(pow(nu,2)-Q2_ep);
    p_prot->Fill(p_p);

    //TVector3 lVpp(-lVep.X(),-lVep.Y(),E_e-lVep.Z());
    //TLorentzVector lLVpp(lVpp, M_p+nu);

    //Each component in the form p_p_z = p_p*cos(phi)
    double phi_p = acos((E_e-hcalt::BBtr_pz[0])/p_p) * 180.0 / PI;
    Phi_p->Fill(phi_p);

    //Get the projection of each event onto the face of HCal - will need proper Jacobian
    //double HCal_z = HCal_d*cos(phi_p);
    //double HCal_x = HCal_d*cos(-phi_e)*sin(phi_p);
    //double HCal_y = HCal_d*sin(-phi_e)*sin(phi_p);

    //Zeroth order approx's
    //int xBlock = std::floor((HCal_d*phi_p - HCal_d*HCal_th)/dBlock)+6;
    //int yBlock = std::floor(-HCal_d*phi_e/dBlock)+12;

    //cout << "xBlock" << xBlock << " yBlock" << yBlock << endl;

    //hcal_col->Fill(xBlock);
    //hcal_row->Fill(yBlock);

    //double dist_z_to_HCal = HCal_d*sin(HCal_th);
    //int yBlock = 

    //W2 should be equal to Mp2 plus some smearing estimated at 30 MeV for now. Need to plot and get an informed value with elastic events.
    //if(W2_ep > pow(M_p,2)-0.03 && W2_ep < pow(M_p,2)+0.03){ 
    /*
    if(W2_ep > 0 && W2_ep < 2){ //Totally arbitrary for testing inelastic data until we have elastic data
      eventlist.push_back( mevent );
      energy_e( E_ep );
    }

    cout << eventlist.size() << endl;
    cout << energy_e.size() << endl;
    */
    //if(mevent<0) cout << mevent << " " << E_ep << endl;
   
    energy_e[mevent]=E_ep;
    cout << E_ep << endl;
  }
  
  cout << mevent << endl;
  /*
  fout->cd();
  W2_p->SetTitle( "W2" );
  W2_p->Write( "W2" );
  W2_p->Draw( "AP" );
  */
  //fout->cd();
  Energy_ep->SetTitle( "Scattered Electron Energy" );
  Energy_ep->Write( "E_ep" );
  Energy_ep->Draw( "AP" );
  /*
  Energy_pp->SetTitle( "Scattered Proton Energy" );
  Energy_pp->Write( "E_pp" );
  Energy_pp->Draw( "AP" );
  */
  cout << "Cut and electron energy histograms written to file." << endl;
  //cout << "Number of events passing global cut = " << eventlist.size() << endl;
  
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
  
  int min_events_per_cell=1000; //require at least 100 events to calculate a new calibration constant for the cell:
  
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

  //int nrows = kNrowsSH;
  //int ncols = kNcolsSH;

  //while( T->GetEntry( eventlist[nevent++] ) ){ 
  for( int nevent=0; nevent<mevent; nevent++ ){

    T->GetEntry(nevent);

    //if(eventlist[nevent]>T->GetEntries()) continue;

    //if(debug==1) cout << "eventlist event: " << eventlist[nevent] << endl;
    if(debug==1) cout << "electron energy: " << energy_e[nevent] << endl;

    //if( nevent%1000 == 0 ) cout << "Main loop: " << nevent << "/" << eventlist.size() << endl;
    
    int r,c;

    double adc_p[kNrowsSH][kNcolsSH] = {0.0};
    double e[kNrowsSH][kNcolsSH] = {0.0};
    double amp[kNrowsSH][kNcolsSH] = {0.0};

    int saturated = 0;

    // Process event with m data
    for( int m = 0; m < hcalt::BBndata; m++ ) {
      // Define row and column
      r = hcalt::BBrow[m];
      c = hcalt::BBcol[m];
      if( r < 0 || c < 0 ) {
	cerr << "Error: row or col negative." << endl;
	continue;
      }
      
      if( r >= kNrowsSH || c >= kNcolsSH ) continue;
      
      // Fill adc and energy arrays
      adc_p[r][c] = hcalt::BBa_p[m]; 
      amp[r][c] = hcalt::BBa_amp[m];
      e[r][c] = hcalt::BBa_c[m];  //Assuming this is the previously calibrated ADC
      
      // Mark saturated array when amplitude meets max RAU
      if( amp[r][c] > 4095 ) {
	saturated = 1;
	cout << "Saturation." << endl;
      }
    }
    
    int best = hcalt::BBcid[0]; //index in the array of the "best" block in the best cluster found in this event. Currently only one cluster is recorded per event

    if(debug==1) cout << "best = " << best << endl;

    int ncellclust = hcalt::BBnblk[0]; //total number of hits in the "best" cluster, default to 9 (3x3) for now
    
    if(debug==1) cout << "nblk = " << ncellclust << endl;

    if( best >= 0 && ncellclust >= 1 && saturated == 0 ){ // require at least 2 hits in this cluster to use as a "good" event for calibration purposes:

      double E_p = energy_e[nevent]; //Energy of e' from first loop

      //Check whether this cluster has its maximum at the edge of the calorimeter
      int rowmax = best/kNcolsSH;
      int colmax = best%kNcolsSH;

      //if(debug==1) cout << "rowmax " << rowmax << " colmax " << colmax << endl; 

      bool edge_max = false;

      //Corrected for HCal geometry
      if( rowmax <= 0 || rowmax >= 26 || colmax <= 0 || colmax >= 6 ){
	edge_max = true;
	cout << "Edgemax" << endl;
      }
     
      if( !edge_max ){  //Only consider clusters with maximum at least one block away from the edge for calibration:
	for(int ihit=0; ihit<ncellclust; ihit++ ){ //outer loop over all cells in the cluster: 

	  int cell_i = hcalt::BBcblkid[ihit]-1;

	  if(debug==1) cout << "cell = " << cell_i << endl;

	  int row_i = cell_i/kNcolsSH;
	  int col_i = cell_i%kNcolsSH;

	  if(debug==1) cout << "r_i:c_i = " << row_i << ":" << col_i << endl;

	  float Ahit_i = adc_p[row_i][col_i]; //Ahit_i is the ADC value for this hit (ped-subtracted)
	  float Ehit_i = e[row_i][col_i]; //Ehit_i is the reconstructed energy for this hit (using previous calibration constants)

	  gIntSpec[row_i][col_i]->Fill(Ahit_i);

	  if(debug==1) cout << "Ahit_i = " << Ahit_i << "; Ehit_i = " << Ehit_i << endl;

	  //Check old calibration constants:
	  oldconstants_sum[cell_i] += Ehit_i/Ahit_i;
	  oldconstants_nevent[cell_i] += 1.0;
	  oldconstants_sum2[cell_i] += pow(Ehit_i/Ahit_i,2);
	  
	  nevents[cell_i] += 1; //increment event counter for cell i
	  if(nevents[cell_i] > 2000) continue; //Added to keep statistics similar across all channels

	  for(int jhit=0; jhit<ncellclust; jhit++ ){ //inner loop over all cells in the cluster: 

	    int cell_j = hcalt::BBcblkid[jhit]-1;

	    int row_j = cell_j/kNcolsSH;
	    int col_j = cell_j%kNcolsSH;

	    float Ahit_j = adc_p[row_j][col_j]; //ADC value of hit j 
	    
	    //increment Matrix element M_{i,j} with Ai * Aj / E_p, where, recall A_i and A_j are ADC values of hit i and hit j, and E_p is the predicted energy of the elastically scattered electron:

	    if(debug==1) cout << "r_j:c_j; cell_i:cell_j = " << row_j << ":" << col_j << "; " << cell_i << ":" << cell_j << endl;


	    if( cell_i >= 288 || cell_j >= 288) {
	      cout << "Error: Cell number greater than expected geometry." << endl;
	      return 0;
	    }
	    
	    M(cell_i,cell_j) += Ahit_i * Ahit_j / E_p;
	    
	  }
	  b(cell_i) += Ahit_i; //increment vector element i with ADC value of hit i. 

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
  ofstream HVOutFile;

  constants_file.open( Form( "/work/halla/sbs/seeds/outFiles/hcalECalParams_%d.txt", run ) );
  HVOutFile.open( Form( "/work/halla/sbs/seeds/outFiles/ECalHVTargets_run%d.txt", run ) );
  HVOutFile << "Module" << '\t' << " HV" << endl;

  TH1D *Constants_chan = new TH1D("Constants_chan","",ncell,0.5,ncell+0.5);
  TH1D *Constants = new TH1D("Constants","",200,0.0,2.0);
  TH2D *Constants_rowcol = new TH2D("Constants_rowcol","",kNcolsSH,0.5,kNcolsSH+0.5,kNrowsSH,0.5,kNrowsSH+0.5);
  
  float ncalibrated_prot = 0.0;
  float average_constant_prot = 0.0;

  TH1D *Old_chan = new TH1D("Old_chan","",ncell, 0.5, ncell+0.5 );
  TH1D *Old_chan_rms = new TH1D("Old_chan_rms","",ncell,0.5,ncell+0.5);
  TH1D *Ratio = new TH1D("Ratio","new/old",200,0.5, 1.5);
  TH1D *Ratio_chan = new TH1D("Ratio_chan","new/old by chan.",ncell,0.5, ncell+0.5);
  TH2D *Ratio_rowcol = new TH2D("Ratio_rowcol","new/old by row/col",kNcolsSH,0.5,kNcolsSH+0.5,kNrowsSH,0.5,kNrowsSH+0.5);

  //loop over all the cells, fill diagnostic histograms, and write out new constants:

  double y[ncell] = {0.0};

  double E_targ = 1000; //Average peak for 14 MeV is 8 pC, so our conversion is 8/0.014 pC/GeV to get GeV. Similarly 14 MeV is 18 mV for maxADC (for reference)

  for(int i=0; i<ncell; i++){

    y[i] = i;

    xErr[i] = E_targ*sqrt(fabs(M(i,i)))/sqrt(nevents[i]); //Covariance error (std dev/root(N))

    NEVvChannel->SetBinContent(i,nevents[i]);

    int r = i/kNcolsSH;
    int c = i%kNcolsSH;

    if( nevents[i] >= min_events_per_cell ){
      
      //If this cell was calibrated:
      //if( badcells[i] == 0 && solution(i)*E_targ >= min_const && solution(i)*E_targ <= max_const ){
	
      //We multiply by 1000 here because the energy units were in GeV and the PMT HVs were typically set to give
      // 1,000 ADC = 1 GeV, so the constants were typically of order 10^-3. Here we multiply them by 1,000 to give a number of O(1).  
      //We have max_adc = 4095 and max_E for GMn = 700 MeV. May need to update as comparitive orders of magnitude are 10 not 10^3
      //Now we have 406 pC and max_E for GMn = 0.7 GeV
	
      double oldconst_avg = oldconstants_sum[i]/oldconstants_nevent[i];
      double oldconst_rms = sqrt(fabs( oldconstants_sum2[i]/oldconstants_nevent[i] - pow(oldconst_avg,2) ) );
	
      //cRatio[r][c] = solution[i]/( oldconst_avg*E_targ );
      cRatio[r][c] = solution[i]/( oldconst_avg );

      CRatiovChannel->SetBinContent(kNcolsSH*r+c+1,cRatio[r][c]);

      if( cRatio[r][c] <= 0 ) cout << "cRatio[" << r << "][" << c << "] = " << solution[i] << "/" << oldconst_avg << "*" << E_targ << endl;

      Old_chan->Fill(i+1, oldconst_avg*E_targ );
      Old_chan_rms->Fill(i+1, oldconst_rms*E_targ);
	
      //compute average calibration constants for both sections of the calorimeter and use the average result for the channels which could not be calibrated:
	
      if( i < kNcolsSH*kNrowsSH ){
	average_constant_prot += solution(i);
	ncalibrated_prot += 1.0;
      } else {
	cout << "Error: looping outside of HCal geometry" << endl;
      }
	
      Constants_chan->Fill( i+1, solution(i)*E_targ );
      Constants_chan->SetBinError( i+1, E_targ*sqrt(fabs(M(i,i)))/sqrt(nevents[i]) );
	

      Constants->Fill( solution(i)*E_targ );
	
      //compute the ratio of new / old calibration constants:
      Ratio->Fill( solution(i)/oldconst_avg );
      Ratio_chan->Fill( i+1, solution(i)/oldconst_avg );
	
      int irow, icol;
	
      if( (i+1) <= kNcolsSH*kNrowsSH ){
	irow = i/kNcolsSH;
	icol = i%kNcolsSH;
      }
	
      Constants_rowcol->Fill( icol, irow, E_targ*solution(i) );
      //Constants_rowcol->SetBinError( icol, irow, E_targ*sqrt(fabs(M(i,i))));
      Ratio_rowcol->Fill( icol, irow, solution(i)/oldconst_avg );
    }
  }
  //}

  average_constant_prot /= ncalibrated_prot;
  
  constants_file << "BB_prot_cfac = " << endl;
  
  for( int i=0; i<ncell; i++ ){
    char cconst[20];
    
    if( badcells[i] == 0 && solution(i)*E_targ >= min_const && solution(i)*E_targ <= max_const ){ //if the cell was calibrated and the calibration result is between 0.1 and 10, write out the new constant to the file:
      sprintf( cconst, "%15.4g", solution(i) );

      CvChannel->SetBinContent( i+1, E_targ*solution(i) );
      
    } else { //Otherwise, use the average calibration constant:
      if( i < kNrowsSH*kNcolsSH ){
	//sprintf( cconst, "%15.4g", average_constant_prot );
	sprintf( cconst, "%15.4g", oldconstants_sum[i]/oldconstants_nevent[i]); //Use old constant instead of average
	//sprintf( cconst, "%15.4g", -10 ); //Make these failures very obvious
      } else {
	cout << "Error: looping outside of normal geometry for cconst" << endl;
      }
    }
    
    badCvChannel->SetBinContent(i+1,badcells[i]);
    
    if(debug==1) cout << "Constant: " << cconst << ".." << endl;
    
    constants_file << cconst;
    if( i+1 != kNrowsSH*kNcolsSH && i+1 != ncell ) constants_file << ",";
    if( (i+1)%kNcolsSH == 0 ) constants_file << endl;
    if( i == kNrowsSH*kNcolsSH ) constants_file << "Error: ncell > kNrowsSH * kNcolsSH" << endl;
  }
  
  constants_file << endl << "BB_prot_gain_cor = " << endl;

  for(int i=0; i<ncell; i++){

    int r = i/kNcolsSH;
    int c = i%kNcolsSH;

    targetHV[r][c] = gHV[r][c]*pow(cRatio[r][c],(1/gAlphas[r][c]));

    HVOutFile << r*12+c+1 << '\t' << -targetHV[r][c] << endl;
    //HVFile << -targetHV[r][c] << endl;

    char cconst[20];
    sprintf( cconst, "%15.4f", targetHV[r][c] );

    if(debug==1) cout << "Target HV: " << cconst << endl;

    constants_file << cconst;
    if( i+1 != kNrowsSH*kNcolsSH && i+1 != ncell ) constants_file << ",";
    if( (i+1)%kNcolsSH == 0 ) constants_file << endl;
    if( i == kNrowsSH*kNcolsSH ) constants_file << "Error: ncell > kNrowsSH * kNcolsSH" << endl;
  }

  double yErr[ncell] = {0.0};
  //TGraphErrors *ccgraph = new TGraphErrors( ncell, &solution[0], y, xErr, yErr );
  TGraphErrors *ccgraph = new TGraphErrors( ncell, y, &solution[0], yErr, xErr );
  ccgraph->GetXaxis()->SetLimits(0.0,288.0);

  //Write out root file for diagnostics and close
  fout->cd();
  CvChannel->SetTitle( "CConst vs Channel" );
  CvChannel->Write( "CvChannel" );
  CvChannel->Draw( "AP" );

  Phi_p->SetTitle( "Proton Scattering Angle Phi" );
  Phi_p->Write( "Phi_p" );
  Phi_p->Draw( "AP" );

  p_prot->SetTitle( "Proton Momentum" );
  p_prot->Write( "p_prot" );
  p_prot->Draw( "AP" );

  Theta_e->SetTitle( "Electron Scattering Angle Theta" );
  Theta_e->Write( "Theta_e" );
  Theta_e->Draw( "AP" );

  Phi_e->SetTitle( "Electron Scattering Angle Phi" );
  Phi_e->Write( "Phi_e" );
  Phi_e->Draw( "AP" );

  //hcal_col->SetTitle( "Expected HCal Col" );
  //hcal_col->Write( "hcal_col" );
  //hcal_col->Draw( "AP" );

  //hcal_row->SetTitle( "Expected HCal Row" );
  //hcal_row->Write( "hcal_row" );
  //hcal_row->Draw( "AP" );

  BBbesttr->SetTitle( "BB Track with Lowest Chi^2" );
  BBbesttr->Write( "BBbesttr" );
  BBbesttr->Draw( "AP" );

  BBtr_chi2->SetTitle( "Chi^2 among all BB tracks" );
  BBtr_chi2->Write( "BBtr_chi2" );
  BBtr_chi2->Draw( "AP" );

  BBtr_bestchi2->SetTitle( "Chi^2 among best BB tracks" );
  BBtr_bestchi2->Write( "BBtr_bestchi2" );
  BBtr_bestchi2->Draw( "AP" );

  ccgraph->SetTitle("Calibration Constants");
  ccgraph->GetXaxis()->SetTitle("Channel");
  ccgraph->GetXaxis()->SetTitle("GeV/pC");
  ccgraph->SetMarkerStyle(21); //Boxes
  ccgraph->Write("AP");
  ccgraph->Draw("AP");

  NEVvChannel->SetTitle( "CConst vs Channel" );
  NEVvChannel->Write( "NEVvChannel" );
  NEVvChannel->Draw( "AP" );
  
  badCvChannel->SetTitle( "CConst vs Channel" );
  badCvChannel->Write( "badCvChannel" );
  badCvChannel->Draw( "AP" );
  
  CRatiovChannel->SetTitle( "Constant Ratio (new/old) vs Channel" );
  CRatiovChannel->Write( "CRatiovChannel" );
  CRatiovChannel->Draw( "AP" );

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

  for( int r=0; r<kNrowsSH; r++ ){
    for( int c=0; c<kNcolsSH; c++ ){
      gIntSpec[r][c]->GetXaxis()->SetRange(-20,150); //For easier viewing
      gIntSpec[r][c]->Draw();
    }
  }

  fout->Close();

  st->Stop();

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

}
  
