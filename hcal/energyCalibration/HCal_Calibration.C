//SSeeds 10.12.21 - Production - One-pass code combining elements from Andrew Puckett's BigCal calibration script and PDatta/EFuchey BB calibration code used to calibrate HCal energy deposition by PMT module. Tested with Eric Fuchey's digitized g4sbs data and used during detector commissioning in GMn run group.

#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <TChain.h>
#include <TSystem.h>
#include <TStopwatch.h>
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TCut.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TLegend.h"
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

double oldGain[kNcell] = {0.0};
double oldRatio[kNcell] = {0.0};

void HCal_Calibration( const char *configfilename, int run = -1 ){
  
  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  // Start the chain for root files passed with config file

  if( !T ) {
    T = new TChain( "T" );
    T->SetBranchStatus( "*", 0 );
    T->SetBranchStatus( "sbs.hcal.*", 1 );
    T->SetBranchStatus( "bb.tr.*", 1 );
    //T->SetBranchStatus( "bb.tdctrig.tdc", 1 );
    //T->SetBranchStatus( "bb.tdctrig.tdcelemID", 1 );
    //T->SetBranchStatus( "Ndata.bb.tdctrig.tdcelemID", 1 );
    T->SetBranchAddress( "bb.tr.chi2", hcalt::BBtr_chi2 );
    T->SetBranchAddress( "bb.tr.n", hcalt::BBtr_n );
    T->SetBranchAddress( "bb.tr.px", hcalt::BBtr_px );
    T->SetBranchAddress( "bb.tr.py", hcalt::BBtr_py );
    T->SetBranchAddress( "bb.tr.pz", hcalt::BBtr_pz );
    T->SetBranchAddress( "bb.tr.p", hcalt::BBtr_p );
    T->SetBranchAddress( "bb.tr.tg_th", hcalt::BBtr_tg_th );
    T->SetBranchAddress( "bb.tr.tg_ph", hcalt::BBtr_tg_ph );
    T->SetBranchAddress( "sbs.hcal.clus.id", hcalt::cid );
    T->SetBranchAddress( "sbs.hcal.clus.row", hcalt::crow );
    T->SetBranchAddress( "sbs.hcal.clus.col", hcalt::ccol );
    T->SetBranchAddress( "sbs.hcal.clus.e", hcalt::ce );
    T->SetBranchAddress( "sbs.hcal.clus.eblk", hcalt::ceblk ); // Highest energy block
    T->SetBranchAddress( "sbs.hcal.clus.nblk", hcalt::cnblk );
    T->SetBranchAddress( "sbs.hcal.nclus", hcalt::nclus ); 
    T->SetBranchAddress( "sbs.hcal.nblk", hcalt::nblk ); // Total number of blocks in cell
    T->SetBranchAddress( "sbs.hcal.clus_blk.id", hcalt::cblkid ); // kNcell-1 index for each block
    T->SetBranchAddress( "sbs.hcal.clus_blk.e", hcalt::cblke ); // Array of block energies
    //T->SetBranchAddress( "bb.tdctrig.tdcelemID", hcalt::TDCT_id );
    //T->SetBranchAddress( "bb.tdctrig.tdc", hcalt::TDCT_tdc );
    //T->SetBranchAddress( "Ndata.bb.tdctrig.tdcelemID", &hcalt::TDCTndata );

    cout << "Opened up tree with nentries=" << T->GetEntries() << endl;
  }

  //for( int iter=0; iter<2; iter++ ){
  double E_e = 1.92; // Energy of beam (incoming electrons from accelerator)
  int minEventPerCell = 500; // Minimum number of scattered p in cell required to calibrate
  int maxEventPerCell = 4000; // Maximum number of scattered p events to contribute
  double ScaleFac = 1.0;
  double highDelta = 0.1; // Minimum M(i,j)/b(i) factor allowed 
  double HCal_d = 14.5; // Distance to HCal from scattering chamber for comm1
  double HCal_th = 35.0; // Angle that the center of HCal is at  
  double opticsCorr = 1.05; // Correction to magnitude of p_e to account for misaligned optics
  double W_mean = 0.93; // Mean of W at current kinematic
  double W_sig = 0.039; // Width of W at current kinematic

  string inConstPath = "/adaqfs/home/a-onl/sbs/HCal_replay/hcal/energyCalibration/const_v2.txt"; //Where db file sees all 0.00175 GeV/pC
  string constPath = "/adaqfs/home/a-onl/sbs/HCal_replay/hcal/energyCalibration/outfiles/newRatio.txt";
  string ratioPath = "/adaqfs/home/a-onl/sbs/HCal_replay/hcal/energyCalibration/outfiles/ratio.txt";
  string oneBlockPath = "/adaqfs/home/a-onl/sbs/HCal_replay/hcal/energyCalibration/outfiles/oneBlock.txt";

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
      if( skey == "HCal_d" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	HCal_d = sval.Atof();
      }
      if( skey == "HCal_th" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	HCal_th = sval.Atof();
      }
      if( skey == "opticsCorr" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	opticsCorr = sval.Atof();
      }
      if( skey == "W_mean" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        W_mean = sval.Atof();
      }
      if( skey == "W_sig" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        W_sig = sval.Atof();
      }
      if( skey == "ScaleFac" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        ScaleFac = sval.Atof();
      }
    }
    delete tokens;
  }
  
  // Declare outfile
  TFile *fout = new TFile( "eCalOut.root", "RECREATE" );
  
  // Initialize vectors and arrays
  double gOldConst[kNcell];
  double gRatio[kNcell] = {0.0}; //Old ratios from first iteration
  double oneBlock[kNcell] = {0.0};
  double GCoeff[kNcell] = {0.0};
  double GCoeff_oneblock[kNcell] = {0.0};
  double GCoeff_divide[kNcell] = {0.0};
  double A[kNcell] = {0.0}; // Array to keep track of ADC values per cell. Memory is reset on each event.

  // Initialize histograms	
  TH1D *hDeltaE = new TH1D( "hDeltaE","1.0-Eclus/p_rec", 100, -1.5, 1.5 );
  TH1D *hClusE = new TH1D( "hClusE","Best Cluster Energy", 100, 0.0, 2.0);
  TH2D *hPAngleCorr = new TH2D( "hPAngCorr","Track p vs Track ang", 100, 30, 60, 100, 0.4, 1.2 );
  TH1D *hW = new TH1D( "W", "W", 250, 0.7, 1.5 );
  hW->GetXaxis()->SetTitle( "GeV" );
  TH1D *hQ2 = new TH1D( "Q2", "Q2", 250, 0.5, 3.0 );
  hQ2->GetXaxis()->SetTitle( "GeV" );
  TH1D *hE_ep = new TH1D( "Scattered Electron Energy","E_ep", 500, 0.0, E_e*1.5 ); 
  hE_ep->GetXaxis()->SetTitle( "GeV" );
  TH1D *hE_pp = new TH1D( "Scattered Proton Energy", "E_pp", 500, 0.0, E_e*1.5 );
  hE_pp->GetXaxis()->SetTitle( "GeV" );
  TH1D *hKE_p = new TH1D( "Scattered Proton Kinetic Energy", "KE_pp", 500, 0.0, E_e*1.5 );
  hKE_p->GetXaxis()->SetTitle( "GeV" );

  TH2D *hADC = new TH2D( "hADC", "HCal Int_ADC Spectra: W Cut", 288, 0, 288, 100., 0., 1. );
  TH2D *hADC_amp = new TH2D( "hADC_amp", "HCal ADC_amp Spectra: W Cut", 288, 0, 288, 100., 0., 10. );

  Long64_t Nevents = T->GetEntries();
  cout << "Opened up tree with nentries: " << Nevents << ".." << endl;

  // Read in previous constants for PMTs
  ifstream inConstFile( inConstPath );
  if( !inConstFile ){
    cerr << "No input constant file present -> setFiles/const.txt expected." << endl;
    return 0;
  }

  cout << "Getting previous calibration constants.." << endl;
  int n1=0;
  double d1=0;
  string line;
  
  while( getline( inConstFile, line ) ){
    if( line.at( 0 )=='#' ) {
      continue;
    }
    istringstream iss( line );
    while( iss >> d1 ){
      gOldConst[n1] = d1;
      cout << "const" << n1 << " = " << d1 << "." << endl;
      n1++;
    }
  }
  cout << "All parameters loaded." << endl;

  // Need to get electron energy and W for cuts on elastic events
  long mevent = 0; 
  cout << "Main loop over all data commencing.." << endl;
    
  TMatrixD Ma(kNcell,kNcell);
  TVectorD ba(kNcell);
  TVectorD bb(kNcell);
  int NEV[kNcell] = {0};
  int NEV_oneblock[kNcell] = {0};

  while( T->GetEntry( mevent++ ) ){ 

    if( mevent%1000 == 0 ) cout << "Main loop: " << mevent << "/" << Nevents << endl;
    
    memset(A, 0, kNcell*sizeof(double));

    // Sort all tracks by lowest Chi^2 (highest confidence track)
    int track_tot = (int)hcalt::BBtr_n[0];
    int track = 0;
    double min = 1000.0; //Set arbitrarily high chi^2 to minimize with loop over tracks
    for( int elem=0; elem<track_tot; elem++ ){            
      if( hcalt::BBtr_chi2[elem] < min )
	track = elem;      
    }

    double E_ep = sqrt( pow(M_e,2) + pow(hcalt::BBtr_p[track],2) ); // Obtain the scattered electron energy
    hE_ep->Fill( E_ep ); // Fill histogram

    double p_ep = opticsCorr*hcalt::BBtr_p[track]; // Obtain the magnitude of scattered electron momentum with correction factor
    double Q2 = 2*E_e*E_ep*( 1-(hcalt::BBtr_pz[track]/p_ep) ); // Obtain Q2 from beam energy, outgoing electron energy, and momenta
    hQ2->Fill( Q2 ); // Fill histogram

    double nu = E_e-E_ep; // Obtain energy transfer
    double W2 = pow( M_p,2 )+2*M_p*nu-Q2; // Obtain W2 from Q2 and nu
    double W = 0.0; // Should be Mp for elastic events. Will use to estimate correction factor while GEM tracking ongoing
    if( W2>0.0 ){
      W = sqrt( W2 ); // Obtain W for real events by throwing out negative W2 solns
      hW->Fill( W );
    }
    
    // Use the electron kinematics to predict the proton momentum assuming elastic scattering on free proton at rest:
    double E_pp = nu+M_p; // Get energy of the proton
    hE_pp->Fill( E_pp ); // Fill histogram

    double p_p = sqrt( pow( E_pp,2 )-W2 ); // Magnitude of the scattered proton momentum

    //double KE_p = pow( p_p,2 )/(2*M_p);
    double KE_p = nu; //For elastics
    hKE_p->Fill( KE_p );

    //Electron polar angle
    double phi_p = TMath::ACos( ( E_e-hcalt::BBtr_pz[track])/p_p ) * 180.0 / PI;    

    //Will need to wait until next mass replay
    /*
    //Cut on BBCal and HCal trigger coincidence
    double bbcal_time=0., hcal_time=0.;
    for(int ihit=0; ihit<hcalt::TDCTndata; ihit++){
      if(hcalt::TDCT_id[ihit]==5) bbcal_time=hcalt::TDCT_tdc[ihit];
      if(hcalt::TDCT_id[ihit]==0) hcal_time=hcalt::TDCT_tdc[ihit];
    }
    double diff = hcal_time - bbcal_time; 
    */

    //if( fabs(W-W_mean)<W_sig && fabs(diff-510.)<20 )

    if( fabs(W-W_mean)<W_sig ){
      

      // Loop over all the clusters first to get the one with the highest energy
      double Etop = -10.0;
      int cluster = 0;
      for(Int_t cl = 0; cl<hcalt::nclus[0]; cl++){
      if(hcalt::ce[cl]>Etop){
        Etop = hcalt::ce[cl];
        cluster = cl;
      	}
      }
      
      //if(cluster>-1)cout << cluster << endl;

      // Reject events where the primary block in the primary cluster is on the edge of the acceptance and throw out events where no clusters exist
      if( hcalt::crow[cluster]==0 ||
	  hcalt::crow[cluster]==23 ||
	  hcalt::ccol[cluster]==0 ||
	  hcalt::ccol[cluster]==11 )
	continue;

      double clusE = 0.0;

      // Get energies with simplest scheme from clusters only
      for( int blk = 0; blk<hcalt::cnblk[cluster]; blk++ ){
	int blkid = int(hcalt::cblkid[blk])-1; //-1 necessary since sbs.hcal.clus_blk.id ranges from 1 to 288 (different than just about everything else)
       
	clusE += hcalt::cblke[blk];
	A[blkid] += hcalt::cblke[blk];

	hADC->Fill( blkid, hcalt::cblke[blk] );
	
	// Simple estimation of the coefficients assuming 100% energy deposition in one block. Will check against the more sophisticated version with chi^2 reduction.
	if(hcalt::nblk[cluster]==1)
	  NEV_oneblock[blkid]++;
	
	NEV[blkid]++;
      }

      // Some of PDatta's favorite histograms..
      hDeltaE->Fill( 1.0-(clusE/KE_p) );
      hClusE->Fill( clusE );
      hPAngleCorr->Fill( phi_p, p_ep );
      
      // Build the matrix as simply as possible
      for(int icol = 0; icol<kNcell; icol++){
	ba(icol)+= A[icol];
	if(hcalt::nblk[cluster]==1){
	  bb(icol)+= A[icol];
	  oneBlock[icol]+= A[icol]*A[icol]/(KE_p*0.0795);
	}
	for(int irow = 0; irow<kNcell; irow++){
	  Ma(icol,irow)+= A[icol]*A[irow]/(KE_p*0.0795); // Must add sampling fraction of the detector (0.06588 GeV/0.8286 GeV) = 0.0795 = 7.95% -> (EFuchey MC E_dep per proton) / (fit to data KE_p)
	  	  
	}
      } 
    }
  }

  // Reject the bad cells and normalize the oneblock check
  int badcell[kNcell];
  int badcell_oneblock[kNcell];
  double y[kNcell] = {0.0}; // For easy TGraphErrors build
  
  for(int i=0; i<kNcell; i++){
    badcell[i] = 0;
    y[i] = i;

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

  for(int i=0; i<kNcell; i++){
    badcell_oneblock[i] = 0;
    if( NEV_oneblock[i] < minEventPerCell || oneBlock[i] < 0.1*bb[i] ){
      
      cout << NEV_oneblock[i] << " " << oneBlock[i] << " " << bb[i] << endl;

      bb(i) = 1.0;
      oneBlock[i] = 1.0;
      badcell_oneblock[i] = 1;
    }
  }
  
  // Invert the matrix, solve for ratios
  TMatrixD M_inv = Ma.Invert();
  TVectorD Coeff = M_inv*ba; // Stayes unmodified for reference
  double oneBlockCoeff[kNcell];

  for( int i=0; i<kNcell; i++ ){
    if(badcell_oneblock[i]==0){
      oneBlockCoeff[i] = gOldConst[i]*(bb[i]/oneBlock[i]);

    }else{
      oneBlockCoeff[i] = gOldConst[i];
    }
  }

  

  for(int i=0; i<kNcell; i++){
    if(badcell[i]==0){
      GCoeff[i]=gOldConst[i]*Coeff[i]; // The new gain coefficient is the old coefficient multiplied by the gain factor that we just solved for. We will take these coefficents as the input for another set of data to estimate the error per channel.
      GCoeff_divide[i]=GCoeff[i]/oneBlockCoeff[i];

      //cout << GCoeff[i] << " / " << oneBlockCoeff[i] << " = " << GCoeff_divide[i] << endl; 

    }else{
      GCoeff[i]=gOldConst[i]; // If the cell is bad, use the old coefficient
      GCoeff_divide[i]=-1.0;
    }
  }

  double yErr[kNcell] = {0.0};
  TGraphErrors *ccgraph = new TGraphErrors( kNcell, y, GCoeff, yErr, yErr ); 
  ccgraph->GetXaxis()->SetLimits(0.0,288.0);  
  ccgraph->GetYaxis()->SetLimits(0.0,0.025);
  ccgraph->SetTitle("Calibration Coefficients");
  ccgraph->GetXaxis()->SetTitle("Channel");
  ccgraph->GetXaxis()->SetTitle("Unitless");
  ccgraph->SetMarkerStyle(20); //idx 20 Circles, idx 21 Boxes
  ccgraph->Write("constants");

  TGraphErrors *ccgraph_oneBlock = new TGraphErrors( kNcell, y, oneBlockCoeff, yErr, yErr ); 
  ccgraph_oneBlock->GetXaxis()->SetLimits(0.0,288.0);  
  ccgraph_oneBlock->GetYaxis()->SetLimits(0.0,0.025);
  ccgraph_oneBlock->SetTitle("Calibration Coefficients One Block");
  ccgraph_oneBlock->GetXaxis()->SetTitle("Channel");
  ccgraph_oneBlock->GetXaxis()->SetTitle("Unitless");
  ccgraph_oneBlock->SetMarkerStyle(21);
  ccgraph_oneBlock->Write("constants_oneblock");

  TGraphErrors *ccgraph_divide = new TGraphErrors( kNcell, y, GCoeff_divide, yErr, yErr ); 
  ccgraph_divide->GetXaxis()->SetLimits(0.0,288.0);  
  ccgraph_divide->SetTitle("Calibration Coefficients / OneBlock Coeff");
  ccgraph_divide->GetXaxis()->SetTitle("Channel");
  ccgraph_divide->GetXaxis()->SetTitle("Unitless");
  ccgraph_divide->SetMarkerStyle(21);
  ccgraph_divide->Write("constants_divide");

  // Write out diagnostic histos and cout coeff
  fout->Write();

  cout << "Gain Coefficients: " << endl << endl;

  // Print to console
  int cell = 0;
  for( int r = 0; r<kNrows; r++){
    for( int c = 0; c<kNcols; c++){
      if( GCoeff[cell]==1 ) GCoeff[cell] = 0.00175;
      cout << GCoeff[cell] << "  ";
      cell++;
    }
    cout << endl;
  }
  cell = 0;

  for( int i=0; i<kNcell; i++){
    GCoeff[i] = GCoeff[i] - 0.00175;
  }

  TGraphErrors *ccgraph_Cdiff = new TGraphErrors( kNcell, y, GCoeff, yErr, yErr ); 
  ccgraph_Cdiff->GetXaxis()->SetLimits(0.0,288.0);  
  ccgraph_Cdiff->GetYaxis()->SetLimits(0.0,0.025);
  ccgraph_Cdiff->SetTitle("Cosmic Delta Calibration Coefficients");
  ccgraph_Cdiff->GetXaxis()->SetTitle("Channel");
  ccgraph_Cdiff->GetXaxis()->SetTitle("Unitless");
  ccgraph_Cdiff->SetMarkerStyle(20); //idx 20 Circles, idx 21 Boxes
  ccgraph_Cdiff->Write("constants_diff");

  cout << endl << "Gain Coefficients diff from cosmic value: " << endl;

  for( int r = 0; r<kNrows; r++){
    for( int c = 0; c<kNcols; c++){
      cout << GCoeff[cell] << "  ";
      cell++;
    }
    cout << endl;
  }
  cell = 0;

  cout << endl << "Number of events available for calibration: " << endl << endl;

  for( int r = 0; r<kNrows; r++){
    for( int c = 0; c<kNcols; c++){
      cout << NEV[cell] << "  ";
      cell++;
    }
    cout << endl;
  }
  cell = 0;

  cout << endl << "One Block:" << endl;

  for( int r = 0; r<kNrows; r++){
    for( int c = 0; c<kNcols; c++){
      cout << oneBlockCoeff[cell] << "  ";
      cell++;
    }
    cout << endl;
  }

  cout << endl << endl;

  // Declare outfiles
  ofstream ratio;
  ofstream GainCoeff;
  ofstream GainCoeff_oneblock;

  ratio.open( ratioPath );
  ratio << "HCal_ratio = " << endl;

  for( int i=0; i<kNcell; i++ ){   
    ratio << Coeff[i] << endl;
  }

  ratio.close();

  GainCoeff.open( constPath );
  GainCoeff << "HCal_gainCoeff = " << endl;

  cell = 0;
  for( int r=0; r<kNrows; r++ ){
    for( int c=0; c<kNcols; c++ ){
      GainCoeff << GCoeff[cell] << "  ";
      cell++;
    }
    GainCoeff << endl;
  }
  cell = 0;

  GainCoeff.close();

  GainCoeff_oneblock.open( oneBlockPath );
  GainCoeff_oneblock << "HCal_gainCoeff_oneBlock = " << endl;

  for( int r=0; r<kNrows; r++ ){
    for( int c=0; c<kNcols; c++ ){
      GainCoeff_oneblock << oneBlockCoeff[cell] << "  ";
      cell++;
    }
    GainCoeff_oneblock << endl;
  }
  cell = 0;

  GainCoeff_oneblock.close();

  cout << "Calibration complete and constants written to file." << endl;

  st->Stop();

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

}

//Sample setup file, comment with #
////////////////////////////////////////////////////////
//setup_HCal_Calibration.txt
////////////////////////////////////////////////////////
//#LH2 full field
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11263*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11265*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11267*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11268*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11252*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11251*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11250*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11249*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11248*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11247*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11246*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11242*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11244*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11245*
//#LH2 quarter field
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11301*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11302*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11303*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11304*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11305*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11306*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11307*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11308*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11309*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11310*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11311*
//#LH2 inverted quarter field
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11276*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11277*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11279*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11280*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11281*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11282*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11283*
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11284*
//#LH2 Zero Field 
//#/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11332*
//endlist
//endcut
//E_e 1.92
//HCal_d 14.5
//HCal_th 35.0
//opticsCorr 1.05
//W_mean 0.93
//W_sig 0.1
//ScaleFac 1.0
/////////////////////////////////////////////////////////
