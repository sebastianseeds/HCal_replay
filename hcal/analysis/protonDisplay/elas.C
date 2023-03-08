//SSeeds 10.19.22 Further simplification to proton spot code. This should just produce elastic yields and show the delta position plots
//SSeeds 10.20.22 Reconstruction of energy calibration code to explain apparent ADC gain discrepancy

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

const Int_t kNcell = 288; // Total number of HCal modules
const Int_t kNrows = 24; // Total number of HCal rows
//const Int_t kNrows = 12; // Total number of HCal rows for top half analysis kin2 GEn
const Int_t kNcols = 12; // Total number of HCal columns
const Int_t kNtrack = 100; // Reasonable max on tracks in BigBite
const Int_t kNtdc = 100; // Reasonable max on tdc trigger channels
const Double_t Xi = -2.20; // Distance from beam center to top of HCal in m
const Double_t Xf = 1.47; // Distance from beam center to bottom of HCal in m
const Double_t Yi = -0.853; // Distance from beam center to opposite-beam side of HCal in m
const Double_t Yf = 0.853; // Distance from beam center to beam side of HCal in m
const Double_t M_p = 0.938272;
const Double_t PI = TMath::Pi();
const Int_t minEventPerCell = 30;
const Double_t tdiffmax = 20;
const Double_t hcalheight = 0.365; //m The height of the center of HCAL above beam
//const Double_t sampFrac = 0.0795; // sampling frac (0.06588 GeV/0.8286 GeV) = 0.0795 = 7.95% -> (MC E_dep per proton) / (fit to data KE_p)
const Double_t sampFrac = 0.077;  //Re-evaluated with MC GEn settings using second to outermost shower column for kin2

void elas( const char *configfilename="k2_elas.cfg", int run = -1 ){

  // Start the chain for root files passed with config file
  TChain *C = new TChain("T");
  
  Double_t E_e = 1.92; // Energy of beam (incoming electrons from accelerator)
  Int_t SBS_field = 0; // Strength in percent of 2100A of SBS magnet
  Double_t W2_mean = M_p; // With perfect optics, this should be true. Will read in until then by fitting W distribution on each run
  Double_t W2_sigma = 0.03; // Reasonable value by default, but will read in from W distribution
  Double_t HCal_d = 17; // HCal distance from face to target chamber
  Double_t HCal_th = 34.7*(TMath::Pi()/180.); // HCal angle theta from downstream beamline
  Double_t tdiff = 510; // Time difference between BBCal trigger and HCal overlapping regions trigger in ns

  // Reading config file
  ifstream configfile(configfilename);
  TString currentline;
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      C->Add(currentline);
      cout << "Loaded file: " << currentline <<  endl;
    }    
  }
  TCut globalcut = "";
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endcut") ){
    if( !currentline.BeginsWith("#") ){
      globalcut += currentline;
      cout << "Applying the following cut to all data: " << globalcut << endl;
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
	cout << "Beam energy (GeV): " << E_e << endl;
      }
      if( skey == "SBS_field" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	SBS_field = sval.Atoi();
	cout << "SBS Magnetic field (%): " << SBS_field << endl;
      }
      if( skey == "W2_mean" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	W2_mean = sval.Atof();
	cout << "W2 mean (GeV): " << W2_mean << endl;
      }
      if( skey == "W2_sigma" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	W2_sigma = sval.Atof();
	cout << "W2 sigma: " << W2_sigma << endl;
      }
      if( skey == "HCal_d" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        HCal_d = sval.Atof();
	cout << "HCal distance (m): " << HCal_d << endl;
      }
      if( skey == "HCal_th" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        HCal_th = sval.Atof();
	cout << "HCal theta (deg): " << HCal_th << endl;
	HCal_th *= (TMath::Pi()/180.); //convert to radians
      }
      if( skey == "tdiff" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        tdiff = sval.Atof();
	cout << "Time difference BBCal/HCal (ns): " << tdiff << endl;
      }
      delete tokens;
    }
  }


  // Path for previous ADC gain constants
  //string inConstPath = "/adaqfs/home/a-onl/sbs/sbs_devel/SBS-replay/DB/db_sbs.hcal.dat";
  //string inConstPath = "/adaqfs/home/a-onl/sbs/HCal_replay/replay/db_sbs.hcal.dat";
  string inConstPath = "/w/halla-scshelf2102/sbs/seeds/SBS-replay/DB/db_sbs.hcal.dat";
  

  // Reading ADC gain parameters from database
  Double_t gOldConst[kNcell];
  
  cout << "Loading previous gain coefficients from file: " << inConstPath << ".." << endl;
  ifstream inConstFile( inConstPath );
  if( !inConstFile ){
    cerr << endl << "ERROR: No input constant file present -> path to db_sbs.hcal.dat expected." << endl;
    return 0;
  }

  Int_t n1=0;
  Double_t d1=0;
  string line;
  bool skip_line = true;
  bool skip_one_line = true;
  bool skip_first_instance = true;
  //bool skip_second_instance = true;
  
  while( getline( inConstFile, line ) ){

    if( n1==( kNcell ) ) break;

    TString Tline = (TString)line;

    if( Tline.BeginsWith("sbs.hcal.adc.gain") && skip_first_instance==true ){
      skip_first_instance = false;
      continue;
    }
    /*
    if( Tline.BeginsWith("sbs.hcal.adc.gain") && skip_first_instance==false && skip_second_instance==true ){
      skip_second_instance = false;
      continue;
    }
    */
    if( Tline.BeginsWith("sbs.hcal.adc.gain") && skip_line==true ) skip_line = false;

    if( skip_line==false && skip_one_line==true ){
      skip_one_line = false;
      continue;
    }

    if( skip_line==false && skip_one_line==false ){
      istringstream iss( line );
      while( iss >> d1 ){
	gOldConst[n1] = d1;
	n1++;
      }
    }
  }
  
  cout << endl << endl << "Old ADC gain params: " << endl;

  for( Int_t r=0; r<kNrows; r++){
    for( Int_t c=0; c<kNcols; c++){
      Int_t i = r*kNcols+c;
      cout << gOldConst[i] << " ";
    }
    cout << endl;
  }
  
  cout << endl << endl << "Setup parameters loaded." << endl;

  //Declare vars
  Double_t atime[kNcell], row[kNcell], col[kNcell], tdctime[kNcell], cblkid[kNcell], cblke[kNcell];
  Double_t nblk, nclus, SHnclus, PSnclus, hcalx, hcaly, hcale;
  UInt_t TBits;

  Double_t BBtr_p[kNtrack], BBtr_px[kNtrack], BBtr_py[kNtrack], BBtr_pz[kNtrack];
  Double_t BBtr_vz[kNtrack];
  Double_t BBtr_n, BBps_x, BBps_y, BBps_e, BBsh_x, BBsh_y, BBsh_e;

  Double_t TDCT_id[kNtdc], TDCT_tdc[kNtdc], hodo_tmean[kNtdc]; 
  Int_t TDCTndata;

  //Setup leaves
  C->SetMakeClass(1); //Allows for viewing of general detector params from Event_Branch
  C->SetBranchStatus( "*", 0 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.atime", 1 );  
  C->SetBranchStatus( "sbs.hcal.clus_blk.row", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.col", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.tdctime", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.id", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.e", 1 );
  C->SetBranchStatus( "sbs.hcal.x", 1 );
  C->SetBranchStatus( "sbs.hcal.y", 1 );
  C->SetBranchStatus( "sbs.hcal.e", 1 );
  C->SetBranchStatus( "sbs.hcal.nblk", 1 );
  C->SetBranchStatus( "sbs.hcal.nclus", 1 );
  C->SetBranchStatus( "bb.sh.nclus", 1 );
  C->SetBranchStatus( "bb.ps.nclus", 1 );
  C->SetBranchStatus( "bb.tr.n", 1 );
  C->SetBranchStatus( "bb.tr.px", 1 );
  C->SetBranchStatus( "bb.tr.py", 1 );
  C->SetBranchStatus( "bb.tr.pz", 1 );
  C->SetBranchStatus( "bb.tr.p", 1 );
  C->SetBranchStatus( "fEvtHdr.fTrigBits", 1 );
  C->SetBranchStatus( "bb.tr.vz", 1 );
  C->SetBranchStatus( "bb.ps.e", 1 );
  C->SetBranchStatus( "bb.ps.x", 1 );
  C->SetBranchStatus( "bb.ps.y", 1 );
  C->SetBranchStatus( "bb.sh.e", 1 );
  C->SetBranchStatus( "bb.sh.x", 1 );
  C->SetBranchStatus( "bb.sh.y", 1 );
  C->SetBranchStatus( "g.trigbits", 1 );
  C->SetBranchStatus( "bb.tdctrig.tdc", 1 );
  C->SetBranchStatus( "bb.tdctrig.tdcelemID", 1 );
  C->SetBranchStatus( "Ndata.bb.tdctrig.tdcelemID", 1 );

  C->SetBranchAddress( "sbs.hcal.clus_blk.atime", atime );
  C->SetBranchAddress( "sbs.hcal.clus_blk.row", row );
  C->SetBranchAddress( "sbs.hcal.clus_blk.col", col );
  C->SetBranchAddress( "sbs.hcal.clus_blk.tdctime", tdctime );
  C->SetBranchAddress( "sbs.hcal.clus_blk.id", cblkid );
  C->SetBranchAddress( "sbs.hcal.clus_blk.e", cblke );
  C->SetBranchAddress( "sbs.hcal.x", &hcalx );
  C->SetBranchAddress( "sbs.hcal.y", &hcaly );
  C->SetBranchAddress( "sbs.hcal.e", &hcale );
  C->SetBranchAddress( "sbs.hcal.nblk", &nblk );
  C->SetBranchAddress( "sbs.hcal.nclus", &nclus );
  C->SetBranchAddress( "bb.sh.nclus", &SHnclus );
  C->SetBranchAddress( "bb.ps.nclus", &PSnclus );
  C->SetBranchAddress( "bb.tr.n", &BBtr_n );
  C->SetBranchAddress( "bb.tr.px", BBtr_px );
  C->SetBranchAddress( "bb.tr.py", BBtr_py );
  C->SetBranchAddress( "bb.tr.pz", BBtr_pz );
  C->SetBranchAddress( "bb.tr.vz", BBtr_vz );
  C->SetBranchAddress( "bb.tr.p", BBtr_p );
  C->SetBranchAddress( "bb.ps.e", &BBps_e );
  C->SetBranchAddress( "bb.ps.x", &BBps_x );
  C->SetBranchAddress( "bb.ps.y", &BBps_y );
  C->SetBranchAddress( "bb.sh.e", &BBsh_e );
  C->SetBranchAddress( "bb.sh.x", &BBsh_x );
  C->SetBranchAddress( "bb.sh.y", &BBsh_y );
  C->SetBranchAddress( "fEvtHdr.fTrigBits", &TBits ); //For GEn, TBits==5 is coin
  C->SetBranchAddress( "bb.tdctrig.tdcelemID", TDCT_id );
  C->SetBranchAddress( "bb.tdctrig.tdc", TDCT_tdc );
  C->SetBranchAddress( "Ndata.bb.tdctrig.tdcelemID", &TDCTndata );

  // Declare outfile
  TFile *fout = new TFile( "k2_precal.root", "RECREATE" );

  TEventList *elist = new TEventList("elist","Elastic Event List");
  C->Draw(">>elist",globalcut);

  // Initialize vectors and arrays
  //Double_t A[kNcell] = {0.0}; // Array to keep track of ADC values per cell. Memory is reset on each event.
  Double_t GCoeff[kNcell] = {0.0};
  Double_t GCoeffDiff[kNcell] = {0.0};

  // Initialize histograms
  TH1D *h_atime = new TH1D( "atime", "HCal ADC Time, All Channels; ns", 160, 0, 160 );
  TH1D *h_coeff = new TH1D( "coeff", "HCal New Raw Coefficients, All Channels; ratio", 100, 0, 4 );
  TH1D *h_E = new TH1D( "E", "HCal Cluster E, All Channels; GeV", 100, 0, 0.5 );
  TH1D *h_E_cut = new TH1D( "E_cut", "HCal Cluster E All Cuts, All Channels; GeV", 100, 0, 0.5 );
  TH1D *h_E_active = new TH1D( "E_active", "HCal Cluster E, Active Area; GeV", 100, 0, 0.5 );
  TH1D *h_E_top = new TH1D( "E_top", "HCal Cluster E, First 12 Rows; GeV", 100, 0, 0.5 );
  TH1D *h_E_topH = new TH1D( "E_topH", "HCal Cluster E, First 12 Rows; GeV", 200, 0, 0.5 );
  TH1D *h_E_topL = new TH1D( "E_topL", "HCal Cluster E, First 12 Rows; GeV", 50, 0, 0.5 );
  TH1D *h_E_exp = new TH1D( "E_exp", "Expected Energy Dep in HCal; GeV", 100, 0, 0.2 );
  TH1D *h_vert = new TH1D( "vert", "Vertex Position; m", 200, -1.0, 1.0 );
  TH2D *h_EvCh = new TH2D( "EvCh", "HCal Cluster E Single Block Clusters; channel, GeV", kNcell, 0, kNcell, 100, 0, 0.5 );
  TH2D *h_CvCh = new TH2D( "CvCh", "HCal Coeff Single Block Clusters; channel, GeV", kNcell, 0, kNcell, 200, 0, 1.0 );
  TH1D *h_W2 = new TH1D( "W2", "W2 No Cuts; GeV", 200, 0.0, 4.0 );
  TH1D *h_W2recon = new TH1D( "W2recon", "W2 No Cuts Reconstructed; GeV", 200, 0.0, 4.0 );
  //TH2D *hdxdy = new TH2D("dxdy",";y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)", 50, -2, 2, 100, -4.0, 0.0 );
  TH2D *hdxdy = new TH2D("dxdy",";y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)", 50, -2, 2, 100, -4.0, 2.0 );
  TH2D *hdxdy_all = new TH2D("dxdy_all",";y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)",125,-2,2,125,-4,6);
  //TH1D *hdx = new TH1D( "dx", "HCal dx; m", 200, -4.0, 0.0 );
  TH1D *hdx = new TH1D( "dx", "HCal dx; m", 200, -4.0, 2.0 );
  TH1D *hdy = new TH1D( "dy", "HCal dy; m", 100, -1.2, 1.2 );
  TH1D *hKE_p = new TH1D( "KE_p", "Scattered Proton Kinetic Energy", 500, 0.0, 5.0 );
  TH2D *hdxVE = new TH2D("dxVE",";x_{HCAL}-x_{expect} (m); GeV", 100, -4.0, 2.0, 100, 0, 0.5 );
  TH1D *hKElow = new TH1D( "KElow", "Lowest Elastic E Sampled in HCal (GeV)", 500, 0.0, 0.2 );
  TH1D *hDiff = new TH1D( "hDiff","HCal time - BBCal time (ns)", 1300, -500, 800 );
  TH1D *hSampFrac = new TH1D( "hSampFrac","HCal Cluster E / Expected KE", 400, 0., 0.05 );
  TH2D *hrowcol = new TH2D( "hrowcol", "HCal Block Position Elastics, HCal; Col; Row", kNcols, 0, kNcols, kNrows, -kNrows, 0 );
  TH1D *hY = new TH1D( "Y", "HCal Y; m", 100, -1.2, 1.2 );


  //Vars to track script
  Long64_t Nevents = elist->GetN();
  Int_t eyield = 0;

  for( Int_t r=0; r<kNrows; r++){
    for ( Int_t c=0; c<kNcols; c++){
      hrowcol->Fill( (c+1), -(r+1) );
    }
  }
  
  //Beam and detector vars
  TLorentzVector Pbeam( 0, 0, E_e, E_e );
  TLorentzVector Ptarg( 0, 0, 0, 0.5*(0.93827+0.93956) );
  TVector3 hcal_origin( -HCal_d*sin(HCal_th), 0, HCal_d*cos(HCal_th) );
  TVector3 hcal_zaxis = hcal_origin.Unit();
  TVector3 hcal_xaxis(0,-1,0);
  TVector3 hcal_yaxis = hcal_zaxis.Cross( hcal_xaxis ).Unit();

  //Declare matrices for chi-square min calibration scheme and keep track of calibrated events
  TMatrixD Ma(kNcell,kNcell);
  TVectorD ba(kNcell);
  TVectorD bb(kNcell);
  Int_t NEV[kNcell] = {0};

  cout << endl << "Proceeding to loop over all events in chain.." << endl;

  for(Long64_t nevent = 0; nevent<Nevents; nevent++){

    if( nevent%10 == 0 ) cout << "Loop: " << nevent << "/" << Nevents << ". Elastic yield: " << eyield << ". \r";
    cout.flush();

    C->GetEntry( elist->GetEntry( nevent ) );

    Double_t A[kNcell] = {0.0}; //Declare double to hold ADC values

    TLorentzVector kprime( BBtr_px[0], BBtr_py[0], BBtr_pz[0], BBtr_p[0] );
    TLorentzVector q = Pbeam - kprime;
    TVector3 vertex( 0, 0, BBtr_vz[0] );
    TVector3 qunit = q.Vect().Unit();

    Double_t sintersect = (hcal_origin-vertex).Dot( hcal_zaxis )/qunit.Dot( hcal_zaxis );
    TVector3 hcal_intersect = vertex + sintersect * qunit;    
    Double_t xexpect = hcal_intersect.Dot( hcal_xaxis );
    Double_t yexpect = hcal_intersect.Dot( hcal_yaxis );
    
    Double_t W2recon = (Ptarg + q).M2();
    Double_t E_ep = BBtr_p[0]; // Obtain the scattered electron energy, neglect mass e
    Double_t p_ep = BBtr_p[0]; // Obtain the magnitude of scattered electron momentum
    Double_t Q2 = 2*E_e*E_ep*( 1-(BBtr_pz[0]/p_ep) );
    Double_t nu = E_e-E_ep; // Obtain energy transfer
    Double_t W2 = pow( M_p,2 )+2*M_p*nu-Q2; // Obtain W2 from Q2 and nu
    Double_t eth = acos( BBtr_pz[0]/BBtr_p[0] );
    Double_t eph = atan2( BBtr_py[0], BBtr_px[0] );
    Double_t phinucleon = eph + TMath::Pi(); //assume coplanarity
    Double_t thetanucleon = acos( (E_e - BBtr_pz[0])/p_ep ); //use elastic constraint on nucleon kinematics
    TVector3 pNhat( sin(thetanucleon)*cos(phinucleon),sin(thetanucleon)*sin(phinucleon),cos(thetanucleon));
    Double_t KE_p = nu; //For elastics
    hKE_p->Fill( KE_p );   

    //Cut on BBCal and HCal trigger coincidence
    Double_t bbcal_time=0., hcal_time=0., coin_time=0., rf_time=0.;
    bool cointrig = false;
    for(Int_t ihit=0; ihit<TDCTndata; ihit++){
      if(TDCT_id[ihit]==5) bbcal_time=TDCT_tdc[ihit];
      if(TDCT_id[ihit]==0) hcal_time=TDCT_tdc[ihit];
      if(TDCT_id[ihit]==1) {
	coin_time=TDCT_tdc[ihit];
	cointrig=true;
      }
      if(TDCT_id[ihit]==4) rf_time=TDCT_tdc[ihit];
    }
    Double_t diff = hcal_time - bbcal_time; 
    hDiff->Fill( diff );

    hdxdy_all->Fill( hcaly - yexpect, hcalx - xexpect );
    h_vert->Fill(BBtr_vz[0]);

    h_W2->Fill(W2);
    h_W2recon->Fill(W2recon);

    ///////////////
    //PRIMARYCUTS//
    //Cut on vertex inside target
    if( abs(BBtr_vz[0])>0.27 ) continue;
    
    //Cut on coin
    if( TBits==1 ) continue; //Presently removes almost everything - will need to look into this further.
    //Atime cut hcal
    //if( atime[0]>60. && atime[0]<90. )
    //Cut on  W2
    if( abs(W2-W2_mean)>W2_sigma ) continue;
    //Cut on BBCal HCal trig diff
    if( abs(diff-tdiff)>tdiffmax ) continue;
    //Bigbite track / HCal angular correlation cut
    //ENDCUTS//
    ///////////
    //eyield++;
    hKElow->Fill( KE_p*sampFrac );


    Double_t SFrac = hcale/KE_p;
    hSampFrac->Fill( SFrac );


    h_E->Fill(hcale);
    hdxVE->Fill(hcalx - xexpect,hcale);

    if( row[0]>0 && row[0]<23 && col[0]>0 && col[0]<11 ){
      h_E_active->Fill(hcale);
      
    }

    if( row[0]<13 ) {
      h_E_top->Fill(hcale);
      h_E_topH->Fill(hcale);
      h_E_topL->Fill(hcale);
    }
    hdxdy->Fill( hcaly-yexpect, hcalx-xexpect );
    hdx->Fill( hcalx-xexpect );
    hdy->Fill( hcaly-yexpect );
    h_atime->Fill( atime[0] );
    hY->Fill( hcaly );

    if( nblk==1 ){
      h_EvCh->Fill( cblkid[0], hcale );
    }
    for( Int_t b=0; b<nblk; b++){      
      //hrowcol->Fill( (Double_t)col[b]+1., -(Double_t)row[b]+1.);
      hrowcol->Fill( (Double_t)col[b], -(Double_t)row[b]-1); 
    }

    //Construct matrix

    ///////////////
    //SECONDARYCUTS

    //Reject events where the primary block in the primary cluster is on the edge of the acceptance
    if( row[0]==0 ||
	row[0]==23 ||
	col[0]==0 ||
	col[0]==11 )
      continue;
    eyield++;

    //ENDCUTS
    /////////

    //if( nblk!=1 ) continue;

    Double_t E_exp = KE_p*sampFrac;
    //Double_t E_exp = (KE_p+M_p)*sampFrac;

    Double_t clusE = 0.0;
    for( Int_t blk = 0; blk<(int)nblk; blk++ ){
      Int_t blkid = int(cblkid[blk])-1; //-1 necessary since sbs.hcal.clus_blk.id ranges from 1 to 288 (different than just about everything else)
      //cout << blkid << endl;
      clusE += cblke[blk];
      A[blkid] += cblke[blk];
      NEV[blkid]++;
      if( nblk==1 ) h_CvCh->Fill(blkid,cblke[blk]/E_exp);
    }
    
    h_E_cut->Fill(clusE);
	
    // Double_t E_exp = KE_p*sampFrac;
    h_E_exp->Fill(E_exp);
    //Build the matrix as simply as possible
    for(Int_t mcol = 0; mcol<kNcell; mcol++){
      ba(mcol)+= A[mcol];
      for(Int_t mrow = 0; mrow<kNcell; mrow++){
	Ma(mcol,mrow) += A[mcol]*A[mrow]/E_exp;
	//Ma(mcol,mrow) += (A[mcol]*A[mrow])/(E_exp*E_exp);
      } 
    }
  }

  cout << endl << "Checking data, inverting matrix, and solving for coefficients.." << endl << endl;

  //Reject the bad cells
  Int_t badcell[kNcell] = {0};
  Double_t y[kNcell] = {0.0}; // For easy TGraphErrors build
  
  for(Int_t i=0; i<kNcell; i++){
    badcell[i] = 0;
    y[i] = i;
    
    //Do not change ADC gain coeff if insufficient events or energy dep in cell
    if( NEV[i] < minEventPerCell || Ma(i,i) < 0.1*ba(i) ){ 

      Double_t elemRatio = Ma(i,i)/ba(i);

      ba(i) = 1.0;  // Set RHS vector for cell i to 1.0 
      Ma(i,i) = 1.0; // Set diagonal element of Matrix M for cell i to 1.0 
      
      for(Int_t j=0; j<kNcell; j++){
	if( j != i ){
	  
	  Ma(i,j) = 0.0; 
	  Ma(j,i) = 0.0;
	}
      }
      badcell[i] = 1;
      //cout << "cell" << i << " bad" << endl;
      //cout << "Number of events in bad cell =" << NEV[i] << endl;
      //cout << "Matrix element/vector element ratio =" << elemRatio << endl;
    }  
  }

  
  //Invert the matrix, solve for ratios
  TMatrixD M_inv = Ma.Invert();
  TVectorD Coeff = M_inv*ba; // Stays unmodified for reference


  for(Int_t i=0; i<kNcell; i++){
    if(badcell[i]==0){
      GCoeff[i]=gOldConst[i]*Coeff[i]; // The new gain coefficient is the old coefficient multiplied by the gain factor that we just solved for. We will take these coefficents as the input for another set of data to estimate the error per channel.
      //cout << GCoeff[i] << endl;
      //GCoeff[i]=gOldConst[i]*Coeff[i]*1.714; // Making certain that these coefficients are doing what I need them to
      //cout << GCoeff[i] << endl;
      h_coeff->Fill(Coeff[i]);
    }else{
      GCoeff[i]=gOldConst[i]; // If the cell is bad, use the old coefficient
      //GCoeff[i]=gOldConst[i]*1.714; // If the cell is bad, use the old coefficient
    }

    //h_coeff->Fill(Coeff[i]);

    GCoeffDiff[i] = (GCoeff[i]-gOldConst[i])/gOldConst[i];
  }

  cout << "Inversion complete. Building histograms and writing out coefficients.." << endl << endl;

  Double_t yErr[kNcell] = {0.0};
  Double_t constErr[kNcell] = {0.0}; //Will need to improve error here

  TGraphErrors *ccgraph = new TGraphErrors( kNcell, y, GCoeff, yErr, constErr ); 
  ccgraph->GetXaxis()->SetLimits(0.0,kNcell);  
  ccgraph->GetYaxis()->SetLimits(0.0,0.25);
  ccgraph->SetTitle("Calibration Coefficients");
  ccgraph->GetXaxis()->SetTitle("Channel");
  ccgraph->GetXaxis()->SetTitle("Unitless");
  ccgraph->SetMarkerStyle(20); // idx 20 Circles, idx 21 Boxes
  ccgraph->Write("constants");

  TGraphErrors *ccdiffgraph = new TGraphErrors( kNcell, y, GCoeffDiff, yErr, constErr ); 
  ccdiffgraph->GetXaxis()->SetLimits(0.0,kNcell);  
  ccdiffgraph->GetYaxis()->SetLimits(0.0,0.25);
  ccdiffgraph->SetTitle("Calibration Coefficient Differences");
  ccdiffgraph->GetXaxis()->SetTitle("Channel");
  ccdiffgraph->GetXaxis()->SetTitle("Unitless");
  ccdiffgraph->SetMarkerStyle(20); // idx 20 Circles, idx 21 Boxes
  ccdiffgraph->Write("constantsDiff");


  Int_t cell = 0;

  cout << endl << "Number of events available for calibration: " << endl << endl;

  for( Int_t r = 0; r<kNrows; r++){
    for( Int_t c = 0; c<kNcols; c++){
      cout << NEV[cell] << "  ";
      cell++;
    }
    cout << endl;
  }
  cell = 0;


  cout << endl << "Raw Gain Coefficients: " << endl << endl;

  for( Int_t r = 0; r<kNrows; r++){
    for( Int_t c = 0; c<kNcols; c++){
      cout << Coeff[cell] << "  ";
      cell++;
    }
    cout << endl;
  }

  cell = 0;


  cout << endl << "Combined Gain Coefficients: " << endl << endl;

  for( Int_t r = 0; r<kNrows; r++){
    for( Int_t c = 0; c<kNcols; c++){
      cout << GCoeff[cell] << "  ";
      cell++;
    }
    cout << endl;
  }

  Double_t posErr[kNcell] = {0.};
  TF1 *f1;
  Double_t X[kNcell];
  Double_t Xval[kNcell];
  Double_t Xerr[kNcell];
  TH1D *Xslice[kNcell];
  for( Int_t x=0; x<kNcell; x++ ){
    X[x] = x;
    Xslice[x] = h_CvCh->ProjectionY(Form("Xslice_%d",x+1),x+1,x+1);
    Xslice[x]->Fit("gaus","Q");
    f1=Xslice[x]->GetFunction("gaus");
    if(Xslice[x]->GetEntries()>0){
      //Xval[x] = f1->GetParameter(1);
      //Xerr[x] = f1->GetParameter(2);
    }
  }

    
  cout << "Total calibration yield for run with current cuts: " << eyield << "." << endl; 

  fout->Write();
  fout->Close();
  
  cout << "Histograms populated and written to file." << endl;
  
}
