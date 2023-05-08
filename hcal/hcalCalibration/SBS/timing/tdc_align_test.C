//SSeeds 5.3.23 - Script adapted from tcal.C at https://github.com/sebastianseeds/HCal_replay
//NOTE: requires $DB_DIR path set correctly in current environment

#include <vector>
#include <iostream>
#include "TCut.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TTreeFormula.h"
#include "TChain.h"
#include "TLatex.h"
#include "TH1D.h"
#include "TStopwatch.h"
#include "../include/hcal.h"

// Set up data structure to hold configuration settings by field setting
typedef struct{
  vector<string> tstamp;
  vector<string> run;
  string gcut;
  Double_t SBSf;
  Double_t E_e;
  Double_t HCal_d;
  Double_t HCal_th;
  Double_t BB_th;
  Double_t W2_mean;
  Double_t W2_sig;
  Double_t dx0_n;
  Double_t dx0_p;
  Double_t dy0;
  Double_t dx_sig_n;
  Double_t dx_sig_p;
  Double_t dy_sig;
  Double_t atime0;
  Double_t atime_sig;
  Int_t useAlshield;
  vector<Double_t> otdcoff;
} PARA;

// Use nested structure to index parameters by target
typedef struct{
  string target;
  PARA para[hcal::gMaxfset];
} TARPARA;

// Create structure to hold configuration files
typedef struct{
  string targ;
  vector<string> config;
} CONF;

//Main <experiment> <configuration> <quasi-replay>; qreplay should only be performed after new offsets obtained
void tdc_align_test( const char *exp = "gmn", Int_t config=4, bool qreplay = false ){

  // Setup qreplay int index for reporting
  Int_t qreplayIdx = 0;
  if( qreplay )
    qreplayIdx = 1;

  // Get index for calibration type
  Int_t calidx;
  for( Int_t s=0; s<hcal::gNstamp; s++ )
    if( hcal::gStamp[s].compare("tdcoff")==0 )
      calidx=s;
  
  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );
  
  // Get the date
  string date = hcal::getDate();

  // Set up output analysis file
  TFile *fout = new TFile( Form("outfiles/tdcalign_%s_conf%d_qr%d.root",exp,config,qreplayIdx), "RECREATE" );

  // Declare path and open log file for reporting
  ofstream report;
  string logPath = Form("logs/log_%s_conf%d_qreplay%d.log",exp,config,qreplayIdx);

  report.open( logPath );
  report << "#Error and performance report from " << exp << " config " << config << " obtained " << date.c_str() << endl << endl;

  cout << "Opened log" << endl;

  // Read in meta-config for multiple field settings within a given SBS config
  string metaconfigname = Form("../config/%s/conf%d/%s_conf%d_meta.cfg",exp,config,exp,config);
  ifstream metaconfigfile(metaconfigname);

  // Set up structure to hold all configuration files by target and record number of field settings per
  CONF targ_config[hcal::gNtar];

  cout << "Made meta config structure" << endl;

  TString mcline;
  for( Int_t f=0; f<hcal::gNtar; f++ ){

    cout << "target from hcal namespace: " << hcal:: gTar[f] << endl;

    // targ_config[f].targ = hcal::gTar[f];
    // string end = "end" + hcal::gTar[f];

    // cout << "end string: " << endl;

    // while( mcline.ReadLine( metaconfigfile ) && !mcline.BeginsWith(end) ){
    //   if( !mcline.BeginsWith("#") ){
    // 	string configpath = Form("../config/%s/conf%d/",exp,config);

    // 	configpath += (string)mcline;

    // 	cout << "configpath: " << configpath << endl;

    // 	targ_config[f].config.push_back(configpath);
    // 	cout << "Recording " << hcal::gTar[f] << " config file: " << configpath << ".." << endl;
    //   }    
    // }
  } 


}
//   cout << "arbitrary config file loaded check: " << targ_config[0].config[0] << endl;

//   // Create array of structures to hold all parameters by target
//   TARPARA allp[hcal::gNtar];

//   // Path to read in old tdc offsets from database
//   string DBpath = gSystem->Getenv("DB_DIR");
//   TString tdcOffsetPath = DBpath + "/db_sbs.hcal.dat";
//   ifstream tdcOffsetFile( tdcOffsetPath );
//   cout << endl << "Opening previous offsets from database file: " << tdcOffsetPath << ".." << endl;
//   if( !tdcOffsetFile ){
//     cerr << endl << "ERROR: No input constant file present -> path to db_sbs.hcal.dat expected." << endl;
//     return 0;
//   }

//   // Loop over potential targets
//   for( Int_t t=0; t<hcal::gNtar; t++ ){
//     string target = hcal::gTar[t];
//     allp[t].target = target;
//     Int_t tarfset = targ_config[t].config.size();

//     // if( tarfset<1 )
//     //   continue;

//     // Loop over all available field settings
//     for( Int_t f=0; f<tarfset; f++ ){

//       cout << "Reading in " << exp << " config " << config << " " << target << " analysis parameters.."  << endl;

//       // Record timestamps
//       ifstream configfile(targ_config[t].config[f]);
//       TString currentline;
//       Int_t stampidx=0;
//       while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endstamp") ){
// 	if( !currentline.BeginsWith("#") ){
// 	  allp[t].para[f].tstamp.push_back((string)currentline);
// 	  cout << "Recording " << hcal::gStamp[stampidx] << " timestamp: " << currentline << ".." << endl;
// 	  stampidx++;
// 	}    
//       }

//       // Record data files
//       while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endlist") ){
// 	if( !currentline.BeginsWith("#") ){
// 	  allp[t].para[f].run.push_back((string)currentline);
// 	  cout << "Recording file: " << currentline << ".." << endl;
// 	}    
//       }

//       // Record globalcuts
//       while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endcut") ){
// 	if( !currentline.BeginsWith("#") ){
// 	  allp[t].para[f].gcut = (string)currentline;
// 	  cout << "Recording globalcut: " << currentline << ".." << endl;
// 	}    
//       }
      
//       // Record physics cuts and parameters
//       while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("#") ){
// 	TObjArray *tokens = currentline.Tokenize(" ");
// 	Int_t ntokens = tokens->GetEntries();
// 	if( ntokens>1 ){
// 	  TString skey = ( (TObjString*)(*tokens)[0] )->GetString();
// 	  if( skey == "SBSf" ){
// 	    TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
// 	    allp[t].para[f].SBSf = sval.Atof();
// 	    cout << "Recording beam energy: " << allp[t].para[f].SBSf << endl;
// 	  }
// 	  if( skey == "E_e" ){
// 	    TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
// 	    allp[t].para[f].E_e = sval.Atof();
// 	    cout << "Loading beam energy: " << allp[t].para[f].E_e << endl;
// 	  }
// 	  if( skey == "HCal_d" ){
// 	    TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
// 	    allp[t].para[f].HCal_d = sval.Atof();
// 	    cout << "Loading HCal distance: " << allp[t].para[f].HCal_d << endl;
// 	  }
// 	  if( skey == "HCal_th" ){
// 	    TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
// 	    allp[t].para[f].HCal_th = sval.Atof() * TMath::DegToRad();	
// 	    cout << "Loading HCal angle: " << allp[t].para[f].HCal_th << endl;
// 	  }
// 	  if( skey == "BB_th" ){
// 	    TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
// 	    allp[t].para[f].BB_th = sval.Atof() * TMath::DegToRad();	
// 	    cout << "Loading BBCal angle: " << allp[t].para[f].BB_th << endl;
// 	  }
// 	  if( skey == "W2_mean" ){
// 	    TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
// 	    allp[t].para[f].W2_mean = sval.Atof();
// 	    cout << "Loading W2 mean cut: " << allp[t].para[f].W2_mean << endl;
// 	  }
// 	  if( skey == "W2_sig" ){
// 	    TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
// 	    allp[t].para[f].W2_sig = sval.Atof();
// 	    cout << "Loading W2 sigma cut: " << allp[t].para[f].W2_sig << endl;
// 	  }
// 	  if( skey == "dx0_n" ){
// 	    TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
// 	    allp[t].para[f].dx0_n = sval.Atof();
// 	    cout << "Loading x position of neutron spot: " << allp[t].para[f].dx0_n << endl;
// 	  }
// 	  if( skey == "dx0_p" ){
// 	    TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
// 	    allp[t].para[f].dx0_p = sval.Atof();
// 	    cout << "Loading y position of proton spot: " << allp[t].para[f].dx0_p << endl;
// 	  }
// 	  if( skey == "dy0" ){
// 	    TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
// 	    allp[t].para[f].dy0 = sval.Atof();
// 	    cout << "Loading y position of both hadron spots: " << allp[t].para[f].dy0 << endl;
// 	  }
// 	  if( skey == "dx_sig_n" ){
// 	    TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
// 	    allp[t].para[f].dx_sig_n = sval.Atof();
// 	    cout << "Loading x sigma of neutron spot: " << allp[t].para[f].dx_sig_n << endl;
// 	  }
// 	  if( skey == "dx_sig_p" ){
// 	    TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
// 	    allp[t].para[f].dx_sig_p = sval.Atof();
// 	    cout << "Loading x sigma of proton spot: " << allp[t].para[f].dx_sig_p << endl;
// 	  }
// 	  if( skey == "dy_sig" ){
// 	    TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
// 	    allp[t].para[f].dy_sig = sval.Atof();
// 	    cout << "Loading y sigma of both hadron spots: " << allp[t].para[f].dy_sig << endl;
// 	  }
// 	  if( skey == "atime0" ){
// 	    TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
// 	    allp[t].para[f].atime0 = sval.Atof();
// 	    cout << "Loading ADC time mean: " << allp[t].para[f].atime0 << endl;
// 	  }
// 	  if( skey == "atime_sig" ){
// 	    TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
// 	    allp[t].para[f].atime_sig = sval.Atof();
// 	    cout << "Loading ADC time sigma: " << allp[t].para[f].atime_sig << endl;
// 	  }
// 	  if( skey == "useAlshield" ){
// 	    TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
// 	    allp[t].para[f].useAlshield = sval.Atoi();
// 	    cout << "Loading Aluminum absorber option: " << allp[t].para[f].useAlshield << endl;
// 	  }
// 	}
// 	delete tokens;
//       } 

//       //Record old tdc offsets by tstamp from database
//       Int_t n0=0;
//       Double_t d0=0;
//       string line;
//       bool read_offset = false;
//       bool found_tstamp = false;
//       if( allp[t].para[f].tstamp[calidx].compare("none")==0 )
// 	found_tstamp = true;
      
//       while( getline( tdcOffsetFile, line ) ){
	
// 	if( n0==hcal::maxHCalChan ) break;
	
// 	TString Tline = (TString)line;
	
// 	if( Tline.BeginsWith(allp[t].para[f].tstamp[calidx]) && !found_tstamp ){
// 	  found_tstamp = true;
// 	  continue;
// 	}
	
// 	if( Tline.BeginsWith("sbs.hcal.tdc.offset") && found_tstamp && !read_offset ){
// 	  read_offset = true;
// 	  continue;
// 	}

// 	if( found_tstamp && read_offset ){
// 	  istringstream iss( line );
// 	  while( iss >> d0 ){
// 	    allp[t].para[f].otdcoff.push_back( d0 );
// 	    n0++;
// 	  }
// 	}
//       }
	
//       cout << endl << endl << "Old TDC offsets: " << endl;
	
//       for( Int_t r=0; r<hcal::maxHCalRows; r++){
// 	for( Int_t c=0; c<hcal::maxHCalCols; c++){
// 	  Int_t i = r*hcal::maxHCalCols+c;
// 	  cout << allp[t].para[f].otdcoff[i] << " ";
// 	}
// 	cout << endl;
//       }
	
//     }//endloop over field settings
//   }//endloop over targets

//   // Declare paths for new tdc offsets
//   string writeParPath = Form("parameters/tdcoffsets_%s_conf%d.txt",exp,config);
//   string readParPath = Form("parameters/tdcoffsets_%s_conf%d.txt",exp,config);
//   Double_t ntdcoff[hcal::maxHCalChan];

//   // Recording fit parameters if quasi-replay selected
//   if( qreplay ){
//     ifstream readParFile( readParPath );
//     if( !readParFile ){
//       cerr << endl << "ERROR: No input constant file present -> parameters/<tdcoffsetfile> expected." << endl;
//       return 0;
//     }

//     cout << "Loading tdc offset parameters.." << endl;

//     // Record new tdc offsets
//     Int_t n0=0;
//     Double_t d0=0;
//     string line;
//     bool read_offset = false;
    
//     while( getline( readParFile, line ) ){
      
//       if( n0==hcal::maxHCalChan ) break;
      
//       TString Tline = (TString)line;
 
//       if( Tline.BeginsWith("sbs.hcal.tdc.offset") && !read_offset ){
// 	  read_offset = true;
// 	  continue;
//       }
      
//       if( read_offset ){
// 	istringstream iss( line );
// 	while( iss >> d0 ){
// 	  ntdcoff[n0]=d0;
// 	  n0++;
// 	}
//       }
//     }
    
//     cout << endl << endl << "New TDC offsets from file: " << endl;
    
//     for( Int_t r=0; r<hcal::maxHCalRows; r++){
//       for( Int_t c=0; c<hcal::maxHCalCols; c++){
// 	Int_t i = r*hcal::maxHCalCols+c;
// 	cout << ntdcoff[i] << " ";
//       }
//       cout << endl;
//     }

//   }//endif qreplay

//   // Re-allocate memory at each run to load different cuts/parameters
//   TChain *C = nullptr;

//   // Create simple output tree
//   TTree *P = new TTree("P","Analysis Tree");

//   // Timing
//   Double_t pblkid_out; //hcal primary cluster, primary block id
//   Double_t tdc_out; //hcal primary cluster tdc time, tree
//   Double_t atime_out; //hcal primary cluster adc time, tree
//   Double_t hodotmean_out; //hodoscope primary cluster mean tdc time

//   // Physics
//   Double_t dx_out; //hcal primary cluster dx
//   Double_t dy_out; //hcal primary cluster dy
//   Double_t W2_out; //W2
//   Double_t Q2_out; //Q2
//   Double_t hcale_out; //hcal primary cluster energy
//   Int_t mag_out; //sbs magnetic field strength (percent)
//   Int_t run_out; //run number
//   Int_t tar_out; //target, LH2 or LD2

//   // Set output tree branches
//   P->Branch( "pblkid", &pblkid_out, "pblkid/D" );
//   P->Branch( "tdc", &tdc_out, "tdc/D" );
//   P->Branch( "atime", &atime_out, "atime/D" );
//   P->Branch( "hodotmean", &hodotmean_out, "hodotmean/D" );
//   P->Branch( "dx", &dx_out, "dx/D" );
//   P->Branch( "dy", &dy_out, "dy/D" );
//   P->Branch( "W2", &W2_out, "W2/D" );
//   P->Branch( "Q2", &Q2_out, "Q2/D" );
//   P->Branch( "hcale", &hcale_out, "hcale/D" );
//   P->Branch( "mag", &mag_out, "mag_out/I" );
//   P->Branch( "run", &run_out, "run_out/I" );
//   P->Branch( "tar", &tar_out, "tar_out/I" );

//   // Begin looping over data

//   // Loop over potential targets
//   for( Int_t t=0; t<hcal::gNtar; t++ ){
    
//     string target = hcal::gTar[t];
//     Int_t tarfset = targ_config[t].config.size();

//     // if( tarfset<1 )
//     //   continue;

//     // Loop over all available field settings
//     for( Int_t f=0; f<tarfset; f++ ){
      
//       Double_t mag = allp[t].para[f].SBSf;

//       cout << "Looping over data in " << exp << " config " << config << " " << target << ".."  << endl;
      
//       Int_t nruns = allp[t].para[f].run.size();

//       for( Int_t r=0; r<nruns; r++ ){
	
// 	cout << "Analyzing run " << r << ".." << endl;
	
// 	C = new TChain("T");
// 	C->Add(allp[t].para[f].run[r].c_str());
	
// 	// Setting up ROOT tree branch addresses
// 	C->SetBranchStatus("*",0);    
// 	Double_t BBtr_p[hcal::maxTracks], BBtr_px[hcal::maxTracks], BBtr_py[hcal::maxTracks], BBtr_pz[hcal::maxTracks];
// 	Double_t BBtr_vz[hcal::maxTracks];
// 	Double_t BBtr_n, BBps_x, BBps_y, BBps_e, BBsh_x, BBsh_y, BBsh_e;	
// 	Double_t HCALx, HCALy, HCALe;
// 	Double_t crow, ccol, nblk;
// 	Double_t cblkid[hcal::maxHCalChan], cblke[hcal::maxHCalChan], cblkatime[hcal::maxHCalChan], cblktime[hcal::maxHCalChan];
// 	Double_t cblkagain[hcal::maxHCalChan];
// 	Double_t HODOtmean;
// 	// Speed up processing by switching on only useful branches
// 	C->SetBranchStatus( "*", 0 );
// 	C->SetBranchStatus( "sbs.hcal.x", 1 );
// 	C->SetBranchStatus( "sbs.hcal.y", 1 );
// 	C->SetBranchStatus( "sbs.hcal.e", 1 );
// 	C->SetBranchStatus( "sbs.hcal.rowblk", 1 );
// 	C->SetBranchStatus( "sbs.hcal.colblk", 1 );
// 	C->SetBranchStatus( "sbs.hcal.nblk", 1 );
// 	C->SetBranchStatus( "sbs.hcal.clus_blk.id", 1 );
// 	C->SetBranchStatus( "sbs.hcal.clus_blk.e", 1 );
// 	C->SetBranchStatus( "sbs.hcal.clus_blk.tdctime", 1 );
// 	C->SetBranchStatus( "sbs.hcal.clus_blk.atime", 1 );
// 	C->SetBranchStatus( "sbs.hcal.clus_blk.again", 1 );
// 	C->SetBranchStatus( "bb.tr.n", 1 );
// 	C->SetBranchStatus( "bb.tr.px", 1 );
// 	C->SetBranchStatus( "bb.tr.py", 1 );
// 	C->SetBranchStatus( "bb.tr.pz", 1 );    
// 	C->SetBranchStatus( "bb.tr.vz", 1 );
// 	C->SetBranchStatus( "bb.tr.p", 1 );
// 	C->SetBranchStatus( "bb.ps.e", 1 );
// 	C->SetBranchStatus( "bb.ps.x", 1 );
// 	C->SetBranchStatus( "bb.ps.y", 1 );
// 	C->SetBranchStatus( "bb.sh.e", 1 );
// 	C->SetBranchStatus( "bb.sh.x", 1 );
// 	C->SetBranchStatus( "bb.sh.y", 1 );
// 	C->SetBranchStatus( "bb.hodotdc.clus.tmean", 1 );
// 	// Linking memory
// 	C->SetBranchAddress( "sbs.hcal.x", &HCALx );
// 	C->SetBranchAddress( "sbs.hcal.y", &HCALy );
// 	C->SetBranchAddress( "sbs.hcal.e", &HCALe );
// 	C->SetBranchAddress( "sbs.hcal.rowblk", &crow );
// 	C->SetBranchAddress( "sbs.hcal.colblk", &ccol );
// 	C->SetBranchAddress( "sbs.hcal.nblk", &nblk ); // Total2 number of blocks in highest E clus
// 	C->SetBranchAddress( "sbs.hcal.clus_blk.id", cblkid ); // kNcell-1 index for each block
// 	C->SetBranchAddress( "sbs.hcal.clus_blk.e", cblke ); // Array of block energies
// 	C->SetBranchAddress( "sbs.hcal.clus_blk.tdctime", cblktime ); // Array of block TDC times
// 	C->SetBranchAddress( "sbs.hcal.clus_blk.atime", cblkatime ); // Array of block ADC times
// 	C->SetBranchAddress( "sbs.hcal.clus_blk.again", cblkagain ); // Array of block ADC gain coeff
// 	C->SetBranchAddress( "bb.tr.n", &BBtr_n );
// 	C->SetBranchAddress( "bb.tr.px", BBtr_px );
// 	C->SetBranchAddress( "bb.tr.py", BBtr_py );
// 	C->SetBranchAddress( "bb.tr.pz", BBtr_pz );
// 	C->SetBranchAddress( "bb.tr.vz", BBtr_vz );
// 	C->SetBranchAddress( "bb.tr.p", BBtr_p );
// 	C->SetBranchAddress( "bb.ps.e", &BBps_e );
// 	C->SetBranchAddress( "bb.ps.x", &BBps_x );
// 	C->SetBranchAddress( "bb.ps.y", &BBps_y );
// 	C->SetBranchAddress( "bb.sh.e", &BBsh_e );
// 	C->SetBranchAddress( "bb.sh.x", &BBsh_x );
// 	C->SetBranchAddress( "bb.sh.y", &BBsh_y ); 
// 	C->SetBranchAddress( "bb.hodotdc.clus.tmean", &HODOtmean );

// 	TCut GCut = allp[t].para[f].gcut.c_str();
// 	TTreeFormula *GlobalCut = new TTreeFormula( "GlobalCut", GCut, C );

// 	// Set up hcal coordinate system with hcal angle wrt exit beamline
// 	vector<TVector3> hcalaxes; hcal::sethcalaxes( allp[t].para[f].HCal_th, hcalaxes );
// 	TVector3 hcalorigin = allp[t].para[f].HCal_d*hcalaxes[2] + hcal::HCalvoff*hcalaxes[0];
// 	Double_t BdL = hcal::maxSBSfield * hcal::sbsdipolegap * (mag/100); //scale crudely by magnetic field
// 	Double_t Eloss_outgoing = hcal::celldiameter/2.0/sin(allp[t].para[f].BB_th) * hcal::tarrho[t] * hcal::dEdx[t];

// 	long nevent = 0, nevents = C->GetEntries(); 
// 	Int_t treenum = 0, currenttreenum = 0;

// 	while (C->GetEntry(nevent++)) {

// 	  cout << nevent << "/" << nevents << " \r";
// 	  cout.flush();

// 	  ///////
// 	  //Single-loop globalcut method. Save pass/fail for output tree.
// 	  currenttreenum = C->GetTreeNumber();
// 	  if( nevent == 1 || currenttreenum != treenum ){
// 	    treenum = currenttreenum; 
// 	    GlobalCut->UpdateFormulaLeaves();
// 	    cout << "Updating formula leaves and switching segment at event: " << nevent << endl;
// 	  }
// 	  bool failedglobal = GlobalCut->EvalInstance(0) == 0;

// 	  if( failedglobal ) continue;
	  
// 	  ///////
// 	  //Physics calculations
// 	  //correct beam energy with vertex information, primary track
// 	  Double_t ebeam_c = allp[t].para[f].E_e - ( (BBtr_vz[0]+hcal::l_tgt[t]/2.0) * hcal::tarrho[t] * hcal::dEdx[t] + hcal::uwallthick[t] * hcal::crho[t] * hcal::cdEdx[t] );

// 	  TVector3 vertex( 0., 0., BBtr_vz[0] );

// 	  //reconstructed momentum, corrected for mean energy loss exiting the target (later we'll also want to include Al shielding on scattering chamber)
// 	  Double_t precon = BBtr_p[0] + Eloss_outgoing;

// 	  //set up four-momenta with some empty for various calculation methods
// 	  TLorentzVector pbeam( 0., 0., ebeam_c, ebeam_c ); //beam momentum
// 	  //TLorentzVector pe( px[0], py[0], pz[0], p[0] ); //e' plvect
// 	  TLorentzVector pe( precon*BBtr_px[0]/BBtr_p[0], precon*BBtr_py[0]/BBtr_p[0], precon*BBtr_pz[0]/BBtr_p[0], precon ); //e' recon plvect
// 	  TLorentzVector ptarg; //target momentum
// 	  ptarg.SetPxPyPzE( 0., 0., 0., hcal::M_t[t] );
// 	  TLorentzVector q = pbeam - pe; //virtual photon mom. or mom. transferred to scatter nucleon (N')
// 	  TVector3 qv = q.Vect();
// 	  TLorentzVector pN; //N' momentum
      
// 	  //simple calculations for e' and N'
// 	  Double_t etheta = acos( pe.Pz() / pe.E() );
// 	  Double_t ephi = atan2( pe.Py(), pe.Px() );
// 	  Double_t pcent = ebeam_c/( 1. + ( ebeam_c/hcal::M_t[t] )*( 1.0 - cos(etheta) ) ); //e' p reconstructed by angles
// 	  Double_t phNexp = ephi + hcal::PI;
// 	  Double_t Q2, W2;

// 	  //e' p reconstruction with track angles (not momentum)
// 	  Double_t nu = pbeam.E() - pcent;
// 	  Double_t pNexp = sqrt( pow(nu, 2.) + 2. * hcal::M_t[t] * nu );
// 	  Double_t thNexp = acos( (pbeam.E()-pcent*cos(etheta)) / pNexp );
// 	  TVector3 pNhat( sin(thNexp) * cos(phNexp), sin(thNexp) * sin(phNexp), cos(thNexp) );
// 	  pN.SetPxPyPzE( pNexp*pNhat.X(), pNexp*pNhat.Y(), pNexp*pNhat.Z(), nu+ptarg.E() );
// 	  Q2 = 2.0 * pbeam.E() * pe.E() * ( 1.0 - cos(etheta) );
// 	  W2 = pow( hcal::M_t[t], 2.0 ) + 2.0*hcal::M_t[t] * (ebeam_c-pe.E()) - Q2;

// 	  //Calculate h-arm quantities
// 	  vector<Double_t> xyhcalexp; hcal::getxyhcalexpect( vertex, pNhat, hcalorigin, hcalaxes, xyhcalexp );
// 	  TVector3 hcalpos = hcalorigin + HCALx*hcalaxes[0] + HCALy*hcalaxes[1]; //from primary blk
// 	  Double_t dx = HCALx - xyhcalexp[0];
// 	  Double_t dy = HCALy - xyhcalexp[1];

// 	  ///////
// 	  //BigBite/SBS Acceptance matching.
// 	  bool failedaccmatch = 
// 	    xyhcalexp[1] > hcal::posHCalYf ||
// 	    xyhcalexp[1] < hcal::posHCalYi ||
// 	    xyhcalexp[0] > hcal::posHCalXf ||
// 	    xyhcalexp[0] < hcal::posHCalXi;
	  
// 	  if( failedaccmatch ) continue;

// 	  pblkid_out = (double)cblkid[0];
// 	  tdc_out = cblktime[0];
// 	  atime_out = cblkatime[0];
// 	  hcale_out = HCALe;
// 	  dx_out = dx;
// 	  dy_out = dy;
// 	  W2_out = W2;
// 	  Q2_out = Q2;
// 	  mag_out = mag;
// 	  run_out = r;
// 	  tar_out = t;
// 	  hodotmean_out = HODOtmean;

// 	  P->Fill();

	  
// 	}//end loop over event

// 	// getting ready for the next run
// 	C->Reset();

//       }//endloop over runs
      
//     }//endloop over fields
  
//   }//endloop over targets
  
//   //Make slices of TDC v ID and fit for offsets


//   fout->Write();

//   st->Stop();

//   // Send time efficiency report to console
//   cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

// }

//   //Declare limits on histograms to avoid fitting mistakes
//   Double_t tllim = -150.;
//   Double_t tulim = 0.;
//   Double_t ttotb = 2*abs(tllim-tulim);
//   Double_t allim = 0.;
//   Double_t aulim = 150.;
//   Double_t atotb = 2*abs(allim-aulim);

//   //Declare physics histograms for match to MC
//   TH1D *hpp_p = new TH1D( "hpp_p", "Expected proton momentum from BB; GeV", 100, 0, 10 );
//   TH1D *hpp_n = new TH1D( "hpp_n", "Expected neutron momentum from BB; GeV", 100, 0, 10 );

//   //Declare diagnostic histograms (as sparse as possible)
//   //Delta plots
//   TH2D *hdx_mag_h = new TH2D( "hdx_mag_h", "Delta X vs Field Setting (LH2); field (percent); x_{HCAL} - x_{exp} (m)", 21, 0, 105, 250, -4, 3 );
//   TH2D *hdy_mag_h = new TH2D( "hdy_mag_h", "Delta Y vs Field Setting (LH2); field (percent); y_{HCAL} - y_{exp} (m)", 21, 0, 105, 250, -1.25, 1.25 );
//   TH2D *hdx_mag_d = new TH2D( "hdx_mag_d", "Delta X vs Field Setting (LD2); field (percent); x_{HCAL} - x_{exp} (m)", 21, 0, 105, 250, -4, 3 );
//   TH2D *hdy_mag_d = new TH2D( "hdy_mag_d", "Delta Y vs Field Setting (LD2); field (percent); y_{HCAL} - y_{exp} (m)", 21, 0, 105, 250, -1.25, 1.25 );

//   //TDC plots
//   //1D over all tdc corrections
//   TH1D *htdc = new TH1D("htdc","TDC (All Channels); ns",ttotb,tllim,tulim);
//   TH1D *htdc_hodocorr = new TH1D("htdc_hodocorr","TDC Hodo Corrected (All Channels); ns",ttotb,tllim,tulim);
//   TH1D *htdc_twcorr = new TH1D("htdc_twcorr","TDC TWalk Corrected (All Channels); ns",ttotb,tllim,tulim);
//   TH1D *htdc_tofpcorr_p = new TH1D("htdc_tofpcorr_p","TDC Proton TOF (vs ep) Corrected (All Channels); ns",ttotb+50,tllim-50,tulim);
//   TH1D *htdc_tofpcorr_n = new TH1D("htdc_tofpcorr_n","TDC Neutron TOF (vs ep) Corrected (All Channels); ns",ttotb+50,tllim-50,tulim);
//   TH1D *htdc_tofidcorr_p = new TH1D("htdc_tofidcorr_p","TDC Proton TOF (vs pos) Corrected (All Channels); ns",ttotb+50,tllim-50,tulim);
//   TH1D *htdc_tofidcorr_n = new TH1D("htdc_tofidcorr_n","TDC Neutron TOF (vs pos) Corrected (All Channels); ns",ttotb+50,tllim-50,tulim);
//   TH1D *htdc_allcorr_p = new TH1D("htdc_allcorr_p","TDC Proton Corrected TW/TOF/Hodo (All Channels); ns",ttotb+50,tllim-50,tulim);
//   TH1D *htdc_allcorr_n = new TH1D("htdc_allcorr_n","TDC Neutron Corrected TW/TOF/Hodo (All Channels); ns",ttotb+50,tllim-50,tulim);
//   TH1D *htp_allcorr_p = new TH1D("htp_allcorr_p","Primary Block TDC Proton Corrected TW/TOF/Hodo (All Channels); ns",ttotb+50,tllim-50,tulim);
//   TH1D *htp_allcorr_n = new TH1D("htp_allcorr_n","Primary Block TDC Neutron Corrected TW/TOF/Hodo (All Channels); ns",ttotb+50,tllim-50,tulim);
//   //2D over all primary blocks and corrections
//   TH2D *htp_ID = new TH2D("htp_ID","TDC Primary Block;Channel;TDC_{HCAL} (ns)",288,0,288,ttotb,tllim,tulim);
//   TH2D *htp_hodocorr_ID = new TH2D("htp_hodocorr_ID","TDC Primary Block - TDC hodo;Channel;TDC_{HCAL}-TDC_{HODO} (ns)",288,0,288,ttotb,tllim,tulim);
//   TH2D *htp_twcorr_ID = new TH2D("htp_twcorr_ID","TDC Primary Block Timewalk Corrected;Channel;TDC_{HCAL}-Corr_{TW} (ns)",288,0,288,ttotb,tllim,tulim);
//   TH2D *htp_tofpcorr_p_ID = new TH2D("htp_tofpcorr_p_ID","TDC Proton Primary Block TOF vs Electron Momentum Corrected;Channel;TDC_{HCAL}-Corr_{TOF} (ns)",288,0,288,ttotb+50,tllim-50,tulim);
//   TH2D *htp_tofpcorr_n_ID = new TH2D("htp_tofpcorr_n_ID","TDC Neutron Primary Block TOF vs Electron Momentum Corrected;Channel;TDC_{HCAL}-Corr_{TOF} (ns)",288,0,288,ttotb+50,tllim-50,tulim);
//   TH2D *htp_tofidcorr_p_ID = new TH2D("htp_tofidcorr_p_ID","TDC Proton Primary Block TOF vs ID Corrected;Channel;TDC_{HCAL}-Corr_{TOF} (ns)",288,0,288,ttotb+50,tllim-50,tulim);
//   TH2D *htp_tofidcorr_n_ID = new TH2D("htp_tofidcorr_n_ID","TDC Neutron Primary Block TOF vs ID Corrected;Channel;TDC_{HCAL}-Corr_{TOF} (ns)",288,0,288,ttotb+50,tllim-50,tulim);
//   TH2D *htp_allcorr_p_ID = new TH2D("htp_allcorr_p_ID","TDC Proton Primary Block TOF/Hodo Corrected;Channel;TDC_{HCAL}-Corr_{TOF} (ns)",288,0,288,ttotb+50,tllim-50,tulim);
//   TH2D *htp_allcorr_n_ID = new TH2D("htp_allcorr_n_ID","TDC Neutron Primary Block TOF/Hodo Corrected;Channel;TDC_{HCAL}-Corr_{TOF} (ns)",288,0,288,ttotb+50,tllim-50,tulim);
//   //2D over all cluster blocks and corrections
//   TH2D *ht_ID = new TH2D("ht_ID","TDC;Channel;TDC_{HCAL} (ns)",288,0,288,ttotb,tllim,tulim);
//   TH2D *ht_hodocorr_ID = new TH2D("ht_hodocorr_ID","TDC Hodo Corrected (all cluster elements);Channel;TDC_{HCAL}-TDC_{HODO} (ns)",288,0,288,ttotb,tllim,tulim);
//   TH2D *ht_twcorr_ID = new TH2D("ht_twcorr_ID","TDC TW Corrected (all cluster elements);Channel;TDC_{HCAL}, TW_{Corr} (ns)",288,0,288,ttotb,tllim,tulim);
//   TH2D *ht_tofpcorr_p_ID = new TH2D("ht_tofpcorr_p_ID","TDC Proton TOF (vs ep) Corrected (all cluster elements);Channel;TDC_{HCAL} TOF v p Corrected (ns)",288,0,288,ttotb+50,tllim-50,tulim);
//   TH2D *ht_tofpcorr_n_ID = new TH2D("ht_tofpcorr_n_ID","TDC Neutron TOF (vs ep) Corrected (all cluster elements);Channel;TDC_{HCAL} TOF v p Corrected (ns)",288,0,288,ttotb+50,tllim-50,tulim);
//   TH2D *ht_tofidcorr_p_ID = new TH2D("ht_tofidcorr_p_ID","TDC Proton TOF (vs ID) Corrected (all cluster elements);Channel;TDC_{HCAL} TOF v ID Corrected (ns)",288,0,288,ttotb+50,tllim-50,tulim);
//   TH2D *ht_tofidcorr_n_ID = new TH2D("ht_tofidcorr_n_ID","TDC Neutron TOF (vs ID) Corrected (all cluster elements);Channel;TDC_{HCAL} TOF v ID Corrected (ns)",288,0,288,ttotb+50,tllim-50,tulim);
//   TH2D *ht_allcorr_p_ID = new TH2D("ht_allcorr_p_ID","TDC Proton TW/TOF/Hodo Corrected (all cluster elements);Channel;TDC_{HCAL}-TDC_{HODO} w/TWalk w/TOF (ns)",288,0,288,ttotb+50,tllim-50,tulim);
//   TH2D *ht_allcorr_n_ID = new TH2D("ht_allcorr_n_ID","TDC Neutron TW/TOF/Hodo Corrected (all cluster elements);Channel;TDC_{HCAL}-TDC_{HODO} w/TWalk w/TOF (ns)",288,0,288,ttotb+50,tllim-50,tulim);

//   //ADCt plots
//   TH1D *hadct = new TH1D("hadct","ADCt (All Channels); ns",atotb,allim,aulim);
//   TH1D *hadct_corr = new TH1D("hadct_corr","ADCt TWalk/Hodo Corrected (All Channels); ns",atotb,allim,aulim);
//   TH2D *ha_ID = new TH2D("ha_ID","ADCt;Channel;ADCtime_{HCAL} (ns)",288,0,288,atotb,allim,aulim);
//   TH2D *hap_ID = new TH2D("hap_ID","ADCt Primary Block;Channel;ADCtime_{HCAL} (ns)",288,0,288,atotb,allim,aulim);
//   TH2D *hapDiff_ID = new TH2D("hapDiff_ID","ADCt Primary Block - TDC hodo;Channel;ADCtime_{HCAL} - TDC_{HODO} (ns)",288,0,288,atotb,allim,aulim);
//   TH2D *hapCorr_ID = new TH2D("hapCorr_ID","ADCt Primary Block, Corrected;Channel;ADCtime_{HCAL}-TDC_{HODO}-Corr_{TW} (ns)",288,0,288,atotb,allim,aulim);
//   TH2D *haDiff_ID = new TH2D("haDiff_ID",";Channel;ADCtime_{HCAL}-TDC_{HODO} (ns)",288,0,288,atotb,allim,aulim);
//   TH2D *haCorr_ID = new TH2D("haCorr_ID",";Channel;ADCtime_{HCAL}-TDC_{HODO}-Corr_{TW} (ns)",288,0,288,atotb,allim,aulim);

//   //Declare self timing histograms
//   TH1D *hclusmean = new TH1D( "hclusmean", "HCal Cluster Mean Time (Primary Block Reftime, nblk>1)", 1000, -50, 50 );
//   TH2D *hclusmeanID = new TH2D( "hclusmeanID", "HCal Cluster Mean Time, (Primary Block Reftime/ID, nblk>1), vs ID", 288, 0, 288, 1000, -50, 50 );
//   TH1D *hclusmean_corrp = new TH1D( "hclusmean_corrp", "HCal Proton Cluster Mean Time TW Corrected (Primary Block Reftime, nblk>1)", 1000, -50, 50 );
//   TH2D *hclusmeanID_corrp = new TH2D( "hclusmeanID_corrp", "HCal Proton Cluster Mean Time TW Corrected (Primary Block Reftime/ID, nblk>1), vs ID", 288, 0, 288, 1000, -50, 50 );
//   TH1D *hclusmean_corrn = new TH1D( "hclusmean_corrn", "HCal Neutron Cluster Mean Time TW Corrected (Primary Block Reftime, nblk>1)", 1000, -50, 50 );
//   TH2D *hclusmeanID_corrn = new TH2D( "hclusmeanID_corrn", "HCal Neutron Cluster Mean Time TW Corrected (Primary Block Reftime/ID, nblk>1), vs ID", 288, 0, 288, 1000, -50, 50 );
//   TH1D *hclusdiff = new TH1D( "hclusdiff", "HCal Cluster Seed - Block Diff, (Primary Block Reftime, nblk>1)", 1000, -50, 50 );  
//   TH2D *hclusdiffID = new TH2D( "hclusdiffID", "HCal Cluster Seed - Block Diff, (Primary Block Reftime, nblk>1), vs ID", 288, 0, 288, 1000, -50, 50 );
//   TH2D *hpvablkdiffID = new TH2D( "hpvablkdiffID", "HCal Cluster Seed - Adjacent Block Diff vs ID", 288, 0, 288, 400, -20, 20 );


//   //Loop over events
//   cout << "Looping over hydrogen data.." << endl;
//   cout << "Looping over hydrogen data.." << endl;
//   if(pass0) gROOT->ProcessLine( "gErrorIgnoreLevel = 6001" ); //Suppress error output to avoid undetectable sbs.hcal.clus_blk.again for pass0

//   //Loop over all hydrogen data
//   for( Int_t f=0; f<nfset_lh2[kIdx]; f++ ){
    
//     Int_t hfieldset = fset_lh2[kIdx][f]; //Check B field setting
//     bool tofready = tofreplay && hfieldset==tofB[kIdx]; //Check if TOF corrections can be applied for these data where only one field setting has been estimated per kinematic for TOF v ep. Will expand to all fields.

//     //Declare energy loss parameters for beam going through the target
//     Double_t pBeam = E_e_h[f]/(1.+E_e_h[f]/M_p*(1.-cos(BB_th_h[f])));
//     Double_t Eloss_outgoing = celldiameter/2.0/sin(BB_th_h[f]) * rho_tgt * dEdx_tgt; //Mean energy loss of the beam prior to the scattering, approximately 1 MeV, could correct further with raster position (likely negligable)
    
//     for( Long64_t nevent = 1; nevent <Nevents_h[f]; nevent++){

//       if ( nevent%10000==0 ) cout << "LH2 kinematic " << kine << " at field " << hfieldset << "% , entry: " << nevent << "/" << Nevents_h[f] << ". Total events gathered for calibration: " << TNEV_h << " \r";
//       cout.flush();

//       if ( nevent%100000==0 ) cout << "LH2 kinematic " << kine << " at field " << hfieldset << "% , entry: " << nevent << "/" << Nevents_h[f] << ". Total events gathered for calibration: " << TNEV_h << endl;
      
//       Ch[f]->GetEntry( elist_h[f]->GetEntry( nevent ) ); 

//       Double_t A[kNcell] = {0.0}; // Array to keep track of ADC values per cell. Outscope on each ev
//       Double_t A_oneblock[kNcell] = {0.0}; // Array to keep track of ADC values per cell for one block clusters only. Outscope on each ev
      
//       //Correct the beam energy with energy loss in target using vertex position
//       Double_t Eloss = (BBtr_vz_h[f][0]+l_tgt/2.0) * rho_tgt * dEdx_tgt + uwallthick_LH2 * rho_Al * dEdx_Al; //approximately 3 MeV
//       Double_t E_corr = E_e_h[f] - Eloss;
//       Double_t p_corr = BBtr_p_h[f][0] - Eloss_outgoing; //Neglecting the mass of e'
//       Double_t etheta = acos( BBtr_pz_h[f][0]/BBtr_p_h[f][0] );
//       Double_t ephi = atan2( BBtr_py_h[f][0], BBtr_px_h[f][0] );

//       TVector3 vertex(0,0,BBtr_vz_h[f][0]); // z location of vertex in hall coordinates
//       TLorentzVector Pbeam(0,0,E_corr,E_corr); //Mass of e negligable
//       TLorentzVector kprime(BBtr_px_h[f][0],BBtr_py_h[f][0],BBtr_pz_h[f][0],BBtr_p_h[f][0]);
//       TLorentzVector Ptarg(0,0,0,M_p); // assume proton for both LH2 and LD2 - can refine where useful with exclusive LH2 data set. Likely better to refine first with dxdy spot cuts on both protons and neutrons.
//       TLorentzVector q = Pbeam - kprime;
//       TLorentzVector PgammaN = Ptarg + q; //(-px, -py, ebeam - pz, Mp + ebeam - p)

//       Double_t pel = E_corr/(1.+E_corr/M_p*(1.-cos(etheta)));
//       Double_t nu = E_corr - BBtr_p_h[f][0];
//       Double_t pp = sqrt(pow(nu,2)+2.*M_p*nu); //momentum of proton
//       Double_t pn = sqrt(pow(nu,2)+2.*M_n*nu); //momentum of neutron
//       Double_t phinucleon = ephi + PI; //assume coplanarity
//       Double_t thetanucleon = acos( (E_corr - BBtr_pz_h[f][0])/pp ); //use elastic constraint on nucleon kinematic
//       TVector3 pNhat( sin(thetanucleon)*cos(phinucleon),sin(thetanucleon)*sin(phinucleon),cos(thetanucleon));

//       //Define HCal coordinate system
//       TVector3 HCAL_zaxis(sin(-HCal_th_h[f]),0,cos(-HCal_th_h[f]));
//       TVector3 HCAL_xaxis(0,-1,0);
//       TVector3 HCAL_yaxis = HCAL_zaxis.Cross(HCAL_xaxis).Unit();	
//       TVector3 HCAL_origin = HCal_d_h[f] * HCAL_zaxis + hcalheight * HCAL_xaxis;

//       //Define intersection points for hadron vector
//       Double_t sintersect = ( HCAL_origin - vertex ).Dot( HCAL_zaxis ) / ( pNhat.Dot( HCAL_zaxis ) );
//       TVector3 HCAL_intersect = vertex + sintersect * pNhat;
//       Double_t yexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_yaxis );
//       Double_t xexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_xaxis );
//       Double_t E_ep = sqrt( pow(M_e,2) + pow(BBtr_p_h[f][0],2) ); // Obtain the scattered electron energy
//       Double_t p_ep = BBtr_p_h[f][0];
//       Double_t Q2 = 2*E_corr*E_ep*( 1-(BBtr_pz_h[f][0]/p_ep) ); // Obtain Q2 from beam energy, outgoing electron energy, and momenta
	
//       //Get some physics
//       Double_t W = PgammaN.M();
//       Double_t W2 = ekineW2_h[f];
//       Double_t E_pp = nu+M_p; // Get energy of the proton
//       Double_t Enucleon = sqrt(pow(pp,2)+pow(M_p,2)); // Check on E_pp, same
//       Double_t KE_p = nu; //For elastics

//       //Money
//       Double_t dx = HCALx_h[f] - xexpect_HCAL;
//       Double_t dy = HCALy_h[f] - yexpect_HCAL;

//       /////////////////////////////////////////////////
//       //Primary W2 cut on elastics and HCal active area
//       if( fabs(W2-W2_mean_h[f])>W2_sig_h[f] ) continue; //Observed mean W2 cut on elastic peak
//       elasYield++; //Events that pass the above cuts constitute elastics
//       /////////////////////////////////////////////////

//       //Calculate/declare new variables for analysis
//       Double_t SFrac = HCALe_h[f]/KE_p;
//       //Double_t E_exp = KE_p*sampFrac[kIdx];
//       Int_t rec_row = ( HCALx_h[f] - HCal_Xmin )/HCal_divx;
//       Int_t rec_col = ( HCALy_h[f] - HCal_Ymin )/HCal_divy;
//       Int_t rec_cell = rec_row*kNcols + rec_col;

//       hdx_mag_h->Fill( hfieldset, dx );
//       hdy_mag_h->Fill( hfieldset, dy );
      
//       //////////////////////////////////////////////////////////////////
//       //Cut on dx and dy.
//       bool pass_y = abs(dy-dy0_h[f])<ifac*dy_sig_h[f];
//       if( !pass_y ) continue;
//       bool pass_p = abs(dx-dx0_p_h[f])<ifac*dx_sig_p_h[f];
//       bool pass_n = abs(dx-dx0_n_h[f])<ifac*dx_sig_n_h[f]; //Redundant for LH2
//       bool isproton = pass_p; //pass_n is meaningless for LH2 (same cut)
//       bool isneutron = false; //No neutrons in LH2
//       bool isamb = pass_p && pass_n; //Will never obtain here
//       if( !pass_p && !pass_n ) continue; //Cut on both n and p spots for each event, cannot know which apriori
//       //////////////////////////////////////////////////////////////////
//       if( isproton ) hpp_p->Fill( pp );
//       if( isneutron ) hpp_n->Fill( pn );
      
//       //Calculate expected block from BB projections to HCAL block for TOF corrections
//       Double_t xdeflect = dx0_n_d[f] - dx0_p_d[f]; //Get SBS field deflection from fits to p/n peaks, ld2
//       Double_t xpexpect = (HCAL_intersect - HCAL_origin).Dot( HCAL_xaxis )-xdeflect;
//       Int_t IDexpect;

//       Int_t rowexpect = (xpexpect-HCal_Xi)/HCal_divx;
//       Int_t colexpect = (yexpect_HCAL-HCal_Yi)/HCal_divy;

//       IDexpect = 12*rowexpect+colexpect;

//       if( HCALy_h[f]<HCal_Yi || HCALy_h[f]>HCal_Yf || HCALx_h[f]<HCal_Xi || HCALx_h[f]>HCal_Xf ){
// 	cout << "Warning: Event passes p/n cut, but is projected off of HCal!" << endl;
// 	continue;
//       }
      
//       //Get time of flight correction by id from fits to MC
//       Double_t tofidcorr = 0;
//       if( isproton && tofready ) tofidcorr = tofcorr_p[f][IDexpect];
//       if( isneutron && tofready ) tofidcorr = tofcorr_n[f][IDexpect];

//       //Get time of flight correction by ep from functional fit to TOF vs mag_p
//       Double_t tofpcorr = 0.; //time of flight correction from momentum fit
//       if( isproton && tofready ){
// 	Double_t protmom = pp;
// 	if( protmom>pulim[kIdx] ) protmom = pulim[kIdx]; //do not pass corrections for momenta > fit limit in MC
// 	tofpcorr = TOFfitp_p[kIdx][0]+TOFfitp_p[kIdx][1]*protmom+TOFfitp_p[kIdx][2]*pow(protmom,2)+TOFfitp_p[kIdx][3]*pow(protmom,3);
//       }
//       if( isneutron && tofready ){
// 	Double_t neutmom = pn;
// 	if( neutmom>pulim[kIdx] ) neutmom = pulim[kIdx]; //do not pass corrections for momenta > fit limit in MC
// 	tofpcorr = TOFfitp_n[kIdx][0]+TOFfitp_n[kIdx][1]*neutmom+TOFfitp_n[kIdx][2]*pow(neutmom,2)+TOFfitp_n[kIdx][3]*pow(neutmom,3);
//       }

//       //Hodo cluster mean time
//       Double_t hodot = HODOtmean_h[f];

//       //Fill some histograms with only the primary block timing
//       Int_t pblkid = int(cblkid_h[f][0])-1;
//       Double_t pTDC = cblktime_h[f][0];
//       Double_t pADCt = cblkatime_h[f][0];
//       Double_t pTWcorr = 0.; //primary block timewalk correction
//       Double_t pADCtcorr = 0.;
//       Double_t pblke = 0.;
//       if( pass0 ){
// 	pblke = cblke_h[f][0]/gOldConst_pass0[pblkid]*gConst_iter1[pblkid];
//       }else{
// 	pblke = cblke_h[f][0]/cblkagain_h[f][0]*gConst_iter1[pblkid];
//       }

//       //Fill TDC histos.
//       if( qreplay ){
// 	pTWcorr = otdcP0[pblkid]*exp(-otdcP1[pblkid]*pblke)+otdcP2[pblkid];
// 	pTDC = pTDC + oldTDCoffsets[pblkid]*TDCCalib - calTDCoffsets[pblkid]*TDCCalib;
//       }
//       htp_ID->Fill( pblkid, pTDC );
//       htp_hodocorr_ID->Fill( pblkid, pTDC-hodot );
//       // Normalize to those corrected to evaluate efficiency
//       if( pTWcorr != 0. ) htp_twcorr_ID->Fill( pblkid, pTDC-pTWcorr );
//       if( isproton ){
// 	if( tofpcorr != 0. ){
// 	  htp_tofpcorr_p_ID->Fill( pblkid, pTDC-tofpcorr );
// 	  htp_allcorr_p->Fill( pTDC-hodot-tofpcorr );
// 	  htp_allcorr_p_ID->Fill( pblkid, pTDC-hodot-tofpcorr );
// 	}
// 	if( tofidcorr != 0. ) htp_tofidcorr_p_ID->Fill( pblkid, pTDC-tofidcorr );
//       }
//       if( isneutron ){
// 	if( tofpcorr != 0. ){
// 	  htp_tofpcorr_n_ID->Fill( pblkid, pTDC-tofpcorr );
// 	  htp_allcorr_n->Fill( pTDC-hodot-tofpcorr );
// 	  htp_allcorr_n_ID->Fill( pblkid, pTDC-hodot-tofpcorr );
// 	}
// 	if( tofidcorr != 0. ) htp_tofidcorr_n_ID->Fill( pblkid, pTDC-tofidcorr );
//       }

//       //Fill ADCt histos
//       if( qreplay ){
// 	pADCtcorr = oadctP0[pblkid]*exp(-oadctP1[pblkid]*pblke);
// 	pADCt = pADCt + oldADCtoffsets[pblkid] - calADCtoffsets[pblkid];
//       }
//       hap_ID->Fill( pblkid, pADCt );
//       hapDiff_ID->Fill( pblkid, pADCt-hodot );
//       hapCorr_ID->Fill( pblkid, pADCt-hodot-pADCtcorr );

//       //Loop over primary cluster
//       Double_t blkseed = 0.;
//       Double_t tavg = 0.;
//       Double_t blkseedcorr = 0.;
//       Double_t tavgcorr = 0.;
//       Int_t Ndiffcorr = 0; //keep track of the number of differences for meantime calculation rejecting blks where no good TW correction exists (TOF passed exclusively by BB and primary projected hadron)
//       Int_t nblk = (int)nblk_h[f];
//       for( Int_t blk = 0; blk<nblk; blk++ ){
// 	Int_t blkid = int(cblkid_h[f][blk])-1; //-1 necessary since sbs.hcal.clus_blk.id ranges from 1 to 288 (different than just about everything else)
// 	Double_t blkatime = cblkatime_h[f][blk];
// 	Double_t blktime = cblktime_h[f][blk];
// 	Double_t blke = 0.;
// 	if( pass0 ){
// 	  blke = cblke_h[f][blk]/gOldConst_pass0[blkid]*gConst_iter1[blkid];
// 	}else{
// 	  blke = cblke_h[f][blk]/cblkagain_h[f][blk]*gConst_iter1[blkid];
// 	}
	
// 	Double_t TWcorr = 0.; //tdc timewalk correction applied by block
// 	Double_t ADCtcorr = 0.; //adct tw corr
// 	if( qreplay ){
// 	  TWcorr = otdcP0[blkid]*exp(-otdcP1[blkid]*blke)+otdcP2[blkid]; //timewalk correction applied
// 	  ADCtcorr = oadctP0[blkid]*exp(-oadctP1[blkid]*blke);
// 	  Double_t blkatime_new = blkatime + oldADCtoffsets[blkid] - calADCtoffsets[blkid]; //reverse replay "+", then qreplay "-"
// 	  Double_t blktime_new = blktime + oldTDCoffsets[blkid]*TDCCalib - calTDCoffsets[blkid]*TDCCalib; //reverse replay "+", then qreplay "-", but first note that both offset values are in tdc units and need conversion to ns with tdccalib
// 	  Double_t blktime_new_corrected = blktime_new - hodot - TWcorr - tofpcorr; //Use p for now

// 	  //Fill ADCt histos
// 	  if( blkatime_new>0 ){
// 	    hadct->Fill( blkatime_new );
// 	    hadct_corr->Fill( blkatime_new - hodot - ADCtcorr );
// 	    ha_ID->Fill( blkid, blkatime_new );
// 	    haDiff_ID->Fill( blkid, blkatime_new - hodot );
// 	    haCorr_ID->Fill( blkid, blkatime_new - hodot - ADCtcorr );
// 	    hadctVe[blkid]->Fill( blke*1000, blkatime_new ); //Convert to MeV for clarity in plots
// 	  }
// 	  //Fill TDC histos
// 	  if( blktime_new<10000 && blktime_new>-500){ //ensure that good tdc times exist on tree
// 	    htdc->Fill( blktime_new );
// 	    ht_ID->Fill( blkid, blktime_new );
// 	    htdc_hodocorr->Fill( blktime_new - hodot );
// 	    ht_hodocorr_ID->Fill( blkid, blktime_new - hodot );
// 	    htdcVe[blkid]->Fill( blke*1000, blktime_new ); //Convert to MeV for clarity in plots
// 	    if( TWcorr != 0. ){
// 	      htdc_twcorr->Fill( blktime_new - TWcorr );
// 	      ht_twcorr_ID->Fill( blkid, blktime_new - TWcorr );
// 	    }
// 	    if( isproton ){
// 	      if( tofpcorr != 0. ){
// 		htdc_tofpcorr_p->Fill( blktime_new - tofpcorr );
// 		ht_tofpcorr_p_ID->Fill( blkid, blktime_new - tofpcorr );
// 	      }
// 	      if( tofidcorr != 0. ){
// 		htdc_tofidcorr_p->Fill( blktime_new - tofidcorr );
// 		ht_tofidcorr_p_ID->Fill( blkid, blktime_new - tofidcorr );
// 	      }
// 	      if( TWcorr != 0. && tofpcorr != 0. ){
// 		htdc_allcorr_p->Fill( blktime_new_corrected );
// 		ht_allcorr_p_ID->Fill( blkid, blktime_new_corrected );
// 	      }
// 	    }
// 	    if( isneutron ){
// 	      if( tofpcorr != 0. ){
// 		htdc_tofpcorr_n->Fill( blktime_new - tofpcorr );
// 		ht_tofpcorr_n_ID->Fill( blkid, blktime_new - tofpcorr );
// 	      }
// 	      if( tofidcorr != 0. ){
// 		htdc_tofidcorr_n->Fill( blktime_new - tofidcorr );
// 		ht_tofidcorr_n_ID->Fill( blkid, blktime_new - tofidcorr );
// 	      }
// 	      if( TWcorr != 0. && tofpcorr != 0. ){
// 		htdc_allcorr_n->Fill( blktime_new_corrected );
// 		ht_allcorr_n_ID->Fill( blkid, blktime_new_corrected );
// 	      }
// 	    }

// 	    //Get self timing histograms
// 	    Double_t blktime_nTW = blktime_new - TWcorr; //Include only correction relevent for self timing
// 	    if( blk==0 ){ //Just set the reference time with the time of the primary block
// 	      blkseed = blktime_new;
// 	      if( TWcorr != 0. ) blkseedcorr = blktime_nTW;
// 	    }else if( blkseed != 0. ){
	      
// 	      Double_t tcdiff = blktime_new - blkseed;
// 	      tavg += tcdiff;

// 	      if( TWcorr != 0. ){
// 		Double_t tcdiffcorr = blktime_nTW - blkseedcorr;
// 		tavgcorr += tcdiffcorr;
// 		Ndiffcorr++;
// 	      }

// 	      if( blkid == pblkid+1 ) hpvablkdiffID->Fill( pblkid, tcdiff );

// 	      hclusdiffID->Fill( blkid, tcdiff );
// 	      hclusdiff->Fill( tcdiff );

// 	    }
	    
// 	  }
// 	}else{
// 	  if( blkatime>0 ){
// 	    hadct->Fill( blkatime );
// 	    hadct_corr->Fill( blkatime - hodot - ADCtcorr );
// 	    ha_ID->Fill( blkid, blkatime );
// 	    haDiff_ID->Fill( blkid, blkatime - hodot );
// 	    haCorr_ID->Fill( blkid, blkatime - hodot - ADCtcorr );
// 	    hadctVe[blkid]->Fill( blke*1000, blkatime ); //Convert to MeV for clarity in plots
// 	  }
// 	  if( blktime<10000 && blktime>-500){ //ensure that good tdc times exist on tree
// 	    htdc->Fill( blktime );
// 	    ht_ID->Fill( blkid, blktime );
// 	    htdc_hodocorr->Fill( blktime - hodot );
// 	    ht_hodocorr_ID->Fill( blkid, blktime - hodot );
// 	    htdcVe[blkid]->Fill( blke*1000, blktime ); //Convert to MeV for clarity in plots
	    
// 	    //Get self timing histograms (no recovery on bad tdc pblks yet)
// 	    if( blk==0 ){ //Just set the reference time with the time of the primary block
// 	      blkseed = blktime;
// 	    }else if( blkseed != 0. ){
	      
// 	      Double_t tcdiff = blktime-blkseed;
	      
// 	      tavg += tcdiff;
// 	      hclusdiffID->Fill( blkid, tcdiff );
// 	      hclusdiff->Fill( tcdiff );
// 	    }
// 	  }
// 	}
	
// 	NEV[blkid]++;
// 	TNEV_h++;

//       } //loop over blocks in primary cluster

//       //Finish getting self timing mean time
//       if( nblk>1 ){
// 	tavg /= (nblk-1);
      
// 	if( tavg!=0 ){
// 	  hclusmeanID->Fill( pblkid, tavg );
// 	  hclusmean->Fill( tavg );
// 	}
//       }
//       if( Ndiffcorr>1 ){
// 	tavgcorr /= Ndiffcorr;
      
// 	if( tavgcorr!=0 ){
// 	  if( isproton ){
// 	    hclusmeanID_corrp->Fill( pblkid, tavgcorr );
// 	    hclusmean_corrp->Fill( tavgcorr );
// 	  }
// 	  if( isneutron ){
// 	    hclusmeanID_corrn->Fill( pblkid, tavgcorr );
// 	    hclusmean_corrn->Fill( tavgcorr );
// 	  }
// 	}

//       }

//     } //loop over hydrogen events
//   } //loop for lh2

//   if( !qreplay ){
//     Int_t cell = 0;

//     cout << endl << "Number of events available for calibration from hydrogen alone: " << endl;
//     cout << endl << "Number of events available for calibration from hydrogen alone: " << endl;

//     for( Int_t r = 0; r<kNrows; r++){
//       for( Int_t c = 0; c<kNcols; c++){
// 	cout << NEV[cell] << "  ";
// 	cout << NEV[cell] << "  ";
// 	cell++;
//       }
//       cout << endl;
//       cout << endl;
//     }
  
//     cout << endl;
//     cout << endl;

//     usleep( 3*second ); //Give some time for review of step
//   }

//   cout << "Looping over deuterium data.." << endl;
//   cout << "Looping over deuterium data.." << endl;

//   //Loop over all deuterium data
//   for( Int_t f=0; f<nfset_ld2[kIdx]; f++ ){

//     Int_t dfieldset = fset_ld2[kIdx][f];
//     bool tofready = tofreplay && dfieldset==tofB[kIdx]; //Check if TOF corrections can be applied for these data

//     //Declare energy loss parameters for beam going through the target
//     Double_t pBeam = E_e_d[f]/(1.+E_e_d[f]/M_p*(1.-cos(BB_th_d[f])));
//     Double_t Eloss_outgoing = celldiameter/2.0/sin(BB_th_d[f]) * rho_tgt * dEdx_tgt; //Mean energy loss of the beam prior to the scattering, approximately 1 MeV, could correct further with raster position (likely negligable)

//     for( Long64_t nevent = 1; nevent <Nevents_d[f]; nevent++){

//       if ( nevent%10000==0 ) cout << "LD2 kinematic " << kine << " at field " << dfieldset << "%, entry: " << nevent << "/" << Nevents_d[f] << ". Total events gathered for calibration: " << TNEV_d << " \r";
//       cout.flush();
      
//       if ( nevent%100000==0 ) cout << "LD2 kinematic " << kine << " at field " << dfieldset << "%, entry: " << nevent << "/" << Nevents_d[f] << ". Total events gathered for calibration: " << TNEV_d << endl;

//       Cd[f]->GetEntry( elist_d[f]->GetEntry( nevent ) ); 

//       Double_t A[kNcell] = {0.0}; // Array to keep track of ADC values per cell. Outscope on each ev
//       Double_t A_oneblock[kNcell] = {0.0}; // Array to keep track of ADC values per cell for one block clusters only. Outscope on each ev
         
//       //Correct the beam energy with energy loss in target using vertex position
//       Double_t Eloss = (BBtr_vz_d[f][0]+l_tgt/2.0) * rho_tgt * dEdx_tgt + uwallthick_LH2 * rho_Al * dEdx_Al; //approximately 3 MeV
//       Double_t E_corr = E_e_d[f] - Eloss;
//       Double_t p_corr = BBtr_p_d[f][0] - Eloss_outgoing; //Neglecting the mass of e'
//       Double_t etheta = acos( BBtr_pz_d[f][0]/BBtr_p_d[f][0] );
//       Double_t ephi = atan2( BBtr_py_d[f][0], BBtr_px_d[f][0] );

//       TVector3 vertex(0,0,BBtr_vz_d[f][0]); // z location of vertex in hall coordinates
//       TLorentzVector Pbeam(0,0,E_corr,E_corr); //Mass of e negligable
//       TLorentzVector kprime(BBtr_px_d[f][0],BBtr_py_d[f][0],BBtr_pz_d[f][0],BBtr_p_d[f][0]);
//       TLorentzVector Ptarg(0,0,0,M_p); // assume proton for both LH2 and LD2 - can refine where useful with exclusive LH2 data set. Likely better to refine first with dxdy spot cuts on both protons and neutrons.
//       TLorentzVector q = Pbeam - kprime;
//       TLorentzVector PgammaN = Ptarg + q; //(-px, -py, ebeam - pz, Mp + ebeam - p)

//       Double_t pel = E_corr/(1.+E_corr/M_p*(1.-cos(etheta)));
//       Double_t nu = E_corr - BBtr_p_d[f][0];
//       Double_t pp = sqrt(pow(nu,2)+2.*M_p*nu); //momentum of proton
//       Double_t pn = sqrt(pow(nu,2)+2.*M_n*nu); //momentum of neutron
//       Double_t phinucleon = ephi + PI; //assume coplanarity
//       Double_t thetanucleon = acos( (E_corr - BBtr_pz_d[f][0])/pp );	
//       TVector3 pNhat( sin(thetanucleon)*cos(phinucleon),sin(thetanucleon)*sin(phinucleon),cos(thetanucleon));

//       //Define HCal coordinate system
//       TVector3 HCAL_zaxis(sin(-HCal_th_d[f]),0,cos(-HCal_th_d[f]));
//       TVector3 HCAL_xaxis(0,-1,0);
//       TVector3 HCAL_yaxis = HCAL_zaxis.Cross(HCAL_xaxis).Unit();	
//       TVector3 HCAL_origin = HCal_d_d[f] * HCAL_zaxis + hcalheight * HCAL_xaxis;

//       //Define intersection points for hadron vector
//       Double_t sintersect = ( HCAL_origin - vertex ).Dot( HCAL_zaxis ) / ( pNhat.Dot( HCAL_zaxis ) );
//       TVector3 HCAL_intersect = vertex + sintersect * pNhat;
//       Double_t yexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_yaxis );
//       Double_t xexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_xaxis );
//       Double_t E_ep = sqrt( pow(M_e,2) + pow(BBtr_p_d[f][0],2) ); // Obtain the scattered electron energy	
//       Double_t p_ep = BBtr_p_d[f][0];
//       Double_t Q2 = 2*E_corr*E_ep*( 1-(BBtr_pz_d[f][0]/p_ep) ); // Obtain Q2 from beam energy, outgoing electron energy, and momenta
	
//       //Get some physics
//       Double_t W = PgammaN.M();
//       Double_t W2 = ekineW2_d[f];
//       Double_t E_pp = nu+M_p; // Get energy of the proton
//       Double_t Enucleon = sqrt(pow(pp,2)+pow(M_p,2)); // Check on E_pp, same
//       Double_t KE_p = nu; //For elastics

//       //Money
//       Double_t dx = HCALx_d[f] - xexpect_HCAL;
//       Double_t dy = HCALy_d[f] - yexpect_HCAL;

//       /////////////////////////////////////////////////
//       //Primary W2 cut on elastics and HCal active area
//       if( fabs(W2-W2_mean_d[f])>W2_sig_d[f] ) continue; //Observed mean W2 cut on elastic peak
//       elasYield++; //Events that pass the above cuts constitute elastics
//       /////////////////////////////////////////////////

//       //Calculate/declare new variables for analysis
//       Double_t SFrac = HCALe_d[f]/KE_p;
//       //Double_t E_exp = KE_p*sampFrac[kIdx];
//       Int_t rec_row = ( HCALx_d[f] - HCal_Xmin )/HCal_divx;
//       Int_t rec_col = ( HCALy_d[f] - HCal_Ymin )/HCal_divy;
//       Int_t rec_cell = rec_row*kNcols + rec_col;
      
//       hdx_mag_d->Fill( dfieldset, dx );
//       hdy_mag_d->Fill( dfieldset, dy );

//       //////////////////////////////////////////////////////////////////
//       //Cut on dx and dy.
//       bool pass_y = abs(dy-dy0_d[f])<ifac*dy_sig_d[f];
//       if( !pass_y ) continue;
//       bool pass_p = abs(dx-dx0_p_d[f])<ifac*dx_sig_p_d[f];
//       bool pass_n = abs(dx-dx0_n_d[f])<ifac*dx_sig_n_d[f];
//       bool isproton = pass_p && !pass_n;
//       bool isneutron = pass_n && !pass_p;
//       bool isamb = pass_p && pass_n;
//       if( !pass_p && !pass_n ) continue; //Cut on both n and p spots for each event, cannot know which apriori      
//       //////////////////////////////////////////////////////////////////
//       if( isproton ) hpp_p->Fill( pp );
//       if( isneutron ) hpp_n->Fill( pn );

//       //Calculate expected block from BB projections to HCAL block for TOF corrections
//       Double_t xdeflect = dx0_n_d[f] - dx0_p_d[f]; //Get SBS field deflection from fits to p/n peaks, ld2
//       Double_t xpexpect = (HCAL_intersect - HCAL_origin).Dot( HCAL_xaxis )-xdeflect;
//       Int_t IDexpect;

//       Int_t rowexpect = (xpexpect-HCal_Xi)/HCal_divx;
//       Int_t colexpect = (yexpect_HCAL-HCal_Yi)/HCal_divy;

//       IDexpect = 12*rowexpect+colexpect;

//       if( HCALy_d[f]<HCal_Yi || HCALy_d[f]>HCal_Yf || HCALx_d[f]<HCal_Xi || HCALx_d[f]>HCal_Xf ){
// 	cout << "Warning: Event passes p/n cut, but is projected off of HCal!" << endl;
// 	continue;
//       }
      
//       //Get time of flight correction by id from fits to MC
//       Double_t tofidcorr = 0;
//       if( isproton && tofready ) tofidcorr = tofcorr_p[f][IDexpect];
//       if( isneutron && tofready ) tofidcorr = tofcorr_n[f][IDexpect];

//       //Get time of flight correction by ep from functional fit to TOF vs mag_p
//       Double_t tofpcorr = 0.; //time of flight correction from momentum fit
//       if( isproton && tofready ){
// 	Double_t protmom = pp;
// 	if( protmom>pulim[kIdx] ) protmom = pulim[kIdx]; //do not pass corrections for momenta > fit limit in MC
// 	tofpcorr = TOFfitp_p[kIdx][0]+TOFfitp_p[kIdx][1]*protmom+TOFfitp_p[kIdx][2]*pow(protmom,2)+TOFfitp_p[kIdx][3]*pow(protmom,3);
//       }
//       if( isneutron && tofready ){
// 	Double_t neutmom = pn;
// 	if( neutmom>pulim[kIdx] ) neutmom = pulim[kIdx]; //do not pass corrections for momenta > fit limit in MC
// 	tofpcorr = TOFfitp_n[kIdx][0]+TOFfitp_n[kIdx][1]*neutmom+TOFfitp_n[kIdx][2]*pow(neutmom,2)+TOFfitp_n[kIdx][3]*pow(neutmom,3);
//       }

//       //Hodo cluster mean time
//       Double_t hodot = HODOtmean_d[f];

//       //Fill some histograms with only the primary block timing
//       Int_t pblkid = int(cblkid_d[f][0])-1;
//       Double_t pTDC = cblktime_d[f][0];
//       Double_t pADCt = cblkatime_d[f][0];
//       Double_t pTWcorr = 0.;
//       Double_t pADCtcorr = 0.;
//       Double_t pblke = 0.;
//       if( pass0 ){
// 	pblke = cblke_d[f][0]/gOldConst_pass0[pblkid]*gConst_iter1[pblkid];
//       }else{
// 	pblke = cblke_d[f][0]/cblkagain_d[f][0]*gConst_iter1[pblkid];
//       }

//       //Fill TDC histos.
//       if( qreplay ){
// 	pTWcorr = otdcP0[pblkid]*exp(-otdcP1[pblkid]*pblke)+otdcP2[pblkid];
// 	pTDC = pTDC + oldTDCoffsets[pblkid]*TDCCalib - calTDCoffsets[pblkid]*TDCCalib;
//       }
//       htp_ID->Fill( pblkid, pTDC );
//       htp_hodocorr_ID->Fill( pblkid, pTDC-hodot );
//       // Normalize to those corrected to evaluate efficiency
//       if( pTWcorr != 0. ) htp_twcorr_ID->Fill( pblkid, pTDC-pTWcorr );
//       if( isproton ){
// 	if( tofpcorr != 0. ){
// 	  htp_tofpcorr_p_ID->Fill( pblkid, pTDC-tofpcorr );
// 	  htp_allcorr_p_ID->Fill( pblkid, pTDC-hodot-tofpcorr );
// 	  htp_allcorr_p->Fill( pTDC-hodot-tofpcorr );
// 	}
// 	if( tofidcorr != 0. ) htp_tofidcorr_p_ID->Fill( pblkid, pTDC-tofidcorr );
//       }
//       if( isneutron ){
// 	if( tofpcorr != 0. ){
// 	  htp_tofpcorr_n_ID->Fill( pblkid, pTDC-tofpcorr );
// 	  htp_allcorr_n_ID->Fill( pblkid, pTDC-hodot-tofpcorr );
// 	  htp_allcorr_n->Fill( pTDC-hodot-tofpcorr );
// 	}
// 	if( tofidcorr != 0. ) htp_tofidcorr_n_ID->Fill( pblkid, pTDC-tofidcorr );
//       }

//       //Fill ADCt histos
//       if( qreplay ){
// 	pADCtcorr = oadctP0[pblkid]*exp(-oadctP1[pblkid]*pblke);
// 	pADCt = pADCt + oldADCtoffsets[pblkid] - calADCtoffsets[pblkid];
//       }
//       hap_ID->Fill( pblkid, pADCt );
//       hapDiff_ID->Fill( pblkid, pADCt-hodot );
//       hapCorr_ID->Fill( pblkid, pADCt-hodot-pADCtcorr );

//       //Loop over primary cluster
//       Double_t blkseed = 0.;
//       Double_t tavg = 0;
//       Double_t blkseedcorr = 0.;
//       Double_t tavgcorr = 0.;
//       Int_t Ndiffcorr = 0; //keep track of the number of differences for meantime calculation rejecting blks where no good TW correction exists (TOF passed exclusively by BB and primary projected hadron)
//       Int_t nblk = (int)nblk_d[f];
//       for( Int_t blk = 0; blk<nblk; blk++ ){
// 	Int_t blkid = int(cblkid_d[f][blk])-1; //-1 necessary since sbs.hcal.clus_blk.id ranges from 1 to 288 (different than just about everything else)
// 	Double_t blkatime = cblkatime_d[f][blk];
// 	Double_t blktime = cblktime_d[f][blk];
// 	Double_t blke = 0.;
// 	if( pass0 ){
// 	  blke = cblke_d[f][blk]/gOldConst_pass0[blkid]*gConst_iter1[blkid];
// 	}else{
// 	  blke = cblke_d[f][blk]/cblkagain_d[f][blk]*gConst_iter1[blkid];
// 	}

// 	Double_t TWcorr = 0.; //tdc timewalk correction applied by block
// 	Double_t ADCtcorr = 0.; //adct tw corr
// 	if( qreplay ){
// 	  TWcorr = otdcP0[blkid]*exp(-otdcP1[blkid]*blke)+otdcP2[blkid]; //timewalk correction applied
// 	  ADCtcorr = oadctP0[blkid]*exp(-oadctP1[blkid]*blke);
// 	  Double_t blkatime_new = blkatime + oldADCtoffsets[blkid] - calADCtoffsets[blkid]; //reverse replay "+", then qreplay "-"
// 	  Double_t blktime_new = blktime + oldTDCoffsets[blkid]*TDCCalib - calTDCoffsets[blkid]*TDCCalib; //reverse replay "+", then qreplay "-", but first note that both offset values are in tdc units and need conversion to ns with tdccalib
// 	  Double_t blktime_new_corrected = blktime_new - hodot - TWcorr - tofpcorr; //Use p for now

// 	  //Fill ADCt histos
// 	  if( blkatime_new>0 ){
// 	    hadct->Fill( blkatime_new );
// 	    hadct_corr->Fill( blkatime_new - hodot - ADCtcorr );
// 	    ha_ID->Fill( blkid, blkatime_new );
// 	    haDiff_ID->Fill( blkid, blkatime_new - hodot );
// 	    haCorr_ID->Fill( blkid, blkatime_new - hodot - ADCtcorr );
// 	    hadctVe[blkid]->Fill( blke*1000, blkatime_new ); //Convert to MeV for clarity in plots
// 	  }
// 	  //Fill TDC histos
// 	  if( blktime_new<10000 && blktime_new>-500){ //ensure that good tdc times exist on tree
// 	    htdc->Fill( blktime_new );
// 	    ht_ID->Fill( blkid, blktime_new );
// 	    htdc_hodocorr->Fill( blktime_new - hodot );
// 	    ht_hodocorr_ID->Fill( blkid, blktime_new - hodot );
// 	    if( TWcorr != 0. ){
// 	      htdc_twcorr->Fill( blktime_new - TWcorr );
// 	      ht_twcorr_ID->Fill( blkid, blktime_new - TWcorr );
// 	    }
// 	    if( isproton ){
// 	      if( tofpcorr != 0. ){
// 		htdc_tofpcorr_p->Fill( blktime_new - tofpcorr );
// 		ht_tofpcorr_p_ID->Fill( blkid, blktime_new - tofpcorr );
// 	      }
// 	      if( tofidcorr != 0. ){
// 		htdc_tofidcorr_p->Fill( blktime_new - tofidcorr );
// 		ht_tofidcorr_p_ID->Fill( blkid, blktime_new - tofidcorr );
// 	      }
// 	      if( TWcorr != 0. && tofpcorr != 0. ){
// 		htdc_allcorr_p->Fill( blktime_new_corrected );
// 		ht_allcorr_p_ID->Fill( blkid, blktime_new_corrected );
// 	      }
// 	    }
// 	    if( isneutron ){
// 	      if( tofpcorr != 0. ){
// 		htdc_tofpcorr_n->Fill( blktime_new - tofpcorr );
// 		ht_tofpcorr_n_ID->Fill( blkid, blktime_new - tofpcorr );
// 	      }
// 	      if( tofidcorr != 0. ){
// 		htdc_tofidcorr_n->Fill( blktime_new - tofidcorr );
// 		ht_tofidcorr_n_ID->Fill( blkid, blktime_new - tofidcorr );
// 	      }
// 	      if( TWcorr != 0. && tofpcorr != 0. ){
// 		htdc_allcorr_n->Fill( blktime_new_corrected );
// 		ht_allcorr_n_ID->Fill( blkid, blktime_new_corrected );
// 	      }
// 	    }

// 	    //Get self timing histograms
// 	    Double_t blktime_nTW = blktime_new - TWcorr; //Include only correction relevent for self timing
// 	    if( blk==0 ){ //Just set the reference time with the time of the primary block
// 	      blkseed = blktime_new;
// 	      if( TWcorr != 0. ) blkseedcorr = blktime_nTW;
// 	    }else if( blkseed != 0. ){
	      
// 	      Double_t tcdiff = blktime_new - blkseed;
// 	      tavg += tcdiff;

// 	      if( TWcorr != 0. && tofpcorr != 0. ){
// 		Double_t tcdiffcorr = blktime_nTW - blkseedcorr;
// 		tavgcorr += tcdiffcorr;
// 		Ndiffcorr++;
// 	      }

// 	      if( blkid == pblkid+1 ) hpvablkdiffID->Fill( pblkid, tcdiff );

// 	      hclusdiffID->Fill( blkid, tcdiff );
// 	      hclusdiff->Fill( tcdiff );

// 	    }
	 
// 	  }	  
// 	}else{
// 	  if( blkatime>0 ){
// 	    hadct->Fill( blkatime );
// 	    hadct_corr->Fill( blkatime - hodot - ADCtcorr );
// 	    ha_ID->Fill( blkid, blkatime );
// 	    haDiff_ID->Fill( blkid, blkatime - hodot );
// 	    haCorr_ID->Fill( blkid, blkatime - hodot - ADCtcorr );
// 	    hadctVe[blkid]->Fill( blke*1000, blkatime ); //Convert to MeV for clarity in plots
// 	  }
// 	  if( blktime<10000 && blktime>-500){ //ensure that good tdc times exist on tree
// 	    htdc->Fill( blktime );
// 	    ht_ID->Fill( blkid, blktime );
// 	    htdc_hodocorr->Fill( blktime - hodot );
// 	    ht_hodocorr_ID->Fill( blkid, blktime - hodot );
// 	    htdcVe[blkid]->Fill( blke*1000, blktime ); //Convert to MeV for clarity in plots
	     
// 	    //Get self timing histograms
// 	    if( blk==0 ){ //Just set the reference time with the time of the primary block
// 	      blkseed = blktime;
// 	    }else if( blkseed != 0. ){
	      
// 	      Double_t tcdiff = blktime-blkseed;

// 	      tavg += tcdiff;
// 	      hclusdiffID->Fill( blkid, tcdiff );
// 	      hclusdiff->Fill( tcdiff );
// 	    }
// 	  }
// 	}

// 	NEV[blkid]++;
// 	TNEV_d++;
	
//       } //loop over blocks in primary cluster

//       //Finish getting self timing mean time
//       if( nblk>1 ){
// 	tavg /= (nblk-1);
      
// 	if( tavg!=0 ){
// 	  hclusmeanID->Fill( pblkid, tavg );
// 	  hclusmean->Fill( tavg );
// 	}
//       }
//       if( Ndiffcorr>1 ){
// 	tavgcorr /= Ndiffcorr;
      
// 	if( tavgcorr!=0 ){
// 	  if( isproton ){
// 	    hclusmeanID_corrp->Fill( pblkid, tavgcorr );
// 	    hclusmean_corrp->Fill( tavgcorr );
// 	  }
// 	  if( isneutron ){
// 	    hclusmeanID_corrn->Fill( pblkid, tavgcorr );
// 	    hclusmean_corrn->Fill( tavgcorr );
// 	  }
// 	}

//       }

//     } //loop over deuterium events
//   } //loop for ld2

//   cout << endl << "Fitting TDC/ADCt/TW distributions" << endl;
//   cout << endl << "Fitting TDC/ADCt/TW distributions" << endl;

//   //Make arrays for TW fits
//   Double_t TDCvseP0[kNcell] = {0.0};
//   Double_t TDCvseP1[kNcell] = {0.0};
//   Double_t TDCvseP2[kNcell] = {0.0};
//   Double_t ADCtvseP0[kNcell] = {0.0};
//   Double_t ADCtvseP1[kNcell] = {0.0};
//   Double_t ADCtvseP2[kNcell] = {0.0};

//   //Set fitting parameters for gaussian fits
//   Double_t asetpar[3];
//   Double_t tsetpar[3];

//   //No error on cell location
//   Double_t cellerr[kNcell] = {0.};

//   //Make arrays for tdc tgraphs
//   Double_t tcell[kNcell] = {0.};
//   Double_t tcval[kNcell] = {0.};
//   Double_t tcvalw[kNcell] = {0.};
//   Double_t tcerr[kNcell] = {0.};
//   TH1D *tcellslice[kNcell];

//   //Make arrays for adc tgraphs
//   Double_t acell[kNcell] = {0.};
//   Double_t acval[kNcell] = {0.};
//   Double_t acvalw[kNcell] = {0.};
//   Double_t acerr[kNcell] = {0.};
//   TH1D *acellslice[kNcell];

//   //Get averages
//   Double_t tcval_avg = 0.;
//   Int_t tcval_Ng = 0;
//   Double_t acval_avg = 0.;
//   Int_t acval_Ng = 0;
  
//   //Fits for timewalk corrections
//   for(Int_t c=0; c<kNcell; c++){
//     //Fit the TDC vs E plots
//     TF1 *fittdcTW = new TF1( "fittdcTW", TW_fit, 0, 300, 3 );
//     fittdcTW->SetParameters(14,0.04,-77);
//     fittdcTW->SetParLimits(0,2,26);
//     fittdcTW->SetParLimits(1,0.01,0.2);
//     fittdcTW->SetParLimits(2,-200,50);

//     if( htdcVe[c]->GetEntries()>tFitMin ){
//       htdcVe[c]->Fit("fittdcTW","Q","",5,300);
//       TDCvseP0[c] = fittdcTW->GetParameter(0);
//       TDCvseP1[c] = fittdcTW->GetParameter(1);
//       TDCvseP2[c] = fittdcTW->GetParameter(2);
//       htdcVe[c]->SetTitle(Form("P0:%f P1:%f P2:%f",TDCvseP0[c],TDCvseP1[c],TDCvseP2[c]));
//     }
//     htdcP0->Fill(TDCvseP0[c]);
//     htdcP1->Fill(TDCvseP1[c]);
//     htdcP2->Fill(TDCvseP2[c]);
//   }

//   //Separate loops to improve fits (avoid parameter recall)
//   for(Int_t c=0; c<kNcell; c++){
//     //Fit the ADCt vs E plots
//     TF1 *fitadctTW = new TF1( "fitadctTW", TW_fit, 0, 300, 3 );
//     fitadctTW->SetParameters(10,0.01,75);
//     fitadctTW->SetParLimits(0,2,9);
//     fitadctTW->SetParLimits(1,0.008,0.05);
//     fitadctTW->SetParLimits(2,-50,200);

//     if( hadctVe[c]->GetEntries()>tFitMin ){
//       hadctVe[c]->Fit("fitadctTW","Q","",5,300);
//       ADCtvseP0[c] = fitadctTW->GetParameter(0);
//       ADCtvseP1[c] = fitadctTW->GetParameter(1);
//       ADCtvseP2[c] = fitadctTW->GetParameter(2);
//       hadctVe[c]->SetTitle(Form("P0:%f P1:%f P2:%f",ADCtvseP0[c],ADCtvseP1[c],ADCtvseP2[c]));
//     }
//     hadctP0->Fill(ADCtvseP0[c]);
//     hadctP1->Fill(ADCtvseP1[c]);
//     hadctP2->Fill(ADCtvseP2[c]);
//   }

//   //Fits for TDC time
//   TCanvas *TDC_top = new TCanvas("TDC_top","TDC_top",1600,1200);
//   TCanvas *TDC_bot = new TCanvas("TDC_bot","TDC_bot",1600,1200);

//   TDC_top->Divide(12,12);
//   TDC_bot->Divide(12,12);

//   gStyle->SetOptStat(0);

//   for(Int_t c=0; c<kNcell; c++){

//     //Index through the canvas
//     TDC_top->cd(c+1);
//     if( c>=144 ){
//       TDC_bot->cd(c-143);
//       gStyle->SetOptStat(0);
//     }

//     //Get slices from htDiff_ID and fit for mean vals
//     Double_t tfitl = 0.;
//     Double_t tfith = 0.;
//     tcell[c] = c;
//     if( !tofreplay ){
//       tcellslice[c] = htp_hodocorr_ID->ProjectionY(Form("tcellslice_%d",c+1),c+1,c+1); //Trying htp_hodocorr_ID from htDiff_ID
//     }else{
//       tcellslice[c] = htp_allcorr_p_ID->ProjectionY(Form("tcellslice_%d",c+1),c+1,c+1);
//     }
//     tcval[c] = oldTDCoffsets[c]; //will overwrite if fit is good.
//     tcvalw[c] = 0.; //leave as zero to evaluate fit values only

//     Int_t sliceN = tcellslice[c]->GetEntries();
//     if( sliceN<tFitMin ){
//       tcellslice[c]->Draw();
//       continue;
//     }

//     Double_t arimean = tcellslice[c]->GetMean();
//     tsetpar[0] = sliceN;
//     tsetpar[1] = arimean;
//     tsetpar[2] = aTDCsig;
//     tfitl = arimean - 4*aTDCsig;
//     tfith = arimean + 4*aTDCsig;
//     TF1 *gausfit = new TF1("gausfit",Gfit,tfitl,tfith,3);
//     gausfit->SetLineWidth(4);
//     gausfit->SetParameter(0,tsetpar[0]);
//     gausfit->SetParameter(1,tsetpar[1]);
//     gausfit->SetParLimits(1,tfitl,tfith);
//     gausfit->SetParameter(2,tsetpar[2]);

//     tcellslice[c]->Fit("gausfit","RBM");
//     tcellslice[c]->Draw();

//     tcval[c] = gausfit->GetParameter(1);
//     tcerr[c] = gausfit->GetParameter(2);
//     tcvalw[c] = gausfit->GetParameter(1);
//     tcellslice[c]->SetTitle(Form("Mean:%f Sigma:%f",tcval[c],tcerr[c]));    

//     tcval_avg += tcval[c];
//     tcval_Ng++;

//   }    
//   TDC_top->Write();
//   TDC_bot->Write();

//   tcval_avg /= tcval_Ng;

//   //Fits for ADC time 
//   TCanvas *ADCt_top = new TCanvas("ADCt_top","ADCt_top",1600,1200);
//   TCanvas *ADCt_bot = new TCanvas("ADCt_bot","ADCt_bot",1600,1200);

//   ADCt_top->Divide(12,12);
//   ADCt_bot->Divide(12,12);

//   gStyle->SetOptStat(0);

//   for(Int_t c=0; c<kNcell; c++){

//     Int_t col = c%kNcols;

//     //Index through the canvas
//     ADCt_top->cd(c+1);
//     if( c>=144 ){
//       ADCt_bot->cd(c-143);
//       gStyle->SetOptStat(0);
//     }

//     //Get slices from haDiff_ID and fit for mean vals
//     Double_t afitl = 0.; //50. empirically
//     Double_t afith = 0.; //75.
//     acell[c] = c;
//     acellslice[c] = haDiff_ID->ProjectionY(Form("acellslice_%d",c+1),c+1,c+1); //Trying hapDiff_ID from haDiff_ID
//     acval[c] = oldADCtoffsets[c]; //will overwrite if fit is good.
//     acvalw[c] = 0.; //will overwrite if fit is good.

//     Int_t sliceN = acellslice[c]->GetEntries();
//     if( sliceN<tFitMin ){
//       acellslice[c]->Draw();
//       continue;
//     }
//     Double_t arimean = acellslice[c]->GetMean();
//     asetpar[0] = sliceN;
//     asetpar[1] = arimean;
//     asetpar[2] = aADCtsig;
//     afitl = arimean - 5*aADCtsig;
//     if( col==11 ) afitl=40.;
//     afith = arimean + 4*aADCtsig;
//     TF1 *gausfit = new TF1("gausfit",Gfit,afitl,afith,3);
//     gausfit->SetLineWidth(4);
//     gausfit->SetParameter(0,asetpar[0]);
//     gausfit->SetParameter(1,asetpar[1]);
//     gausfit->SetParLimits(1,afitl,afith);
//     gausfit->SetParameter(2,asetpar[2]);

//     acellslice[c]->Fit("gausfit","RBM");
//     acellslice[c]->Draw();

//     acval[c] = gausfit->GetParameter(1);
//     acerr[c] = gausfit->GetParameter(2);
//     acvalw[c] = gausfit->GetParameter(1);
//     acellslice[c]->SetTitle(Form("Mean:%f Sigma:%f",acval[c],acerr[c]));    

//     acval_avg += acval[c];
//     acval_Ng++;

//   }
//   ADCt_top->Write();
//   ADCt_bot->Write();

//   acval_avg /= acval_Ng;

  
//   //Make graphs with errors for reporting
//   TGraphErrors *gtdc_c = new TGraphErrors( kNcell, tcell, tcvalw, cellerr, tcerr );
//   gtdc_c->GetXaxis()->SetLimits(-10,290);  
//   gtdc_c->GetYaxis()->SetLimits(tllim,tulim);
//   gtdc_c->SetTitle("TDC_{hcal}-TDCmean_{hodo} vs Cell");
//   gtdc_c->GetXaxis()->SetTitle("Cell");
//   gtdc_c->GetYaxis()->SetTitle("TDC_{HCAL}-TDCMEAN_{HODO}");
//   gtdc_c->SetMarkerStyle(20); // idx 20 Circles, idx 21 Boxes
//   gtdc_c->Draw();
//   gtdc_c->Write("gtdc_c");


//   TGraphErrors *gadct_c = new TGraphErrors( kNcell, acell, acvalw, cellerr, acerr );
//   gadct_c->GetXaxis()->SetLimits(-10,290);  
//   gadct_c->GetYaxis()->SetLimits(allim,aulim);
//   gadct_c->SetTitle("ADCt_{hcal}-TDCmean_{hodo} vs Cell");
//   gadct_c->GetXaxis()->SetTitle("Cell");
//   gadct_c->GetYaxis()->SetTitle("ADCt_{HCAL}-TDCMEAN_{HODO}");
//   gadct_c->SetMarkerStyle(20); // idx 20 Circles, idx 21 Boxes
//   gadct_c->Draw();
//   gadct_c->Write("gadct_c");

//   // Create parameter files
//   ofstream tdcP0;
//   ofstream tdcP1;
//   ofstream tdcP2;
//   ofstream adctP0;
//   ofstream adctP1;
//   ofstream adctP2;
//   ofstream tdcoff;
//   ofstream adctoff;

//   ofstream tdcfinaloff; //Add to align to proton TOF/TW/Hodo corrected times

//   //Write to outfiles if !qreplay
//   if( !qreplay ){
//     //Write first TDCvE parameter
//     tdcP0.open( tdctwP0path );
//     tdcP0 << "#HCal tdc vs E fit parameter P0 obtained " << date.c_str() << endl;
//     tdcP0 << "#Exponential fit for i PMTs -> y = P0_i*exp(-P1_i*x)+P2_i, P0_i = " << endl;
//     cout << "#HCal tdc vs E fit parameter P0 obtained " << date.c_str() << endl;
//     cout << "#Exponential fit for i PMTs -> y = P0_i*exp(-P1_i*x)+P2_i, P0_i = " << endl;

//     for( Int_t i=0; i<kNcell; i++ ){   
//       tdcP0 << TDCvseP0[i] << endl;
//       cout << TDCvseP0[i] << endl;
//     }

//     tdcP0.close();
//     cout << endl << endl;

//     //Write second TDCvE parameter
//     tdcP1.open( tdctwP1path );
//     tdcP1 << "#HCal tdc vs E fit parameter P1 obtained " << date.c_str() << endl;
//     tdcP1 << "#Exponential fit for i PMTs -> y = P0_i*exp(-P1_i*x)+P2_i, P1_i = " << endl;
//     cout << "#HCal tdc vs E fit parameter P1 obtained " << date.c_str() << endl;
//     cout << "#Exponential fit for i PMTs -> y = P0_i*exp(-P1_i*x)+P2_i, P1_i = " << endl;

//     for( Int_t i=0; i<kNcell; i++ ){   
//       tdcP1 << TDCvseP1[i] << endl;
//       cout << TDCvseP1[i] << endl;
//     }

//     tdcP1.close();
//     cout << endl << endl;

//     //Write third TDCvE parameter
//     tdcP2.open( tdctwP2path );
//     tdcP2 << "#HCal tdc vs E fit parameter P2 obtained " << date.c_str() << endl;
//     tdcP2 << "#Exponential fit for i PMTs -> y = P0_i*exp(-P1_i*x)+P2_i, P2_i = " << endl;
//     cout << "#HCal tdc vs E fit parameter P2 obtained " << date.c_str() << endl;
//     cout << "#Exponential fit for i PMTs -> y = P0_i*exp(-P1_i*x)+P2_i, P2_i = " << endl;

//     for( Int_t i=0; i<kNcell; i++ ){   
//       tdcP2 << TDCvseP2[i] << endl;
//       cout << TDCvseP2[i] << endl;
//     }

//     tdcP2.close();
//     cout << endl << endl;

//     //Write first ADCtvE parameter
//     adctP0.open( adcttwP0path );
//     adctP0 << "#HCal adct vs E fit parameter P0 obtained " << date.c_str() << endl;
//     adctP0 << "#Exponential fit for i PMTs -> y = P0_i*exp(-P1_i*x)+P2_i, P0_i = " << endl;
//     cout << "#HCal adct vs E fit parameter P0 obtained " << date.c_str() << endl;
//     cout << "#Exponential fit for i PMTs -> y = P0_i*exp(-P1_i*x)+P2_i, P0_i = " << endl;

//     for( Int_t i=0; i<kNcell; i++ ){   
//       adctP0 << ADCtvseP0[i] << endl;
//       cout << ADCtvseP0[i] << endl;
//     }

//     adctP0.close();
//     cout << endl << endl;

//     //Write second ADCtvE parameter
//     adctP1.open( adcttwP1path );
//     adctP1 << "#HCal adct vs E fit parameter P1 obtained " << date.c_str() << endl;
//     adctP1 << "#Exponential fit for i PMTs -> y = P0_i*exp(-P1_i*x)+P2_i, P1_i = " << endl;
//     cout << "#HCal adct vs E fit parameter P1 obtained " << date.c_str() << endl;
//     cout << "#Exponential fit for i PMTs -> y = P0_i*exp(-P1_i*x)+P2_i, P1_i = " << endl;

//     for( Int_t i=0; i<kNcell; i++ ){   
//       adctP1 << ADCtvseP1[i] << endl;
//       cout << ADCtvseP1[i] << endl;
//     }

//     adctP1.close();
//     cout << endl << endl;

//     //Write third ADCtvE parameter
//     adctP2.open( adcttwP2path );
//     adctP2 << "#HCal adct vs E fit parameter P2 obtained " << date.c_str() << endl;
//     adctP2 << "#Exponential fit for i PMTs -> y = P0_i*exp(-P1_i*x)+P2_i, P2_i = " << endl;
//     cout << "#HCal adct vs E fit parameter P2 obtained " << date.c_str() << endl;
//     cout << "#Exponential fit for i PMTs -> y = P0_i*exp(-P1_i*x)+P2_i, P2_i = " << endl;

//     for( Int_t i=0; i<kNcell; i++ ){   
//       adctP2 << ADCtvseP2[i] << endl;
//       cout << ADCtvseP2[i] << endl;
//     }

//     adctP2.close();
//     cout << endl << endl;

//     //Write TDC offsets
//     tdcoff.open( tdcconstPath );
//     tdcoff << "#HCal TDC offsets obtained " << date.c_str() << endl;
//     tdcoff << "#Offsets obtained from fits over TDC distributions = " << endl;
//     cout << "#HCal TDC offsets obtained " << date.c_str() << endl;
//     cout << "#Offsets obtained from fits over TDC distributions = " << endl;
  
//     for( Int_t i=0; i<kNcell; i++ ){
//       if( oldTDCoffsets[i]>0 && ( tcval[i]==0 || abs(tcerr[i])<0.05 ) ){ //This will only work if most channels are closely aligned already
// 	tdcoff << tcval_avg << endl;
// 	cout << tcval_avg << " with adjustment" << endl;
//       }else{
// 	tdcoff << tcval[i]/TDCCalib + oldTDCoffsets[i] - TDC_target/TDCCalib << endl;
// 	cout << tcval[i]/TDCCalib + oldTDCoffsets[i] - TDC_target/TDCCalib << endl;
//       }
//     }  

//     tdcoff.close();

//     cout << endl << endl;

//     //Write ADCt offsets
//     adctoff.open( adctconstPath );
//     adctoff << "#HCal ADCT offsets obtained " << date.c_str() << endl;
//     adctoff << "#Offsets obtained from fits over ADCT distributions = " << endl;
//     cout << "#HCal ADCT offsets obtained " << date.c_str() << endl;
//     cout << "#Offsets obtained from fits over ADCT distributions = " << endl;
  
//     for( Int_t i=0; i<kNcell; i++ ){   
//       adctoff << acval[i] + oldADCtoffsets[i] - ADCt_target << endl;
//       cout << acval[i] + oldADCtoffsets[i] - ADCt_target << endl;
//     }
  
//     adctoff.close();
//   }

//   if( tofreplay && !tofqreplay ){
//     //Write TDC offsets after all corrections are applied. Align to protons (more data!)
//     tdcfinaloff.open( tdcfinalconstPath );
//     tdcfinaloff << "#HCal Fully Corrected (TW/TOF/Hodo) TDC offsets obtained " << date.c_str() << " on iteration " << iter << endl;
//     tdcfinaloff << "#Offsets obtained from fits over TDC distributions = " << endl;
//     cout << "#HCal Fully Corrected (TW/TOF/Hodo) TDC offsets obtained " << date.c_str() << " on iteration " << iter << endl;
//     cout << "#Offsets obtained from fits over TDC distributions = " << endl;
  
//     for( Int_t i=0; i<kNcell; i++ ){
//       if( oldTDCoffsets[i]>0 && ( tcval[i]==0 || abs(tcerr[i])<0.05 ) ){ //This will only work if most channels are closely aligned already - sets channel to average value over all channels if fit and old offset are bad
// 	tdcfinaloff << tcval_avg << endl;
// 	cout << tcval_avg << " with adjustment" << endl;
//       }else{
// 	tdcfinaloff << tcval[i]/TDCCalib + oldTDCoffsets[i] - TDC_target/TDCCalib << endl;
// 	cout << tcval[i]/TDCCalib + oldTDCoffsets[i] - TDC_target/TDCCalib << endl;
//       }
//     }  

//     tdcfinaloff.close();

//   }

//   fout->Write();

//   cout << endl << endl << "Elastic yield for analyzed runs: " << elasYield << ". Total event bins available for calibration = LH2 + LD2: " << TNEV_h << " + " << TNEV_d << " = " << TNEV_h + TNEV_d << ". Total number of events analyzed: " << lh2Events + ld2Events << "." << endl << endl;
//   cout << endl << endl << "Elastic yield for analyzed runs: " << elasYield << ". Total event bins available for calibration = LH2 + LD2: " << TNEV_h << " + " << TNEV_d << " = " << TNEV_h + TNEV_d << ". Total number of events analyzed: " << lh2Events + ld2Events << "." << endl << endl;
  
//   if( tofqreplay ){
//     cout << "Final alignment with all corrections iteration complete. Histograms written to file." << endl;
//     cout << "Final alignment with all corrections iteration complete. Histograms written to file." << endl << endl;
//   }else if( tofreplay ){
//     cout << "All corrections diagnostic iteration complete. Histograms written to file." << endl;
//     cout << "All corrections diagnostic iteration complete. Histograms written to file." << endl << endl;
//   }else if( qreplay ){
//     cout << "Simple diagnostic iteration complete. Histograms written to file." << endl;
//     cout << "Simple diagnostic iteration complete. Histograms written to file." << endl << endl;
//   }else{
//     cout << "Timing analysis complete. Constants and histograms extracted from data and written to file." << endl;
//     cout << "Timing analysis complete. Constants and histograms extracted from data and written to file." << endl << endl;
//   }

//   st->Stop();

//   // Send time efficiency report to console
//   cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;
//   cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

//   report.close();

// }
