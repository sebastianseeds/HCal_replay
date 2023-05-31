//SSeeds 5.8.23 - Script adapted from ecal.C at https://github.com/sebastianseeds/HCal_replay to calibrate energy with elastic events from all available data. Assumes tdc/adct alignment scripts have been written and new alignments are available.
//NOTE: requires $DB_DIR path set correctly in current environment. Script assumes hardware differences on adctoffset timestamps and outputs alignment constants for all timestamps within the configuration provided.
//NOTE: written in anticipation of additional cluster information which includes all block ID per cluster. As of pass0/1, this information is not on the tree.

#include <vector>
#include <iostream>
#include <algorithm>
#include "TCut.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TTreeFormula.h"
#include "TChain.h"
#include "TLatex.h"
#include "TH1D.h"
#include "TStopwatch.h"
#include "../include/hcal.h"

//ADC gain coefficient general constraints
const Int_t first_hcal_chan = 0;
const Int_t total_bins = 240;
const Int_t lower_lim = 0;
const Int_t upper_lim = 140;
const Int_t fit_event_min = 50;
const Double_t observed_adct_sigma = 4.0; //rough estimate
const Double_t ADCt_target = 50.;
const Double_t minEventPerCell = 100;
const Double_t highDelta = 0.01;

// Set up data structure to hold configuration settings by field setting
typedef struct{
  vector<string> tstamp; //timestamp in database where {0,gain;1,tdcoffsets;2,tdccalib;3,adcoffsets}
  vector<string> run; //all runs
  string gcut; //globalcut
  Double_t SBSf;
  Double_t E_e;
  Double_t HCal_d;
  Double_t HCal_th;
  Double_t HCal_sf;
  Double_t HCal_es;
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
  vector<Double_t> oadcgain;
  vector<Double_t> oadctoff;
  vector<Double_t> otdcoff;
} PARA;

// Use nested structure to index parameters by target
typedef struct{
  string targ;
  vector<string> config;
  PARA para[hcal::gMaxfset];
} TARPARA;

// Create structure to define calibration sets for categorical analysis
typedef struct{
  vector<Int_t> target_index;
  vector<Int_t> field_index;
  vector<string> timestamp;
} OFFSET;

typedef struct{
  Double_t GCoeff[hcal::maxHCalChan] = {0.0};
  Double_t GCoeff_oneblock[hcal::maxHCalChan] = {0.0};
  Double_t GCoeff_divide[hcal::maxHCalChan] = {0.0};
  TMatrixD Ma(hcal::maxHCalChan,hcal::maxHCalChan);
  TMatrixD Ma_oneblock(hcal::maxHCalChan,hcal::maxHCalChan);
  TVectorD ba(hcal::maxHCalChan);
  TVectorD ba_oneblock(hcal::maxHCalChan);
  Int_t NEV[hcal::maxHCalChan] = {0};
  Int_t NEV_oneblock[hcal::maxHCalChan] = {0};
  Double_t err[hcal::maxHCalChan] = {0.};
  Double_t err_ev[hcal::maxHCalChan] = {0.};
  Double_t err_oneblock[hcal::maxHCalChan] = {0.};
  Double_t err_ev_oneblock[hcal::maxHCalChan] = {0.};
  Int_t deadclus = 0;
} GCAL;

//Main <experiment> <configuration> <best-cluster-option><quasi-replay-option>; qreplay should only be performed after new offsets obtained
void ecal( const char *exp = "gmn", Int_t config=4, bool bclus = false; bool qreplay = true ){

  // Setup qreplay int index for reporting
  Int_t qreplayIdx = 0;
  if( qreplay )
    qreplayIdx = 1;

  // Setup bclus int index for reporting
  Int_t bclusIdx = 0;
  if( bclus )
    bclusIdx = 1;

  // Get index for timestamp type
  Int_t ecalidx;
  Int_t toffidx;
  Int_t tcalidx;
  Int_t aoffidx;
  for( Int_t s=0; s<hcal::gNstamp; s++ ){
    if( hcal::gStamp[s].compare("gain")==0 )
      ecalidx=s;
    if( hcal::gStamp[s].compare("tdcoff")==0 )
      toffidx=s;
    if( hcal::gStamp[s].compare("tdccalib")==0 )
      tcalidx=s;
    if( hcal::gStamp[s].compare("adctoff")==0 )
      aoffidx=s;
  }

  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );
  
  // Get the date
  string date = hcal::getDate();

  // Set up output analysis file
  TFile *fout = new TFile( Form("outfiles/ecal_%s_conf%d_qr%d_bclus%d.root",exp,config,qreplayIdx,bclusIdx), "RECREATE" );

  // Declare path and open log file for reporting
  ofstream report;
  string logPath = Form("logs/log_ecal_%s_conf%d_qreplay%d_bclus%d.log",exp,config,qreplayIdx,bclusIdx);

  report.open( logPath );
  report << "#Error and performance report from " << exp << " config " << config << " obtained " << date.c_str() << endl << endl;

  // Read in meta-config for multiple field settings within a given SBS config
  string metaconfigname = Form("../config/%s/conf%d/%s_conf%d_meta.cfg",exp,config,exp,config);
  ifstream metaconfigfile(metaconfigname);

  // Create array of structures to hold all parameters by target
  TARPARA allp[hcal::gNtar];
  // Create structure to define calibration sets
  OFFSET oset;

  TString mcline;
  for( Int_t t=0; t<hcal::gNtar; t++ ){

    //targ_config[t].targ = hcal::gTar[t];
    allp[t].targ = hcal::gTar[t];
    string end = "end" + hcal::gTar[t];

    while( mcline.ReadLine( metaconfigfile ) && !mcline.BeginsWith(end) ){
      if( !mcline.BeginsWith("#") ){
	string configpath = Form("../config/%s/conf%d/",exp,config);

	configpath += (string)mcline;

	allp[t].config.push_back(configpath);
	report << "Recording " << hcal::gTar[t] << " config file: " << configpath << ".." << endl;
      }    
    }
  } 

  cout << endl << "Configuration " << config << " " << exp << " energy calibration parameters loaded for " << endl;

  // Loop over potential targets
  for( Int_t t=0; t<hcal::gNtar; t++ ){
    string target = hcal::gTar[t];
    Int_t tarfset = allp[t].config.size();

    report << endl << endl << "Total field settings for target " << target << ": " << tarfset << endl;

    // Loop over all available field settings
    for( Int_t f=0; f<tarfset; f++ ){

      report << endl << endl << "Reading in " << exp << " config " << config << " " << target << " field set " << f << " analysis parameters.." << endl;

      // Record timestamps
      ifstream configfile(allp[t].config[f]);
      TString currentline;
      Int_t stampidx=0;
      while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endstamp") ){
	if( !currentline.BeginsWith("#") ){
	  allp[t].para[f].tstamp.push_back((string)currentline);
	  report << "Recording " << hcal::gStamp[stampidx] << " timestamp: " << currentline << ".." << endl;
	  stampidx++;
	}    
      }

      // Record data files
      while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endlist") ){
	if( !currentline.BeginsWith("#") ){
	  allp[t].para[f].run.push_back((string)currentline);
	  report << "Recording file: " << currentline << ".." << endl;
	}    
      }

      // Record globalcuts
      while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endcut") ){
	if( !currentline.BeginsWith("#") ){
	  allp[t].para[f].gcut += (string)currentline;
	  report << "Recording globalcut: " << currentline << ".." << endl;
	}    
      }
      
      // Record physics cuts and parameters
      while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("#") ){
	TObjArray *tokens = currentline.Tokenize(" ");
	Int_t ntokens = tokens->GetEntries();
	if( ntokens>1 ){
	  TString skey = ( (TObjString*)(*tokens)[0] )->GetString();
	  if( skey == "SBSf" ){
	    TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	    allp[t].para[f].SBSf = sval.Atof();
	    report << "Recording beam energy: " << allp[t].para[f].SBSf << endl;
	  }
	  if( skey == "E_e" ){
	    TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	    allp[t].para[f].E_e = sval.Atof();
	    report << "Loading beam energy: " << allp[t].para[f].E_e << endl;
	  }
	  if( skey == "HCal_d" ){
	    TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	    allp[t].para[f].HCal_d = sval.Atof();
	    report << "Loading HCal distance: " << allp[t].para[f].HCal_d << endl;
	  }
	  if( skey == "HCal_th" ){
	    TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	    allp[t].para[f].HCal_th = sval.Atof() * TMath::DegToRad();	
	    report << "Loading HCal angle: " << allp[t].para[f].HCal_th << endl;
	  }
	  if( skey == "HCal_sf" ){
	    TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	    allp[t].para[f].HCal_sf = sval.Atof();	
	    report << "Loading HCal sampling fraction: " << allp[t].para[f].HCal_sf << endl;
	  }
	  if( skey == "HCal_es" ){
	    TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	    allp[t].para[f].HCal_es = sval.Atof();	
	    report << "Loading HCal E/sig factor: " << allp[t].para[f].HCal_es << endl;
	  }
	  if( skey == "BB_th" ){
	    TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	    allp[t].para[f].BB_th = sval.Atof() * TMath::DegToRad();	
	    report << "Loading BBCal angle: " << allp[t].para[f].BB_th << endl;
	  }
	  if( skey == "W2_mean" ){
	    TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	    allp[t].para[f].W2_mean = sval.Atof();
	    report << "Loading W2 mean cut: " << allp[t].para[f].W2_mean << endl;
	  }
	  if( skey == "W2_sig" ){
	    TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	    allp[t].para[f].W2_sig = sval.Atof();
	    report << "Loading W2 sigma cut: " << allp[t].para[f].W2_sig << endl;
	  }
	  if( skey == "dx0_n" ){
	    TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	    allp[t].para[f].dx0_n = sval.Atof();
	    report << "Loading x position of neutron spot: " << allp[t].para[f].dx0_n << endl;
	  }
	  if( skey == "dx0_p" ){
	    TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	    allp[t].para[f].dx0_p = sval.Atof();
	    report << "Loading y position of proton spot: " << allp[t].para[f].dx0_p << endl;
	  }
	  if( skey == "dy0" ){
	    TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	    allp[t].para[f].dy0 = sval.Atof();
	    report << "Loading y position of both hadron spots: " << allp[t].para[f].dy0 << endl;
	  }
	  if( skey == "dx_sig_n" ){
	    TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	    allp[t].para[f].dx_sig_n = sval.Atof();
	    report << "Loading x sigma of neutron spot: " << allp[t].para[f].dx_sig_n << endl;
	  }
	  if( skey == "dx_sig_p" ){
	    TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	    allp[t].para[f].dx_sig_p = sval.Atof();
	    report << "Loading x sigma of proton spot: " << allp[t].para[f].dx_sig_p << endl;
	  }
	  if( skey == "dy_sig" ){
	    TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	    allp[t].para[f].dy_sig = sval.Atof();
	    report << "Loading y sigma of both hadron spots: " << allp[t].para[f].dy_sig << endl;
	  }
	  if( skey == "atime0" ){
	    TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	    allp[t].para[f].atime0 = sval.Atof();
	    report << "Loading ADC time mean: " << allp[t].para[f].atime0 << endl;
	  }
	  if( skey == "atime_sig" ){
	    TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	    allp[t].para[f].atime_sig = sval.Atof();
	    report << "Loading ADC time sigma: " << allp[t].para[f].atime_sig << endl;
	  }
	  if( skey == "useAlshield" ){
	    TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	    allp[t].para[f].useAlshield = sval.Atoi();
	    report << "Loading Aluminum absorber option: " << allp[t].para[f].useAlshield << endl;
	  }
	}
	delete tokens;
      } 

      //Record old gain and offset parameters by tstamp from database (assumes one file for all timstamps)
      string DBpath = gSystem->Getenv("DB_DIR");
      TString adcGainPath = DBpath + "/db_sbs.hcal.dat";
      ifstream adcGainFile( adcGainPath );
      report << endl << "Recording previous gain coefficients from database file " << adcGainPath << " with timestamp " << allp[t].para[f].tstamp[ecalidx] << ".." << endl;
      if( !adcGainFile ){
	cerr << endl << "ERROR: No input constant file present -> path to db_sbs.hcal.dat expected." << endl;
	return 0;
      }

      Int_t n0=0;
      Double_t d0=0;
      string line;
      bool read_offset = false;
      bool found_tstamp = false;
      if( allp[t].para[f].tstamp[ecalidx].compare("none")==0 )
	found_tstamp = true;
      
      while( getline( adcGainFile, line ) ){
	
	if( n0==hcal::maxHCalChan ) break;
	
	TString Tline = (TString)line;
	
	if( Tline.BeginsWith(allp[t].para[f].tstamp[ecalidx]) && !found_tstamp ){
	  found_tstamp = true;
	  continue;
	}
	
	if( Tline.BeginsWith("sbs.hcal.adc.gain") && found_tstamp && !read_offset ){
	  read_offset = true;
	  continue;
	}

	if( found_tstamp && read_offset ){
	  istringstream iss( line );
	  while( iss >> d0 ){
	    allp[t].para[f].oadcgain.push_back( d0 );
	    n0++;
	  }
	}
      }

      for( Int_t r=0; r<hcal::maxHCalRows; r++){
	for( Int_t c=0; c<hcal::maxHCalCols; c++){
	  Int_t i = r*hcal::maxHCalCols+c;
	  report << allp[t].para[f].oadcgain[i] << " ";
	}
	report << endl;
      }
	
      //Get old adct offsets for replacement with new ones
      n0=0;
      d0=0;
      string line;
      read_offset = false;
      found_tstamp = false;
      if( allp[t].para[f].tstamp[aoffidx].compare("none")==0 )
	found_tstamp = true;
      
      while( getline( adcGainFile, line ) ){
	
	if( n0==hcal::maxHCalChan ) break;
	
	TString Tline = (TString)line;
	
	if( Tline.BeginsWith(allp[t].para[f].tstamp[aoffidx]) && !found_tstamp ){
	  found_tstamp = true;
	  continue;
	}
	
	if( Tline.BeginsWith("sbs.hcal.adc.timeoffset") && found_tstamp && !read_offset ){
	  read_offset = true;
	  continue;
	}

	if( found_tstamp && read_offset ){
	  istringstream iss( line );
	  while( iss >> d0 ){
	    allp[t].para[f].oadctoff.push_back( d0 );
	    n0++;
	  }
	}
      }
	
      report << endl << endl << "Old ADCt offsets: " << endl;
      for( Int_t r=0; r<hcal::maxHCalRows; r++){
	for( Int_t c=0; c<hcal::maxHCalCols; c++){
	  Int_t i = r*hcal::maxHCalCols+c;
	  report << allp[t].para[f].oadctoff[i] << " ";
	}
	report << endl;
      }

      //Get old tdc offsets for replacement with new ones
      n0=0;
      d0=0;
      string line;
      read_offset = false;
      found_tstamp = false;
      if( allp[t].para[f].tstamp[toffidx].compare("none")==0 )
	found_tstamp = true;
      
      while( getline( adcGainFile, line ) ){
	
	if( n0==hcal::maxHCalChan ) break;
	
	TString Tline = (TString)line;
	
	if( Tline.BeginsWith(allp[t].para[f].tstamp[toffidx]) && !found_tstamp ){
	  found_tstamp = true;
	  continue;
	}
	
	if( Tline.BeginsWith("sbs.hcal.tdc.offset") && found_tstamp && !read_offset ){
	  read_offset = true;
	  continue;
	}

	if( found_tstamp && read_offset ){
	  istringstream iss( line );
	  while( iss >> d0 ){
	    allp[t].para[f].otdcoff.push_back( d0 );
	    n0++;
	  }
	}
      }

      report << endl << endl << "Old TDC offsets: " << endl;
      for( Int_t r=0; r<hcal::maxHCalRows; r++){
	for( Int_t c=0; c<hcal::maxHCalCols; c++){
	  Int_t i = r*hcal::maxHCalCols+c;
	  report << allp[t].para[f].otdcoff[i] << " ";
	}
	report << endl;
      }

      //Currently set to produce gain coefficients by adc gain timestamp only
      if( std::find(oset.timestamp.begin(),oset.timestamp.end(),allp[t].para[f].tstamp[ecalidx]) == oset.timestamp.end() ){
	oset.timestamp.push_back( allp[t].para[f].tstamp[ecalidx] );
	oset.target_index.push_back( t );
	oset.field_index.push_back( f );
	report << "FOR t " << t << " f " << f << " ADDING timestamp " << oset.timestamp.back();
      }

      cout << "  target " << hcal::gTar[t] << ", field set " << allp[t].para[f].SBSf << "%.. " << endl;
      usleep( hcal::second );

    }//endloop over field settings
  }//endloop over targets
  
  Int_t calsetsize = oset.timestamp.size();

  // Declare paths for new gain coefficients
  string writeCoeffPath = Form("parameters/adcgaincoeff_%s_conf%d_bclus%d.txt",exp,config,bclusIdx);
  string readCoeffPath = Form("parameters/adcgaincoeff_%s_conf%d_bclus%d.txt",exp,config,bclusIdx);

  // Declare path to new adct alignments for better elastic selection
  string readadctParPath = Form("../timing/parameters/adctoffsets_%s_conf%d.txt",exp,config);

  // Declare path to new tdc alignments
  string readtdcParPath = Form("../timing/parameters/tdcoffsets_%s_conf%d.txt",exp,config);  

  if( qreplay )
    writeCoeffPath = "/dev/null"; //Safety to prevent overwriting constants on quasi-replay

  Double_t ntdcoff[hcal::maxHCalChan];
  Double_t nadctoff[hcal::maxHCalChan];
  Double_t nadcgaincoeff[hcal::maxHCalChan];

  ifstream readadctParFile( readadctParPath );
  if( !readadctParFile ){
    cerr << endl << "ERROR: No input constant file present -> parameters/<adctoffsetfile> expected." << endl;
    return 0;
  }

  report << "Loading newly calculated adct offset parameters.." << endl;

  // Record new adct offsets
  Int_t n0=0;
  Double_t d0=0;
  string adctline;
  bool read_offset = false;
    
  while( getline( readadctParFile, adctline ) ){
      
    if( n0==hcal::maxHCalChan ) break;
      
    TString Tline = (TString)adctline;
 
    if( Tline.BeginsWith("sbs.hcal.adc.timeoffset") && !read_offset ){
      read_offset = true;
      continue;
    }
      
    if( read_offset ){
      istringstream iss( adctline );
      while( iss >> d0 ){
	nadctoff[n0]=d0;
	n0++;
      }
    }
  }
  
  report << endl << endl << "New ADCT offsets from file: " << endl;
  cout << endl << endl << "New ADCT offsets from file: " << endl;
    
  for( Int_t r=0; r<hcal::maxHCalRows; r++){
    for( Int_t c=0; c<hcal::maxHCalCols; c++){
      Int_t i = r*hcal::maxHCalCols+c;
      report << nadctoff[i] << " ";
      cout << nadctoff[i] << " ";
    }
    report << endl;
    cout << endl;
  }

  cout << "New ADCT offsets loaded for quasi-replay.." << endl;

  ifstream readtdcParFile( readtdcParPath );
  if( !readtdcParFile ){
    cerr << endl << "ERROR: No input constant file present -> parameters/<adctoffsetfile> expected." << endl;
    return 0;
  }

  report << "Loading newly calculated adct offset parameters.." << endl;

  // Record new tdc offsets
  n0=0;
  d0=0;
  string tdcline;
  read_offset = false;
    
  while( getline( readtdcParFile, tdcline ) ){
      
    if( n0==hcal::maxHCalChan ) break;
      
    TString Tline = (TString)tdcline;
 
    if( Tline.BeginsWith("sbs.hcal.tdc.offset") && !read_offset ){
      read_offset = true;
      continue;
    }
      
    if( read_offset ){
      istringstream iss( tdcline );
      while( iss >> d0 ){
	ntdcoff[n0]=d0;
	n0++;
      }
    }
  }
  
  report << endl << endl << "New TDC offsets from file: " << endl;
  cout << endl << endl << "New TDC offsets from file: " << endl;
    
  for( Int_t r=0; r<hcal::maxHCalRows; r++){
    for( Int_t c=0; c<hcal::maxHCalCols; c++){
      Int_t i = r*hcal::maxHCalCols+c;
      report << ntdcoff[i] << " ";
      cout << ntdcoff[i] << " ";
    }
    report << endl;
    cout << endl;
  }

  cout << "New TDC offsets loaded for quasi-replay.." << endl;

  //  if quasi-replay selected
  if( qreplay ){
    ifstream readCoeffFile( readCoeffPath );
    if( !readCoeffFile ){
      cerr << endl << "ERROR: No input constant file present -> parameters/<adcgaincoefffile> expected." << endl;
      return 0;
    }

    report << "Loading newly calculated adc gain coefficients.." << endl;

    // Record new adc gain offsets
    n0=0;
    d0=0;
    string cline;
    read_offset = false;
    
    while( getline( readCoeffFile, cline ) ){
      
      if( n0==hcal::maxHCalChan ) break;
      
      TString Tline = (TString)cline;
 
      if( Tline.BeginsWith("sbs.hcal.adc.gain") && !read_offset ){
	  read_offset = true;
	  continue;
      }
      
      if( read_offset ){
	istringstream iss( cline );
	while( iss >> d0 ){
	  nadcgaincoeff[n0]=d0;
	  n0++;
	}
      }
    }
    
    report << endl << endl << "New ADC gain coefficients from file: " << endl;
    cout << endl << endl << "New ADC gain coefficients from file: " << endl;
    
    for( Int_t r=0; r<hcal::maxHCalRows; r++){
      for( Int_t c=0; c<hcal::maxHCalCols; c++){
	Int_t i = r*hcal::maxHCalCols+c;
	report << nadcgaincoeff[i] << " ";
	cout << nadcgaincoeff[i] << " ";
      }
      report << endl;
      cout << endl;
    }

    cout << "New ADC gain coefficients loaded for quasi-replay.." << endl;

  }//endif qreplay

  // Re-allocate memory at each run to load different cuts/parameters
  TChain *C = nullptr;

  // Create simple output tree
  TTree *P = new TTree("P","Analysis Tree");

  // Timing
  Double_t pblkid_out; //hcal cluster, primary block id
  Double_t tdc_out; //hcal cluster tdc time, tree
  Double_t atime_out; //hcal cluster adc time, tree
  Double_t hodotmean_out; //hodoscope cluster mean tdc time

  // Physics
  Double_t dx_out; //hcal cluster dx
  Double_t dy_out; //hcal cluster dy
  Double_t W2_out; //Invariant mass squared W2
  Double_t Q2_out; //Inverse momentum squared Q2
  Double_t hcale_out; //hcal cluster energy
  Double_t sfrac_out; //hcal cluster energy
  Double_t thetapq_pout; //proton thetapq
  Double_t thetapq_nout; //neutron thetapq
  Double_t nu_out; //hadron KE
  Double_t ep_out; //track reconstructed e' momentum
  Int_t nblk_out; //total blocks in primary cluster
  Int_t pid_out; //particle id (0:ambiguous,1:proton,2:neutron)
  Int_t mag_out; //sbs magnetic field strength (percent)
  Int_t run_out; //run number
  Int_t seg_out; //segment number
  Int_t tar_out; //target, LH2 or LD2
  Int_t badclus; //0: all blocks pass coin, 1: at least one block fails coin

  //cluster block variables
  Double_t cblkid_out[hcal::maxHCalBlk];
  Double_t cblkatime_out[hcal::maxHCalBlk];
  Double_t cblktime_out[hcal::maxHCalBlk];
  Double_t cblke_out[hcal::maxHCalBlk];

  // Set output tree branches
  P->Branch( "pblkid", &pblkid_out, "pblkid/D" );
  P->Branch( "tdc", &tdc_out, "tdc/D" );
  P->Branch( "atime", &atime_out, "atime/D" );
  P->Branch( "hodotmean", &hodotmean_out, "hodotmean/D" );
  P->Branch( "dx", &dx_out, "dx/D" );
  P->Branch( "dx", &dx_out, "dx/D" );
  P->Branch( "dy", &dy_out, "dy/D" );
  P->Branch( "W2", &W2_out, "W2/D" );
  P->Branch( "Q2", &Q2_out, "Q2/D" );
  P->Branch( "nu", &nu_out, "nu/D" );
  P->Branch( "hcale", &hcale_out, "hcale/D" );
  P->Branch( "sfrac", &sfrac_out, "sfrac/D" );
  P->Branch( "thetapq_p", &thetapq_pout, "thetapq_p/D" );
  P->Branch( "thetapq_n", &thetapq_nout, "thetapq_n/D" );
  P->Branch( "nblk", &nblk_out, "nblk_out/I" );
  P->Branch( "pid", &pid_out, "pid_out/I" );
  P->Branch( "mag", &mag_out, "mag_out/I" );
  P->Branch( "run", &run_out, "run_out/I" );
  P->Branch( "seg", &seg_out, "seg_out/I" );
  P->Branch( "tar", &tar_out, "tar_out/I" );
  P->Branch( "badclus", &badclus_out, "badclus_out/I" );

  P->Branch( "cblkid", &cblkid_out, Form("cblkid[%d]/D",hcal::maxHCalBlk) );
  P->Branch( "cblkatime", &cblkatime_out, Form("cblkatime[%d]/D",hcal::maxHCalBlk) );
  P->Branch( "cblktime", &cblktime_out, Form("cblktime[%d]/D",hcal::maxHCalBlk) );
  P->Branch( "cblke", &cblke_out, Form("cblke[%d]/D",hcal::maxHCalBlk) );

  cout << "Analysis output tree \"P\" created.." << endl;
  report << "Analysis output tree \"P\" created.." << endl;

  // //Set up matrices, vectors, and arrays to minimize chisquare over many elastic events
  // Double_t GCoeff[hcal::maxHCalChan] = {0.0};
  // Double_t GCoeff_oneblock[hcal::maxHCalChan] = {0.0};
  // Double_t GCoeff_divide[hcal::maxHCalChan] = {0.0};
  // TMatrixD Ma(hcal::maxHCalChan,hcal::maxHCalChan);
  // TMatrixD Ma_oneblock(hcal::maxHCalChan,hcal::maxHCalChan);
  // TMatrixD Ma_err(hcal::maxHCalChan,hcal::maxHCalChan);
  // TVectorD ba(hcal::maxHCalChan);
  // TVectorD ba_oneblock(hcal::maxHCalChan);
  // TVectorD ba_err(hcal::maxHCalChan);
  // Int_t NEV[hcal::maxHCalChan] = {0};
  // Int_t NEV_oneblock[hcal::maxHCalChan] = {0};
  // Double_t err[hcal::maxHCalChan] = {0.};
  // Double_t err_ev[hcal::maxHCalChan] = {0.};
  // Double_t err_oneblock[hcal::maxHCalChan] = {0.};
  // Double_t err_ev_oneblock[hcal::maxHCalChan] = {0.};
  // Int_t deadclus = 0;
  
  GCAL gset[calsetsize];

  //////////////////////////
  // Begin looping over data
  
  // Loop over potential targets
  for( Int_t t=0; t<hcal::gNtar; t++ ){
    
    string target = hcal::gTar[t];
    //Int_t tarfset = targ_config[t].config.size();
    Int_t tarfset = allp[t].config.size();
    
    // Loop over all available field settings
    for( Int_t f=0; f<tarfset; f++ ){
      
      Double_t mag = allp[t].para[f].SBSf;
      
      cout << "Looping over data in " << exp << " config " << config << ", target " << target << ", field " << mag << ".."  << endl;
      report << "Looping over data in " << exp << " config " << config << ", target " << target << ", field " << mag << ".."  << endl;
      
      Int_t nruns = allp[t].para[f].run.size();
      
      for( Int_t r=0; r<nruns; r++ ){
	
	//Add the run(s) from config structure
	C = new TChain("T");
	C->Add(allp[t].para[f].run[r].c_str());
	
	//Get run from string
	string runpathstr = allp[t].para[f].run[r];
	string runstr = "";
	string segstr = "";
	Int_t skipN = 0;
	Int_t rnct = 5; //run number digits
	
	for( Int_t l=0; l<runpathstr.length(); l++ ){
	  if( runpathstr[l]=='_' && skipN==1 ){ //read 5 digits from the second underscore
	    for( Int_t d=1; d<=rnct; d++ ){
	      runstr += runpathstr[l+d];
	    }
	  }
	  if( runpathstr[l]=='_' && skipN==4 ){ //read all remaining digits from last underscore
	    Int_t b = l;
	    while( runpathstr[b+1]!='.' ){
	      segstr += runpathstr[b+1];
	      b++;
	    }
	  }
	  if( runpathstr[l]=='_' ){ //skip first underscore
	    skipN++;
	  }
	}
	
	Int_t run = stoi( runstr );
	Int_t seg = stoi( segstr );
	
	report << "Analyzing run " << run << " segment " << seg << ".." << endl;
	
	// Establish set index to build calibration sets with TARPARA object
	Int_t setidx=-1;
	
	for( Int_t set=0; set<tssize; set++ ){
	  report << "Set number " << set << " oset: " << oset.timestamp[set] << ", allp: " << allp[t].para[f].tstamp[offidx] << endl;
	  if( oset.timestamp[set].compare(allp[t].para[f].tstamp[offidx]) == 0 )
	    setidx = set;
	}
	if( setidx==-1 ){
	  cout << "ERROR: run DB timestamp not in calibration set on set index " << setidx << ", with parameter timestamp: " << allp[t].para[f].tstamp[offidx] << ".." << endl;
	  return;
	}
	
	// Setting up ROOT tree branch addresses
	C->SetBranchStatus("*",0);    
	Double_t BBtr_p[hcal::maxTracks], BBtr_px[hcal::maxTracks], BBtr_py[hcal::maxTracks], BBtr_pz[hcal::maxTracks];
	Double_t BBtr_vz[hcal::maxTracks];
	Double_t BBtr_n, BBps_x, BBps_y, BBps_e, BBsh_x, BBsh_y, BBsh_e;	
	Double_t HCALx, HCALy, HCALe;
	Double_t pblkrow, pblkcol, nblk, nclus;
	Int_t Ncid;
	Double_t cblkid[hcal::maxHCalChan], cblke[hcal::maxHCalChan], cblkatime[hcal::maxHCalChan], cblktime[hcal::maxHCalChan], cblkagain[hcal::maxHCalChan];
	Double_t cid[hcal::maxHCalClus], crow[hcal::maxHCalClus], ccol[hcal::maxHCalClus], ce[hcal::maxHCalClus], cx[hcal::maxHCalClus], cy[hcal::maxHCalClus], catime[hcal::maxHCalClus];
	Double_t HODOtmean;
	// Speed up processing by switching on only useful branches
	C->SetBranchStatus( "*", 0 );
	C->SetBranchStatus( "sbs.hcal.x", 1 );
	C->SetBranchStatus( "sbs.hcal.y", 1 );
	C->SetBranchStatus( "sbs.hcal.e", 1 );
	C->SetBranchStatus( "sbs.hcal.rowblk", 1 );
	C->SetBranchStatus( "sbs.hcal.colblk", 1 );
	C->SetBranchStatus( "sbs.hcal.nblk", 1 );
	C->SetBranchStatus( "sbs.hcal.nclus", 1 );
	C->SetBranchStatus( "sbs.hcal.clus_blk.id", 1 );
	C->SetBranchStatus( "sbs.hcal.clus_blk.e", 1 );
	C->SetBranchStatus( "sbs.hcal.clus_blk.tdctime", 1 );
	C->SetBranchStatus( "sbs.hcal.clus_blk.atime", 1 );
	//C->SetBranchStatus( "sbs.hcal.clus_blk.again", 1 );
	C->SetBranchStatus( "sbs.hcal.clus.id", 1 );
	C->SetBranchStatus( "sbs.hcal.clus.row", 1 );
	C->SetBranchStatus( "sbs.hcal.clus.col", 1 );
	C->SetBranchStatus( "sbs.hcal.clus.e", 1 );
	C->SetBranchStatus( "sbs.hcal.clus.x", 1 );
	C->SetBranchStatus( "sbs.hcal.clus.y", 1 );
	C->SetBranchStatus( "sbs.hcal.clus.atime", 1 );
	C->SetBranchStatus( "bb.tr.n", 1 );
	C->SetBranchStatus( "bb.tr.px", 1 );
	C->SetBranchStatus( "bb.tr.py", 1 );
	C->SetBranchStatus( "bb.tr.pz", 1 );    
	C->SetBranchStatus( "bb.tr.vz", 1 );
	C->SetBranchStatus( "bb.tr.p", 1 );
	C->SetBranchStatus( "bb.ps.e", 1 );
	C->SetBranchStatus( "bb.ps.x", 1 );
	C->SetBranchStatus( "bb.ps.y", 1 );
	C->SetBranchStatus( "bb.sh.e", 1 );
	C->SetBranchStatus( "bb.sh.x", 1 );
	C->SetBranchStatus( "bb.sh.y", 1 );
	C->SetBranchStatus( "bb.hodotdc.clus.tmean", 1 );
	C->SetBranchStatus( "bb.hodotdc.clus.tmean", 1 );
	C->SetBranchStatus( "bb.gem.track.nhits", 1 );
	C->SetBranchStatus( "bb.etot_over_p", 1 );
	C->SetBranchStatus( "Ndata.sbs.hcal.clus.id", 1 ); //Odd maxing out at 10 clusters on all cluster Ndata objects, so this is needed in addition to sbs.hcal.nclus

	// Linking memory
	C->SetBranchAddress( "sbs.hcal.x", &HCALx );
	C->SetBranchAddress( "sbs.hcal.y", &HCALy );
	C->SetBranchAddress( "sbs.hcal.e", &HCALe );
	C->SetBranchAddress( "sbs.hcal.rowblk", &pblkrow );
	C->SetBranchAddress( "sbs.hcal.colblk", &pblkcol );
	C->SetBranchAddress( "sbs.hcal.nblk", &nblk ); // Total number of blocks in highest E clus
	C->SetBranchAddress( "sbs.hcal.nclus", &nclus ); // Total number of clusters
	C->SetBranchAddress( "sbs.hcal.clus_blk.id", cblkid ); // kNcell-1 index for each block
	C->SetBranchAddress( "sbs.hcal.clus_blk.e", cblke ); // Array of block energies
	C->SetBranchAddress( "sbs.hcal.clus_blk.tdctime", cblktime ); // Array of block TDC times
	C->SetBranchAddress( "sbs.hcal.clus_blk.atime", cblkatime ); // Array of block ADC times
	//C->SetBranchAddress( "sbs.hcal.clus_blk.again", cblkagain ); // Array of block ADC gain coeff
	C->SetBranchAddress( "sbs.hcal.clus.id", cid );
	C->SetBranchAddress( "sbs.hcal.clus.row", crow );
	C->SetBranchAddress( "sbs.hcal.clus.col", ccol );
	C->SetBranchAddress( "sbs.hcal.clus.e", ce );
	C->SetBranchAddress( "sbs.hcal.clus.x", cx );
	C->SetBranchAddress( "sbs.hcal.clus.y", cy );
	C->SetBranchAddress( "sbs.hcal.clus.atime", catime );
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
	C->SetBranchAddress( "bb.hodotdc.clus.tmean", &HODOtmean );
	C->SetBranchAddress( "Ndata.sbs.hcal.clus.id", &Ncid ); //Odd maxing out at 10 clusters on all cluster Ndata objects, so this is needed in addition to sbs.hcal.nclus

	TCut GCut = allp[t].para[f].gcut.c_str();
	TTreeFormula *GlobalCut = new TTreeFormula( "GlobalCut", GCut, C );

	// Set up hcal coordinate system with hcal angle wrt exit beamline
	vector<TVector3> hcalaxes; hcal::sethcalaxes( allp[t].para[f].HCal_th, hcalaxes );
	TVector3 hcalorigin = allp[t].para[f].HCal_d*hcalaxes[2] + hcal::HCalvoff*hcalaxes[0];
	Double_t BdL = hcal::maxSBSfield * hcal::sbsdipolegap * (mag/100); //scale crudely by magnetic field
	Double_t Eloss_outgoing = hcal::celldiameter/2.0/sin(allp[t].para[f].BB_th) * hcal::tarrho[t] * hcal::dEdx[t];

	long nevent = 0, nevents = C->GetEntries(); 
	Int_t treenum = 0, currenttreenum = 0;

	while (C->GetEntry(nevent++)) {
	  
	  cout << "Analyzing run " << run << ", segment " << seg << ": " << nevent << "/" << nevents << " \r";
	  cout.flush();
	  
	  ///////
	  //Single-loop elastic globalcut method. Save pass/fail for output tree.
	  currenttreenum = C->GetTreeNumber();
	  if( nevent == 1 || currenttreenum != treenum ){
	    treenum = currenttreenum; 
	    GlobalCut->UpdateFormulaLeaves();
	    report << "Updating formula leaves and switching segment at event: " << nevent << endl;
	  }
	  bool failedglobal = GlobalCut->EvalInstance(0) == 0;
	  
	  if( failedglobal ) continue;
	  
	  ///////
	  //HCal Active Area Cut
	  bool failedactivearea = 
	    pblkrow==0 || 
	    pblkrow==23 || 
	    pblkcol==0 || 
	    pblkcol==11;

	  if( !bclus && failedactivearea ) continue; //All events with primary cluster element on edge blocks cut

	  ///////
	  //HCal coincidence time cut (using adctime while hcal tdc suspect, new offsets)
	  Int_t pblkid = cblkid[0]; //define primary block, primary cluster ID
	  Double_t natime = cblkatime[0]+allp[t].para[f].oadctoff[pblkid]-nadctoff[pblkid]; //new atime
	  Double_t ntime = cblktime[0]+allp[t].para[f].otdcoff[pblkid]-ntdcoff[pblkid]; //new tdc
	  Double_t atime0 = allp[t].para[f].atime0; //observed elastic peak in adc time
	  Double_t atimesig = allp[t].para[f].atime_sig; //observed width of elastic peak in adc time

	  bool failedcoin = abs(natime-atime0)>3*atimesig;
	  
	  if( failedcoin ) continue; //All events where adctime outside of reasonable window cut
	
	  //ADC arrays reset per event 
	  Double_t A[hcal::maxHCalChan] = {0.0}; // Array to keep track of ADC values per cell. Outscope on each ev
	  Double_t A_oneblock[hcal::maxHCalChan] = {0.0}; // Array to keep track of ADC values per cell for one block clusters only. Outscope on each ev


	  ///////
	  //Physics calculations
	  //correct beam energy with vertex information, primary track
	  Double_t ebeam_c = allp[t].para[f].E_e - ( (BBtr_vz[0]+hcal::l_tgt[t]/2.0) * hcal::tarrho[t] * hcal::dEdx[t] + hcal::uwallthick[t] * hcal::crho[t] * hcal::cdEdx[t] );

	  TVector3 vertex( 0., 0., BBtr_vz[0] );

	  //reconstructed momentum, corrected for mean energy loss exiting the target (later we'll also want to include Al shielding on scattering chamber)
	  Double_t precon = BBtr_p[0] + Eloss_outgoing;

	  //set up four-momenta with some empty for various calculation methods
	  TLorentzVector pbeam( 0., 0., ebeam_c, ebeam_c ); //beam momentum
	  //TLorentzVector pe( px[0], py[0], pz[0], p[0] ); //e' plvect
	  TLorentzVector pe( precon*BBtr_px[0]/BBtr_p[0], precon*BBtr_py[0]/BBtr_p[0], precon*BBtr_pz[0]/BBtr_p[0], precon ); //e' recon plvect
	  TLorentzVector ptarg; //target momentum
	  ptarg.SetPxPyPzE( 0., 0., 0., hcal::M_t[t] );
	  TLorentzVector q = pbeam - pe; //virtual photon mom. or mom. transferred to scatter nucleon (N')
	  TVector3 qv = q.Vect();
	  TLorentzVector pN; //N' momentum
      
	  //simple calculations for e' and N'
	  Double_t etheta = acos( pe.Pz() / pe.E() );
	  Double_t ephi = atan2( pe.Py(), pe.Px() );
	  Double_t pcent = ebeam_c/( 1. + ( ebeam_c/hcal::M_t[t] )*( 1.0 - cos(etheta) ) ); //e' p reconstructed by angles
	  Double_t phNexp = ephi + hcal::PI;
	  Double_t Q2, W2;

	  //e' p reconstruction with track angles (not momentum)
	  Double_t nu = pbeam.E() - pcent;
	  Double_t pNexp = sqrt( pow(nu, 2.) + 2. * hcal::M_t[t] * nu );
	  Double_t thNexp = acos( (pbeam.E()-pcent*cos(etheta)) / pNexp );
	  TVector3 pNhat( sin(thNexp) * cos(phNexp), sin(thNexp) * sin(phNexp), cos(thNexp) );
	  pN.SetPxPyPzE( pNexp*pNhat.X(), pNexp*pNhat.Y(), pNexp*pNhat.Z(), nu+ptarg.E() );
	  Q2 = 2.0 * pbeam.E() * pe.E() * ( 1.0 - cos(etheta) );
	  W2 = pow( hcal::M_t[t], 2.0 ) + 2.0*hcal::M_t[t] * (ebeam_c-pe.E()) - Q2;

	  //Calculate h-arm quantities
	  vector<Double_t> xyhcalexp; hcal::getxyhcalexpect( vertex, pNhat, hcalorigin, hcalaxes, xyhcalexp );
	  TVector3 hcalpos = hcalorigin + HCALx*hcalaxes[0] + HCALy*hcalaxes[1]; //from primary blk
	  Double_t KE_p = nu; //For elastics total predicted by earm
	  Double_t SFrac = HCALe/KE_p; //Measured
	  Double_t dx = HCALx - xyhcalexp[0];
	  Double_t dy = HCALy - xyhcalexp[1];
	  TVector3 neutdir = (hcalpos - vertex).Unit();
	  Double_t protdeflect = tan( 0.3 * BdL / qv.Mag() ) * (allp[t].para[f].HCal_d - (hcal::sbsdist + hcal::sbsdipolegap/2.0) );
	  TVector3 protdir = ( hcalpos + protdeflect * hcalaxes[0] - vertex ).Unit();
	  Double_t thetapq_p = acos( protdir.Dot( qv.Unit() ) );
	  Double_t thetapq_n = acos( neutdir.Dot( qv.Unit() ) );
	  
	  ///////
	  //Target Energy in HCal for calibrations
	  Double_t HCalSF = allp[t].para[f].HCal_sf;
	  Double_t HCalES = allp[t].para[f].HCal_es; 
	  Double_t KE_exp = KE_p*HCalSF/HCalES; //Expected E in HCal

	  ///////
	  //BigBite/SBS Acceptance matching.
	  bool failedaccmatch = 
	    xyhcalexp[1] > hcal::posHCalYf ||
	    xyhcalexp[1] < hcal::posHCalYi ||
	    xyhcalexp[0] > hcal::posHCalXf ||
	    xyhcalexp[0] < hcal::posHCalXi;
	  
	  if( failedaccmatch ) continue;

	  ///////
	  //dy elastic cut
	  Double_t dy0 = allp[t].para[f].dy0;
	  Double_t dysig = allp[t].para[f].dy_sig;
	  bool faileddy = abs(dy-dy0)<3*dysig;

	  if( faileddy ) continue;

	  ///////
	  //PID (evaluated where possible)
	  Double_t pid = -1;
	  Double_t dx0p = allp[t].para[f].dx0_p;
	  Double_t dx0n = allp[t].para[f].dx0_n;
	  Double_t dxsigp = allp[t].para[f].dx_sig_p;
	  Double_t dxsign = allp[t].para[f].dx_sig_n;

	  bool isproton = abs(dx-dx0p)<3*dxsigp;
	  bool isneutron = abs(dx-dx0n)<3*dxsign;
	  bool isambiguous = isproton && isneutron;

	  if( isproton ) 
	    pid = 1;
	  if( isneutron )
	    pid = 2;
	  if( isambiguous )
	    pid = 0;

	  ///////
	  //dx elastic cut
	  if( !isproton && !isneutron ) continue;

	  ///////
	  //Get corrected primary cluster energy
	  Double_t clusE = 0.0;
	  Double_t clusblkE = 0.0;
	  Int_t badclus = 0;
	  for( Int_t blk = 0; blk<Ncid; blk++ ){
	    Int_t blkid = int(cblkid[blk])-1; //-1 necessary since sbs.hcal.clus_blk.id ranges from 1 to 288 (different than just about everything else)
	    Double_t blke = cblke[blk];

	    ////////
	    //Cluster block coincidence adc time check, 5sig estimate
	    if( abs(cblkatime[blk]-atime0)>5*atimesig )
	      badclus = 1;

	    //Correct the block energy with new gain coeff on quasi-replay
	    if( qreplay )
	      blke = cblke[blk]/allp[t].para[f].oadcgain[blkid]*nadcgaincoeff[blkid];

	    clusE += blke;
	    if( blk==0 ) clusblkE += blke;
	    
	    ////////
	    //Fill output tree HCal cluster block variables
	    cblkid_out[blk] = blkid;
	    cblkatime_out[blk] = cblkatime[blk]+allp[t].para[f].oadctoff[blkid]-nadctoff[blkid];
	    cblktime_out[blk] = cblktime[blk]+allp[t].para[f].otdcoff[blkid]-ntdcoff[blkid];
	    cblke_out[blk] = blke;
	  }
	  
	  //Check for inconsistencies
	  cout << clusE << " " << HCALe << endl;
	  cout << clusblkE << " " << cblke[0] << endl;
	  
	  ///////
	  //Fill output tree general variables
	  pblkid_out = cblkid[0];
	  tdc_out = ntime;
	  atime_out = natime;
	  hcale_out = HCALe;
	  dx_out = dx;
	  dy_out = dy;
	  W2_out = W2;
	  Q2_out = Q2;
	  nu_out = nu;
	  hcale_out = HCALe;
	  sfrac_out = SFrac;
	  thetapq_pout = thetapq_p;
	  thetapq_nout = thetapq_n;
	  nblk_out = Ncid;
	  pid_out = pid;
	  mag_out = mag;
	  run_out = run;
	  seg_out = seg;
	  tar_out = t;
	  badclus_out = badclus;
	  hodotmean_out = HODOtmean;

	  P->Fill();

	  ////////
	  //Proceed only if calibrating
	  if( qreplay ) continue;


	  bool badblock = false;
	  clusE = 0.0; //reset cluster energy
	  Double_t cluspC = 0.0;
	  for( Int_t blk = 0; blk<Ncid; blk++ ){
	    badblock = false;
	    Int_t blkid = int(cblkid[blk])-1;
	    Double_t blke = cblke[blk];
	    Double_t blkpC = cblke[blk]/allp[t].para[f].oadcgain[blkid];

	    ////////
	    //Cluster block coincidence time cut (using adct for now)
	    if( abs(cblkatime[blk]-atime0)>5*atimesig )
	      badblock=true;
	    
	    //if( badblock ) continue;

	    clusE += blke;
	    cluspC += blkpC;
	    A[blkid] += blkpC;
	    if( Ncid==1 ){
	      gset[setidx].A_oneblock[blkid] += blkpC;
	      gset[setidx].err_ev_oneblock[blkid] += pow( ((1.+0.5*blkpC)/blkpC+(0.05*KE_p)/KE_p) , 2 ); //Probably better approximation of the error since all the KE_p should be deposited in this block.
	      gset[setidx].NEV_oneblock[blkid]++;
	    }
	    gset[setidx].err_ev[blkid] += pow( ((1.+blkpC)/blkpC+(0.05*KE_p)/KE_p) , 2 ); //Getting ready for std. dev of mean
	    gset[setidx].NEV[blkid]++;
	    
	  }

	  //Build the matrix as simply as possible
	  for(Int_t icol = 0; icol<hcal::maxHCalChan; icol++){
	    gset[setidx].ba(icol) += A[icol];

	    if( Ncid==1 ) 
	      gset[setidx].ba_oneblock(icol) += A_oneblock[icol];

	    for(Int_t irow = 0; irow<hcal::maxHCalChan; irow++){
	      gset[setidx].Ma(icol,irow) += A[icol]*A[irow]/KE_exp;
	      if( Ncid==1 ){	    
		gset[setidx].Ma_oneblock(icol,irow) += A_oneblock[icol]*A_oneblock[irow]/KE_exp;
	      } //endloop oneblock
	    } //endloop over matrix element rows
	  } //endloop over matrix element cols

	}//endloop over event

	// getting ready for the next run
	C->Reset();

      }//endloop over runs
      
    }//endloop over fields
  
  }//endloop over targets
  
  if( !qreplay ){
    cout << endl << "Checking data, inverting matrix, and solving for coefficients.." << endl << endl;
    report << endl << "Checking data, inverting matrix, and solving for coefficients.." << endl << endl;
  }

  //Loop over all independent data sets
  for( Int_t s=0; s<calsetsize; s++ ){
  
    //Reject the bad cells and normalize the oneblock check
    Int_t badcell[hcal::maxHCalChan];
    Int_t badcell_oneblock[hcal::maxHCalChan];
    Int_t cellBad = 0;
    Double_t y[hcal::maxHCalChan] = {0.0}; // For easy TGraphErrors build
  
    for(Int_t i=0; i<hcal::maxHCalChan; i++){
      if( qreplay ){
	report << "Skipping matrix element checks for iteration 1 (quasi-replay).." << endl << endl;
	break;
      }else{
	report << "Proceeding to energy matrix quality check analysis.." << endl << endl;
      }

      badcell[i] = 0;
      y[i] = i;
    
      //Do not change ADC gain coeff if insufficient events or energy dep in cell
      if( gset[s].NEV[i] < minEventPerCell || gset[s].Ma(i,i) < highDelta*gset[s].ba(i) ){ 

	cellBad = 1;

	Double_t elemRatio = gset[s].Ma(i,i)/gset[s].ba(i);

	gset[s].ba(i) = 1.0;  // Set RHS vector for cell i to 1.0 
	gset[s].Ma(i,i) = 1.0; // Set diagonal element of Matrix M for cell i to 1.0 
      
	for(Int_t j=0; j<hcal::maxHCalChan; j++){
	  if( j != i ){
	    gset[s].Ma(i,j) = 0.0;
	    gset[s].Ma(j,i) = 0.0;
	  }
	}

	badcell[i] = 1;
	report << "Cell " << i << " bad" << endl;
	report << "Number of events in bad cell =" << gset[s].NEV[i] << endl;
	report << "Matrix element/vector element ratio =" << elemRatio << endl << endl;
      } 

      badcell_oneblock[i] = 0;

      //Do not change ADC gain coeff if insufficient events or energy dep in cell, oneblock
      if( gset[s].NEV_oneblock[i] < minEventPerCell || gset[s].Ma_oneblock(i,i) < highDelta*gset[s].ba_oneblock(i) ){ 

	cellBad = 1;

	Double_t elemRatio = gset[s].Ma_oneblock(i,i)/gset[s].ba_oneblock(i);

	gset[s].ba_oneblock(i) = 1.0;  // Set RHS vector for cell i to 1.0 
	gset[s].Ma_oneblock(i,i) = 1.0; // Set diagonal element of Matrix M for cell i to 1.0 
      
	for(Int_t j=0; j<hcal::maxHCalChan; j++){
	  if( j != i ){
	  
	    gset[s].Ma_oneblock(i,j) = 0.0; 
	    gset[s].Ma_oneblock(j,i) = 0.0;
	  }
	}
	badcell_oneblock[i] = 1;
	report << "Oneblock: cell" << i << " bad" << endl;
	report << "Oneblock: number of events in bad cell =" << gset[s].NEV_oneblock[i] << endl;
	report << "Oneblock: matrix element/vector element ratio =" << elemRatio << endl << endl;
      }   

      //Calculate error per block (element of A)
      gset[s].err[i] = sqrt(gset[s].err_ev[i]/gset[s].NEV[i]);
      gset[s].err_oneblock[i] = sqrt(gset[s].err_ev_oneblock[i]/gset[s].NEV_oneblock[i]);
    } //endloop over channels
  
  

    if( cellBad!=0 && !qreplay ) cout << "Bad cells detected on set " >> oset.timestamp[s] >> ". See report for details." << endl << endl;
    if( cellBad==0 && !qreplay ) report << "No bad cells detected on set " >> oset.timestamp[s] >> "." << endl << endl;

    //If iteration == 0, invert the matrix, solve for ratios
    if( !qreplay ){
      TMatrixD M_inv = gset[s].Ma.Invert();
      TMatrixD M_inv_oneblock = gset[s].Ma_oneblock.Invert();
      TVectorD Coeff = gset[s].M_inv*gset[s].ba; // Stays unmodified for reference
      TVectorD Coeff_oneblock = gset[s].M_inv_oneblock*gset[s].ba_oneblock; // Stays unmodified for reference

      for(Int_t i=0; i<hcal::maxHCalChan; i++){

	if(badcell_oneblock[i]==0){
	  gset[s].GCoeff_oneblock[i]=Coeff_oneblock[i];
	}else{
	  gset[s].GCoeff_oneblock[i]=allp[t].para[f].oadcgain[i];
	}

	if(badcell[i]==0&&Coeff[i]>0){ //Only update the coefficent if coefficient passes quality checks
	  gset[s].GCoeff[i]=Coeff[i];
	  gset[s].GCoeff_divide[i]=gset[s].GCoeff[i]/gset[s].GCoeff_oneblock[i];
	}else{
	  gset[s].GCoeff[i]=allp[t].para[f].oadcgain[i]; // If the cell calibration is bad, use the old coefficient
	  gset[s].GCoeff_divide[i]=-1.0;
	}
      }

      Double_t yErr[hcal::maxHCalChan] = {0.0};

      TGraphErrors *ccgraph = new TGraphErrors( hcal::maxHCalChan, y, GCoeff, yErr, err ); 
      ccgraph->GetXaxis()->SetLimits(0.0,288.0);  
      ccgraph->GetYaxis()->SetLimits(0.0,0.25);
      ccgraph->SetTitle(Form("Calibration Coefficients, %s",oset.timestamp[s]));
      ccgraph->GetXaxis()->SetTitle("Channel");
      ccgraph->GetYaxis()->SetTitle("Unitless");
      ccgraph->SetMarkerStyle(20); // idx 20 Circles, idx 21 Boxes
      ccgraph->Write(Form("constants_%d",s));

      TGraphErrors *ccgraph_oneBlock = new TGraphErrors( hcal::maxHCalChan, y, GCoeff_oneblock, yErr, err_oneblock ); 
      ccgraph_oneBlock->GetXaxis()->SetLimits(0.0,288.0);  
      ccgraph_oneBlock->GetYaxis()->SetLimits(0.0,0.25);
      ccgraph_oneBlock->SetTitle(Form("Calibration Coefficients One Block, %s",oset.timestamp[s]));
      ccgraph_oneBlock->GetXaxis()->SetTitle("Channel");
      ccgraph_oneBlock->GetYaxis()->SetTitle("Unitless");
      ccgraph_oneBlock->SetMarkerStyle(21);
      ccgraph_oneBlock->Write(Form("constants_oneblock_%d",s));

      TGraphErrors *ccgraph_divide = new TGraphErrors( hcal::maxHCalChan, y, GCoeff_divide, yErr, yErr ); 
      ccgraph_divide->GetXaxis()->SetLimits(0.0,288.0);  
      ccgraph_divide->SetTitle(Form("Calibration Coefficients / OneBlock Coeff, %s",oset.timestamp[s]));
      ccgraph_divide->GetXaxis()->SetTitle("Channel");
      ccgraph_divide->GetYaxis()->SetTitle("Unitless");
      ccgraph_divide->SetMarkerStyle(21);
      ccgraph_divide->Write(Form("constants_divide_%d",s));

      //Declare coeff outfile
      ofstream GainCoeff;

      GainCoeff.open( writeCoeffPath );

      //Console/txt outs
      Int_t cell = 0;

      GainCoeff << "#HCal gain coefficients from SBS-" << kine << " obtained " << date.c_str() << " timestamp set " << oset.timestamp[s] << endl;
      if( oset.timestamp[s].compare("none")!=0 )
	GainCoeff << oset.timestamp[s] << endl;
      GainCoeff << "sbs.hcal.adc.gain =" << endl;

      cout << "Gain Coefficients" << endl;
      report << "Gain Coefficients" << endl;
      for( Int_t r=0; r<kNrows; r++ ){
	for( Int_t c=0; c<kNcols; c++ ){
	  GainCoeff << gset[s].GCoeff[cell] << "  ";
	  cout << gset[s].GCoeff[cell] << "  ";
	  report << gset[s].GCoeff[cell] << "  ";
	  cell++;
	}
	GainCoeff << endl;
	cout << endl;
	report << endl;
      }

      cell = 0;

      cout << endl << "One Block Std:" << endl;
      report << endl << "One Block Std:" << endl;
      GainCoeff << endl << "#One Block Std = " << endl;

      for( Int_t r = 0; r<kNrows; r++){
	for( Int_t c = 0; c<kNcols; c++){
	  GainCoeff << gset[s].GCoeff_oneblock[cell] << "  ";
	  cout << gset[s].GCoeff_oneblock[cell] << "  ";
	  report << gset[s].GCoeff_oneblock[cell] << "  ";
	  cell++;
	}
	GainCoeff << endl;
	cout << endl;
      }

      cell = 0;

      cout << endl << "Total Number of events available for calibration: " << endl << endl;
      report << endl << "Total Number of events available for calibration: " << endl << endl;
      GainCoeff << endl << "#Number of events available for calibration = " << endl;

      for( Int_t r = 0; r<kNrows; r++){
	for( Int_t c = 0; c<kNcols; c++){
	  GainCoeff << gset[s].NEV[cell] << "  ";
	  cout << gset[s].NEV[cell] << "  ";
	  report << gset[s].NEV[cell] << "  ";
	  cell++;
	}
	GainCoeff << endl;
	cout << endl;
	report << endl;
      }

      GainCoeff.close();
    } //end if iter==0

  }//endloop over calibration sets

  fout->Write();

  st->Stop();

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;
  report << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

  report.close();

  cout << endl << endl << "Detailed processing report can be found here: " << logPath << "." << endl;

}
