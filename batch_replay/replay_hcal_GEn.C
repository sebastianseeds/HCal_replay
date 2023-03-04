// sseeds 11.13.22 Rework of GEn Replay to fully replay all segments by stream
// Should be run with shell script run_hcal_gen_replay.sh
// See shell script for usage
#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <iostream>

#include "TSystem.h"
#include "THaGlobals.h"
#include "TString.h"
#include "TFile.h"
#include "TList.h"
#include "TObject.h"
#include "TClonesArray.h"

#include "THaEvData.h"
#include "THaEvent.h"
#include "THaRun.h"
#include "THaAnalyzer.h"
#include "THaVarList.h"
#include "THaInterface.h"
#include "THaGoldenTrack.h"
#include "THaDecData.h"
#include "THaPrimaryKine.h"

#include "SBSBigBite.h"
#include "SBSBBShower.h"
#include "SBSBBTotalShower.h"
#include "SBSGRINCH.h"
#include "SBSEArm.h"
#include "SBSHCal.h"
#include "SBSGEMStand.h"
#include "SBSTimingHodoscope.h"
#include "SBSGEMSpectrometerTracker.h"
#include "SBSGEMTrackerBase.h"
#include "SBSRasteredBeam.h"
#endif

void replay_hcal_GEn(int run_number = 124, uint nev = -1, uint nseg = 0, Int_t pedestalmode=0, Int_t stream=0)
{
  //Add BB information for energy calibration
  SBSBigBite* bigbite = new SBSBigBite("bb", "BigBite spectrometer" );
  SBSBBTotalShower* ts= new SBSBBTotalShower("ts", "sh", "ps", "BigBite shower");
  ts->SetDataOutputLevel(0);
  bigbite->AddDetector( ts );
  ts->SetStoreEmptyElements(kFALSE);
  ts->GetShower()->SetStoreEmptyElements(kFALSE);
  ts->GetPreShower()->SetStoreEmptyElements(kFALSE);

  SBSGenericDetector* bbtrig= new SBSGenericDetector("bbtrig","BigBite shower ADC trig");
  bbtrig->SetModeADC(SBSModeADC::kADC);
  bbtrig->SetModeTDC(SBSModeTDC::kTDC);
  bbtrig->SetStoreEmptyElements(kFALSE);
  bigbite->AddDetector( bbtrig );

  SBSGenericDetector* tdctrig= new SBSGenericDetector("tdctrig","BigBite shower TDC trig");
  tdctrig->SetModeADC(SBSModeADC::kNone);
  tdctrig->SetModeTDC(SBSModeTDC::kTDC);
  tdctrig->SetStoreEmptyElements(kFALSE);
  bigbite->AddDetector( tdctrig );

  SBSTimingHodoscope* hodotdc = new  SBSTimingHodoscope("hodotdc", "BigBite hodo");
  hodotdc->SetModeTDC(SBSModeTDC::kTDC);
  hodotdc->SetModeADC(SBSModeADC::kNone);
  hodotdc->SetStoreEmptyElements(kFALSE);
  hodotdc->SetDataOutputLevel(1);// => this adds in the output the elements belonging to the "main" cluster.
  
  SBSTimingHodoscope* hodoadc = new  SBSTimingHodoscope("hodoadc", "BigBite hodo");
  hodoadc->SetModeTDC(SBSModeTDC::kNone);
  hodoadc->SetModeADC(SBSModeADC::kADC);
  hodoadc->SetStoreEmptyElements(kFALSE);
  hodoadc->SetStoreRawHits(kFALSE);
  hodotdc->SetDataOutputLevel(0);
  bigbite->AddDetector(hodotdc);
  bigbite->AddDetector(hodoadc);
  SBSGEMSpectrometerTracker *bbgem = new SBSGEMSpectrometerTracker("gem", "BigBite Hall A GEM data");
  bool pm =  ( pedestalmode != 0 );
  //this will override the database setting:
  ( static_cast<SBSGEMTrackerBase *> (bbgem) )->SetPedestalMode( pm );
  bigbite->AddDetector(bbgem);
  gHaApps->Add(bigbite);

  SBSEArm *harm = new SBSEArm("sbs","Hadron Arm with HCal");
  SBSHCal* hcal =  new SBSHCal("hcal","HCAL");
  hcal->SetStoreRawHits(kTRUE);
  hcal->SetStoreEmptyElements(kTRUE);
  hcal->SetDataOutputLevel(1);
  harm->AddDetector(hcal);

  SBSGenericDetector* sbstrig= new SBSGenericDetector("trig","HCal trigs");
  sbstrig->SetModeADC(SBSModeADC::kWaveform);
  sbstrig->SetStoreRawHits(kTRUE);
  //trig->SetDataOutputLevel(1);
  harm->AddDetector( sbstrig );
  
  gHaApps->Add( harm );

  gHaPhysics->Add( new THaGoldenTrack( "BB.gold", "BigBite golden track", "bb" ));
  gHaPhysics->Add( new THaPrimaryKine( "e.kine", "electron kinematics", "bb", 0.0, 0.938272 ));
  //--- Set up the run we want to replay ---
  // This often requires a bit of coding to search directories, test
  // for non-existent files, etc.
  // Variables for searching for split data files.
  bool seg_ok = true;
  bool path_ok = true;
  TString exp = "hcal_gen_pass2";
  // Create file name patterns.
  string firstname = "/e1209016_%d";
  TString prefix = gSystem->Getenv("DATA_DIR");

  //std::cout<<"!!!!!!!!!!!!!!!!!"<<prefix<<std::endl;

  THaAnalyzer* analyzer = new THaAnalyzer;

  //Look for datafile
  string data_path;
  string searchname = Form(".evio.%d.0",stream);
  string searchcombined(string(firstname)+searchname);
  const char* RunFileNamePattern = searchcombined.c_str();
  vector<TString> pathlist;
  pathlist.push_back(prefix);
  
  if( prefix != "/adaqeb1/data1" )
    pathlist.push_back( "/adaqeb1/data1" );
  
  if( prefix != "/adaqeb2/data1" )
    pathlist.push_back( "/adaqeb2/data1" );
  
  if( prefix != "/adaqeb3/data1" )
    pathlist.push_back( "/adaqeb3/data1" );
  
  if( prefix != "/adaq1/data1/sbs" )
    pathlist.push_back( "/adaq1/data1/sbs" );
  
  if( prefix != "/cache/halla/sbs/raw" )
    pathlist.push_back( "/cache/halla/sbs/GEnII/raw" );

  for(int i=0; i<pathlist.size(); i++){
    TString run_name = Form(RunFileNamePattern, run_number);
    TString path = pathlist[i];
    TString seg_path_str = path+run_name;

    string segment_path = (string)seg_path_str;

    //Convert the std::string to const char * pointer expected by gSystem->AccessPathName.
    //const char * seg_path = seg_path_str.c_str();
    const char * seg_path = segment_path.c_str();
    std::cout << "Searching path " << seg_path << ".." << std::endl;
    
    if( gSystem->AccessPathName(seg_path)) {
      std::cout << "No data found in " << pathlist[i] << "." << std::endl;
    }
    else{
      std::cout<< "Found " << Form(RunFileNamePattern,run_number) << "!" << std::endl;   
      std::cout<< "Setting default path.." <<std::endl;
      data_path = pathlist[i];
      break;
    }
    
  }
  while(seg_ok) {
    string endname = Form(".evio.%d.%d",stream,nseg);
    string combined(string(firstname)+endname);
    const char* RunFileNamePattern = combined.c_str();
    // Check if segment exits
    string run_name = Form(RunFileNamePattern, run_number);
    string seg_path_str = data_path+run_name;
    const char * seg_path = seg_path_str.c_str();
    std::cout<<seg_path<<std::endl;
    
    if( gSystem->AccessPathName(seg_path)) {
      seg_ok = false;
      std::cout << "Segment " << nseg << " not found. Exiting.." << std::endl;
      continue;
    }
    else{
      std::cout<<"Found "<<Form(RunFileNamePattern,run_number)<<"!"<<std::endl;
    }
    
    vector<TString> path;
    path.push_back(".");
    path.push_back(data_path);

    THaRun* run = new THaRun( path, Form(RunFileNamePattern, run_number) );
    
    run->SetDataRequired(7);//for the time being
    
    cout << "Number of events to replay (-1=all)? ";
    if( nev > 0 )
      run->SetLastEvent(nev);
    
    //TString out_dir = gSystem->Getenv("OUT_DIR");
    TString out_dir ="/adaqfs/home/a-onl/sbs/Rootfiles";
    if( out_dir.IsNull() )  out_dir = ".";
    TString out_file = out_dir + "/" + exp + Form("_r%d_nev%d_stream%d.root", run_number,nev,stream);
    
    analyzer->SetOutFile( out_file );

    analyzer->SetOdefFile( "replay_hcal.odef" );
    analyzer->SetCutFile( "replay_hcal.cdef" );
    // Set EPICS event type
    analyzer->SetEpicsEvtType(150);
    analyzer->SetEpicsEvtType(151);

    analyzer->SetVerbosity(2);  // write cut summary to stdout
    analyzer->EnableBenchmarks();
    
    run->Print();
       
    analyzer->Process(run);
    nseg++;
    }

  // Clean up
  analyzer->Close();
  delete analyzer;
  //gHaVars->Clear();
  gHaPhysics->Delete();
  gHaApps->Delete();
}
