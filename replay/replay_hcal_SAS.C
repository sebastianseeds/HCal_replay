// sseeds 11.18.21 HCAL expert replay designed to run more quickly and to perform HCAL standalone analysis.
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

void replay_hcal_SAS(int run_number = 124, uint nev = -1, uint nseg = 0)
{
  //load SBS-offline
  gSystem->Load("libsbs.so");
  
  //Add BB information for energy calibration
  SBSBigBite* bigbite = new SBSBigBite("bb", "BigBite spectrometer" );
  SBSBBTotalShower* ts= new SBSBBTotalShower("ts", "sh", "ps", "BigBite shower");
  ts->SetDataOutputLevel(0);
  bigbite->AddDetector( ts );
  ts->SetStoreEmptyElements(kFALSE);
  ts->GetShower()->SetStoreEmptyElements(kFALSE);
  ts->GetPreShower()->SetStoreEmptyElements(kFALSE);

  SBSGenericDetector* tdctrig= new SBSGenericDetector("tdctrig","BigBite shower TDC trig");
  tdctrig->SetModeADC(SBSModeADC::kNone);
  tdctrig->SetModeTDC(SBSModeTDC::kTDC);
  bigbite->AddDetector( tdctrig );
  gHaApps->Add(bigbite);

  SBSEArm *harm = new SBSEArm("sbs","Hadron Arm with HCal");
  SBSHCal* hcal =  new SBSHCal("hcal","HCAL");
  hcal->SetStoreRawHits(kTRUE);
  hcal->SetStoreEmptyElements(kTRUE);
  //hcal->SetDataOutputLevel(1);
  harm->AddDetector(hcal);

  SBSGenericDetector* sbstrig= new SBSGenericDetector("trig","HCal trigs");
  sbstrig->SetModeADC(SBSModeADC::kWaveform);
  sbstrig->SetStoreRawHits(kTRUE);
  //trig->SetDataOutputLevel(1);
  harm->AddDetector( sbstrig );
  
  gHaApps->Add( harm );

  gHaPhysics->Add( new THaPrimaryKine( "e.kine", "electron kinematics", "bb", 0.0, 0.938272 ));
  //--- Set up the run we want to replay ---
  
  // This often requires a bit of coding to search directories, test
  // for non-existent files, etc.
  //Variables for searching for split data files.
  //int seg = 0;
  bool seg_ok = true;
  TString exp = "hcal_trim";
  // Create file name patterns.
  //string firstname = "hcal_trigtest_%d";
  //string firstname = "e1209019_trigtest_%d";
  string firstname = "e1209019_%d";
  //string firstname = "hcal_adc_tdc_%d";
  
  THaAnalyzer* analyzer = new THaAnalyzer;
  
  while(seg_ok) {
    //string endname = Form(".evio.%d",nseg);
    string endname = Form(".evio.0.%d",nseg);
    //string endname = Form(".dat.%d",nseg);
    //string endname = Form(".evio");
    string combined(string(firstname)+endname);
    const char* RunFileNamePattern = combined.c_str();
    vector<TString> pathList;
    pathList.push_back(".");
    //string data_path = "/adaqfs/home/a-onl/sbs/data/";
    string data_path = "/adaqeb1/data1/";
    //pathList.push_back(Form("%s/data","/adaqfs/home/a-onl/sbs"));
    pathList.push_back(Form("%s/data1","/adaqeb1"));
    //pathList.push_back(Form("%s/data","/adaqfs/home/a-onl/skbarcus"));
    
    // Check if segment exits
    string run_name = Form(RunFileNamePattern, run_number);
    string seg_path_str = data_path+run_name;
    //Convert the std::string to const char * pointer expected by gSystem->AccessPathName.
    const char * seg_path = seg_path_str.c_str();
    std::cout<<seg_path<<std::endl;
    //if( gSystem->AccessPathName("/adaqfs/home/a-onl/sbs/data/hcal_trigtest_183.evio.0")) {
    
    if( gSystem->AccessPathName(seg_path)) {
      seg_ok = false;
      std::cout << "Segment " << nseg << " not found. Exiting" << std::endl;
      continue;
    }
    else{
      std::cout<<"Found "<<Form(RunFileNamePattern,run_number)<<"."<<std::endl;
    }
    
    THaRun* run = new THaRun( pathList, Form(RunFileNamePattern, run_number) );
    
    run->SetDataRequired(7);//for the time being
    
    cout << "Number of events to replay (-1=all)? ";
    if( nev > 0 )
      run->SetLastEvent(nev);
    
    //--- Set up any physics calculations we want to do ---
    
    // Extract the reconstructed target quantities of the golden track
    // Not really a physics calculation, but a convenience function.
    // It effectively converts the L.tr.* variables, which are arrays, 
    // to scalers L.gold.*
    
    //gHaPhysics->Add( new THaGoldenTrack( "R.gold", "RHRS golden track", "R" ));
    
    // Single-arm electron kinematics for the one spectrometer we have set up.
    // We assume a carbon-12 target (12 AMU)
    //gHaPhysics->Add( new THaPrimaryKine( "R.ekine", "RHRS electron kinematics",
    //"R", 0.511e-3, 12*0.9315 ));
    
    // Vertex position calculated from RHRS golden track and ideal beam
    // (will poor resolution if raster is on)
    //gHaPhysics->Add( new THaReactionPoint( "R.vx", "Vertex R", "R", "IB" ));
    
    //--- Define what we want the analyzer to do ---
    // The only mandatory items are the output definition and output file names
    
    //THaAnalyzer* analyzer = new THaAnalyzer;
    
    //TString out_dir = gSystem->Getenv("OUT_DIR");
    TString out_dir ="/adaqfs/home/a-onl/sbs/Rootfiles";
    if( out_dir.IsNull() )  out_dir = ".";
    TString out_file = out_dir + "/" + exp + Form("_%d_%d.root", run_number,nev);
    
    analyzer->SetOutFile( out_file );

    analyzer->SetOdefFile( "/adaqfs/home/a-onl/sbs/HCal_replay/replay/replay_hcal.odef" );
    analyzer->SetCutFile( "/adaqfs/home/a-onl/sbs/HCal_replay/replay/replay_hcal.cdef" );

    analyzer->SetVerbosity(2);  // write cut summary to stdout
    analyzer->EnableBenchmarks();
    
    //--- Analyze the run ---
    // Here, one could replay more than one run with
    // a succession of Process calls. The results would all go into the
    // same ROOT output file
    
    run->Print();
    
    
    analyzer->Process(run);
    nseg++;
    }

  // Clean up

  analyzer->Close();
  delete analyzer;
  //gHaCuts->Clear();
  gHaVars->Clear();
  gHaPhysics->Delete();
  gHaApps->Delete();

  // Open the ROOT file so that a user doing interactive analysis can 
  // immediately look at the results. Of course, don't do this in batch jobs!
  // To close the file later, simply type "delete rootfile" or just exit

  //TFile* rootfile = new TFile( out_file, "READ" );
}
