#include <iostream>

#include "TSystem.h"
#include "TString.h"
#include "TFile.h"

#include "THaShower.h"
#include "THaEvent.h"
#include "THaEvData.h"
#include "THaRun.h"
#include "THaAnalyzer.h"
#include "THaVarList.h"

#include "SBSBigBite.h"
#include "SBSBBTotalShower.h"
#include "SBSBBShower.h"
#include "SBSTimingHodoscope.h"

void replay_BBHodoRaw(int run_number = 10746, uint nev = 50000, TString exp = /*"bbhodo"*/ "e1209019")
{
  uint nseg = 0;
  
  SBSBigBite* bigbite = new SBSBigBite("bb", "BigBite spectrometer" );
  SBSTimingHodoscope* hodotdc = new  SBSTimingHodoscope("hodotdc", "BigBite hodo");
  hodotdc->SetModeTDC(SBSModeTDC::kTDC);
  hodotdc->SetModeADC(SBSModeADC::kNone);
  hodotdc->SetStoreRawHits(true);
  SBSTimingHodoscope* hodoadc = new  SBSTimingHodoscope("hodoadc", "BigBite hodo");
  hodoadc->SetModeTDC(SBSModeTDC::kNone);
  hodoadc->SetModeADC(SBSModeADC::kADC);
  bigbite->AddDetector(hodotdc);
  bigbite->AddDetector(hodoadc);
  gHaApps->Add(bigbite);

  //--- Set up the run we want to replay ---

  // This often requires a bit of coding to search directories, test
  // for non-existent files, etc.
  // Create file name patterns.
  string firstname = /*"bbhodo_%d"*/ "e1209019__%d";
  string endname = Form(".evio.%d",nseg);
  string combined(string(firstname)+endname);
  const char* RunFileNamePattern = combined.c_str();
  vector<TString> pathList;
  pathList.push_back(".");
  pathList.push_back(Form("%s/data1/","/adaqeb1"));


  THaAnalyzer* analyzer = new THaAnalyzer;
  THaEvent* event = new THaEvent;
  THaRun* run = 0;
  int seg = 0;
  bool seg_ok = true;
  while(seg_ok) {
    TString data_fname;
    data_fname = TString::Format("%s/%s_%d.evio.0.%d",getenv("DATA_DIR"),exp.Data(),run_number,seg);
    std::cout << "Looking for segment " << seg << " file " << data_fname.Data() << std::endl;
    if( gSystem->AccessPathName(data_fname)) {
      seg_ok = false;
      std::cout << "Segment " << seg << " not found. Exiting" << std::endl;
      continue;
    }
    run = new THaRun(data_fname);
    run->SetLastEvent(nev);

    run->SetDataRequired(0);//for the time being
    run->SetDate(TDatime());

    analyzer->SetEvent( event );
    TString out_dir = gSystem->Getenv("OUT_DIR");
    if( out_dir.IsNull() )  out_dir = ".";
    TString out_file = out_dir + "/" + exp + Form("_%d_%d.root", run_number,nev);

    analyzer->SetOutFile( out_file );
    
    analyzer->SetCutFile( "replay_hodo.cdef" );
    analyzer->SetOdefFile( "replay_hodo.odef" );
    analyzer->SetVerbosity(2);  // write cut summary to stdout
    analyzer->EnableBenchmarks();

    //--- Analyze the run ---
    // Here, one could replay more than one run with
    // a succession of Process calls. The results would all go into the
    // same ROOT output file
    
    run->Print();
    analyzer->Process(run);
    // Cleanup this run segment
    delete run;
    
    seg++; // Increment for next search
  }
  


  analyzer->Close();
  delete analyzer;
  gHaVars->Clear();
  gHaPhysics->Delete();
  gHaApps->Delete();

}
