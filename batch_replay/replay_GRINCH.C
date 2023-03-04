#if !defined(__CLING__) || defined(__ROOTCLING__)
R__ADD_INCLUDE_PATH($SBS/include)
R__ADD_LIBRARY_PATH($SBS/lib64)
R__ADD_LIBRARY_PATH($SBS/lib)
R__LOAD_LIBRARY(libsbs)

#include <iostream>

#include "TSystem.h"
#include "TString.h"
#include "TFile.h"
#include "TList.h"
#include "TObject.h"

#include "THaEvData.h"
#include "THaRun.h"
#include "THaEvent.h"
#include "THaAnalyzer.h"
#include "THaVarList.h"
#include "THaInterface.h"
#include "THaApparatus.h"

#include "SBSBigBite.h"
#include "SBSGRINCH.h"
#include "SBSGRINCH_ClusterList.h"

#include "SBSEArm.h"
#include "SBSHCal.h"

#endif

using namespace std;

//void replay_GRINCH(Int_t runnum = 10460, Int_t lastEvent = -1){
void replay_GRINCH(Int_t runnum = 62, Int_t lastEvent = -1){
  //gSystem->Load(Form("%s/lib64/libsbs.so",getenv("SBSOFFLINE")));
  SBSGenericDetector *grinch_tdc = new SBSGenericDetector("grinch_tdc","GRINCH TDC data");
  SBSGenericDetector *grinch_adc = new SBSGenericDetector("grinch_adc","GRINCH ADC data");
  //grinch2->SetDisableRefADC(true);
  grinch_adc->SetModeADC(SBSModeADC::kWaveform);
  grinch_adc->SetModeTDC(SBSModeTDC::kNone);

  grinch_tdc->SetModeTDC(SBSModeTDC::kTDC);
  //grinch_tdc->SetModeTDC(SBSModeTDC::kCommonStartTDC);
  grinch_tdc->SetModeADC(SBSModeADC::kNone);
  grinch_tdc->SetDisableRefTDC(true);

   //grinch_tdc->SetWithTDC(true);

 SBSBigBite* bb = new SBSBigBite("bb", "BigBite spectrometer" );
  bb->AddDetector(grinch_adc);
  bb->AddDetector(grinch_tdc);

  // Set up the analyzer - we use the standard one,
  // but this could be an experiment-specific one as well.
  // The Analyzer controls the reading of the data, executes
  // tests/cuts, loops over Apparatus's and PhysicsModules,
  // and executes the output routines.
  THaAnalyzer* analyzer = new THaAnalyzer;

  gHaApps->Add(bb);

  // A simple event class to be output to the resulting tree.
  // Creating your own descendant of THaEvent is one way of
  // defining and controlling the output.
  THaEvent* event = new THaEvent;

  // Define the run(s) that we want to analyze.
  // We just set up one, but this could be many.
  //  THaRun* run = new THaRun( "prod12_4100V_TrigRate25_4.dat" );
//  THaRun* run = new THaRun(TString::Format("/cache/halla/sbs/raw/grinch_%d.evio.0",runnum));
  THaRun* run = new THaRun(TString::Format("/adaq1/data1/sbs/grinch_%d.evio.0",runnum));
//  THaRun* run = new THaRun(TString::Format("/adaq1/data1/sbs/e1209019_trigtest_%d.evio.0",runnum));
  run->SetLastEvent(lastEvent);

  run->SetDataRequired(0);
//  run->SetDataRequired(7);
  run->SetDate(TDatime());

  analyzer->SetVerbosity(2);
  analyzer->SetCutFile( "replay_grinch.cdef" );
  //analyzer->SetOdefFile("output_grinch.def");
  analyzer->SetOdefFile("replay_grinch.odef");
  analyzer->SetMarkInterval(500);

  // Define the analysis parameters
  analyzer->SetEvent( event );
  //analyzer->SetOutFile( TString::Format("rootfiles/fadcAERO_%d.root",runnum));
  //analyzer->SetOutFile( TString::Format("../ReplayedFiles/replayed_bbgrinch_%d.evio.0.root",runnum));
  analyzer->SetOutFile( TString::Format("../Rootfiles/grinch_%d.evio.0.root",runnum));
//analyzer->SetOutFile( TString::Format("../Rootfiles/e1209019_trigtest_%d.evio.0.root",runnum));
  // File to record cuts accounting information
  //analyzer->SetSummaryFile(TString::Format("sbs_hcal_led_%d.log",runnum)); // optional

  //analyzer->SetCompressionLevel(0); // turn off compression
  analyzer->Process(run);     // start the actual analysis
}
