// sseedso 1.13.23 HCAL expert replay designed to run more quickly and to perform HCAL standalone analysis.
//#if !defined(__CLING__) || defined(__ROOTCLING__)
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
//#endif

void replay_hcal_GMn(UInt_t runnum=10491, Long_t nevents=-1, Long_t firstevent=0, const char *fname_prefix="e1209019", UInt_t firstsegment=0, UInt_t maxsegments=1, Int_t pedestalmode=0, Int_t cmplots=1)
{

  THaAnalyzer* analyzer = new THaAnalyzer;

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
  hodotdc->SetDisableRefTDC(true);
  hodotdc->SetStoreRawHits(kTRUE);
  hodotdc->SetStoreEmptyElements(kFALSE);
  
  //hodotdc->SetDataOutputLevel(1);// => this adds in the output the elements belonging to the "main" cluster.
  
  SBSTimingHodoscope* hodoadc = new  SBSTimingHodoscope("hodoadc", "BigBite hodo");
  hodoadc->SetModeTDC(SBSModeTDC::kNone);
  hodoadc->SetModeADC(SBSModeADC::kADC);
  hodoadc->SetStoreEmptyElements(kFALSE);
  hodoadc->SetStoreRawHits(kFALSE);
  hodotdc->SetDataOutputLevel(0);
  //bigbite->AddDetector( new THaShower("ps", "BigBite preshower") );
  bigbite->AddDetector(hodotdc);
  bigbite->AddDetector(hodoadc);

  SBSGenericDetector *grinch_adc = new SBSGenericDetector("grinch_adc","GRINCH ADC data");
  //grinch2->SetDisableRefADC(true);
  grinch_adc->SetModeADC(SBSModeADC::kADC);
  grinch_adc->SetModeTDC(SBSModeTDC::kNone);
  grinch_adc->SetStoreEmptyElements(kFALSE);
  grinch_adc->SetStoreRawHits(kTRUE); //kTRUE to get multiple hits from the ADC

  SBSGenericDetector *grinch_tdc = new SBSGenericDetector("grinch_tdc","GRINCH TDC data");
  grinch_tdc->SetModeTDC(SBSModeTDC::kTDC);
  //grinch_tdc->SetModeTDC(SBSModeTDC::kCommonStartTDC);
  grinch_tdc->SetModeADC(SBSModeADC::kNone);
  grinch_tdc->SetDisableRefTDC(true);
  grinch_tdc->SetStoreEmptyElements(kFALSE);
  grinch_tdc->SetStoreRawHits(kTRUE); //kTRUE to get multiple hits from the TDC
  bigbite->AddDetector(grinch_adc);
  bigbite->AddDetector(grinch_tdc);

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
  //hcal->SetDataOutputLevel(1);
  harm->AddDetector(hcal);

  SBSGenericDetector* sbstrig= new SBSGenericDetector("trig","HCal trigs");
  sbstrig->SetModeADC(SBSModeADC::kWaveform);
  sbstrig->SetStoreRawHits(kTRUE);
  //trig->SetDataOutputLevel(1);
  harm->AddDetector( sbstrig );
  
  gHaApps->Add( harm );

  // add decoder
  THaApparatus* decL = new THaDecData("DL","Misc. Decoder Data");
  gHaApps->Add( decL );
  
  // add *rastered* beam
  THaApparatus* Lrb = new SBSRasteredBeam("Lrb","Raster Beamline for FADC");
  gHaApps->Add(Lrb);
  
  THaApparatus* sbs = new SBSRasteredBeam("SBSrb","Raster Beamline for FADC");
  gHaApps->Add(sbs);

  gHaPhysics->Add( new THaGoldenTrack( "BB.gold", "BigBite golden track", "bb" ));
  gHaPhysics->Add( new THaPrimaryKine( "e.kine", "electron kinematics", "bb", 0.0, 0.938272 ));

  LHRSScalerEvtHandler *lScaler = new LHRSScalerEvtHandler("Left","HA scaler event type 140");
  // lScaler->SetDebugFile(&debugFile);
  gHaEvtHandlers->Add(lScaler);

  SBSScalerEvtHandler *sbsScaler = new SBSScalerEvtHandler("sbs","SBS Scaler Bank event type 1");
  //sbsScaler->AddEvtType(1);             // Repeat for each event type with scaler banks
  sbsScaler->SetUseFirstEvent(kTRUE);
  gHaEvtHandlers->Add(sbsScaler);
   
  THaEvent* event = new THaEvent;

  TString prefix = gSystem->Getenv("DATA_DIR");
  
  bool segmentexists = true;
  int segment=firstsegment; 

  int lastsegment=firstsegment;
  
  TDatime now = TDatime();
  
  int stream = 0;
 
  TClonesArray *filelist = new TClonesArray("THaRun",10);

  vector<TString> pathlist;
  pathlist.push_back( prefix );

  if( prefix != "/adaqeb1/data1" )
    pathlist.push_back( "/adaqeb1/data1" );

  if( prefix != "/adaqeb2/data1" )
    pathlist.push_back( "/adaqeb2/data1" );

  if( prefix != "/adaqeb3/data1" )
    pathlist.push_back( "/adaqeb3/data1" );

  if( prefix != "/adaq1/data1/sbs" )
    pathlist.push_back( "/adaq1/data1/sbs" );

  if( prefix != "/cache/mss/halla/sbs/raw" )
    pathlist.push_back( "/cache/mss/halla/sbs/raw" );

  for( int i=0; i<pathlist.size(); i++ ){
    cout << "search paths = " << pathlist[i] << endl;
  }

  TDatime RunDate = TDatime(); 

  int max1 = maxsegments;

  int segcounter=0;
  
  //This loop adds all file segments found to the list of THaRuns to process:
  while( segcounter < max1 && segment - firstsegment < maxsegments ){
    
    TString codafilename;
    codafilename.Form("%s_%d.evio.%d.%d", fname_prefix, runnum, stream, segment );

    TString ftest(fname_prefix);
    if( ftest == "bbgem" || ftest == "e1209019_trigtest" ){
      codafilename.Form("%s_%d.evio.%d", fname_prefix, runnum, segment );
    }

    segmentexists = false;

    cout << "codafilename = " << codafilename << endl;

    for( int ipath=0; ipath<pathlist.size(); ipath++ ){
      TString searchname;
      searchname.Form( "%s/%s", pathlist[ipath].Data(), codafilename.Data() );
      
      if( !gSystem->AccessPathName( searchname.Data() ) ){
	segmentexists = true;
	break;
      }
    }
   
    if( segmentexists ){
      new( (*filelist)[segcounter] ) THaRun( pathlist, codafilename.Data(), "GMN run" );
      cout << "Added segment " << segment << ", CODA file name = " << codafilename << endl;

    }

    if( segmentexists ){
      segcounter++;
      lastsegment = segment;
    }
    segment++;
  }

  cout << "n segments to analyze = " << segcounter << endl;

  prefix = gSystem->Getenv("OUT_DIR");

  TString outfilename;

  if( nevents > 0 ){ 

    outfilename.Form( "%s/hcal_gmn_replayed_%d_stream%d_seg%d_%d_firstevent%d_nevent%d.root", prefix.Data(), runnum,
		      stream, firstsegment, lastsegment, firstevent, nevents );
  } else {
    outfilename.Form( "%s/hcal_gmn_fullreplay_%d_stream%d_seg%d_%d.root", prefix.Data(), runnum,
		      stream, firstsegment, lastsegment );
  }
 

  analyzer->SetVerbosity(2);
  analyzer->SetMarkInterval(100);

  analyzer->EnableBenchmarks();
  
  // Define the analysis parameters
  analyzer->SetEvent( event );
  analyzer->SetOutFile( outfilename.Data() );
  // File to record cuts accounting information
  
  prefix = gSystem->Getenv("LOG_DIR");
  analyzer->SetSummaryFile(Form("%s/replay_hcal_gmn_%d_stream%d_seg%d_%d.log", prefix.Data(), runnum, 
				stream, firstsegment, lastsegment));
  
  prefix = gSystem->Getenv("SBS_REPLAY");
  prefix += "/replay/";

  analyzer->SetOdefFile( "/w/halla-scshelf2102/sbs/seeds/HCal_replay/replay/replay_hcal.odef" );
  analyzer->SetCutFile( "/w/halla-scshelf2102/sbs/seeds/HCal_replay/replay/replay_hcal.cdef" );

  filelist->Compress();

  for( int iseg=0; iseg<filelist->GetEntries(); iseg++ ){
    THaRun *run = ( (THaRun*) (*filelist)[iseg] );
    if( nevents > 0 ) run->SetLastEvent(firstevent+nevents-1); //not sure if this will work as we want it to for multiple file segments chained together

    run->SetFirstEvent( firstevent );
    
    run->SetDataRequired(THaRunBase::kDate|THaRunBase::kRunNumber);

    run->Init();
    
    if( run->GetSegment() >= firstsegment && run->GetSegment() - firstsegment < maxsegments ){
      analyzer->Process(run);     // start the actual analysis
    }
  }
}




  bool seg_ok = true;
  bool path_ok = true;
  TString exp = "hcal_general";
  // Create file name patterns.
  string firstname = "/e1209019_%d";
  TString prefix = gSystem->Getenv("DATA_DIR");
  //THaAnalyzer* analyzer = new THaAnalyzer;

  //Look for datafile
  string data_path;
  string searchname = ".evio.0.0";
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
    pathlist.push_back( "/cache/halla/sbs/raw" );

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
    string endname = Form(".evio.0.%d",nseg);
    string combined(string(firstname)+endname);
    const char* RunFileNamePattern = combined.c_str();
    /*
    vector<TString> pathlist;
    pathlist.push_back("prefix");

    if( prefix != "/adaqeb1/data1" )
      pathlist.push_back( "/adaqeb1/data1" );

    if( prefix != "/adaqeb2/data1" )
      pathlist.push_back( "/adaqeb2/data1" );

    if( prefix != "/adaqeb3/data1" )
      pathlist.push_back( "/adaqeb3/data1" );

    if( prefix != "/adaq1/data1/sbs" )
      pathlist.push_back( "/adaq1/data1/sbs" );

    if( prefix != "/cache/halla/sbs/raw" )
      pathlist.push_back( "/cache/halla/sbs/raw" );

    //string data_path = "/adaqfs/home/a-onl/sbs/data/";
    //string data_path = "/adaqeb1/data1/";
    //pathlist.push_back(Form("%s/data","/adaqfs/home/a-onl/sbs"));
    //pathlist.push_back(Form("%s/data1","/adaqeb1"));
    //pathlist.push_back(Form("%s/data","/adaqfs/home/a-onl/skbarcus"));
    */
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
    
    TString out_dir = gSystem->Getenv("OUT_DIR");
    //TString out_dir ="/adaqfs/home/a-onl/sbs/Rootfiles";
    if( out_dir.IsNull() )  out_dir = ".";
    TString out_file = out_dir + "/" + exp + Form("_%d_%d.root", run_number,nev);
    
    analyzer->SetOutFile( out_file );

    analyzer->SetOdefFile( "/w/halla-scshelf2102/sbs/seeds/HCal_replay/replay/replay_hcal.odef" );
    analyzer->SetCutFile( "/w/halla-scshelf2102/sbs/seeds/HCal_replay/replay/replay_hcal.cdef" );

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
