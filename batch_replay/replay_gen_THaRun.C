// Replay script for SBS-GEn
// October 2022
//
// Requires Podd 1.7.4

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
#include "MultiFileRun.h"
#include "THaAnalyzer.h"
#include "THaVarList.h"
#include "THaInterface.h"
#include "THaGoldenTrack.h"
#include "THaPrimaryKine.h"
#include "THaDecData.h"

#include "SBSBigBite.h"
#include "SBSBBShower.h"
#include "SBSBBTotalShower.h"
#include "SBSGRINCH.h"
//#include "SBSCherenkovDetector.h"
#include "SBSEArm.h"
#include "SBSHCal.h"
#include "SBSGEMStand.h"
#include "SBSTimingHodoscope.h"
#include "SBSGEMSpectrometerTracker.h"
#include "SBSGEMTrackerBase.h"
#include "SBSRasteredBeam.h"
#include "LHRSScalerEvtHandler.h"
#include "SBSScalerEvtHandler.h"
//#endif

using namespace std;

void replay_gen(UInt_t runnum=10491, Long_t nevents=-1, Long_t firstevent=1, const char *fname_prefix="e1209016", UInt_t firstsegment=0, UInt_t maxsegments=1, Int_t maxstream=2, Int_t pedestalmode=0, Int_t cmplots=1)
{

  THaAnalyzer* analyzer = new THaAnalyzer;

  SBSBigBite* bigbite = new SBSBigBite("bb", "BigBite spectrometer" );
  //bigbite->AddDetector( new SBSBBShower("ps", "BigBite preshower") );
  //bigbite->AddDetector( new SBSBBShower("sh", "BigBite shower") );
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

  //SBSGenericDetector *grinch_tdc = new SBSGenericDetector("grinch_tdc","GRINCH TDC data");
  //SBSCherenkovDetector *grinch_tdc = new SBSCherenkovDetector("grinch_tdc","GRINCH TDC data");
  SBSGRINCH *grinch_tdc = new SBSGRINCH("grinch_tdc","GRINCH TDC data");
  SBSGenericDetector *grinch_adc = new SBSGenericDetector("grinch_adc","GRINCH ADC data");
  grinch_adc->SetModeADC(SBSModeADC::kWaveform);
  grinch_adc->SetModeTDC(SBSModeTDC::kNone);
  grinch_adc->SetStoreEmptyElements(kFALSE);
  grinch_adc->SetStoreRawHits(kFALSE);

  grinch_tdc->SetModeTDC(SBSModeTDC::kTDC);
  //grinch_tdc->SetModeTDC(SBSModeTDC::kCommonStartTDC);
  grinch_tdc->SetModeADC(SBSModeADC::kNone);
  grinch_tdc->SetStoreEmptyElements(kFALSE);
  grinch_tdc->SetStoreRawHits(kFALSE);
  grinch_tdc->SetDisableRefTDC(true);
  bigbite->AddDetector(grinch_adc);
  bigbite->AddDetector(grinch_tdc);

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
  //bigbite->AddDetector( new THaShower("ps", "BigBite preshower") );
  bigbite->AddDetector(hodotdc);
  bigbite->AddDetector(hodoadc);
  //bigbite->AddDetector( new SBSGEMSpectrometerTracker("gem", "GEM tracker") );
  SBSGEMSpectrometerTracker *bbgem = new SBSGEMSpectrometerTracker("gem", "BigBite Hall A GEM data");
  bool pm =  ( pedestalmode != 0 );
  //this will override the database setting:
  bbgem->SetPedestalMode( pm );
  bbgem->SetMakeCommonModePlots( cmplots );
  bigbite->AddDetector(bbgem);
  gHaApps->Add(bigbite);

  SBSEArm *harm = new SBSEArm("sbs","Hadron Arm with HCal");
  SBSHCal* hcal =  new SBSHCal("hcal","HCAL");
  //hcal->SetStoreRawHits(kTRUE);
  hcal->SetStoreEmptyElements(kFALSE);
  harm->AddDetector(hcal);

  SBSGenericDetector* sbstrig= new SBSGenericDetector("trig","HCal trigs");
  sbstrig->SetModeADC(SBSModeADC::kWaveform);
  //sbstrig->SetStoreRawHits(kTRUE);
  sbstrig->SetStoreEmptyElements(kFALSE);
  harm->AddDetector( sbstrig );

  SBSGEMSpectrometerTracker *sbsgem = new SBSGEMSpectrometerTracker("gem", "Super BigBite Hall A GEM data");
  sbsgem->SetPedestalMode( pm );
  sbsgem->SetMakeCommonModePlots( cmplots );
  harm->AddDetector(sbsgem);

  gHaApps->Add(harm);

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

  //gHaEvtHandlers->Add (new THaScalerEvtHandler("Left","HA scaler event type 140"));
  //gHaEvtHandlers->Add (new THaScalerEvtHandler("SBS","HA scaler event type 141"));

  //bigbite->SetDebug(2);
  //harm->SetDebug(2);

  LHRSScalerEvtHandler *lScaler = new LHRSScalerEvtHandler("Left","HA scaler event type 140");
  // lScaler->SetDebugFile(&debugFile);
  gHaEvtHandlers->Add(lScaler);

  SBSScalerEvtHandler *sbsScaler = new SBSScalerEvtHandler("sbs","SBS Scaler Bank event type 1");
  //sbsScaler->AddEvtType(1);             // Repeat for each event type with scaler banks
  sbsScaler->SetUseFirstEvent(kTRUE);
  gHaEvtHandlers->Add(sbsScaler);

  //THaInterface::SetDecoder( SBSSimDecoder::Class() );

  TString prefix = gSystem->Getenv("DATA_DIR");

  //Here we need to revert back to the old behavior based on THaRun for the time being:

  vector<TString> pathlist;
  if( !prefix.IsNull() )
    pathlist.push_back( prefix.Data() );

  if( prefix != "/adaqeb1/data1" )
    pathlist.push_back( "/adaqeb1/data1" );

  if( prefix != "/adaqeb2/data1" )
    pathlist.push_back( "/adaqeb2/data1" );
  
  if( prefix != "/adaqeb3/data1" )
    pathlist.push_back( "/adaqeb3/data1" );

  if( prefix != "/adaq1/data1/sbs" )
    pathlist.push_back( "/adaq1/data1/sbs" );

  //This may take forever to scan
  if( prefix != "/cache/mss/halla/sbs/raw" )
    pathlist.push_back( "/cache/mss/halla/sbs/raw" );

  for( const auto& path: pathlist ) {
    cout << "search paths = " << path << endl;
  }

  int nstream=maxstream+1;

  int stream = 0, segcounter=0;
 
  int segment = firstsegment; 
  int lastsegment = firstsegment;

  TClonesArray *filelist = new TClonesArray("THaRun",10);

  while( segcounter < maxsegments*nstream && segment-firstsegment < maxsegments && stream < nstream ){

    TString codafilename;
    //codafilename.Form("%s_%u.evio.[0-9]+.[0-9]+", fname_prefix, runnum);
  
    codafilename.Form( "%s_%d.evio.%d.%d", fname_prefix, runnum, stream, segment );
    // TString ftest(fname_prefix);
    // if( ftest == "bbgem" || ftest == "e1209019_trigtest" ) {
    //   codafilename.Form("%s_%u.evio.[0-9]+", fname_prefix, runnum);
    // } these lines obsolete:
    
    bool segmentexists=false;

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
      new( (*filelist)[segcounter] ) THaRun( pathlist, codafilename.Data(), "GEN run" );
      cout << "Added segment " << segment << ", CODA file name = " << codafilename << endl;
      segcounter++;
    }

    if( stream < maxstream ){ 
      stream++; 
    } else if( segment - firstsegment < maxsegments ){
      stream = 0;
      segment++;
    }
  }

  cout << "n segments to analyze = " << segcounter << endl;
  //auto* run = new Podd::MultiFileRun(pathlist, codafilename.Data());

  //run->SetFirstSegment(firstsegment);  // Starting segment number
  //run->SetMaxSegments(maxsegments);    // Number of segments
  //run->SetMaxStreams(maxstream+1);  // Number of streams

  // Count physics events sequentially, always starting at 1. This count will be compared to the run's event range.
  // The run's event range can then be specified independently of absolute event numbers.
  analyzer->SetCountMode(THaAnalyzer::kCountPhysics);

  // if( nevents > 0 ) run->SetLastEvent(firstevent+nevents-1);
  // run->SetFirstEvent( firstevent );
  // run->SetDataRequired(THaRunBase::kDate|THaRunBase::kRunNumber);
  // Int_t st = run->Init();
  // if( st != THaRunBase::READ_OK ) {
  //   cerr << "========= Error initializing run" << endl;
  //   return;
  // }

  //  cout << "n segments to analyze = " << run->GetNFiles() << endl;

  prefix = gSystem->Getenv("OUT_DIR");

  TString outfilename;
  //  UInt_t lastsegment = run->GetLastSegment();
  if( nevents > 0 ){
    outfilename.Form( "%s/%s_replayed_%u_stream%d_%d_seg%u_%u_firstevent%ld_nevent%ld.root", prefix.Data(), fname_prefix, runnum,
		      0, maxstream, firstsegment, lastsegment, firstevent, nevents );
  } else {
    outfilename.Form( "%s/%s_fullreplay_%u_stream%d_%d_seg%u_%u.root", prefix.Data(), fname_prefix, runnum,
		      0, maxstream, firstsegment, lastsegment );
  }

  analyzer->SetVerbosity(2);
  analyzer->SetMarkInterval(100);

  analyzer->EnableBenchmarks();

  // Define the analysis parameters
  analyzer->SetOutFile( outfilename.Data() );

  // File to record cuts accounting information
  prefix = gSystem->Getenv("LOG_DIR");
  analyzer->SetSummaryFile(Form("%s/replay_gen_%d_stream%d_%d_seg%d_%d.log", prefix.Data(), runnum,
				0, maxstream, firstsegment, lastsegment));

  prefix = gSystem->Getenv("SBS_REPLAY");
  prefix += "/replay/";

  TString odef_filename = "replay_gen.odef";

  odef_filename.Prepend( prefix );

  analyzer->SetOdefFile( odef_filename );

  //added cut list in order to have
  //TString cdef_filename = "replay_gen_farm.cdef";
  TString cdef_filename = "replay_gen.cdef";

  //the above cdef file includes almost no significant cuts (basically just any hit in BigBite shower or preshower
  cdef_filename.Prepend( prefix );

  analyzer->SetCutFile( cdef_filename );
  //analyzer->SetCompressionLevel(0); // turn off compression

  filelist->Compress();

  //reverting to the behavior of the old GMN script for the time being:
  for( int iseg=0; iseg<filelist->GetEntries(); iseg++ ){
    THaRun *run = ( (THaRun*) (*filelist)[iseg] );
    if( nevents > 0 ) run->SetLastEvent( firstevent + nevents - 1 ); 
    run->SetFirstEvent( firstevent );
    
    run->SetDataRequired(THaRunBase::kDate|THaRunBase::kRunNumber);

    run->Init();
   
    Int_t st = analyzer->Init(run);
    if( st != 0 ) {
      cerr << "========= Error initializing analyzer" << endl;
      return;
    }
    analyzer->Process(run);     // start the actual analysis
  }
}


// int main(int argc, char *argv[])
// {
//   //does this int main actually get called when we run our script?
//   new THaInterface( "The Hall A analyzer", &argc, argv, 0, 0, 1 );
//   uint runnum = 220;
//   uint firstsegment = 0;
//   uint maxsegments = 0;
//   string name_prefix = "e1209019";
//   long firstevent=0;
//   long nevents=-1;
//   int pedestalmode=0;
//   uint nev = -1;
//   if(argc<2 || argc>8){
//     cout << "Usage: replay_gmn runnum(uint) nevents(ulong) firstevent(ulong)"
// 	 << endl
// 	 << "                  name_prefix(string) firstsegment(uint) maxsegments(uint) "
// 	 << endl
// 	 << "                  pedestalmode(int)"
// 	 << endl;
//     return -1;
//   }
//   if(argc>=2) runnum = atoi(argv[1]);
//   if(argc>=3) firstsegment = atoi(argv[2]);
//   if(argc>=4) maxsegments = atoi(argv[3]);
//   if(argc>=5) name_prefix = argv[4];
//   if(argc>=6) firstevent = atoi(argv[5]);
//   if(argc>=7) nevents = atoi(argv[6]);
//   if(argc>=8) pedestalmode = atoi(argv[7]);

//   replay_gmn(runnum, nevents, firstevent,
// 	     name_prefix.c_str(), firstsegment, maxsegments,
// 	     pedestalmode);
//   return 0;
// }
