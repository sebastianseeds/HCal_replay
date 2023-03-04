
#include <cstdlib>
#include <iostream>
#include <string> 

#include "TStyle.h"
#include "TString.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TSystem.h"

#include "THaRun.h"
#include "THaEvent.h"
#include "THaAnalyzer.h"
#include "THaApparatus.h"
#include "THaVarList.h"
#include "THaInterface.h"
#include "THaHRS.h"
#include "THaDecData.h"

#include "SBSRasteredBeam.h"
#include "SBSScalerEvtHandler.h"
#include "LHRSScalerEvtHandler.h"
#include "SBSScalerHelicity.h"
#include "SBSScalerHelicityReader.h"

int MAX_SEG = 700; 

const std::string SCRIPT = "[replay_beam]: "; 

void replay_beam(int runNum,Long_t firstevent=0,Long_t nevents=-1,int maxsegments=1){

   std::cout << "=============================" << std::endl;
   std::cout << SCRIPT << "Input parameters: "  << std::endl;
   std::cout << "Run:          " << runNum      << std::endl;
   std::cout << "First event:  " << firstevent  << std::endl;
   std::cout << "N events:     " << nevents     << std::endl;
   std::cout << "max segments: " << maxsegments << std::endl;
   std::cout << "=============================" << std::endl;

   // Set up the analyzer
   THaAnalyzer* analyzer = THaAnalyzer::GetInstance();
   if(analyzer){
      analyzer->Close();
   }else{
      analyzer = new THaAnalyzer;
   }
 
   THaEvent* event = new THaEvent;

   // add the LHRS (where the beam signals are)  
   THaHRS* HRSL = new THaHRS("L","Left arm HRS");
   HRSL->AutoStandardDetectors(kFALSE);
   gHaApps->Add( HRSL );

   // add decoder
   THaApparatus* decL = new THaDecData("DL","Misc. Decoder Data");
   gHaApps->Add( decL );

   SBSScalerHelicity *scalerHelicity = new SBSScalerHelicity("scalHel","Scaler Helicity Info"); 
   // scalerHelicity->SetVerbosity(2); 
  
   // add *rastered* beam (LHRS)
   THaApparatus* Lrb = new SBSRasteredBeam("Lrb","Raster Beamline for FADC");
   // add scaler helicity 
   Lrb->AddDetector(scalerHelicity);
   gHaApps->Add(Lrb);

   // add *rastered* beam (SBS)
   THaApparatus* sbs = new SBSRasteredBeam("SBSrb","Raster Beamline for FADC");
   gHaApps->Add(sbs);

   // std::ofstream debugFile; 
   // debugFile.open("lhrs-scaler-dump.txt");

   // LHRS scaler data (ROC10)  
   LHRSScalerEvtHandler *lScaler = new LHRSScalerEvtHandler("Left","HA scaler event type 140");
   // lScaler->SetDebugFile(&debugFile);
   gHaEvtHandlers->Add(lScaler);

   // SBS scaler data (sbsvme29) 
   SBSScalerEvtHandler *sbsScaler = new SBSScalerEvtHandler("sbs","SBS Scaler Bank event type 1");
   // sbsScaler->AddEvType(1);             // Repeat for each event type with scaler banks
   sbsScaler->SetUseFirstEvent(kTRUE);
   gHaEvtHandlers->Add(sbsScaler);

   analyzer->SetEvent(event);
   analyzer->SetCompressionLevel(1);

   // File to record cuts accounting information
   analyzer->SetSummaryFile("sbs_beam_test.log"); // optional

   analyzer->SetVerbosity(10);   

   // find all files associated with the run 
  
   TString prefix = gSystem->Getenv("DATA_DIR");

   char fname_prefix[200]; 
   sprintf(fname_prefix,"e1209019"); 

   int stream=0;

   TClonesArray *filelist = new TClonesArray("THaRun",10);

   // search paths to find the raw data
   std::vector<TString> pathlist;
   pathlist.push_back( prefix );

   if(prefix!="/adaqeb1/data1" )  pathlist.push_back( "/adaqeb1/data1" );
   if(prefix!="/adaq1/data1/sbs") pathlist.push_back( "/adaq1/data1/sbs" );
   if(prefix!="/cache/mss/halla/sbs/raw") pathlist.push_back( "/cache/mss/halla/sbs/raw" );

   for( int i=0; i<pathlist.size(); i++ ){
      std::cout << "search paths = " << pathlist[i] << std::endl;
   }

   if(nevents<0){
      std::cout << SCRIPT << "Requested all events!  Changing maxsegments to " << MAX_SEG << std::endl;
      maxsegments = MAX_SEG;
   }

   int max1 = maxsegments;
   int segcounter=0,segment=0,firstsegment=0,lastsegment=0; 

   bool segmentexists=false;

   //This loop adds all file segments found to the list of THaRuns to process:
   while( segcounter < max1 && segment - firstsegment < maxsegments ){
      std::cout << "segcounter = " << segcounter << ", segment = " << segment << std::endl;
      TString codafilename;
      //codafilename.Form( "%s/bbgem_%d.evio.%d", prefix.Data(), runNum, segment );
      codafilename.Form("%s_%d.evio.%d.%d", fname_prefix, runNum, stream, segment );
      
      TString ftest(fname_prefix);
      if( ftest == "bbgem" || ftest == "e1209019_trigtest" ){
         codafilename.Form("%s_%d.evio.%d", fname_prefix, runNum, segment );
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
         cout << "Added segment " << segment << ", CODA file name = " << codafilename << ", segcounter = " << segcounter << std::endl;
      }
      if( segmentexists ){
         segcounter++;
         lastsegment = segment;
      }
      segment++;
   }

   std::cout << SCRIPT << "n segments to analyze = " << segcounter << endl;
   
   // set up file paths
   prefix = gSystem->Getenv("OUT_DIR");
   
   TString outfilename = Form("%s/gmn_replayed-beam_%d_stream%d_seg%d_%d.root",prefix.Data(),runNum,
		              stream,firstsegment,lastsegment); 

   prefix = gSystem->Getenv("LOG_DIR");
   analyzer->SetSummaryFile(Form("%s/replay_beam_gmn.log", prefix.Data()));

   prefix = gSystem->Getenv("SBS_REPLAY");
   prefix += "/replay/";

   TString odef_filename = Form("%s/replay_beam.odef",prefix.Data());
   TString cdef_filename = Form("%s/replay_beam.cdef",prefix.Data());

   analyzer->SetOdefFile(odef_filename);
   analyzer->SetCutFile(cdef_filename);      
   analyzer->SetOutFile(outfilename);
   std::cout << SCRIPT << "Output file " << outfilename.Data() << " set up " << std::endl; 
  
   filelist->Compress();

   int runSegment=-1;
  
   bool replaySuccess=false;
   std::vector<int> failSegment; 

   for(int iseg=0; iseg<filelist->GetEntries(); iseg++ ){
      THaRun *run = ( (THaRun*) (*filelist)[iseg] );
      if(nevents>0) run->SetLastEvent(nevents); //not sure if this will work as we want it to for multiple file segments chained together
      run->SetFirstEvent( firstevent );
      run->SetDataRequired(THaRunBase::kDate|THaRunBase::kRunNumber);
      // FIXME: to address prescale factor issue (temporary) 
      // TDatime now; 
      // run->SetDate(now); 
      // run->SetDataRequired(0);
      run->Init(); 
      runSegment = run->GetSegment(); 
      if( runSegment>=firstsegment && (runSegment-firstsegment)<maxsegments ){
         analyzer->Process(run);     // start the actual analysis
         std::cout << SCRIPT << "--> Done!" << std::endl;
	 replaySuccess = true;
      }else{
	 std::cout << SCRIPT << "ERROR! Cannot analyze file segment = " << runSegment << std::endl;
	 failSegment.push_back(runSegment);
      }
   }

   if(replaySuccess) std::cout << SCRIPT << "Replay of run " << runNum << " COMPLETE!" << std::endl;

   int NFAIL = failSegment.size();
   if(NFAIL!=0){
      for(int i=0;i<NFAIL;i++){
         std::cout << SCRIPT << "WARNING: Failed replaying file segment " << failSegment[i] << std::endl;
      }
   }

   // clean up  
   delete analyzer;
   gHaVars->Clear();
   gHaPhysics->Delete();
   gHaApps->Delete();
}
