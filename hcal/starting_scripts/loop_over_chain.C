//seeds 6.16.23 - short script to illustrate a method to loop over events from a single run via TChain
//Usage: root 'loop_over_chain.C(<run number>)' 
//Where <run number> is some viable option from /lustre19/expphy/volatile/halla/sbs/seeds/rootfiles/

// #include <sstream>
// #include <fstream>
// #include "TMath.h"
// #include "TH1F.h"
// #include <TH2.h>
// #include <TStyle.h>
// #include <TGraph.h>
// #include <TROOT.h>
// #include <TMath.h>
// #include <TLegend.h>
// #include <TPaveLabel.h>
// #include <TProfile.h>
// #include <TPolyLine.h>
// #include <TObjArray.h>
// #include <cmath>
// #include <cstdio>
// #include <cstdlib>
// #include <math.h>
// #include <stack>
// #include "TLorentzVector.h"
// #include "TCut.h"
// #include "TLatex.h"
// using namespace std;

#include <TSystem.h>
#include <TChain.h>
#include <TString.h>
#include "TFile.h"
#include "TTree.h"
#include <TNtuple.h>
#include "TCanvas.h"
#include "TLatex.h"
#include <iostream>

//global variables
const Int_t maxTracks = 16; 
const Double_t PI = TMath::Pi();
const Double_t M_e = 0.0051;
const Double_t M_p = 0.938272;
const Double_t M_n = 0.939565;

const Double_t hcalheight = -0.2897;

//Arguments: <runnum> = runnumber (only used if manyRuns = false); <manyRuns> = if enabled, uses SBS 4 hydrogen at 50% data, all runs (11589,11590,11592)
void loop_over_chain( Int_t runnum = 13747, bool manyRuns = true ){//main
   
  //create a new TChain to add many root files to
  TChain *C = new TChain("T");

  //One can add many runs to the TChain (NOTE: these files from the general SBS replay data don't contain the HCAL samples (waveform) branches)
  if( manyRuns ){
    //Just add a single run, single segment
    //C->Add("/work/halla/sbs/sbs-gmn/pass0/SBS4/LH2/rootfiles/e1209019_fullreplay_11589_stream0_seg10_13.root");
    //Just add a single run
    //C->Add("/work/halla/sbs/sbs-gmn/pass0/SBS4/LH2/rootfiles/e1209019_fullreplay_11589*");
    //Add all SBS-4 50% field data
    C->Add("/work/halla/sbs/sbs-gmn/pass0/SBS4/LH2/rootfiles/e1209019_fullreplay_11589*");
    C->Add("/work/halla/sbs/sbs-gmn/pass0/SBS4/LH2/rootfiles/e1209019_fullreplay_11590*");
    C->Add("/work/halla/sbs/sbs-gmn/pass0/SBS4/LH2/rootfiles/e1209019_fullreplay_11592*");
  }else{
    //alternatively, adds the expert HCal root files by run number (runnum). note the syntax for Form(), runnum is inserted at %d. The C->Add() command adds all files which match the criteria provided, so the wildcard * gets all file segments which exist at this path for runnum
    C->Add(Form("/lustre19/expphy/volatile/halla/sbs/seeds/rootfiles/hcal_gmn_fullreplay_%d*",runnum));
  }


  //declare hcal variables. All of these have a single value per event corresponding to the primary cluster primary block (HCALtdc, HCALatime), or just primary cluster (HCALx, HCALy, HCALe).
  Double_t HCALx, HCALy, HCALe, HCALtdc, HCALatime, HCALid;

  //By default, switches OFF (second argument) all branches (first argument) to speed up processing
  C ->SetBranchStatus("*",0);

  //Switch ON the branches to look at
  C->SetBranchStatus("sbs.hcal.x",1);
  C->SetBranchStatus("sbs.hcal.y",1);
  C->SetBranchStatus("sbs.hcal.e",1);
  C->SetBranchStatus("sbs.hcal.tdctimeblk",1);
  C->SetBranchStatus("sbs.hcal.atimeblk",1);
  C->SetBranchStatus("sbs.hcal.idblk",1);

  //Link the branch addresses to the declared variables. note the & dereference needed for single valued variables
  C->SetBranchAddress("sbs.hcal.x", &HCALx);
  C->SetBranchAddress("sbs.hcal.y", &HCALy);
  C->SetBranchAddress("sbs.hcal.e", &HCALe);
  C->SetBranchAddress("sbs.hcal.tdctimeblk", &HCALtdc);
  C->SetBranchAddress("sbs.hcal.atimeblk", &HCALatime);
  C->SetBranchAddress("sbs.hcal.idblk", &HCALid);

  //declare a simple output root file for analysis. 
  std::string outputfilename = "out.root";
  TFile *fout = new TFile(outputfilename.c_str(),"RECREATE");

  //declare 1D histograms to fill per event
  TH1D *hTDC = new TH1D("hTDC","HCal TDC, all channels; ns", 800, -400, 400);
  TH1D *hADCtime = new TH1D("hADCtime","HCal ADC time, all channels; ns", 150, 1, 150); //start at 1 to avoid underflow

  //declare a 2D histogram to check ADC time vs HCal channel
  TH2D *hADCtime_ID = new TH2D("hADCtime_ID","HCal ADC time vs channel; channel id; ns", 288, 0, 288, 150, 0, 150);

  //make a long int just in case the number of events requires more bits to store than a simple int can do
  Long64_t nevents = C->GetEntries();
  
  //send the total number of entries after adding to the chain to terminal
  cout << "Total entries from files: " << nevents << endl;

  //for loop over all events
  for( Long64_t nevent = 1; nevent <nevents; nevent++ ){

    //Write out to console on each event to keep track of progress
    if( nevent%10000 == 0 )
      cout << " Entry = " << nevent << " / " << nevents << endl;
    
    //Get the event from the chain
    C->GetEntry(nevent); 

    //Fill the histograms with the tdc and adc time values for this event
    hTDC->Fill(HCALtdc);
    hADCtime->Fill(HCALatime);
    hADCtime_ID->Fill(HCALid,HCALatime);
    
  }//endloop over events

  //write the output file to outputfilename (current directory by default)
  fout->Write();

}// end main







