//seeds 6.16.23 - short script to illustrate a method to loop over events from a single run via TChain
//Usage: root 'loop_over_chain.C(<run number>)' 
//Where <run number> is some viable option from /lustre19/expphy/volatile/halla/sbs/seeds/rootfiles/
//
//white/seeds 6.27.23 - some updates which include:
//// - added HCal-specific detector parameters
//// - added additional hcal and bigbite variables for general analysis
//// - added loop over primary cluster blocks for time analysis
//// - added loop over all cluster primary block for time analysis
//// - added fits to 1D histograms after main event loop
//// - default main argument changed (manyRuns = true)
//// - comments and cleanup

#include <sstream>
#include <fstream>
#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <math.h>

#include "TMath.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "TCut.h"
#include "TSystem.h"
#include "TChain.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TVector3.h"

//global variables
const Double_t PI = TMath::Pi();
const Double_t M_e = 0.0051;
const Double_t M_p = 0.938272;
const Double_t M_n = 0.939565;
//HCal-specific detector parameters
const Int_t kNcell = 288; //total number of channels/blocks
const Int_t kNrows = 24; //total number of hcal rows
const Int_t kNcols = 12; //total nubmer of hcal columns
const Int_t maxclus = 30; //reasonable upper limit on total number of clusters in hcal on any given event
const Int_t maxblk = 20; //reasonable upper limit on total number of blocks in primary cluster in hcal on any given event
const Int_t maxtracks = 15; //reasonable upper limit on surviving e' tracks in bigbite
const Double_t hcalheight = -0.2897; //hcal vertical offset correction given common hcal database file db_sbs.hcal.dat, 6.27.23
const Double_t atime_approx_FWHM = 10; //HCal adctime approximate full-width-half-max of elastic peak

//Arguments: <runnum> = runnumber (only used if manyRuns = false); <manyRuns> = if enabled, uses SBS 4 hydrogen at 50% data, all runs (11589,11590,11592)
void loop_over_chain_v2( Int_t runnum = 13747, bool manyRuns = true ){//main
   
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
    //alternatively, adds the expert HCal root files by run number (runnum). note the syntax for Form(), runnum is inserted at %d. The C->Add() command adds all files which match the criteria provided, so the wildcard "*" gets all file segments which exist at this path for runnum
    C->Add(Form("/lustre19/expphy/volatile/halla/sbs/seeds/rootfiles/hcal_gmn_fullreplay_%d*",runnum));

  }

  //declare hcal variables. All of these have a single value per event corresponding to the primary cluster primary block (HCALtdc, HCALatime), or just primary cluster (HCALx, HCALy, HCALe).
  Double_t HCALx, HCALy, HCALe, HCALtdc, HCALatime, HCALid, HCALnblk, HCALnclus; //declare single valued variables for single valued branches
  Double_t HCALid_clusblk[maxblk], HCALe_clusblk[maxblk], HCALatime_clusblk[maxblk], HCALtdc_clusblk[maxblk]; //declare arrays of size = hcal total blocks for clus_blk branches
  Double_t HCALid_clus[maxclus], HCALe_clus[maxclus], HCALatime_clus[maxclus], HCALtdc_clus[maxclus], HCALnblk_clus[maxclus]; //declare arrays size = hcal total clusters for clus branches

  // declare BigBite variables
  Double_t BBps_e, BBps_x, BBps_y, BBsh_e, BBsh_x, BBsh_y, BBtr_p, BBps_atime;

  //declare kinematic variables from e' track
  Double_t ekineW2;
  
  //By default, switches OFF (second argument) all branches (first argument) to speed up processing
  C ->SetBranchStatus( "*", 0 );

  ////Only switch ON the branches to look at
  //HCal single valued branches ('primary' is highest energy)
  C->SetBranchStatus( "sbs.hcal.x", 1 ); //location (m) of center (energy weighted centroid) of primary cluster in dispersive (vertical) direction hcal
  C->SetBranchStatus( "sbs.hcal.y", 1 ); //location (m) of center (energy weighted centroid) of primary cluster in transverse (horizontal) direction hcal
  C->SetBranchStatus( "sbs.hcal.e", 1 ); //total energy (GeV) primary cluster
  C->SetBranchStatus( "sbs.hcal.tdctimeblk", 1 ); //tdc time value (ns) primary block in primary cluster
  C->SetBranchStatus( "sbs.hcal.atimeblk", 1 ); //adc time value (ns) primary block in primary cluster
  C->SetBranchStatus( "sbs.hcal.idblk", 1 ); //Channel number (1-288) of primary block in primary cluster
  //HCal primary cluster block branches (sbs.hcal.nblk defines array size per event, each element (block) sorted by highest to lowest energy)
  C->SetBranchStatus( "sbs.hcal.nblk", 1 ); //Total number of blocks in primary cluster (only single valued branch here)
  C->SetBranchStatus( "sbs.hcal.clus_blk.id", 1 ); //id of nth block element in primary cluster (1-288)
  C->SetBranchStatus( "sbs.hcal.clus_blk.e", 1 ); //energy of nth block element in primary cluster (GeV)
  C->SetBranchStatus( "sbs.hcal.clus_blk.atime", 1 ); //adc time of nth block element in primary cluster (ns)
  C->SetBranchStatus( "sbs.hcal.clus_blk.tdctime", 1 ); //tdc time of nth block element in primary cluster (ns)
  //HCal all cluster primary block branches (sbs.hcal.nclus defines array size per event, each element (cluster) sorted by highest to lowest energy)
  C->SetBranchStatus( "sbs.hcal.nclus", 1 ); //Total number of blocks in primary cluster (only single valued branch here)
  C->SetBranchStatus( "sbs.hcal.clus.id", 1 ); //id of primary block element in nth cluster (1-288)
  C->SetBranchStatus( "sbs.hcal.clus.e", 1 ); //energy of primary block element in nth cluster (GeV)
  C->SetBranchStatus( "sbs.hcal.clus.atime", 1 ); //adc time of primary block element in nth cluster (ns)
  C->SetBranchStatus( "sbs.hcal.clus.tdctime", 1 ); //tdc time of primary block element in nth cluster (ns)
  C->SetBranchStatus( "sbs.hcal.clus.nblk", 1 ); //number of blocks in nth cluster
  //BigBite branches
  C->SetBranchStatus( "bb.ps.atimeblk", 1 ); //adc time of primary block primary cluster in BigBite Calorimeter (BBCal) in e-arm
  C->SetBranchStatus( "bb.tr.p", 1 ); //best track momentum in e-arm
  C->SetBranchStatus( "bb.ps.e", 1 ); //BBCal preshower primary cluster total energy (GeV)
  C->SetBranchStatus( "bb.ps.x", 1 ); //BBCal preshower location (m) of center (energy weighted centroid) of primary cluster in dispersive (vertical) direction
  C->SetBranchStatus( "bb.ps.y", 1 ); //BBCal preshower location (m) of center (energy weighted centroid) of primary cluster in transverse (horizontal) direction
  C->SetBranchStatus( "bb.sh.e", 1 ); //BBCal shower primary cluster total energy (GeV)
  C->SetBranchStatus( "bb.sh.x", 1 ); //BBCal shower location (m) of center (energy weighted centroid) of primary cluster in dispersive (vertical) direction
  C->SetBranchStatus( "bb.sh.y", 1 ); //BBCal shower location (m) of center (energy weighted centroid) of primary cluster in transverse (horizontal) direction
  //Kinematic variables calculated from e' tracks
  C->SetBranchStatus( "e.kine.W2", 1 ); //Invariant mass (W^2) transfer (GeV). Should be the mass of the scattered nucleon for elastic scattering.

  ////Link the branch addresses to the declared variables. note the & dereference needed for single valued variables
  //HCal single valued branch links
  C->SetBranchAddress( "sbs.hcal.x", &HCALx );
  C->SetBranchAddress( "sbs.hcal.y", &HCALy );
  C->SetBranchAddress( "sbs.hcal.e", &HCALe );
  C->SetBranchAddress( "sbs.hcal.tdctimeblk", &HCALtdc );
  C->SetBranchAddress( "sbs.hcal.atimeblk", &HCALatime );
  C->SetBranchAddress( "sbs.hcal.idblk", &HCALid );
  //HCal primary cluster block branch links
  C->SetBranchAddress( "sbs.hcal.nblk", &HCALnblk );
  C->SetBranchAddress( "sbs.hcal.clus_blk.id", HCALid_clusblk ); //remove dereference & for array with size of cluster per event
  C->SetBranchAddress( "sbs.hcal.clus_blk.e", HCALe_clusblk ); //remove dereference & for array with size of cluster per event
  C->SetBranchAddress( "sbs.hcal.clus_blk.atime", HCALatime_clusblk ); //remove dereference & for array with size of cluster per event
  C->SetBranchAddress( "sbs.hcal.clus_blk.tdctime", HCALtdc_clusblk ); //remove dereference & for array with size of cluster per event
  //HCal all cluster primary block branches (sbs.hcal.nclus defines array size per event, each element (cluster) sorted by highest to lowest energy)
  C->SetBranchAddress( "sbs.hcal.nclus", &HCALnclus );
  C->SetBranchAddress( "sbs.hcal.clus.id", HCALid_clus ); //remove dereference & for array with size of cluster per event
  C->SetBranchAddress( "sbs.hcal.clus.e", HCALe_clus ); //remove dereference & for array with size of cluster per event
  C->SetBranchAddress( "sbs.hcal.clus.atime", HCALatime_clus ); //remove dereference & for array with size of cluster per event
  C->SetBranchAddress( "sbs.hcal.clus.tdctime", HCALtdc_clus ); //remove dereference & for array with size of cluster per event
  C->SetBranchAddress( "sbs.hcal.clus.nblk", HCALnblk_clus ); //remove dereference & for array with size of cluster per event
  //Bigbite branch links
  C->SetBranchAddress( "bb.ps.atimeblk", &BBps_atime );
  C->SetBranchAddress( "bb.tr.p", &BBtr_p );
  C->SetBranchAddress( "bb.ps.e", &BBps_e );
  C->SetBranchAddress( "bb.ps.x", &BBps_x );
  C->SetBranchAddress( "bb.ps.y", &BBps_y );
  C->SetBranchAddress( "bb.sh.e", &BBsh_e );
  C->SetBranchAddress( "bb.sh.e", &BBsh_e );
  C->SetBranchAddress( "bb.sh.e", &BBsh_e );

  C->SetBranchAddress( "e.kine.W2", &ekineW2 );

  //declare a simple output root file for analysis. Toggle the output file name depending on whether the analysis is over a single run or many runs.
  std::string outputfilename;
  if( manyRuns ){
    outputfilename = "manyRuns_out.root";
  }else{
    outputfilename = Form("run%d_out.root",runnum);
  }
  TFile *fout = new TFile(outputfilename.c_str(),"RECREATE");

  ////declare 1D histograms to fill per event
  //No cut histograms to compare with histograms with cuts
  TH1D *hTDC = new TH1D("hTDC","HCal TDC, all channels; ns", 800, -400, 400); //RAW, no cuts
  TH1D *hHCALe = new TH1D("hHCALe","HCal Energy, all channels; GeV", 100, 0, 0.2); //RAW, no cuts
  TH1D *hADCtime = new TH1D("hADCtime","HCal ADC time, all channels; ns", 150, 1, 150); //RAW, no cuts
  TH1D *hBBps_e =  new TH1D("hBBps_e","BigBite preshower energy;GeV", 100, 0, 2); //RAW, no cuts
  //Elastic cut histograms
  TH1D *hTDC_ecut = new TH1D("hTDC_ecut","HCal TDC elastic cut, all channels; ns", 800, -400, 400); //RAW, no cuts
  TH1D *hHCALe_ecut = new TH1D("hHCALe_ecut","HCal Energy elastic cut, all channels; GeV", 100, 0, 0.2); //RAW, no cuts
  TH1D *hADCtime_ecut = new TH1D("hADCtime_ecut","HCal ADC time elastic cut, all channels; ns", 150, 1, 150); //RAW, no cuts
  //kinematic variable histograms
  TH1D *hekineW2 =  new TH1D("hekineW2","W2;GeV^2", 100, 0, 3);

  ////declare a 2D histogram to check ADC time vs HCal channel
  //No cut histograms
  TH2D *hADCtime_ID = new TH2D("hADCtime_ID","HCal ADC time vs channel; channel id; ns", kNcell, 0, kNcell, 150, 0, 150);

  //Elastic cut histograms
  TH2D *hADCtime_ID_ecut = new TH2D("hADCtime_ID_ecut","HCal ADC time vs channel elastic cut; channel id; ns", kNcell, 0, kNcell, 150, 0, 150);
  TH2D *hADCtime_clus_ecut = new TH2D("hADCtime_clus_ecut","HCal ADC time vs cluster, primary block, elastic cut; cluster id; ns", maxclus, 0, maxclus, 150, 0, 150);
  TH2D *hADCtime_clusblk_ecut = new TH2D("hADCtime_clusblk_ecut","HCal ADC time vs primary cluster block, elastic cut; cluster id; ns", maxblk, 0, maxblk, 150, 0, 150);

  //make a long int just in case the number of events requires more bits to store than a simple int can do
  Long64_t nevents = C->GetEntries();
  
  //send the total number of entries after adding to the chain to terminal
  cout << endl << "Total entries from files: " << nevents << endl << endl;

  //for loop over all events
  for( Long64_t nevent = 1; nevent < nevents; nevent++ ){

    //Write out to console on each event to keep track of progress. Added .flush() with "\r" to ease viewing
    cout << " Entry = " << nevent << " / " << nevents << "\r";
    cout.flush();

    //Get the event from the chain
    C->GetEntry(nevent); 

    /////////////////////////////////////
    //BEGIN CUTS AND FILLING HISTOGRAMS//
    /////////////////////////////////////

    //fill no cut histograms
    hHCALe->Fill(HCALe); //hcal primary cluster energy
    hBBps_e->Fill(BBps_e); //bbcal preshower primary cluster energy
    hTDC->Fill(HCALtdc); //hcal primary cluster primary block tdc time
    hADCtime->Fill(HCALatime); //hcal primary cluster primary block adc time
    hADCtime_ID->Fill(HCALid,HCALatime); //hcal primary cluster primary block adc time vs hcal channel
    hekineW2->Fill(ekineW2); //invariant mass from e'

    ///CUTS
    //make cuts to get better distribution now, want hcal.e > 10gev
    if( HCALe < 0.02 ){
      continue; //Do not process anything else for this event if the energy of the primary cluster in hcal is < 10.5 MeV
    }
    //cut on preshower energy for particle id (cut out pions from electrons in e-arm)
    if( BBps_e < 0.2 ){
      continue;
    }
    //cut on invariant mass to isolate elastic peak
    if( ekineW2 > 1.2 ){
      continue;
    }
    //cut with HCALid attempt(does nothing to graph, if reverse to <, deletes data):
    //    if (HCALid > 3000) {
    //  continue;
    //   } 

    //Fill the histograms with variables that survive the above cuts. These histograms have much better signal to noise.
    hHCALe_ecut->Fill(HCALe);
    hTDC_ecut->Fill(HCALtdc);
    hADCtime_ecut->Fill(HCALatime);
    hADCtime_ID_ecut->Fill(HCALid,HCALatime);

    //Look at each block in the primary custer that survives the elastic cuts by looping over all of its blocks
    for( Int_t b=0; b<HCALnblk; b++ ){
      Double_t clusblk_id = HCALid_clusblk[b]; //gives the hcal block id of the cluster block element (not yet used)
      Double_t clusblk_idx = (Double_t)b; //gives the block index in order of decending energy for this event. Cast the Int_t into a Double_t for clarity.
      Double_t clusblk_atime = HCALatime_clusblk[b]; //gives the cluster block adc time

      hADCtime_clusblk_ecut->Fill( clusblk_idx, clusblk_atime );
    }

    //Look at each cluster that survives the elastic cuts by looping over all of its clusters
    for( Int_t c=0; c<HCALnclus; c++ ){
      Double_t cluster_id = HCALid_clus[c]; //gives the hcal block id of the cluster element (not yet used)
      Double_t cluster_idx = (Double_t)c; //gives the cluster index in order of decending energy for this event. Cast the Int_t into a Double_t for clarity.
      Double_t cluster_atime = HCALatime_clus[c]; //gives the cluster adc time

      hADCtime_clus_ecut->Fill( cluster_idx, cluster_atime );
    }
    
  }//endloop over events

  ////////////////////////////////////////////////////
  //FIT ADC TIME HISTOGRAMS (NO CUT AND ELASTIC CUT)//
  ////////////////////////////////////////////////////

  //declare a canvas to hold our fits
  TCanvas *c1 = new TCanvas("atime fits","ADC time with/without elastic cut",1600,1200);
  c1->Divide(1,2); //divide the canvas with 2 rows and 1 column
  gStyle->SetOptStat(0); //remove the stats box in the top right from the histograms on the canvas
  
  //switch to the first sub-canvas
  c1->cd(1);

  //declare a 1D fit function to fit no cut adc time histogram
  TF1 *f1;

  //declare some dynamic fit variables
  Int_t atimeBinMax = hADCtime->GetMaximumBin(); //Get the bin correspoding to the maximum y value in the histogram
  Double_t atimeBinCenter = hADCtime->GetBinCenter( atimeBinMax ); //Get the time in ns associated with that bin
  Double_t atime_fitLowerLim = atimeBinCenter - atime_approx_FWHM; //Define a lower limit on the gaussian fit with the atimeBinCenter and the approximate FWHM of the elastic peak
  Double_t atime_fitUpperLim = atimeBinCenter + atime_approx_FWHM; //Define an upper limit on the gaussian fit with the atimeBinCenter and the approximate FWHM of the elastic peak

  //Do the fit
  hADCtime->Fit( "gaus", "Q", "", atime_fitLowerLim, atime_fitUpperLim ); //fit the no cut adc time with a root predifined gaussian with the dynamic limits defined above
  f1 = hADCtime->GetFunction( "gaus" ); //set the fit function equal to this fit

  //extract fit parameters. gaussian has three parameters, f(x) = p0 * exp(-0.5 * ((x - p1) / p2)^2)
  Double_t atime_mean = f1->GetParameter(1); //where p1 is the mean
  Double_t atime_sig = f1->GetParameter(2); //where p2 is the std dev

  hADCtime->SetTitle(Form("ADC time, no cut. Mean:%f Sigma:%f",atime_mean,atime_sig)); //Write the relevant fit parameters to the title of the histogram.
  hADCtime->Draw();

  //switch to the second sub-canvas
  c1->cd(2);

  //declare a 1D fit function to fit elastic cut adc time histogram
  TF1 *f2;

  //declare some dynamic fit variables same as before
  Int_t atimeBinMax_ecut = hADCtime_ecut->GetMaximumBin();
  Double_t atimeBinCenter_ecut = hADCtime_ecut->GetBinCenter( atimeBinMax_ecut );
  Double_t atime_fitLowerLim_ecut = atimeBinCenter_ecut - atime_approx_FWHM;
  Double_t atime_fitUpperLim_ecut = atimeBinCenter_ecut + atime_approx_FWHM;

  //Do the fit same as before
  hADCtime_ecut->Fit( "gaus", "Q", "", atime_fitLowerLim_ecut, atime_fitUpperLim_ecut );
  f2 = hADCtime_ecut->GetFunction( "gaus" );

  //extract fit parameters same as before
  Double_t atime_mean_ecut = f2->GetParameter(1);
  Double_t atime_sig_ecut = f2->GetParameter(2);

  hADCtime_ecut->SetTitle(Form("ADC time, elastic cut. Mean:%f Sigma:%f",atime_mean_ecut,atime_sig_ecut));
  hADCtime_ecut->Draw();

  //write the canvas to the output file
  c1->Write();
  
  //write the output file to outputfilename (current directory by default)
  fout->Write();

  cout << endl << "Ended analysis without errors" << endl << endl; //final reporting

}// end main







