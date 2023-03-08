//SSeeds 7.13.22 - Short script to extract the boundary crossing (BC) time from the sensitive detector (SD) track information from g4sbs simulations

#include <ctime>
#include <iostream>
#include <vector>
#include "TStopwatch.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TROOT.h"
#include "TObjArray.h"
#include "TChainElement.h"
#include "TClonesArray.h"
#include "TLegend.h"
#include "TMath.h"
#include "gmn_tree.C"
#include <TROOT.h>

//#include "G4SBSRunData.hh"
//#include "/work/halla/sbs/seeds/G4SBS/g4sbs_install/root_macros/gmn_tree.C"

//Static Detector Parameters
const int maxTracks = 1000; // Reasonable limit on tracks to be stored per event
const int maxTdcChan = 10; // Set to accomodate original 5 TDCTrig channels with buffer
const double hcalheight = -0.2897; // Height of HCal above beamline

//Constants
const double PI = TMath::Pi();
const double M_e = 0.00051; //Mass electron (GeV)
const double M_p = 0.938272; //Mass proton (GeV)
const double M_n = 0.939565; //Mass neutron (GeV)

//Geant Parameters
//PDG PID, 2112=Neutron, 2212=Proton
//pdg.lbl.gov/2020/reviews/rpp2020-rev-monte-carlo-numbering.pdf
const int prot = 2212;
const int neut = 2112;

// Get today's date
string getDate(){
  time_t now = time(0);
  tm ltm = *localtime(&now);
  
  string yyyy = to_string(1900 + ltm.tm_year);
  string mm = to_string(1 + ltm.tm_mon);
  string dd = to_string(ltm.tm_mday);
  string date = mm + '_' + dd + '_' + yyyy;
  
  return date;
}

void simTOF(){
  
  // Define a clock to check overall time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );
  
  // Get the date
  string date = getDate();

  //Set processing params
  //TTree *T = 0;
  TChain *C = new TChain("T");
  TFile *fout = new TFile( "simTOFout.root", "RECREATE" );
  int nevent=1;

  //Load in data
  //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/w/halla-scshelf2102/sbs/seeds/G4SBS/g4sbs_install/run_g4sbs_here/gmn_3.5GeV2_elastic_TOF.root");
  //C->Add("/w/halla-scshelf2102/sbs/seeds/G4SBS/g4sbs_install/run_g4sbs_here/gmn_3.5GeV2_elastic_TOF.root");
  //C->Add("/w/halla-scshelf2102/sbs/seeds/G4SBS/g4sbs_install/run_g4sbs_here/gmn_13.3GeV2.root");
  C->Add("/lustre19/expphy/volatile/halla/sbs/pdbforce/g4sbs_output/sbs4_elas_sbs30p/gmn_SBS4_elastic.root");
  //C->Add("gmn_3.5GeV2_elastic_TOF.root");

  gmn_tree *T = new gmn_tree("C");

  /*
  //Set tree variables
  Int_t Harm_HCalScint_hit_nhits;
  vector<int> Harm_HCalScint_hit_sdtridx(100,0);
  //Harm_HCalScint_hit_sdtridx=0;
  vector<int> Harm_HCalScint_hit_cell(100,0);
  //int Harm_HCalScint_hit_cell[100]={0};
  //int Harm_HCalScint_hit_cell;
  //Harm_HCalScint_hit_cell=0;
  vector<double> SDTrack_T(100,0.);
  //SDTrack_T=0;
  vector<int> SDTrack_PID(100,0);
  //SDTrack_PID=0;

  //cout << "cell[0]: " << Harm_HCalScint_hit_cell[0] << endl;


  //Set up the output tree in fout
  //Track number hits
  //C->SetBranchStatus( "Harm.HCalScint.hit.nhits", 1 );
  TBranch *b_Harm_HCalScint_hit_nhits;
  C->SetBranchAddress( "Harm.HCalScint.hit.nhits", &Harm_HCalScint_hit_nhits, &b_Harm_HCalScint_hit_nhits );
  //Track ID
  //C->SetBranchStatus( "Harm.HCalScint.hit.sdtridx", 1 );
  TBranch *b_Harm_HCalScint_hit_sdtridx;
  C->SetBranchAddress( "Harm.HCalScint.hit.sdtridx", &Harm_HCalScint_hit_sdtridx, &b_Harm_HCalScint_hit_sdtridx );
  //Track cell
  //C->SetBranchStatus( "Harm.HCalScint.hit.cell", 1 );
  TBranch *b_Harm_HCalScint_hit_cell;
  C->SetBranchAddress( "Harm.HCalScint.hit.cell", &Harm_HCalScint_hit_cell, &b_Harm_HCalScint_hit_cell );
  //Track PID
  //C->SetBranchStatus( "SDTrack.PID", 1 );
  TBranch *b_SDTrack_PID;
  C->SetBranchAddress( "SDTrack.PID", &SDTrack_PID, &b_SDTrack_PID );
  //Track TOF (Sensitive detector crossing time)-(Scattering time)
  //C->SetBranchStatus( "SDTrack.T", 1 );
  TBranch *b_SDTrack_T;
  C->SetBranchAddress( "SDTrack.T", &SDTrack_T, &b_SDTrack_T );
  */

  //f->GetObject("T",T);

  if( T->GetEntries()!=0 ){
    cout << "Opened file successfully." << endl;
  }else{
    cout << "Error: No file found." << endl;
    return;
  }
  
  
  //Define histograms
  TH1D *TOF_all = new TH1D("TOF_all","Time of Flight over All Modules", 100, 0, 50);
  TH1D *TOF_all_p = new TH1D("TOF_all_p","Time of Flight over All Modules: Protons", 100, 0, 50);
  TH1D *TOF_all_n = new TH1D("TOF_all_n","Time of Flight over All Modules: Neutrons", 100, 0, 50);
  TH2D *TOF_vs_ID = new TH2D("TOF_vs_ID","Time of Flight vs Channel; channel; ns",288,0,288,100,0,50);
  TH2D *TOF_vs_ID_p = new TH2D("TOF_vs_ID_p","Time of Flight vs Channel: Protons; channel; ns",288,0,288,100,0,50);
  TH2D *TOF_vs_ID_n = new TH2D("TOF_vs_ID_n","Time of Flight vs Channel: Neutrons; channel; ns",288,0,288,100,0,50);
  

  //Define loop parameters
  Long64_t Nevents = T->GetEntries();

  cout << "Opened tree with " << Nevents << " simulated events." << endl;

  while( T->GetEntry( nevent++ ) ){ 
    
    if( nevent%1000 == 0 ) cout << "Processing event: " << nevent << "/" << Nevents << endl;


    /*
    //Generate TOF from SD track information
    for( int ihit=0; ihit<Harm_HCalScint_hit_nhits; ihit++ ){
    double_t BCtime = ( *(SDTrack_T) )[( *(Harm_HCalScint_hit_sdtridx) )[ ihit ]];
    int_t PID = ( *(SDTrack_PID) )[( *(Harm_HCalScint_hit_sdtridx) )[ ihit ]];
    int_t cell = Harm_HCalScint_hit_cell;
    TOF_all->Fill( BCtime );
    TOF_vs_ID->Fill( cell, BCtime );
    if( SDTrack == prot ){
    TOF_all_p->Fill( BCtime );
    TOF_vs_ID_p->Fill( cell, BCtime );
    }else{
    TOF_all_n->Fill( BCtime );
    TOF_vs_ID_n->Fill( cell, BCtime );
    }
    }
    */
    /*
      for( int ihit=0; ihit<Harm_HCalScint_hit_nhits; ihit++ ){
      double BCtime = SDTrack_T[ Harm_HCalScint_hit_sdtridx[ ihit ] ];
      int PID = SDTrack_PID[ Harm_HCalScint_hit_sdtridx[ ihit ] ];
      TOF_all->Fill( BCtime );
      //TOF_vs_ID->Fill( Harm_HCalScint_hit_cell, BCtime );
      if( SDTrack_PID[ ihit ] == prot ){
      TOF_all_p->Fill( BCtime );
      //TOF_vs_ID_p->Fill( Harm_HCalScint_hit_cell, BCtime );
      }else if( SDTrack_PID[ ihit ] == neut ){
      TOF_all_n->Fill( BCtime );
      //TOF_vs_ID_n->Fill( Harm_HCalScint_hit_cell, BCtime );
      }
      }
    */
    /*
    for( int ihit=0; ihit<C->Harm_HCalScint_hit_nhits; ihit++ ){
      double_t BCtime = ( *(C->SDTrack_T) )[( *(C->Harm_HCalScint_hit_sdtridx) )[ ihit ]];
      cout << BCtime << endl;
    }
    */

    //cout << "nhits: " << Harm_HCalScint_hit_nhits << endl;
    //cout << "cell.size(): " << Harm_HCalScint_hit_cell.size() << endl;
    cout << "cell[0]: " << Harm_HCalScint_hit_cell[0] << endl;
    //cout << "cell[1]: " << Harm_HCalScint_hit_cell[1] << endl;

  }
  fout->Write();

  st->Stop();

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

}
