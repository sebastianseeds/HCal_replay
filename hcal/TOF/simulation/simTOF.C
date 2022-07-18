//SSeeds 7.13.22 - Short script to extract the boundary crossing (BC) time from the sensitive detector (SD) track information from g4sbs simulations

#include "TMatrixD.h"
#include "TVectorD.h"
#include "TFile.h"
#include "TChain.h"
#include <iostream>
#include <fstream>
#include <vector>
#include "TString.h"
#include "TDecompSVD.h"
#include "TCut.h"
#include "TEventList.h"
#include "gmn_tree.C"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TVector3.h"
#include "TMath.h"
#include "TF2.h"
#include "TChainElement.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include "TObjArray.h"
#include "TObjString.h"

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
  TChain *C = new TChain("T");
  TFile *fout = new TFile( "simTOFout.root", "RECREATE" );
  int nevent=1;

  //Load in data
  C->Add("/work/halla/sbs/seeds/G4SBS/g4sbs_install/run_g4sbs_here/gmn_3.5GeV2_elastic_TOF.root");

  gmn_tree *T = new gmn_tree(C);


  if( C->GetEntries()!=0 ){
    cout << "Opened file successfully." << endl;
  }else{
    cout << "Error: No file found." << endl;
    return;
  }
  
  
  //Define TOF histograms
  TH1D *TOF_all = new TH1D("TOF_all","Time of Flight over All Modules", 100, 35, 45);
  TH1D *TOF_all_p = new TH1D("TOF_all_p","Time of Flight over All Modules: Protons", 100, 35, 45);
  TH1D *TOF_all_n = new TH1D("TOF_all_n","Time of Flight over All Modules: Neutrons", 100, 35, 45);
  TH2D *TOF_vs_ID = new TH2D("TOF_vs_ID","Time of Flight vs Channel; channel; ns",288,0,288,100,35,45);
  TH2D *TOF_vs_ID_p = new TH2D("TOF_vs_ID_p","Time of Flight vs Channel: Protons; channel; ns",288,0,288,100,35,45);
  TH2D *TOF_vs_ID_n = new TH2D("TOF_vs_ID_n","Time of Flight vs Channel: Neutrons; channel; ns",288,0,288,100,35,45);
  //Define supplementary histograms
  TH1D *hvX = new TH1D("hvX","Number of Detections vs Transverse X; X (m)",100,-2,3);
  TH1D *pvX = new TH1D("pvX","Number of Protons vs Transverse X; X (m)",100,-2,3);
  TH1D *nvX = new TH1D("nvX","Number of Neutrons vs Transverse X; X (m)",100,-2,3);
  TH1D *hvY = new TH1D("hvY","Number of Detections vs Dispersive Y; Y (m)",100,-8,-4);
  TH1D *pvY = new TH1D("pvY","Number of Protons vs Dispersive Y; Y (m)",100,-8,-4);
  TH1D *nvY = new TH1D("nvY","Number of Neutrons vs Dispersive Y; Y (m)",100,-8,-4);
  
  //Define loop parameters
  Long64_t Nevents = C->GetEntries();

  cout << "Opened tree with " << Nevents << " simulated events." << endl;

  while( T->GetEntry( nevent++ ) ){ 
    
    if( nevent%1000 == 0 ){
      cout << "Processing event: " << nevent << "/" << Nevents << "\r";
      cout.flush();
    }
    
    //Generate TOF from SD track information
    for( int ihit=0; ihit<T->Harm_HCalScint_hit_nhits; ihit++ ){
      Double_t BCtime = ( *(T->SDTrack_T) )[( *(T->Harm_HCalScint_hit_sdtridx) )[ ihit ]];
      Int_t PID = ( *(T->SDTrack_PID) )[( *(T->Harm_HCalScint_hit_sdtridx) )[ ihit ]];
      Double_t hadPosY = ( *(T->SDTrack_posx) )[( *(T->Harm_HCalScint_hit_sdtridx) )[ ihit ]];
      Double_t hadPosX = ( *(T->SDTrack_posy) )[( *(T->Harm_HCalScint_hit_sdtridx) )[ ihit ]];

      //cout << hadPosY << endl;

      Int_t cell = (*(T->Harm_HCalScint_hit_cell))[ihit];
      hvY->Fill( hadPosY );
      hvX->Fill( hadPosX );
      TOF_all->Fill( BCtime );
      TOF_vs_ID->Fill( cell, BCtime );
      if( PID == prot ){
	pvY->Fill( hadPosY );
	pvX->Fill( hadPosX );
	TOF_all_p->Fill( BCtime );
	TOF_vs_ID_p->Fill( cell, BCtime );
      }else if( PID==neut ){
	nvY->Fill( hadPosY );
	nvX->Fill( hadPosX );
	TOF_all_n->Fill( BCtime );
	TOF_vs_ID_n->Fill( cell, BCtime );
      }
    }
  }
  fout->Write();

  st->Stop();

cout << "Analysis complete. Results written to file: simTOFout.root" << endl;

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

}
