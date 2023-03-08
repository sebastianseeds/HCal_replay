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
const int Nchan = 288;
const int Nrows = 24;
const int Ncols = 12;

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
  //C->Add("/work/halla/sbs/seeds/G4SBS/g4sbs_install/run_g4sbs_here/gmn_3.5GeV2_elastic_TOF.root");
  C->Add("/work/halla/sbs/ilsmythe/G4SBS/install/run_g4sbs_here/gmn_4.5GeV2.root");

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
  TH2D *TOF_vs_row_p = new TH2D("TOF_vs_row_p","Time of Flight vs Row Proton; row; ns",Nrows,0,Nrows,100,35,45);  
  TH2D *TOF_vs_col_p = new TH2D("TOF_vs_col_p","Time of Flight vs Column Proton; col; ns",Ncols,0,Ncols,100,35,45);
  TH2D *TOF_vs_row_n = new TH2D("TOF_vs_row_n","Time of Flight vs Row Neutron; row; ns",Nrows,0,Nrows,100,35,45);  
  TH2D *TOF_vs_col_n = new TH2D("TOF_vs_col_n","Time of Flight vs Column Neutron; col; ns",Ncols,0,Ncols,100,35,45);

  /*
  TH1D *TOFp[Nchan];
  TH1D *TOFn[Nchan];
  for( int i=0; i<Nchan; i++ ){
    TOFp[i] = new TH1D(Form("TOF Proton Channel %d",i),";TOF (ns)", 40, 35, 45);
    TOFn[i] = new TH1D(Form("TOF Neutron Channel %d",i),";TOF (ns)", 40, 35, 45);
  }
  */
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

      //cout << (int)cell/Ncols << " " << cell%Ncols << endl;

      hvY->Fill( hadPosY );
      hvX->Fill( hadPosX );
      TOF_all->Fill( BCtime );
      TOF_vs_ID->Fill( cell, BCtime );
      if( PID == prot ){
	pvY->Fill( hadPosY );
	pvX->Fill( hadPosX );
	TOF_all_p->Fill( BCtime );
	TOF_vs_ID_p->Fill( cell, BCtime );
	TOF_vs_row_p->Fill( (int)cell/Ncols, BCtime );
	TOF_vs_col_p->Fill( cell%Ncols, BCtime );
      }else if( PID==neut ){
	nvY->Fill( hadPosY );
	nvX->Fill( hadPosX );
	TOF_all_n->Fill( BCtime );
	TOF_vs_ID_n->Fill( cell, BCtime );
	TOF_vs_row_n->Fill( (int)cell/Ncols, BCtime );
	TOF_vs_col_n->Fill( cell%Ncols, BCtime );
      }
    }
  }
  
  //Construct graphs
  Double_t posErr[Nchan] = {0.};
  TF1 *f1;
  //For Protons
  Double_t X[Nchan];
  Double_t Xval[Nchan];
  Double_t Xerr[Nchan];
  TH1D *slice_p[Nchan];
  Double_t TOFmin_p = 100;
  Double_t corrErr[Nchan] = {0.};
  for( Int_t x=1; x<=Nchan; x++ ){
    X[x-1] = x-1;
    slice_p[x-1] = TOF_vs_ID_p->ProjectionY(Form("slice_p_%d",x),x,x);
    slice_p[x-1]->Fit("gaus","Q","",35.,45.);
    f1=slice_p[x-1]->GetFunction("gaus");
    if(slice_p[x-1]->GetEntries()>10){
      Xval[x-1] = f1->GetParameter(1);
      Xerr[x-1] = f1->GetParameter(2);
      if( Xval[x-1]>38&&Xval[x-1]<TOFmin_p ) TOFmin_p=Xval[x-1];
    }
  }
  TGraphErrors *cTOF_X = new TGraphErrors( Nchan, X, Xval, posErr, Xerr );
  cTOF_X->GetXaxis()->SetLimits(0,Nchan);  
  cTOF_X->GetYaxis()->SetLimits(35.,45.);
  cTOF_X->SetTitle("Time of Flight Proton");
  cTOF_X->GetXaxis()->SetTitle("Channel");
  cTOF_X->GetYaxis()->SetTitle("Time of Flight (ns)");
  cTOF_X->SetMarkerStyle(20); // idx 20 Circles, idx 21 Boxes
  cTOF_X->Write("cTOF_P");

  cout << TOFmin_p << endl << endl;
  for( Int_t v=0; v<Nchan; v++ ){
    Xval[v]-=TOFmin_p;
    cout << Xval[v] << endl;
  }

  TGraphErrors *cTOF_X_corr = new TGraphErrors( Nchan, X, Xval, posErr, corrErr );
  cTOF_X_corr->GetXaxis()->SetLimits(0,Nchan);  
  cTOF_X_corr->GetYaxis()->SetLimits(35.,45.);
  cTOF_X_corr->SetTitle("TOF Corrections Proton");
  cTOF_X_corr->GetXaxis()->SetTitle("Channel");
  cTOF_X_corr->GetYaxis()->SetTitle("Correction (ns)");
  cTOF_X_corr->SetMarkerStyle(20); // idx 20 Circles, idx 21 Boxes
  cTOF_X_corr->Write("cTOF_P_corr");
  
  //For Neutrons
  Double_t Y[Nchan];
  Double_t Yval[Nchan];
  Double_t Yerr[Nchan];
  Double_t Rval[Nrows];
  Double_t Rerr[Nrows];
  TH1D *slice_n[Nchan];
  Double_t TOFmin_n=100;
  for( Int_t x=1; x<=Nchan; x++ ){
    Y[x-1] = x-1;
    slice_n[x-1] = TOF_vs_ID_n->ProjectionY(Form("slice_n_%d",x),x,x);
    slice_n[x-1]->Fit("gaus","Q","",35.,45.);
    f1=slice_n[x-1]->GetFunction("gaus");
    if(slice_n[x-1]->GetEntries()>10){
      Yval[x-1] = f1->GetParameter(1);
      Yerr[x-1] = f1->GetParameter(2);
      if( Yval[x-1]>5&&Yval[x-1]<TOFmin_n ) TOFmin_n=Yval[x-1];
    }
  }
  TGraphErrors *cTOF_Y = new TGraphErrors( Nchan, Y, Yval, posErr, Yerr );
  cTOF_Y->GetXaxis()->SetLimits(0,Nchan);  
  cTOF_Y->GetYaxis()->SetLimits(30.,45.);
  cTOF_Y->SetTitle("Time of Flight Neutron");
  cTOF_Y->GetXaxis()->SetTitle("Channel");
  cTOF_Y->GetYaxis()->SetTitle("Time of Flight (ns)");
  cTOF_Y->SetMarkerStyle(20); // idx 20 Circles, idx 21 Boxes
  cTOF_Y->Write("cTOF_N");

  cout << TOFmin_n << endl << endl;
  for( Int_t v=0; v<Nchan; v++ ){ 
    Yval[v]-=TOFmin_n;
    cout << Yval[v] << endl;
  }

  TGraphErrors *cTOF_Y_corr = new TGraphErrors( Nchan, Y, Yval, posErr, corrErr );
  cTOF_Y_corr->GetXaxis()->SetLimits(0,Nchan);  
  cTOF_Y_corr->GetYaxis()->SetLimits(35.,45.);
  cTOF_Y_corr->SetTitle("TOF Corrections Neutron");
  cTOF_Y_corr->GetXaxis()->SetTitle("Channel");
  cTOF_Y_corr->GetYaxis()->SetTitle("Correction (ns)");
  cTOF_Y_corr->SetMarkerStyle(20); // idx 20 Circles, idx 21 Boxes
  cTOF_Y_corr->Write("cTOF_N_corr");  

  fout->Write();

  st->Stop();

  cout << "Analysis complete. Results written to file: simTOFout.root" << endl;

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

}
