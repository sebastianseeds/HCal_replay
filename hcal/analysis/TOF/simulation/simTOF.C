//SSeeds 7.13.22 - Short script to extract the boundary crossing (BC) time from the sensitive detector (SD) track information from g4sbs simulations
//SSeeds 3.8.23 - Added TOFvR and TOFvp plots for better HCal TOF correction extraction. Iteration 0 skips fits for diagnosis. Iteration 1 fits the slices of TOF vs p for evaluation of functional forms. 
//SSeeds 3.9.23 - Diagnosis complete and const arrays filled. Shouldn't need iter==0 moving forward.

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
const Int_t maxTracks = 1000; // Reasonable limit on tracks to be stored per event
const Int_t maxTdcChan = 10; // Set to accomodate original 5 TDCTrig channels with buffer
const Double_t hcalheight = -0.2897; // Height of HCal above beamline
const Int_t Nchan = 288;
const Int_t Nrows = 24;
const Int_t Ncols = 12;
const Double_t hcalycenter = -5.8;

//Constants
const Double_t PI = TMath::Pi();
const Double_t M_e = 0.00051; //Mass electron (GeV)
const Double_t M_p = 0.938272; //Mass proton (GeV)
const Double_t M_n = 0.939565; //Mass neutron (GeV)

//Geant Parameters
//PDG PID, 2112=Neutron, 2212=Proton
//pdg.lbl.gov/2020/reviews/rpp2020-rev-monte-carlo-numbering.pdf
const Int_t prot = 2212;
const Int_t neut = 2112;


//For Proton Momentum and TOF over all kinematics
const Double_t pcompr = 2.5; //range of px, py, pz on either size of zero
const Double_t drllim = 4.5; //diagnostic radius lower limit
const Double_t drulim = 7.5; //diagnostic radius upper limit
const Int_t drbins = 120; //diagnostic radius number of bins
const Double_t dpllim = 0.0; //same scheme: momentum
const Double_t dpulim = 15.0; 
const Int_t dpbins = 150;
const Double_t dtllim = 0.0; //same scheme: time of flight (TOF)
const Double_t dtulim = 60.0;
const Int_t dtbins = 480;
const Double_t pllim[6] = {1.5,3.5,3.0,1.5,1.5,1.5}; //tuned momentum lower limit
const Double_t pulim[6] = {2.8,6.7,8.9,5.4,3.9,3.7}; //tuned momentum upper limit
const Double_t tllim[6] = {35,46,48,46,35,35}; //tuned TOF lower limit
const Double_t tulim[6] = {45,50,52,58,50,50}; //tuned TOF upper limit
const Double_t TOF_sig[6] = {0.174,0.190,0.172,0.173,0.274,0.208}; //measured TOF sigma
const Double_t alin_h[6] = {39.5,47.5,48.75,47.5,38.0,38.0}; //measured TOF value at upper bin
const Double_t alin_l[6] = {43.5,48.5,50.75,55.0,43.0,43.0}; //measured TOF value at lower bin
const Int_t Npbins[6] = {39,96,177,117,72,66}; //10*prange*3
const Int_t Ntbins[6] = {100,40,40,120,150,150}; //10*trange
//const Int_t Npbins = 39;
//const Int_t Ntbins = 100;

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

//Main, kine->kinematic, iter->0:no fitting 1:fitting
void simTOF( Int_t kine = 4, Int_t iter = 0 ){
  
  // Define a clock to check overall time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );
  
  // Get the date
  string date = getDate();

  //Indexing kinematics for processing
  Int_t kIdx;
  if( kine == 4 ) kIdx=0;
  if( kine == 7 ) kIdx=1;
  if( kine == 11 ) kIdx=2;
  if( kine == 14 ) kIdx=3;
  if( kine == 8 ) kIdx=4;
  if( kine == 9 ) kIdx=5;

  //Set processing params
  TChain *C = new TChain("T");
  TFile *fout = new TFile( Form("output/simTOFout_sbs%d.root",kine), "RECREATE" );
  Int_t nevent=1;

  C->Add(Form("/lustre19/expphy/volatile/halla/sbs/seeds/simulation/gmn_sbs%d_v2_job*",kine));

  //Use root makeclass() functionality, no rush
  gmn_tree *T = new gmn_tree(C);

  if( C->GetEntries()!=0 ){
    cout << "Opened file successfully." << endl;
  }else{
    cout << "Error: No file found." << endl;
    return;
  }
  
  //Define some limits
  Double_t prange = pulim[kIdx]-pllim[kIdx];
  Double_t trange = tulim[kIdx]-tllim[kIdx];
  //Int_t Npbins = 10*prange*3; //empirical
  //Int_t Ntbins = 10*trange; //empirical

  //Define diagnostic histograms
  TH1D *TOF_all = new TH1D("TOF_all","Time of Flight over All Modules", dtbins, dtllim, dtulim);
  TH2D *TOF_vs_ID = new TH2D("TOF_vs_ID","Time of Flight vs Channel; channel; ns", Nchan, 0, Nchan, dtbins, dtllim, dtulim);
  TH2D *TOF_vs_R = new TH2D("TOF_vs_R","Time of Flight vs Radius from HCal Origin; R (m); TOF (ns)", drbins, drllim, drulim, dtbins, dtllim, dtulim);
  TH2D *R_vs_p = new TH2D("R_vs_p","Radius from HCal Origin vs Hadron Momentum; p_{had} (GeV); R (m)", drbins, drllim, drulim, dpbins, dpllim, dpulim);
  TH2D *TOF_vs_p = new TH2D("TOF_vs_p","Time of Flight vs Hadron Momentum", dpbins, dpllim, dpulim, dtbins, dtllim, dtulim);
  TH2D *cvy = new TH2D("cvy","col vs y",Ncols,0,Ncols,180,-6.7,-4.9); //empirial y limits
  TH1D *hp = new TH1D("hp","Hadron p; GeV",dpbins,dpllim,dpulim);
  TH2D *TOF_vs_px = new TH2D("TOF_vs_px","Time of Flight vs Hadron x Momentum", 100, -pcompr, pcompr, dtbins, dtllim, dtulim);
  TH2D *TOF_vs_py = new TH2D("TOF_vs_py","Time of Flight vs Hadron y Momentum", 100, -pcompr, pcompr, dtbins, dtllim, dtulim);
  TH2D *TOF_vs_pz = new TH2D("TOF_vs_pz","Time of Flight vs Hadron z Momentum", 100, -pcompr, pcompr, dtbins, dtllim, dtulim);

  //Define TOF histograms
  TH1D *TOF_all_p = new TH1D("TOF_all_p","Time of Flight over All Modules: Protons", Ntbins[kIdx], tllim[kIdx], tulim[kIdx]);
  TH1D *TOF_all_n = new TH1D("TOF_all_n","Time of Flight over All Modules: Neutrons", Ntbins[kIdx], tllim[kIdx], tulim[kIdx]);
  TH2D *TOF_vs_ID_p = new TH2D("TOF_vs_ID_p","Time of Flight vs Channel: Protons; channel; ns",Nchan,0,Nchan,Ntbins[kIdx],tllim[kIdx],tulim[kIdx]);
  TH2D *TOF_vs_ID_n = new TH2D("TOF_vs_ID_n","Time of Flight vs Channel: Neutrons; channel; ns",Nchan,0,Nchan,Ntbins[kIdx],tllim[kIdx],tulim[kIdx]);
  TH2D *TOF_vs_row_p = new TH2D("TOF_vs_row_p","Time of Flight vs Row Proton; row; ns",Nrows,0,Nrows,Ntbins[kIdx],tllim[kIdx],tulim[kIdx]);  
  TH2D *TOF_vs_col_p = new TH2D("TOF_vs_col_p","Time of Flight vs Column Proton; col; ns",Ncols,0,Ncols,Ntbins[kIdx],tllim[kIdx],tulim[kIdx]);
  TH2D *TOF_vs_row_n = new TH2D("TOF_vs_row_n","Time of Flight vs Row Neutron; row; ns",Nrows,0,Nrows,Ntbins[kIdx],tllim[kIdx],tulim[kIdx]);  
  TH2D *TOF_vs_col_n = new TH2D("TOF_vs_col_n","Time of Flight vs Column Neutron; col; ns",Ncols,0,Ncols,Ntbins[kIdx],tllim[kIdx],tulim[kIdx]);

  TH2D *TOF_vs_R_p = new TH2D("TOF_vs_R_p","Proton Time of Flight vs Radius from HCal Origin; R (m); TOF (ns)", drbins, 4.5, 7.5, Ntbins[kIdx], tllim[kIdx], tulim[kIdx]);
  TH2D *TOF_vs_R_n = new TH2D("TOF_vs_R_n","Neutron Time of Flight vs Radius from HCal Origin; R (m); TOF (ns)", drbins, 4.5, 7.5, Ntbins[kIdx], tllim[kIdx], tulim[kIdx]);

  TH2D *TOF_vs_p_zoom = new TH2D("TOF_vs_p_zoom","Time of Flight vs Hadron Momentum", Npbins[kIdx], pllim[kIdx], pulim[kIdx], Ntbins[kIdx], tllim[kIdx], tulim[kIdx]);
  TH2D *TOF_vs_pp = new TH2D("TOF_vs_pp","Time of Flight vs Proton Momentum", Npbins[kIdx], pllim[kIdx], pulim[kIdx], Ntbins[kIdx], tllim[kIdx], tulim[kIdx]);
  TH2D *TOF_vs_pn = new TH2D("TOF_vs_pn","Time of Flight vs Neutron Momentum", Npbins[kIdx], pllim[kIdx], pulim[kIdx], Ntbins[kIdx], tllim[kIdx], tulim[kIdx]);
  TH2D *TOF_vs_px_p = new TH2D("TOF_vs_px_p","Time of Flight vs Proton x Momentum", 100, -pcompr, pcompr, Ntbins[kIdx], tllim[kIdx], tulim[kIdx]);
  TH2D *TOF_vs_px_n = new TH2D("TOF_vs_px_n","Time of Flight vs Neutron x Momentum", 100, -pcompr, pcompr, Ntbins[kIdx], tllim[kIdx], tulim[kIdx]);
  TH2D *TOF_vs_py_p = new TH2D("TOF_vs_py_p","Time of Flight vs Proton y Momentum", 100,-pcompr, pcompr, Ntbins[kIdx], tllim[kIdx], tulim[kIdx]);
  TH2D *TOF_vs_py_n = new TH2D("TOF_vs_py_n","Time of Flight vs Neutron y Momentum", 100, -pcompr, pcompr, Ntbins[kIdx], tllim[kIdx], tulim[kIdx]);
  TH2D *TOF_vs_pz_p = new TH2D("TOF_vs_pz_p","Time of Flight vs Proton z Momentum", 100,-pcompr, pcompr, Ntbins[kIdx], tllim[kIdx], tulim[kIdx]);
  TH2D *TOF_vs_pz_n = new TH2D("TOF_vs_pz_n","Time of Flight vs Neutron z Momentum", 100, -pcompr, pcompr, Ntbins[kIdx], tllim[kIdx], tulim[kIdx]);

  //Momentum component slices
  TH2D *TOF_vs_px_cent = new TH2D("TOF_vs_px_cent","Time of Flight vs Hadron x Momentum, Central py/pz", 100, -pcompr, pcompr, Ntbins[kIdx], tllim[kIdx], tulim[kIdx]);
  TH2D *TOF_vs_py_cent = new TH2D("TOF_vs_py_cent","Time of Flight vs Hadron y Momentum, Central px/pz", 100, -pcompr, pcompr, Ntbins[kIdx], tllim[kIdx], tulim[kIdx]);
  TH2D *TOF_vs_pz_cent = new TH2D("TOF_vs_pz_cent","Time of Flight vs Hadron z Momentum, Central px/py", 100, -pcompr, pcompr, Ntbins[kIdx], tllim[kIdx], tulim[kIdx]);


  //Define supplementary histograms
  TH1D *hpX = new TH1D("hpX","Hadron px; GeV",100,-pcompr,pcompr);
  TH1D *hpY = new TH1D("hpY","Hadron py; GeV",100,-pcompr,pcompr);
  TH1D *hpZ = new TH1D("hpZ","Hadron pz; GeV",100,-pcompr,pcompr);
  TH1D *hvX = new TH1D("hvX","Number of Detections vs X; X (m)",100,-2,3);
  TH1D *pvX = new TH1D("pvX","Number of Protons vs X; X (m)",100,-2,3);
  TH1D *nvX = new TH1D("nvX","Number of Neutrons vs X; X (m)",100,-2,3);
  TH1D *hvY = new TH1D("hvY","Number of Detections vs Y; Y (m)",100,-8,-4);
  TH1D *pvY = new TH1D("pvY","Number of Protons vs Y; Y (m)",100,-8,-4);
  TH1D *nvY = new TH1D("nvY","Number of Neutrons vs Y; Y (m)",100,-8,-4);
  
  //Define loop parameters
  Long64_t Nevents = C->GetEntries();

  cout << "Opened tree with " << Nevents << " simulated events." << endl;

  while( T->GetEntry( nevent++ ) ){ 
    
    if( nevent%1000 == 0 ){
      cout << "Processing event: " << nevent << "/" << Nevents << "\r";
      cout.flush();
    }
    
    //Generate TOF from SD track information
    for( Int_t ihit=0; ihit<T->Harm_HCalScint_hit_nhits; ihit++ ){ //loop over all hits in HCal
      Double_t BCtime = ( *(T->SDTrack_T) )[( *(T->Harm_HCalScint_hit_sdtridx) )[ ihit ]]; //boundary crossing time
      Int_t PID = ( *(T->SDTrack_PID) )[( *(T->Harm_HCalScint_hit_sdtridx) )[ ihit ]]; //particle ID
      Double_t hadPosY = ( *(T->SDTrack_posx) )[( *(T->Harm_HCalScint_hit_sdtridx) )[ ihit ]]; //hadron position X
      Double_t hadPosX = ( *(T->SDTrack_posy) )[( *(T->Harm_HCalScint_hit_sdtridx) )[ ihit ]]; //hadron position Y
      Double_t hadpx = ( *(T->SDTrack_momx) )[( *(T->Harm_HCalScint_hit_sdtridx) )[ ihit ]];
      Double_t hadpy = ( *(T->SDTrack_momy) )[( *(T->Harm_HCalScint_hit_sdtridx) )[ ihit ]];
      Double_t hadpz = ( *(T->SDTrack_momz) )[( *(T->Harm_HCalScint_hit_sdtridx) )[ ihit ]];
      Double_t hadPosR = sqrt(pow(hadPosX,2)+pow(hadPosY,2));
      Double_t hadp = sqrt(pow(hadpx,2)+pow(hadpy,2)+pow(hadpz,2));

      Int_t cell = (*(T->Harm_HCalScint_hit_cell))[ihit];

      R_vs_p->Fill(hadPosR,hadp);

      hp->Fill(hadp);
      hpX->Fill(hadpx);
      hpY->Fill(hadpy);
      hpZ->Fill(hadpz);
      cvy->Fill( cell, hadPosY );
      hvY->Fill( hadPosY );
      hvX->Fill( hadPosX );
      TOF_all->Fill( BCtime );
      TOF_vs_ID->Fill( cell, BCtime );
      TOF_vs_R->Fill( hadPosR, BCtime );
      TOF_vs_p->Fill( hadp, BCtime );
      TOF_vs_p_zoom->Fill( hadp, BCtime );
      TOF_vs_px->Fill( hadpx, BCtime );
      TOF_vs_py->Fill( hadpy, BCtime );
      TOF_vs_pz->Fill( hadpz, BCtime );
      if( PID == prot ){
	pvY->Fill( hadPosY );
	pvX->Fill( hadPosX );
	TOF_all_p->Fill( BCtime );
	TOF_vs_ID_p->Fill( cell, BCtime );
	TOF_vs_row_p->Fill( (Int_t)cell/Ncols, BCtime );
	TOF_vs_col_p->Fill( cell%Ncols, BCtime );
	TOF_vs_R_p->Fill( hadPosR, BCtime );
	TOF_vs_pp->Fill( hadp, BCtime );
	TOF_vs_px_p->Fill( hadpx, BCtime );
	TOF_vs_py_p->Fill( hadpy, BCtime );
	TOF_vs_pz_p->Fill( hadpz, BCtime );
      }else if( PID == neut ){
	nvY->Fill( hadPosY );
	nvX->Fill( hadPosX );
	TOF_all_n->Fill( BCtime );
	TOF_vs_ID_n->Fill( cell, BCtime );
	TOF_vs_row_n->Fill( (Int_t)cell/Ncols, BCtime );
	TOF_vs_col_n->Fill( cell%Ncols, BCtime );
	TOF_vs_R_n->Fill( hadPosR, BCtime );
	TOF_vs_pn->Fill( hadp, BCtime );
	TOF_vs_px_n->Fill( hadpx, BCtime );
	TOF_vs_py_n->Fill( hadpy, BCtime );
	TOF_vs_pz_n->Fill( hadpz, BCtime );
      }
    }
  }
  
  if( iter==0 ){
    fout->Write();

    st->Stop();

    cout << "Prelimary analysis complete. Results written to file: simTOFout_sbs" << kine << ".root" << endl;

    // Send time efficiency report to console
    cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;
    return;
  }
    

  //Construct graphs
  Double_t posErr[Npbins[kIdx]];
  TF1 *f1;
  Double_t alin_div = (alin_l[kIdx]-alin_h[kIdx])/Npbins[kIdx];
  Double_t corrErr[Npbins[kIdx]];
  Int_t lNsig = 17;
  Int_t uNsig = 2;

  //For Proton Momentum
  Double_t X[Npbins[kIdx]];
  Double_t Xval[Npbins[kIdx]];
  Double_t Xerr[Npbins[kIdx]];
  TH1D *slice_p[Npbins[kIdx]];
  Double_t TOFmin_p = 100;
  for( Int_t x=0; x<=Npbins[kIdx]; x++ ){
    posErr[x] = 0.;
    corrErr[x] = 0.;
    Double_t afitm = alin_l[kIdx]-x*alin_div;
    Double_t afitl = afitm-lNsig*TOF_sig[kIdx];
    Double_t afith = afitm+uNsig*TOF_sig[kIdx];
    X[x] = x*(prange/Npbins[kIdx])+pllim[kIdx];
    slice_p[x] = TOF_vs_pp->ProjectionY(Form("slice_p_%d",x+1),x+1,x+1);
    slice_p[x]->Fit("gaus","Q","",afitl,afith);
    f1=slice_p[x]->GetFunction("gaus");
    if(slice_p[x]->GetEntries()>10){
      Xval[x] = f1->GetParameter(1);
      Xerr[x] = f1->GetParameter(2);
      if( Xval[x]>38&&Xval[x]<TOFmin_p ) TOFmin_p=Xval[x];
    }
  }
  TGraphErrors *cTOF_X = new TGraphErrors( Npbins[kIdx], X, Xval, posErr, Xerr );
  cTOF_X->GetXaxis()->SetLimits(0,Npbins[kIdx]);  
  cTOF_X->GetYaxis()->SetLimits(pllim[kIdx],pulim[kIdx]);
  cTOF_X->SetTitle("Time of Flight vs Proton Momentum");
  cTOF_X->GetXaxis()->SetTitle("p (GeV)");
  cTOF_X->GetYaxis()->SetTitle("Time of Flight (ns)");
  cTOF_X->SetMarkerStyle(20); // idx 20 Circles, idx 21 Boxes
  cTOF_X->Write("cTOF_P");

  cout << TOFmin_p << endl << endl;
  for( Int_t v=0; v<Npbins[kIdx]; v++ ){
    Xval[v]-=TOFmin_p;
    cout << Xval[v] << endl;
  }

  TGraphErrors *cTOF_X_corr = new TGraphErrors( Npbins[kIdx], X, Xval, posErr, corrErr );
  cTOF_X_corr->GetXaxis()->SetLimits(0,Npbins[kIdx]);  
  cTOF_X_corr->GetYaxis()->SetLimits(pllim[kIdx],pulim[kIdx]);
  cTOF_X_corr->SetTitle("TOF vs P Corrections Proton");
  cTOF_X_corr->GetXaxis()->SetTitle("p (GeV)");
  cTOF_X_corr->GetYaxis()->SetTitle("Correction (ns)");
  cTOF_X_corr->SetMarkerStyle(21); // idx 20 Circles, idx 21 Boxes
  cTOF_X_corr->Write("cTOF_P_corr");
  
  //For Neutrons
  Double_t Y[Npbins[kIdx]];
  Double_t Yval[Npbins[kIdx]];
  Double_t Yerr[Npbins[kIdx]];
  TH1D *slice_n[Npbins[kIdx]];
  Double_t TOFmin_n=100;
  for( Int_t x=0; x<Npbins[kIdx]; x++ ){
    Double_t afitm = alin_l[kIdx]-x*alin_div;
    Double_t afitl = afitm-lNsig*TOF_sig[kIdx];
    Double_t afith = afitm+uNsig*TOF_sig[kIdx];
    Y[x] = x*(prange/Npbins[kIdx])+pllim[kIdx];
    slice_n[x] = TOF_vs_pn->ProjectionY(Form("slice_n_%d",x+1),x+1,x+1);
    slice_n[x]->Fit("gaus","Q","",afitl,afith);
    f1=slice_n[x]->GetFunction("gaus");
    if(slice_n[x]->GetEntries()>10){
      Yval[x] = f1->GetParameter(1);
      Yerr[x] = f1->GetParameter(2);
      if( Yval[x]>5&&Yval[x]<TOFmin_n ) TOFmin_n=Yval[x];
    }
  }
  TGraphErrors *cTOF_Y = new TGraphErrors( Npbins[kIdx], Y, Yval, posErr, Yerr );
  cTOF_Y->GetXaxis()->SetLimits(0,Npbins[kIdx]);  
  cTOF_Y->GetYaxis()->SetLimits(pllim[kIdx],pulim[kIdx]);
  cTOF_Y->SetTitle("Time of Flight vs Neutron Momentum");
  cTOF_Y->GetXaxis()->SetTitle("p (GeV)");
  cTOF_Y->GetYaxis()->SetTitle("Time of Flight (ns)");
  cTOF_Y->SetMarkerStyle(20); // idx 20 Circles, idx 21 Boxes
  cTOF_Y->Write("cTOF_N");

  cout << TOFmin_n << endl << endl;
  for( Int_t v=0; v<Npbins[kIdx]; v++ ){ 
    Yval[v]-=TOFmin_n;
    cout << Yval[v] << endl;
  }

  TGraphErrors *cTOF_Y_corr = new TGraphErrors( Npbins[kIdx], Y, Yval, posErr, corrErr );
  cTOF_Y_corr->GetXaxis()->SetLimits(0,Npbins[kIdx]);  
  cTOF_Y_corr->GetYaxis()->SetLimits(pllim[kIdx],pulim[kIdx]);
  cTOF_Y_corr->SetTitle("TOF vs P Corrections Neutron");
  cTOF_Y_corr->GetXaxis()->SetTitle("p (GeV)");
  cTOF_Y_corr->GetYaxis()->SetTitle("Correction (ns)");
  cTOF_Y_corr->SetMarkerStyle(20); // idx 20 Circles, idx 21 Boxes
  cTOF_Y_corr->Write("cTOF_N_corr");  

  fout->Write();

  st->Stop();

  cout << "Analysis complete. Results written to file: simTOFout_sbs" << kine << ".root" << endl;

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

}
