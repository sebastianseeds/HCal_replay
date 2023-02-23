//SSeeds 10.7.22 Check on energy of each participating member of the primary cluster over cuts on trigger sums
//sseeds 10.14.22 - Cutoff around 100 MeV supported by lowest n sampled E at 97 MeV
//Judicious choice for threshold should be 90-95 MeV

#include <ctime>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <TH2F.h>
#include <TChain.h>
#include <TCanvas.h>
#include <iostream>
#include <TSystem.h>
#include <TStopwatch.h>
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TString.h"
#include "TCut.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TVector3.h"
#include "TGraphErrors.h"
#include "TLorentzVector.h"

const Int_t kNcell = 288; //Total number of HCal modules
const Int_t kNrows = 24; //Total number of HCal rows
const Int_t kNcols = 12; //Total number of HCal columns
const Int_t kNtrack = 100; //Reasonable number of tracks in BigBite
const Int_t kNsum = 11; //Total number of HCal 2x2 sector sums (10) + totalsum (1)
const Int_t kNsect = 18; //Total number of HCal 4x4 block sectors
const Int_t kNsrows = 6; //Total number of sector rows
const Int_t kNscols = 3; //Total number of sector columns
const Double_t Xi = -2.20; //Distance from beam center to top of HCal in m
const Double_t Xf = 1.47; //Distance from beam center to bottom of HCal in m
const Double_t Yi = -0.853; //Distance from beam center to opposite-beam side of HCal in m
const Double_t Yf = 0.853; //Distance from beam center to beam side of HCal in m
const Double_t M_p = 0.938272; //Mass proton
const Double_t PI = TMath::Pi(); //TMath pi
	
/////////////////////
///HCal Sector Map///
//Each sector 4x4bl//
//bl (0,0) topright//
//Front side f/beam//
/////////////////////
/*  3  |  2  |  1  */
/* ----2-----1---- */
/*  6  |  5  |  4  */
/* ----4-----3---- */
/*  9  |  8  |  7  */
/* ----6-----5---- */
/* 12  | 11  | 10  */
/* ----8-----7---- */
/* 15  | 14  | 13  */
/* ---10-----9---- */
/* 18  | 17  | 16  */
/////////////////////
///Each # between ///
///  is a trigID  ///
/////////////////////
  
//Declare sector map locations in 2d array
const int sMap[6][3] = {{3,2,1},{6,5,4},{9,8,7},{12,11,10},{15,14,13},{18,17,16}};
  
//Main
void sumE( const char *configfilename="ssumE.cfg", int run = -1 ){

  //Start the chain for root files passed with config file
  TChain *C = new TChain("T");
  
  Double_t E_e = 1.92; //Energy of beam (incoming electrons from accelerator)
  Int_t SBS_field = 0; //Strength in percent of 2100A of SBS magnet
  Double_t W_mean = M_p; //With perfect optics, this should be true. Will read in until then by fitting W distribution on each run
  Double_t W_sigma = 0.03; //Reasonable value by default, but will read in from W distribution

  //Reading config file
  ifstream configfile(configfilename);
  TString currentline;
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      C->Add(currentline);
      cout << "Loaded file: " << currentline <<  endl;
    }    
  }
  TCut globalcut = "";
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endcut") ){
    if( !currentline.BeginsWith("#") ){
      globalcut += currentline;
      cout << "Applying the following cut to all data: " << globalcut << endl;
    }    
  }
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("#") ){
    TObjArray *tokens = currentline.Tokenize(" ");
    int ntokens = tokens->GetEntries();
    if( ntokens>1 ){
      TString skey = ( (TObjString*)(*tokens)[0] )->GetString();
      if( skey == "E_e" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	E_e = sval.Atof();
      }
      if( skey == "SBS_field" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	SBS_field = sval.Atoi();
      }
      if( skey == "W_mean" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	W_mean = sval.Atof();
      }
      if( skey == "W_sigma" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	W_sigma = sval.Atof();
      }
      delete tokens;
    }
  }

  //Declare vars
  Double_t atime[kNcell], row[kNcell], col[kNcell], id[kNcell], tdctime[kNcell], en[kNcell], Tr_ap[kNcell], TrID[kNcell];
  Double_t nblk, nclus, SHnclus, PSnclus, hcalx, hcaly, hcale;
  Double_t Tn;
  Double_t Tpx[kNtrack], Tpy[kNtrack], Tpz[kNtrack], Tp[kNtrack];
  Int_t Trdata;
  Double_t eW2;

  //Switch on branches
  C->SetBranchStatus( "*", 0 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.atime", 1 );  
  C->SetBranchStatus( "sbs.hcal.clus_blk.row", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.col", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.id", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.tdctime", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.e", 1 ); //check
  C->SetBranchStatus( "sbs.hcal.x", 1 );
  C->SetBranchStatus( "sbs.hcal.y", 1 );
  C->SetBranchStatus( "sbs.hcal.e", 1 );
  C->SetBranchStatus( "sbs.hcal.nblk", 1 ); //check
  C->SetBranchStatus( "sbs.hcal.nclus", 1 );
  C->SetBranchStatus( "sbs.trig.a_p", 1 ); //check
  C->SetBranchStatus( "sbs.trig.adcelemID", 1 ); //check
  C->SetBranchStatus( "Ndata.sbs.trig.adcelemID", 1 ); //check
  C->SetBranchStatus( "bb.sh.nclus", 1 );
  C->SetBranchStatus( "bb.ps.nclus", 1 );
  C->SetBranchStatus( "bb.tr.n", 1 );
  C->SetBranchStatus( "bb.ps.e", 1 );
  C->SetBranchStatus( "bb.sh.e", 1 );
  C->SetBranchStatus( "bb.tr.px", 1 );
  C->SetBranchStatus( "bb.tr.py", 1 );
  C->SetBranchStatus( "bb.tr.pz", 1 );
  C->SetBranchStatus( "bb.tr.vz", 1 );
  C->SetBranchStatus( "bb.tr.p", 1 );
  C->SetBranchStatus( "e.kine.W2", 1);

  //Map branches to vars
  C->SetBranchAddress( "sbs.hcal.clus_blk.atime", atime );
  C->SetBranchAddress( "sbs.hcal.clus_blk.row", row );
  C->SetBranchAddress( "sbs.hcal.clus_blk.col", col);
  C->SetBranchAddress( "sbs.hcal.clus_blk.id", id);
  C->SetBranchAddress( "sbs.hcal.clus_blk.tdctime", tdctime );
  C->SetBranchAddress( "sbs.hcal.clus_blk.e", en ); //check
  C->SetBranchAddress( "sbs.hcal.x", &hcalx );
  C->SetBranchAddress( "sbs.hcal.y", &hcaly );
  C->SetBranchAddress( "sbs.hcal.e", &hcale ); 
  C->SetBranchAddress( "sbs.hcal.nblk", &nblk ); //check
  C->SetBranchAddress( "sbs.hcal.nclus", &nclus );
  C->SetBranchAddress( "sbs.trig.a_p", Tr_ap ); //check
  C->SetBranchAddress( "sbs.trig.adcelemID", TrID ); //check
  C->SetBranchAddress( "Ndata.sbs.trig.adcelemID", &Trdata ); //check
  C->SetBranchAddress( "bb.sh.nclus", &SHnclus );
  C->SetBranchAddress( "bb.ps.nclus", &PSnclus );
  C->SetBranchAddress( "bb.tr.n", &Tn );
  C->SetBranchAddress( "bb.tr.px", Tpx );
  C->SetBranchAddress( "bb.tr.py", Tpy );
  C->SetBranchAddress( "bb.tr.pz", Tpz );
  C->SetBranchAddress( "bb.tr.p", Tp );
  C->SetBranchAddress( "e.kine.W2", &eW2 );

  //Declare outfile
  TFile *fout = new TFile( "sumEout.root", "RECREATE" );

  //Declare event list to apply globalcut
  TEventList *elist = new TEventList("elist","Elastic Event List");
  C->Draw(">>elist",globalcut);

  //Initialize histograms
  TH1D *h_W = new TH1D( "h_W", "W", 50, 0, 5 );

  TH1D *h_Etot_s[kNsum];
  TH1D *h_NCEtot_s[kNsum];
  TH1D *h_Eblk_s[kNsum];
  TH1D *h_NCEblk_s[kNsum];
  TH2D *h_eVsId_s[kNsum];
  TH2D *h_NCeVsId_s[kNsum];

  for( int s=0; s<kNsum; s++){
    h_Etot_s[s] = new TH1D( Form("h_Etot_s%d",s), Form("HCal Clus Tot E - Sum %d; GeV",s), 200, 0.001, 0.5 );
    h_NCEtot_s[s] = new TH1D( Form("h_NCEtot_s%d",s), Form("NoCut HCal Clus Tot E - Sum %d; GeV",s), 200, 0.001, 0.5 );
    h_Eblk_s[s] = new TH1D( Form("h_Eblk_s%d",s), Form("HCal Clus Blk E - Sum %d; GeV",s), 200, 0., 0.5 );
    h_NCEblk_s[s] = new TH1D( Form("h_NCEblk_s%d",s), Form("NoCut HCal Clus Blk E - Sum %d; GeV",s), 200, 0., 0.5 );
    h_eVsId_s[s] = new TH2D( Form("h_eVsId_s%d",s), Form("HCal Cluster E vs Block - Sum %d; GeV",s), 288, 0.5, 288.5, 200, 0, 0.5 );
    h_NCeVsId_s[s] = new TH2D( Form("h_NCeVsId_s%d",s), Form("NoCut HCal Cluster E vs Block - Sum %d; GeV",s), 288, 0.5, 288.5, 200, 0, 0.5 );
  }

  //Declare loop counters for limits and monitoring
  Long64_t Nevents = elist->GetN();
  Int_t yield = 0;

  cout << endl << "Proceeding to loop over all events in chain.." << endl;

  //Main loop over events
  for(Long64_t nevent = 0; nevent<Nevents; nevent++){

    //Monitoring
    if( nevent%10 == 0 ) cout << "Loop: " << nevent << "/" << Nevents << ". Elastics: " << yield << ". \r";
    cout.flush();

    //Get event
    C->GetEntry( elist->GetEntry( nevent ) ); 

    //Declare vars to be redefined on event
    bool elas = false;
    double W = eW2;
    h_W->Fill( W );

    //Liberal cut on W
    if( fabs( W-W_mean ) < W_sigma ){
      elas=true;
      yield++;
    }

    //Loop over all sum triggers on this event
    for( Int_t b=0; b<Trdata; b++){
      
      int trigID = (int)TrID[b];

      Double_t totE = 0.;
      Double_t NCtotE = 0.;

      //Continue to next trigger if the integrated adc value isn't positive
      if( Tr_ap[b]<=0. ){
	continue;
      }

      //Loop through all blocks in the primary HCal cluster for this event
      for( Int_t i=0; i<nblk; i++ ){
	
	//Declare cluster specific vars
	Double_t clusID = id[i];
	Double_t clusEN = en[i];
	Int_t sector = -1;

	//Cut on coincidence time peak in HCal ADC
	bool coin = atime[b]>40. && atime[b]<65.;
		
	//Reverse the columns to match the map (blocks start at top right)
	int tcol = (int)col[i]/4;
	if( tcol==0 ){
	  tcol=2;
	}else if( tcol==2 ){
	  tcol=0;
	}

	//Define the sector the cluster block is in
	sector = sMap[(int)(row[i]/4)][tcol];

	if( sector<0 || sector>kNsect ){
	  cerr << "ERROR: Cluster member out of acceptance. Exiting..." << endl;
	  break;
	}

	//Fill histograms for total sum trigger
	if( trigID==0 ) {
	  h_NCeVsId_s[0]->Fill( clusID, clusEN );
	  h_NCEblk_s[0]->Fill( clusEN );
	  NCtotE += clusEN;
	  if( elas==true && coin==true ){
	    h_eVsId_s[0]->Fill( clusID, clusEN );
	    h_Eblk_s[0]->Fill( clusEN );
	    totE += clusEN;
	  }
	}

	//Fill histograms for each individual trigger sum
	for( int s=1; s<kNsum; s++){
	  if( trigID==s ){ 
	    int cStart = (trigID+1)%2;
	    if( cStart == 0 ) cStart = 2;
	    int rStart = (trigID-1)/2;

	    //Check to see if block is in trigger sum region then fill histos and sum up energy
	    if( sector==sMap[rStart][cStart] || 
		sector==sMap[rStart][cStart-1] || 
		sector==sMap[rStart+1][cStart] || 
		sector==sMap[rStart+1][cStart-1] ){
	      h_NCEblk_s[s]->Fill( clusEN );
	      NCtotE += clusEN;
	      if( elas==true && coin==true ){
		h_Eblk_s[s]->Fill( clusEN );
		totE += clusEN;
	      }
	    }
	    h_NCeVsId_s[s]->Fill( clusID, clusEN );
	    if( elas==true && coin==true )
	      h_eVsId_s[s]->Fill( clusID, clusEN );
	  } 
	}
      }	  
      //Fill histos with summed energy
      for( int s=0; s<kNsum; s++){
	if( trigID==s ){
	  h_Etot_s[s]->Fill(totE);
	  h_NCEtot_s[s]->Fill(NCtotE);
	}
      }
    } 
  }
   
  //Make a few canvases to hold all important plots
  TCanvas *sum_Etot = new TCanvas("sum_Etot","Energy Total by Sum",1600,1200);
  sum_Etot->Divide(2,5);
  gStyle->SetOptStat(0);
  for( int i=0; i<kNsum-1; i++ ){
    sum_Etot->cd(i+1);
    h_Etot_s[i+1]->Draw();    
  }
  sum_Etot->Write();
  //gPad->Update();
  TCanvas *sum_NCEtot = new TCanvas("sum_NCEtot","NoCut Energy Total by Sum",1600,1200);
  sum_NCEtot->Divide(2,5);
  for( int i=0; i<kNsum-1; i++ ){   
    sum_NCEtot->cd(i+1);
    h_NCEtot_s[i+1]->Draw();
    h_Etot_s[i+1]->SetLineColor(kGreen);
    h_Etot_s[i+1]->Draw("same");
  }
  sum_NCEtot->Write();
  //gPad->Update();
  TCanvas *Etot = new TCanvas("Etot","Energy Total",1600,1200);
  Etot->cd();
  h_Etot_s[0]->Draw();
  Etot->Write();
  //gPad->Update();
  TCanvas *NCEtot = new TCanvas("NCEtot","NoCut Energy Total",1600,1200);
  NCEtot->cd();
  h_NCEtot_s[0]->Draw();
  h_Etot_s[0]->SetLineColor(kGreen);
  h_Etot_s[0]->Draw("same");
  NCEtot->Write();

  //gPad->Update();

  fout->Write();
  fout->Close();
  
  cout << "Histograms populated and written to file: sumEout.root" << endl;
  
}
