#include <TH2F.h>
#include <TChain.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <TSystem.h>
#include "hcal.h"
#include <vector>
#include <TStopwatch.h>
#include <TMath.h>
#include <TF1.h>
#include<TStyle.h>
#include<TList.h>
#include<TChain.h>
using namespace std;

#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <TChain.h>
#include <TSystem.h>
#include <TStopwatch.h>
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TCut.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TLegend.h"
#include "TVector3.h"
#include "TGraphErrors.h"
//#include "hcal.h"

void HCal_BBCal_Time_Diff( const char *configfilename, int run = -1 )
{
  //Define a new stopwatch.
  TStopwatch *st=new TStopwatch();
  st->Start(kTRUE);

  // Start the chain for root files passed with config file
  TChain *T = 0;

  if( !T ) {
    T = new TChain( "T" );
    T->SetBranchStatus( "*", 0 );
    T->SetBranchStatus( "sbs.hcal.*", 1 );
    T->SetBranchStatus( "bb.tr.*", 1 );
    //T->SetBranchStatus( "bb.tdctrig.tdc", 1 );
    //T->SetBranchStatus( "bb.tdctrig.tdcelemID", 1 );
    //T->SetBranchStatus( "Ndata.bb.tdctrig.tdcelemID", 1 );
    T->SetBranchAddress( "bb.tr.chi2", hcalt::BBtr_chi2 );
    T->SetBranchAddress( "bb.tr.n", hcalt::BBtr_n );
    T->SetBranchAddress( "bb.tr.px", hcalt::BBtr_px );
    T->SetBranchAddress( "bb.tr.py", hcalt::BBtr_py );
    T->SetBranchAddress( "bb.tr.pz", hcalt::BBtr_pz );
    T->SetBranchAddress( "bb.tr.p", hcalt::BBtr_p );
    T->SetBranchAddress( "bb.tr.tg_th", hcalt::BBtr_tg_th );
    T->SetBranchAddress( "bb.tr.tg_ph", hcalt::BBtr_tg_ph );
    T->SetBranchAddress( "sbs.hcal.clus.id", hcalt::cid );
    T->SetBranchAddress( "sbs.hcal.clus.row", hcalt::crow );
    T->SetBranchAddress( "sbs.hcal.clus.col", hcalt::ccol );
    T->SetBranchAddress( "sbs.hcal.clus.e", hcalt::ce );
    T->SetBranchAddress( "sbs.hcal.clus.eblk", hcalt::ceblk ); // Highest energy block
    T->SetBranchAddress( "sbs.hcal.clus.nblk", hcalt::cnblk );
    T->SetBranchAddress( "sbs.hcal.nclus", hcalt::nclus ); 
    T->SetBranchAddress( "sbs.hcal.nblk", hcalt::nblk ); // Total number of blocks in cell
    T->SetBranchAddress( "sbs.hcal.clus_blk.id", hcalt::cblkid ); // kNcell-1 index for each block
    T->SetBranchAddress( "sbs.hcal.clus_blk.e", hcalt::cblke ); // Array of block energies
    //T->SetBranchAddress( "bb.tdctrig.tdcelemID", hcalt::TDCT_id );
    //T->SetBranchAddress( "bb.tdctrig.tdc", hcalt::TDCT_tdc );
    //T->SetBranchAddress( "Ndata.bb.tdctrig.tdcelemID", &hcalt::TDCTndata );

    cout << "Opened up tree with nentries=" << T->GetEntries() << endl;
  }

  // Reading config file
  ifstream configfile(configfilename);
  TString currentline;
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      T->Add(currentline);
    }    
  }/*
  TCut globalcut = "";
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endcut") ){
    if( !currentline.BeginsWith("#") ){
      globalcut += currentline;
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
      if( skey == "HCal_d" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	HCal_d = sval.Atof();
      }
      if( skey == "HCal_th" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	HCal_th = sval.Atof();
      }
      if( skey == "opticsCorr" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	opticsCorr = sval.Atof();
      }
      if( skey == "W_mean" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        W_mean = sval.Atof();
      }
      if( skey == "W_sig" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        W_sig = sval.Atof();
      }
      if( skey == "ScaleFac" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        ScaleFac = sval.Atof();
      }
    }
    delete tokens;
    }*/


  //Stop and print out stopwatch information.
  st->Stop();
  cout<<"CPU time = "<<st->CpuTime()<<" s = "<<st->CpuTime()/60.<<" min   Real time = "<<st->RealTime()<<" s = "<<st->RealTime()/60.<<" min"<<endl;
}
