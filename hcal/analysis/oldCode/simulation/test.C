//SSeeds 7.13.22 - test to make sure I'm not sane
#include "TTree.h"
#include "TChain.h"
#include <iostream>

void test(){

//Set processing params
TChain *C = new TChain("T");
TFile *fout = new TFile( "simTOFout.root", "RECREATE" );
long nevent=0;

//Load in data
C->Add("/w/halla-scshelf2102/sbs/seeds/G4SBS/g4sbs_install/run_g4sbs_here/gmn_3.5GeV2_elastic_TOF.root");
//C->Add("/w/halla-scshelf2102/sbs/seeds/G4SBS/g4sbs_install/run_g4sbs_here/gmn_13.3GeV2.root");
//C->Add("/lustre19/expphy/volatile/halla/sbs/pdbforce/g4sbs_output/sbs4_elas_sbs30p/gmn_SBS4_elastic.root");

  //TFile *f = new TFile("/w/halla-scshelf2102/sbs/seeds/G4SBS/g4sbs_install/run_g4sbs_here/gmn_13.3GeV2.root");

  //Set tree variables
  Int_t Harm_HCalScint_hit_nhits;
  int Harm_HCalScint_hit_sdtridx[100];
  int Harm_HCalScint_hit_cell[100];
  double SDTrack_T[100];
  int SDTrack_PID[100];


  TBranch *b_Harm_HCalScint_hit_nhits;
  C->SetBranchAddress( "Harm.HCalScint.hit.nhits", &Harm_HCalScint_hit_nhits );
  C->SetBranchAddress( "Harm.HCalScint.hit.cell", Harm_HCalScint_hit_nhits );
  

  if( C->GetEntries()!=0 ){
    cout << "Opened file successfully." << endl;
  }else{
    cout << "Error: No file found." << endl;
    //return;
  }

  while( C->GetEntry( nevent++ ) ){ 
    
    if( nevent%1000 == 0 ) cout << "Processing event: " << nevent << "/" << C->GetEntries() << endl;
    /*
    for( int ihit=0; ihit<C->Harm_HCalScint_hit_nhits; ihit++ ){
      double_t BCtime = ( *(C->SDTrack_T) )[( *(C->Harm_HCalScint_hit_sdtridx) )[ ihit ]];
      cout << BCtime << endl;
    }
    */

    cout << Harm_HCalScint_hit_nhits << endl;
    cout << Harm_HCalScint_hit_cell << endl;


  }
}
