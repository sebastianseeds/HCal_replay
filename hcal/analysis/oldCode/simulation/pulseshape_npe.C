#include <iostream>
#include <vector>
#include "TChain.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TROOT.h"
#include "TObjArray.h"
#include "TChainElement.h"
#include "TClonesArray.h"
#include "TLegend.h"
#include "TF1.h"
#include "TF2.h"
#include "TMath.h"
#include "G4SBSRunData.hh"
//#include "../run_g4sbs_here/plsp_tree.C"
#include "gmn_tree.C"

void pulseshape_npe(){

  TChain *C = new TChain("T");
  C->Add("lustre19/expphy/volatile/halla/sbs/pdbforce/g4sbs_output/sbs4_elas_sbs30p/gmn_SBS4_elastic.root");
  //C->Add("../run_g4sbs_here/plsp_gmn13_5_HCal_100_1.root");
  //C->Add("../run_g4sbs_here/plsp_gmn13_5_HCal_100_2.root");
  //C->Add("../run_g4sbs_here/plsp_gmn13_5_HCal_100_3.root");
  string SDname = "Harm_HCal";

  G4SBSRunData *rd;
  long ngen = 0;
  float hcaldist = 0.0;
  int nfiles = 0;
  TObjArray *Flist = C->GetListOfFiles();
  TIter next(Flist);
  TChainElement *chEl = 0;
  set<TString> bad_file_list;

  while( (chEl=(TChainElement*)next()) ){
    TFile newfile(chEl->GetTitle());
    newfile.GetObject("run_data",rd);
    if( rd ){
      ngen +=rd->fNtries;
      nfiles++;
    } else{
      bad_file_list.insert( chEl->GetTitle());
    }    
  }

  double Ibeam = 30e-6; //A
  double weight = Ibeam/double(ngen)/1.602e-19;
 
  cout << "Total number of generated events = " << ngen << endl;

  gmn_tree *G = new gmn_tree(C);
  
  // TCanvas *c1[10000];
  // TGraph *g1[10000];
  // gROOT->SetBatch(); // This prevents canvas pop-ups

  // std::vector<double> time;
  // std::vector<double> edep;
  double time[25];
  double edep[25];
  char gtitle[100], cantitle[50];

  for ( int i=0; i<25; i++ ){
    time[ i ] = 0.0;
    edep[ i ] = 0.0;
  }

  Long64_t nevent = C->GetEntries();
 
  int ntimebins = 0;
  int caniter = 0, lasthit = 0;
  int iter = 0;
  for (Long64_t i=0; i<nevent; i++){
    C->GetEntry(i);     
    double timewindow = G->Harm_HCal_gatewidth;

    for (int ihit=0; ihit<G->Harm_HCal_hit_nhits; ihit++) { 

      int MID = (*(G->SDTrack_MID))[(*(G->Harm_HCalScint_hit_sdtridx))[ihit]];
      int PrTID = (*(G->PTrack_TID))[(*(G->Harm_HCalScint_hit_ptridx))[ihit]];
      int PID = (*(G->SDTrack_PID))[(*(G->Harm_HCalScint_hit_sdtridx))[ihit]];
      int PrPID = (*(G->PTrack_PID))[(*(G->Harm_HCalScint_hit_ptridx))[ihit]];
      
      if( (*(G->Harm_HCal_hit_NumPhotoelectrons))[ihit]>60 && PID == 2112 ){

	ntimebins = (*(G->Harm_HCal_hit_NPE_vs_time))[0].size();
	double wbin = timewindow/double(ntimebins);

	for (int tbin=0; tbin<ntimebins; tbin++) {

	  time[ tbin ] = tbin*wbin + 0.5*wbin ;
	  // edep.push_back( (*(G->Harm_HCal_hit_edep_vs_time))[ihit][tbin] );

	  edep[ tbin ] += (*(G->Harm_HCal_hit_NPE_vs_time))[ihit][tbin];	  
	  // if( tbin == 0) cout << " edep " << (*(G->Harm_HCal_hit_NPE_vs_time))[ihit][tbin] << " edep_tot " << edep[ tbin ] << endl;
	  
	}

	iter++;

	// sprintf(cantitle,"c%d",caniter);
	// c1[ihit] = new TCanvas(cantitle,cantitle,1000,800);
	// c1[ihit]->cd();
	// g1[ihit] = new TGraph( ntimebins, &(time[0]), &(edep[0]) );
	// sprintf(gtitle,"Event # %lld  |  Edep vs Time (%s)  |  hit # %d",i,SDname.c_str(),ihit);
	// g1[ihit]->SetTitle(gtitle);
	// g1[ihit]->SetLineColor(2);
	// g1[ihit]->SetLineWidth(2);
	// g1[ihit]->SetMarkerColor(1);
	// g1[ihit]->SetMarkerStyle(20);
	// g1[ihit]->GetXaxis()->SetTitle("Time (ns)");
	// g1[ihit]->GetYaxis()->SetTitle("edep (MeV)");
	// g1[ihit]->GetYaxis()->SetMaxDigits(2);
	// g1[ihit]->Draw("AP");

	// if ( caniter == 0 ) c1[ihit]->SaveAs(TString::Format("plsp_gmn_13_5_%s_10k.pdf[",SDname.c_str())); 
	// c1[ihit]->SaveAs(TString::Format("plsp_gmn_13_5_%s_10k.pdf[",SDname.c_str()));

	caniter += 1;
	lasthit = ihit;
	
      }
    }
    // time.clear();
    // edep.clear();
  }


  for ( int i=0; i<25; i++ ){
    edep[ i ] = edep[ i ]/iter;
  }

  // for ( int i=0; i<25; i++ ){
  //   cout << " time [ " << i << " ]= " << time[i] << " edep [ " << i << " ]= " << edep[i] << " iter " << caniter  << endl;

  // }

  cout << iter << endl;

  TCanvas *c1 = new TCanvas("c1","c1",1000,800);
  c1->cd();
  TGraph *g1 = new TGraph( ntimebins, time, edep );
  // sprintf(gtitle,"Event # %lld  |  NPE vs Time (%s)  |  hit # %d",i,SDname.c_str(),ihit);
  g1->SetTitle("NPE vs Time");
  g1->SetLineColor(2);
  g1->SetLineWidth(2);
  g1->SetMarkerColor(1);
  g1->SetMarkerStyle(20);
  g1->GetXaxis()->SetTitle("Time (ns)");
  g1->GetYaxis()->SetTitle("NPE");
  //g1->GetYaxis()->SetMaxDigits(2);
  g1->Draw("AP");
  
  
  // time.clear();
  // edep.clear();

  // c1[lasthit]->SaveAs(TString::Format("plsp_gmn_13_5_%s_10k.pdf[",SDname.c_str()));
  // cout << "No. of plots generated " << caniter << endl;

}


