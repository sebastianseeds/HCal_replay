//sseeds 5.12.23 - script to consolidate many configuration files into one .csv file

#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <TChain.h>
#include <TSystem.h>
#include <TStopwatch.h>
#include <unistd.h>
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
#include "TLorentzVector.h"
#include "TString.h"
#include "TF1.h"

const Int_t nkine = 6;
const Int_t ntar = 2;
const Int_t nfields_lh2 = 13;
const Int_t nfields_ld2 = 12;
const Int_t kIdx[6] = {4, 7, 11, 14, 8, 9}; //Index of SBS run
const Double_t sampFrac[6] = {0.0641,0.0718,0.0734,0.0710,0.0686,0.0682};
const Double_t pCsigfactor[6] = {0.58,0.82,0.83,0.78,0.66,0.72};
const Int_t nfset_lh2[6] = {3,1,2,2,4,1}; //Number of LH2 magnetic field settings in order SBS{4,7,11,14,8,9}
const Int_t nfset_ld2[6] = {3,1,2,1,4,1}; //Number of LD2 magnetic field settings in order SBS{4,7,11,14,8,9}
const Int_t fset_lh2[6][4] = {{0,30,50,-1},
			      {85,-1,-1,-1},
			      {0,100,-1,-1},
			      {0,70,-1,-1},
			      {0,50,70,100},
			      {70,-1,-1,-1}};
const Int_t fset_ld2[6][4] = {{0,30,50,-1},
			      {85,-1,-1,-1},
			      {0,100,-1,-1},
			      {70,-1,-1,-1},
			      {0,50,70,100},
			      {70,-1,-1,-1}};


typedef struct{
  Int_t configure;
  std::string target;
  Int_t field;
  Double_t hcal_sf;
  Double_t hcal_es;
  Double_t W2_mean; // Mean of W at current kinematic
  Double_t W2_sig; // Width of W at current kinematic
  Double_t dx0_n; // Position of neutron spot, x-x_expected
  Double_t dx0_p; // Position of proton spot, x-x_expected
  Double_t dy0; // Position of hadron spot, y-y_expected
  Double_t dx_sig_n; // Max spread of neutron spot, x-x_expected
  Double_t dx_sig_p; // Max spread of proton spot, x-x_expected
  Double_t dy_sig; // Max spread of hadron spot, y-y_expected
  Double_t atime0; // Expected location in ADC time of signal
  Double_t atime_sig; // 1 sig of atime distribution
  Int_t useAlshield;
  std::string gcut;

} CON;

void csvGen( ){
  
  Int_t Nconfig = nfields_lh2+nfields_ld2;

  CON config[Nconfig];
  Int_t cidx = 0;

  std::string outputfilename = "gmncuts_pass1.csv";

  for( Int_t k=0; k<nkine; k++ ){
    
    Int_t kine = kIdx[k];
    
    for( Int_t t=0; t<ntar; t++ ){

      std::string tar;
      Int_t nfields;
      if( t==0 ){ 
	tar = "lh2";
	nfields = nfset_lh2[k];
      }else{
	tar = "ld2";
	nfields = nfset_ld2[k];
      }

      for( Int_t f=0; f<nfields; f++ ){

	Int_t fieldsetting;
	if( t==0 )
	  fieldsetting = fset_lh2[k][f];
	else
	  fieldsetting = fset_ld2[k][f];

	if( fieldsetting==-1 ) 
	  cout << "ERROR: -1 FIELD" << endl;

	std::string configfilename = Form("SBS%d/secal_%s_sbs%d_f%d.cfg",kine,tar.c_str(),kine,fieldsetting);
	
	config[cidx].hcal_sf = sampFrac[k];
	config[cidx].hcal_es = pCsigfactor[k];
	config[cidx].field = fieldsetting;
	config[cidx].target = tar;
	config[cidx].configure = kine;
  
	ifstream configfile(configfilename);
	TString currentline;
	while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endlist") ){
	  cout << "" << endl;
	}
	while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endcut") ){
	  if( !currentline.BeginsWith("#") ){
	    config[cidx].gcut += currentline;
	  }    
	}
	while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("#") ){
	  TObjArray *tokens = currentline.Tokenize(" ");
	  Int_t ntokens = tokens->GetEntries();
	  if( ntokens>1 ){
	    TString skey = ( (TObjString*)(*tokens)[0] )->GetString();
	    if( skey == "W2_mean" ){
	      TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	      config[cidx].W2_mean = sval.Atof();
	    }
	    if( skey == "W2_sig" ){
	      TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	      config[cidx].W2_sig = sval.Atof();
	    }
	    if( skey == "dx0_n" ){
	      TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	      config[cidx].dx0_n = sval.Atof();
	    }
	    if( skey == "dx0_p" ){
	      TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	      config[cidx].dx0_p = sval.Atof();
	    }
	    if( skey == "dy0" ){
	      TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	      config[cidx].dy0 = sval.Atof();
	    }
	    if( skey == "dx_sig_n" ){
	      TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	      config[cidx].dx_sig_n = sval.Atof();
	    }
	    if( skey == "dx_sig_p" ){
	      TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	      config[cidx].dx_sig_p = sval.Atof();
	    }
	    if( skey == "dy_sig" ){
	      TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	      config[cidx].dy_sig = sval.Atof();
	    }
	    if( skey == "atime0" ){
	      TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	      config[cidx].atime0 = sval.Atof();
	    }
	    if( skey == "atime_sig" ){
	      TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	      config[cidx].atime_sig = sval.Atof();
	    }
	    if( skey == "useAlshield" ){
	      TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	      config[cidx].useAlshield = sval.Atoi();
	    }
	  }
	  delete tokens;
	}
	cidx++;

      }//endloop fields
    }//endloop target
  }//endloop kinematic

  //Declare coeff outfile
  ofstream gmncutsfile;
  
  gmncutsfile.open( outputfilename );
  
  gmncutsfile << "sbsconfig,target,field,hcal_sf,hcal_es,W2_mean,W2_sig,dx0_n,dx0_p,dy0,dx_sig_n,dx_sig_p,dy_sig,atime0,atime_sig,useAlshield,gcut" << endl;

  for( Int_t n=0; n<Nconfig; n++ ){
    gmncutsfile << config[n].configure << "," << config[n].target << "," << config[n].field << "," << config[n].hcal_sf << "," << config[n].hcal_es << "," << config[n].W2_mean << "," << config[n].W2_sig << "," << config[n].dx0_n << "," << config[n].dx0_p << "," << config[n].dy0 << "," << config[n].dx_sig_n << "," << config[n].dx_sig_p << "," << config[n].dy_sig << "," << config[n].atime0 << "," << config[n].atime_sig << "," << config[n].useAlshield << "," << config[n].gcut << endl;
  }

  gmncutsfile.close();

  

}
