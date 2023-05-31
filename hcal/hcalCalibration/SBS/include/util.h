#ifndef UTIL_H
#define UTIL_H

#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <dirent.h>
#include "TVector3.h"
#include "TMatrixD.h"
#include "TLorentzVector.h"
#include "TLatex.h"
#include "TString.h"

namespace util {

  //Console
  string getDate();

  //General
  void readParam(std::string const_path,vector<Double_t> &param); //reads parameters from txt file
  Int_t countSets(std::string const_path,std::string type); //Counts the number of calibration sets for plotting
  void readDB(std::string const_path,std::string timestamp,std::string type,Double_t param[]); 
  void readDB(std::string const_path,std::string timestamp,std::string type,Double_t &param); //overload 1
  void readDB(std::string const_path,std::string timestamp,std::string type,vector<Double_t> &param); //overload 2
  void tsCompare(std::string tsA,std::string tsB,std::string &tsLate);

  //Functions to read csv files
  void ReadRunList(std::string runsheet_dir,    // Dir. path containing CSV files with run info
		   std::string experiment,      // SBS experiment {gmn, gen, etc.}
		   Int_t &nruns,                // No. of runs to analyze
		   Int_t sbsconf,               // SBS configuration
		   Int_t replay_pass,           // replay pass
		   Int_t verbose,               // verbosity
		   vector<calrun> &crun);       // Output: Vector of crun structs

  void ReadCutList(std::string cutsheet_dir,    // Dir. path containing CSV files with cut info
		   std::string experiment,      // experiment {gmn,gen,genrp,gep,etc.}
		   Int_t sbsconfig,             // SBS configuration
		   Int_t replay_pass,           // replay pass
		   std::string target,          // target
		   Int_t field,                 // sbs field in percent
		   Int_t verbose,               // verbosity
		   vector<calcut> &ccut);       // Output: Vector of ccut structs

  //Geometry/physics
  void sethcalaxes( Double_t sbstheta_rad,                 //SBS arm theta in radians
		    vector<TVector3> &hcal_axes );         //Pass 3-vector (unit) for modification in HCal coordinates

  void getxyhcalexpect( TVector3 vertex,                   //Unit 3-vector from the z-component of the vertex
			TVector3 pNhat,                    //Unit direction of momentum transfer to recoil nucleon
			TVector3 hcal_origin,              //Location of hcal origin in hall coordinates
			vector<TVector3> hcal_axes,        //Define the hcal axes (unit)
			vector<Double_t> &xyhcalexpect );  //Expected location of elastic recoil nucleon on hcal

  void checkPID( std::string target,                       //Target can sometimes disambiguate p/n peaks
		 Double_t dx0p,                            //Distribution dx location of proton peak
		 Double_t dx0n,                            //Distribution dx location of neutron peak
		 Double_t dxsigp,                          //Distribution dx proton sigma
		 Double_t dxsign,                          //Distribution dx neutron sigma
		 Double_t dx,                              //Event dx
		 Int_t &pid );                             //Modified pid(-1:neither,0:ambiguous,1:proton,2:neutron)
  
  //Fit functions
  Double_t g_gfit(Double_t *x, Double_t *par);
  Double_t g_gfit_bg(Double_t *x, Double_t *par);
  Double_t g_sgfit(Double_t *x, Double_t *par);
  Double_t g_sgfit_bg(Double_t *x, Double_t *par);
  Double_t g_expofit(Double_t *x, Double_t *par);
  Double_t g_twfit(Double_t *x, Double_t *par);
  Double_t g_scexpofit(Double_t *x, Double_t *par);
  Double_t g_lfit(Double_t *x, Double_t *par);
  Double_t g_p3fit(Double_t *x, Double_t *par);
  Double_t g_p4fit(Double_t *x, Double_t *par);
  Double_t g_p6fit(Double_t *x, Double_t *par);

}

#endif
