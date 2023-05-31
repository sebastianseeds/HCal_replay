#include "../include/sbs.h"
#include <vector>
#include <iostream>
#include <algorithm>

void ctsettest(){
  
  std::string DBpath = gSystem->Getenv("DB_DIR");
  std::string tdcOffsetPath = DBpath + "/db_sbs.hcal.dat";
  std::string tdcCalibPath = DBpath + "/db_sbs.hcal.dat";
  std::string newtwPath = "../timing/parameters/tdctw_class_gmn_conf4_pass0.txt";
  std::string testPath = "../timing/parameters/test.txt";

  cout << "Looking for constants in file: " << newtwPath << endl;

  Double_t calib_const = -1;
  Double_t tdcoffsets[hcal::maxHCalChan] = {0.};
  vector<Double_t> tdctw;

  //std::string current_timestamp = "none";
  std::string current_timestamp = "-------[ 2021-10-24 04:30:00 ]";
  std::string tw_timestamp = "-------[ 2021-10-24 04:30:00 ]";
  //std::string current_timestamp = "--------[ 2024-07-26 00:00:00 ]";
  std::string tdc_calib_type = "sbs.hcal.tdc.calib";
  std::string tdc_offset_type = "sbs.hcal.tdc.offset";
  std::string timewalkA_type = "sbs.hcal.tdc.chan_tw_a";
  std::string timewalkB_type = "sbs.hcal.tdc.chan_tw_b";
  std::string timewalk_type = "sbs.hcal.tdc.tw";


  //util::readDB(tdcCalibPath,current_timestamp,tdc_calib_type,calib_const);
  //util::readDB(tdcOffsetPath,current_timestamp,tdc_offset_type,tdcoffsets);
  //util::readDB(newtwPath,current_timestamp,timewalkB_type,tdcoffsets);
  util::readDB(newtwPath,tw_timestamp,timewalk_type,tdctw);

  Int_t nset = util::countSets( newtwPath, timewalk_type );
  
  std::cout << "Number of calibration sets for " << timewalk_type << ": " << nset << std::endl;

}
