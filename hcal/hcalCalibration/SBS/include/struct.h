#ifndef STRUCT_H
#define STRUCT_H

#include "../include/hcal.h"
#include "TVector3.h"
#include "TMatrixD.h"
#include "TVectorD.h"

typedef struct calrun {

  Int_t runnum;
  Int_t sbsconf;
  std::string target;
  Int_t sbsmag;           // SBS magnet current (A)
  Int_t bbmag;            // BB magnet current (A)
  Double_t ebeam;         // GeV, avg. over entire run 
  Double_t ebeam_std;     // GeV, std. dev. over entire run 
  Double_t charge;        // C, total charge collected by the run
  Double_t DAQltime;      // %
  std::string tdc_ts;     // hcal tdc offset timestamp
  std::string tdcc_ts;    // hcal tdc calibration constant timestamp
  std::string adct_ts;    // hcal adc time offset timestamp
  std::string adcg_ts;    // hcal adc gain coefficient timestamp

  // constructor 
  calrun(): 
  runnum(0),sbsconf(0),target("NONE"),sbsmag(0),bbmag(0),ebeam(0),ebeam_std(0),charge(0),DAQltime(0),tdc_ts("NONE"),tdcc_ts("NONE"),adct_ts("NONE"),adcg_ts("NONE")
  {}

  // define an ostream operator to print to screen conveniently
  friend ostream &operator <<(ostream &out, const calrun &crun) {
    out << " ------------------------------------------------------" << std::endl;
    out << " Run number                          : " << crun.runnum << std::endl;
    out << " SBS config                          : " << crun.sbsconf << std::endl;
    out << " Target                              : " << crun.target << std::endl;
    out << " SBS mag. cur. (A)                   : " << crun.sbsmag << std::endl;
    out << " BB mag. cur. (A)                    : " << crun.bbmag << std::endl;
    out << " Avg. ebeam (GeV)                    : " << crun.ebeam << std::endl;
    out << " Total Charge                        : " << crun.charge << std::endl;
    out << " HCal TDC offset timestamp           : " << crun.tdc_ts << std::endl;
    out << " HCal TDC cal const timestamp        : " << crun.tdcc_ts << std::endl;
    out << " HCal ADCt offset timestamp          : " << crun.adct_ts << std::endl;
    out << " HCal ADC gain coefficent timestamp  : " << crun.adcg_ts << std::endl;
    out << " ------------------------------------------------------" << std::endl << std::endl;
    return out;
  }

  // sets data by reading runsheet (exclusively for util::ReadRunList functions)
  void SetDataRunSheet(std::vector<std::string> data) {
    sbsconf = stoi(data[0]);
    runnum = stoi(data[1]);
    target = data[2];
    sbsmag = stoi(data[3]);
    bbmag = stoi(data[4]);
    ebeam = stod(data[5]);
    ebeam_std = stod(data[6]);
    charge = stod(data[7]);
    tdc_ts = data[8];
    tdcc_ts = data[9];
    adct_ts = data[10];
    adcg_ts = data[11];
  }

} calrun_t;

typedef struct calcut {

  Int_t sbsconfig;
  std::string target;
  Int_t field;           // SBS magnet current (A)
  Double_t hcal_sf;
  Double_t hcal_es;
  Double_t W2_mean;
  Double_t W2_sig;
  Double_t dx0_n;
  Double_t dx0_p;
  Double_t dy0;
  Double_t dx_sig_n;
  Double_t dx_sig_p;
  Double_t dy_sig;
  Double_t atime0;
  Double_t atime_sig;
  Int_t useAlshield;
  std::string gcut;

  // constructor 
  calcut(): 
  sbsconfig(0),target("NONE"),field(0),hcal_sf(0),hcal_es(0),W2_mean(0),W2_sig(0),dx0_n(0),dx0_p(0),dy0(0),dx_sig_n(0),dx_sig_p(0),dy_sig(0),atime0(0),atime_sig(0),useAlshield(0),gcut("NONE")
  {}

  // define an ostream operator to print to screen conveniently
  friend ostream& operator <<(ostream &out, const calcut& ccut) {
    out << " -------------------------------------------------" << std::endl;
    out << " SBS config                          : " << ccut.sbsconfig << std::endl;
    out << " Target                              : " << ccut.target << std::endl;
    out << " SBS field strength                  : " << ccut.field << std::endl;
    out << " HCal Samp Frac                      : " << ccut.hcal_sf << std::endl;
    out << " HCal E/sig ratio                    : " << ccut.hcal_es << std::endl;
    out << " W2 Mean                             : " << ccut.W2_mean << std::endl;
    out << " W2 sigma                            : " << ccut.W2_sig << std::endl;
    out << " dx mean neutron                     : " << ccut.dx0_n << std::endl;
    out << " dx mean proton                      : " << ccut.dx0_p << std::endl;
    out << " dy mean                             : " << ccut.dy0 << std::endl;
    out << " dx sigma neutron                    : " << ccut.dx_sig_n << std::endl;
    out << " dx sigma proton                     : " << ccut.dx_sig_p << std::endl;
    out << " dy sigma                            : " << ccut.dy_sig << std::endl;
    out << " ADC time mean                       : " << ccut.atime0 << std::endl;
    out << " ADC time sigma                      : " << ccut.atime_sig << std::endl;
    out << " Use Al Shield (0:out,1:in)          : " << ccut.useAlshield << std::endl;
    out << " Global Cut                          : " << ccut.gcut << std::endl;
    out << " -------------------------------------------------" << std::endl << std::endl;
    return out;
  }

  // sets data by reading runsheet (exclusively for util::ReadRunList functions)
  void SetDataCutSheet(std::vector<std::string> data) {
    sbsconfig = stoi(data[0]);
    target = data[1];
    field = stoi(data[2]);
    hcal_sf = stod(data[3]);
    hcal_es = stod(data[4]);
    W2_mean = stod(data[5]);
    W2_sig = stod(data[6]);
    dx0_n = stod(data[7]);
    dx0_p = stod(data[8]);
    dy0 = stod(data[9]);
    dx_sig_n = stod(data[10]);
    dx_sig_p = stod(data[11]);
    dy_sig = stod(data[12]);
    atime0 = stod(data[13]);
    atime_sig = stod(data[14]);
    useAlshield = stoi(data[15]);
    gcut = data[16];
  }

} calcut_t;  

typedef struct calset {
  
  //general members
  std::string timestamp;
  Double_t old_param[hcal::maxHCalChan];
  Double_t old_paramB[hcal::maxHCalChan];
  Double_t old_paramC[hcal::maxHCalChan];
  Double_t new_param[hcal::maxHCalChan];
  Double_t new_paramB[hcal::maxHCalChan]; //can be used for timewalk fit extraction
  //member for tdc alignment
  std::string calib_ts;
  Double_t tdc_calib;
  //members for adc gain
  Double_t new_param_oneblock[hcal::maxHCalChan];
  Double_t new_param_divide[hcal::maxHCalChan];
  TMatrixD Ma;
  TMatrixD Ma_oneblock;
  TVectorD ba;
  TVectorD ba_oneblock;
  Int_t NEV[hcal::maxHCalChan];
  Int_t NEV_oneblock[hcal::maxHCalChan];
  Double_t err[hcal::maxHCalChan];
  Double_t err_ev[hcal::maxHCalChan];
  Double_t err_oneblock[hcal::maxHCalChan];
  Double_t err_ev_oneblock[hcal::maxHCalChan];

  // constructor
  calset():
  timestamp("NONE"),old_param{},old_paramB{},old_paramC{},new_param{},new_paramB{},calib_ts("NONE"),tdc_calib(0),new_param_oneblock{},new_param_divide{},Ma(),Ma_oneblock(),ba(),ba_oneblock(),NEV{},NEV_oneblock{},err{},err_ev{},err_oneblock{},err_ev_oneblock{}
  {}

} calset_t;

#endif
