#ifndef HCAL_H
#define HCAL_H
#include <TPaveStats.h>

const Int_t MAX_FADC_SAMPLES = 250;
const Int_t MAX_HCAL_MODULES = 288;
const Int_t MAX_TRIGGER_MODULES = 16;
const Int_t MAX_HCAL_TDC_MODULES = 288;
const Int_t MAX_BB_MODULES = 241;

namespace hcalt {
  //HCAL core vars
  Double_t samps[MAX_HCAL_MODULES*MAX_FADC_SAMPLES+1000] = {0.};
  Double_t nsamps[MAX_HCAL_MODULES+1000] = {0.};
  Double_t row[MAX_HCAL_MODULES+1000] = {0.};
  Double_t col[MAX_HCAL_MODULES+1000] = {0.};
  Double_t samps_idx[MAX_HCAL_MODULES+1000] = {0.};
  Double_t a[MAX_HCAL_MODULES+1000] = {0.};
  Double_t a_p[MAX_HCAL_MODULES+1000] = {0.};
  Double_t a_c[MAX_HCAL_MODULES+1000] = {0.};
  Double_t ped[MAX_HCAL_MODULES+1000] = {0.};
  Double_t a_amp[MAX_HCAL_MODULES+1000] = {0.};
  Double_t a_amp_p[MAX_HCAL_MODULES+1000] = {0.};
  Int_t ndata = 0;
  Double_t ledbit = -1;
  Double_t ledcount = 0.0;
  Double_t atimeblk = 0.0;
  //Double_t ce = 0.0;
  Double_t tdc[MAX_HCAL_TDC_MODULES+100] = {0.};
  Double_t a_time[MAX_HCAL_MODULES+1000] = {0.};
  Double_t e[MAX_HCAL_MODULES+1000] = {0.};
  //Cluster vars - TODO, check memory allocation for these 
  Double_t nclus[MAX_HCAL_MODULES+1000] = {0.}; 
  Double_t nblk[MAX_HCAL_MODULES+1000] = {0.};
  Double_t cid[MAX_HCAL_MODULES+1000] = {0.};
  Double_t cblkid[MAX_HCAL_MODULES+1000] = {0.};
  Double_t crow[MAX_HCAL_MODULES+1000] = {0.};
  Double_t ccol[MAX_HCAL_MODULES+1000] = {0.};
  Double_t cbrow[MAX_HCAL_MODULES+1000] = {0.};
  Double_t cbcol[MAX_HCAL_MODULES+1000] = {0.};
  Double_t ce[MAX_HCAL_MODULES+1000] = {0.};
  Double_t ceblk[MAX_HCAL_MODULES+1000] = {0.};
  Double_t cnblk[MAX_HCAL_MODULES+1000] = {0.};
  Double_t cblke[MAX_HCAL_MODULES+1000] = {0.};
  Double_t cblkatime[MAX_HCAL_MODULES+1000] = {0.};
  Double_t catime[MAX_HCAL_MODULES+1000] = {0.};
  Double_t cx[MAX_HCAL_MODULES+1000] = {0.};
  Double_t cy[MAX_HCAL_MODULES+1000] = {0.};
  //Trigger generic detector vars
  Double_t Tsamps[MAX_TRIGGER_MODULES*MAX_FADC_SAMPLES+1000];
  Double_t Tnsamps[MAX_TRIGGER_MODULES+1000] = {0.};
  Double_t Tsamps_idx[MAX_TRIGGER_MODULES+1000] = {0.};
  Double_t Ta_p[MAX_TRIGGER_MODULES+1000] = {0.};
  Double_t Ta_amp_p[MAX_TRIGGER_MODULES+1000] = {0.};
  Double_t Ta_amp[MAX_TRIGGER_MODULES+1000] = {0.};
  Double_t TelemID[MAX_TRIGGER_MODULES+1000] = {0.};
  Double_t Tcol[MAX_TRIGGER_MODULES+1000] = {0.};
  Int_t Tndata = 0;
  //BB Trigger generic detector vars
  Int_t TDCTndata = 0;
  Double_t TDCT_id[MAX_TRIGGER_MODULES+1000] = {0.};
  Double_t TDCT_tdc[MAX_TRIGGER_MODULES+1000] = {0.};
  //BB track vars
  Double_t BBtr_n = 0;
  Double_t BBtr_px[MAX_BB_MODULES+1000] = {0.};
  Double_t BBtr_py[MAX_BB_MODULES+1000] = {0.};
  Double_t BBtr_pz[MAX_BB_MODULES+1000] = {0.};  
  Double_t BBtr_vx[MAX_BB_MODULES+1000] = {0.};
  Double_t BBtr_vy[MAX_BB_MODULES+1000] = {0.};
  Double_t BBtr_vz[MAX_BB_MODULES+1000] = {0.};
  Double_t BBtr_p[MAX_BB_MODULES+1000] = {0.};
  Double_t BBtr_chi2[MAX_BB_MODULES+1000] = {0.};
  Double_t BBtr_tg_th[MAX_BB_MODULES+1000] = {0.};
  Double_t BBtr_tg_ph[MAX_BB_MODULES+1000] = {0.};
  Double_t BBrow[MAX_BB_MODULES+1000] = {0.};
  Double_t BBcol[MAX_BB_MODULES+1000] = {0.};
  Int_t BBndata = 0;
  Double_t BBa_p[MAX_BB_MODULES+1000] = {0.};
  Double_t BBa_c[MAX_BB_MODULES+1000] = {0.};
  Double_t BBa_amp[MAX_BB_MODULES+1000] = {0.};
  Double_t BBceblk[MAX_BB_MODULES+1000] = {0.};
  Double_t BBcnblk[MAX_BB_MODULES+1000] = {0.};
  Double_t BBcblke[MAX_BB_MODULES+1000] = {0.};
  Double_t BBcblkid[MAX_BB_MODULES+1000] = {0.};
  Double_t BBcid[MAX_BB_MODULES+1000] = {0.};
  Double_t BBnblk[MAX_BB_MODULES+1000] = {0.};
  Double_t BBnclus[MAX_BB_MODULES+1000] = {0.}; 
  Double_t BBsh_x[MAX_BB_MODULES+1000] = {0.};
  Double_t BBsh_y[MAX_BB_MODULES+1000] = {0.};
  Double_t BBsh_e[MAX_BB_MODULES+1000] = {0.};
  Double_t BBps_x[MAX_BB_MODULES+1000] = {0.};
  Double_t BBps_y[MAX_BB_MODULES+1000] = {0.};
  Double_t BBps_e[MAX_BB_MODULES+1000] = {0.};

}

void fixStats()
{
  gPad->Update();
  TPaveStats *ps = (TPaveStats*)gPad->GetPrimitive("stats");
  if(ps) {
    ps->SetX1NDC(0.6);
    ps->SetY1NDC(0.55);
  }
}



#endif // HCAL_H
