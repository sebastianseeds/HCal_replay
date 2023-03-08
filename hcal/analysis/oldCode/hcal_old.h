#ifndef HCAL_H
#define HCAL_H
#include <TPaveStats.h>

const Int_t MAX_FADC_SAMPLES = 250;
const Int_t MAX_HCAL_MODULES = 288;
const Int_t MAX_TRIGGER_MODULES = 16;
const Int_t MAX_HCAL_TDC_MODULES = 288;

namespace hcalt {
  //HCAL core vars
  Double_t samps[MAX_HCAL_MODULES*MAX_FADC_SAMPLES+1000];
  Double_t nsamps[MAX_HCAL_MODULES+1000] = {0};
  Double_t row[MAX_HCAL_MODULES+1000] = {0};
  Double_t col[MAX_HCAL_MODULES+1000] = {0};
  Double_t samps_idx[MAX_HCAL_MODULES+1000] = {0};
  Double_t a[MAX_HCAL_MODULES+1000] = {0};
  Double_t a_p[MAX_HCAL_MODULES+1000] = {0};
  Double_t ped[MAX_HCAL_MODULES+1000] = {0};
  Double_t a_amp[MAX_HCAL_MODULES+1000] = {0};
  Double_t a_amp_p[MAX_HCAL_MODULES+1000] = {0};
  Int_t ndata = 0;
  Double_t ledbit = -1;
  Double_t ledcount = 0;
  Double_t tdc[MAX_HCAL_TDC_MODULES+100];
  //Cluster vars
  Double_t nclus[MAX_HCAL_MODULES+1000] = {0};
  Double_t nblk[MAX_HCAL_MODULES+1000] = {0};
  Double_t cid[MAX_HCAL_MODULES+1000] = {0};
  Double_t cblkid[MAX_HCAL_MODULES+1000] = {0};
  Double_t crow[MAX_HCAL_MODULES+1000] = {0};
  Double_t ccol[MAX_HCAL_MODULES+1000] = {0};
  Double_t ce[MAX_HCAL_MODULES+1000] = {0};
  Double_t ceblk[MAX_HCAL_MODULES+1000] = {0};
  Double_t cnblk[MAX_HCAL_MODULES+1000] = {0};
  Double_t cblke[MAX_HCAL_MODULES+1000] = {0};
  //Trigger generic detector vars
  Double_t Tsamps[MAX_TRIGGER_MODULES*MAX_FADC_SAMPLES+1000];
  Double_t Tnsamps[MAX_TRIGGER_MODULES+1000] = {0};
  Double_t Tsamps_idx[MAX_TRIGGER_MODULES+1000] = {0};
  Double_t Ta_p[MAX_TRIGGER_MODULES+1000] = {0};
  Double_t Ta_amp_p[MAX_TRIGGER_MODULES+1000] = {0};
  Double_t Ta_amp[MAX_TRIGGER_MODULES+1000] = {0};
  Double_t TelemID[MAX_TRIGGER_MODULES+1000] = {0};
  Double_t Tcol[MAX_TRIGGER_MODULES+1000] = {0};
  Int_t Tndata = 0;
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
