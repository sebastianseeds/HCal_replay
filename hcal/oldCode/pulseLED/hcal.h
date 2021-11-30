#ifndef HCAL_H
#define HCAL_H
#include <TPaveStats.h>

const Int_t MAX_FADC_SAMPLES = 250;
const Int_t MAX_HCAL_MODULES = 288;
const Int_t MAX_HCAL_TDC_MODULES = 288;
namespace hcalt {
  Double_t samps[MAX_HCAL_MODULES*MAX_FADC_SAMPLES+1000];
  Double_t nsamps[MAX_HCAL_MODULES+1000] = {0};
  Double_t row[MAX_HCAL_MODULES+1000] = {0};
  Double_t col[MAX_HCAL_MODULES+1000] = {0};
  Double_t samps_idx[MAX_HCAL_MODULES+1000] = {0};
  Double_t a[MAX_HCAL_MODULES+1000] = {0};
  Int_t ndata = 0;
  Double_t ledbit = -1;
  Double_t ledcount = 0;
  Double_t tdc[MAX_HCAL_TDC_MODULES+100];

};

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
