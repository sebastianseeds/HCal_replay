#include "hcal.h"
const Int_t kNrows = 12;
const Int_t kNcols = 12;
const Int_t kN = kNrows*kNcols;
const Int_t kMinRows = kNrows*0+5;
const Int_t kBinsHisto = 100;

const Bool_t kCutsApplyToAllInVert = true;

const Int_t kMaxNPECount = 2500;

const Int_t gCanvSizeX = 250*1.5;
const Int_t gCanvSizeY = 150*1.5;
Float_t gBinDivider = 50.0;

const Float_t kConv = 2000./4095.; // ADC to mV

Bool_t kPlotNPEGraph = false;
Bool_t kPlotMaxADC = false;
Bool_t kPlotVs = false;
Bool_t kPlotPos = false;
Bool_t kPlotPeakVsEntry = false;

Bool_t kUseCodedRefMean = true;

// Values come from Module 2-1 Run 430
//const Double_t kRefMean   = 708.4;
//const Double_t kRefStdDev = 453.955;
//const Double_t kRefNPE    = 2.43518;
// New values from Module 2-1 PMT 190123 Run 438
//const Double_t kRefMean   = 2521.23;
//const Double_t kRefStdDev = 1292.66;
//const Double_t kRefNPE    = 3.80417;
// New values from Module 2-1 PMT 190123 Run 438
const Double_t kRefMean   = 0;
const Double_t kRefStdDev = 0;
const Double_t kRefNPE    = 4.04037;


const Double_t kQuickPed[kN] = {
  28630.0, 27564.0, 21824.0, 28146.0,
  27070.0, 26316.0, 27525.0, 24491.0,
  27198.0, 27429.0, 23833.0, 25223.0,
  23157.0, 27580.0, 27996.0, 18100.0,
  26665.0, 28228.0, 23733.0, 29258.0,
  23501.0, 24824.0, 25250.0, 30725.0
};


const Float_t kMaxTrigSlopeDeviation = 150.;
//Float_t kLimits[kN] = {
//   20000, 20000, 20000, 20000,
//   10000, 10000, 1500, 1500,
//   10000, 10000, 1500, 1500,
//   10000, 10000, 1500, 1500,
//   10000, 10000, 1500, 1500,
//   10000, 10000, 1500, 1500
//  };

Float_t kLimits[kN] = {
   20000, 20000, 20000, 20000,
   7000, 7000, 7000, 7000,
   7000, 7000, 7000, 6000,
   7000, 7000, 7000, 7000,
   7000, 7000, 7000, 7000,
   7000, 7000, 7000, 7000
  };


//Float_t kLimitsADC[kN] = {
//   2500, 2500, 2500, 2500,
//   1100, 1024, 120, 100,
//   1000, 1024,  60, 200,
//    950, 1024, 200, 120,
//    800, 1024, 120, 140,
//   1024, 1024, 180, 100
//};

Float_t kLimitsADC[kN] = {
   1220, 1220, 1220, 1220,
   430, 430, 430, 430,
   430, 430, 430, 430,
   430, 430, 430, 430,
   430, 430, 430, 430,
   430, 430, 430, 430
};


//const Float_t kMinStrig[kN] = {
//  //400., 400., 400., 400.,
//  300., 300., 300., 300.,
//  30., 30., 5., 5.,
//  15., 15., 5., 5.,
//  15., 15., 5., 5.,
//  15., 15., 5., 5.,
//  15., 15., 5., 5.
//};

const Float_t kMinStrig[kN] = {
  //400., 400., 400., 400.,
  145., 145., 145., 145.,
  24.5, 24.5, 24.5, 24.5,
  24.5, 24.5, 24.5, 24.5,
  24.5, 24.5, 24.5, 24.5,
  24.5, 24.5, 24.5, 24.5,
  24.5, 24.5, 24.5, 24.5
};


const Float_t kMinStrigPlot[kN] = {
  7., 7., 7., 7.,
  7., 7., 7., 7.,
  7., 7., 7., 7.,
  7., 7., 7., 7.,
  7., 7., 7., 7.,
  7., 7., 7., 7.
};

const Float_t kMinCut[kN] = {
  500., 500., 500., 500.,
  200., 200., 200., 200.,
  200., 200., 200., 200.,
  200., 200., 200., 200.,
  200., 200., 200., 200.,
  200., 200., 200., 200.
};

const Float_t kMinAmp[kN] = {
  2000., 2000., 2000., 2000.,
  700., 700., 700., 700.,
  700., 700., 700., 700.,
  700., 700., 700., 700.,
  700., 700., 700., 700.,
  700., 700., 700., 700.
};


//const Float_t kMinStrigPartner[kN] = {
//  200., 200., 200., 200.,
//  50., 50., 10., 10.,
//  50., 50., 10., 10.,
//  50., 50., 10., 10.,
//  50., 50., 10., 10.,
//  50., 50., 10., 10.
//};

const Float_t kMinStrigPartner[kN] = {
  200., 200., 200., 200.,
  100., 100., 100., 100.,
  100., 100., 100., 100.,
  100., 100., 100., 100.,
  100., 100., 100., 100.,
  100., 100., 100., 100.
};

const Float_t kFitLowADC[kN] = {
  0.020, 0.020, 0.020, 0.020,
  0.200, 0.200, 0.200, 0.200,
  0.200, 0.200, 0.200, 0.200,
  0.200, 0.200, 0.200, 0.200,
  0.200, 0.200, 0.200, 0.200,
  0.200, 0.500, 0.200, 0.200
};

namespace Flags {
  const UInt_t AbsFitLimit    = 1<<0;
  const UInt_t FitWithGaus = 1<<1;
};

const Float_t kFitLow[kN] = {
  0.150, 0.150, 0.150, 0.150,
  0.150, 0.150, 0.150, 0.150,
  0.150, 0.150, 0.150, 0.150,
  0.150, 0.150, 0.150, 0.150,
  0.150, 0.150, 0.150, 0.150,
  0.150, 0.150, 0.150, 0.150
};


const Float_t kGain[kN] = {
  1.00, 1.00, 1.00, 1.00,
  1.00, 1.00, 1.00, 1.00,
  1.00, 1.00, 1.00, 1.00,
  1.00, 1.00, 1.00, 1.00,
  1.00, 1.00, 1.00, 1.00,
  1.00, 1.00, 1.00, 1.00
};



const Float_t kCutsMaxPartners[kN] = {
    10000., 10000., 10000., 10000.,
    10000., 10000., 10000., 10000.,
    10000., 10000., 10000., 10000.,
    10000., 10000., 10000., 10000.,
    10000., 10000., 10000., 10000.,
    10000., 10000., 10000., 10000.
  };
const Float_t kCutsMin[kN] = {
  1000., 2750., 2750., 1750., 
  100., 100., 10., 10., 
  100., 100., 10., 10., 
  100., 100., 10., 10., 
  100., 100., 10., 10., 
  100., 100., 10., 10.
  };

const Float_t kMaxTrig[kNcols] = {
  900.e4, 900.e4, 900.e4, 900e4
};
Float_t amps[kNrows][kNcols];
Float_t peds[kNrows][kNcols];
Float_t peds_err[kNrows][kNcols];
Float_t raws[kNrows][kNcols];
Float_t smax[kNrows][kNcols];
Int_t   spos[kNrows][kNcols];
Float_t sped[kNrows][kNcols];
Float_t sadc[kNrows][kNcols];

Double_t gMaxY[kNrows][kNcols];
Double_t gXMaxY[kNrows][kNcols];
Double_t gXHMaxY[kNrows][kNcols];
Bool_t   gMinCutPassed[kNrows][kNcols];
Double_t gPed[kNrows][kNcols];


TMarker *gMark[kNrows][kNcols];

Double_t running_vals[kNrows][kNcols][kMaxNPECount];
Int_t    running_idx[kNrows][kNcols];
std::vector<Double_t> vNPE[kNrows][kNcols];
std::vector<Double_t> vNPEE[kNrows][kNcols];
std::vector<Double_t> vAvg[kNrows][kNcols];
std::vector<Double_t> vSig[kNrows][kNcols];
std::vector<Double_t> vCount[kNrows][kNcols];
std::vector<Double_t> vZero[kNrows][kNcols];
Double_t running_valsNC[kNrows][kNcols][kMaxNPECount];
Int_t    running_idxNC[kNrows][kNcols];
std::vector<Double_t> vNPENC[kNrows][kNcols];
std::vector<Double_t> vAvgNC[kNrows][kNcols];
std::vector<Double_t> vSigNC[kNrows][kNcols];
std::vector<Double_t> vCountNC[kNrows][kNcols];
std::vector<Double_t> vZeroNC[kNrows][kNcols];



Float_t trigSlopes[kNcols] = {0.106429, 0.108255, 0.105606, 0.101965 };
Float_t trigInterc[kNcols] = {43.2741, 2.2662, 14.9244, 60.0468 };
Bool_t checkTrigSlope(Int_t c)
{
  return 1;
  return (TMath::Abs(sadc[0][c]-amps[0][c]*trigSlopes[c]+trigInterc[c])<kMaxTrigSlopeDeviation);
}

void fixAmpStats()
{
  gPad->Update();
  TPaveStats *ps = (TPaveStats*)gPad->GetPrimitive("stats");
  if(ps) {
    ps->SetX1NDC(0.50);
    ps->SetY1NDC(0.50);
    ps->SetY2NDC(1.0);
    // Now find the lines with conv and sigma and remove them
    TList *list = ps->GetListOfLines();
    TText *conv = ps->GetLineWith("Constant");
    list->Remove(conv);
  }
  gPad->Modified();
  gPad->Update();
}



Bool_t checkPeakPos(Int_t pos, Int_t r, Int_t c)
{
  return 1;
  if(r==0) {
    return pos>60&&pos<90;
  }
  Int_t diff = spos[0][c]-pos;
  //if(c>1)
  //  return diff>4&&diff<9;
  //return diff>0&&diff<5;
  return diff>0&&diff<10;
}

void FindPeak(TH1F *h, Float_t &maxy, Float_t &maxx, Float_t &fwhmxlow, Float_t &fwhmxhigh, Bool_t fit = kFALSE)
{
  maxy = maxx = -1e5;
  Int_t nbins = h->GetNbinsX();
  Float_t y;
  Int_t maxbin = 0;
  Int_t b;
  for(b = 1; b <= nbins; b++) {
    y = h->GetBinContent(b);
    if(y>maxy) {
      maxy = y;
      maxx = h->GetBinLowEdge(b);
      maxbin = b;
    }
  }
  // Now find the FWHM (low end)
  b = maxbin;
  Float_t halfmaxy = maxy/2.;
  y = h->GetBinContent(b);
  while(y > halfmaxy && b-- > 0) {
    y = h->GetBinContent(b);
    fwhmxlow = h->GetBinLowEdge(b);
  }
  b = maxbin;
  y = h->GetBinContent(b);
  //std::cout << "y: " << y << ", halfmaxy: " << halfmaxy << std::endl;
  if(fit) {
    h->Fit("landau","","",h->GetBinLowEdge(1),h->GetBinLowEdge(nbins));
  } else {
    while(y > halfmaxy && b++ <= nbins) {
      y = h->GetBinContent(b);
      fwhmxhigh = h->GetBinLowEdge(b);
    }
  }
}

void FindPeak(TH1F *h, Double_t &maxy, Double_t &xmaxy, Int_t startBin = 1)
{
  xmaxy = h->GetBinLowEdge(0);
  maxy = h->GetBinContent(0)-1e4;
  Int_t nbins = h->GetNbinsX();
  Double_t y =  0;
  for(Int_t b = startBin; b <= nbins; b++) {
    y = h->GetBinContent(b);
    if(y>maxy) {
      maxy = y;
      xmaxy = h->GetBinLowEdge(b);
    }
  }
}

Bool_t CheckAmp(Int_t r, Int_t c) {
  return amps[r][c] > kMinAmp[c+r*kNcols]*.75;
}

Bool_t CheckSides(Int_t r, Int_t c)
{
  Int_t m = c+r*kNcols;
  return sadc[r][c]>kMinStrigPartner[m]||CheckAmp(r,c);
}



void FitPeak(TH1F *h, Double_t range, Double_t flow, Double_t fhigh,
    Double_t &maxy, Double_t &xmaxy, Double_t &xhmaxy, Bool_t useLandau,
    UInt_t flags = 0)
{
  Double_t temp;
  FindPeak(h,maxy,xmaxy);

    h->Scale(1.0/maxy);
    maxy = 1.0;
  Double_t dx = range*0.005;
  if(useLandau) {
    TF1 *f1 = new TF1("f1","landau");
    Double_t fit_xmin, fit_xmax;
    if( flags & Flags::AbsFitLimit) {
      fit_xmin = range*flow;
      fit_xmax = range*fhigh;
    } else {
      fit_xmin = xmaxy-range*flow;
      fit_xmax = h->GetBinLowEdge(h->GetNbinsX())*fhigh;
    }
    TFitResultPtr res = h->Fit(f1,"SQME","",fit_xmin,fit_xmax);
    if(true||res!=-1) {
      Double_t con = res->Value(0);
      Double_t mu  = res->Value(1);
      Double_t sig = res->Value(2);
      xmaxy = res->Value(1);
      xmaxy = f1->GetMaximumX(xmaxy*0.5,xmaxy*1.5);
      maxy  = TMath::Landau(xmaxy,mu,sig)*con;
      // Now find out when it drops below half maximum
      Double_t hmaxy = maxy/2.0;
      Double_t x = xmaxy;
      Double_t y = maxy;
      if(flags&Flags::FitWithGaus&&false) {
        // If we need to fit with a gaussian in the end, lets do that now
        res = h->Fit("gaus","SQ","",xmaxy-0.175*range,xmaxy+0.075*range);
        xmaxy  = res->Value(1);
        // Instead of FWHM just return the sigma of the gaussian
        xhmaxy = res->Value(2);
      } else {
        //while(x < range && y > hmaxy) {
        //  y = TMath::Landau(x,mu,sig)*con;
        //  xhmaxy = x;
        //  x += dx;
        //}
        // Cheat and return the mpv (mu) with its error
        maxy = mu;
        xmaxy = res->Error(1);
        xhmaxy = 100.*xmaxy/mu;
     }
      
    }
  } else {
    TFitResultPtr res = h->Fit("expo","QS","",xmaxy,
        h->GetBinLowEdge(h->GetNbinsX())*fhigh);
    if(res!=-1) {
      Double_t p0 = res->Value(0);
      Double_t p1 = res->Value(1);
      xmaxy = (TMath::Log(maxy)-p0)/p1;
      Double_t x = xmaxy;
      Double_t y = maxy;
      Double_t hmaxy = maxy/2.0;
      while(x < range && y > hmaxy) {
        y = TMath::Exp(p0+p1*x);
        xhmaxy = x;
        x += dx;
      }
    }
  }
  
  //TMarker *mark = new TMarker(xmaxy,1.1,23);
  //mark->SetMarkerColor(kGreen+1);
  //mark->Draw();
}

TFitResultPtr FitPeak(TH1F *h, Float_t range)
{
  Float_t temp;
  Float_t maxx = -1e4;
  Float_t maxy = -1e4;
  //Float_t width = range*0.1;
  FindPeak(h,maxy,maxx,temp,temp);

  //return h->Fit("gaus","QS","",maxx-width,maxx+width);
  return h->Fit("landau","S","",maxx-range*0.075,h->GetBinLowEdge(h->GetNbinsX())*0.9);
}

void FindMinMax(TH1F *h, Float_t mu, Float_t sig, Float_t &maxx, Float_t &xhalfmaxy, Float_t range)
{
  Float_t dx = range*0.01;
  Float_t x = dx;
  Float_t maxy = -1e4;
  maxx = -1e4;
  xhalfmaxy = -1e4;
  Float_t y1 = TMath::Landau(x,mu,sig);
  Float_t y0 = y1-1.0;
  Int_t count = 0;
  while(y0 < y1 && x < range) {
    count++;
    if(count%10==0) {
      std::cout << TString::Format(
          "Checking x=%5.2f, y0=%5.2f, y1=%5.2f",x,y0,y1) << std::endl;
    }

    y0 = y1;
    if(y0>maxy) {
      maxy = y0;
      maxx = x;
    }
    x += dx;
    y1 = TMath::Landau(x,mu,sig);
      }
}

void MakeBins(Float_t min, Float_t max, Float_t dx, Int_t &nbins, Float_t* &bins)
{
  nbins = 0;
  Float_t x = min;
  while(x < max) {
    x+=dx;
    nbins++;
  }
  bins = new Float_t[nbins+1];
  x = min;
  for(Int_t b = 0; b < nbins; b++) {
    bins[b] = x;
    x+=dx;
  }
  bins[nbins] = max; // Define the overflow bin
}

void MakeIntBins(Int_t min, Int_t max, Bool_t inclusive, Int_t &nbins,
    Float_t* &bins)
{
  Float_t dx = Float_t( (max-min)+inclusive )/gBinDivider;
  if(dx < 1.25 ) {
    dx = 1.0;
  } else {
    Int_t base = Int_t(dx);
    Float_t ddx = (dx-base);
    dx = base;
    if(ddx >= 0.75) {
      dx += 1.0;
    } else if ( ddx >= 0.25 ) {
      dx += 0.5;
    }
  }
  MakeBins(min,max+inclusive,dx,nbins,bins);
}


void ProcessRunningNPE(Int_t r, Int_t c, Bool_t nc = false)
{
  Int_t count = running_idx[r][c];
  running_idx[r][c] = 0;
  if(count < 2)
    return; // Do nothing!

  Double_t avg = 0.;
  Double_t sig2 = 0.;
  Double_t npe = 0;
  for(Int_t i = 0; i < count; i++) {
    avg += running_vals[r][c][i];
  }
  avg /= Double_t(count);
  for(Int_t i = 0; i < count; i++) {
    sig2 += TMath::Power(running_vals[r][c][i] - avg,2.0);
    running_vals[r][c][i] = 0.0;
  }
  sig2 /= Double_t(count-1); 
  Double_t sig = TMath::Sqrt(sig2);
  vNPE[r][c].push_back(avg*avg/sig2);
  vNPEE[r][c].push_back(1./TMath::Sqrt(count));
  vAvg[r][c].push_back(avg);
  vSig[r][c].push_back(sig);
  vCount[r][c].push_back(vNPE[r][c].size());
  vZero[r][c].push_back(0.0);
}

void ProcessRunningNPENC(Int_t r, Int_t c, Bool_t nc = false)
{
  Int_t count = running_idxNC[r][c];
  running_idxNC[r][c] = 0;
  if(count < 2)
    return; // Do nothing!

  Double_t avg = 0.;
  Double_t sig2 = 0.;
  Double_t npe = 0;
  for(Int_t i = 0; i < count; i++) {
    avg += running_valsNC[r][c][i];
  }
  avg /= Double_t(count);
  for(Int_t i = 0; i < count; i++) {
    sig2 += TMath::Power(running_valsNC[r][c][i] - avg,2.0);
    running_valsNC[r][c][i] = 0.0;
  }
  sig2 /= Double_t(count-1); 
  vNPENC[r][c].push_back(avg*avg/sig2);
  vAvgNC[r][c].push_back(avg);
  vSigNC[r][c].push_back(TMath::Sqrt(sig2));
  vCountNC[r][c].push_back(vNPENC[r][c].size());
  vZeroNC[r][c].push_back(0.0);
}

void PrintFitInfo(Int_t r, Int_t c)
{
  std::cout << TString::Format("[[%02d,%02d]] %7.2f +/- %7.2f (%7.2f%%)",r+1,c+1,gMaxY[r][c],gXMaxY[r][c],gXHMaxY[r][c]) << std::endl;
}

void printMaxY(Int_t r, Int_t c)
{
  //std::cout << TString::Format("[[%02d,%02d]] %7.2f  %7.2f  %7.2f %7.2f mV",r+1,c+1,gMaxY[r][c],gXMaxY[r][c],gXHMaxY[r][c],kConv*gXMaxY[r][c]/kGain[c+r*kNcols]) << std::endl;
  std::cout << TString::Format("[[%02d,%02d]] %7.2f  %7.2f mV  %7.2f mV",r+1,c+1,gMaxY[r][c],gXMaxY[r][c],gXHMaxY[r][c]) << std::endl;
}
TH1F *histos[kNrows][kNcols];
TH1F *histosNoCut[kNrows][kNcols];
TH1F *histosPeak[kNrows][kNcols];
TH1F *histosPeakNoCut[kNrows][kNcols];
TH1F *histosPos[kNrows][kNcols];
TH2F *histosVs[kNrows][kNcols];

void FillHistos(Int_t r, Int_t c)
{
  if(sadc[r][c]>kMinStrigPlot[c+r*kNcols]) {
    histosPeak[r][c]->Fill(sadc[r][c]);
  }
  histosPos[r][c]->Fill(spos[0][c]-(r>0?spos[r][c]:0));
  histos[r][c]->Fill(amps[r][c]);
  histosVs[r][c]->Fill(amps[r][c],sadc[r][c]);
  running_vals[r][c][running_idx[r][c]] = amps[r][c];
  running_idx[r][c]++;

  if(running_idx[r][c]==kMaxNPECount) {
    ProcessRunningNPE(r,c);
  }
}

void FindPedestals(TChain *T, TCanvas *canvas)
{
  TH1F *h;
  TH1F *h2;
  Int_t m = 0;
  for(Int_t r = 0; r < kNrows; r++) {
    for(Int_t c = 0; c < kNcols; c++) { 
      m++;
      canvas->cd(m);
      TString var = Form("raw%d",m);
      TString hname = Form("h%02d%02d",r+1,c+1);
      T->Draw(Form("%s>>%s",var.Data(),hname.Data()),
          //Form("%s>0&&%s<30e3&&smax%d<4095",var.Data(),var.Data(),m));
          Form("%s>0&&smax%d<4095",var.Data(),m));
      h = (TH1F*)gDirectory->Get(hname.Data());
      TSpectrum *spec = new TSpectrum(5,1);
      spec->Search(h);
      if(spec->GetNPeaks()>0){
         Double_t ped = spec->GetPositionX()[0];
         //T->Draw(Form("(%s-%g)>>h2%02d%02d(100,0,10000)",
         //    var.Data(),ped,r+1,c+1),
         //    Form("%s-%g>100",var.Data(),ped));
         Double_t pm = 1000.;
         T->Draw(Form("%s>>h2%02d%02d(100,%g,%g)",var.Data(),r+1,c+1,
           ped-pm,ped+pm));
         h2 = (TH1F*)gDirectory->Get(Form("h2%02d%02d",r+1,c+1));
         TFitResultPtr res = h2->Fit("gaus","QS");
         ped = res->Value(1);
         gPed[r][c] = ped;
         std::cout << Form("Pedestal [[%02d,%02d]] %7.1f",r+1,c+1,ped) << std::endl;
      }
    }
  }
}

Int_t analyze_cosmics(Int_t run = 457)
{
  //hcalgui::SetupGUI();
  gStyle->SetLabelSize(0.05,"XY");
  gStyle->SetTitleFontSize(0.09);
  gStyle->SetOptFit(111);
  gStyle->SetOptStat(10);
  gStyle->SetTitleAlign(33);

  TChain *T = new TChain("TInt");
  if( T->Add(TString::Format("%s/cosmics/hcal_%d.root",
       getenv("HCAL_ROOTFILES"),run) ) <= 0) {
     std::cout << "ROOT file now found!" << std::endl;
     return -1;
   }

  TString canvasName = TString::Format("canvas_Run%d",run);
  TCanvas *canvas = new TCanvas(canvasName.Data(),
    canvasName.Data(),kNcols*gCanvSizeX,kNrows*gCanvSizeY);
  canvas->Divide(kNcols,kNrows);
  // First, find the pedestal
  FindPedestals(T,canvas);

  Int_t m = 0;
  Int_t lposmin,lposmax;
  Int_t lnbinspeak;
  Float_t *lbinspeak;
  for(Int_t r = 0; r < kNrows; r++) {
    for(Int_t c = 0; c < kNcols; c++) {
      m++;
      T->SetBranchAddress(TString::Format("c%d",m),&(amps[r][c]));
      //T->SetBranchAddress(TString::Format("raw%d",m),&(amps[r][c]));
      T->SetBranchAddress(TString::Format("ped%d",m),&(peds[r][c]));
      T->SetBranchAddress(TString::Format("ped_err%d",m),&(peds_err[r][c]));
      T->SetBranchAddress(TString::Format("raw%d",m),&(raws[r][c]));
      T->SetBranchAddress(TString::Format("smax%d",m),&(smax[r][c]));
      T->SetBranchAddress(TString::Format("spos%d",m),&(spos[r][c]));
      T->SetBranchAddress(TString::Format("sped%d",m),&(sped[r][c]));
      TString name = TString::Format("hR%02dC%02d",r+1,c+1);
      histos[r][c] = new TH1F(name.Data(),
        TString::Format("ADC-Ped %02d-%02d",r+1,c+1),kBinsHisto,0,kLimits[m-1]);
      lnbinspeak = kLimitsADC[m-1] <= 150 ? kLimitsADC[m-1] : 150;
      histosNoCut[r][c] = new TH1F(TString::Format("hNocutR%02dC%02d",r+1,c+1),
        TString::Format("ADC-Ped %02d-%02d",r+1,c+1),kBinsHisto,0,kLimits[m-1]);
      MakeIntBins(0,kLimitsADC[m-1],true,lnbinspeak,lbinspeak);
      histosPeak[r][c] = new TH1F(TString::Format("hPeakR%02d%02d",r+1,c+1),
      //  TString::Format("Peak %02d-%02d",r+1,c+1),lnbinspeak,0.0,kLimitsADC[m-1]);
        TString::Format("Peak %02d-%02d",r+1,c+1),lnbinspeak,lbinspeak);
      lposmin=r>0?0:55;
      lposmax=r>0?10:75;
      histosPos[r][c] = new TH1F(TString::Format("hPosR%02d%02d",r+1,c+1),
        TString::Format("PeakPos %02d-%02d",r+1,c+1),lposmax-lposmin,lposmin,lposmax);
      histosPeakNoCut[r][c] = new TH1F(
        TString::Format("hPeakNocutR%02dC%02d",r+1,c+1),
        TString::Format("Peak %02d-%02d",r+1,c+1),lnbinspeak,lbinspeak);

      histosVs[r][c] = new TH2F(TString::Format("hVsR%02dc%02d",r+1,c+1),
        TString::Format("Peak vs ADC-Ped %02d-%02d",r+1,c+1),
          500,0,kLimits[m-1],500,0,kLimitsADC[m-1]);
      histos[r][c]->SetLineColor(kBlue+1);
      histos[r][c]->SetLineWidth(2);
      histosPeak[r][c]->SetLineColor(kBlue+1);
      histosPeak[r][c]->SetLineWidth(2);
      histosPos[r][c]->SetLineColor(kBlue+1);
      histosPos[r][c]->SetLineWidth(2);
      histosNoCut[r][c]->SetLineColor(kRed+1);
      histosPeakNoCut[r][c]->SetLineColor(kRed+1);
      histosVs[r][c]->SetMarkerColor(kBlue+1);
      //histosNoCut[r][c]->SetLineWidth(2);
      for(Int_t k = 0; k < kMaxNPECount; k++) {
        running_vals[r][c][k] = 0.0;
      }
      running_idx[r][c] = 0;
    }
  }

  Int_t entries = T->GetEntries();
  if(run==430) {
    entries-=2000; // Just because last minute I was rotating PMTs
  }
  if(entries>148000 && (run == 442 | run==438))
    entries=148000;
  Int_t cc,mm;
  Int_t colCount[kNcols];
  Bool_t trig[kNcols];
  Int_t trigCount = 0;
  Int_t trigCol = -1;
  Bool_t no_sat = true;
  Int_t mp, mn;
  for(Int_t entry = 0; entry < entries; entry++) {
    T->GetEntry(entry);

    m = 0;
    for(Int_t r = 0; r < kNrows; r++) {
      for(Int_t c = 0; c < kNcols; c++) {
        m = c+r*kNcols;
        amps[r][c] = raws[r][c] - gPed[r][c];
        sadc[r][c] = (smax[r][c]-sped[r][c])*kGain[m]*kConv;
        if(sadc[r][c]>kMinStrigPlot[m])
          histosPeakNoCut[r][c]->Fill(sadc[r][c]);
        //histosPeakNoCut[r][c]->Fill(spos[r][c]);
        if(amps[r][c]>500)
        histosNoCut[r][c]->Fill(amps[r][c]);
	//if(sadc[r][c]>kMinStrigPartner[c+r*kNcols]) {
        gMinCutPassed[r][c] = (amps[r][c]>kMinCut[m]);	
        if(sadc[r][c]>kMinStrig[c+r*kNcols]) {
          running_valsNC[r][c][running_idx[r][c]] = amps[r][c];
          running_idxNC[r][c]++;
          if(running_idxNC[r][c]==kMaxNPECount) {
            ProcessRunningNPENC(r,c);
          }
        }
      }
    }
    // First, determine which module on top row triggered
    trigCount = 0;
    no_sat = true;
    trigCol = -1;
    for(Int_t c = 0; c < kNcols; c++) {
      //trig[c] = (amps[0][c]>kCutsMaxPartners[c]?true:false);
      trig[c] = (sadc[0][c]>kMinStrig[c]?true:false);
      if(smax[0][c]>4095)
        no_sat = false;
      if(trig[c]) {
        trigCount++;
        trigCol = c;
      }
    }
    Bool_t passedCuts[kNrows];
    Bool_t goodTrig = (trigCount==1);
    if(false&&goodTrig) {
      // Quickly check the other columns to ensure they failed the minimum cut
      for(Int_t c = 0; c < kNcols; c++) {
        if(c != trigCol) {
          for(Int_t r = 0; r < kNrows; r++) {
            if(gMinCutPassed[r][c])
              goodTrig = false;
          }
        }
      }
    }
    if(no_sat&&goodTrig&&checkTrigSlope(trigCol)) {
      Int_t passedRows = 0;
      Int_t c = trigCol;
      for(Int_t r = 0; r < kNrows; r++) {
        m = c + r*kNcols;
        passedCuts[r] = (sadc[r][c]>kMinStrig[m]&&smax[r][c]<4095);
        if(!trig[c]||sadc[0][c]>kMaxTrig[c]||!checkPeakPos(spos[r][c],r,c))
          passedCuts[r] = false;
        if(r>0) {
          // Check that both sides did not trigger
          if( (c>0 && CheckSides(r,c-1)) ||
              (c+1<kNcols&&CheckSides(r,c+1)) ) {
            passedCuts[r] = false;
          }
          // Check that the top module triggered
          if(!CheckAmp(r-1,c)) {
            passedCuts[r]=false;
          }
        }
        // Finally, check that the bottom module triggered
        if(r+1<kNrows && !CheckAmp(r+1,c)) {
          passedCuts[r] = false;
        }
        if(passedCuts[r]&&r<kMinRows)
          passedRows++;
        if(!kCutsApplyToAllInVert&&passedCuts[r]) {
          FillHistos(r,c);
        }
      }
      if(kCutsApplyToAllInVert&&passedRows>=kMinRows) {
        for(Int_t r = 0; r < kNrows; r++) {
          FillHistos(r,c);
          //if(sadc[r][c]>kMinStrigPlot[c+r*kNcols]) {
          //  histosPeak[r][c]->Fill(sadc[r][c]);
          //}
          //histosPos[r][c]->Fill(spos[0][c]-(r>0?spos[r][c]:0));
          //histos[r][c]->Fill(amps[r][c]);
          //histosVs[r][c]->Fill(amps[r][c],sadc[r][c]);
          //running_vals[r][c][running_idx[r][c]] = amps[r][c];
          //running_idx[r][c]++;
          //if(running_idx[r][c]==kMaxNPECount) {
          //  ProcessRunningNPE(r,c);
          //}
        }
      }
    }
  }


  m = 0;
  fstream outSummary(TString::Format("summary/cosmics_%d.dat",run),std::ios::out);
  Double_t ref_mean, ref_stdev, ref_npe;
  if(kUseCodedRefMean) {
    ref_mean = kRefMean;
    ref_stdev = kRefStdDev;
    ref_npe = kRefNPE;
  } else {
    ref_mean = histos[1][0]->GetMean();
    ref_stdev = histos[1][0]->GetRMS();
    ref_npe = (ref_mean*ref_mean)/(ref_stdev*ref_stdev);
  }
  Double_t mean,stdev,ncount,npe,npeNC,npeHisto;
  for(Int_t r = 0; r < kNrows; r++) {
    for(Int_t c = 0; c < kNcols; c++) {
      m++;
      canvas->cd(m);
      // Make the graphs
      if(vNPE[r][c].size() == 0 && running_idx[r][c] > 2) {
        ProcessRunningNPE(r,c);
      }
      if(vNPENC[r][c].size() == 0 && running_idxNC[r][c]>2) {
        ProcessRunningNPENC(r,c);
      }
      TGraphErrors *graph = new TGraphErrors(vCount[r][c].size(),
        vCount[r][c].data(),vNPE[r][c].data(),
        vZero[r][c].data(),vNPEE[r][c].data());
      //TGraph *graph = new TGraph(vCount[r][c].size(), vCount[r][c].data(),
      //  vNPE[r][c].data());
      //graph->SetMarkerColor(kBlue+1);
      graph->SetMarkerStyle(20);
      TGraph *graphNC = new TGraph(vCountNC[r][c].size(), vCountNC[r][c].data(),
        vNPENC[r][c].data());
      graphNC->SetMarkerColor(kBlue+1);
      graphNC->SetMarkerStyle(20);

      if(kPlotPeakVsEntry&&r==0) {
        Int_t mm = c+r*kNcols+1;
        TH2F *h = new TH2F(TString::Format("hTest%02d",mm),"",1000,0,entries,
            1000, 390, 1000);
        T->Draw(TString::Format("smax%d-sped%d:Entry$>>hTest%02d",mm,mm,mm),
            TString::Format("smax%d<1000&&smax%d-sped%d>390",mm,mm,mm),"COLZ");
        h->Fit("pol1");
      } else
      if(kPlotNPEGraph) {
        graph->Draw("AP");
        //graphNC->Draw("AP");
        gPad->Update();
      } else if(kPlotMaxADC) {
        histosPeak[r][c]->Draw();
        // Find the peak and fit to that
        if(true||r>0) {
          gPad->SetLogy(kFALSE);
          //Double_t lmaxy,lxmaxy,lxhmaxy;
          FitPeak(histosPeak[r][c],kLimitsADC[c+r*kNcols],
              //(r==0?0.02:0.075),0.9,lmaxy,lxmaxy,lxhmaxy,r>0);
              kFitLowADC[c+r*kNcols],0.9,gMaxY[r][c],gXMaxY[r][c],gXHMaxY[r][c],r>0);
          //gMark[r][c] = new TMarker(gXMaxY[r][c],1.1,23);
          //gMark[r][c]->SetMarkerColor(kGreen+1);
          //gMark[r][c]->Draw();
          //if(rpeak!=-1) {
            //std::cout << TString::Format("[[%02d,%02d]] = %5.2f  %5.2f  %5.2f",r+1,c+1,rpeak->Value(0),rpeak->Value(1), rpeak->Value(2)) << std::endl;
            //Float_t lmaxx,lxhalfmaxy=-1e4;
            //FindMinMax(rpeak->Value(0),rpeak->Value(1),lmaxx,lxhalfmaxy,kLimitsADC[c+r*kNcols]);
            //std::cout << TString::Format("[[%02d,%02d]] mu= %5.2f  sig= %5.2f",r+1,c+1,rpeak->Value(0),rpeak->Value(1)) << std::endl;
            //std::cout << TString::Format("[[%02d,%02d]] %7.2f  %7.2f  %7.2f %7.2fmV",r+1,c+1,gMaxY[r][c],gXMaxY[r][c],gXHMaxY[r][c],kConv*gXMaxY[r][c]/kGain[c+r*kNcols]) << std::endl;
          //}
        } else {
          gPad->SetLogy(kTRUE*0);
          Float_t maxpeaky,maxpeakx,fwhmxlow,fwhmxhigh = -1;
          FindPeak(histosPeak[r][c],maxpeaky,maxpeakx,fwhmxlow,fwhmxhigh,true);
          std::cout <<   TString::Format("[[%02d,%02d]] = %5.2f  %5.2f %5.2f",r+1,c+1,maxpeaky,maxpeakx,fwhmxhigh) << std::endl;
        }
        fixStats();
        //histosPeak[r][c]->Fit("gaus","","",40,80);
        //Double_t tmp1,tmp2;
        //FindPeak(histosPeakNoCut[r][c],tmp1,tmp2,10);
        //histosPeakNoCut[r][c]->Scale(1.0/tmp1);
        //histosPeakNoCut[r][c]->Draw("SAME");
        //histosPeakNoCut[r][c]->Draw();
        //gPad->SetLogy(kTRUE);
      } else if (kPlotVs) {
        histosVs[r][c]->Draw("COLZ");
        fixStats();
        histosVs[r][c]->Fit("pol1");
      } else if (kPlotPos) {
        histosPos[r][c]->Draw();
      } else {
        FitPeak(histos[r][c],kLimits[c+r*kNcols],
          kFitLow[c+r*kNcols],0.9,gMaxY[r][c],gXMaxY[r][c],gXHMaxY[r][c],
          true,Flags::AbsFitLimit|Flags::FitWithGaus);
        //histos[r][c]->SetTitle(Form("Amplitude %02d-%02d [ %7.2f +/- %7.2f ]",r+1,c+1,gMaxY[r][c],gXMaxY[r][c]));
        histos[r][c]->Draw("");

        fixAmpStats();
        //histosNoCut[r][c]->Draw("SAME");
        //histosNoCut[r][c]->Draw("");
        //gPad->SetLogy(kTRUE);
        //fixStats();
      }
      mean = histos[r][c]->GetMean();
      stdev = histos[r][c]->GetRMS();
      ncount = histos[r][c]->GetEntries();
      npeHisto = (mean*mean)/(stdev*stdev);
      TFitResultPtr res = graph->Fit("pol0","QS");
      TFitResultPtr resNC = graphNC->Fit("pol0","QS");
      npe = npeNC = -1;
      if(res==0)
        npe = res->Value(0);
      if(resNC==0)
        npeNC = resNC->Value(0);

      /*std::cout << "[" << r+1 << ", " << c+1 << "] npe = " << npe
        << " (mean: " << mean << ", stdev: " << stdev << ", n: " << ncount << ", qe: " << (npe/ref_npe) << ")" << std::endl;

      */
      //printf("[%02d,%02d] %7.5f %7.5f %7.5f\n",r+1,c+1,npe,npeNC,npeHisto);
      outSummary << TString::Format("[%02d,%02d] %7.5f %7.5f %7.5f",r+1,c+1,npe,npeNC,npeHisto) << std::endl;
      //std::cout << "[ " << r+1 << ", " << c+1 << " ] " << npe << " " << npeNC << " " << npeHisto << std::endl;
      //histosNoCut[r][c]->Draw();
      //gPad->SetLogy(true);
    }
  }
  canvas->SaveAs(TString::Format("plots/cosmics_%d.png",run));
  outSummary.close();


  std::cout << "19-pin PMTs:" << std::endl;
  for(Int_t r = 1; r < kNrows; r++) {
    for(Int_t c = 0; c < 2; c++) {
      PrintFitInfo(r,c);
    }
  }
  std::cout << "20-pin PMTs:" << std::endl;
  for(Int_t r = 1; r < kNrows; r++) {
    for(Int_t c = 2; c < kNcols; c++) {
      PrintFitInfo(r,c);
    }
  }

  return 0;
}
