#include <TH2F.h>
#include <TChain.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <TSystem.h>
#include "hcal.h"

const Int_t kNrows = 12;
const Int_t kNcols = 12;

const Int_t kNumModules = kNrows*kNcols;
const Int_t DISP_MIN_SAMPLE = 85*0;
const Int_t DISP_MAX_SAMPLE = 135*0+200*0+150*0+20;
const Int_t DISP_FADC_SAMPLES = (DISP_MAX_SAMPLE-DISP_MIN_SAMPLE);
const Int_t numSamples = 50;

const Int_t minADC = 0;
const Int_t maxADC = 4000;
const Int_t kCanvSize = 100;

std::string user_input;

Int_t gCurrentEntry = -1;

TChain *T = 0;
Int_t foundModules = 0;

Double_t nhit = 0;
TH1F *histos[kNrows][kNcols];
Bool_t gSaturated[kNrows][kNcols];

//TCanvas *C1 = new TCanvas();
//C1->cd();

TH1F* MakeHisto(Int_t row, Int_t col, Int_t bins)
{
  TH1F *h = new TH1F(TString::Format("h%02d%02d",row,col),
		     TString::Format("%d-%d",row+1,col+1),bins,DISP_MIN_SAMPLE,DISP_MAX_SAMPLE);
  h->SetStats(0);
  h->SetLineWidth(2);
  return h;
}

void displayEvent(string file, Int_t entry = 1)
{
  
  //TFile *infile = new TFile(file);
  //TTree *D = (TTree*)infile->Get("T");  
  
  gCurrentEntry = entry;
  
  T->GetEntry(gCurrentEntry);
  std::cout << "Displaying event " << gCurrentEntry << std::endl;
  
  Int_t r,c,idx,n,sub;
  // Clear old histograms, just in case modules are not in the tree
  for(r = 0; r < kNrows; r++) {
    for(c = 0; c < kNcols; c++) {
      histos[r][c]->Reset("ICES M");
      gSaturated[r][c] = false;
    }
  }

  Float_t peak[kNrows][kNcols];
  Double_t adc[kNrows][kNcols];
  Double_t tdc[kNrows][kNcols];
  //Initialize all arrays to zero
  for(r  = 0; r < kNrows; r++) {
    for(c  = 0; c < kNcols; c++) {
      peak[r][c] = 0.0;
      adc[r][c] = 0.0;
      tdc[r][c] = 0.0;
    }
  }
  
  //
  for(Int_t m = 0; m < hcalt::ndata; m++) {
    r = hcalt::row[m]-1;
    c = hcalt::col[m]-1;
    if(r < 0 || c < 0) {
      std::cerr << "Why is row negative? Or col?" << std::endl;
      continue;
    }
    
    std::cout << "row value: " << r << endl;
    std::cout << "column value: " << c << endl;
    
    if(r>= kNrows || c >= kNcols)
      continue;
    idx = hcalt::samps_idx[m];
    n = hcalt::nsamps[m];
    adc[r][c] = hcalt::a[m];
    tdc[r][c] = hcalt::tdc[m];
    //std::cout << "n=" << hcalt::nsamps[m] << std::endl;
    bool displayed = false;
    for(Int_t s = DISP_MIN_SAMPLE; s < DISP_MAX_SAMPLE && s < n; s++) {
      displayed = true;
      histos[r][c]->SetBinContent(s+1-DISP_MIN_SAMPLE,hcalt::samps[idx+s]);
      if(peak[r][c]<hcalt::samps[idx+s])
	peak[r][c]=hcalt::samps[idx+s];
      if(peak[r][c]>4095) {
	gSaturated[r][c] = true;
      }
      //std::cout << "setting bin content: [" << r+1 << ", " << c+1 << ", " << s << "] = " << hcalt::samps[idx+s] << std::endl;
    }
    if(!displayed) {
      std::cerr << "Skipping empty module: " << m << std::endl;
      for(Int_t s = 0;  s < DISP_FADC_SAMPLES; s++) {
	histos[r][c]->SetBinContent(s+1,-404);
      }
    }
  }
  
  for(r = 0; r < kNrows; r++) {
    for(c = 0; c < kNcols; c++) {
      sub = r/6;
      //subCanv[sub]->cd(c*kNrows + r + 1);
      subCanv[sub]->cd((r%6)*kNcols + c + 1);
      histos[r][c]->SetTitle(TString::Format("%d-%d (ADC=%g,TDC=%g)",r+1,c+1,adc[r][c],tdc[r][c]));
      if(gSaturated[r][c])
	histos[r][c]->SetLineColor(kRed+1);
      else
	histos[r][c]->SetLineColor(kBlue+1);
      if(tdc[r][c]!=0)
	histos[r][c]->SetLineColor(kGreen+1);
      
      histos[r][c]->Draw();
      gPad->Update();
      std::cout << " [" << r << ", " << c << "]=" << peak[r][c];
    }
  }
  std::cout << std::endl;
  //gSystem->mkdir("images/",kTRUE);
  //std::cerr << "Saving canvas!" << std::endl;
  //canvas->SaveAs("images/display_hcal.png");
  //  canvVector[i]->SaveAs(TString::Format("images/canvas_%d.png",int(i)));
  
}

Int_t display(Int_t run = 989, Int_t event = -1)
{
  hcalgui::SetupGUI();
  gStyle->SetLabelSize(0.05,"XY");
  gStyle->SetTitleFontSize(0.08);
  
  if(!T) { 
    T = new TChain("T");
    T->Add(TString::Format("fadc_f1tdc_%d.root",run));
    T->SetBranchStatus("*",0);
    T->SetBranchStatus("sbs.hcal.*",1);
    T->SetBranchAddress("sbs.hcal.nsamps",hcalt::nsamps);
    T->SetBranchAddress("sbs.hcal.a",hcalt::a);
    T->SetBranchAddress("sbs.hcal.tdc",hcalt::tdc);
    //T->SetBranchAddress("sbs.hcal.ledbit",&hcalt::ledbit);
    //T->SetBranchAddress("sbs.hcal.ledcount",&hcalt::ledcount);
    T->SetBranchAddress("sbs.hcal.samps",hcalt::samps);
    T->SetBranchAddress("sbs.hcal.samps_idx",hcalt::samps_idx);
    T->SetBranchAddress("sbs.hcal.row",hcalt::row);
    T->SetBranchAddress("sbs.hcal.col",hcalt::col);
    T->SetBranchStatus("Ndata.sbs.hcal.row",1);
    T->SetBranchAddress("Ndata.sbs.hcal.row",&hcalt::ndata);
    std::cerr << "Opened up tree with nentries=" << T->GetEntries() << std::endl;
    for(Int_t r = 0; r < kNrows; r++) {
      for(Int_t c = 0; c < kNcols; c++) {
	histos[r][c] = MakeHisto(r,c,DISP_FADC_SAMPLES);
	gSaturated[r][c] = false;
      }
    }
    //canvas = new TCanvas("canvas","canvas",kCanvSize*kNcols,kCanvSize*kNrows);
    //canvas->Divide(kNcols,kNrows);
  }
  //return 0;
  gCurrentEntry = event;
  while( user_input != "q" ) {
    if(is_number(user_input)) {
      gCurrentEntry = std::stoi(user_input);
    } else {
      gCurrentEntry++;
    }
    //dis1(event);
    displayEvent(gCurrentEntry);
    std::cout << "Display options: <enter> == next event, or q to stop." << std::endl;
    //std::cin >> user_input;
    getline(std::cin,user_input);
  }
  return 0;
}


void singleEvent(Int_t run = 989, Int_t event = 1)
{
  
  display(run,event);
  
}
