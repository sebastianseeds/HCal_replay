
#include <TH2F.h>
#include <TChain.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <TSystem.h>
#include "hcal.h"

const Int_t kNrows = 24;
const Int_t kNcols = 12;

const Int_t kNumModules = kNrows*kNcols;
const Int_t DISP_MIN_SAMPLE = 0;
const Int_t DISP_MAX_SAMPLE = 50;
//const Int_t DISP_MAX_SAMPLE = 30;
//const Int_t DISP_FADC_SAMPLES = 200;
const Int_t DISP_FADC_SAMPLES = (DISP_MAX_SAMPLE-DISP_MIN_SAMPLE);
const Int_t numSamples = 50;

const Int_t minADC = 0;
const Int_t maxADC = 4000;
const Int_t kCanvSize = 100;

std::string user_input;


Int_t gCurrentEntry = -1;

TChain *T = 0;
Int_t foundModules = 0;
TCanvas *canvas = 0;

TCanvas *subCanv[4];

void clicked_displayEntryButton();
void clicked_displayNextButton();
namespace hcalgui {
  TGMainFrame *main = 0;
  TGHorizontalFrame *frame1 = 0;
  TGTab *fTab;
  TGLayoutHints *fL3;
  TGCompositeFrame *tf;
  TGTextButton *exitButton;
  TGTextButton *displayEntryButton;
  TGTextButton *displayNextButton;
  TGNumberEntry *entryInput;
  TGLabel *ledLabel;

  TRootEmbeddedCanvas *canv[4];

  TGCompositeFrame* AddTabSub(Int_t sub) {
    tf = fTab->AddTab(Form("HCAL Sub%d",sub+1));

    TGCompositeFrame *fF5 = new TGCompositeFrame(tf, (12+1)*kCanvSize,(6+1)*kCanvSize , kHorizontalFrame);
    TGLayoutHints *fL4 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX |
        kLHintsExpandY, 5, 5, 5, 5);
    TRootEmbeddedCanvas *fEc1 = new TRootEmbeddedCanvas(Form("hcalSubCanv%d",sub), fF5, 6*kCanvSize,8*kCanvSize);
    //TRootEmbeddedCanvas *fEc1 = new TRootEmbeddedCanvas(0, fF5, 600, 600);
    //Int_t wid = fEc1->GetCanvasWindowId();
    //subCanv[sub] = new TCanvas(Form("subCanv%d",sub),10,10,wid);
    //subCanv[sub]->Divide(12,6);
    //fEc1->AdoptCanvas(subCanv[sub]);
    canv[sub] = fEc1;
    fF5->AddFrame(fEc1,fL4);
    tf->AddFrame(fF5,fL4);
    return tf;
  }

  void SetupGUI() {
    if(!main) {
      main = new TGMainFrame(gClient->GetRoot(), 1000, 900);
      frame1 = new TGHorizontalFrame(main, 150, 20, kFixedWidth);
      ledLabel = new TGLabel(frame1,"LED Bit:    , Count:      ");
      displayEntryButton = new TGTextButton(frame1,"&Display Entry","clicked_displayEntryButton()");
      entryInput = new TGNumberEntry(frame1,0,5,-1,TGNumberFormat::kNESInteger);
      displayNextButton = new TGTextButton(frame1,"&Next Entry","clicked_displayNextButton()");
      exitButton = new TGTextButton(frame1, "&Exit", 
          "gApplication->Terminate(0)");
      TGLayoutHints *frame1LH = new TGLayoutHints(kLHintsTop|kLHintsLeft|
          kLHintsExpandX,2,2,2,2);
      frame1->AddFrame(ledLabel,frame1LH);
      frame1->AddFrame(displayEntryButton,frame1LH);
      frame1->AddFrame(entryInput,frame1LH);
      frame1->AddFrame(displayNextButton,frame1LH);
      frame1->AddFrame(exitButton,frame1LH);
      frame1->Resize(800, displayNextButton->GetDefaultHeight());
      main->AddFrame(frame1, new TGLayoutHints(kLHintsBottom | kLHintsRight, 2, 2, 5, 1));

      // Create the tab widget
      fTab = new TGTab(main, 300, 300);
      fL3 = new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 5);

      // Create Tab1 (HCAL Sub1)
      for(Int_t i = 0; i < 4; i++) {
        tf = AddTabSub(i);
      }
      main->AddFrame(fTab, new TGLayoutHints(kLHintsBottom | kLHintsExpandX |
                                          kLHintsExpandY, 2, 2, 5, 1));
      main->MapSubwindows();
      main->Resize();   // resize to default size
      main->MapWindow();

      for(Int_t i = 0; i < 4; i++) {
        subCanv[i] = canv[i]->GetCanvas();
	if( kNrows<12 || kNcols<12) {
	  //subCanv[i]->Divide(kNrows,kNcols,0.001,0.001);
	  subCanv[i]->Divide(kNcols,kNrows,0.001,0.001);
        } else {
	  subCanv[i]->Divide(12,6,0.001,0.001);
        }
      }
    }
  }
};


Double_t nhit = 0;
TH1F *histos[kNrows][kNcols];
Bool_t gSaturated[kNrows][kNcols];

TH1F* MakeHisto(Int_t row, Int_t col, Int_t bins)
{
  TH1F *h = new TH1F(TString::Format("h%02d%02d",row,col),
      TString::Format("%d-%d",row+1,col+1),bins,DISP_MIN_SAMPLE,DISP_MAX_SAMPLE);
  h->SetStats(0);
  h->SetLineWidth(2);
  return h;
}


bool is_number(const std::string& mystring)
{
  std::string::const_iterator it = mystring.begin();
  while (it != mystring.end() && std::isdigit(*it)) ++it;
  return !mystring.empty() && it == mystring.end();
}

void displayEvent(Int_t entry = -1)
{
  if(entry == -1) {
    gCurrentEntry++;
  } else {
    gCurrentEntry = entry;
  }

  if(gCurrentEntry<0) {
    gCurrentEntry = 0;
  }

  T->GetEntry(gCurrentEntry);
  std::cout << "Displaying event " << gCurrentEntry << std::endl;
  hcalgui::ledLabel->SetText(TString::Format("LED Bit: %02d, Count: %5d",Int_t(hcalt::ledbit),Int_t(hcalt::ledcount)));

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
  for(r  = 0; r < kNrows; r++) {
    for(c  = 0; c < kNcols; c++) {
      peak[r][c] = 0.0;
      adc[r][c] = 0.0;
      tdc[r][c] = 0.0;
    }
  }
  for(Int_t m = 0; m < hcalt::ndata; m++) {
    r = hcalt::row[m];
    c = hcalt::col[m];
    if(r < 0 || c < 0) {
      std::cerr << "Why is row negative? Or col?" << std::endl;
      continue;
    }
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

void clicked_displayNextButton()
{
  //if(gCurrentEntry>gMaxEntries);
  hcalgui::entryInput->SetIntNumber(++gCurrentEntry);
  displayEvent(gCurrentEntry);
}

void clicked_displayEntryButton()
{
  gCurrentEntry = hcalgui::entryInput->GetIntNumber();
  displayEvent(gCurrentEntry);
}


Int_t hcal_display(Int_t run = 1198, Int_t event = -1)
{
  hcalgui::SetupGUI();
  gStyle->SetLabelSize(0.05,"XY");
  gStyle->SetTitleFontSize(0.08);

  if(!T) { 
    T = new TChain("T");
    T->Add(TString::Format("/volatile/halla/sbs/seeds/rootfiles/hcal_%d*.root",run));
    //T->Add(TString::Format("/adaqfs/home/a-onl/sbs/Rootfiles/hcal_%d*.root",run));
    //T->Add(TString::Format("rootFiles/LED/fadc_f1tdc_%d.root",run));
    T->SetBranchStatus("*",0);
    T->SetBranchStatus("sbs.hcal.*",1);
    T->SetBranchAddress("sbs.hcal.nsamps",hcalt::nsamps);
    T->SetBranchAddress("sbs.hcal.a",hcalt::a);
    T->SetBranchAddress("sbs.hcal.tdc",hcalt::tdc);
    T->SetBranchAddress("sbs.hcal.ledbit",&hcalt::ledbit);
    T->SetBranchAddress("sbs.hcal.ledcount",&hcalt::ledcount);
    T->SetBranchAddress("sbs.hcal.samps",hcalt::samps);
    T->SetBranchAddress("sbs.hcal.samps_idx",hcalt::samps_idx);
    T->SetBranchAddress("sbs.hcal.adcrow",hcalt::row);
    T->SetBranchAddress("sbs.hcal.adccol",hcalt::col);
    T->SetBranchStatus("Ndata.sbs.hcal.adcrow",1);
    T->SetBranchAddress("Ndata.sbs.hcal.adcrow",&hcalt::ndata);
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

