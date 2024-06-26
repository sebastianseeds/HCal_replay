#include <TH2F.h>
#include <TChain.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <TSystem.h>
#include "hcal.h"

const Int_t kNcols = 16;

const Int_t kNumModules = kNcols;

const Int_t DISP_MIN_SAMPLE = 0;
const Int_t DISP_MAX_SAMPLE = 50;

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
    tf = fTab->AddTab(Form("HCAL Trigger fADC Sub%d",sub+1));

    TGCompositeFrame *fF5 = new TGCompositeFrame(tf, (12+1)*kCanvSize,(6+1)*kCanvSize , kHorizontalFrame);
    TGLayoutHints *fL4 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX |
        kLHintsExpandY, 5, 5, 5, 5);
    TRootEmbeddedCanvas *fEc1 = new TRootEmbeddedCanvas(Form("hcalSubCanv%d",sub), fF5, 6*kCanvSize,8*kCanvSize);
    canv[sub] = fEc1;
    fF5->AddFrame(fEc1,fL4);
    tf->AddFrame(fF5,fL4);
    return tf;
  }

  void SetupTGUI() {
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
      tf = AddTabSub(0);
      
      main->AddFrame(fTab, new TGLayoutHints(kLHintsBottom | kLHintsExpandX |
                                          kLHintsExpandY, 2, 2, 5, 1));
      main->MapSubwindows();
      main->Resize();   // resize to default size
      main->MapWindow();

      subCanv[0] = canv[0]->GetCanvas();
      subCanv[0]->Divide(4,4,0.001,0.001);
    }
  }
};


Double_t nhit = 0;
TH1F *histos[kNcols];
Bool_t gSaturated[kNcols];


TH1F* MakeHisto(Int_t col, Int_t bins)
{
  TH1F *h = new TH1F(TString::Format("h%02d",col),
      TString::Format("%d",col+1),bins,DISP_MIN_SAMPLE,DISP_MAX_SAMPLE);
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
  
  Int_t c,Tidx,Tn;
  // Clear old histograms, just in case modules are not in the tree
  
  for(c = 0; c < kNcols; c++) {
    histos[c]->Reset("ICES M");
    gSaturated[c] = false;
  }
  
  Double_t peak[kNcols];
  Double_t Tadc[kNcols];
  
  for(c  = 0; c < kNcols; c++) {
    peak[c] = 0.0;
    Tadc[c] = 0.0;
  }
  
  for(Int_t m = 0; m < hcalt::Tndata; m++) {
    c = hcalt::Tcol[m];
    peak[c] = hcalt::Ta_amp[m];
    if(c < 0) {
      std::cerr << "Error: column negative." << std::endl;
      continue;
    }
    
    if(c >= kNcols)
      continue;
    
    // Fill adc, tdc, and cluster arrays
    Tidx = hcalt::Tsamps_idx[m];
    Tn = hcalt::Tnsamps[m];
    Tadc[c] = hcalt::Ta_p[m];
    bool displayed = false;
    for(Int_t s = DISP_MIN_SAMPLE; s < DISP_MAX_SAMPLE && s < Tn; s++) {
      displayed = true;
      histos[c]->SetBinContent(s+1-DISP_MIN_SAMPLE,hcalt::Tsamps[Tidx+s]); //Prepared for trigger samples to be added to the root tree
      if(peak[c]<hcalt::Tsamps[Tidx+s])
	peak[c]=hcalt::Tsamps[Tidx+s];
      if(peak[c]>4095) {
	gSaturated[c] = true;
      }
    }
    if(!displayed) {
      std::cerr << "Skipping empty module: " << m << std::endl;
      for(Int_t s = 0;  s < DISP_FADC_SAMPLES; s++) {
	histos[c]->SetBinContent(s+1,-404);
      }
    }
  }
  
  for(c = 0; c < kNcols; c++) {
    subCanv[0]->cd( c + 1 );
    histos[c]->SetTitle(TString::Format("%d (ADC=%g)",c+1,Tadc[c]));
    histos[c]->SetMaximum(60);
    if(gSaturated[c])
      histos[c]->SetLineColor(kRed+1);
    else
      histos[c]->SetLineColor(kBlue+1);
    histos[c]->Draw();
    gPad->Update();
    std::cout << " [" << c << "]=" << peak[c];
  }
}

void clicked_displayNextButton()
{
  hcalgui::entryInput->SetIntNumber(++gCurrentEntry);
  displayEvent(gCurrentEntry);
}

void clicked_displayEntryButton()
{
  gCurrentEntry = hcalgui::entryInput->GetIntNumber();
  displayEvent(gCurrentEntry);
}

Int_t dissBB(Int_t run = 1198, Int_t event = -1)
{

  // Initialize function with user commands

  cout << "Enter run number for analysis." << endl;
  cin >> run;

  hcalgui::SetupTGUI();
  gStyle->SetLabelSize(0.05,"XY");
  gStyle->SetTitleFontSize(0.08);  

  if(!T) { 
    T = new TChain("T");
    //T->Add(TString::Format("/volatile/halla/sbs/seeds/rootfiles/hcal_%d*.root",run));
    T->Add(TString::Format("../Rootfiles/bbshower_%d*.root",run));
    T->SetBranchStatus("*",0);
    
    T->SetBranchStatus("bb.ps.*",1);
    T->SetBranchStatus("bb.sh.*",2);

    /*
    T->SetBranchStatus("sbs.hcal.*",1);
    T->SetBranchStatus("sbs.trig.*",2);
    T->SetBranchAddress("sbs.hcal.nsamps",hcalt::nsamps);
    T->SetBranchAddress("sbs.hcal.a",hcalt::a);
    T->SetBranchAddress("sbs.hcal.a_p",hcalt::a_p);
    T->SetBranchAddress("sbs.hcal.a_amp",hcalt::a_amp);
    T->SetBranchAddress("sbs.hcal.a_amp_p",hcalt::a_amp_p);
    T->SetBranchAddress("sbs.hcal.ped",hcalt::ped);
    T->SetBranchAddress("sbs.hcal.tdc",hcalt::tdc);
    T->SetBranchAddress("sbs.hcal.ledbit",&hcalt::ledbit);
    T->SetBranchAddress("sbs.hcal.ledcount",&hcalt::ledcount);
    T->SetBranchAddress("sbs.hcal.samps",hcalt::samps);
    T->SetBranchAddress("sbs.hcal.samps_idx",hcalt::samps_idx);
    T->SetBranchAddress("sbs.hcal.adcrow",hcalt::row);
    T->SetBranchAddress("sbs.hcal.adccol",hcalt::col);
    
    // Add trigger generic detector branches
    T->SetBranchAddress("sbs.trig.a_p",hcalt::Ta_p);
    T->SetBranchAddress("sbs.trig.adcelemID",hcalt::TelemID);
    T->SetBranchAddress("sbs.trig.nsamps",hcalt::Tnsamps);
    T->SetBranchAddress("sbs.trig.samps",hcalt::Tsamps);
    T->SetBranchAddress("sbs.trig.samps_idx",hcalt::Tsamps_idx);
    T->SetBranchAddress("sbs.trig.adccol",hcalt::Tcol);
    T->SetBranchAddress("sbs.trig.a_amp_p",hcalt::Ta_amp_p);
    T->SetBranchAddress("sbs.trig.a_amp",hcalt::Ta_amp);
    
    // Add clustering branches
    T->SetBranchAddress("sbs.hcal.clus.id",hcalt::cid);
    T->SetBranchAddress("sbs.hcal.clus.row",hcalt::crow);
    T->SetBranchAddress("sbs.hcal.clus.col",hcalt::ccol);
    T->SetBranchAddress("sbs.hcal.clus.e",hcalt::ce);
    T->SetBranchAddress("sbs.hcal.clus.eblk",hcalt::ceblk);

    T->SetBranchStatus("Ndata.sbs.hcal.adcrow",1);
    T->SetBranchAddress("Ndata.sbs.hcal.adcrow",&hcalt::ndata);

    T->SetBranchStatus("Ndata.sbs.trig.adccol",2);
    T->SetBranchAddress("Ndata.sbs.trig.adccol",&hcalt::Tndata);
    */


    std::cerr << "Opened up tree with nentries=" << T->GetEntries() << std::endl;
    for(Int_t c = 0; c < kNcols; c++) {
      histos[c] = MakeHisto(c,DISP_FADC_SAMPLES);
      gSaturated[c] = false;
    }
  }

  if( T->GetEntries()<=0 ){
    cout << "Error: no data." << endl;
    return 0;
  }

  gCurrentEntry = event;
  while( user_input != "q" ) {
    if(is_number(user_input)) {
      gCurrentEntry = std::stoi(user_input);
    } else {
      gCurrentEntry++;
    }
    displayEvent(gCurrentEntry);
    std::cout << "Display options: <enter> == next event, or q to stop." << std::endl;
    getline(std::cin,user_input);
  }

  return 0;
}

