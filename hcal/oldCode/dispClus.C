#include <TH2F.h>
#include <TChain.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <TSystem.h>
#include "hcal.h"

const Int_t kNrows = 24;
const Int_t kNcols = 12;

Int_t Nrows;
Int_t Ncols;

const Int_t DISP_MIN_SAMPLE = 0;
const Int_t DISP_MAX_SAMPLE = 50;

const Int_t DISP_FADC_SAMPLES = (DISP_MAX_SAMPLE-DISP_MIN_SAMPLE);

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
    canv[sub] = fEc1;
    fF5->AddFrame(fEc1,fL4);
    tf->AddFrame(fF5,fL4);
    return tf;
  }

  void SetupHMGUI() {
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
      for(Int_t i = 0; i < 1; i++) {
        tf = AddTabSub(i);
      }
      main->AddFrame(fTab, new TGLayoutHints(kLHintsBottom | kLHintsExpandX |
                                          kLHintsExpandY, 2, 2, 5, 1));
      main->MapSubwindows();
      main->Resize();   // resize to default size
      main->MapWindow();

      subCanv[0] = canv[0]->GetCanvas();
      subCanv[0]->Divide(2,0,0.001,0.001);

    }
  }
};


Double_t nhit = 0;
TH2D *heatMapHisto;
TH2D *clustHisto;
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
    
    // Check event increment and increment
    if( entry == -1 ) {
      gCurrentEntry++;
    } else {
      gCurrentEntry = entry;
    }

    if( gCurrentEntry < 0 ) {
      gCurrentEntry = 0;
    }

    // Get the event from the TTree
    T->GetEntry( gCurrentEntry );
    std::cout << "Displaying event " << gCurrentEntry << std::endl;
    hcalgui::ledLabel->SetText(TString::Format("LED Bit: %02d, Count: %5d",Int_t(hcalt::ledbit),Int_t(hcalt::ledcount)));
  
    int r,c;
    int nblk = hcalt::nblk[0];
    int cblkid[nblk];
    int cid = hcalt::cid[0];

    // Clear old histogram
    heatMapHisto->Reset("ICES");
    clustHisto->Reset("ICES");

    // Declare signal amplitude, adc, and tdc arrays for this event
    double adc_p[kNrows][kNcols] = {0.0};
    //double blk_clus[kNrows][kNcols] = {0.0};

    // Check the cluster and set array values to 1 for hits
    for( int el = 0; el < kNrows*kNcols; el++ ){
      r = el/kNcols;
      c = el%kNcols;
      
      clustHisto->SetBinContent( c+1, r+1, 0.1 );
      
      for( int b = 0; b < nblk; b++ ){
	int cblkid = hcalt::cblkid[b];
	//if( el = hcalt::cblkid[b] ) blk_clus[r][c]=1.0;
	if( el == cblkid-1 ) 
	  clustHisto->SetBinContent( c+1, r+1, 1.0 );
	  
	if( el == cid-1 )
	  clustHisto->SetBinContent( c+1, r+1, 0.5 );
      }
	
    }
      
    // Process event with m data
      for( int m = 0; m < hcalt::ndata; m++ ) {
      // Define row and column
      r = hcalt::row[m];
      c = hcalt::col[m];
      if( r < 0 || c < 0 ) {
	cerr << "Error: row or col negative." << endl;
	continue;
      }
    
      if( r >= kNrows || c >= kNcols ) continue;

      // Fill adc array
      adc_p[r][c] = hcalt::a_p[m];

      heatMapHisto->SetBinContent( c+1, r+1, adc_p[r][c] );
      

    }

    subCanv[0]->cd(1); //Heatmap 2d histogram
    subCanv[0]->SetGrid();
    gStyle->SetOptStat(0);
    gStyle->SetPalette(53);
    heatMapHisto->SetTitle("HCal Pedestal Subtracted ADC (sRAU)");
    heatMapHisto->Draw("colz");
    gPad->Update();

    subCanv[0]->cd(2); //Cluster 2d histogram
    subCanv[0]->SetGrid();
    gStyle->SetOptStat(0);
    clustHisto->SetTitle("HCal Cluster Elements");
    clustHisto->Draw("colz");
    gPad->Update();
    
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

Int_t dispClus(Int_t run = 1198, Int_t event = -1)
{

  // Initialize function with user commands

  cout << "Enter run number for analysis." << endl;
  cin >> run;

  Nrows=1;
  Ncols=1;
  heatMapHisto = new TH2D( "Hits by Module", "", kNcols, 0., kNcols, kNrows, 0., kNrows );
  clustHisto = new TH2D( "Cluster Elements", "", kNcols, 0., kNcols, kNrows, 0., kNrows );


  hcalgui::SetupHMGUI();
  gStyle->SetLabelSize(0.05,"XY");
  gStyle->SetTitleFontSize(0.08);  

  if(!T) { 
    T = new TChain("T");
    T->Add(TString::Format("/volatile/halla/sbs/seeds/rootfiles/hcal_%d*.root",run));
    T->SetBranchStatus("*",0);
    T->SetBranchStatus("sbs.hcal.*",1);
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

    // Add clustering branches
    T->SetBranchAddress("sbs.hcal.clus.id",hcalt::cid);
    T->SetBranchAddress("sbs.hcal.clus.row",hcalt::crow);
    T->SetBranchAddress("sbs.hcal.clus.col",hcalt::ccol);
    T->SetBranchAddress("sbs.hcal.clus.e",hcalt::ce);
    T->SetBranchAddress("sbs.hcal.clus.eblk",hcalt::ceblk);
    T->SetBranchAddress("sbs.hcal.clus.nblk",hcalt::cnblk);
    T->SetBranchAddress("sbs.hcal.nclus",hcalt::nclus);
    T->SetBranchAddress("sbs.hcal.nblk",hcalt::nblk);
    T->SetBranchAddress("sbs.hcal.clus_blk.id",hcalt::cblkid);
    T->SetBranchAddress("sbs.hcal.clus.eblk",hcalt::ceblk);

    
    T->SetBranchStatus("Ndata.sbs.hcal.adcrow",1);
    T->SetBranchAddress("Ndata.sbs.hcal.adcrow",&hcalt::ndata);
    std::cerr << "Opened up tree with nentries=" << T->GetEntries() << std::endl;

  }

  if(T->GetEntries()<=0){
    cout << "Error: no data for this run available." << endl;
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

