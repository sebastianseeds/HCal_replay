#include <TH2F.h>
#include <TChain.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <TSystem.h>
#include "hcal.h"

const Int_t kNrows = 12;
const Int_t kNcols = 12;

const Int_t DISP_MIN_SAMPLE = 0;
const Int_t DISP_MAX_SAMPLE = 30;  //For extended fADC window LED pulse setup
const Int_t DISP_FADC_SAMPLES = (DISP_MAX_SAMPLE-DISP_MIN_SAMPLE);

const Int_t minADC = 0;
const Int_t maxADC = 4000;

Int_t gCurrentEntry = -1;

TChain *T = 0;
TCanvas *C1 = new TCanvas();

//Augment fn's for cosmic Calibration
void goodHistoTest(TH1F*, double, int, int);
void getPedestal(Int_t);
void pulseSpec(TH1F*, TF1*, double, int, int);

//Augmenting with parameters necessary to obtain integrated landau fit values from good pulses
Double_t nhit = 0;
TH1F *histos[kNrows][kNcols];
TH1F *PMTIntSpec[kNrows][kNcols]; //=new TH1F("","",20,0,5000); //Empirical limits
TH1F *PMTMaxSpec[kNrows][kNcols]; //=new TH1F("","",20,0,1000); //Empirical limits
TH1F *passHisto = new TH1F("passHisto","passHisto", 19, 0, 19);
vector<double> intValues;
vector<double> maxValues;
vector<int> intValues_row;
vector<int> intValues_col;
int signalTotal = 1;
TFile *HistosFile = new TFile("HistosFile.root","RECREATE");  //File for checking histogram fits and integrity
int checkHistos = 0;  //Iterator to save only some histograms for a check
TF1 *fFit;  //Declaration for landau fit check
double pedestals[kNrows][kNcols];
Bool_t gSaturated[kNrows][kNcols];
int badFit = 0; //Counts up per landau fit failing chi^2 test
int chisqr = 0; //Counts up per chi sqr <1 fit
int DRAWNHISTO = 0; //Counts up per histo drawn
double con;
double mu;
double sig;
double xmaxy;
double maxy;

TH1F* MakeHisto(Int_t row, Int_t col, Int_t bins){
  TH1F *h = new TH1F(TString::Format("h%02d%02d",row,col),
      TString::Format("%d-%d",row+1,col+1),bins,DISP_MIN_SAMPLE,DISP_MAX_SAMPLE);
  h->SetStats(0);
  h->SetLineWidth(2);
  return h;
}

//Rows and Columns start at zero and go to kNrows-1 and kNcols-1
void displayEvent(Int_t entry = -1){
  if(entry == -1) {
    gCurrentEntry++;
  } else {
    gCurrentEntry = entry;
  }

  if(gCurrentEntry<0) {
    gCurrentEntry = 0;
  }

  T->GetEntry(gCurrentEntry);

  Int_t r,c,idx,n,sub;
  // Clear old histograms, just in case modules are not in the tree
  for(r = 0; r < kNrows; r++) {
    for(c = 0; c < kNcols; c++) {
      histos[r][c]->Reset("ICES M");  //Resets integrals, contents, errors, stats, min/max (all)
      gSaturated[r][c] = false;  //Sets all saturation values to zero
    }
  }

  // Populate arrays to store values by pmt. Set all to zero initially.
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

  //
  for(Int_t m = 0; m < hcalt::ndata; m++) {
    r = hcalt::row[m]-1; //Why the -1 here?
    c = hcalt::col[m]-1; //Why the -1 here?
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
    bool displayed = false;
    for(Int_t s = DISP_MIN_SAMPLE; s < DISP_MAX_SAMPLE && s < n; s++) {
      displayed = true;
      histos[r][c]->SetBinContent(s+1-DISP_MIN_SAMPLE,hcalt::samps[idx+s]);
      if(peak[r][c]<hcalt::samps[idx+s])
        peak[r][c]=hcalt::samps[idx+s];
      if(peak[r][c]>4095) {
        gSaturated[r][c] = true;
      }
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
      histos[r][c]->SetTitle(TString::Format("%d-%d (ADC=%g,TDC=%g)",r+1,c+1,adc[r][c],tdc[r][c]));
      if(gSaturated[r][c])
        histos[r][c]->SetLineColor(kRed+1);
      else
        histos[r][c]->SetLineColor(kBlue+1);
      if(tdc[r][c]!=0){
        histos[r][c]->SetLineColor(kGreen+1);
      }
      //if(tdc[r][c]!=0 && !gSaturated[r][c]){ //Commented per recommendation to cut on tdc last.
      if(!gSaturated[r][c]){
	goodHistoTest(histos[r][c],tdc[r][c],r,c);
	checkHistos++;
      }
    }
  }
}

//Open pedestal text file and read out all pedestal values by HV run
void getPedestal(int run = 989){  //is 0-12 in file along the row or column?---------------------------------------*
  int lineNumber = 0;
  ifstream pedFile(Form("Pedestals_%i.txt",run));

  cout << "Getting pedestal values from run " << run << ".." << endl;
  
  int n1;
  double d1, d2;

  int rval, cval;
  string line;
  
  while (getline(pedFile,line))
    {
      if(line.at(0) == '#'){
	continue;
      }

      stringstream ss(line);
      ss >> n1 >> d1 >> d2;
      
      rval = floor(n1/kNrows);
      cval = n1 % kNcols;
      
      pedestals[rval][cval] = d1; 
      
      //cout << "row = " << rval << ", col = " << cval << ", array val = " << pedestals[rval][cval] << endl;
    }

  cout << "Pedestal values extracted." << endl;
  
}

//Saturation and TDC test already applied in display(). Test now for chi^2 goodness of fit to landau.
void goodHistoTest(TH1F *testHisto, double tdcVal, int row, int col){
  
  //Get pedestal value for the run and pmt 
  double pedVal = pedestals[row][col];
  
  //Subtract pedestal value from each bin and cut off bin 20 (overflow)
  for(int b=1 ; b<=testHisto->GetNbinsX() ; b++) {
    testHisto->SetBinContent(b,testHisto->GetBinContent(b)-pedVal);
    //if(testHisto->GetBinContent(b) < 0) testHisto->SetBinContent(b,0.0); //romoved to check pedestals
  }

  
  //Pass threshold from max value as cut
  
  if(testHisto->GetBinContent(testHisto->GetMaximumBin())>50 && testHisto->GetMaximumBin()>2 && testHisto->GetMaximumBin()<17){  //Cut out noise by ensuring all pulses are greater than the std dev avg of the ped and that the max lies within the center of the histogram

    //cout << "Threshold values are: Max value " << testHisto->GetBinContent(testHisto->GetMaximumBin()) << " at bin " << testHisto->GetMaximumBin() << endl;
    
    testHisto->Fit("landau", "Q", 0, 19);  //Fit histogram with landau, limits window minus overflow on 20
    fFit=testHisto->GetFunction("landau"); //Set error of bin content
    
    int dof = 18;  //20 bins -1 per histo, removing last bin to erase overflow -> 19 - 1 = 18
    double rChiSqr = (fFit->GetChisquare())/dof;
    
    if (rChiSqr > 500.0){  //Keep updating this empirically. Set very high to cut out very bad fits
      if(badFit % 20 == 0){
	HistosFile->cd();
	testHisto->SetTitle(Form("Bad Fit R%d C%d",row,col));
	testHisto->Write(Form("badFitHisto%d_rChiSqr_gt5_r%d_c%d",badFit,row,col));
	testHisto->Draw("AP");
	DRAWNHISTO++;
      }
      
      badFit++;
      
    } else {
      if(rChiSqr < 1 && chisqr%10==0){
	cout << "Chi square <1 event detected" << endl;
	HistosFile->cd();
	testHisto->SetTitle(Form("Chi Square <1 Fit R%d C%d",row,col));
	testHisto->Write(Form("FitHisto_rChiSqr_lt1_r%d_c%d",row,col));
	testHisto->Draw("AP");
	chisqr++;
      }
      //if(row==0 && col==0 && testHisto->GetBinContent(testHisto->GetMaximumBin())<50){
      //HistosFile->cd();
      //testHisto->Write(Form("BadFitHisto_NoiseFit_r%d_c%d",row,col));
      //testHisto->Draw("AP");
      //}
      pulseSpec(testHisto, fFit, tdcVal, row, col); //Call integratePulse() for histograms that pass all tests
      
    }
  }  
}

//integrate good signals, get max from good signals, return spectrums and averages
void pulseSpec(TH1F *intHisto, TF1 *fFit, double tdcVal, int row, int col){

  
  //intValues.push_back(fFit->Integral(0,1000));  //Pushing integral way out to approximate indefinite int
  intValues.push_back(intHisto->Integral(0,19));  //Just getting total bin content
  //Using TMath and obtaining parameters from fit function - similar to method in analyze_cosmics.C l:377
  con = fFit->GetParameter(0);
  mu = fFit->GetParameter(1);
  sig = fFit->GetParameter(2);
  xmaxy = fFit->GetMaximumX(mu*0.5,mu*1.5);
  maxy  = TMath::Landau(xmaxy,mu,sig)*con;

  maxValues.push_back(maxy);
  intValues_row.push_back(row);
  intValues_col.push_back(col);
    
  //PMTIntSpec[row][col]->Fill(fFit->Integral(0,1000));
  PMTIntSpec[row][col]->Fill(intHisto->Integral(0,19)); //Just getting total bin content
  
  PMTMaxSpec[row][col]->Fill(maxy);
  
  //if(row == 0 && col == 0) cout << "ROW AND COL ZERO" << endl; //Not sure why r=0,c=0 is always empty

  if (signalTotal % 10000 == 0){
    cout << "Writing a reference good fit histogram to file.." << endl;
    cout << "Maximum y value for this histogram = " << maxy << "." << endl;
    cout << "Pedestal for same histogram = " << pedestals[row][col] << "." << endl;
    cout << "TDC value for same histogram = " << tdcVal << "." << endl << endl;
    HistosFile->cd();
    intHisto->SetTitle(Form("Good Fit R%d C%d",row,col));
    intHisto->Write(Form("goodFitHisto%d_r%d_c%d",signalTotal,row,col));
    intHisto->Draw("AP");
    DRAWNHISTO++;
  }
  signalTotal++;
}

Int_t LEDCal(Int_t run = 989, Int_t event = -1)
{
  //Build spectrum histograms. Empirical limits.

  cout << "Building spectrum histograms.." << endl;
  for(int r=0; r<kNrows; r++){
    for(int c=0; c<kNcols; c++){  
      PMTIntSpec[r][c] = new TH1F(Form("Int ADC Spect R%d C%d",r,c),Form("Int ADC Spect R%d C%d",r,c),50,0,5000); 
      PMTMaxSpec[r][c] = new TH1F(Form("Max ADC Spect R%d C%d",r,c),Form("Max ADC Spect R%d C%d",r,c),50,0,1000); 
    }
  }
  
  getPedestal(run);
  
  if(!T) { 
    T = new TChain("T");
    T->Add(TString::Format("fadc_f1tdc_%d.root",run));
    T->SetBranchStatus("*",0);
    T->SetBranchStatus("sbs.hcal.*",1);
    T->SetBranchAddress("sbs.hcal.nsamps",hcalt::nsamps);
    T->SetBranchAddress("sbs.hcal.a",hcalt::a);
    T->SetBranchAddress("sbs.hcal.tdc",hcalt::tdc);
    T->SetBranchAddress("sbs.hcal.ledbit",&hcalt::ledbit);
    T->SetBranchAddress("sbs.hcal.ledcount",&hcalt::ledcount);
    T->SetBranchAddress("sbs.hcal.samps",hcalt::samps);
    T->SetBranchAddress("sbs.hcal.samps_idx",hcalt::samps_idx);
    T->SetBranchAddress("sbs.hcal.row",hcalt::row);
    T->SetBranchAddress("sbs.hcal.col",hcalt::col);
    T->SetBranchStatus("Ndata.sbs.hcal.row",1);
    T->SetBranchAddress("Ndata.sbs.hcal.row",&hcalt::ndata);
    cout << "Opened up tree with nentries=" << T->GetEntries() << endl;
    for(Int_t r = 0; r < kNrows; r++) {
      for(Int_t c = 0; c < kNcols; c++) {
        histos[r][c] = MakeHisto(r,c,DISP_FADC_SAMPLES);
        gSaturated[r][c] = false;
      }
    }
  }

  gCurrentEntry = event;

  cout << "Total events to process " << T->GetEntries() << endl;
  
  for (int i = gCurrentEntry; i < T->GetEntries(); i++){ //Not sure about GetEntries() here
    displayEvent(gCurrentEntry);
    gCurrentEntry++;
    
    //Keep count of events processed for monitoring
    if (gCurrentEntry%10000 == 0){
      cout << "Current event: " << gCurrentEntry << endl;
    }
  }
  cout << "Finished loop over run " << run << "." << endl;
  cout << "Total good signals = " << signalTotal << "." << endl;

  cout << "Obtaining averages for max value and integrated pulses.." << endl;
  //Now take in three vectors and get average integrated and max value pulse per pmt (unique r c value)
  double pmtIntTotal[kNrows*kNcols]={0};
  double pmtEvTotal[kNrows*kNcols]={0};
  double pmtMaxTotal[kNrows*kNcols]={0};

  for(int ir=0; ir<kNrows; ir++){
    
    for(int ic=0; ic<kNcols; ic++){
      
      for(int iev=0; iev<intValues.size(); iev++){
	
	int evRow = intValues_row[iev];
	int evCol = intValues_col[iev];

	if(evRow == ir && evCol == ic){
	  
	  pmtIntTotal[ir*kNrows+ic] += intValues[iev]; //build value for int tot, one val per pmt
	  pmtMaxTotal[ir*kNrows+ic] += maxValues[iev];
	  pmtEvTotal[ir*kNrows+ic] += 1.0;
	  
	}
      }
    }
  }

  //cout << "pmtIntTotal[0] = " <<  pmtIntTotal[0]  << "." << endl;

  //cout << "pmtIntTotal size = " << pmtIntTotal.size() << "." << endl;

  //cout << "pmtMaxTotal size = " << pmtMaxTotal.size() << "." << endl;

  //cout << "pmtEvTotal[0] = " << pmtEvTotal[0] << "." << endl;

  //cout << "pmtEvTotal.size = " << pmtEvTotal.size() << "." << endl;

  cout << "Writing spectrum histograms to file.." << endl;
  
  for(int r=0; r<kNrows; r++){
    for(int c=0; c<kNcols; c++){
      HistosFile->cd();
      PMTIntSpec[r][c]->Write(Form("Int ADC Spect R%d C%d",r,c));
      PMTIntSpec[r][c]->Draw("AP");

      HistosFile->cd();
      PMTMaxSpec[r][c]->Write(Form("Max ADC Spect R%d C%d",r,c));
      PMTMaxSpec[r][c]->Draw("AP");
    }
  }
  
  ofstream outFile;

  outFile.open(Form("AvePulseRun_%i.txt",run));

  outFile << "#Averaged int pulse and pulse max data for run " << run << "." << endl;
  outFile << "#Total number of pulse events is " << intValues.size() << "." << endl;
  outFile << "#Row  Column  pmtIntValue  pmtMaxValue  EvTotal" << endl;
  
  for(int e=0; e<(kNrows*kNcols); e++){
    if(pmtEvTotal[e] == 0){
      pmtIntTotal[e]=-1.0; //Make empty bins obvious
      pmtMaxTotal[e]=-1.0; //Make empty bins obvious
    } else {
    pmtIntTotal[e] = pmtIntTotal[e]/pmtEvTotal[e];
    pmtMaxTotal[e] = pmtMaxTotal[e]/pmtEvTotal[e];
    }
    cout << "PMT row = " << floor(e/kNrows) << ", col = " << e % kNcols << ": Int Pulse Avg = " << pmtIntTotal[e] << " and Max Pulse Avg = " << pmtMaxTotal[e] << " with " << pmtEvTotal[e] << " events averaged over." << endl;
    outFile << floor(e/kNrows) << "  " << e % kNcols << "  " << pmtIntTotal[e] << "  " << pmtMaxTotal[e] << "  " << pmtEvTotal[e] << endl;
  }

  cout << "Integrated averages, maximum averages, and total events per PMT written to file AvePulseRun_" << run << ".txt." << endl;

  cout << "ADC spectrum histograms written to file HistosFile.root" << endl;

  cout << "Checked histograms = " << checkHistos << " written to file HistosFile.root." << endl;
  
  cout << "Total bad fits = " << badFit << " rejected." << endl;

  cout << "Total fitted histograms drawn to file HistosFile.root = " << DRAWNHISTO << "." << endl;
  
  return 0;
}

