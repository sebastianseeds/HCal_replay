#include <TH2F.h>
#include <TChain.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <TSystem.h>
#include "hcal.h"
#include <vector>
#include <TStopwatch.h>
#include <TMath.h>
#include <TF1.h>
#include<TStyle.h>
#include<TList.h>
#include<TChain.h>
using namespace std;

const Int_t nleds = 5;      //Number of LEDs used in voltage scan.
const Int_t kNrows = 12;
const Int_t kNcols = 12;
const Int_t channels = kNrows * kNcols;
Float_t pedestal[channels], pedestal_int[channels], ped_std_dev[channels], ped_int_std_dev[channels];
Int_t gCurrentEntry = 0;
Float_t led_means[channels][nleds],led_int_means[channels][nleds],led_std_dev[channels][nleds],led_int_std_dev[channels][nleds],npe[channels][nleds];
const Int_t nscans = 8;
Int_t runs[nscans] = {1198,1200,1201,1202,1203,1204,1205,1206};
Int_t hv_cmu[nscans] = {1200,1250,1300,1350,1400,1450,1425,1375};
Int_t hv_jlab[nscans] = {1200,1250,1300,1350,1400,1450,1425,1375};
Int_t hvmin = hv_cmu[0]-50, hvmax = hv_cmu[5]+50;
Int_t r,c;

TChain *T = 0;

//Create power function for fitting.
Double_t fit_pow(Double_t *X,Double_t *par) 
{
  //Double_t fitval = par[0]*pow(X[0],par[1]);
  Double_t fitval = par[0]*pow((X[0]-par[1]),par[2]) + par[3];
  return fitval;
}

//Create Gaussian for fitting.
Double_t fit_gaus(Double_t *X,Double_t *par) 
{
  Double_t fitval = par[0]*TMath::Exp(-0.5*pow(((X[0]-par[1])/par[2]),2));
  return fitval;
}

void Plot_LED_Voltage_Scans(Int_t run = 1203)
{
  //Define a new stopwatch.
  TStopwatch *st=new TStopwatch();
  st->Start(kTRUE);

  gROOT->SetBatch(kTRUE);//Always run in batch mode otherwise 100s of plot will appear.

    if(!T) 
    { 
      T = new TChain("T");
      
      //============  Reading the Rootfile =======================//
      
      for(Int_t i=0;i<nscans;i++)
	{
	  const TString rootfilePath = "/home/daq/test_fadc/Voltage_Scans/LEDs/";
	  std::ostringstream str;
	  str << rootfilePath<<"LED_Voltage_Scan_Run"<<runs[i];
	  TString basename = str.str().c_str();
	  TString rootfile = basename + ".root";
	  cout<<basename<<endl;
	  
	  Long_t split=0;
	  char* file = 0;
	  
	  //====adding splits rootfiles =======================//
	  
	  Long_t u=0;
	  while ( !gSystem->AccessPathName(rootfile.Data()) ) 
	    {
	      T->Add(rootfile.Data());
	      cout << "ROOT file " << rootfile << " added to TChain." << endl;
	      u++;
	      rootfile = basename + "__" + u + ".root";
	    }
	}
      if(!T->GetEntries())
	{
	  cerr<< "No root file was found" << endl;
	  return;
	}
      //==finish adding splits rootfiles=====================//
      
      T->SetBranchStatus("*",1);
      T->SetBranchAddress("pedestal",pedestal);
      T->SetBranchAddress("ped_std_dev",ped_std_dev);
      T->SetBranchAddress("led_means",led_means);
      T->SetBranchAddress("led_std_dev",led_std_dev);
      T->SetBranchAddress("pedestal_int",pedestal_int);
      T->SetBranchAddress("ped_int_std_dev",ped_int_std_dev);
      T->SetBranchAddress("led_int_means",led_int_means);
      T->SetBranchAddress("led_int_std_dev",led_int_std_dev);
      T->SetBranchAddress("npe",npe);
      /* 
      T->SetBranchAddress("sbs.hcal.nsamps",hcalt::nsamps);       //Number of samples for given row-col
      T->SetBranchAddress("sbs.hcal.a",hcalt::a);                 //Raw ADC amplitudes
      T->SetBranchAddress("sbs.hcal.tdc",hcalt::tdc);             //Raw TDC value
      //T->SetBranchAddress("sbs.hcal.ledbit",&hcalt::ledbit);
      //T->SetBranchAddress("sbs.hcal.ledcount",&hcalt::ledcount);
      T->SetBranchAddress("sbs.hcal.samps",hcalt::samps);         //RAW ADC samples
      T->SetBranchAddress("sbs.hcal.samps_idx",hcalt::samps_idx); //Index in samples vector for given row-col module 
      T->SetBranchAddress("sbs.hcal.row",hcalt::row);             //Row for block in data vectors
      T->SetBranchAddress("sbs.hcal.col",hcalt::col);             //Col for block in data vectors
      T->SetBranchAddress("sbs.hcal.ledbit",&hcalt::ledbit);       //Stores the LED bit.
      T->SetBranchStatus("Ndata.sbs.hcal.row",1);
      T->SetBranchAddress("Ndata.sbs.hcal.row",&hcalt::ndata);    //Number of fADC samples per PMT per event.
      */
      std::cerr << "Opened up tree with nentries=" << T->GetEntries() << std::endl;
    }

    TFile* file = new TFile(Form("/home/daq/test_fadc/Voltage_Scans/LEDs/LED_Voltage_Scan_Results.root"),"RECREATE");
    gStyle->SetOptFit(1111);
    //Create histograms to hold the voltage scan results.
    TH2F *hhv_scan[channels][nleds];
    for(Int_t i = 0; i<channels; i++)
      {
	for(Int_t j = 0; j<nleds; j++)
	  {
	    //gStyle->SetOptFit(1111);
	    hhv_scan[i][j] = new TH2F(Form("hpmt%d_led%d",i,j+1),Form("HV Scan for Module %d and LED %d",i,j+1),100,hvmin,hvmax,1000,0,8200);
	    hhv_scan[i][j]->GetXaxis()->SetTitle("HV (-V)");
	    hhv_scan[i][j]->GetXaxis()->CenterTitle();
	    hhv_scan[i][j]->GetYaxis()->SetTitle("Average fADC Amplitude from LED");
	    hhv_scan[i][j]->GetYaxis()->CenterTitle();
	    hhv_scan[i][j]->SetMarkerStyle(8);
	    hhv_scan[i][j]->SetMarkerSize(1);
	  }
      }

    //Create histograms to hold the voltage scan fADC integral results.
    TH2F *hhv_scan_int[channels][nleds];
    for(Int_t i = 0; i<channels; i++)
      {
	for(Int_t j = 0; j<nleds; j++)
	  {
	    //gStyle->SetOptFit(1111);
	    hhv_scan_int[i][j] = new TH2F(Form("hpmt%d_led%d_int",i,j+1),Form("HV Scan fADC Integrals for Module %d and LED %d",i,j+1),100,hvmin,hvmax,10000,0,100000);
	    hhv_scan_int[i][j]->GetXaxis()->SetTitle("HV (-V)");
	    hhv_scan_int[i][j]->GetXaxis()->CenterTitle();
	    hhv_scan_int[i][j]->GetYaxis()->SetTitle("Average fADC Integral from LED");
	    hhv_scan_int[i][j]->GetYaxis()->CenterTitle();
	    hhv_scan_int[i][j]->SetMarkerStyle(8);
	    hhv_scan_int[i][j]->SetMarkerSize(1);
	  }
      }

    T->GetEntry(0);

    for(Int_t i=0;i<channels;i++)
      {
	//cout<<Form("Pedestal for Module %d = %.1f with standard deviation of %.2f.",i,pedestal[i],ped_std_dev[i])<<endl;
      }

    for(Int_t i=0;i<channels;i++)
      {
	for(Int_t j=0;j<nleds;j++)
	  {
	    //cout<<Form("Mean for Module %d LED %d = %.1f with standard deviation of %.2f. %.2f photoelectrons were detected.",i,j+1,led_means[i][j],led_std_dev[i][j],npe[i][j])<<endl;
	  }
      }
    //Rootfile to hold voltage scan results plots.
    //TFile* file = new TFile(Form("/home/daq/test_fadc/Voltage_Scans/LEDs/LED_Voltage_Scan_Results.root"),"RECREATE");
    //TCanvas* c1=new TCanvas("c1");
    //TF1 *func_pow_fit[channels][nleds];
    TCanvas* c1[channels][nleds];
    TCanvas* c1_int[channels][nleds];
    for(Int_t i=0;i<channels;i++)
      {
	for(Int_t j=0;j<nleds;j++)
	  {
	    c1[i][j] = new TCanvas(Form("PMT%d_LED%d",i,j));
	    c1_int[i][j] = new TCanvas(Form("PMT%d_LED%d_int",i,j));
	  }
      }

    //Loop over all rootfiles to access the LED fADC fit results for each data run.
    for(Int_t i=0;i<nscans;i++)
      {
	T->GetEntry(gCurrentEntry);
	
	//Loop over all PMTs.
	for(Int_t j=0;j<channels;j++)
	  {
	    //Loop over all LEDs in scan.
	    for(Int_t k=0;k<nleds;k++)
	      {
		hhv_scan[j][k]->Fill(hv_cmu[i],led_means[j][k]-pedestal[j]);
		hhv_scan_int[j][k]->Fill(hv_cmu[i],led_int_means[j][k]-pedestal_int[j]);
	      }
	  }
	gCurrentEntry++;
      }
 
    //Store the results histograms after fitting them with a power function.
    TF1 *func_pow_fit[channels][nleds];
    TF1 *func_pow_fit_int[channels][nleds];
    //gStyle->SetOptFit(1111);

    //Loop over all PMTs.
    for(Int_t i=0;i<channels;i++)
      {
	//Assign HCal rows and cols (counting from 1) to each PMT module.
	r=(i/12)+1;
	c=i-12*(r-1)+1;
	//cout<<Form("PMT %d is r%d-c%d.",i,r,c)<<endl;
	//Loop over all LEDs in scan.
	for(Int_t j=0;j<nleds;j++)
	  {
	    func_pow_fit[i][j] = new TF1("func_pow_fit",fit_pow, hvmin+50, hvmax+50, 4);
	    func_pow_fit[i][j]->SetLineColor(2);
	    func_pow_fit[i][j]->SetNpx(4000);

	    func_pow_fit_int[i][j] = new TF1("func_pow_fit_int",fit_pow, hvmin+50, hvmax+50, 4);
	    func_pow_fit_int[i][j]->SetLineColor(2);
	    func_pow_fit_int[i][j]->SetNpx(4000);

	    /*
	    func_pow_fit[i][j]->SetParameter(0,1);
	    func_pow_fit[i][j]->SetParameter(1,12);
	    */
	    
	    //Set Fit parameters for JLab and CMU PMTs separately. (JLab are cols 5-8.)
	    if(c==5||c==6||c==7||c==8)
	      {
		func_pow_fit[i][j]->SetParameter(0,0);
		//func_pow_fit[i][j]->SetParLimits(0,0,1);
		func_pow_fit[i][j]->SetParameter(1,10);
		func_pow_fit[i][j]->SetParLimits(1,0,1500);
		func_pow_fit[i][j]->SetParameter(2,8);
		func_pow_fit[i][j]->SetParLimits(2,1,100);
		func_pow_fit[i][j]->SetParameter(3,0);
		//func_pow_fit[i][j]->SetParLimits(3,-200,200);
	      }
	    else
	      {
		func_pow_fit[i][j]->SetParameter(0,0);
		//func_pow_fit[i][j]->SetParLimits(0,0,1);
		func_pow_fit[i][j]->SetParameter(1,1200);
		func_pow_fit[i][j]->SetParLimits(1,0,1500);
		func_pow_fit[i][j]->SetParameter(2,12);
		func_pow_fit[i][j]->SetParLimits(2,1,100);
		func_pow_fit[i][j]->SetParameter(3,0);
		//func_pow_fit[i][j]->SetParLimits(3,-200,200);
	      }
	    
	    /*
	    func_pow_fit[i][j]->SetParameter(0,1);
	    //func_pow_fit[i][j]->SetParLimits(0,1e-60,1e-5);
	    func_pow_fit[i][j]->SetParameter(1,10);
	    //func_pow_fit[i][j]->SetParLimits(1,1,20);
	    func_pow_fit[i][j]->SetParameter(2,200);
	    func_pow_fit[i][j]->SetParLimits(2,100,400);
	    */
	    cout<<Form("Fit for PMT %d and LED %d.",i,j)<<endl;
	    c1[i][j]->cd();
	    hhv_scan[i][j]->Fit(func_pow_fit[i][j],"rm+");

	    //gPad->Update();
	    TPaveText *pt = new TPaveText(0.25,0.7,0.4,0.85,"NDC NB");
	    pt->SetTextSize(0.04);
	    pt->SetFillColor(0);
	    pt->SetTextAlign(12);
	    pt->AddText(Form("%.2e(HV - %.2f)^%.2f + %.2f",func_pow_fit[i][j]->GetParameter(0),func_pow_fit[i][j]->GetParameter(1),func_pow_fit[i][j]->GetParameter(2),func_pow_fit[i][j]->GetParameter(3)));
	    pt->Draw("");
	    //gPad->Update();
	    //pt->DrawFile("/home/daq/test_fadc/Voltage_Scans/LEDs/LED_Voltage_Scan_Results.root");
	    //c1->Update();
	    //gStyle->SetOptFit(1111);
	    hhv_scan[i][j]->Draw("sames");
	    //pt->Draw("");
	    //gStyle->SetOptFit(1111);
	    //gPad->Update();
	    //pt->Write();//Saves only the TPaveText object.
	    //func_pow_fit[i][j]->Write();//Saves only the fit line.
	    hhv_scan[i][j]->Write();
	    c1[i][j]->Write();
	    //c1->Write();//This will save the whole canvas with the TPaveText but still doesn't save stats panel. Also pops out new window from root file which is very irritating.
	    //gPad->Update();
	  }
      }
    
    Double_t pe = 0,pe_int = 0;
    for(Int_t i=0;i<channels;i++)
      {
	for(Int_t j=0;j<nleds;j++)
	  {
	    pe = pow( (led_means[i][j]-pedestal[i]) / pow( pow(led_std_dev[i][j],2.0)-pow(ped_std_dev[i],2.0) ,0.5)  ,2.0);
	    pe_int = pow( (led_int_means[i][j]-pedestal_int[i]) / pow( pow(led_int_std_dev[i][j],2.0)-pow(ped_int_std_dev[i],2.0) ,0.5)  ,2.0);
	    cout<<Form("PMT %d LED %d has Landau fit amplitude NPE = %f and fADC integral NPE = %f.",i,j,pe,pe_int)<<endl;
	  }
      }
   
    file->Close();

    //Access a stored histgram.
    TFile *f = new TFile("/home/daq/test_fadc/Voltage_Scans/LEDs/LED_Voltage_Scan_Run1204.root");
    TH1F *h1 = (TH1F*)f->Get("hpmt1_led4");
    Int_t binx,biny,binz;
    //h1->GetBinXYZ(200,binx,biny,binz);
    //cout<<Form("h1->GetMaximum() = %d.",h1->GetMaximum())<<endl;
    /*
      cout<<Form("h1->GetXaxis()->FindBin(pedestal[1]) = %d.",h1->GetXaxis()->FindBin(pedestal[1]))<<endl;
      cout<<Form("Maximum bin number. h1->GetMaximumBin() = %d.",h1->GetMaximumBin())<<endl;
      cout<<Form("X position of max bin. h1->GetXaxis()->GetBinCenter(h1->GetMaximumBin()) = %d.",h1->GetXaxis()->GetBinCenter(h1->GetMaximumBin()))<<endl;
      cout<<Form("h1->GetBinContent(h1->GetMaximumBin()) = %f.",h1->GetBinContent(h1->GetMaximumBin()))<<endl;
    cout<<h1->GetBinContent(355)<<endl;
    */
    TCanvas* cSingle=new TCanvas("cSingle");
    cSingle->SetGrid();
    hhv_scan[143][4]->Draw();
    //h1->Draw();

    file->Close();

    //Stop and print out stopwatch information.
    st->Stop();
    cout<<"CPU time = "<<st->CpuTime()<<" s = "<<st->CpuTime()/60.<<" min   Real time = "<<st->RealTime()<<" s = "<<st->RealTime()/60.<<" min"<<endl;
}
