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
#include <fstream>

using namespace std;


// compile me by setting COMPILE to 1, then
// g++ `root-config --libs --cflags` led3.C -o led3.exe -Wall
// you can execute me via
// ./led3.exe runnumber
// where runnumber is the run number to analyze (see variable "run" for current default)

//#define COMPILE 1
// otherwise (interactive mode) just comment the line above

// you might want to change here...
Int_t limit_evts = 0;    //0-> analyze all events. 1-> only analyze events up until max_evts.
Int_t max_evts = 1000;   //Maximum number of events to analyze if limit_evts = 1.
// verbosity, 0 or 1
Int_t verbosity = 0;

//const TString rootfilePath = "/Users/vanessabrio/rootfiles/";
const TString rootfilePath = "/Users/sebastianseeds/HCAL/rootFiles/LED/";

//const TString rootfilePath = "./";

const Int_t NLEDS = 5; // total numer of LEDs
const Int_t NEVTS_LED = 1000; // # of events dedicated to each led
const Int_t DISP_MIN_SAMPLE = 0;
const Int_t DISP_MAX_SAMPLE = 30;
const Int_t DISP_FADC_SAMPLES = (DISP_MAX_SAMPLE-DISP_MIN_SAMPLE);
const Int_t kNrows = 24;
const Int_t kNcols = 12;
const Int_t channels = kNrows * kNcols;

Int_t loop_max = 0;   //Dummy variable set equal to max_evts or nevt for the loop.
Double_t entry;
Int_t gCurrentEntry = 0;
Int_t r,c,n,idx;
Int_t nevt;
TChain *T = 0;
Int_t skip = 3;            //Gives number of lines to skip at top of data file.
Int_t nlines = 0;          //Counts number of lines in the data file.
Int_t ncols;               //Set how many columns of data we have in the data file.
char str[1000];            //Variable to read lines of the data file.
Float_t pedestal[channels];

#ifndef COMPILE
void Vanessa_led3(Int_t run = 1288)
{
#else
  int main(int argc, char **argv)
  {
    Int_t run = 1288;
    if (argc == 2){
      run = atoi(argv[1]);
      cout << "run "<< run << endl;
    }
#endif
    //Define a new stopwatch.
    TStopwatch *st=new TStopwatch();
    st->Start(kTRUE);

    if(!T)
    {
      T = new TChain("T");

      //============  Reading the Rootfile =======================//

      std::ostringstream str;
      str << rootfilePath<<"fadc_f1tdc_"<<run;
      TString basename = str.str().c_str();
      TString rootfile = basename + ".root";
      cout<<basename<<endl;

      //====adding splits rootfiles =======================//

      Long_t u=0;
      while ( !gSystem->AccessPathName(rootfile.Data()) )
      {
        T->Add(rootfile.Data());
        cout << "ROOT file " << rootfile << " added to TChain." << endl;
        u++;
        //rootfile = basename + "__" + u + ".root";
        rootfile = basename + "_" + u + ".root";
      }
      //==finish adding splits rootfiles=====================//

      nevt = T->GetEntries();

      if( nevt < 0)
        //if(!T->GetEntries())
      {
        cerr<< "No root file was found" << endl;
#ifndef COMPILE
        return;
#else
        return -1;
#endif
      }

      //T->Add(TString::Format("rootfiles/fadc_f1tdc_%d.root",run));
      T->SetBranchStatus("*",0);
      T->SetBranchStatus("sbs.hcal.*",1);
      T->SetBranchAddress("sbs.hcal.nsamps",hcalt::nsamps);       //Number of samples for given row-col
      T->SetBranchAddress("sbs.hcal.a",hcalt::a);                 //Raw ADC amplitudes
      T->SetBranchAddress("sbs.hcal.tdc",hcalt::tdc);             //Raw TDC value
      T->SetBranchAddress("sbs.hcal.ledbit",&hcalt::ledbit);
      T->SetBranchAddress("sbs.hcal.ledcount",&hcalt::ledcount);
      T->SetBranchAddress("sbs.hcal.samps",hcalt::samps);         //RAW ADC samples
      T->SetBranchAddress("sbs.hcal.samps_idx",hcalt::samps_idx); //Index in samples vector for given row-col module
      T->SetBranchAddress("sbs.hcal.row",hcalt::row);             //Row for block in data vectors
      T->SetBranchAddress("sbs.hcal.col",hcalt::col);             //Col for block in data vectors
      T->SetBranchStatus("Ndata.sbs.hcal.row",1);
      T->SetBranchAddress("Ndata.sbs.hcal.row",&hcalt::ndata);    //???
      std::cerr << "Opened up tree with nentries=" << nevt << std::endl;
    }


    cout<<"Hello"<<endl;
    T->GetEntry(gCurrentEntry);
    cout<<"There!"<<endl;
    cout << "Displaying event " << gCurrentEntry << endl;
    cout<<"ndata = "<<hcalt::ndata<<endl;                                  //Number of channels.
    cout<<"row[] = "<<hcalt::row[16]<<"   col[] = "<<hcalt::col[16]<<endl;
    cout<<"a[] = "<<hcalt::a[1]<<endl;
    cout<<"tdc[] = "<<hcalt::tdc[1]<<endl;

    n = hcalt::nsamps[0];
    cout<<"n = "<<n<<endl;

    idx = hcalt::samps_idx[0];
    cout<<"idx = "<<idx<<endl;

    //nevt = T->GetEntries();
    cout<<"Number of events = "<<nevt<<endl;

    Double_t adc[kNrows][kNcols];
    int ledbitArray[nevt];

    //Loop over all events.
    if(limit_evts==1)
    {
      loop_max = max_evts;
    }
    else
    {
      loop_max = nevt;
    }

    
    //THIS SEEMS QUESTIONABLE. WHY THIS WAY? WHY NOT - Double_t npe[channels][1+NLEDS] = {0} ?
    Double_t npe[channels][1+NLEDS];
    Double_t avg[channels][1+NLEDS];
    Double_t rms[channels][1+NLEDS];
    Int_t neventsled[channels][1+NLEDS];
    memset(npe, 0, channels * (1+NLEDS) * sizeof(Double_t));
    memset(avg, 0, channels * (1+NLEDS) * sizeof(Double_t));
    memset(rms, 0, channels * (1+NLEDS) * sizeof(Double_t));
    memset(pedestal, 0, channels * sizeof(Float_t));
    memset(neventsled, 0, channels * (1+NLEDS) * sizeof(Int_t));
    
    /////////////////////////////////////////////////////////

    // calculates pedestal from first NEVTS_LED events
    for (int entrynum = 0; entrynum<NEVTS_LED ; entrynum++)
    {
      T->GetEntry(entrynum);

      for(Int_t j=0; j<hcalt::ndata; j++)
      {
        idx = hcalt::samps_idx[j];

        for(Int_t s = DISP_MIN_SAMPLE; s < DISP_MAX_SAMPLE && s < n; s++)
        {
          Double_t value = hcalt::samps[idx+s];
          pedestal[ j ] += value;
        }
      }
    }
    for(Int_t j=0; j<hcalt::ndata; j++)
    {
        pedestal[ j ] /= NEVTS_LED;
    }

    TH1D *hsbshcala_noped[channels][NLEDS+1];
    for (int j=0; j < channels; j++){
      for( int led=0; led<NLEDS+1; led++){
     hsbshcala_noped[j][led] = new TH1D(
         Form("hsbshcala_noped_mod_%d_led_%d",j,led),
         Form("hsbshcala_noped_mod_%d_led_%d",j,led),
         10000,-10000,10000);
      }
    }

    // calculate, for every led, the average of adc - pedestal
    for (int entrynum = NEVTS_LED; entrynum<loop_max ; entrynum++)
    {
      T->GetEntry(entrynum);
      memset(adc,      0, channels * sizeof(Double_t));

      if ( hcalt::ledbit > 0 )
      {
        ledbitArray[entrynum] = 1 + (int)( log( hcalt::ledbit ) / log( 2 ) );
      } else {
        ledbitArray[entrynum] = 0;
      }

      if (entrynum % 1000 == 0) cout << "Event "<< entrynum <<endl;

      for(Int_t j=0; j<hcalt::ndata; j++)
      {
        r = hcalt::row[j]-1;
        c = hcalt::col[j]-1;
        adc[r][c] = hcalt::a[j];

        avg[ j ][ ledbitArray[entrynum] ] += ( adc[r][c] - pedestal[j] );
        neventsled[ j ][ ledbitArray[entrynum] ]++ ;

        hsbshcala_noped[j][ledbitArray[entrynum]]->Fill( adc[r][c] - pedestal[j] );
      }
    }
    
    for (Int_t module = 0; module < channels; module++)
    {
      cout <<" Module: "<< module <<" <pedestal>: "<< pedestal[module] / DISP_FADC_SAMPLES <<endl;

      for (Int_t led=0; led<NLEDS+1; led++)
      {
        avg[ module ][ led ] /= neventsled[ module ][ led ];
      }
    }


    // calculate <rms> and npe
    for (int entrynum = NEVTS_LED; entrynum<nevt ; entrynum++)
    {
      T->GetEntry(entrynum);
      int led = ledbitArray[ entrynum ];

      for(Int_t j=0; j<hcalt::ndata; j++)
      {
        rms[ j ][ led ] += pow( hcalt::a[j] - pedestal[j] - avg[ j ][ led ], 2);
      }
    }
    for (Int_t module = 0; module < channels; module++)
    {
      for (Int_t led=0; led<NLEDS+1; led++)
      {
        rms[ module ][ led ] = sqrt( rms[ module ][ led ] / neventsled[ module ][ led ]);
        npe[ module ][ led ] = pow( avg[ module ][ led ] / rms[ module ][ led ], 2);
      }
    }

    //===============// SAS //================//
    ofstream outFile;
    outFile.open(Form("Comparisons/VPed_%d.txt",run));
    //===============// SAS //================//

    // output

    for (Int_t module = 0; module < channels; module++)
    {
      r = hcalt::row[module]-1;
      c = hcalt::col[module]-1;
      printf("[ %2d, %2d ] led:", r, c);

      //===============// SAS //================//
      outFile << pedestal[module]/DISP_FADC_SAMPLES << endl;
      //===============// SAS //================//

      for (int led=0; led<NLEDS + 1;led++)
      {
        printf(" %8d", led);
      }
      cout << endl; 
      cout << "<sbs.hcal.a>   ";
      for (int led=0; led<1+NLEDS; led++)
      {
        printf(" %8.2f", avg[ module ][ led ]);
      }
      cout << endl;
      cout << "<rms>          ";
      for (int led=0; led<1+NLEDS; led++)
      {
        printf(" %8.2f", rms[ module ][ led ]);
      }
      cout << endl;
      cout << "npe            ";
      for (int led=0; led<1+NLEDS; led++)
      {
        printf(" %8.2f", npe[ module ][ led ]);
      }
      cout << endl;
      cout << "nevents        ";
      for (int led=0; led<1+NLEDS; led++)
      {
        printf(" %8.2d", neventsled[ module ][ led ]);
      }
      cout << endl;
      cout << "RMS histogram  ";
      for (int led=0; led<1+NLEDS; led++)
      {
        printf(" %8.2f", hsbshcala_noped[module][led]->GetRMS() );
      }
      cout << endl;

    //hsbshcala_noped->Draw();
    //cout << "RMS histogram: "<< hsbshcala_noped[module]->GetRMS() << endl;
    }
    

    Int_t HV_CMU, HV_JLAB, HV;
    cout << "HV CMU? "<<endl;
    cin >> HV_CMU;
    cout << "HV JLAB? "<<endl;
    cin >> HV_JLAB;

    char *basename = Form("run_%d.csv",run);
    ofstream outfile(basename);

    // calculate and show signal averages
    cout << "channel, HV, pedestal, npe no led, npe led 1, ... npe last led"<< endl;
    outfile << "module, HV, pedestal, npe_led0, npe_led1, npe_led2, npe_led3, npe_led4, npe_led5"<< endl;
    for (Int_t i=0; i<channels; i++)
    {
      HV = ( i % 12 < 4 || i % 12 > 7 ) ? HV_CMU : HV_JLAB;

      printf("%3d %5d", i, HV);
      outfile << i <<"," << HV << "," << pedestal[i];

      for (int led=0; led<NLEDS+1; led++)
      {
        printf(" %12.2f", npe[i][led]);
        outfile << "," << npe[i][led];
      }
      cout<<endl;
      outfile<< endl;
    }

    outfile.close();

    /////////////////////////////////////////////////////////
#ifndef COMPILE
    return; 
#else
    return 0;
#endif


  }
