#include <iostream>
#include <fstream>
#include <string>

#include "TGraph.h"

using namespace std;

int plotHV3()
{
  const int header = 1; // header lines to skip
  const int npmt = 12 * 24; // number of pmt's
  const float peakAmp = 600; // peak amplitude we want
  const float nevents_threshold = 600; // min number of events to consider the current point
   

  // files to read
  const char* files[]={
   // "run_1196.csv",
    // "run_1197.csv",
    //"run_1198.csv",
    //"run_1200.csv",
   // "run_1201.csv",
    //"run_1203.csv",
   // "run_1204.csv",
    //"run_1205.csv",
   // "run_1206.csv",
      //"run_1286.csv",
      //"run_1288.csv",
      "run_1289.csv",
      "run_1290.csv",
      "run_1291.csv",
      "run_1292.csv",
      "run_1293.csv",
      "run_1294.csv",
      "run_1297.csv",
      "run_1300.csv",
      "run_1302.csv",
      "run_1303.csv",
  };
  const int nfiles = sizeof( files ) / sizeof( *files );

  // variables to read lines of the data file.
  string str;      
  float HVArray  [ nfiles ][ npmt ];
  float pedestal [ nfiles ][ npmt ];
  int channel;
  int ndata[ npmt ];
  const int nleds = 6; // 0..5
  Float_t npe_led0 [nfiles][ npmt ];
  Float_t npe_led1 [nfiles][ npmt ];
  Float_t npe_led2 [nfiles][ npmt ];
  Float_t npe_led3 [nfiles][ npmt ];
  Float_t npe_led4 [nfiles][ npmt ];
  Float_t npe_led5 [nfiles][ npmt ];

  // variables for plots
  float HV      [ nfiles ];
  float npe_0   [ nfiles ];
  float npe_1   [ nfiles ];
  float npe_2   [ nfiles ];
  float npe_3   [ nfiles ];
  float npe_4   [ nfiles ];
  float npe_5   [ nfiles ];
  TGraph *plots  [ npmt ][ nleds ];
  TString title;
  TCanvas *canvas [ npmt ][ nleds ];
  int npoints;
  float minx = 10000;
  float maxx = 0;

  // variables for fit
  TF1 *powerlaw;
  float y[ npmt ]; // HV to set the pmt's in order to have amplitude = peakAmp
  gStyle->SetOptFit(1111); // things to show in the box
  gStyle->SetStatX(.5); // box position on x axis, normalized in [0,1] 
  gStyle->SetStatY(.87); // box position on y axis, normalized in [0,1]  

  cout << "Going to open "<< nfiles <<" files" << endl;

  string newfiles[nfiles];
  string buffer, ch;
  for (int ifile = 0; ifile < nfiles; ifile++)
  {
    ifstream fp ( files[ ifile ] );
    newfiles[ ifile ] = files[ ifile ];
    newfiles[ ifile ] += "_tmp";
    ofstream fop( newfiles[ ifile ] );

    for (int h = 0; h < header + npmt; h++)
    {
      getline(fp, str);
      buffer = "";
      for (int i=0; i < str.length(); i++)
      {
        ch = str.substr( i, 1 );
        if ( ch == "," ) ch = " ";
        buffer += ch;
      }
      fop << buffer << endl;
    }

    fp.close();
    fop.close();
  }

  for (int ifile = 0; ifile < nfiles; ifile++)
  {
    ifstream fp ( newfiles[ ifile ] );
    cout << "Opened file "<< files[ ifile ]<<endl;

    // skip header
    for (int h = 0; h < header; h++)
      getline(fp, str);

    for (int pmt = 0; pmt < npmt; pmt++)
    {
      fp >> channel 
        >> HVArray[ ifile ][ pmt ] 
        >> pedestal[ ifile ][ pmt ]
        >> npe_led0 [ifile][ pmt ]
        >> npe_led1 [ifile][ pmt ]
        >> npe_led2 [ifile][ pmt ]
        >> npe_led3 [ifile][ pmt ]
        >> npe_led4 [ifile][ pmt ]
        >> npe_led5 [ifile][ pmt ];

      /*
      cout << pmt <<" "<< channel <<" "
        << HVArray[ ifile ][ pmt ] <<" "
        << pedestal[ ifile ][ pmt ] <<" "
        << npe_led0 [ifile][ pmt ] <<" "
        << npe_led1 [ifile][ pmt ] <<" "
        << npe_led2 [ifile][ pmt ] <<" "
        << npe_led3 [ifile][ pmt ] <<" "
        << npe_led4 [ifile][ pmt ] <<" "
        << npe_led5 [ifile][ pmt ];
        */
    }

    fp.close();
  }
  system("rm *_tmp");
  system("mkdir -p fit");

  powerlaw = new TF1("powerlaw","[0] * x**[1]", minx - 100, maxx + 100);
  powerlaw->SetParameter(0,100);
  powerlaw->SetParameter(1, 1);
  powerlaw->SetParName(0,"Normalization");
  powerlaw->SetParName(1,"Exponent");

  for ( int pmt = 0; pmt < npmt; pmt++)
  {
    ndata [ pmt ] = 0;
    for ( int index = 0; index < nfiles; index++)
    {
      HV [ ndata [pmt] ] = HVArray[ index ][ pmt ];
      npe_0 [ ndata [pmt] ] = npe_led0[ index ][ pmt ];
      npe_1 [ ndata [pmt] ] = npe_led1[ index ][ pmt ];
      npe_2 [ ndata [pmt] ] = npe_led2[ index ][ pmt ];
      npe_3 [ ndata [pmt] ] = npe_led3[ index ][ pmt ];
      npe_4 [ ndata [pmt] ] = npe_led4[ index ][ pmt ];
      npe_5 [ ndata [pmt] ] = npe_led5[ index ][ pmt ];
      ndata[pmt]++;
    }
    title.Form("PMT_%d_led_0", pmt);
    canvas [ pmt ][0] = new TCanvas(title);

    // log scale?
    //gPad->SetLogx();
    //gPad->SetLogy();

    plots [ pmt ][0] = new TGraph( ndata[pmt], HV, npe_0);
    plots [ pmt ][0]->SetTitle( title );
    plots [ pmt ][0]->GetYaxis()->SetTitle( "# of photoelectrons" );
    plots [ pmt ][0]->GetYaxis()->SetLimits(0,150);
    plots [ pmt ][0]->GetYaxis()->SetMoreLogLabels(kTRUE); // try and add some more labels on log axis
    plots [ pmt ][0]->GetXaxis()->SetTitle( "HV (V)");
    plots [ pmt ][0]->SetMarkerStyle(3);
    plots [ pmt ][0]->SetMarkerSize(3);
    plots [ pmt ][0]->Draw("APL");
    //plots [ pmt ][0]->Fit( "powerlaw","Q" );
    canvas[ pmt ][0]->SaveAs(Form("fit/npeled3/Fit_PMT_%d_led_0.png", pmt));

    title.Form("PMT_%d_led_1", pmt);
    canvas [ pmt ][1] = new TCanvas(title);
    plots [ pmt ][1] = new TGraph( ndata[pmt], HV, npe_1);
    plots [ pmt ][1]->SetTitle( title );
    plots [ pmt ][1]->GetYaxis()->SetTitle( "# of photoelectrons" );
   // plots [ pmt ][1]->GetYaxis()->SetLimits(0,400);
    plots [ pmt ][1]->GetYaxis()->SetMoreLogLabels(kTRUE); // try and add some more labels on log axis
    plots [ pmt ][1]->GetXaxis()->SetTitle( "HV (V)");
    plots [ pmt ][1]->SetMarkerStyle(3);
    plots [ pmt ][1]->SetMarkerSize(3);
    plots [ pmt ][1]->Draw("AP");
    //plots [ pmt ][1]->Fit( "powerlaw","Q" );
    canvas[ pmt ][1]->SaveAs(Form("fit/npeled3/Fit_PMT_%d_led_1.png", pmt));

    title.Form("PMT_%d_led_2", pmt);
    canvas [ pmt ][2] = new TCanvas(title);
    plots [ pmt ][2] = new TGraph( ndata[pmt], HV, npe_2);
    plots [ pmt ][2]->SetTitle( title );
    plots [ pmt ][2]->GetYaxis()->SetTitle( "# of photoelectrons" );
     // plots [ pmt ][2]->GetYaxis()->SetLimits(0,1000);
    plots [ pmt ][2]->GetYaxis()->SetMoreLogLabels(kTRUE); // try and add some more labels on log axis
    plots [ pmt ][2]->GetXaxis()->SetTitle( "HV (V)");
    plots [ pmt ][2]->SetMarkerStyle(3);
    plots [ pmt ][2]->SetMarkerSize(3);
    plots [ pmt ][2]->Draw("AP");
   // plots [ pmt ][2]->Fit( "powerlaw","Q" );
    canvas[ pmt ][2]->SaveAs(Form("fit/npeled3/Fit_PMT_%d_led_2.png", pmt));

    title.Form("PMT_%d_led_3", pmt);
    canvas [ pmt ][3] = new TCanvas(title);
    plots [ pmt ][3] = new TGraph( ndata[pmt], HV, npe_3);
    plots [ pmt ][3]->SetTitle( title );
    plots [ pmt ][3]->GetYaxis()->SetTitle( "# of photoelectrons" );
     // plots [ pmt ][3]->GetYaxis()->SetLimits(0,2000);
    plots [ pmt ][3]->GetYaxis()->SetMoreLogLabels(kTRUE); // try and add some more labels on log axis
    plots [ pmt ][3]->GetXaxis()->SetTitle( "HV (V)");
    plots [ pmt ][3]->SetMarkerStyle(3);
    plots [ pmt ][3]->SetMarkerSize(3);
    plots [ pmt ][3]->Draw("AP");
    //plots [ pmt ][3]->Fit( "powerlaw","Q" );
    canvas[ pmt ][3]->SaveAs(Form("fit/npeled3/Fit_PMT_%d_led_3.png", pmt));

    title.Form("PMT_%d_led_4", pmt);
    canvas [ pmt ][4] = new TCanvas(title);
    plots [ pmt ][4] = new TGraph( ndata[pmt], HV, npe_4);
    plots [ pmt ][4]->SetTitle( title );
    plots [ pmt ][4]->GetYaxis()->SetTitle( "# of photoelectrons" );
     // plots [ pmt ][4]->GetYaxis()->SetLimits(0,3000);
    plots [ pmt ][4]->GetYaxis()->SetMoreLogLabels(kTRUE); // try and add some more labels on log axis
    plots [ pmt ][4]->GetXaxis()->SetTitle( "HV (V)");
    plots [ pmt ][4]->SetMarkerStyle(3);
    plots [ pmt ][4]->SetMarkerSize(3);
    plots [ pmt ][4]->Draw("AP");
    //plots [ pmt ][4]->Fit( "powerlaw","Q" );
    canvas[ pmt ][4]->SaveAs(Form("fit/npeled3/Fit_PMT_%d_led_4.png", pmt));

    title.Form("PMT_%d_led_5", pmt);
    canvas [ pmt ][5] = new TCanvas(title);
    plots [ pmt ][5] = new TGraph( ndata[pmt], HV, npe_5);
    plots [ pmt ][5]->SetTitle( title );
    plots [ pmt ][5]->GetYaxis()->SetTitle( "# of photoelectrons" );
    //  plots [ pmt ][5]->GetYaxis()->SetLimits(0,5000);
    plots [ pmt ][5]->GetYaxis()->SetMoreLogLabels(kTRUE); // try and add some more labels on log axis
    plots [ pmt ][5]->GetXaxis()->SetTitle( "HV (V)");
    plots [ pmt ][5]->SetMarkerStyle(3);
    plots [ pmt ][5]->SetMarkerSize(3);
    plots [ pmt ][5]->Draw("AP");
    //plots [ pmt ][5]->Fit( "powerlaw","Q" );
    canvas[ pmt ][5]->SaveAs(Form("fit/npeled3/Fit_PMT_%d_led_5.png", pmt));
  }

  return 0;
}
