//Code to draw pmt tgraphs, average int charge vs HV setting sseeds

#include <TH2F.h>
#include <TChain.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <TSystem.h>
#include "hcal.h"

const Int_t kNrows = 12; //Number of pmt rows
const Int_t kNcols = 12; //Number of pmt columns

TFile *CosmicCalPlots = new TFile("CosmicCalPlots.root","RECREATE");
TGraph *G1 = new TGraph(); //Single graph is filled then overwritten after writing to canvas for each pmt

TCanvas *C1 = new TCanvas("C1","PMT 0-71 intPulse",900,700); //First 12 x 6 grid
TCanvas *C2 = new TCanvas("C2","PMT 72-143 intPulse",900,700); //Second 12 x 6 grid

void plotResults(int run1 = 978, int run2 = 980, int run3 = 987, int run4 = 988, int run5 = 989){ //Can pass 5 HV runs at a time.
  
  int runs[5] = {run1, run2, run3, run4, run5}; //Configure runs for loop
  
  G1->GetXaxis()->SetTitle("HV Setting (units)"); //Same for all
  G1->GetYaxis()->SetTitle("Average Integrated Charge (units)"); //Same for all
  
  C1->Divide(kNcols,ceil(kNrows/2)); //Divide up the canvas into two
  C2->Divide(kNcols,floor(kNrows/2));

  cout << ceil(kNrows/2) << "  " << floor(kNrows/2) << endl;

  //Load in run/HV correspondance
  ifstream inFile1("runHV.txt"); //Read in file for first run
  cout << "Reading in runHV.txt..." << endl;

  

  
  for (int R=0; R<5; R++){

    


    
    ifstream inFile1(Form("intPulseRun_%d.txt",runs[R])); //Read in file for first run

    cout << "Reading in " << Form("intPulseRun_%d.txt",runs[R]) << "..." << endl;
    
    //intPulseRun must be of the format <row  col  intVal  evTot> after lines which start with #

    int row, col, evTot;
    double intVal;
    string line;
    
    while (getline(inFile1,line))
      {
	
	if(line.at(0) == '#'){
	  continue;
	}
	
	stringstream ss(line);
	ss >> row >> col >> intVal >> evTot;
	
	if(row<ceil(kNrows/2)) {
	  C1->cd(kNrows*row+col+1);
	  cout << "Draw at position " << kNrows*row+col+1 << "." << endl;
	} else {
	  C2->cd(kNrows*(row-ceil(kNrows/2))+col+1);
	  cout << "Draw at position " << kNrows*(row-ceil(kNrows/2))+col+1 << "." << endl;	  
	}
	
	G1->SetTitle(Form("IntAve Run%d Row%d Col%d",runs[R],row,col));
	cout << Form("IntAve Run%d Row%d Col%d",runs[R],row,col) << endl;
	G1->SetMinimum(0.0);
	G1->SetMaximum(5000.0);
	G1->SetPoint(R,runs[R],intVal);
	//cout << R << "  " << runs[R] << "  " << intVal << endl;
	G1->SetMarkerStyle(119);
	G1->SetMarkerSize(3);
	G1->SetMarkerColor(2);
	G1->Draw("AP same");

	CosmicCalPlots->cd();
	G1->Write(Form("Run%d,Row%d,Col%d",runs[R],row,col));	
	G1->Draw("AP");
	
	//cout << "row = " << row << ", col = " << col << ", intVal = " << intVal << endl;
      }

    
    //Plot the distribution of pedestals by HV run
    ifstream inFile1(Form("Pedestals_%d.txt",runs[R])); //Read in file for first run

    
  }
}
