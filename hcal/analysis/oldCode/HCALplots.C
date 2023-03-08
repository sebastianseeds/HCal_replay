//Second attempt at plotting ave int charge vs HV setting - sseeds 11.12.20

#include <TH2F.h>
#include <TChain.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <TSystem.h>
#include "hcal.h"

const int kNrows = 12; //Number of pmt rows
const int kNcols = 12; //Number of pmt columns
const int HVt = 5;
const int PMTt = 144;

TFile *CosmicCalPlots = new TFile("CosmicCalPlots.root","RECREATE");
TGraph *G1 = new TGraph(); //Single graph is filled then overwritten after writing to canvas for each pmt for int val
TGraph *G2 = new TGraph(); //For max val

TCanvas *C1 = new TCanvas("C1","PMT intPulse",900,700); //First 12 x 6 grid
TCanvas *C2 = new TCanvas("C2","PMT intPulse",900,700); //Temporary for max val graphs

TF1 *pFit = new TF1("pFit","[0]*pow(x,[1])",0,2000); //Produces remarkably bad fits
TF1 *testFit;
//TF1 *pFit;

TF1 *pFitted;

void HCALplots(int run1 = 978, int run2 = 980, int run3 = 987, int run4 = 988, int run5 = 989){ //Can pass 5 HV runs at a time.

  //int HVt = 5;
  
  int runs[5] = {run1, run2, run3, run4, run5}; //Configure runs for loop
  
  G1->GetXaxis()->SetTitle("HV Setting (units)"); //Same for all
  G1->GetYaxis()->SetTitle("Average Integrated Charge (units)"); //Same for all

  //Load in run/HV correspondance
  int lineNumberA = 0;
  ifstream inFile1("runHV.txt"); //Read in file for first run
  cout << "Reading in runHV.txt.." << endl;

  int run_N;
  double HV;
  string lineA;
  
  vector<int> run_Nv = {0};
  vector<double> HVv = {0};
  
  while (getline(inFile1,lineA))
    {
      //Ignore lines that start with #
      if(lineA.at(0) == '#'){
	continue;
      }
      
      lineNumberA++;
      stringstream ss(lineA);
      ss >> run_N >> HV;
      
      run_Nv.push_back(run_N);
      HVv.push_back(HV);
      
    }  

  //Build array of HV values corresponding to run number
  double HVs[HVt];
  
  for(int i=0; i<HVt; i++){
    for(int ii=0; ii<HVv.size(); ii++){
      if(run_Nv[ii] == runs[i]){
	HVs[i] = HVv[ii]; 
	cout << "HV = " << HVs[i] << " for run = " << run_Nv[ii] << endl; //Each HV i value per run number i
      }
    }
  }

  cout << "HV settings correlated to run number.." << endl;
  
  //Build tgraphs of pedestal values by HV setting
  double pedestals[HVt][kNrows][kNcols];
  int pedNumber = 0;
  
  for (int i=0; i<HVt; i++){

    ifstream pedFile(Form("Pedestals_%d.txt",runs[i]));

    cout << "Getting pedestal values from run " << runs[i] << " at HV " << HVs[i] << "." << endl;
  
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
	
	pedestals[i][rval][cval] = d1;
	//cout << "run = " << i << "row = " << rval << ", col = " << cval << ", array val = " << pedestals[rval][cval] << endl;

	pedNumber++;
      }
  }

  cout << "Pedestals collected.." << endl;

  //Loop to write pedestal tgraphs to file by PMT
  vector<TGraph*> pGv;
  vector<TGraph*> pGv2;
  double ped[HVt];
  double ped2[PMTt];
  
  for(int R=0; R<kNrows; R++){
    for(int C=0; C<kNcols; C++){
      for(int HV=0; HV<HVt; HV++){
	ped[HV]=pedestals[HV][R][C];
      }
      pGv.push_back(new TGraph(HVt,HVs,ped));
      CosmicCalPlots->cd();
      pGv[R*12+C]->SetMarkerStyle(8);
      pGv[R*12+C]->SetMarkerSize(2);
      pGv[R*12+C]->SetLineColor(kWhite);

      pGv[R*12+C]->SetTitle(Form("Pedestals R%d C%d",R,C));
      pGv[R*12+C]->SetTitle(Form("Pedestals R%d C%d",R,C));
      
      pGv[R*12+C]->Write(Form("Pedestals Row %d, Col %d",R,C));
      pGv[R*12+C]->Draw("AP");
    }
  }

  //Loop to write pedestal tgraphs to file by HV value
  int index = 0;
  double PMTindex[PMTt];

  for(int HV=0; HV<HVt; HV++){
    for(int R=0; R<kNrows; R++){
      for(int C=0; C<kNcols; C++){
      	PMTindex[R*12+C]=index;
	index++;
	ped2[R*12+C]=pedestals[HV][R][C];

	cout << "HV " << HVs[HV] << " pedestal value for index " << PMTindex[R*12+C] << " is " << pedestals[HV][R][C] << endl;
	
      }

    }
    pGv2.push_back(new TGraph(PMTt,PMTindex,ped2));
    CosmicCalPlots->cd();
    pGv2[HV]->SetMarkerStyle(8);
    pGv2[HV]->SetMarkerSize(2);
    pGv2[HV]->SetLineColor(kWhite);

    pGv2[HV]->SetTitle(Form("Pedestals HV setting -%f",HVs[HV]));
      
    pGv2[HV]->Write(Form("Pedestals HV setting -%f",HVs[HV]));
    pGv2[HV]->Draw("AP");

    index = 0;
  }
  
  cout << "Pedestal TGraphs constructed.." << endl;

  //Loop to extract integrated pulse values from txt file
  int row, col, evTot;
  double intVal;
  double maxVal;
  string lineB;

  vector<double> intVals = {0};
  vector<double> maxVals = {0};
  vector<int> rowVals = {0};
  vector<int> colVals = {0};
  vector<int> HVVals = {0};
  
  for (int R=0; R<5; R++){

    int lineNumberB = 0;
    ifstream inFile1(Form("avePulseRun_%d.txt",runs[R])); //Read in file for first run (old version reads in intPulseRun_%d.txt)

    cout << "Reading in " << Form("avePulseRun_%d.txt",runs[R]) << "..." << endl;
    
    //avePulseRun must be of the format <row  col  intVal  maxVal  evTot> after lines which start with #
    while (getline(inFile1,lineB))
      {
	
	if(lineB.at(0) == '#'){
	  continue;
	}

	lineNumberB++;
	stringstream ss(lineB);
	ss >> row >> col >> intVal >> maxVal >> evTot;

	intVals.push_back(intVal);
	maxVals.push_back(maxVal);
       	rowVals.push_back(row);
       	colVals.push_back(col);
       	HVVals.push_back(HVs[R]);

      }            
  }

  cout << "Integrated pulse values collected.." << endl;

  int HVct = 0;
  vector<double> intFitC = {0};
  vector<double> intFitP = {0};
  vector<double> maxFitC = {0};
  vector<double> maxFitP = {0};

  //Loop to write all integrated pulse tgraphs to file

  for (int cL=0; cL<(kNcols*kNrows); cL++){
    C1->Clear();
    C1->cd();
    HVct = 0;
    
    for (int r=0; r<kNrows; r++){
      for (int c=0; c<kNcols; c++){
	for (int i=1; i<intVals.size()+1; i++){
	  if (rowVals[i]*kNcols+colVals[i]==cL){
    
	    G1->SetPoint(HVct, HVVals[i], intVals[i]);
	    G1->SetTitle(Form("IntAve Row%d Col%d",r,c));
	    
	    G1->SetMinimum(0.0);
	    G1->SetMarkerStyle(8);
	    G1->SetMarkerSize(2);
	    
	    HVct++;
	    
	  }
	  
	}
      }
    }
    
    pFit->SetParameter(1,9);
    
    G1->Fit("pFit", "Q"); 
    
    intFitC.push_back(pFit->GetParameter(0));
    intFitP.push_back(pFit->GetParameter(1));
    
    int rw = floor(cL/kNcols);
    int cl = cL%kNrows;
    
    CosmicCalPlots->cd();
    G1->SetMinimum(0.0);
    G1->SetMarkerStyle(8);
    G1->SetMarkerSize(2);
    G1->SetLineColor(kWhite);
    G1->SetTitle(Form("IntAve Row%d Col%d",rw,cl));
    G1->Write(Form("Row %d,Col %d",rw,cl));
    G1->Draw("AP");
    
  }
  
  cout << "Integrated pulse tgraphs written to file CosmicCalPlots.root." << endl;
  
  for (int cL=0; cL<(kNcols*kNrows); cL++){
    C2->Clear();
    C2->cd();
    HVct = 0;
    
    for (int r=0; r<kNrows; r++){
      for (int c=0; c<kNcols; c++){
	for (int i=1; i<maxVals.size()+1; i++){
	  if (rowVals[i]*kNcols+colVals[i]==cL){
    
	    G2->SetPoint(HVct, HVVals[i], maxVals[i]);
	    G2->SetTitle(Form("MaxAve Row%d Col%d",r,c));
	    
	    G2->SetMinimum(0.0);
	    G2->SetMarkerStyle(8);
	    G2->SetMarkerSize(2);
	    
	    HVct++;
	    
	  }
	  
	}
      }
    }
    
    pFit->SetParameter(1,9);
    
    G2->Fit("pFit", "Q"); 
    
    maxFitC.push_back(pFit->GetParameter(0));
    maxFitP.push_back(pFit->GetParameter(1));
    
    int rw = floor(cL/kNcols);
    int cl = cL%kNrows;
    
    CosmicCalPlots->cd();
    G2->SetMinimum(0.0);
    G2->SetMarkerStyle(8);
    G2->SetMarkerSize(2);
    G2->SetLineColor(kWhite);
    G2->SetTitle(Form("MaxAve Row%d Col%d",rw,cl));
    G2->Write(Form("Row %d,Col %d",rw,cl));
    G2->Draw("AP");
    
  }

  cout << "Pulse Maximum tgraphs written to file CosmicCalPlots.root." << endl << endl;

  
  cout << "Fit parameters collected.." << endl;

  //Loop to write fit parameters to txt file
  ofstream outFile;
  outFile.open("PMTFitParams.txt");
  
  outFile << "#Fit parameters for PMT by row and col. Fit Function C[HV]^P." << endl;
  outFile << "# Row Column intC intP maxC maxP" << endl;
  for(int i=0; i<intFitC.size(); i++){
    outFile << floor(i/kNrows) << "  " << i % kNcols << "  " << intFitC[i] << "  " << intFitP[i] << "  " << maxFitC[i] << "  " << maxFitP[i] << endl;
  }
  
  cout << "Extraction complete. Results written to CosmicCalPlots.root. Fit parameters written to PMTFitParams.txt." << endl;

}

