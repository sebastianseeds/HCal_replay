//SSeeds 5.6.24 - Production update to use alphas to predict new HV settings based on current HV settings and desired factor

#include <ctime>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <TH2F.h>
#include <TChain.h>
#include <TCanvas.h>
#include <iostream>
#include <TSystem.h>
#include <TStopwatch.h>
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TString.h"
#include "TCut.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TVector3.h"
#include "TGraphErrors.h"
#include "hcal.h"

const int kNcell = 288; // Total number of HCal modules
const int kNrows = 24; // Total number of HCal rows
const int kNcols = 12; // Total number of HCal columns

const double HV_hard_lim = 2550;

//main. gain_factor is desired factor increase in ADC gain
int HV_targ( double gain_factor=2. ){
  
  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  /////
  //CONFIGURE PATHS
  string HVfilePath = "oldHV.txt"; //Current HV settins (one per line in cell order)
  string alphasFilePath = "alphas.txt"; //alpha parameters (one per line in cell order)
  string inConstPath = "hcal_gainCoeff_102121.txt"; //not necessary for relative HV targets
  string HVoutPath = "HVTargets.txt"; //Output path
  string HVoutPath_safe= "HVTargets_safe.txt"; //Output path
  /////

  // Declare outfile
  TFile *fout = new TFile( "HVTarg.root", "RECREATE" );
  
  // Initialize vectors and arrays
  double gHV[kNrows][kNcols];
  double gAlphas[kNrows][kNcols];
  double gTargetHV[kNrows][kNcols];
  //double gOldConst[kNcell];
  std::vector<double> gOldConst; //unnecessary for relative increases using gain_factor
  std::vector<double> maxHV;
  std::vector<double> minHV;

  // Read in HV settings for PMTs
  ifstream HVFile( HVfilePath );
  
  if( !HVFile ){
    cerr << "No HV settings present -> setFiles/HV_run" << endl;
    return 0;
  }
  
  cout << "Getting HV settings for each pmt.." << endl;
  int n1=0;
  double d1;
  int rval, cval;
  string line;
  
  while( getline( HVFile, line ) ){
    if( line.at(0) == '#' ) {
      continue;
    }
    stringstream ss( line );
    ss >> d1;
    rval = floor( n1/kNcols );
    cval = n1 % kNcols;
    gHV[rval][cval] = -d1; 

    cout << "HV" << n1 << " = " << -d1 << "." << endl;

    n1++;
  }
  
  // Read in alpha parameters for PMTs
  ifstream alphaFile( alphasFilePath );
  if( !alphaFile ){
    cerr << "No PMT alphas file present -> setFiles/alphas.txt expected." << endl;
    return 0;
  }

  cout << "Getting alpha parameters for each pmt.." << endl;
  n1=0;
  d1=0;
  string line2;
  
  while( getline( alphaFile, line2 ) ){
    if( line2.at( 0 )=='#' ) {
      continue;
    }
    stringstream ss( line2 );
    ss >> d1;
    rval = floor( n1/kNcols );
    cval = n1 % kNcols;
    gAlphas[rval][cval] = d1;

    cout << "alpha" << n1 << " = " << d1 << "." << endl;

    n1++;
  }
 
  // Read in previous constants for PMTs
  ifstream inConstFile( inConstPath );
  bool constfile = true;
  if( !inConstFile ){
    cerr << "No input constant file present -> not necessary for analysis." << endl;
    constfile = false;
  }

  cout << "Getting previous calibration constants.." << endl;
  n1=0;
  d1=0;
  string line4;
  
  while (std::getline(inConstFile, line)) {
    // Skip lines starting with '#'
    if (line.empty() || line.at(0) == '#' || !constfile) {
      continue;
    }

    // Read values from the line
    std::stringstream ss(line);
    double d;
    while (ss >> d) {
      gOldConst.push_back(d);
      cout << "const" << gOldConst.size() << " = " << d << "." << endl;
    }
  }

  // Read max HV values from plateau analysis
  std::ifstream file("setFiles/Vplat_l1.txt");
  if (!file) {
    std::cerr << "Error: Could not open file Vplat_l1.txt." << std::endl;
    return 0;
  }

  std::string line5;
  while (std::getline(file, line5)) {
    std::stringstream ss(line5);
    int first, second;
    char comma;

    // Read the line5 in the format "first,second"
    if (ss >> first >> comma >> second) {
      minHV.push_back(first);
      maxHV.push_back(second);

      //cout << "max plateau HV cell" << maxHV.size() << " = " << second << "." << endl;
    }
  }

  // cout << "!!!!!!!!!!!!!!!!!!!!" << endl;
  // for( size_t i=0; i<maxHV.size(); ++i )
  //   cout << maxHV[i] << endl;

  cout << "All parameters loaded." << endl;

  cout << "Printing old HVs to console..." << endl;

  // Print old HVs to console
  int cell = 0;
  for( int r = 0; r<kNrows; r++){
    for( int c = 0; c<kNcols; c++){
      if(gHV[r][c]>maxHV[cell])
	cout << -gHV[r][c] << "!!  ";
      else
	cout << -gHV[r][c] << "  ";      
      cell++;
    }
    cout << endl;
  }

  
  cout << endl << "Printing new HVs to console..." << endl;
  // Print to console
  cell = 0;
  for( int r = 0; r<kNrows; r++){
    for( int c = 0; c<kNcols; c++){
      //cout << gHV[r][c]*pow(gain_factor*gOldConst[cell],(1/gAlphas[r][c])) << "  ";
      double projHV = gHV[r][c]*pow(gain_factor,(1/gAlphas[r][c]));

      if(projHV>HV_hard_lim)
	projHV=HV_hard_lim; //prevent damage to PMTs

      if(projHV>maxHV[cell])
	cout << -projHV << "!!  ";
      else
	cout << -projHV << "  ";

      cell++;
    }
    cout << endl;
  }

  vector<double> targHV;
  vector<double> targHV_outPlat;
  vector<double> oldHV;
  
  vector<double> safeHV;

  // Declare outfiles
  ofstream HCal_HV;

  HCal_HV.open( HVoutPath );
  HCal_HV << "#Generated from " << HVfilePath.c_str() << " and alphas.txt with max HV limited by " << HV_hard_lim << endl;
  HCal_HV << "#Module" << '\t' << " HV" << endl;

  //Adding minus signs to match HV setting input and proper voltage setting
  cell = 0;
  for( int i=0; i<kNcell; i++ ){   
    int r = i/kNcols;
    int c = i%kNcols;

    //gTargetHV[r][c] = -gHV[r][c]*pow(gOldConst[i],(1/gAlphas[r][c]));
    gTargetHV[r][c] = gHV[r][c]*pow(gain_factor,(1/gAlphas[r][c]));

    if(gTargetHV[r][c]>HV_hard_lim)
      gTargetHV[r][c]=HV_hard_lim; //prevent damage to PMTs

    //Build vector of out of plateau HV region projections for plot
    if(gTargetHV[r][c]>maxHV[cell]){
      targHV_outPlat.push_back(-gTargetHV[r][c]);
      safeHV.push_back(-maxHV[cell]);
    }else{
      targHV_outPlat.push_back(100); //arbitrarily out of window
      safeHV.push_back(-gTargetHV[r][c]);
    }

    cout << safeHV.back() << " = " << -gTargetHV[r][c] << " or " << -maxHV[cell] << endl;

    targHV.push_back(-gTargetHV[r][c]);
    oldHV.push_back(-gHV[r][c]);

    HCal_HV << r*kNcols+c+1 << '\t' << -gTargetHV[r][c] << endl;
  }

  TGraph *g1 = new TGraph(targHV.size());
  TGraph *g2 = new TGraph(oldHV.size());
  TGraph *g3 = new TGraph(targHV_outPlat.size());
  TGraph *g4 = new TGraph(safeHV.size());
  for (int i = 0; i < targHV.size(); ++i) {
    g1->SetPoint(i, i, targHV[i]);
    g2->SetPoint(i, i, oldHV[i]);
    g3->SetPoint(i, i, targHV_outPlat[i]);
    g4->SetPoint(i, i, safeHV[i]);

    //cout << safeHV[i] << " ";

  }

  g1->SetMarkerStyle(20);
  g1->SetMarkerColor(kBlue);
  g1->SetLineColor(kBlue);
    
  g2->SetMarkerStyle(21);
  g2->SetMarkerColor(kBlack);
  g2->SetLineColor(kBlack);
    
  g3->SetMarkerStyle(29);
  g3->SetMarkerColor(kMagenta);
  g3->SetMarkerSize(2);
  g3->SetLineColor(kMagenta);

  g1->GetXaxis()->SetTitle("Cell");

  TCanvas *canvas = new TCanvas("canvas", "HV Targets", 2100, 800);
  g1->SetTitle(Form("Target HVs/cell %0.1fx gain relative to %s",gain_factor,HVfilePath.c_str()));
  g1->GetYaxis()->SetRangeUser(-2800,0);
  g1->Draw("AP");
  g2->Draw("P SAME");
  g3->Draw("P SAME");

  // Create a horizontal dotted line at HV_hard_lim
  TLine *hline = new TLine(0, -HV_hard_lim, targHV.size(), -HV_hard_lim);
  hline->SetLineStyle(2);  // Dotted line
  hline->SetLineWidth(2);
  hline->SetLineColor(kRed);
  hline->Draw("same");

  // Create a legend
  TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
  legend->AddEntry(g1, "target HV", "p");
  legend->AddEntry(g3, "target HV, Out of Plateau Region", "p");
  legend->AddEntry(g2, "old HV", "p");
  legend->AddEntry(hline, "HV Hard Limit", "l");
  legend->Draw("same");

  canvas->Write();
  HCal_HV.close();

  //Write out new canvas for safe HVs
  g4->SetMarkerStyle(20);
  g4->SetMarkerColor(kGreen);
  g4->SetLineColor(kGreen);

  TCanvas *c2 = new TCanvas("c2", "HV Safe Targets", 2100, 800);
  c2->cd();
  g4->SetTitle(Form("Target HVs/cell %0.1fx gain relative to %s (in plateau)",gain_factor,HVfilePath.c_str()));
  g4->GetYaxis()->SetRangeUser(-2800,0);
  g4->Draw("AP");
  g2->Draw("P SAME");

  hline->Draw("same");

  // Create a legend
  TLegend *leg2 = new TLegend(0.7, 0.7, 0.9, 0.9);
  leg2->AddEntry(g4, "safe target HV", "p");
  leg2->AddEntry(g2, "old HV", "p");
  leg2->AddEntry(hline, "HV Hard Limit", "l");
  leg2->Draw("same");

  c2->Write();

  // Declare safe outfile
  ofstream HCal_HV_safe;

  HCal_HV_safe.open( HVoutPath_safe );
  HCal_HV_safe << "#Generated from " << HVfilePath.c_str() << " and alphas.txt with max HV limits by maximum in plateau region and max HV hard limit " << HV_hard_lim << endl;
  HCal_HV_safe << "#Module" << '\t' << " HV" << endl;

  for( int i=0; i<kNcell; i++ ){   
    int r = i/kNcols;
    int c = i%kNcols;

    HCal_HV_safe << r*kNcols+c+1 << '\t' << safeHV[i] << endl;
  }

  HCal_HV_safe.close();

  fout->Write();

  cout << "General HV targets written to file " << HVoutPath << "." << endl;
  cout << "Safe HV targets written to file " << HVoutPath_safe << "." << endl;

  st->Stop();

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

  return 0;

}


