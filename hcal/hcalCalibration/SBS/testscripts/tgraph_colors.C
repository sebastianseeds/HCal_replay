//sseeds 8.5.23 Test script to overlay a TH2D with two colors of tgraph points representing fit means per bin

#include "../include/sbs.h"
#include <vector>
#include <iostream>
#include <algorithm>

const Double_t E_approx_FWHM = 20.;
const Double_t pC_approx_FWHM = 10.;

///BASE
// void overlayWithGaussianFits(TH2D* hist, TCanvas* canvas, bool pCopt, const std::vector<int>& redBins) {

//   if (!hist || !canvas) {
//     std::cerr << "Null histogram or canvas pointer!" << std::endl;
//     return;
//   }

//   // Use the provided canvas
//   canvas->cd();

//   // Create TGraphErrors for standard and red points
//   TGraphErrors* graph = new TGraphErrors();
//   TGraphErrors* redGraph = new TGraphErrors();

//   for (int binX = 1; binX <= hist->GetNbinsX(); ++binX) {
//     // Project the 2D histogram onto a 1D histogram along the y-axis for the current x-bin
//     TH1D* projY = hist->ProjectionY("_py", binX, binX);

//     // Fit the 1D histogram with a Gaussian
//     projY->Fit("gaus", "Q"); // Q for Quiet mode

//     // Get the fitted function and retrieve mean and std dev
//     TF1* fitFunction = projY->GetFunction("gaus");
//     if (fitFunction) { // Ensure the fit was successful
//       double mean = fitFunction->GetParameter(1);
//       double sigma = fitFunction->GetParameter(2);

//       // Get the errors on the mean and sigma
//       double meanError = fitFunction->GetParError(1);
//       double sigmaError = fitFunction->GetParError(2);

//       // Determine which graph to fill based on whether the bin is in redBins
//       TGraphErrors* currentGraph = graph;
//       if (std::find(redBins.begin(), redBins.end(), binX) != redBins.end()) {
// 	currentGraph = redGraph;
//       }

//       int pointIdx = currentGraph->GetN();
//       currentGraph->SetPoint(pointIdx, binX, mean);
//       currentGraph->SetPointError(pointIdx, 0, sigma);
//     }
//     delete projY; // Clean up

//   }


//   // Draw the original histogram, standard graph, and red graph
//   hist->Draw("COLZ");
//   graph->SetMarkerStyle(20);
//   graph->SetMarkerColor(kBlack);
//   graph->SetLineColor(kBlack);
//   graph->Draw("P SAME");

//   redGraph->SetMarkerStyle(20);
//   redGraph->SetMarkerColor(kRed);
//   redGraph->SetLineColor(kRed);
//   redGraph->Draw("P SAME");

//   // Update the canvas
//   canvas->Update();
// }

////FROM CURRENT SCRIPT
void overlayWithGaussianFits(TH2D* hist, TCanvas* canvas, bool pCopt, const std::vector<int>& redBins) {
  if (!hist || !canvas) {
    std::cerr << "Null histogram or canvas pointer!" << std::endl;
    return;
  }

  // Use the provided canvas
  canvas->cd();

  // Create TGraphErrors for standard and red points
  TGraphErrors* graph = new TGraphErrors();
  TGraphErrors* redGraph = new TGraphErrors();

  for (int binX = 1; binX <= hist->GetNbinsX(); ++binX) {
    // Project the 2D histogram onto a 1D histogram along the y-axis for the current x-bin
    TH1D* projY = hist->ProjectionY("_py", binX, binX);

    //Skip the fit if less than 100 entries exist in the bin
    // if( projY->GetEntries()<100 )
    //   continue;
	
    //declare some dynamic fit variables
    Double_t FWHM = E_approx_FWHM;
    if(pCopt) FWHM = pC_approx_FWHM;
    Int_t binMax = projY->GetMaximumBin();
    Double_t binCenter = projY->GetBinCenter( binMax );
    Double_t fitLowerLim = binCenter - FWHM;
    Double_t fitUpperLim = binCenter + FWHM;
	
    // Fit a Gaussian to this projection
    //TF1 *gaussFit = new TF1("gausFit", "gaus");
    //projY->Fit(gaussFit, "Q", "", fitLowerLim, fitUpperLim ); // "Q" for quiet mode
    //projY->Fit(gaussFit, "Q"); // "Q" for quiet mode

    // // Fit the 1D histogram with a Gaussian
     projY->Fit("gaus", "Q"); // Q for Quiet mode

    // // Get the fitted function and retrieve mean and std dev
     TF1* gaussFit = projY->GetFunction("gaus");


    if (gaussFit) { // Ensure the fit was successful

      // Determine the x-value as the center of the current bin
      double xCenter = hist->GetXaxis()->GetBinCenter(binX);

      double mean = gaussFit->GetParameter(1);
      double sigma = gaussFit->GetParameter(2);

      // Get the errors on the mean and sigma
      double meanError = gaussFit->GetParError(1);
      double sigmaError = gaussFit->GetParError(2);

      // Determine which graph to fill based on whether the bin is in redBins
      TGraphErrors* currentGraph = graph;
      if (std::find(redBins.begin(), redBins.end(), binX) != redBins.end()) {
	currentGraph = redGraph;
      }

      int pointIdx = currentGraph->GetN();
      currentGraph->SetPoint(pointIdx, xCenter, mean);
      currentGraph->SetPointError(pointIdx, 0, sigma);
    }
    delete projY; // Clean up
  }

  // Draw the original histogram, standard graph, and red graph
  hist->Draw("COLZ");
  graph->SetMarkerStyle(20);
  graph->SetMarkerColor(kBlack);
  graph->SetLineColor(kBlack);
  graph->Draw("P SAME");

  redGraph->SetMarkerStyle(20);
  redGraph->SetMarkerColor(kRed);
  redGraph->SetLineColor(kRed);
  redGraph->Draw("P SAME");

  // Update the canvas
  canvas->Update();
}


// void overlayWithGaussianFits(TH2D* h2, TCanvas* canvas, std::vector<double>& specialValues) {
//     // Create TGraphErrors to store mean and std dev
//     TGraphErrors* graph = new TGraphErrors();
//     TGraphErrors* highlightedGraph = new TGraphErrors();

//     int nBinsX = h2->GetNbinsX();
//     for (int binx = 1; binx <= nBinsX; ++binx) {
//         TH1D* projY = h2->ProjectionY("py", binx, binx);
//         // TF1* gaussFit = new TF1("gaussFit", "gaus", projY->GetXaxis()->GetXmin(), projY->GetXaxis()->GetXmax());
        
// 	//Skip the fit if less than 100 entries exist in the bin
// 	if( projY->GetEntries()<100 ){
// 	  //graphErrors->SetPoint(i - 1, h2d->GetXaxis()->GetBinCenter(i), 0);
// 	  //graphErrors->SetPointError(i - 1, 0, 0); // Assuming no error on x
// 	  continue;
// 	}
	
// 	//declare some dynamic fit variables
// 	Double_t FWHM = E_approx_FWHM;
// 	Int_t binMax = projY->GetMaximumBin();
// 	Double_t binCenter = projY->GetBinCenter( binMax );
// 	Double_t fitLowerLim = binCenter - FWHM;
// 	Double_t fitUpperLim = binCenter + FWHM;
	
// 	// Fit a Gaussian to this projection
// 	TF1 *gaussFit = new TF1("gausFit", "gaus");
// 	projY->Fit(gaussFit, "Q", "", fitLowerLim, fitUpperLim ); // "Q" for quiet mode


//         if (projY->GetEntries() > 0) {
//             projY->Fit(gaussFit, "Q");  // Quiet mode
//             double mean = gaussFit->GetParameter(1);
//             double sigma = gaussFit->GetParameter(2);

//             if (std::find(specialValues.begin(), specialValues.end(), h2->GetXaxis()->GetBinCenter(binx)) != specialValues.end()) {
//                 // This x-value is special and should be colored differently
//                 int pointIndex = highlightedGraph->GetN();
//                 highlightedGraph->SetPoint(pointIndex, h2->GetXaxis()->GetBinCenter(binx), mean);
//                 highlightedGraph->SetPointError(pointIndex, 0, sigma);
//             } else {
//                 int pointIndex = graph->GetN();
//                 graph->SetPoint(pointIndex, h2->GetXaxis()->GetBinCenter(binx), mean);
//                 graph->SetPointError(pointIndex, 0, sigma);
//             }
//         }
//         delete gaussFit;
//     }

//     canvas->cd();

//     // Draw the histogram
//     h2->Draw("COLZ");

//     // Draw the graphs
//     graph->SetMarkerStyle(20);
//     graph->SetMarkerColor(kBlack);
//     graph->SetLineColor(kBlack);
//     graph->Draw("P SAME");

//     highlightedGraph->SetMarkerStyle(27);
//     highlightedGraph->SetMarkerColor(kRed);
//     highlightedGraph->SetLineColor(kRed);
//     highlightedGraph->Draw("P SAME");

//     canvas->Update();
// }

// void overlayWithGaussianFits(TH2D* hist, TCanvas* canvas, const std::vector<int>& redBins) {
//   if (!hist || !canvas) {
//     std::cerr << "Null histogram or canvas pointer!" << std::endl;
//     return;
//   }

//   // Use the provided canvas
//   canvas->cd();

//   // Create TGraphErrors for standard and red points
//   TGraphErrors* graph = new TGraphErrors();
//   TGraphErrors* redGraph = new TGraphErrors();

//   for (int binX = 1; binX <= hist->GetNbinsX(); ++binX) {
//     // Project the 2D histogram onto a 1D histogram along the y-axis for the current x-bin
//     TH1D* projY = hist->ProjectionY("_py", binX, binX);

//     // Fit the 1D histogram with a Gaussian
//     projY->Fit("gaus", "Q"); // Q for Quiet mode

//     // Get the fitted function and retrieve mean and std dev
//     TF1* fitFunction = projY->GetFunction("gaus");
//     if (fitFunction) { // Ensure the fit was successful

//       // Determine the x-value as the center of the current bin
//       double xCenter = hist->GetXaxis()->GetBinCenter(binX);

//       double mean = fitFunction->GetParameter(1);
//       double sigma = fitFunction->GetParameter(2);

//       // Get the errors on the mean and sigma
//       double meanError = fitFunction->GetParError(1);
//       double sigmaError = fitFunction->GetParError(2);

//       // Determine which graph to fill based on whether the bin is in redBins
//       TGraphErrors* currentGraph = graph;
//       if (std::find(redBins.begin(), redBins.end(), binX) != redBins.end()) {
// 	currentGraph = redGraph;
//       }

//       int pointIdx = currentGraph->GetN();
//       currentGraph->SetPoint(pointIdx, xCenter, mean);
//       currentGraph->SetPointError(pointIdx, 0, sigma);
//     }
//     delete projY; // Clean up
//   }

//   // Draw the original histogram, standard graph, and red graph
//   hist->Draw("COLZ");
//   graph->SetMarkerStyle(20);
//   graph->SetMarkerColor(kBlack);
//   graph->SetLineColor(kBlack);
//   graph->Draw("P SAME");

//   redGraph->SetMarkerStyle(20);
//   redGraph->SetMarkerColor(kRed);
//   redGraph->SetLineColor(kRed);
//   redGraph->Draw("P SAME");

//   // Update the canvas
//   canvas->Update();
// }

void tgraph_colors(){
  
  // Initialize the random number generator
  TRandom3 randGen;

  // Create a 2D histogram
  TH2D *h2 = new TH2D("h2", "TH2D with Gaussian Distributions in Y;X;Y", 100, 0, 100, 100, 0, 100);

  // For each x-bin, determine a Gaussian distribution in y
  for (int binX = 1; binX <= 100; ++binX) {
    double meanY = randGen.Uniform(1, 100);

    // Populate the y-values for this x-bin using the Gaussian distribution
    for (int i = 0; i < 1000; ++i) { // 1000 entries per bin to make the Gaussian shape clearer
      double valX = binX; // center of the bin in x
      double valY = randGen.Gaus(meanY, 5); // Assuming a fixed standard deviation of 5 for y
      h2->Fill(valX, valY);
    }
  }
  
  std::vector<Int_t> lh2runs = {5,25,88,89,90};

  //overlay energy vs run histo
  TCanvas *c1 = new TCanvas("c1","E vs Run",1200,1200);
  c1->cd();
  overlayWithGaussianFits(h2,c1,false,lh2runs);
  //c1->Write();

}
