#include "TH2D.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TEllipse.h"
#include "TF2.h"
#include "TMath.h"

int testEllipseInt() {
    int nBinsX = 100;
    int nBinsY = 100;
    double xMin = -3.0;
    double xMax = 3.0;
    double yMin = -3.0;
    double yMax = 3.0;

    // Create the TH2D histogram
    TH2D* histogram = new TH2D("histogram", "Gaussian Distribution", nBinsX, xMin, xMax, nBinsY, yMin, yMax);

    // Generate random numbers following a Gaussian distribution
    TRandom3 random;
    double meanX = 0.5;
    double sigmaX = 0.5;
    double meanY = -0.5;
    double sigmaY = 1.0;
    int nEvents = 10000;

    for (int i = 0; i < nEvents; ++i) {
        double x = random.Gaus(meanX, sigmaX);
        double y = random.Gaus(meanY, sigmaY);
        histogram->Fill(x, y);
    }

    // Adjust histogram limits based on filled data range
    histogram->GetXaxis()->SetRangeUser(histogram->GetXaxis()->GetXmin(), histogram->GetXaxis()->GetXmax());
    histogram->GetYaxis()->SetRangeUser(histogram->GetYaxis()->GetXmin(), histogram->GetYaxis()->GetXmax());

    // Fit an ellipse to the 2-sigma contour
    TF2* ellipseFit = new TF2("ellipseFit", "[0]*TMath::Power((x-[1]), 2) + [2]*TMath::Power((y-[3]), 2) - 1.0", -3.0, 3.0, -3.0, 3.0);
    ellipseFit->SetParameters(1.0, 0.0, 1.0, 0.0);
    histogram->Fit(ellipseFit, "QN");

    // Create a canvas and draw the histogram
    TCanvas* canvas = new TCanvas("canvas", "Gaussian Distribution Canvas", 800, 600);
    histogram->Draw("colz");

    // Draw the fitted ellipse
    TEllipse* ellipse = new TEllipse(ellipseFit->GetParameter(1), ellipseFit->GetParameter(3), TMath::Sqrt(1.0 / ellipseFit->GetParameter(0)), TMath::Sqrt(1.0 / ellipseFit->GetParameter(2)), 0, 360, 0);
    ellipse->SetLineColor(kRed);
    ellipse->SetLineWidth(2);
    ellipse->SetFillStyle(0); // Set fill style to transparent
    ellipse->Draw();

    canvas->SaveAs("output.png");

    return 0;
}
