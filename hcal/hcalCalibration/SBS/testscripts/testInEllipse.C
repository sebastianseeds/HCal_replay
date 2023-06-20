#include <iostream>
#include "TCanvas.h"
#include "TH2D.h"
#include "TEllipse.h"

bool isPointInEllipse(double x, double y, double centerX, double centerY, double semiMajorAxis, double semiMinorAxis, double rotationAngle) {
  // Apply rotation angle
  double cosAngle = std::cos(rotationAngle);
  double sinAngle = std::sin(rotationAngle);
  double xRot = (x - centerX) * cosAngle + (y - centerY) * sinAngle;
  double yRot = (y - centerY) * cosAngle - (x - centerX) * sinAngle;

  // Check if point is within the ellipse equation
  double result = ((xRot * xRot) / (semiMajorAxis * semiMajorAxis)) + ((yRot * yRot) / (semiMinorAxis * semiMinorAxis));
  return result <= 1.0;
}

int testInEllipse( double pointX, double pointY ) {
  
  gStyle->SetPalette(53);

  int nPoints = 100000; // Number of points in the histogram

  double centerX = 0.0;
  double centerY = 0.0;
  double semiMajorAxis = 1.0;
  double semiMinorAxis = 2.0;
  double rotationAngle = 0.0 * M_PI / 180.0; // Convert angle to radians

  // Create a TH2D histogram with a Gaussian distribution in both x and y
  TH2D *histogram = new TH2D("histogram", "Random Gaussian Distribution", 100, -10, 10, 100, -10, 10);
  histogram->GetXaxis()->SetTitle("X");
  histogram->GetYaxis()->SetTitle("Y");

  // Fill the histogram with random values
  for (int i = 0; i < nPoints; ++i) {
    double x = gRandom->Gaus(0, 1);
    double y = gRandom->Gaus(0, 1);
    histogram->Fill(x, y);
  }

  // Check if the arbitrary point is within the ellipse
  bool isInsideEllipse = isPointInEllipse(pointX, pointY, centerX, centerY, semiMajorAxis, semiMinorAxis, rotationAngle);

  // Create a canvas to display the histogram
  TCanvas *canvas = new TCanvas("canvas", "Histogram Canvas", 800, 600);

  // Draw the histogram
  histogram->Draw("colz");

  // Draw the ellipse on the histogram
  TEllipse *ellipse = new TEllipse(centerX, centerY, semiMajorAxis, semiMinorAxis, 0.0, 360.0, rotationAngle * 180.0 / M_PI);
  ellipse->SetFillStyle(0);
  ellipse->SetLineColor(kRed);
  ellipse->SetLineWidth(2);
  ellipse->Draw();

  // Draw the arbitrary point
  TMarker *point = new TMarker(pointX, pointY, 20);
  point->SetMarkerSize(2);
  if (isInsideEllipse) {
    point->SetMarkerColor(kBlue);
  } else {
    point->SetMarkerColor(kRed);
  }
  point->Draw();

  return 0;

}
