#include <TCanvas.h>
#include <TText.h>
#include "../include/sbs.h"

const int linecount = 22;

void testReport(const char *experiment = "gmn", Int_t config=4, Int_t pass = 0, bool h2only = false) {

  // Get the date
  string date = util::getDate();

  Int_t Ncal_set_size = 1;
  calset gain_coeff[hcal::gNstamp];

  //Add output report canvas
  TCanvas *c1[Ncal_set_size];

  for( Int_t s=0; s<Ncal_set_size; s++ ){

    //c1[s] = new TCanvas(Form("c%d",s), Form("Configuration/Cut Information, Calibration Set %d",s), 200, 10, 1900, 1000);
    c1[s] = new TCanvas(Form("c%d",s), Form("Configuration/Cut Information, Calibration Set %d",s), 200, 10, 900, 500);
  
    // Set margin.
    c1[s]->SetLeftMargin(0.05);

    // Create a TText object.
    TText *t = new TText();

    // Set text align to the left (horizontal alignment = 1).
    t->SetTextAlign(11);

    //make an array of strings
    std::string target_option = "All Available";
    if( h2only )
      target_option = "LH2";
    std::string report[linecount] = {
      "General Info",
      Form("Experiment: %s, Configuration: %d, Pass: %d", experiment, config, pass),
      Form("Creation Date: %s", date.c_str() ),
      Form("Target(s) Used: %s", target_option.c_str() ),
      Form("Calibration Set: %s", gain_coeff[s].timestamp.c_str() ),
      "",
      "Elastic Cuts",
      Form("Global Elastic Cuts: %s", gain_coeff[s].gcut.c_str() ),
      Form("W2 mean: %f", gain_coeff[s].W2mean),
      Form("W2 sigma: %f", gain_coeff[s].W2sigma),
      Form("dx mean, neutron : %f", gain_coeff[s].dxmean_n),
      Form("dx mean, proton : %f", gain_coeff[s].dxmean_p),
      Form("dx sigma, neutron : %f", gain_coeff[s].dxsigma_p),
      Form("dx sigma, proton : %f", gain_coeff[s].dxsigma_p),
      Form("dy mean : %f", gain_coeff[s].dymean),
      Form("dy sigma : %f", gain_coeff[s].dysigma),
      Form("adc time mean : %f", gain_coeff[s].atimemean),
      Form("adc time sigma : %f", gain_coeff[s].atimesigma),
      "",
      "Other Cuts",
      Form("Minimum Ev per Cell : %d", gain_coeff[s].minEv),
      Form("Minimum Energy Deposited in Cell (factor, vs expectation) : %f", gain_coeff[s].highdelta)
    };
    // Loop to write the lines to the canvas.
    Double_t spacefactor = 1/(double)linecount;
    for(int i = 0; i < linecount; i++) {
      // Vertical position adjusted according to line number.
      double verticalPosition = 0.9 - i * 0.04;
      t->DrawTextNDC(0.1, verticalPosition, report[i].c_str());
    }

  }

}
