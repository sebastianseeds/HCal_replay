//sseeds 6.2.23 Function to write a canvas to pdf from file. Should be obsolete moving forward as these canvases will be saved automatically by calibration scripts.

#include <vector>
#include <iostream>
#include <algorithm>
#include "TCut.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TTreeFormula.h"
#include "TChain.h"
#include "TLatex.h"
#include "TH1D.h"
#include "TStopwatch.h"
#include "../include/sbs.h"

//Saves all canvas from a file with names as they appear in the file with iterator appended
void save_all_canvas( std::string file_path, std::string out_path ){
  TFile *file = new TFile(file_path.c_str(), "READ"); // Open the ROOT file in read mode

  TCanvas *canvas = nullptr; // Declare a pointer to TCanvas

  TIter next(file->GetListOfKeys());
  TKey *key;
  Int_t failSafe = 0;

  while ((key = static_cast<TKey *>(next()))) {
    TObject *obj = key->ReadObj();

    if (obj->IsA() == TCanvas::Class()) {
      canvas = static_cast<TCanvas *>(obj);
      // Get the name of the canvas
      std::string canvasName = canvas->GetName();
      std::string canvasTrunc = "";

      for( Int_t i=0; i<canvasName.size(); i++ ){
	if( canvasName[i]==',' ) 
	  break;
	canvasTrunc+=canvasName[i];
      }

      std::string outputName = out_path + canvasTrunc;

      std::string canvasPath = Form("%s_%d.pdf", outputName.c_str(), failSafe);

      // Save the canvas as a PDF file
      canvas->SaveAs(canvasPath.c_str());

      if(failSafe>hcal::gNstamp){
	std::cout << "Error: Canvas saves over limit set by maximum number of calibration sets. Check output calibration file." << std::endl;
	return;
      }
      failSafe++;

    }
    delete obj; // Clean up the memory

  }

  file->Close(); // Close the ROOT file
  delete file; // Clean up the memory

}

void quality_canv( const char *experiment = "gmn", Int_t config = 4, Int_t pass = 0, bool h2only = true ){
  
  Int_t hcal_energy_Nsets = 1;

  std::string h2opt = "";
  if( h2only )
    h2opt = "_lh2only";
  std::string plotdir = Form("../quality_plots/%s/conf%d%s/",experiment,config,h2opt.c_str());
  std::string outdir_path = gSystem->Getenv("OUT_DIR");
  std::string adct_type = "adct";
  std::string tdc_type = "tdc";

  std::string adct_align_qr0_path = outdir_path + Form("/hcal_calibrations/pass%d/timing/%salign_class_%s_conf%d_qr0_pass%d.root",pass,adct_type.c_str(),experiment,config,pass);

  std::string tdc_align_qr0_path = outdir_path + Form("/hcal_calibrations/pass%d/timing/%salign_class_%s_conf%d_qr0_pass%d.root",pass,tdc_type.c_str(),experiment,config,pass);

  std::string ecal_qr0_path = outdir_path + Form("/hcal_calibrations/pass%d/energy/ecal%s_%s_conf%d_qr0_pass%d.root",pass,h2opt.c_str(),experiment,config,pass);

  //ADCt canvases
  save_all_canvas(adct_align_qr0_path,plotdir);

  //TDC canvases
  save_all_canvas(tdc_align_qr0_path,plotdir);

  //Energy canvases
  save_all_canvas(ecal_qr0_path,plotdir);

}
