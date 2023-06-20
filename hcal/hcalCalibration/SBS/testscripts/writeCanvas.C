//sseeds 6.2.23 Test function to write a canvas to png from file

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


// //Saves all canvas from a file with names as they appear in the file
// void save_all_canvas( std::string file_path, std::string out_path, Int_t set ){
//   TFile *file = new TFile(file_path.c_str(), "READ"); // Open the ROOT file in read mode

//   TCanvas *canvas = nullptr; // Declare a pointer to TCanvas

//   TIter next(file->GetListOfKeys());
//   TKey *key;
//   Int_t failSafe = 0;

//   while ((key = static_cast<TKey *>(next()))) {
//     TObject *obj = key->ReadObj();

//     if (obj->IsA() == TCanvas::Class()) {
//       canvas = static_cast<TCanvas *>(obj);
//       // Get the name of the canvas
//       std::string canvasName = canvas->GetName();
//       std::string canvasTrunc = "";

//       for( Int_t i=0; i<canvasName.size(); i++ ){
// 	if( canvasName[i]==' ' )
// 	  continue;
// 	if( canvasName[i]==',' ) 
// 	  break;
// 	canvasTrunc+=canvasName[i];
//       }

//       std::string outputName = out_path + canvasTrunc;

//       std::string canvasPath = Form("%s_set%d.pdf", outputName.c_str(), set);

//       canvas->SaveAs(canvasPath.c_str());

//       if(failSafe>hcal::gNstamp){
// 	std::cout << "Error: Canvas saves over limit set by maximum number of calibration sets. Check output calibration file." << std::endl;
// 	return;
//       }
//       failSafe++;

//     }
//     delete obj; // Clean up the memory
//   }

//   file->Close(); // Close the ROOT file
//   delete file; // Clean up the memory

// }

void writeCanvas( const char *experiment = "gmn", Int_t config = 4, Int_t pass = 0 ){
  
  Int_t hcal_energy_Nsets = 1;

  std::string plotdir = "test_";
  std::string outdir_path = gSystem->Getenv("OUT_DIR");
  std::string adct_type = "adct";
  std::string tdc_type = "tdc";

  std::string adct_align_qr0_path = outdir_path + Form("/hcal_calibrations/pass%d/timing/%salign_class_%s_conf%d_qr0_pass%d.root",pass,adct_type.c_str(),experiment,config,pass);

  std::string tdc_align_qr0_path = outdir_path + Form("/hcal_calibrations/pass%d/timing/%salign_class_%s_conf%d_qr0_pass%d.root",pass,tdc_type.c_str(),experiment,config,pass);

  //ADCt canvases
  //for( Int_t set=0; set<hcal_energy_Nsets; set++ )
  save_all_canvas(adct_align_qr0_path,plotdir);

  save_all_canvas(tdc_align_qr0_path,plotdir);

}


// void writeCanvas(){

//   TFile *file = new TFile("/volatile/halla/sbs/seeds/hcal_calibrations/pass0/timing/adctalign_class_gmn_conf4_qr0_pass0.root", "READ"); // Open the ROOT file in read mode
  
//   TCanvas *canvas = nullptr; // Declare a pointer to TCanvas

//   TIter next(file->GetListOfKeys());
//   TKey *key;
//   while ((key = static_cast<TKey *>(next()))) {
//     TObject *obj = key->ReadObj();
//     if (obj->IsA() == TCanvas::Class()) {
//       canvas = static_cast<TCanvas *>(obj);
//       // Do whatever you need with the canvas
//       canvas->Draw();
//     }
//     //delete obj; // Clean up the memory
//   }

//   //file->Close(); // Close the ROOT file
//   //delete file; // Clean up the memory
// }
