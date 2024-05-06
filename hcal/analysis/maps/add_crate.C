//script to add crate to channel map
//#Ch Amp ADCpp ADCpp ADC sp disc TDCpp TDCpp disc TDC sum HV PMT_r-c adc_crate_b-a tdc_crate

#include <fstream>
#include <iostream>
#include <sstream>

int add_crate() {
  std::ifstream inputFile("textfiles/channel_map.txt");
  std::ofstream outputFile("textfiles/channel_map_updated.txt");

  std::string line;
  while (std::getline(inputFile, line)) {
    std::istringstream iss(line);
    std::string token;
    int fieldCounter = 0;
    int fValue = 0;
    while (iss >> token) {
      ++fieldCounter;
      // Check the 5th field for the "f<x>" value
      if (fieldCounter == 5) {
	// Assuming the format is always "f<x>-<y>"
	size_t dashPos = token.find('-');
	if (dashPos != std::string::npos) {
	  std::string fNumStr = token.substr(1, dashPos - 1);
	  fValue = std::stoi(fNumStr);
	}
	break;
      }
    }

    int X = (fValue < 17) ? 16 : 17;
    int Y = (fValue < 17) ? 28 : 29;
    int Z = 17;

    outputFile << line << " " << X << "-" << Y << " " << Z << std::endl;
  }

  inputFile.close();
  outputFile.close();

  return 0;
}
