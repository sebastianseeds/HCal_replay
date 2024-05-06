//converts channel map output from add_crate.C to latex table friendly rows
#include <fstream>
#include <iostream>
#include <sstream>

int map_to_latex() {
    // Input and output file names
    std::string inputFileName = "textfiles/channel_map_updated.txt";
    std::string outputFileName = "textfiles/latex_rows.txt";

    // Open the input file for reading
    std::ifstream inputFile(inputFileName);
    if (!inputFile.is_open()) {
        std::cerr << "Could not open the input file." << std::endl;
        return 1;
    }

    // Open the output file for writing
    std::ofstream outputFile(outputFileName);
    if (!outputFile.is_open()) {
        std::cerr << "Could not open the output file." << std::endl;
        return 1;
    }

    std::string line;
    int lineNumber = 0;
    // Read each line from the input file
    while (std::getline(inputFile, line)) {
        // Skip lines that start with '#'
        if (!line.empty() && line[0] == '#') {
            continue;
        }

        std::istringstream iss(line);
        std::string item;
        bool first = true;

        // Add \rowcolor{lightgray} for every other row starting from the second row
        if (lineNumber % 2 == 1) {
            outputFile << "\\rowcolor{lightgray} ";
        }
        lineNumber++;

        // Split the line by spaces and write to the output file
        while (iss >> item) {
            if (!first) {
                outputFile << " & ";
            }
            outputFile << item;
            first = false;
        }

        // Finish the LaTeX table row
        outputFile << " \\\\" << std::endl;
    }

    // Close the files
    inputFile.close();
    outputFile.close();

    std::cout << "Processing completed. Latex rows located " << outputFileName << std::endl;
    return 0;
}
