//sseeds 4.24.24 script to aggregate led plateau and alpha parameter results into landscape latex table
//alphas exist for all cells. Due to fiber and led errors in the test lab, data exist only for the bottom half of hcal for plateau regions. Averages for the CMU (col 1-4 + 9-12) and JLab (col 5-8) type PMTs included here for each LED setting on top half.
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <iomanip>
#include <cmath>  // For round function

using namespace std;

struct Plateau {
  double start;
  double end;
};

map<int, Plateau> readLEDFile(const string& filename) {
  ifstream file(filename);
  map<int, Plateau> data;
  string line;
  int i = 145;  // Start from cell 145
  while (getline(file, line)) {
    stringstream ss(line);
    double start, end;
    char comma;
    ss >> start >> comma >> end;
    data[i++] = {start, end};
  }
  file.close();
  return data;
}

Plateau computeAverage(const map<int, Plateau>& data, bool center_region) {
  double sumStart = 0, sumEnd = 0;
  int count = 0;

  for (int i = 145; i <= 288; ++i) {
    int col = (i-1)%12 + 1;
    bool CMU = col<5 || col>8;
    bool JLab = col>4 && col<9;

    if( (center_region && CMU) || (!center_region && JLab) )
      continue;

    sumStart += data.at(i).start;
    sumEnd += data.at(i).end;
    count++;
  }
  // Round to the nearest whole number
  return {round(sumStart / count), round(sumEnd / count)};
}

void writeLatexTable(ofstream& out, int cell, double alpha, const vector<Plateau>& plateaus) {
  if( (cell-1) % 24 == 0 )
    out << "\\rowcolor{lightgray} ";

  out << fixed << setprecision(2); // Set the number format to fixed point with two decimals
  out << "cell: " << cell << ", alpha: " << alpha << ", plateaus: ";
  int iter = 1;
  for (auto& plateau : plateaus) {
    out << "led" << iter << " " << (int)plateau.start << "-" << (int)plateau.end;
    if(iter<4)
      out << ", ";
    iter++;
  }

  if( cell % 12 == 0 )
    out << " \\\\\n";
  else
    out << " &";

}

int agg_alphas_plat() {
  vector<double> alphas;
  ifstream alphaFile("alphas.txt");
  string line;
  while (getline(alphaFile, line)) {
    alphas.push_back(stod(line));
  }
  alphaFile.close();

  map<int, Plateau> led1 = readLEDFile("Vplat_l1.txt");
  map<int, Plateau> led2 = readLEDFile("Vplat_l2.txt");
  map<int, Plateau> led3 = readLEDFile("Vplat_l3.txt");
  map<int, Plateau> led4 = readLEDFile("Vplat_l4.txt");

  Plateau avgCMU_l1 = computeAverage(led1, false);
  Plateau avgCMU_l2 = computeAverage(led2, false);
  Plateau avgCMU_l3 = computeAverage(led3, false);
  Plateau avgCMU_l4 = computeAverage(led4, false);
  Plateau avgJLab_l1 = computeAverage(led1, true);
  Plateau avgJLab_l2 = computeAverage(led2, true);
  Plateau avgJLab_l3 = computeAverage(led3, true);
  Plateau avgJLab_l4 = computeAverage(led4, true);

  ofstream outFile("agg_alphas_plat_table.tex");

  outFile << "\\documentclass{article}\n"
	  << "\\usepackage{lscape}\n"
	  << "\\usepackage{longtable}\n"
	  << "\\usepackage{fancyhdr} % For custom headers/footers\n"
	  << "\\usepackage{geometry} % For adjusting margins\n"
	  << "\\usepackage[table]{xcolor} % For table row colors\n"
	  << "\\usepackage{pdflscape}  % Provides the landscape environment for PDF-friendly rotation\n"
	  << "\n"
	  << "% Adjust the page layout to remove page numbers and potentially change margins\n"
	  << "\\geometry{\n"
	  << "  left=1in,\n"
	  << "  right=1in,\n"
	  << "  top=1in,\n"
	  << "  bottom=1in,\n"
	  << "  headheight=15pt,\n"
	  << "  includeheadfoot\n"
	  << "}\n"
	  << "\n"
	  << "% Define a page style for landscape pages without page numbers\n"
	  << "\\fancypagestyle{lscape}{\n"
	  << "  \\fancyhf{} % clear all header and footer fields\n"
	  << "  \\renewcommand{\\headrulewidth}{0pt} % remove the header rule\n"
	  << "  \\renewcommand{\\footrulewidth}{0pt} % remove the footer rule\n"
	  << "}\n"
	  << "\n"
	  << "\\pagestyle{plain} % Use the plain page style for normal pages\n"
	  << "\n"
	  << "\\begin{document}\n"
	  << "\\definecolor{darkblue}{rgb}{0.0, 0.0, 0.55}\n"
	  << "\\definecolor{capri}{rgb}{0.0, 0.75, 1.0}\n"
	  << "\\renewcommand{\\thesubsection}{\\thesection.\\alph{subsection}}\n"
	  << "\\newcommand{\\scripty}[1]{\\ensuremath{\\mathcalligra{#1}}\\;}\n"
	  << "\\newcommand*\\VF[1]{\\mathbf{#1}}\n"
	  << "\\newcommand*\\dif{\\mathop{}\\!\\mathrm{d}}\n"
	  << "\\newcommand{\\Lagr}{\\mathcal{L}}\n"
	  << "\\newcommand{\\Hami}{\\mathcal{H}}\n"
	  << "\\newcommand{\\elec}{4\\pi\\epsilon_0}\n"
	  << "\n"
	  << "\\begin{landscape}\n"
	  << "\\scriptsize\n"
	  << "\\begin{longtable}{|p{1.5cm}|p{1.5cm}|p{1.5cm}|p{1.5cm}|p{1.5cm}|p{1.5cm}|p{1.5cm}|p{1.5cm}|p{1.5cm}|p{1.5cm}|p{1.5cm}|p{1.5cm}|}\n"
	  << "\\caption{HCal LED Plateaus and Alphas}\\label{tab:hcalalphas} \\\\\n"
	  << "\\hline\n"
	  << "\\endfirsthead\n"
	  << "\n"
	  << "\\multicolumn{12}{c}%\n"
	  << "{{\\bfseries \\tablename\\ \\thetable{} -- continued from previous page}} \\\\\n"
	  << "\\hline\n"
	  << "\\endhead\n"
	  << "\n"
	  << "\\hline\n"
	  << "\\multicolumn{12}{|r|}{{Continued on next page}} \\\\ \\hline\n"
	  << "\\endfoot\n"
	  << "\n"
	  << "\\hline \\hline\n"
	  << "\\endlastfoot\n";

  for (int i = 1; i <= 288; ++i) {
    double alpha = alphas[i-1];
    vector<Plateau> plateaus;

    int col = (i-1)%12 + 1;
    bool CMU = col<5 || col>8;
    bool JLab = col>4 && col<9;

    if (i <= 144) {
      if (CMU) { 
	plateaus = {avgCMU_l1, avgCMU_l2, avgCMU_l3, avgCMU_l4};
      } else {
	plateaus = {avgJLab_l1, avgJLab_l2, avgJLab_l3, avgJLab_l4};
      }
    } else {
      plateaus = {led1[i], led2[i], led3[i], led4[i]};  // Corrected to access directly without adding 145
    }

    // Print details to console
    cout << "For cell " << i << ", alpha " << alpha 
	 << ", led1 range " << plateaus[0].start << "-" << plateaus[0].end
	 << ", led2 range " << plateaus[1].start << "-" << plateaus[1].end
	 << ", led3 range " << plateaus[2].start << "-" << plateaus[2].end
	 << ", led4 range " << plateaus[3].start << "-" << plateaus[3].end << endl;

    writeLatexTable(outFile, i, alpha, plateaus);
  }

  outFile << "\\end{longtable}\n"
	  << "\\end{landscape}\n"
	  << "\\end{document}";

  outFile.close();

  return 0;
}
