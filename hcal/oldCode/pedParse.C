#include <fstream>
#include <iostream>
#include <stdio.h>
#include <math.h>

void pedParse(){

  std::ifstream infile("Pedestals_989.txt");

  int n1;
  double d1, d2;
  
  int n_rows = 12;
  int n_cols = 12;
  
  double Pedestals[12][12] = {0};

  //int count;
  //infile>>count;
  infile.ignore(3, '\n');
  
  int rval, cval;
  
  while (infile >> n1 >> d1 >> d2)
    {
      rval = floor(n1/n_rows);
      cval = n1 % n_cols;
      
      Pedestals[rval][cval] = d1;

      std::cout << "row = " << rval << ", col = " << cval << ", array val = " << Pedestals[rval][cval] << std::endl;
    }
}
