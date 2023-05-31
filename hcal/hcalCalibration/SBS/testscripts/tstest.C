#include "../include/sbs.h"
#include <vector>
#include <iostream>
#include <algorithm>

void tstest(){

  std::string tsA = "-------[ 2021-10-24 04:30:00 ]";
  std::string tsB = "-------[ 2020-10-24 04:30:00 ]";
  std::string tsC = "-------[ 2021-05-24 04:30:00 ]";
  std::string tsD = "-------[ 2021-05-12 04:30:00 ]";
  std::string tsE = "-------[ 2021-05-12 02:30:00 ]";
  std::string tsF = "-------[ 2021-05-12 02:30:00  ]";
  std::string tsG = "-----------[ 2021-05-12 02:30:00 ]";
  std::string tsH = "-------[ 2021-04-12 02:30:00 ]";
  std::string tsI = "none";
  
  std::string testA = tsC;
  std::string testB = tsB;
  std::string testC = tsA;

  std::string temp;

  util::tsCompare(testA,testB,temp);
  util::tsCompare(temp,testC,temp);

  cout << temp << endl;

  // if( testA.compare(testB)!=0 ){
    
  //   if( testA.size() != testB.size() ){
  //     cout << "ERROR: string sizes don't match" << endl;
  //     return;
  //   }

  //   for( Int_t l=0; l<testA.size(); l++ ){
	
  // 	char a = testA[l];
  // 	char b = testB[l];
  // 	bool decision = false;

  // 	if( a!=b ){

  // 	  if( (Int_t)a>(Int_t)b ){
  // 	    //cout << testA << " > " << testB << endl;
  // 	    temp = testA;
  // 	    decision = true;
  // 	    //return;
  // 	  }else if( (Int_t)a<(Int_t)b ){
  // 	    //cout << testA << " < " << testB << endl;
  // 	    temp = testB;
  // 	    decision = true;
  // 	    //return;
  // 	  }else{
  // 	    //cout << testA << " = " << testB << endl;
  // 	    //cout << "ERROR: strings are the same" << endl;
  // 	    //return;
  //        continue;
  // 	  }
	
  // 	  if( decision ){
  // 	    cout << temp << endl;
  // 	    return;

  // 	  }
	  
  // 	}
      

  //     }
  //   }

}
