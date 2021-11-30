#include <TH2F.h>
#include <TChain.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <TSystem.h>
#include <vector>
#include <TStopwatch.h>
#include <TMath.h>
#include <TF1.h>
#include<TStyle.h>
#include<TList.h>
#include<TChain.h>
#include <TDatime.h>
#include <cmath>
using namespace std;

const Int_t kNrows = 24;
const Int_t kNcols = 12;
const Int_t channels = kNrows * kNcols;
const Int_t hv_card_chs = 12;
const Int_t hv_cards = 12;
Int_t bad_CMU_offset = 0;    //Voltage added to bad CMU modules.
Int_t bad_JLab_offset = 0;  //Voltage added to bad JLab modules.
Int_t CMU_offset = 300;    //Voltage added to CMU modules.
Int_t JLab_offset = 1200;  //Voltage added to JLab modules.

std::string user_input;
Int_t skip = 1;                            //Gives number of lines to skip at top of data file. 
Int_t nlines = 0;                          //Counts number of lines in the data file. 
Int_t ncols;                               //Set how many columns of data we have in the data file.
char str[1000];                           //Variable to read lines of the data file.
Int_t mod[channels], hv[channels];
Float_t mod_temp, hv_temp;
//Float_t rpi20_S0[hv_card_chs],rpi20_S1[hv_card_chs],rpi20_S3[hv_card_chs],rpi20_S4[hv_card_chs],rpi20_S6[hv_card_chs],rpi20_S7[hv_card_chs],rpi20_S9[hv_card_chs],rpi20_S10[hv_card_chs],rpi20_S11[hv_card_chs],rpi20_S13[hv_card_chs],rpi20_S14[hv_card_chs],rpi20_S15[hv_card_chs];
//Float_t rpi21_S0[hv_card_chs],rpi21_S1[hv_card_chs],rpi21_S3[hv_card_chs],rpi21_S4[hv_card_chs],rpi21_S6[hv_card_chs],rpi21_S7[hv_card_chs],rpi21_S9[hv_card_chs],rpi21_S10[hv_card_chs],rpi21_S11[hv_card_chs],rpi21_S13[hv_card_chs],rpi21_S14[hv_card_chs],rpi21_S15[hv_card_chs];
//rpi21 S2 is cosmic scintillator paddles.
//Int_t rpi21_S2[hv_card_chs] = {0, -1200, 0, -1850, 0, -1000, 0, -1100, 0, 0, 0, 0};
Int_t hv_card_first_modules_rpi20[12] = {73,74,77,78,81,82,1,2,5,6,9,10};//Each entry is the module number of the first PMT HV in the first channel of that HV card. HV cards are in ascending order. i.e. hv_card_first_modules_rpi20[6] = 1 which means for the sixth HV card in rpi20 the first PMT HV in the first channel is PMT module 1.
Int_t hv_card_first_modules_rpi21[12] = {145,146,149,150,153,154,217,218,221,222,225,226};//Each entry is the module number of the first PMT HV in the first channel of that HV card. HV cards are in ascending order.i.e. hv_card_first_modules_rpi20[6] = 1 which means for the sixth HV card in rpi20 the first PMT HV in the first channel is PMT module 1. Note for rpi21 there is the cosmic paddle card which is filled separately. 
std::pair <Int_t,Int_t> hvs[channels];

void Generate_HV_GUI_Settings(Int_t run = 1271)
{
  //Define a new stopwatch.
  TStopwatch *st=new TStopwatch();
  st->Start(kTRUE);

  //Create a date object.
  TDatime time;

 FILE *fp;
 fp = fopen(Form("/home/daq/test_fadc/Voltage_Scans/HV_Settings/HV_Calibrated_from_Run_%d.txt",run),"r");
 //Read in data.
  while (1) 
    {
      //Skips the first skip lines of the file. 
      if (nlines < skip)
	{
	  fgets(str,1000,fp);
	  nlines++;
	}
      //Reads the two columns of data into x and y.
      else
	{
	  //Read in the number of columns of data in your data file. 
	  ncols = fscanf(fp,"%f %f", &mod_temp, &hv_temp);
	  
	  if (ncols < 0) break;    
	  
	  mod[nlines-skip] = (int) mod_temp;
	  hv[nlines-skip] = (int) hv_temp;
	  hv[nlines-skip] = hv[nlines-skip];

	  //Create a vector of pairs holding (PMT module #, that PMT's HV setting).
	  hvs[nlines-skip].first = (int) mod_temp;                  // the type of first is string
	  hvs[nlines-skip].second = (int) hv_temp; 
	  
	  nlines++;
	}
    }

  for(Int_t i=0; i<channels; i++)
    {
      //cout<<Form("PMT module %d has HV = %d.",hvs[i].first,hvs[i].second)<<endl;
    }

  Int_t rpi20_S0[hv_card_chs] = {hv[72]+bad_CMU_offset,hv[74]+bad_CMU_offset,hv[84],hv[86],hv[96],hv[98],hv[108]+bad_CMU_offset,hv[110],hv[120],hv[122],hv[132],hv[134]};
  Int_t rpi20_S1[hv_card_chs] = {hv[73]+bad_CMU_offset,hv[75]+bad_CMU_offset,hv[85],hv[87],hv[97],hv[99],hv[109],hv[111],hv[121],hv[123],hv[133],hv[135]};
  Int_t rpi20_S3[hv_card_chs] = {hv[76]+bad_JLab_offset,hv[78],hv[88],hv[90],hv[100],hv[102]+bad_JLab_offset,hv[112],hv[114],hv[124],hv[126],hv[136],hv[138]};
  Int_t rpi20_S4[hv_card_chs] = {hv[77]+bad_JLab_offset,hv[79],hv[89],hv[91],hv[101],hv[103],hv[113],hv[115],hv[125],hv[127],hv[137],hv[139]};
  Int_t rpi20_S6[hv_card_chs] = {hv[80],hv[82],hv[92],hv[94],hv[104],hv[106],hv[116],hv[118],hv[128],hv[130]+bad_CMU_offset,hv[140],hv[142]};
  Int_t rpi20_S7[hv_card_chs] = {hv[81],hv[83]+bad_CMU_offset,hv[93],hv[95],hv[105],hv[107],hv[117],hv[119],hv[129],hv[131]+bad_CMU_offset,hv[141],hv[143]};
  Int_t rpi20_S9[hv_card_chs] = {hv[0],hv[2],hv[12],hv[14],hv[24],hv[26]+bad_CMU_offset,hv[36],hv[38],hv[48],hv[50],hv[60],hv[62]};
  Int_t rpi20_S10[hv_card_chs] = {hv[1],hv[3],hv[13],hv[15],hv[25]+bad_CMU_offset,hv[27],hv[37],hv[39],hv[49],hv[51],hv[61],hv[63]};
  Int_t rpi20_S11[hv_card_chs] = {hv[4]+bad_JLab_offset,hv[6],hv[16],hv[18]+bad_JLab_offset,hv[28],hv[30]+bad_JLab_offset,hv[40],hv[42],hv[52],hv[54],hv[64],hv[66]};
  Int_t rpi20_S13[hv_card_chs] = {hv[5],hv[7],hv[17],hv[19]+bad_JLab_offset,hv[29],hv[31]+bad_JLab_offset,hv[41],hv[43],hv[53],hv[55],hv[65],hv[67]};
  Int_t rpi20_S14[hv_card_chs] = {hv[8],hv[10],hv[20]+bad_CMU_offset,hv[22]+bad_CMU_offset,hv[32]+bad_CMU_offset,hv[34]+bad_CMU_offset,hv[44],hv[46],hv[56],hv[58],hv[68],hv[70]};
  Int_t rpi20_S15[hv_card_chs] = {hv[9],hv[11],hv[21]+bad_CMU_offset,hv[23]+bad_CMU_offset,hv[33],hv[35],hv[45],hv[47],hv[57],hv[59],hv[69],hv[71]};

  Int_t rpi21_S0[hv_card_chs] = {hv[144]+bad_CMU_offset,hv[146]+bad_CMU_offset,hv[156],hv[158],hv[168]+bad_CMU_offset,hv[170],hv[180],hv[182],hv[192],hv[194],hv[204],hv[206]};
  Int_t rpi21_S1[hv_card_chs] = {hv[145]+bad_CMU_offset,hv[147]+bad_CMU_offset,hv[157]+bad_CMU_offset,hv[159],hv[169],hv[171],hv[181],hv[183],hv[193],hv[195],hv[205],hv[207]};
  Int_t rpi21_S2[hv_card_chs] = {0, -1110, 0, -1850, 0, -1280, 0, -1123, 0, 0, 0, 0};
  Int_t rpi21_S3[hv_card_chs] = {hv[148]+bad_JLab_offset,hv[150],hv[160],hv[162],hv[172],hv[174],hv[184],hv[186],hv[196],hv[198],hv[208],hv[210]};
  Int_t rpi21_S4[hv_card_chs] = {hv[149]+bad_JLab_offset,hv[151],hv[161],hv[163],hv[173],hv[175],hv[185],hv[187],hv[197],hv[199],hv[209],hv[211]};
  Int_t rpi21_S6[hv_card_chs] = {hv[152],hv[154],hv[164],hv[166],hv[176],hv[178],hv[188],hv[190],hv[200],hv[202],hv[212],hv[214]};
  Int_t rpi21_S7[hv_card_chs] = {hv[153]+bad_CMU_offset,hv[155],hv[165],hv[167],hv[177]+bad_CMU_offset,hv[179],hv[189],hv[191],hv[201],hv[203]+bad_CMU_offset,hv[213],hv[215]};
  Int_t rpi21_S9[hv_card_chs] = {hv[216],hv[218],hv[228],hv[230],hv[240],hv[242],hv[252],hv[254],hv[264],hv[266],hv[276],hv[278]};
  Int_t rpi21_S10[hv_card_chs] = {hv[217],hv[219],hv[229],hv[231],hv[241],hv[243],hv[253],hv[255],hv[265],hv[267],hv[277],hv[279]};
  Int_t rpi21_S11[hv_card_chs] = {hv[220],hv[222],hv[232],hv[234],hv[244],hv[246],hv[256],hv[258],hv[268],hv[270],hv[280],hv[282]};
  Int_t rpi21_S13[hv_card_chs] = {hv[221],hv[223],hv[233],hv[235],hv[245],hv[247],hv[257],hv[259],hv[269],hv[271],hv[281],hv[283]+bad_JLab_offset};
  Int_t rpi21_S14[hv_card_chs] = {hv[224],hv[226],hv[236],hv[238],hv[248],hv[250],hv[260],hv[262],hv[272],hv[274],hv[284],hv[286]};
  Int_t rpi21_S15[hv_card_chs] = {hv[225],hv[227],hv[237],hv[239],hv[249],hv[251],hv[261],hv[263],hv[273],hv[275],hv[285],hv[287]};
  
  for(Int_t i=0; i<hv_card_chs; i++) 
    {
      rpi20_S0[i] = rpi20_S0[i]+CMU_offset;
      rpi20_S1[i] = rpi20_S1[i]+CMU_offset;
      rpi20_S3[i] = rpi20_S3[i]+JLab_offset;
      rpi20_S4[i] = rpi20_S4[i]+JLab_offset;
      rpi20_S6[i] = rpi20_S6[i]+CMU_offset;
      rpi20_S7[i] = rpi20_S7[i]+CMU_offset;
      rpi20_S9[i] = rpi20_S9[i]+CMU_offset;
      rpi20_S10[i] = rpi20_S10[i]+CMU_offset;
      rpi20_S11[i] = rpi20_S11[i]+JLab_offset;
      rpi20_S13[i] = rpi20_S13[i]+JLab_offset;
      rpi20_S14[i] = rpi20_S14[i]+CMU_offset;
      rpi20_S15[i] = rpi20_S15[i]+CMU_offset;
      
      rpi21_S0[i] = rpi21_S0[i]+CMU_offset;
      rpi21_S1[i] = rpi21_S1[i]+CMU_offset;
      rpi21_S3[i] = rpi21_S3[i]+JLab_offset;
      rpi21_S4[i] = rpi21_S4[i]+JLab_offset;
      rpi21_S6[i] = rpi21_S6[i]+CMU_offset;
      rpi21_S7[i] = rpi21_S7[i]+CMU_offset;
      rpi21_S9[i] = rpi21_S9[i]+CMU_offset;
      rpi21_S10[i] = rpi21_S10[i]+CMU_offset;
      rpi21_S11[i] = rpi21_S11[i]+JLab_offset;
      rpi21_S13[i] = rpi21_S13[i]+JLab_offset;
      rpi21_S14[i] = rpi21_S14[i]+CMU_offset;
      rpi21_S15[i] = rpi21_S15[i]+CMU_offset;
    }
  

  //Int_t cmu_volt = -1800;
  //Int_t jlab_volt = -2200;
  /*
  for(Int_t i=0; i<hv_card_chs; i++) 
    {
      rpi20_S0[i] = cmu_volt;
      rpi20_S1[i] = cmu_volt;
      rpi20_S3[i] = jlab_volt;
      rpi20_S4[i] = jlab_volt;
      rpi20_S6[i] = cmu_volt;
      rpi20_S7[i] = cmu_volt;
      rpi20_S9[i] = cmu_volt;
      rpi20_S10[i] = cmu_volt;
      rpi20_S11[i] = jlab_volt;
      rpi20_S13[i] = jlab_volt;
      rpi20_S14[i] = cmu_volt;
      rpi20_S15[i] = cmu_volt;
      
      rpi21_S0[i] = cmu_volt;
      rpi21_S1[i] = cmu_volt;
      rpi21_S3[i] = jlab_volt;
      rpi21_S4[i] = jlab_volt;
      rpi21_S6[i] = cmu_volt;
      rpi21_S7[i] = cmu_volt;
      rpi21_S9[i] = cmu_volt;
      rpi21_S10[i] = cmu_volt;
      rpi21_S11[i] = jlab_volt;
      rpi21_S13[i] = jlab_volt;
      rpi21_S14[i] = cmu_volt;
      rpi21_S15[i] = cmu_volt;
    }
*/

  /*
  //Loop over all HV cards (ignoring cosmic paddle card) for the PMTs.
  for(Int_t i=0; i<hv_card; i++)
    {
      Int_t rpi20_start = hv_card_first_modules_rpi20[i];
      Int_t rpi21_start = hv_card_first_modules_rpi21[i];

      //Loop over and fill each hv card channel.
      for(Int_t j=0; j<hv_card; j++)
	{
	  if(j%2!=0)
	    {
	    }
	}
    }
  */

  /*
  for(Int_t i=0; i<hv_card_chs; i++)
    {
      //
      if((int) mod[i]%2!=0)
	{
	  rpi20_S6[i] = hvs[hv_card_first_modules_rpi20[6]+(i*12)-1].second;
	}
      if((int) mod[i]%2==0)
	{
	  rpi20_S6[i] = hvs[hv_card_first_modules_rpi20[6]+(i*12)+1].second;
	}
    }
  */

  cout<<"rpi20_S0 = ";
  for(Int_t i=0; i<hv_card_chs; i++)
    {
      cout<<rpi20_S0[i]<<", ";
    }

  //Open file to write output to.
  //std::ofstream output1 (Form("/home/daq/slowc/UPPER/hv_set/Voltages_Calibrated_from_Run_%d.hv",run), std::ofstream::out);
  //std::ofstream output2 (Form("/home/daq/slowc/LOWER/hv_set/Voltages_Calibrated_from_Run_%d.hv",run), std::ofstream::out);

  std::ofstream output1 (Form("/home/daq/slowc/UPPER/hv_set/Voltages_Calibrated_from_Run_%d_CMU_+%d_JLab_+%d.hv",run,CMU_offset,JLab_offset), std::ofstream::out);
  std::ofstream output2 (Form("/home/daq/slowc/LOWER/hv_set/Voltages_Calibrated_from_Run_%d_CMU_+%d_JLab_+%d.hv",run,CMU_offset,JLab_offset), std::ofstream::out);

  //std::ofstream output1 (Form("/home/daq/slowc/UPPER/hv_set/CMU_%d_JLab_%d.hv",cmu_volt,jlab_volt), std::ofstream::out);
  //std::ofstream output2 (Form("/home/daq/slowc/LOWER/hv_set/CMU_%d_JLab_%d.hv",cmu_volt,jlab_volt), std::ofstream::out);

  output1<<Form("# This file was generated at %d:%d:%d on %d/%d/%d. HVFrame-0(rpi20:2001) demand voltage set(DV)",time.GetHour(),time.GetMinute(),time.GetSecond(),time.GetMonth(),time.GetDay(),time.GetYear())<<endl;
  output1<<Form("rpi20:2001 S0 DV %d %d %d %d %d %d %d %d %d %d %d %d",rpi20_S0[0],rpi20_S0[1],rpi20_S0[2],rpi20_S0[3],rpi20_S0[4],rpi20_S0[5],rpi20_S0[6],rpi20_S0[7],rpi20_S0[8],rpi20_S0[9],rpi20_S0[10],rpi20_S0[11])<<endl;
  output1<<Form("rpi20:2001 S1 DV %d %d %d %d %d %d %d %d %d %d %d %d",rpi20_S1[0],rpi20_S1[1],rpi20_S1[2],rpi20_S1[3],rpi20_S1[4],rpi20_S1[5],rpi20_S1[6],rpi20_S1[7],rpi20_S1[8],rpi20_S1[9],rpi20_S1[10],rpi20_S1[11])<<endl;
  output1<<Form("rpi20:2001 S3 DV %d %d %d %d %d %d %d %d %d %d %d %d",rpi20_S3[0],rpi20_S3[1],rpi20_S3[2],rpi20_S3[3],rpi20_S3[4],rpi20_S3[5],rpi20_S3[6],rpi20_S3[7],rpi20_S3[8],rpi20_S3[9],rpi20_S3[10],rpi20_S3[11])<<endl;
  output1<<Form("rpi20:2001 S4 DV %d %d %d %d %d %d %d %d %d %d %d %d",rpi20_S4[0],rpi20_S4[1],rpi20_S4[2],rpi20_S4[3],rpi20_S4[4],rpi20_S4[5],rpi20_S4[6],rpi20_S4[7],rpi20_S4[8],rpi20_S4[9],rpi20_S4[10],rpi20_S4[11])<<endl;
  output1<<Form("rpi20:2001 S6 DV %d %d %d %d %d %d %d %d %d %d %d %d",rpi20_S6[0],rpi20_S6[1],rpi20_S6[2],rpi20_S6[3],rpi20_S6[4],rpi20_S6[5],rpi20_S6[6],rpi20_S6[7],rpi20_S6[8],rpi20_S6[9],rpi20_S6[10],rpi20_S6[11])<<endl;
  output1<<Form("rpi20:2001 S7 DV %d %d %d %d %d %d %d %d %d %d %d %d",rpi20_S7[0],rpi20_S7[1],rpi20_S7[2],rpi20_S7[3],rpi20_S7[4],rpi20_S7[5],rpi20_S7[6],rpi20_S7[7],rpi20_S7[8],rpi20_S7[9],rpi20_S7[10],rpi20_S7[11])<<endl;
  output1<<Form("rpi20:2001 S9 DV %d %d %d %d %d %d %d %d %d %d %d %d",rpi20_S9[0],rpi20_S9[1],rpi20_S9[2],rpi20_S9[3],rpi20_S9[4],rpi20_S9[5],rpi20_S9[6],rpi20_S9[7],rpi20_S9[8],rpi20_S9[9],rpi20_S9[10],rpi20_S9[11])<<endl;
  output1<<Form("rpi20:2001 S10 DV %d %d %d %d %d %d %d %d %d %d %d %d",rpi20_S10[0],rpi20_S10[1],rpi20_S10[2],rpi20_S10[3],rpi20_S10[4],rpi20_S10[5],rpi20_S10[6],rpi20_S10[7],rpi20_S10[8],rpi20_S10[9],rpi20_S10[10],rpi20_S10[11])<<endl;
  output1<<Form("rpi20:2001 S11 DV %d %d %d %d %d %d %d %d %d %d %d %d",rpi20_S11[0],rpi20_S11[1],rpi20_S11[2],rpi20_S11[3],rpi20_S11[4],rpi20_S11[5],rpi20_S11[6],rpi20_S11[7],rpi20_S11[8],rpi20_S11[9],rpi20_S11[10],rpi20_S11[11])<<endl;
  output1<<Form("rpi20:2001 S13 DV %d %d %d %d %d %d %d %d %d %d %d %d",rpi20_S13[0],rpi20_S13[1],rpi20_S13[2],rpi20_S13[3],rpi20_S13[4],rpi20_S13[5],rpi20_S13[6],rpi20_S13[7],rpi20_S13[8],rpi20_S13[9],rpi20_S13[10],rpi20_S13[11])<<endl;
  output1<<Form("rpi20:2001 S14 DV %d %d %d %d %d %d %d %d %d %d %d %d",rpi20_S14[0],rpi20_S14[1],rpi20_S14[2],rpi20_S14[3],rpi20_S14[4],rpi20_S14[5],rpi20_S14[6],rpi20_S14[7],rpi20_S14[8],rpi20_S14[9],rpi20_S14[10],rpi20_S14[11])<<endl;
  output1<<Form("rpi20:2001 S15 DV %d %d %d %d %d %d %d %d %d %d %d %d",rpi20_S15[0],rpi20_S15[1],rpi20_S15[2],rpi20_S15[3],rpi20_S15[4],rpi20_S15[5],rpi20_S15[6],rpi20_S15[7],rpi20_S15[8],rpi20_S15[9],rpi20_S15[10],rpi20_S15[11])<<endl;

  output2<<Form("# This file was generated at %d:%d:%d on %d/%d/%d. HVFrame-0(rpi21:2001) demand voltage set(DV)",time.GetHour(),time.GetMinute(),time.GetSecond(),time.GetMonth(),time.GetDay(),time.GetYear())<<endl;
  output2<<Form("rpi21:2001 S0 DV %d %d %d %d %d %d %d %d %d %d %d %d",rpi21_S0[0],rpi21_S0[1],rpi21_S0[2],rpi21_S0[3],rpi21_S0[4],rpi21_S0[5],rpi21_S0[6],rpi21_S0[7],rpi21_S0[8],rpi21_S0[9],rpi21_S0[10],rpi21_S0[11])<<endl;
  output2<<Form("rpi21:2001 S1 DV %d %d %d %d %d %d %d %d %d %d %d %d",rpi21_S1[0],rpi21_S1[1],rpi21_S1[2],rpi21_S1[3],rpi21_S1[4],rpi21_S1[5],rpi21_S1[6],rpi21_S1[7],rpi21_S1[8],rpi21_S1[9],rpi21_S1[10],rpi21_S1[11])<<endl;
  output2<<Form("rpi21:2001 S2 DV %d %d %d %d %d %d %d %d %d %d %d %d",rpi21_S2[0],rpi21_S2[1],rpi21_S2[2],rpi21_S2[3],rpi21_S2[4],rpi21_S2[5],rpi21_S2[6],rpi21_S2[7],rpi21_S2[8],rpi21_S2[9],rpi21_S2[10],rpi21_S2[11])<<endl;
  output2<<Form("rpi21:2001 S3 DV %d %d %d %d %d %d %d %d %d %d %d %d",rpi21_S3[0],rpi21_S3[1],rpi21_S3[2],rpi21_S3[3],rpi21_S3[4],rpi21_S3[5],rpi21_S3[6],rpi21_S3[7],rpi21_S3[8],rpi21_S3[9],rpi21_S3[10],rpi21_S3[11])<<endl;
  output2<<Form("rpi21:2001 S4 DV %d %d %d %d %d %d %d %d %d %d %d %d",rpi21_S4[0],rpi21_S4[1],rpi21_S4[2],rpi21_S4[3],rpi21_S4[4],rpi21_S4[5],rpi21_S4[6],rpi21_S4[7],rpi21_S4[8],rpi21_S4[9],rpi21_S4[10],rpi21_S4[11])<<endl;
  output2<<Form("rpi21:2001 S6 DV %d %d %d %d %d %d %d %d %d %d %d %d",rpi21_S6[0],rpi21_S6[1],rpi21_S6[2],rpi21_S6[3],rpi21_S6[4],rpi21_S6[5],rpi21_S6[6],rpi21_S6[7],rpi21_S6[8],rpi21_S6[9],rpi21_S6[10],rpi21_S6[11])<<endl;
  output2<<Form("rpi21:2001 S7 DV %d %d %d %d %d %d %d %d %d %d %d %d",rpi21_S7[0],rpi21_S7[1],rpi21_S7[2],rpi21_S7[3],rpi21_S7[4],rpi21_S7[5],rpi21_S7[6],rpi21_S7[7],rpi21_S7[8],rpi21_S7[9],rpi21_S7[10],rpi21_S7[11])<<endl;
  output2<<Form("rpi21:2001 S9 DV %d %d %d %d %d %d %d %d %d %d %d %d",rpi21_S9[0],rpi21_S9[1],rpi21_S9[2],rpi21_S9[3],rpi21_S9[4],rpi21_S9[5],rpi21_S9[6],rpi21_S9[7],rpi21_S9[8],rpi21_S9[9],rpi21_S9[10],rpi21_S9[11])<<endl;
  output2<<Form("rpi21:2001 S10 DV %d %d %d %d %d %d %d %d %d %d %d %d",rpi21_S10[0],rpi21_S10[1],rpi21_S10[2],rpi21_S10[3],rpi21_S10[4],rpi21_S10[5],rpi21_S10[6],rpi21_S10[7],rpi21_S10[8],rpi21_S10[9],rpi21_S10[10],rpi21_S10[11])<<endl;
  output2<<Form("rpi21:2001 S11 DV %d %d %d %d %d %d %d %d %d %d %d %d",rpi21_S11[0],rpi21_S11[1],rpi21_S11[2],rpi21_S11[3],rpi21_S11[4],rpi21_S11[5],rpi21_S11[6],rpi21_S11[7],rpi21_S11[8],rpi21_S11[9],rpi21_S11[10],rpi21_S11[11])<<endl;
  output2<<Form("rpi21:2001 S13 DV %d %d %d %d %d %d %d %d %d %d %d %d",rpi21_S13[0],rpi21_S13[1],rpi21_S13[2],rpi21_S13[3],rpi21_S13[4],rpi21_S13[5],rpi21_S13[6],rpi21_S13[7],rpi21_S13[8],rpi21_S13[9],rpi21_S13[10],rpi21_S13[11])<<endl;
  output2<<Form("rpi21:2001 S14 DV %d %d %d %d %d %d %d %d %d %d %d %d",rpi21_S14[0],rpi21_S14[1],rpi21_S14[2],rpi21_S14[3],rpi21_S14[4],rpi21_S14[5],rpi21_S14[6],rpi21_S14[7],rpi21_S14[8],rpi21_S14[9],rpi21_S14[10],rpi21_S14[11])<<endl;
  output2<<Form("rpi21:2001 S15 DV %d %d %d %d %d %d %d %d %d %d %d %d",rpi21_S15[0],rpi21_S15[1],rpi21_S15[2],rpi21_S15[3],rpi21_S15[4],rpi21_S15[5],rpi21_S15[6],rpi21_S15[7],rpi21_S15[8],rpi21_S15[9],rpi21_S15[10],rpi21_S15[11])<<endl;

  //Stop and print out stopwatch information.
  st->Stop();
  cout<<"CPU time = "<<st->CpuTime()<<" s = "<<st->CpuTime()/60.<<" min   Real time = "<<st->RealTime()<<" s = "<<st->RealTime()/60.<<" min"<<endl;
}
