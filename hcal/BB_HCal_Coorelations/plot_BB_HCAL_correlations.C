#include "TChain.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include <iostream>
#include <fstream>

const double Mp = 0.938272;
const double Mn = 0.939565;
const int ncells_hcal = 288;

double expoFit(double *x, double *par) {
  double offset = par[0];
  double amp = par[1];
  double coeff = par[2];
  double e = x[0];
  return offset+TMath::Exp(amp+coeff*e);
}
/*
void plot_BB_HCAL_correlations(const char *rootfilename, 
			       const char *outfilename="outfiles/BBHCALCOOR_SBS4.root", 
			       double ebeam=3.7278, 
			       double bbtheta=36.0, 
			       double sbstheta=31.9, 
			       double hcaldist=11.0, 
			       double dx0=0.7, 
			       double dy0=0.12, 
			       double dxsigma=0.06, 
			       double dysigma=0.1, 
			       double Wmin=0.65, 
			       double Wmax=0.85, 
			       double dpel_min=-0.06, 
			       double dpel_max=0.06 ){
*/
/*
void plot_BB_HCAL_correlations(const char *rootfilename, 
			       const char *outfilename="outfiles/BBHCALCOOR_SBS14.root", 
			       double ebeam=5.965, 
			       double bbtheta=46.5, 
			       double sbstheta=17.5, 
			       double hcaldist=14.0, 
			       double dx0=0.0, 
			       double dy0=0.0, 
			       double dxsigma=0.06, 
			       double dysigma=0.1, 
			       double Wmin=0.6, 
			       double Wmax=1.2, 
			       double dpel_min=-0.06, 
			       double dpel_max=0.06 ){
  */


void plot_BB_HCAL_correlations(const char *rootfilename, 
			       const char *outfilename="outfiles/BBHCALCOOR_SBS8_LD2.root", 
			       double ebeam=5.965, 
			       double bbtheta=26.5, 
			       double sbstheta=29.4, 
			       double hcaldist=11., 
			       double dx0=0.9, 
			       double dy0=0.62, 
			       double dxsigma=0.09, 
			       double dysigma=0.15, 
			       double Wmin=1.00, 
			       double Wmax=1.3, 
			       double dpel_min=-0.06, 
			       double dpel_max=0.06 ){


  TChain *C = new TChain("T");

  //SBS-14 LH2 CH
  /*
  C->Add("/adaqfs/home/a-onl/sbs/Rootfiles/hcal_general_13241_-1*");
  C->Add("/adaqfs/home/a-onl/sbs/Rootfiles/hcal_general_13242_-1*");
  C->Add("/adaqfs/home/a-onl/sbs/Rootfiles/hcal_general_13243_-1*");
  C->Add("/adaqfs/home/a-onl/sbs/Rootfiles/hcal_general_13312_-1*");
  C->Add("/adaqfs/home/a-onl/sbs/Rootfiles/hcal_general_13313_-1*");
  C->Add("/adaqfs/home/a-onl/sbs/Rootfiles/hcal_general_13320_-1*");
  C->Add("/adaqfs/home/a-onl/sbs/Rootfiles/hcal_general_13321_-1*");
  C->Add("/adaqfs/home/a-onl/sbs/Rootfiles/e1209019_fullreplay_13349*");
  C->Add("/adaqfs/home/a-onl/sbs/Rootfiles/e1209019_fullreplay_13396*");
  */

  /*
  //SBS-4 LH2 CH
  C->Add("/adaqfs/home/a-onl/sbs/Rootfiles/hcal_gmn_11588*");
  C->Add("/adaqfs/home/a-onl/sbs/Rootfiles/hcal_gmn_11587*");
  */
  /*
  //SBS-4 LH2 ifarm
  C->Add("/volatile/halla/sbs/seeds/GMN_Replays/e1209019_fullreplay_11547*");
  C->Add("/volatile/halla/sbs/seeds/GMN_Replays/e1209019_fullreplay_11548*");
  //C->Add("/volatile/halla/sbs/seeds/GMN_Replays/e1209019_fullreplay_11587*");
  //C->Add("/volatile/halla/sbs/seeds/GMN_Replays/e1209019_fullreplay_11588*");
  //C->Add("/volatile/halla/sbs/seeds/GMN_Replays/e1209019_fullreplay_11573*");
  //C->Add("/volatile/halla/sbs/seeds/GMN_Replays/e1209019_fullreplay_11589*");
  //C->Add("/volatile/halla/sbs/seeds/GMN_Replays/e1209019_fullreplay_11590*");
  //C->Add("/volatile/halla/sbs/seeds/GMN_Replays/e1209019_fullreplay_11592*");
  C->Add("/volatile/halla/sbs/seeds/GMN_Replays/e1209019_fullreplay_11616*");
  */
  /*
  //SBS-8 LH2 ifarm
  C->Add("/volatile/halla/sbs/puckett/GMN_REPLAYS/rootfiles/e1209019_fullreplay_13483*");
  C->Add("/volatile/halla/sbs/puckett/GMN_REPLAYS/rootfiles/e1209019_fullreplay_13484*");
  C->Add("/volatile/halla/sbs/puckett/GMN_REPLAYS/rootfiles/e1209019_fullreplay_13485*");
  C->Add("/volatile/halla/sbs/puckett/GMN_REPLAYS/rootfiles/e1209019_fullreplay_13486*");
  C->Add("/volatile/halla/sbs/puckett/GMN_REPLAYS/rootfiles/e1209019_fullreplay_13487*");
  C->Add("/volatile/halla/sbs/puckett/GMN_REPLAYS/rootfiles/e1209019_fullreplay_13488*");
  C->Add("/volatile/halla/sbs/puckett/GMN_REPLAYS/rootfiles/e1209019_fullreplay_13489*");
  C->Add("/volatile/halla/sbs/puckett/GMN_REPLAYS/rootfiles/e1209019_fullreplay_13490*");
  */

  //SBS-8 LD2 ifarm
  C->Add("/volatile/halla/sbs/seeds/GMN_Replays/e1209019_fullreplay_13581*");
  C->Add("/volatile/halla/sbs/seeds/GMN_Replays/e1209019_fullreplay_13582*");
  C->Add("/volatile/halla/sbs/seeds/GMN_Replays/e1209019_fullreplay_13583*");
  C->Add("/volatile/halla/sbs/seeds/GMN_Replays/e1209019_fullreplay_13584*");
  C->Add("/volatile/halla/sbs/seeds/GMN_Replays/e1209019_fullreplay_13585*");
  C->Add("/volatile/halla/sbs/seeds/GMN_Replays/e1209019_fullreplay_13586*");
  C->Add("/volatile/halla/sbs/seeds/GMN_Replays/e1209019_fullreplay_13587*");
  C->Add("/volatile/halla/sbs/seeds/GMN_Replays/e1209019_fullreplay_13588*");
  C->Add("/volatile/halla/sbs/seeds/GMN_Replays/e1209019_fullreplay_13589*");
  C->Add("/volatile/halla/sbs/seeds/GMN_Replays/e1209019_fullreplay_13590*");
  C->Add("/volatile/halla/sbs/seeds/GMN_Replays/e1209019_fullreplay_13591*");
  C->Add("/volatile/halla/sbs/seeds/GMN_Replays/e1209019_fullreplay_13592*");
  C->Add("/volatile/halla/sbs/seeds/GMN_Replays/e1209019_fullreplay_13593*");

  //double uA_hr_SBS8_tot[8] = {7.13,3.08,2.92,6.32,7.44,6.78,5.89,6.90};
  double uA_hr_SBS8_tot[13] = {2.98,5.2,5.25,5.25,5.5,2.06,5.35,5.2,5.4,0.64,5.45,0.16,1.6};

  int totEvents = C->GetEntries();

  cout << "Total number of events to analyze: " << totEvents << "." << endl;

  bbtheta *= TMath::DegToRad();
  sbstheta *= TMath::DegToRad();

  double hcalheight = 0.365; //m The height of the center of HCAL above beam

  //The following are the positions of the "first" row and column from HCAL database (top right block as viewed from upstream)
  double xoff_hcal = 0.92835;
  double yoff_hcal = 0.47305;
  
  double blockspace_hcal = 0.15254;

  int nrows_hcal=24;
  int ncols_hcal=12;

  //int ncells_hcal=nrows_hcal*ncols_hcal;

  //For monitoring
  int posCutCount = 0;
  int BBCutCount = 0;
  int WCutCount = 0;
  int tracks = 0;
  int inTime = 0;
  
  double Wmin_elastic = Wmin;
  double Wmax_elastic = Wmax;
  
  //HCAL variables:
  double xHCAL, yHCAL, EHCAL, tHCAL;
  double tdcID[288];
  double aHCAL[288];
  double tHCAL_blk[288];

  //BigBite track variables:
  double ntrack;
  int MAXNTRACKS=1000;

  double px[MAXNTRACKS], py[MAXNTRACKS], pz[MAXNTRACKS], p[MAXNTRACKS];  
  double vx[MAXNTRACKS], vy[MAXNTRACKS], vz[MAXNTRACKS];

  //BigBite shower/preshower variables:
  double Eps_BB, Esh_BB;
  double xps_BB, yps_BB, xsh_BB, ysh_BB;
  Double_t bb_tdctrig_tdc[6] = {0.}, bb_tdctrig_tdcelemID[6] = {0.};
  Int_t Ndata_bb_tdctrig_tdcelemID = 0;

  //BigBite hodoscope variables:
  double tHODO;
  double nClusHODO;

  //Event variables
  UInt_t trig;
  
  C->SetBranchStatus("*",0);
  C->SetBranchStatus("sbs.hcal.x",1);
  C->SetBranchStatus("sbs.hcal.y",1);
  C->SetBranchStatus("sbs.hcal.e",1);
  C->SetBranchStatus("sbs.hcal.clus_blk.tdctime",1);
  C->SetBranchStatus("sbs.hcal.clus_blk.id",1);
  C->SetBranchStatus("sbs.hcal.clus_blk.e",1);
  //C->SetBranchStatus("sbs.hcal.tdctimeblk",1);
  C->SetBranchStatus("bb.tr.n",1);
  C->SetBranchStatus("bb.tr.px",1);
  C->SetBranchStatus("bb.tr.py",1);
  C->SetBranchStatus("bb.tr.pz",1);
  C->SetBranchStatus("bb.tr.vx",1);
  C->SetBranchStatus("bb.tr.vy",1);
  C->SetBranchStatus("bb.tr.vz",1);
  C->SetBranchStatus("bb.ps.e",1);
  C->SetBranchStatus("bb.ps.x",1);
  C->SetBranchStatus("bb.ps.y",1);
  C->SetBranchStatus("bb.sh.e",1);
  C->SetBranchStatus("bb.sh.x",1);
  C->SetBranchStatus("bb.sh.y",1);
  C->SetBranchStatus("bb.hodotdc.clus.tmean",1);
  C->SetBranchStatus("bb.hodotdc.nclus",1);
  C->SetBranchStatus("fEvtHdr.fTrigBits",1);
  C->SetBranchStatus("bb.tdctrig.tdc",1);
  C->SetBranchStatus("Ndata.bb.tdctrig.tdcelemID",1);
  C->SetBranchStatus("bb.tdctrig.tdcelemID",1);

  C->SetBranchAddress("sbs.hcal.x",&xHCAL);
  C->SetBranchAddress("sbs.hcal.y",&yHCAL);
  C->SetBranchAddress("sbs.hcal.e",&EHCAL);
  C->SetBranchAddress("sbs.hcal.clus_blk.tdctime",&tHCAL_blk);
  C->SetBranchAddress("sbs.hcal.clus_blk.id",&tdcID);
  C->SetBranchAddress("sbs.hcal.clus_blk.e",&aHCAL);
  //C->SetBranchAddress("sbs.hcal.tdctimeblk",&tHCAL);

  C->SetBranchAddress("bb.tr.n",&ntrack);
  C->SetBranchAddress("bb.tr.px",px);
  C->SetBranchAddress("bb.tr.py",py);
  C->SetBranchAddress("bb.tr.pz",pz);
  C->SetBranchAddress("bb.tr.p",p);
  C->SetBranchAddress("bb.tr.vx",vx);
  C->SetBranchAddress("bb.tr.vy",vy);
  C->SetBranchAddress("bb.tr.vz",vz);
  C->SetBranchAddress("bb.ps.e",&Eps_BB);
  C->SetBranchAddress("bb.sh.e",&Esh_BB);
  C->SetBranchAddress("bb.ps.x",&xps_BB);
  C->SetBranchAddress("bb.ps.y",&yps_BB);
  C->SetBranchAddress("bb.sh.x",&xsh_BB);
  C->SetBranchAddress("bb.sh.y",&ysh_BB);
  C->SetBranchAddress("bb.hodotdc.clus.tmean",&tHODO);
  C->SetBranchAddress("bb.hodotdc.nclus",&nClusHODO);
  C->SetBranchAddress("fEvtHdr.fTrigBits",&trig);
  C->SetBranchAddress("bb.tdctrig.tdc", &bb_tdctrig_tdc);
  C->SetBranchAddress("bb.tdctrig.tdcelemID", &bb_tdctrig_tdcelemID);
  C->SetBranchAddress("Ndata.bb.tdctrig.tdcelemID", &Ndata_bb_tdctrig_tdcelemID);

  // Output file
  TFile *fout = new TFile(outfilename,"RECREATE");

  // Declare histograms
  TH1D *hdpel = new TH1D("hdpel",";p/p_{elastic}(#theta)-1;", 250, -1.0, 0.5);
  TH1D *hW = new TH1D("hW",";W (GeV);", 400,0.0,4.0);

  TH1D *hdpel_cutBBCAL = new TH1D("hdpel_cutBBCAL","",250, -1.0,0.5);
  TH1D *hW_cutBBCAL = new TH1D("hW_cutBBCAL","",400,0.0,4.0);
  TH1D *htDiff_HODO_HCAL = new TH1D("htDiff_HODO_HCAL","",150,-150,0);
  
  TH1D *hdx_HCAL = new TH1D("hdx_HCAL",";x_{HCAL}-x_{expect} (m);", 500, -2.5, 2.5);
  TH1D *hdy_HCAL = new TH1D("hdy_HCAL",";y_{HCAL}-y_{expect} (m);", 500, -1.25, 1.25);
  TH2D *hdxdy_HCAL = new TH2D("hdxdy_HCAL",";y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)", 250, -1.25, 1.25, 250, -2.5, 2.5 );
  TH2D *hxcorr_HCAL = new TH2D("hxcorr_HCAL",";x_{expect} (m);x_{HCAL} (m)", 250, -2.5, 2.5, 250, -2.5, 2.5 );
  TH2D *hycorr_HCAL = new TH2D("hycorr_HCAL",";y_{expect} (m);y_{HCAL} (m)", 250, -1.25, 1.25, 250, -1.25, 1.25);
  TH2D *htDiff_vs_HCALID = new TH2D("htDiff_vs_HCALID",";Channel;TDC_{HCAL}-TDC_{HODO} (ns)",288,0,288,300,-150,0);
  TH2D *htBBHDiff_vs_HCALID = new TH2D("htBBHDiff_vs_HCALID",";Channel;TDC_{HCAL}-TDC_{BBTrig} (ns)",288,0,288,300,350,500);
  TH2D *htDiff_vs_HCALID_corr = new TH2D("htDiff_vs_HCALID_corr",";Channel;TDC_{HCAL}-TDC_{HODO} (ns)",288,0,288,300,-150,0);
  TH2D *htDiff_vs_HCALID_corrByChannel = new TH2D("htDiff_vs_HCALID_corrByChannel",";Channel;TDC_{HCAL}-TDC_{HODO} (ns)",288,0,288,300,-150,0);  

  TH2D *htDiff_vs_ADCint[ncells_hcal+1];

  for( int i=0; i<ncells_hcal+1; i++ ){
    htDiff_vs_ADCint[i] = new TH2D(Form("htDiff_vs_ADCint_bl%d",i),Form(";E_{bl%d} (GeV);TDC_{HCAL}-TDC_{HODO} (ns)",i),270,0.0,0.9,600,-250,50);
  }

  //TH2D *htDiff_vs_ADCint_bl175 = new TH2D("htDiff_vs_ADCint_bl175",";E_{bl175} (GeV);TDC_{HCAL}-TDC_{HODO} (ns)",90,0.0,0.3,300,-150,0); //
  //TH2D *htDiff_vs_ADCint_bl224 = new TH2D("htDiff_vs_ADCint_bl224",";E_{bl224} (GeV);TDC_{HCAL}-TDC_{HODO} (ns)",90,0.0,0.3,300,-150,0); //
  //TH2D *htDiff_vs_ADCint_bl92 = new TH2D("htDiff_vs_ADCint_bl92",";E_{bl92} (GeV);TDC_{HCAL}-TDC_{HODO} (ns)",90,0.0,0.3,300,-150,0); //

  TH2D *htcorr_HCAL_HODO = new TH2D("htcorr_HCAL_HODO",";TDC_{HCAL} (ns);TDC_{HODO} (ns)", 150, -150, 0, 150, -15, 15);
  TH2D *htcorr_HCAL_BBTrig = new TH2D("htcorr_HCAL_BBTrig",";TDC_{HCAL} (ns);TDC_{BBTrig} (ns)", 150, -150, 0, 100, 300, 400);

  TH1D *hvz = new TH1D("hvz","",250,-0.15,0.15);

  TH2D *hdy_HCAL_vs_z = new TH2D("hdy_HCAL_vs_z","",250,-0.15,0.15,250,-1.25,1.25);
  TH2D *hdy_HCAL_vs_ptheta = new TH2D("hdy_HCAL_vs_ptheta","",250,sbstheta-0.3,sbstheta+0.3,250,-1.25,1.25);
  
  TH1D *hE_HCAL = new TH1D("hE_HCAL",";HCAL E (GeV);",250,0.0,1.0);
  TH1D *hE_HCAL_cut = new TH1D("hE_HCAL_cut",";HCAL E (GeV);",250,0.0,1.0);
  
  TH1D *hW_cut_HCAL = new TH1D("hW_cut_HCAL",";W (GeV);", 400,0.0,4.0);
  TH1D *hdpel_cut_HCAL = new TH1D("hdpel_cut_HCAL",";p/p_{elastic}(#theta)-1;", 250,-0.5,0.5);

  TH1D *hE_preshower = new TH1D("hE_preshower",";E_{PS} (GeV);",250,0.0,1.25);
  TH1D *hEoverP = new TH1D("hEoverP",";E/p;",250,0.0,2.0);
  TH2D *hEoverP_vs_preshower = new TH2D("hEoverP_vs_preshower",";E_{PS} (GeV);E/p",125,0.0,1.25,125,0.0,2.0);

  TH1D *hE_preshower_cut = new TH1D("hE_preshower_cut",";E_{PS} (GeV);",250,0.0,1.25);
  TH1D *hEoverP_cut = new TH1D("hEoverP_cut",";E/p;",250,0.0,2.0);
  TH2D *hEoverP_vs_preshower_cut = new TH2D("hEoverP_vs_preshower_cut",";E_{PS} (GeV);E/p",125,0.0,1.25,125,0.0,2.0);

  TH2D *hdxdy_HCAL_cut = new TH2D("hdxdy_HCAL_cut",";y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)", 250, -1.25, 1.25, 250, -2.5, 2.5 );

  TH1D *hdx_HCAL_cut = new TH1D("hdx_HCAL_cut",";x_{HCAL}-x_{expect} (m);", 500, -2.5, 2.5);
  TH1D *hdy_HCAL_cut = new TH1D("hdy_HCAL_cut",";y_{HCAL}-y_{expect} (m);", 500, -1.25, 1.25);

  TH1D *hdeltaphi = new TH1D("hdeltaphi",";#phi_{p}-#phi_{e}-#pi;", 250, -0.25,0.25 );
  TH1D *hthetapq = new TH1D("hthetapq",";#theta_{pq};", 250,0.0,0.5);
  TH1D *hpmiss_perp = new TH1D("hpmiss_perp",";p_{miss,#perp} (GeV);", 250,0.0,0.5);
  TH1D *hdeltaptheta = new TH1D("hdeltaptheta",";#theta_{p}-#theta_{p,expect} (rad);", 250, -0.25,0.25);
  
  //For these only apply W and preshower cuts:
  TH1D *hdeltaphi_cut = new TH1D("hdeltaphi_cut",";#phi_{p}-#phi_{e}-#pi;", 250, -0.5,0.5 );
  TH1D *hthetapq_cut = new TH1D("hthetapq_cut",";#theta_{pq};", 250,0.0,0.5);
  TH1D *hpmiss_perp_cut = new TH1D("hpmiss_perp_cut",";p_{miss,#perp} (GeV);", 250,0.0,0.5);
  TH1D *hdeltaptheta_cut = new TH1D("hdeltaptheta_cut",";#theta_{p}-#theta_{p,expect} (rad);", 250, -0.25,0.25);

  TH1D *hetheta_cut = new TH1D("hetheta_cut","#theta_{e} (deg)",250,bbtheta*TMath::RadToDeg()-10.0,bbtheta*TMath::RadToDeg()+10.0);
  TH1D *hephi_cut = new TH1D("hephi_cut","#phi_{e} (deg)",250,-45,45);
  TH1D *hpphi_cut = new TH1D("hpphi_cut","#phi_{p} (deg)",250,135,225);

  double pel_central = ebeam/(1.+ebeam/Mp*(1.-cos(bbtheta)));
  double Q2_central = 2.*ebeam*pel_central*(1.-cos(bbtheta));
  
  TH1D *hep_cut = new TH1D("hep_cut",";p_{e} (GeV);",250,0.5*pel_central,1.5*pel_central);
  TH1D *hQ2_cut = new TH1D("hQ2_cut",";Q^{2} (GeV^{2});",250,0.5*Q2_central, 1.5*Q2_central );
  
  TH1D *hptheta_cut = new TH1D("hptheta_cut",";#theta_{p} (deg);",250,sbstheta*TMath::RadToDeg()-15.0,sbstheta*TMath::RadToDeg()+15.0);

  double pp_central = sqrt(Q2_central*(1.+Q2_central/(4.*Mp*Mp)));
  
  TH1D *hpp_cut = new TH1D("hpp_cut","p_{N,expect} (GeV)",250,0.5*pp_central,1.5*pp_central);

  double Tcentral = sqrt(pow(pp_central,2)+pow(Mp,2))-Mp;
  
  TH1D *hpEkin_cut = new TH1D("hpEkin_cut","E_{kin,expect} (GeV)",250,0.5*Tcentral,1.5*Tcentral);

  TH1D *hEHCALoverEkin = new TH1D("hEHCALoverEkin","E_{HCAL}/E_{kin,expect}",250,0.0,0.5);

  TH2D *hEHCALoverEkin_x = new TH2D("hEHCALoverEkin_x","E_{HCAL}/E_{kin,expect}",50,-1.42,2.16,250,0.0,1.0);

  TH2D *hEHCALoverEkin_y = new TH2D("hEHCALoverEkin_y","E_{HCAL}/E_{kin,expect}",50,-0.84,0.8,250,0.0,1.0);

  TH1D *hEHCALoverEkin_ySlice[ncols_hcal];
  TH1D *hEHCALoverEkin_xSlice[nrows_hcal];
  /*
  for( int i=0; i<nrows_hcal; i++ ){
    if(i%2==0) hEHCALoverEkin_ySlice[i/2] = new TH1D(Form("hEHCALoverEkin_ySlice_c%d",i/2),Form("E_{c%d}/E_{kin,expect} (GeV)",i/2),250,0.0,1.0);
    hEHCALoverEkin_xSlice[i] = new TH1D(Form("hEHCALoverEkin_xSlice_r%d",i),Form("E_{r%d}/E_{kin,expect} (GeV)",i),250,0.0,1.0);
  }
  */
  TH1D *hvz_cut = new TH1D("hvz_cut",";vertex z (m);", 250,-0.125,0.125);

  vector<long> elasEvent;
  long nevent = 0;

  while( C->GetEntry( nevent++ ) ){
    if( nevent % 100000 == 0 ) cout << "NEvents: " << nevent << "/" << totEvents << ", tracks: " << tracks+1 << ", PosCutCount: " << posCutCount << ", BBCutCount " << BBCutCount << ", WCutCount " << WCutCount << ", trigger: " << trig << ", inTime: " << inTime << "." << endl;

    //if( trig!=0 ) cout << "TRIG: " << trig << endl;

    if( ntrack > 0 ){

      tracks++;

      double etheta = acos( pz[0]/p[0] );
      double ephi = atan2( py[0], px[0] );

      TVector3 vertex(0,0,vz[0]);

      TLorentzVector Pbeam(0,0,ebeam,ebeam); //Mass of e negligable
      TLorentzVector kprime(px[0],py[0],pz[0],p[0]);
      TLorentzVector Ptarg(0,0,0,Mp);

      TLorentzVector q = Pbeam - kprime;

      TLorentzVector PgammaN = Ptarg + q; //(-px, -py, ebeam - pz, Mp + ebeam - p)
      
      double pel = ebeam/(1.+ebeam/Mp*(1.-cos(etheta)));

      hdpel->Fill( p[0]/pel - 1.0 );

      hW->Fill( PgammaN.M() );

      hvz->Fill( vertex.Z() );

      double nu = ebeam - p[0];
      double pp = sqrt(pow(nu,2)+2.*Mp*nu);
      double phinucleon = ephi + TMath::Pi(); //assume coplanarity
      double thetanucleon = acos( (ebeam - p[0]*cos(etheta))/pp ); //use elastic constraint on nucleon kinematics
      
      TVector3 pNhat( sin(thetanucleon)*cos(phinucleon),sin(thetanucleon)*sin(phinucleon),cos(thetanucleon));
      
      //Define HCal coordinate system
      TVector3 HCAL_zaxis(-sin(sbstheta),0,cos(sbstheta));
      TVector3 HCAL_xaxis(0,1,0);
      TVector3 HCAL_yaxis = HCAL_zaxis.Cross(HCAL_xaxis).Unit();
      
      TVector3 HCAL_origin = hcaldist * HCAL_zaxis + hcalheight * HCAL_xaxis;
      
      TVector3 TopRightBlockPos_DB(xoff_hcal,yoff_hcal,0);
      
      TVector3 TopRightBlockPos_Hall( hcalheight + (nrows_hcal/2-0.5)*blockspace_hcal,
				      (ncols_hcal/2-0.5)*blockspace_hcal, 0 );
      
      
      //Assume that HCAL origin is at the vertical and horizontal midpoint of HCAL and locate cluster position
      xHCAL += TopRightBlockPos_Hall.X() - TopRightBlockPos_DB.X();
      yHCAL += TopRightBlockPos_Hall.Y() - TopRightBlockPos_DB.Y();
      
      double sintersect = ( HCAL_origin - vertex ).Dot( HCAL_zaxis ) / (pNhat.Dot( HCAL_zaxis ) );
      
      TVector3 HCAL_intersect = vertex + sintersect * pNhat;
      
      double yexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_yaxis );
      double xexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_xaxis );
      
      hdx_HCAL->Fill( xHCAL - xexpect_HCAL );
      hdy_HCAL->Fill( yHCAL - yexpect_HCAL );
      
      //Calculate the proton spot
      hdxdy_HCAL->Fill( yHCAL - yexpect_HCAL, xHCAL - xexpect_HCAL );
      
      //Plot correlation between obs and exp pos of nucleon
      hxcorr_HCAL->Fill( xexpect_HCAL, xHCAL );
      hycorr_HCAL->Fill( yexpect_HCAL, yHCAL );
      
      //Plot other correlations with horizontal pos of nucleon
      hdy_HCAL_vs_z->Fill( vertex.Z(), yHCAL - yexpect_HCAL );
      hdy_HCAL_vs_ptheta->Fill( thetanucleon, yHCAL - yexpect_HCAL );

      //Construct exp nucleon vectors
      TVector3 HCALpos(xHCAL,yHCAL,0);
      TVector3 HCALpos_global = HCAL_origin + HCALpos.X() * HCAL_xaxis + HCALpos.Y() * HCAL_yaxis;

      TVector3 HCAL_ray = HCALpos_global - vertex;
      HCAL_ray = HCAL_ray.Unit();

      //Get angles of nucleon vector wrt beamline
      double ptheta_recon = HCAL_ray.Theta();
      double pphi_recon = HCAL_ray.Phi();

      if( pphi_recon < 0.0 ) pphi_recon += 2.0*TMath::Pi();
      if( pphi_recon > 2.0*TMath::Pi() ) pphi_recon -= 2.0*TMath::Pi();

      //Plot exp - obs angles
      hdeltaphi->Fill( pphi_recon - phinucleon );
      hdeltaptheta->Fill( ptheta_recon - thetanucleon );
      
      double thetapq = acos( HCAL_ray.Dot( pNhat ) );
      hthetapq->Fill( thetapq );

      TVector3 pNrecon = pp*HCAL_ray;

      //Get energy of nucleon
      double Enucleon = sqrt(pow(pp,2)+pow(Mp,2));

      TLorentzVector PNrecon( pNrecon,Enucleon );

      TLorentzVector Pmiss = PgammaN - PNrecon;

      double pmiss_perp = (Pmiss.Vect() - Pmiss.Vect().Dot(q.Vect().Unit())*q.Vect().Unit()).Mag();

      hpmiss_perp->Fill( pmiss_perp );
      
      hE_HCAL->Fill( EHCAL );

      //Nucleon position cut
      if( pow( (xHCAL-xexpect_HCAL - dx0)/dxsigma,2) + pow( (yHCAL-yexpect_HCAL - dy0)/dysigma,2) <= pow(2.5,2) ){

	posCutCount++;

	//BB sh and ps cuts (E/p)
	if( Eps_BB >= 0.15 && abs( (Eps_BB+Esh_BB)/p[0] - 1. ) < 0.25 ){

	  BBCutCount++;

	  Double_t bbcal_time=0., hcal_time=0.;
	  for(Int_t ihit=0; ihit<Ndata_bb_tdctrig_tdcelemID; ihit++){
	    if(bb_tdctrig_tdcelemID[ihit]==5) bbcal_time=bb_tdctrig_tdc[ihit];
	    if(bb_tdctrig_tdcelemID[ihit]==0) hcal_time=bb_tdctrig_tdc[ihit];
	  }
	  Double_t BB_H_diff = hcal_time - bbcal_time; 
	  /*
	  if(fabs(diff-510.)<20.){
	    h2_bbcal_hcal_corr->Fill( bb_sh_rowblk+1, sbs_hcal_rowblk+1);
	  }
	  */


	  hW_cut_HCAL->Fill( PgammaN.M() );
	  hdpel_cut_HCAL->Fill( p[0]/pel - 1.0 );
	  
	  
	  hE_HCAL_cut->Fill( EHCAL );

	  //Cuts on elastics via W and p/p'
	  if( Wmin < PgammaN.M() && PgammaN.M() < Wmax && dpel_min < p[0]/pel-1.&&p[0]/pel-1. < dpel_max ){
	      
	    WCutCount++;
	    elasEvent.push_back(nevent);

	    hep_cut->Fill( p[0] );
	    hQ2_cut->Fill( 2.*ebeam*p[0]*(1.-cos(etheta)) );
	    hptheta_cut->Fill( thetanucleon * TMath::RadToDeg() );
	    hetheta_cut->Fill( etheta* TMath::RadToDeg() );
	    hpp_cut->Fill( pp );
	    hpEkin_cut->Fill( nu );

	    //Uniformity plots
	    hEHCALoverEkin->Fill( EHCAL/nu );
	    hEHCALoverEkin_x->Fill( xHCAL, EHCAL/nu );
	    hEHCALoverEkin_y->Fill( yHCAL, EHCAL/nu );

	    //HCal/Hodo timing correlation
	    if( EHCAL>0.02 && tHCAL_blk[0]>-400 && nClusHODO==1 ){
	      htcorr_HCAL_HODO->Fill( tHCAL_blk[0], tHODO );
	      htcorr_HCAL_BBTrig->Fill( tHCAL_blk[0], bbcal_time );

	      //cout << bbcal_time << endl;

	      double HHdiff = tHCAL_blk[0]-tHODO;
	      double BBHdiff = bbcal_time-tHCAL_blk[0];
	      htDiff_HODO_HCAL->Fill( HHdiff );
	      htDiff_vs_HCALID->Fill( tdcID[0], HHdiff );

	      //cout << BBHdiff << endl;

	      htBBHDiff_vs_HCALID->Fill( tdcID[0], BBHdiff );

	      int tdc = (int)tdcID[0];

	      //cout << (int)tdcID[0] << " " << tdcID[0] << endl;

	      if( tdc<0 || tdc>288 ) cout << "ERROR: TDC out of bounds at tdcID = " << tdc << endl;

	      //if( aHCAL[0]<0 || aHCAL[0]>0.3 || HHdiff<-150 || HHdiff>0 ) cout << "E=" << aHCAL[0] << " tDiff=" << HHdiff << endl;

	      // Cut out afterpulse with if statement
	      if( HHdiff<-65. && HHdiff>-80. ) htDiff_vs_ADCint[tdc]->Fill( aHCAL[0], HHdiff );

	      //if( tdcID[0]==175 ) htDiff_vs_ADCint_bl175->Fill( aHCAL[0], HHdiff );
	      //if( tdcID[0]==224 ) htDiff_vs_ADCint_bl224->Fill( aHCAL[0], HHdiff );
	      //if( tdcID[0]==92 ) htDiff_vs_ADCint_bl92->Fill( aHCAL[0], HHdiff );

	      if( HHdiff<-65. && HHdiff>-80. ) htDiff_vs_HCALID_corr->Fill( tdcID[0], HHdiff+22.0*aHCAL[0] ); //Timewalk corrected from linear fit to time diff vs clus e

	      inTime++;
	    }
	    hephi_cut->Fill( ephi * TMath::RadToDeg() );
	    hpphi_cut->Fill( pphi_recon * TMath::RadToDeg() );

	    hvz_cut->Fill( vz[0] );

	  }
	}

	//cut on elastic peak for EoverP and preshower cut plots:
	if( Wmin < PgammaN.M() && PgammaN.M() < Wmax && dpel_min < p[0]/pel-1.&&p[0]/pel-1. < dpel_max ){ 
	  hE_preshower_cut->Fill( Eps_BB );
	  hEoverP_cut->Fill( (Eps_BB+Esh_BB)/p[0] );
	  hEoverP_vs_preshower_cut->Fill( Eps_BB,  (Eps_BB+Esh_BB)/p[0] );
	  
	}

      }

      if( Eps_BB >= 0.15 && abs( (Eps_BB+Esh_BB)/p[0] - 1.0 ) <= 0.3 ){
	hdpel_cutBBCAL->Fill( p[0]/pel - 1.0 );
	hW_cutBBCAL->Fill( PgammaN.M() );
	if( Wmin <= PgammaN.M() && PgammaN.M() <= Wmax && 
	    dpel_min <= p[0]/pel - 1. && p[0]/pel - 1. < dpel_max ){
	  hdx_HCAL_cut->Fill( xHCAL - xexpect_HCAL );
	  hdy_HCAL_cut->Fill( yHCAL - yexpect_HCAL );
	  hdxdy_HCAL_cut->Fill( yHCAL - yexpect_HCAL, xHCAL - xexpect_HCAL );

	  hdeltaphi_cut->Fill( pphi_recon - phinucleon );
	  hthetapq_cut->Fill( thetapq );
	  hpmiss_perp_cut->Fill( pmiss_perp );
	  hdeltaptheta_cut->Fill( ptheta_recon - thetanucleon );
	}
      }

      hE_preshower->Fill( Eps_BB );
      hEoverP->Fill( (Eps_BB+Esh_BB)/p[0] );
      hEoverP_vs_preshower->Fill( Eps_BB,  (Eps_BB+Esh_BB)/p[0] );
    }
  }

  TF1 *expo_fit[ncells_hcal] = {};
  double slope[ncells_hcal] = {0.0};

  for( int c=0; c<ncells_hcal; c++ ){
    TF1 *f1;

    if( htDiff_vs_ADCint[c]->GetEntries()>30 ){
      //expo_fit[c] = new TF1(Form("expo_fit cell%d",c), expoFit, 0, 1000, 3);
      //expo_fit[c]->SetLineColor(4);
      //expo_fit[c]->SetNpx(1000);
      //htDiff_vs_ADCint[c]->Fit(expoFit[c],"+RQ");
      //htDiff_vs_ADCint[c]->SetTitle(Form("htDiff_vs_ADCint_bl%d p0%d p1%d p2%d",c,expoFit[c]->GetParameter(0),expoFit[c]->GetParameter(1),expoFit[c]->GetParameter(2)));

      htDiff_vs_ADCint[c]->Fit("pol1","Q","",0.0,0.2);
      f1=htDiff_vs_ADCint[c]->GetFunction("pol1");
      htDiff_vs_ADCint[c]->SetName(Form("FITTED_htDiff_vs_ADCint_bl%d",c));
      htDiff_vs_ADCint[c]->SetTitle(Form("bl%d offs%f sl%f",c,f1->GetParameter(0),f1->GetParameter(1)));

      slope[c] = f1->GetParameter(1);
    }

  }

  int totalElastics = elasEvent.size();

  cout << endl << "Looping back over elastics to apply time-walk corrections with " << totalElastics << " elastics events." << endl << endl;

  nevent = 0;

  while( C->GetEntry( elasEvent[nevent]-1 ) ){

    if( elasEvent[nevent] < 1 ){
      cout << "Last good event: " << elasEvent[nevent-1] << " at event in loop: " << nevent << "." << endl;
      break;
    }

    int id = (int)tdcID[0];

    if( nevent % 10000 == 0 ) cout << "ElasEvent: " << nevent << "/" << totalElastics << " at overall event: " << elasEvent[nevent] << "/" << totEvents << ". tdcID: " << tdcID[0] << ". id: " << id << "." << endl;

    if( EHCAL>0.02 && tHCAL_blk[0]>-400 && nClusHODO==1 ){
      
      double HHdiff = tHCAL_blk[0]-tHODO;
      
      if( HHdiff<-65. && HHdiff>-80. ) htDiff_vs_HCALID_corrByChannel->Fill( tdcID[0], HHdiff+slope[id]*aHCAL[0] ); //Timewalk corrected from linear fit to time diff vs clus e
      //if( HHdiff<-65. && HHdiff>-80. ) htDiff_vs_HCALID_corrByChannel->Fill( tdcID[0], HHdiff+22.0*aHCAL[0] ); 
    }
    nevent++;
  }


  for( int i=0; i<nrows_hcal; i++ ){
    TF1 *f2;

    if(i%2==0) {
      hEHCALoverEkin_ySlice[i/2] = new TH1D(Form("hEHCALoverEkin_ySlice_c%d",i/2),Form("E_{c%d}/E_{kin,expect} (GeV)",i/2),250,0.0,1.0);
      for(int j=0; j<250; j++){
	hEHCALoverEkin_ySlice[i/2]->SetBinContent(j, hEHCALoverEkin_y->GetBinContent(i/2+1,j+1));
      }
      hEHCALoverEkin_ySlice[i/2]->Fit("gaus","Q","",0.01,0.1);
      f2=hEHCALoverEkin_ySlice[i/2]->GetFunction("gaus");
      //hEHCALoverEkin_ySlice[i/2]->SetName(Form("FITTED_hEHCALoverEkin_ySlice_c%d",i/2));
      //hEHCALoverEkin_ySlice[i/2]->SetTitle(Form("bl%d amp%f pos%f stddev%f",i,f2->GetParameter(0),f2->GetParameter(1),f2->GetParameter(2)));
    }
    TF1 *f3;
    hEHCALoverEkin_xSlice[i] = new TH1D(Form("hEHCALoverEkin_xSlice_r%d",i),Form("E_{r%d}/E_{kin,expect} (GeV)",i),250,0.0,1.0);
    for(int j=0; j<250; j++){
      hEHCALoverEkin_xSlice[i]->SetBinContent(j, hEHCALoverEkin_x->GetBinContent(i+1,j+1));
    }
    hEHCALoverEkin_xSlice[i]->Fit("gaus","Q","",0.01,0.1);
    f3=hEHCALoverEkin_xSlice[i]->GetFunction("gaus");
    //hEHCALoverEkin_xSlice[i]->SetName(Form("FITTED_hEHCALoverEkin_xSlice_c%d",i));
    //hEHCALoverEkin_xSlice[i]->SetTitle(Form("bl%d amp%f pos%f stddev%f",i,f3->GetParameter(0),f3->GetParameter(1),f3->GetParameter(2)));
  }
  
  double charge = 0.0;
  for( int i=0; i<8; i++){
    charge += uA_hr_SBS8_tot[i];
  }

  cout << "Charge normalized elastic yield: " << WCutCount/charge << " elastics/uA-hr." << endl;

  fout->Write();

}
