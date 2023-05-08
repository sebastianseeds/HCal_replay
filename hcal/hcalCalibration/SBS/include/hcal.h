#ifndef HCAL_H
#define HCAL_H

#include <cmath>
#include "TVector3.h"
#include "TLorentzVector.h"

const Int_t NTARGET = 3;
const Int_t NSTAMP = 4;
const Int_t NTSTAMP = 10;

namespace hcal{
  //HCal Detector
  const Int_t maxTracks = 1000; // Reasonable limit on tracks to be stored per event
  //HCal - Note that actual measurement of vertical length is 381.6cm, indicating that the MC figures are correct
  const Int_t maxHCalChan = 288; // Total HCal channels
  const Int_t maxHCalRows = 24; // Total HCal rows
  const Int_t maxHCalCols = 12; // Total HCal cols
  const Int_t maxClusters = 10; // Total HCal clusters with information saved
  const Double_t HCalvoff = -0.2897; // Height of HCal above beamline in m
  const Double_t HCalblk_l_h = 0.15494; // Horizontal length of all HCAL blocks in m from MC database
  const Double_t HCalblk_l_v = 0.15875; // Vertical length of all HCAL blocks in m from MC database
  const Double_t posHCalXi = -2.3531; // Distance from beam center to top of HCal in m from MC database
  const Double_t posHCalXf = 1.45309; // Distance from beam center to bottom of HCal in m from MC database
  const Double_t posHCalYi = -0.93155; // Distance from beam center to opposite-beam side of HCal in m from MC database
  const Double_t posHCalYf = 0.93155; // Distance from beam center to beam side of HCal in m from MC database
  const Double_t HCalSampFrac = 0.077;  //Re-evaluated with MC GEn settings using second to outermost shower column for kin2
  //Trigger TDC
  const Int_t maxTDCTrigChan = 10; // Set to accomodate original 5 TDCTrig channels with buffer
  const Double_t tdiffwidecut = 50; // Set at 50ns from nominal 510ns (passed by user) from GMn
  const Int_t gNstamp = NSTAMP;
  const string gStamp[NSTAMP] = { "gain", "tdcoff", "tdccalib", "adcoff" }; //Ordered list timestamp types

  //BBCal
  const Int_t maxBBCalShChan = 189; // Total BBCal Shower Channels
  const Int_t maxBBCalShRows = 27;
  const Int_t maxBBCalShCols = 7;
  const Int_t maxBBCalPSChan = 52; // Total BBCal Preshower Channels
  const Int_t maxBBCalPSRows = 26;
  const Int_t maxBBCalPSCols = 2;

  //Beamline and target
  const Int_t chargeConvert = 3318; // See D.Flay Doc DB sbs.jlab.org/DocDB/0001/000164/002/dflay_bcm-ana-update_02-21-22.pdf p.8
  const Int_t clockActual = 103700; // Needed to convert the 104kHz clock to the actual counting rate
  const Int_t gMaxfset = 5; //Reasonable maximum on total SBS field settings per configuration
  const Int_t gNtar = 3; //Total number of possible targets in SBS (LH2, LD2, Helium-3)
  const string gTar[NTARGET] = { "lh2", "ld2", "he3" }; //Ordered list of available targets
  
  //SBS Magnet
  const Double_t Dgap = 48.0*2.54/100.0; //about 1.22 m
  const Double_t maxSBSfield = 1.26; //Tesla
  const Double_t SBSfield = 1.0; //fraction of max field. TODO: should be variable per run
  const Double_t SBSdist = 2.25; //m
  const Double_t sbsdipolegap = 48.0*2.54/100.;  // ~1.22 m
  const Double_t sbsmaxfield = 3.1 * atan( 0.85/(11.0 - 2.25 - 1.22/2.0 ))/0.3/1.22/0.7;
  //GEMs
  const Double_t GEMpitch = 10*TMath::DegToRad();

  ///////////////
  ///Physics/Math
  const Double_t PI = TMath::Pi();
  const Double_t M_t[NTARGET] = {0.938272,(0.938272+0.939565)/2,(0.938272+0.939565+0.939565)/3};
  const Double_t M_e = 0.00051;
  const Double_t M_p = 0.938272;
  const Double_t M_n = 0.939565;
  const UInt_t us = 1000000; //For conversion to seconds used by reporting time delays
  const Int_t protonMCidx = 2212; //pdg.lbl.gov/2020/reviews/rpp2020-rev-monte-carlo-numbering.pdf
  const Int_t neutronMCidx = 2112; //pdg.lbl.gov/2020/reviews/rpp2020-rev-monte-carlo-numbering.pdf
  const Int_t electronMCidx = 11; //pdg.lbl.gov/2020/reviews/rpp2020-rev-monte-carlo-numbering.pdf
  const Int_t photonMCidx = 22; //pdg.lbl.gov/2020/reviews/rpp2020-rev-monte-carlo-numbering.pdf

  ////////////////////////////
  ////Static Target/Scattering Chamber Parameters (LH2,LD2,He3)
  const Double_t l_tgt[NTARGET] = {0.15,0.15,0.3}; // Length of the target (m)
  const Double_t tarrho[NTARGET] = {0.0723,0.169,0.26}; //g/cc, target density, He3 is a guess
  const Double_t crho[NTARGET] = {2.71,2.71,2.7}; // Density of {Al,Al,glass}
  const Double_t cdiam[NTARGET] = {0.04064,0.04064,0.04064}; //m
  const Double_t cdEdx[NTARGET] = {0.0021,0.0021,0.0021}; // guess
  const Double_t cthick[NTARGET] = {0.02,0.02,0.5};  //cm, target cell thickness, He3 is a guess
  const Double_t uwallthick[NTARGET] = {0.0145,0.0145,0.5}; //cm, upstream wall thickness
  const Double_t dwallthick[NTARGET] = {0.015,0.015,0.5};  //cm, downstream wall thickness
  const Double_t dEdx[NTARGET] = {0.00574,0.01148,0.01722}; //According to NIST ESTAR, the collisional stopping power of hydrogen is about 5.74 MeV*cm2/g at 2 GeV energy, else guess

  const Double_t celldiameter = 1.6*2.54/100; //m, guess
  const Double_t dEdx_Al = 0.0021; //According to NIST ESTAR, the collisional stopping power of Aluminum is about 2.1 MeV*cm2/g between 1-4 GeV
  const Double_t dEdx_glass = 0.0021; // guess
  const Double_t Alshieldthick = 2.54/8.0; //= 1/8 inch * 2.54 cm/inch 

  //////////////////////////////
  ////GMn Experimental params
  Int_t gIdx[6] = {4, 7, 11, 14, 8, 9}; //Index of SBS run
  TCut gglobalcut[6] = {"bb.tr.n==1&&bb.ps.e>0.2&&abs(bb.tr.vz[0])<0.08&&bb.gem.track.nhits>3&&abs(bb.etot_over_p-0.92)<0.2&&sbs.hcal.e>0.01&&bb.ps.e+bb.sh.e>1.7", 
			"bb.tr.n>0&&abs(bb.tr.vz[0])<0.08&&bb.gem.track.nhits>2&&bb.tr.p[0]>2.0&&sbs.hcal.e>0.01&&bb.ps.e>0.2", 
			"bb.tr.n>0&&abs(bb.tr.vz[0])<0.08&&bb.gem.track.nhits>2&&bb.tr.p[0]>2.0&&sbs.hcal.e>0.01&&bb.ps.e>0.2", 
			"bb.tr.n>0&&abs(bb.tr.vz[0])<0.08&&bb.gem.track.nhits>2&&bb.tr.p[0]>1.6&&sbs.hcal.e>0.01&&bb.ps.e>0.2", 
			"bb.tr.n==1&&bb.ps.e>0.2&&abs(bb.tr.vz[0])<0.08&&bb.gem.track.nhits>3&&sbs.hcal.e>0.01&&abs(bb.tr.tg_th[0])<0.15&&abs(bb.tr.tg_ph[0])<0.3", 
			"bb.tr.n==1&&bb.ps.e>0.2&&abs(bb.tr.vz[0])<0.08&&bb.gem.track.nhits>3&&sbs.hcal.e>0.01&&abs(bb.tr.tg_th[0])<0.15&&abs(bb.tr.tg_ph[0])<0.3"};
  TCut gglobcut_earm[6] = {"bb.tr.n==1&&bb.ps.e>0.2&&abs(bb.tr.vz[0])<0.08&&bb.gem.track.nhits>3&&abs(bb.etot_over_p-0.92)<0.2&&bb.ps.e+bb.sh.e>1.7", 
			   "bb.tr.n>0&&abs(bb.tr.vz[0])<0.08&&bb.gem.track.nhits>2&&bb.tr.p[0]>2.0&&bb.ps.e>0.2", 
			   "bb.tr.n>0&&abs(bb.tr.vz[0])<0.08&&bb.gem.track.nhits>2&&bb.tr.p[0]>2.0&&bb.ps.e>0.2", 
			   "bb.tr.n>0&&abs(bb.tr.vz[0])<0.08&&bb.gem.track.nhits>2&&bb.tr.p[0]>1.6&&bb.ps.e>0.2", 
			   "bb.tr.n==1&&bb.ps.e>0.2&&abs(bb.tr.vz[0])<0.08&&bb.gem.track.nhits>3&&abs(bb.tr.tg_th[0])<0.15&&abs(bb.tr.tg_ph[0])<0.3", 
			   "bb.tr.n==1&&bb.ps.e>0.2&&abs(bb.tr.vz[0])<0.08&&bb.gem.track.nhits>3&&abs(bb.tr.tg_th[0])<0.15&&abs(bb.tr.tg_ph[0])<0.3"};
  Double_t gE_e[6] = {3.7278, 7.9072, 9.8594, 5.9649, 5.9648, 4.0148}; //GeV
  Double_t gHCal_d[6] = {11., 14., 14.5, 14., 11., 11.}; //m
  Double_t gHCal_th[6] = {31.9, 16.1, 13.3, 17.3, 29.4, 22.}; //degrees
  Double_t gBB_th[6] = {36., 40., 42., 46.5, 26.5, 49.}; //degrees
  //The remaining parameters are good starting points, but should be configured by field. dx positions are left out due to strong dependence on field.
  Double_t gW2_mean[6] = {0.917994, 0.88, 0.925, 0.870819, 0.91, 0.91}; //GeV
  Double_t gW2_sig[6] = {0.167922, 0.5, 0.325, 0.19443, 0.21, 0.17}; //GeV
  Double_t gdy0[6] = {-0.0270143, 0.0158817, -0.0119, 0.00220912, -0.0119, -0.0175041}; //m
  Double_t gdy_sig[6] = {0.0869443, 0.155968, 0.08675, 0.0824202, 0.08675, 0.134702}; //m
  Double_t gatime0[6] = {51.466, 63.1685, 50.36, 58.7825, 50.36, 50.6158}; //ns
  Double_t gatime_sig[6] = {3.67744, 3.69617, 3.73, 3.68101, 3.73, 3.52261}; //ns
  Double_t guseAlshield[6] = {0, 0, 0, 0, 0, 0}; //0==out, 1==in

  //Utility
  //Get today's date
  string getDate(){
    time_t now = time(0);
    tm ltm = *localtime(&now);
  
    string yyyy = to_string(1900 + ltm.tm_year);
    string mm = to_string(1 + ltm.tm_mon);
    string dd = to_string(ltm.tm_mday);
    string date = mm + '_' + dd + '_' + yyyy;
  
    return date;
  }
  //sets hcal axes by kinematic (from hall coordinates)
  void sethcalaxes( Double_t sbstheta_rad,                           // SBS angle in radian 
		    vector<TVector3> &hcal_axes ) {
    TVector3 hcal_zaxis( sin(-sbstheta_rad), 0, cos(-sbstheta_rad) ); // Clock-wise rotation about Y axis
    TVector3 hcal_xaxis( 0, -1, 0 );                                  // -Y axis of Hall CoS = X axis of hcal CoS
    TVector3 hcal_yaxis = hcal_zaxis.Cross(hcal_xaxis).Unit();
    hcal_axes.push_back(hcal_xaxis);
    hcal_axes.push_back(hcal_yaxis);
    hcal_axes.push_back(hcal_zaxis);
  }
  //gets expected location of scattered nucleon assuming straight line projections from BB track
  void getxyhcalexpect( TVector3 vertex, TVector3 pNhat, TVector3 hcal_origin, 
			vector<TVector3> hcal_axes, vector<Double_t> &xyhcalexpect ) {
    // intersection of a ray with a plane in hcal coordinates
    Double_t sintersect = ( hcal_origin - vertex).Dot(hcal_axes[2] ) / ( pNhat.Dot(hcal_axes[2]) );
    // ray from Hall origin onto the face of hcal where the nucleon hit
    TVector3 hcal_intersect = vertex + sintersect*pNhat; 

    Double_t xexpect_hcal = ( hcal_intersect - hcal_origin ).Dot( hcal_axes[0] );
    Double_t yexpect_hcal = ( hcal_intersect - hcal_origin ).Dot( hcal_axes[1] );

    xyhcalexpect.push_back( xexpect_hcal );
    xyhcalexpect.push_back( yexpect_hcal );
  }


  //Fit functions
  //gaussian fit
  Double_t g_gfit(Double_t *x, Double_t *par){
    Double_t amp = par[0];
    Double_t offset = par[1];
    Double_t sigma = par[2];
    return amp*exp(-0.5*pow((x[0]-offset)/sigma,2.));
  }

  //expo fit
  Double_t g_expofit(Double_t *x, Double_t *par){
    Double_t amp = par[0];
    Double_t str = par[1];
    return exp(amp+str*x[0]);
  }

  //expo fit with offset
  Double_t g_scexpofit(Double_t *x, Double_t *par){
    Double_t yint = par[0];
    Double_t amp = par[1];
    Double_t str = par[2];
    return yint+exp(amp+str*x[0]);
  }

  //linear fit
  Double_t g_lfit(Double_t *x, Double_t *par){
    Double_t yint = par[0];
    Double_t p1 = par[1];;
    return yint+p1*x[0];
  }

  //2nd order poly fit
  Double_t g_p2fit(Double_t *x, Double_t *par){
    Double_t yint = par[0];
    Double_t p1 = par[1];
    Double_t p2 = par[2];
    return yint+p1*x[0]+p2*pow(x[0],2);
  }

  //3rd order poly fit
  Double_t g_p3fit(Double_t *x, Double_t *par){
    Double_t yint = par[0];
    Double_t p1 = par[1];
    Double_t p2 = par[2];
    Double_t p3 = par[3];
    return yint+p1*x[0]+p2*pow(x[0],2)+p3*pow(x[0],3);
  }

  //4th order poly fit
  Double_t g_p4fit(Double_t *x, Double_t *par){
    Double_t yint = par[0];
    Double_t p1 = par[1];
    Double_t p2 = par[2];
    Double_t p3 = par[3];
    Double_t p4 = par[4];
    return yint+p1*x[0]+p2*pow(x[0],2)+p3*pow(x[0],3)+p4*pow(x[0],4);
  }

  //5th order poly fit
  Double_t g_p5fit(Double_t *x, Double_t *par){
    Double_t yint = par[0];
    Double_t p1 = par[1];
    Double_t p2 = par[2];
    Double_t p3 = par[3];
    Double_t p4 = par[4];
    Double_t p5 = par[5];
    return yint+p1*x[0]+p2*pow(x[0],2)+p3*pow(x[0],3)+p4*pow(x[0],4)+p5*pow(x[0],5);
  }

  //6th order poly fit
  Double_t g_p6fit(Double_t *x, Double_t *par){
    Double_t yint = par[0];
    Double_t p1 = par[1];
    Double_t p2 = par[2];
    Double_t p3 = par[3];
    Double_t p4 = par[4];
    Double_t p5 = par[5];
    Double_t p6 = par[6];
    return yint+p1*x[0]+p2*pow(x[0],2)+p3*pow(x[0],3)+p4*pow(x[0],4)+p5*pow(x[0],5)+p6*pow(x[0],6);
  }

  //7th order poly fit
  Double_t g_p7fit(Double_t *x, Double_t *par){
    Double_t yint = par[0];
    Double_t p1 = par[1];
    Double_t p2 = par[2];
    Double_t p3 = par[3];
    Double_t p4 = par[4];
    Double_t p5 = par[5];
    Double_t p6 = par[6];
    Double_t p7 = par[7];
    return yint+p1*x[0]+p2*pow(x[0],2)+p3*pow(x[0],3)+p4*pow(x[0],4)+p5*pow(x[0],5)+p6*pow(x[0],6)+p7*pow(x[0],7);
  }

  //8th order poly fit
  Double_t g_p8fit(Double_t *x, Double_t *par){
    Double_t yint = par[0];
    Double_t p1 = par[1];
    Double_t p2 = par[2];
    Double_t p3 = par[3];
    Double_t p4 = par[4];
    Double_t p5 = par[5];
    Double_t p6 = par[6];
    Double_t p7 = par[7];
    Double_t p8 = par[8];
    return yint+p1*x[0]+p2*pow(x[0],2)+p3*pow(x[0],3)+p4*pow(x[0],4)+p5*pow(x[0],5)+p6*pow(x[0],6)+p7*pow(x[0],7)+p8*pow(x[0],8);
  }

}

#endif // HCAL_H
