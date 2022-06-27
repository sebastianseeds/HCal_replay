// P.Datta, S.Seeds. 10.12.21 This script will generate energy of the elastically scattered electron and proton (e-p)
// taking into consideration the angular acceptance of BBSH for a stationary proton target. Takes three inputs: Beam
// energy, BB angle, & BB distance

#include <TMath.h>
#include <iostream>
#include <vector>

using namespace std;

static const Double_t Mp = 0.938; // GeV

int GMnElasPeak_p(){

  Double_t E_beam=0., BB_angle=0., BB_distance=0.;
  cout << " E_beam(GeV)? BB_angle(deg)? BB_distance(m)?" << endl;
  cin >> E_beam >> BB_angle >> BB_distance;

  Double_t sh_faceDist = 1. + BB_distance;
  
  Double_t sh_ypos[7] = {-0.2565, -0.171, -.0855, 0.0, 0.0855, 0.171, 0.2565};
  Double_t effective_BBang[7] = {0.};

  Double_t deltaBBang = 0.;
  for(int shcol=0; shcol<7; shcol++){
    effective_BBang[shcol] = (sh_ypos[shcol]/sh_faceDist) + BB_angle*(TMath::Pi()/180.); //rad
    cout << "SH Col= " << shcol+1 << " Effective BB angle= " << effective_BBang[shcol]*(180./TMath::Pi()) << " deg" << endl;
  }

  cout << endl;
  
  Double_t Q2 = 4.0;

  Double_t elasPeak[7] = {0.};
  Double_t elasPeak_p[7] = {0.};
  Double_t epsilon[7] = {0.};
  for(int shcol=0; shcol<7; shcol++){
    elasPeak[shcol] = E_beam/( 1. + (2.*E_beam/Mp)*pow(TMath::Sin(effective_BBang[shcol]/2.), 2.) );
    epsilon[shcol] = pow(1+2*(1+Q2/(4*Mp))*pow(TMath::Tan(effective_BBang[shcol]/2),2),-1);
    Double_t nu = E_beam-elasPeak[shcol];
    elasPeak_p[shcol] = nu+Mp;
    cout << "SH Col= " << shcol+1 << " Elastic Peak= " << elasPeak[shcol] << " GeV. Epsilon= " << epsilon[shcol] << "." << endl;
  }

  cout << endl;

  for(int shcol=0; shcol<7; shcol++){
    cout << "SH Col= " << shcol+1 << ". Elastic proton Peak= " << elasPeak_p[shcol] << " GeV" << endl;
  }

  return 0;
}

