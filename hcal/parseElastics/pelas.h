#ifndef PELAS_H
#define PELAS_H

//Static Detector Parameters
const Int_t maxTracks = 1000; // Reasonable limit on tracks to be stored per event
const Int_t maxTdcChan = 10; // Set to accomodate original 5 TDCTrig channels with buffer
const Int_t maxHCalChan = 288; // Total HCal channels
const Int_t maxBBCalShChan = 189; // Total BBCal Shower Channels
const Int_t maxBBCalPSChan = 52; // Total BBCal Preshower Channels
const Double_t hcalheight = 0.365; // Height of HCal above beamline

//Physics/Math
const Double_t PI = TMath::Pi();
const Double_t M_e = 0.00051;
const Double_t M_p = 0.938272;
const Double_t M_n = 0.939565;

//Static Target Parameters
const Double_t l_tgt = 0.15; // Length of the target (m)
const Double_t rho_tgt = 0.0723; // Density of target (g/cc)
const Double_t rho_Al = 2.7; // Density of aluminum windows (g/cc)
const Double_t celldiameter = 1.6*2.54; //cm, right now this is a guess
const Double_t Ztgt = 1.0;
const Double_t Atgt = 1.0;
const Double_t Mmol_tgt = 1.008; //g/mol

//For energy-loss correction to beam energy:
const Double_t dEdx_tgt=0.00574; //According to NIST ESTAR, the collisional stopping power of hydrogen is about 5.74 MeV*cm2/g at 2 GeV energy
const Double_t dEdx_Al = 0.0021; //According to NIST ESTAR, the collisional stopping power of Aluminum is about 2.1 MeV*cm2/g between 1-4 GeV
const Double_t uwallthick_LH2 = 0.0145; //cm
const Double_t dwallthick_LH2 = 0.015; //cm
const Double_t cellthick_LH2 = 0.02; //cm, this is a guess;
const Double_t Alshieldthick = 2.54/8.0; //= 1/8 inch * 2.54 cm/inch 

#endif // PELAS_H
