//SSeeds 7.26.23 Script to cut on elastics and build hcal energy/sampling-fraction distributions and timing offsets from local SBS-offline (including new gain coefficients and timing offset parameters) for final quality checks.

#include <vector>
#include <iostream>
#include <algorithm>
#include "TCut.h"
#include "TLorentzVector.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TVector3.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TTreeFormula.h"
#include "TChain.h"
#include "TLatex.h"
#include "TH1D.h"
#include "TStopwatch.h"
#include "../../include/sbs.h"

const Double_t bins_SFE = 400.;
const Double_t start_SFE = 0.;
const Double_t end_SFE = 1.;
const Double_t bins_W2 = 200.;
const Double_t start_W2 = 0.;
const Double_t end_W2 = 2.;
const Double_t bins_dy = 100.;
const Double_t start_dy = -2.;
const Double_t end_dy = 2.;
const Int_t atimeNsig = 6;

//Main <experiment> <configuration>
void hcalE_digmc( const char *experiment = "gmn", Int_t config=4, Int_t pass=0, Int_t run=13747 ){
  
  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );
    
  // Get the date
  string date = util::getDate();

  // outfile path
  std::string outdir_path = gSystem->Getenv("OUT_DIR");
  std::string qcheck_path = outdir_path + Form("/hcal_calibrations/pass%d/energy/hcalE_qcheck_%s_sbs%d_p%d.root",pass,experiment,config,pass);

  //Set up path variables and output files
  TFile *fout = new TFile( qcheck_path.c_str(), "RECREATE" );
  
  // Get information from .csv files
  std::string struct_dir = Form("../../config/%s/",experiment); //unique to my environment for now
  Int_t nruns = -1; //Analyze all available runs for this configuration
  
  // re-allocate memory at each run to load different cuts/parameters
  TChain *C = nullptr;

  //Get experimental configuration parameters
  SBSconfig config_parameters(experiment,config);    
  cout << config_parameters;
  Double_t ebeam = config_parameters.GetEbeam();
  Double_t hcaltheta_rad = config_parameters.GetHCALtheta_rad();
  Double_t hcaldist = config_parameters.GetHCALdist();
  Double_t sbsdist = config_parameters.GetSBSdist();
  Double_t bbtheta_rad = config_parameters.GetBBtheta_rad(); //in radians
  std::string sbs_timestamp = config_parameters.GetSBSTimestamp();

  //Add quality plots
    
  TH1D *hE = new TH1D("hE_mc",
		      "HCal E, MC; GeV",
		      bins_SFE,
		      start_SFE,
		      end_SFE);
  
  TH1D *hSF = new TH1D("hSF_mc",
		       "HCal Sampling Fraction; E/KE_{p}",
		       bins_SFE,
		       start_SFE,
		       end_SFE);
  
  TH1D *hW2 = new TH1D("hW2_mc",
		       "W^2; GeV^2",
		       bins_W2,
		       start_W2,
		       end_W2);
  
  TH1D *hdy = new TH1D("hdy_mc",
		       "dy; m",
		       bins_dy,
		       start_dy,
		       end_dy);

  TH1D *hfgl = new TH1D("hfgl_mc",
			"failedglobal bool",
			4,
			0,
			2);



  //Get replayed mc data file
  std::string rootfile_dir = "/lustre19/expphy/volatile/halla/sbs/seeds/sim042123/";
  std::string rootfile_path = rootfile_dir + Form("replayed_gmn_sbs%d*.root",config);

  //Get available cuts for quality check. Use first element (0) of cut
  vector<calcut> cut;
  Int_t Ncal_set_idx = 0; //ignoring additional calibration sets for diagnostic only
  Int_t mag = 0;
  Int_t verb = 0;
  std::string current_target = "lh2";
  Int_t pass = 0;

  //reconfigure cuts to use available data
  if( config==14 || config==8 || config==9 || config==11 )
    pass = 1;
  if( config==7 )
    mag=85;
  if( config==9 )
    mag=70;

  util::ReadCutList(struct_dir,experiment,config,Ncal_set_idx,pass,current_target,mag,verb,cut);
  std::cout << cut[0];

  //Get target configuration
  SBStarget target_parameters(current_target);
  Int_t target_index = target_parameters.GetTargIndex();  //Target index (1:lh2,2:ld2,3:he3)
  Double_t target_length = target_parameters.GetTargLength();
  Double_t target_rho = target_parameters.GetTargRho();
  Double_t cell_rho = target_parameters.GetCellRho();
  Double_t cell_diam = target_parameters.GetCellDiam();
  Double_t cell_dEdx = target_parameters.GetCelldEdx();
  Double_t upstream_wthick = target_parameters.GetUpstreamWallThick();
  Double_t target_dEdx = target_parameters.GetTargdEdx();
  Double_t M_avg = target_parameters.GetAvgMass();

  // Setting up chain and branch addresses
  C = new TChain("T");
  C->Add(rootfile_path.c_str());
  
  C->SetBranchStatus("*",0);    
  Double_t BBtr_p[hcal::maxTracks], BBtr_px[hcal::maxTracks], BBtr_py[hcal::maxTracks], BBtr_pz[hcal::maxTracks];
  Double_t BBtr_vz[hcal::maxTracks];
  Double_t BBtr_n, BBps_x, BBps_y, BBps_e, BBsh_x, BBsh_y, BBsh_e;	
  Double_t HCALx, HCALy, HCALe;
  Double_t pblkrow, pblkcol, nblk, nclus;
  Int_t Ncid;
  Double_t cblkid[hcal::maxHCalChan], cblke[hcal::maxHCalChan], cblkatime[hcal::maxHCalChan], cblktime[hcal::maxHCalChan], cblkagain[hcal::maxHCalChan];
  Double_t cid[hcal::maxHCalClus], crow[hcal::maxHCalClus], ccol[hcal::maxHCalClus], ce[hcal::maxHCalClus], cx[hcal::maxHCalClus], cy[hcal::maxHCalClus], catime[hcal::maxHCalClus];
  
  // Speed up processing by switching on only useful branches
  C->SetBranchStatus( "*", 0 );
  C->SetBranchStatus( "sbs.hcal.x", 1 );
  C->SetBranchStatus( "sbs.hcal.y", 1 );
  C->SetBranchStatus( "sbs.hcal.e", 1 );
  C->SetBranchStatus( "sbs.hcal.rowblk", 1 );
  C->SetBranchStatus( "sbs.hcal.colblk", 1 );
  C->SetBranchStatus( "sbs.hcal.nblk", 1 );
  C->SetBranchStatus( "sbs.hcal.nclus", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.id", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.e", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.tdctime", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.atime", 1 );
  //C->SetBranchStatus( "sbs.hcal.clus_blk.again", 1 );
  C->SetBranchStatus( "sbs.hcal.clus.id", 1 );
  C->SetBranchStatus( "sbs.hcal.clus.row", 1 );
  C->SetBranchStatus( "sbs.hcal.clus.col", 1 );
  C->SetBranchStatus( "sbs.hcal.clus.e", 1 );
  C->SetBranchStatus( "sbs.hcal.clus.x", 1 );
  C->SetBranchStatus( "sbs.hcal.clus.y", 1 );
  C->SetBranchStatus( "sbs.hcal.clus.atime", 1 );
  C->SetBranchStatus( "bb.tr.n", 1 );
  C->SetBranchStatus( "bb.tr.px", 1 );
  C->SetBranchStatus( "bb.tr.py", 1 );
  C->SetBranchStatus( "bb.tr.pz", 1 );    
  C->SetBranchStatus( "bb.tr.vz", 1 );
  C->SetBranchStatus( "bb.tr.p", 1 );
  C->SetBranchStatus( "bb.ps.e", 1 );
  C->SetBranchStatus( "bb.ps.x", 1 );
  C->SetBranchStatus( "bb.ps.y", 1 );
  C->SetBranchStatus( "bb.sh.e", 1 );
  C->SetBranchStatus( "bb.sh.x", 1 );
  C->SetBranchStatus( "bb.sh.y", 1 );
  C->SetBranchStatus( "bb.gem.track.nhits", 1 );
  C->SetBranchStatus( "bb.etot_over_p", 1 );
  C->SetBranchStatus( "Ndata.sbs.hcal.clus.id", 1 ); //Odd maxing out at 10 clusters on all cluster Ndata objects, so this is needed in addition to sbs.hcal.nclus
  
  // Linking memory
  C->SetBranchAddress( "sbs.hcal.x", &HCALx );
  C->SetBranchAddress( "sbs.hcal.y", &HCALy );
  C->SetBranchAddress( "sbs.hcal.e", &HCALe );
  C->SetBranchAddress( "sbs.hcal.rowblk", &pblkrow );
  C->SetBranchAddress( "sbs.hcal.colblk", &pblkcol );
  C->SetBranchAddress( "sbs.hcal.nblk", &nblk ); // Total number of blocks in highest E clus
  C->SetBranchAddress( "sbs.hcal.nclus", &nclus ); // Total number of clusters
  C->SetBranchAddress( "sbs.hcal.clus_blk.id", cblkid ); // kNcell-1 index for each block
  C->SetBranchAddress( "sbs.hcal.clus_blk.e", cblke ); // Array of block energies
  C->SetBranchAddress( "sbs.hcal.clus_blk.tdctime", cblktime ); // Array of block TDC times
  C->SetBranchAddress( "sbs.hcal.clus_blk.atime", cblkatime ); // Array of block ADC times
  //C->SetBranchAddress( "sbs.hcal.clus_blk.again", cblkagain ); // Array of block ADC gain coeff
  C->SetBranchAddress( "sbs.hcal.clus.id", cid );
  C->SetBranchAddress( "sbs.hcal.clus.row", crow );
  C->SetBranchAddress( "sbs.hcal.clus.col", ccol );
  C->SetBranchAddress( "sbs.hcal.clus.e", ce );
  C->SetBranchAddress( "sbs.hcal.clus.x", cx );
  C->SetBranchAddress( "sbs.hcal.clus.y", cy );
  C->SetBranchAddress( "sbs.hcal.clus.atime", catime );
  C->SetBranchAddress( "bb.tr.n", &BBtr_n );
  C->SetBranchAddress( "bb.tr.px", BBtr_px );
  C->SetBranchAddress( "bb.tr.py", BBtr_py );
  C->SetBranchAddress( "bb.tr.pz", BBtr_pz );
  C->SetBranchAddress( "bb.tr.vz", BBtr_vz );
  C->SetBranchAddress( "bb.tr.p", BBtr_p );
  C->SetBranchAddress( "bb.ps.e", &BBps_e );
  C->SetBranchAddress( "bb.ps.x", &BBps_x );
  C->SetBranchAddress( "bb.ps.y", &BBps_y );
  C->SetBranchAddress( "bb.sh.e", &BBsh_e );
  C->SetBranchAddress( "bb.sh.x", &BBsh_x );
  C->SetBranchAddress( "bb.sh.y", &BBsh_y ); 
  C->SetBranchAddress( "Ndata.sbs.hcal.clus.id", &Ncid ); //Odd maxing out at 10 clusters on all cluster Ndata objects, so this is needed in addition to sbs.hcal.nclus
  
  //globalcut enables
  C->SetBranchStatus( "bb.tr.tg_th", 1 );
  C->SetBranchStatus( "bb.tr.tg_ph", 1 );
  
  //Use TTreeFormula to avoid looping over data an additional time
  TCut GCut = cut[0].gcut.c_str();

  TTreeFormula *GlobalCut = new TTreeFormula( "GlobalCut", GCut, C );
    
  // Set up hcal coordinate system with hcal angle wrt exit beamline
  vector<TVector3> hcalaxes; util::sethcalaxes( hcaltheta_rad, hcalaxes );
  TVector3 hcalorigin = hcaldist*hcalaxes[2] + hcal::HCalvoff*hcalaxes[0];
  Double_t BdL = hcal::maxSBSfield * hcal::sbsdipolegap * (mag/100); //scale crudely by magnetic field
  Double_t Eloss_outgoing = cell_diam/2.0/sin(bbtheta_rad) * target_rho * target_dEdx;
  
  long nevent = 0, nevents = C->GetEntries(); 
  Int_t treenum = 0, currenttreenum = 0, elastic_yield = 0;

  //Main loop over events in run
  while (C->GetEntry(nevent++)) {
    
    cout << "Analyzing " <<  nevent << "/" << nevents << " with total elastic yield with data cuts: " << elastic_yield << " \r";
    cout.flush();

    ///////
    //Single-loop elastic globalcut method. Save pass/fail for output tree.
    currenttreenum = C->GetTreeNumber();
    if( nevent == 1 || currenttreenum != treenum ){
      treenum = currenttreenum; 
      GlobalCut->UpdateFormulaLeaves();
    }
    bool failedglobal = GlobalCut->EvalInstance(0) == 0;

    ///////
    //HCal Active Area Cut
    bool failedactivearea = 
      pblkrow==0 || 
      pblkrow==23 || 
      pblkcol==0 || 
      pblkcol==11;
    
    if( failedactivearea ) 
      continue; //All events with primary cluster element on edge blocks cut
    
    ///////
    //HCal primary cluster coincidence time cut (using adctime while hcal tdc suspect, new offsets)
    Int_t pblkid = cblkid[0]-1; //define primary block, primary cluster ID
    
    Double_t natime = cblkatime[0];
    Double_t atime0 = cut[0].atime0; //observed elastic peak in adc time
    Double_t atimesig = cut[0].atime_sig; //observed width of elastic peak in adc time
    
    bool failedcoin = abs(natime-atime0)>atimeNsig*atimesig;
    
    ///////
    //Physics calculations
    //correct beam energy with vertex information, primary track
    Double_t ebeam_c = ebeam - ( (BBtr_vz[0]+target_length/2.0) * target_rho * target_dEdx + upstream_wthick * cell_rho * cell_dEdx );
    
    TVector3 vertex( 0., 0., BBtr_vz[0] );
    
    //reconstructed momentum, corrected for mean energy loss exiting the target (later we'll also want to include Al shielding on scattering chamber)
    Double_t precon = BBtr_p[0] + Eloss_outgoing;
    
    //set up four-momenta with some empty for various calculation methods
    TLorentzVector pbeam( 0., 0., ebeam_c, ebeam_c ); //beam momentum
    TLorentzVector pe( precon*BBtr_px[0]/BBtr_p[0], precon*BBtr_py[0]/BBtr_p[0], precon*BBtr_pz[0]/BBtr_p[0], precon ); //e' recon plvect
    TLorentzVector ptarg; //target momentum
    ptarg.SetPxPyPzE( 0., 0., 0., M_avg );
    TLorentzVector q = pbeam - pe; //virtual photon mom. or mom. transferred to scatter nucleon (N')
    TVector3 qv = q.Vect();
    TLorentzVector pN; //N' momentum
    
    //simple calculations for e' and N'
    Double_t etheta = acos( pe.Pz() / pe.E() );
    Double_t ephi = atan2( pe.Py(), pe.Px() );
    Double_t pcent = ebeam_c/( 1. + ( ebeam_c/M_avg )*( 1.0 - cos(etheta) ) ); //e' p reconstructed by angles
    Double_t phNexp = ephi + hcal::PI;
    Double_t Q2, W2;
    
    //e' p reconstruction with track angles (not momentum)
    Double_t nu = pbeam.E() - pcent;
    Double_t pNexp = sqrt( pow(nu, 2.) + 2. * M_avg * nu );
    Double_t thNexp = acos( (pbeam.E()-pcent*cos(etheta)) / pNexp );
    TVector3 pNhat( sin(thNexp) * cos(phNexp), sin(thNexp) * sin(phNexp), cos(thNexp) );
    pN.SetPxPyPzE( pNexp*pNhat.X(), pNexp*pNhat.Y(), pNexp*pNhat.Z(), nu+ptarg.E() );
    Q2 = 2.0 * pbeam.E() * pe.E() * ( 1.0 - cos(etheta) );
    W2 = pow( M_avg, 2.0 ) + 2.0*M_avg * (ebeam_c-pe.E()) - Q2;
    
    //Calculate h-arm quantities
    vector<Double_t> xyhcalexp; util::getxyhcalexpect( vertex, pNhat, hcalorigin, hcalaxes, xyhcalexp );
    TVector3 hcalpos = hcalorigin + HCALx*hcalaxes[0] + HCALy*hcalaxes[1]; //from primary blk
    Double_t KE_p = nu; //For elastics total predicted by earm
    Double_t SFrac = HCALe/KE_p; //Measured
    Double_t dx = HCALx - xyhcalexp[0];
    Double_t dy = HCALy - xyhcalexp[1];
    
    //Target Energy in HCal for calibrations
    Double_t hcal_samp_frac = cut[0].hcal_sf;
    Double_t hcal_esratio = cut[0].hcal_es;
    
    ///////
    //dy elastic cut
    Double_t dy0 = cut[0].dy0;
    Double_t dysig = cut[0].dy_sig;
    bool faileddy = abs(dy-dy0)>3*dysig;
    
    ////////////////////////////
    //Primary W2 cut on elastics
    Double_t W2_mean = cut[0].W2_mean;
    Double_t W2_sig = cut[0].W2_sig;
    bool failedW2 = fabs(W2-W2_mean)>W2_sig;
    
    //Track elastics with cuts used on data. Will track global cuts seperately.
    if( !faileddy && !failedW2 ) 
      elastic_yield++;

    ///////
    //Get corrected primary cluster energy
    Double_t clusE = 0.0;
    //for( Int_t blk = 0; blk<Ncid; blk++ ){
    for( Int_t blk = 0; blk<nblk; blk++ ){
      Int_t blkid = int(cblkid[blk])-1; //-1 necessary since sbs.hcal.clus_blk.id ranges from 1 to 288 (different than just about everything else)
      Double_t blke = cblke[blk];
      
      clusE += blke;
      
    }
    
    Double_t sampling_fraction = clusE/KE_p;
    
    //Fill histograms
    hE->Fill( clusE );
    hSF->Fill( sampling_fraction );
    hW2->Fill( W2 );
    hdy->Fill( dy );
    hfgl->Fill( failedglobal );

  }//endloop over events
  
  cout << "Ended loop over MC data. " << endl;
  
  fout->Write();
  
  st->Stop();

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;    

}
