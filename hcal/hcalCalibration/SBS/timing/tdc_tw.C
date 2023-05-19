//SSeeds 5.19.23 Script to extract timewalk exponential fit parameters from tdc distributions
//NOTE: requires $DB_DIR path set correctly in current environment. Script assumes hardware differences on tdcoffset timestamps and outputs alignment constants for all timestamps within the configuration provided.
//ADDITIONAL NOTE: this script expects that tdc_align.C, adct_align.C, and ecal.C are performed first and will look for new offset parameters to update tdc data set, adct parameters to improve elastic selection, and ecal parameters to update energies

#include <vector>
#include <iostream>
#include <algorithm>
#include "TCut.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TTreeFormula.h"
#include "TChain.h"
#include "TLatex.h"
#include "TH1D.h"
#include "TStopwatch.h"
#include "../include/sbs.h"

// TDC vs E extraction constraints
const Int_t first_hcal_chan = 0;
const Int_t tdc_bins = 500;
const Double_t tdc_lower_lim = -200.;
const Double_t tdc_upper_lim = 50.;
const Int_t E_bins = 1000;
const Double_t E_lower_lim = 0.;
const Double_t E_upper_lim = 500.;
const Int_t fit_event_min = 100;

//Main <experiment> <configuration> <quasi-replay> <replay-pass> <number-of-calibration-sets>; qreplay should only be performed after new offsets obtained
void tdc_tw( const char *experiment = "gmn", Int_t config = 4, bool qreplay = false, Int_t pass = 0, Int_t Nset = 2 ){

  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );
  
  // Get the date
  string date = util::getDate();

  // Set up universal paths and variables
  string db_path = gSystem->Getenv("DB_DIR");
  std::string new_tdcoffset_path = Form("../timing/parameters/tdcoffsets_class_%s_conf%d_pass%d.txt",experiment,config,pass);
  std::string new_adcgain_path = Form("../energy/parameters/adcgaincoeff_%s_conf%d.txt",experiment,config);
  std::string new_adctoffset_path = Form("../timing/parameters/adctoffsets_class_%s_conf%d_pass%d.txt",experiment,config,pass);
  std::string new_tdctw_path = Form("../timing/parameters/tdctw_class_%s_conf%d_pass%d.txt",experiment,config,pass);
  if( qreplay )
    new_tdctw_path = "/dev/null";  //Safety to prevent overwriting constants on quasi-replay
  std::string old_db_path = db_path + "/db_sbs.hcal.dat";
  std::string db_tdcoffset_variable = "sbs.hcal.tdc.offset";
  std::string db_tdccalib_variable = "sbs.hcal.tdc.calib";
  std::string db_tdctwA_variable = "sbs.hcal.tdc.tw_a";
  std::string db_tdctwB_variable = "sbs.hcal.tdc.tw_b";
  std::string db_gain_variable = "sbs.hcal.adc.gain";
  std::string db_adctoffset_variable = "sbs.hcal.adc.timeoffset";

  // Declare output analysis file
  TFile *fout = new TFile( Form("outfiles/tdctw_class_%s_conf%d_qr%d_pass%d.root",experiment,config,(Int_t)qreplay,pass), "RECREATE" );

  // Get information from .csv files
  std::string struct_dir = Form("../config/%s/",experiment); //unique to my environment for now
  Int_t nruns = -1; //Analyze all available runs for this configuration
  Int_t verb = 0; //Don't print diagnostic info by default

  //read the run list to parse run numbers associated to input parameters and set up for loop over runs
  vector<calrun> runs; 
  util::ReadRunList(struct_dir,experiment,nruns,config,pass,verb,runs); //nruns=-1 modifies nruns to loop over all available
  
  //Set up structure to hold calibration parameters and initialize size/index
  calset tdctw_cal[Nset];
  Int_t Ncal_set_size = 0;
  Int_t Ncal_set_idx = 0;
  
  // re-allocate memory at each run to load different cuts/parameters
  TChain *C = nullptr;
  
  // create output tree
  TTree *P = new TTree("P","Analysis Tree");

  // Timing
  Double_t pblkid_out; //hcal cluster, primary block id
  Double_t tdc_out; //hcal cluster tdc time, tree
  Double_t atime_out; //hcal cluster adc time, tree
  Double_t hodotmean_out; //hodoscope cluster mean tdc time

  // Physics
  Double_t dx_out; //hcal cluster dx
  Double_t dy_out; //hcal cluster dy
  Double_t W2_out; //Invariant mass squared W2
  Double_t Q2_out; //Inverse momentum squared Q2
  Double_t hcale_out; //hcal cluster energy
  Double_t ep_out; //track reconstructed e' momentum
  Int_t nblk_out; //total blocks in primary cluster
  Int_t mag_out; //sbs magnetic field strength (percent)
  Int_t run_out; //run number
  Int_t tar_out; //target, LH2 or LD2

  //cluster block variables
  Double_t cblkid_out[hcal::maxHCalBlk];
  Double_t cblkatime_out[hcal::maxHCalBlk];
  Double_t cblktime_out[hcal::maxHCalBlk];
  Double_t cblke_out[hcal::maxHCalBlk];

  // Set output tree branches
  P->Branch( "pblkid", &pblkid_out, "pblkid/D" );
  P->Branch( "tdc", &tdc_out, "tdc/D" );
  P->Branch( "atime", &atime_out, "atime/D" );
  P->Branch( "hodotmean", &hodotmean_out, "hodotmean/D" );
  P->Branch( "dx", &dx_out, "dx/D" );
  P->Branch( "dy", &dy_out, "dy/D" );
  P->Branch( "W2", &W2_out, "W2/D" );
  P->Branch( "Q2", &Q2_out, "Q2/D" );
  P->Branch( "hcale", &hcale_out, "hcale/D" );
  P->Branch( "nblk", &nblk_out, "nblk_out/I" );
  P->Branch( "mag", &mag_out, "mag_out/I" );
  P->Branch( "run", &run_out, "run_out/I" );
  P->Branch( "tar", &tar_out, "tar_out/I" );

  P->Branch( "cblkid", &cblkid_out, Form("cblkid[%d]/D",hcal::maxHCalBlk) );
  P->Branch( "cblkatime", &cblkatime_out, Form("cblkatime[%d]/D",hcal::maxHCalBlk) );
  P->Branch( "cblktime", &cblktime_out, Form("cblktime[%d]/D",hcal::maxHCalBlk) );
  P->Branch( "cblke", &cblke_out, Form("cblke[%d]/D",hcal::maxHCalBlk) );

  //Get experimental configuration parameters
  SBSconfig config_parameters(experiment,config);    
  cout << config_parameters;
  Double_t ebeam = config_parameters.GetEbeam();
  Double_t hcaltheta_rad = config_parameters.GetHCALtheta_rad();
  Double_t hcaldist = config_parameters.GetHCALdist();
  Double_t sbsdist = config_parameters.GetSBSdist();
  Double_t bbtheta_rad = config_parameters.GetBBtheta_rad(); //in radians
  std::string sbs_timestamp = config_parameters.GetSBSTimestamp();

  //Declare histograms for modification later where additional calibration sets are necessary
  TH2D *htdcvE[Nset][hcal::maxHCalChan];
  TH2D *htdcvE_pblk[Nset][hcal::maxHCalChan];

  for( Int_t set=0; set<Nset; set++ ){
    for( Int_t chan=0; chan<hcal::maxHCalChan; chan++ ){
      
      htdcvE[set][chan] = new TH2D(Form("htdcvE_set%d_chan%d",set,chan),
				   "null",
				   E_bins,
				   E_lower_lim,
				   E_upper_lim,
				   tdc_bins,
				   tdc_lower_lim,
				   tdc_upper_lim);
      
      htdcvE_pblk[set][chan] = new TH2D(Form("htdcvE_pblk_set%d_chan%d",set,chan),
				   "null",
				   E_bins,
				   E_lower_lim,
				   E_upper_lim,
				   tdc_bins,
				   tdc_lower_lim,
				   tdc_upper_lim);
      
    }
  }

  //minimize reporting
  Int_t target_change_index;
  Int_t total_analyzed_events;
  Int_t missing_tdc_events;

  //Main loop over runs
  for (Int_t r=0; r<nruns; r++) {

    //Get run experimental parameters
    std::string current_gain_timestamp = runs[r].adcg_ts;
    std::string current_adctoffset_timestamp = runs[r].adcg_ts;
    std::string current_tdcoffset_timestamp = runs[r].tdc_ts;
    std::string current_tdccalib_timestamp = runs[r].tdcc_ts;
    Int_t current_runnumber = runs[r].runnum;

    std::string current_target = runs[r].target;
    std::string targ_uppercase = current_target; transform(targ_uppercase.begin(), targ_uppercase.end(), targ_uppercase.begin(), ::toupper );
    Int_t mag = runs[r].sbsmag / 21; //convert to percent where max field is at 2100A

    //Get run paths
    std::string rootfile_dir = Form("/w/halla-scshelf2102/sbs/sbs-%s/pass%d/SBS%d/%s/rootfiles/",experiment,pass,config,targ_uppercase.c_str());
    std::string rootfile_path = rootfile_dir + Form("*%d*",current_runnumber);

    //Get target configuration
    SBStarget target_parameters(current_target);
    Int_t target_index = target_parameters.GetTargIndex();
    Double_t target_length = target_parameters.GetTargLength();
    Double_t target_rho = target_parameters.GetTargRho();
    Double_t cell_rho = target_parameters.GetCellRho();
    Double_t cell_diam = target_parameters.GetCellDiam();
    Double_t cell_dEdx = target_parameters.GetCelldEdx();
    Double_t upstream_wthick = target_parameters.GetUpstreamWallThick();
    Double_t target_dEdx = target_parameters.GetTargdEdx();
    Double_t M_avg = target_parameters.GetAvgMass();
    
    //send to console if the target changes
    if( target_change_index != target_index ){
      cout << target_parameters;
      target_change_index = target_index;
    }

    ////Timestamp juggling. There is probably a more efficient and clear way to do this.

    //Record new offset/coefficient parameters. Repopulate array on each run for changing timestamps
    Double_t new_tdc_offsets[hcal::maxHCalChan] = {0.};
    Double_t new_adct_offsets[hcal::maxHCalChan] = {0.};
    Double_t new_adcg_coeff[hcal::maxHCalChan] = {1.};
    Double_t new_tdctwA_par[hcal::maxHCalChan] = {0.};
    Double_t new_tdctwB_par[hcal::maxHCalChan] = {0.};

    //Get timestamp for new tdc offsets for this run
    std::string tdc_active_timestamp = "-------[ 0000-00-00 00:00:00 ]";
    util::tsCompare( current_tdcoffset_timestamp, config_parameters.GetSBSTimestamp(), tdc_active_timestamp );
    util::tsCompare( tdc_active_timestamp, current_tdccalib_timestamp, tdc_active_timestamp );
    util::readDB( new_tdcoffset_path, tdc_active_timestamp, db_tdcoffset_variable, new_tdc_offsets );
    
    cout << "1" << endl;

    //get timestamp for new adct offsets for this run
    std::string adct_active_timestamp = "-------[ 0000-00-00 00:00:00 ]";
    util::tsCompare( current_adctoffset_timestamp, config_parameters.GetSBSTimestamp(), adct_active_timestamp );
    util::readDB( new_adctoffset_path, adct_active_timestamp, db_adctoffset_variable, new_adct_offsets );

    cout << "2" << endl;


    //get timestamp for new adc gain coefficients for this run
    std::string active_timestamp = "-------[ 0000-00-00 00:00:00 ]";
    util::tsCompare( current_gain_timestamp, config_parameters.GetSBSTimestamp(), active_timestamp );
    util::readDB( new_adcgain_path, active_timestamp, db_gain_variable, new_adcg_coeff );
    //quasi-replay, get new timewalk parameters
    if( qreplay ){
      util::readDB( new_tdctw_path, active_timestamp, db_tdctwA_variable, new_tdctwA_par );
      util::readDB( new_tdctw_path, active_timestamp, db_tdctwB_variable, new_tdctwB_par );
    }

    cout << "3" << endl;


    //////////
    //Select correct calibration set. Timewalk corrections will go with energy, active_timestamp left unaltered.
    bool new_ts = true;
    bool new_offset_ts = true;
    //Look through available sets for duplicate timestamps
    for( Int_t cs=0; cs<Ncal_set_size; cs++ ){ 
      if( tdctw_cal[cs].timestamp == active_timestamp ){
	new_ts = false;
	Ncal_set_idx = cs;
      }
    }
    
    //if a new timestamp exists, expand the set and update index/old-offsets/old-calib
    if( new_ts || Ncal_set_size==0){

      cout << "Adding new calibration set. Active (energy calibration) timestamp:" << active_timestamp << endl << endl;

      Ncal_set_idx = Ncal_set_size;

      if( Ncal_set_size == hcal::gNstamp-1 ){
	cout << "Error: Allowable calibration sets exceed maximum of " << hcal::gNstamp << ". Check timestamps and adjust limits as necessary." << endl;
	return;
      }

      //update timestamp and old coefficients for corrections. Don't use active timestamp, use db timestamps for old
      tdctw_cal[Ncal_set_size].timestamp = active_timestamp.c_str();
      util::readDB( old_db_path, current_gain_timestamp, db_gain_variable, tdctw_cal[Ncal_set_size].old_param ); //gain 
      util::readDB( old_db_path, current_tdcoffset_timestamp, db_tdcoffset_variable, tdctw_cal[Ncal_set_size].old_paramB );  //tdcoffset
      util::readDB( old_db_path, current_adctoffset_timestamp, db_adctoffset_variable, tdctw_cal[Ncal_set_size].old_paramC );  //adctoffset   
      util::readDB( old_db_path, current_tdccalib_timestamp, db_tdccalib_variable, tdctw_cal[Ncal_set_size].tdc_calib );  //tdcoffset


      //build analysis histograms
      for( Int_t chan=0; chan<hcal::maxHCalChan; chan++ ){
	htdcvE[Ncal_set_size][chan]->SetTitle(Form("TDC vs E, set: %d, chan: %d;MeV;ns",Ncal_set_size,chan));
	htdcvE_pblk[Ncal_set_size][chan]->SetTitle(Form("TDC vs E primary block only, set: %d, chan: %d;MeV;ns",Ncal_set_size,chan));
      }

      Ncal_set_size++;
    }     

    //Get available cuts for current config/target/field combination. Use first element (0) of cut
    vector<calcut> cut;
    util::ReadCutList(struct_dir,experiment,config,pass,current_target,mag,verb,cut);

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
    Double_t HODOtmean;

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
    C->SetBranchStatus( "bb.hodotdc.clus.tmean", 1 );
    C->SetBranchStatus( "bb.hodotdc.clus.tmean", 1 );
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
    C->SetBranchAddress( "bb.hodotdc.clus.tmean", &HODOtmean );
    C->SetBranchAddress( "Ndata.sbs.hcal.clus.id", &Ncid ); //Odd maxing out at 10 clusters on all cluster Ndata objects, so this is needed in addition to sbs.hcal.nclus

    //Use TTreeFormula to avoid looping over data an additional time
    TCut GCut = cut[0].gcut.c_str();
    TTreeFormula *GlobalCut = new TTreeFormula( "GlobalCut", GCut, C );

    // Set up hcal coordinate system with hcal angle wrt exit beamline
    vector<TVector3> hcalaxes; util::sethcalaxes( hcaltheta_rad, hcalaxes );
    TVector3 hcalorigin = hcaldist*hcalaxes[2] + hcal::HCalvoff*hcalaxes[0];
    Double_t BdL = hcal::maxSBSfield * hcal::sbsdipolegap * (mag/100); //scale crudely by magnetic field
    Double_t Eloss_outgoing = cell_diam/2.0/sin(bbtheta_rad) * target_rho * target_dEdx;

    long nevent = 0, nevents = C->GetEntries(); 
    Int_t treenum = 0, currenttreenum = 0;

    //Main loop over events in run
    while (C->GetEntry(nevent++)) {

      cout << "Analyzing run " << current_runnumber << ": " <<  nevent << "/" << nevents << " \r";
      cout.flush();

      ///////
      //Single-loop elastic globalcut method. Save pass/fail for output tree.
      currenttreenum = C->GetTreeNumber();
      if( nevent == 1 || currenttreenum != treenum ){
	treenum = currenttreenum; 
	GlobalCut->UpdateFormulaLeaves();
      }
      bool failedglobal = GlobalCut->EvalInstance(0) == 0;
	  
      if( failedglobal ) continue;

      ///////
      //HCal coincidence time cut (using adctime while hcal tdc suspect, new offsets)
      Int_t pblkid = cblkid[0]-1; //define primary block, primary cluster ID
      Double_t hcaladct = cblkatime[0]+tdctw_cal[Ncal_set_idx].old_paramC[pblkid]-new_adct_offsets[pblkid]; //new atime
      Double_t atime0 = cut[0].atime0; //observed elastic peak in adc time
      Double_t atimesig = cut[0].atime_sig; //observed width of elastic peak in adc time

      bool failedcoin = abs(hcaladct-atime0)>3*atimesig;
	  
      if( failedcoin ) continue; //All events where adctime outside of reasonable window cut

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
      Double_t dx = HCALx - xyhcalexp[0];
      Double_t dy = HCALy - xyhcalexp[1];

      ///////
      //BigBite/SBS Acceptance matching.
      bool failedaccmatch = 
	xyhcalexp[1] > hcal::posHCalYf ||
	xyhcalexp[1] < hcal::posHCalYi ||
	xyhcalexp[0] > hcal::posHCalXf ||
	xyhcalexp[0] < hcal::posHCalXi;
	  
      if( failedaccmatch ) continue;

      //////
      //Cut on bad tdc and keep track of missing primary cluster events
      total_analyzed_events++;
      bool missedtdc = cblktime[0]==-1000 ||
		       cblktime[0]>10000;
      if( missedtdc ){
	missing_tdc_events++;
	continue;
      }
	
      ///////
      //Address timewalk with primary cluster block tdc time vs energy
      Double_t pblke = 0;
      Double_t tdc_tc = -1000;
      for( Int_t blk = 0; blk<nblk; blk++ ){
	
	//get the channel id for the blk-th block in the cluster
	Int_t blkid = (Int_t)cblkid[blk]-1;

	//get updated primary cluster block tdc time
	Double_t blktime = cblktime[blk] + 
	  tdctw_cal[Ncal_set_idx].old_param[pblkid] * tdctw_cal[Ncal_set_idx].tdc_calib - 
	  new_tdc_offsets[pblkid] * tdctw_cal[Ncal_set_idx].tdc_calib;
	  
	//add trigger correction from hodo
	Double_t blktime_tc = blktime - HODOtmean;

	//get updated primary cluster block energy
	Double_t blke = cblke[blk]/tdctw_cal[Ncal_set_idx].old_param[blkid]*new_adcg_coeff[blkid];

	//on quasi-replay, apply timewalk correction by block as well
	if( qreplay )
	  blktime_tc -= new_tdctwA_par[blkid]*exp(-new_tdctwB_par[blkid]*blke);

	//get primary cluster primary block information
	if( blk==0 ){
	  pblke = blke;
	  tdc_tc = blktime_tc;
	}

	//build tdc vs E over all cluster blocks
	htdcvE[Ncal_set_idx][blk]->Fill( blke, blktime_tc );

      }

      //build tdc vs E for primary block only
      htdcvE_pblk[Ncal_set_idx][pblkid]->Fill( pblke, tdc_tc );
      
      //fill analysis tree
      pblkid_out = pblkid;
      tdc_out = cblktime[0];
      atime_out = cblkatime[0];
      hcale_out = HCALe;
      dx_out = dx;
      dy_out = dy;
      W2_out = W2;
      Q2_out = Q2;
      nblk_out = nblk;
      mag_out = mag;
      run_out = current_runnumber;
      tar_out = target_index;
      hodotmean_out = HODOtmean;
      
      P->Fill();
      
    }//end loop over event
    
    // getting ready for the next run
    C->Reset();
    
  }//endloop over runs
  
  // Declare canvas and graph for plots
  TCanvas *TDCtw_top[Ncal_set_size];
  TCanvas *TDCtw_bot[Ncal_set_size];

  // set up arrays for timewalk fits
  Double_t TDCvseP0[Ncal_set_size][hcal::maxHCalChan];
  Double_t TDCvseP1[Ncal_set_size][hcal::maxHCalChan];
  Double_t TDCvseP2[Ncal_set_size][hcal::maxHCalChan];

  //Loop over all independent data sets
  for( Int_t s=0; s<Ncal_set_size; s++ ){

    TDCtw_top[s] = new TCanvas(Form("TDC_top_s%d",s),Form("TDC_top_s%d",s),1600,1200);
    TDCtw_bot[s] = new TCanvas(Form("TDC_bot_s%d",s),Form("TDC_bot_s%d",s),1600,1200);

    //Fits for timewalk corrections
    for(Int_t c=0; c<hcal::maxHCalChan; c++){
      
      //initialize the arrays
      TDCvseP0[s][c] = 0.;
      TDCvseP1[s][c] = 0.;
      TDCvseP2[s][c] = 0.;

      //Index through the canvas
      TDCtw_top[s]->cd(c+1);
      if( c>=144 ){
	TDCtw_bot[s]->cd(c-143);
	gStyle->SetOptStat(0);
      }

      //Fit the TDC vs E plots
      TF1 *fittdcTW = new TF1( "fittdcTW", util::g_twfit, 0, 300, 3 );
      fittdcTW->SetParameters(14,0.04,-77);
      fittdcTW->SetParLimits(0,2,26);
      fittdcTW->SetParLimits(1,0.01,0.2);
      fittdcTW->SetParLimits(2,-200,50);

      if( htdcvE_pblk[s][c]->GetEntries()>fit_event_min ){
	htdcvE_pblk[s][c]->Fit("fittdcTW","Q","",5,300);
	TDCvseP0[s][c] = fittdcTW->GetParameter(0);
	TDCvseP1[s][c] = fittdcTW->GetParameter(1);
	TDCvseP2[s][c] = fittdcTW->GetParameter(2);
	htdcvE_pblk[s][c]->SetTitle(Form("P0:%f P1:%f P2:%f",TDCvseP0[s][c],TDCvseP1[s][c],TDCvseP2[s][c]));
	htdcvE_pblk[s][c]->Draw();
      }
    }

    TDCtw_top[s]->Write();
    TDCtw_bot[s]->Write();

  }

  // Create timewalk fit parameter files
  ofstream tdctw;
  tdctw.open( new_tdcoffset_path );

  if( !qreplay ){
    for( Int_t p=0; p<2; p++ ){
      if( p==0 )
	tdctw.open( new_tdctw_path );
      else
	tdctw.open( new_tdctw_path, fstream::app );
      
      tdctw << endl << endl << "#HCal tdc vs E fit parameter P" << p << " obtained " << date.c_str() << endl;
      tdctw << "#Exponential fit for i PMTs -> y = P0_i*exp(-P1_i*x)+P2_i" << endl;

      if( p==0 )
	tdctw << db_tdctwA_variable << " =" << endl;
      else
	tdctw << db_tdctwB_variable << " =" << endl;

      Int_t cell = 0;
      for( Int_t r=0; r<hcal::maxHCalRows; r++ ){
	for( Int_t c=0; c<hcal::maxHCalCols; c++ ){
	  if( p==0 )
	    tdctw << TDCvseP0[cell] << "  ";
	  else
	    tdctw << TDCvseP1[cell] << "  ";
	  cell++;
	}
	tdctw << endl;
      }
    }
  }

  tdctw.close();

  //clean up
  for( Int_t unused_set_index=Ncal_set_size; unused_set_index<Nset; unused_set_index++ ){
    for( Int_t c=0; c<hcal::maxHCalChan; c++ ){
      delete htdcvE[unused_set_index][c];
      delete htdcvE_pblk[unused_set_index][c];
    }
  }

  fout->Write();

  st->Stop();

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;    

}

