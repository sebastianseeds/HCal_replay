void bbcal_hcal_clust_corr(){

  cout << "Processing macro.." << endl;

  TStopwatch *st = new TStopwatch();
  st->Start(kTRUE);

  Int_t Ndata_bb_tdctrig_tdcelemID = 0;
  Double_t bb_sh_rowblk = 0., sbs_hcal_rowblk = 0.; 
  Double_t bb_sh_nclus = 0., bb_sh_e = 0., sbs_hcal_nclus = 0.;
  Double_t bb_tdctrig_tdc[6] = {0.}, bb_tdctrig_tdcelemID[6] = {0.};

  TH2F *h2_bbcal_hcal_corr = new TH2F("h2_bbh_corr","BBCal-HCal cluster correlation; BB Shower Rows; HCal Rows",27,1,28,24,1,25);

  // Declare branches
  TTree *T = (TTree*) gDirectory->Get("T");
  T->SetBranchStatus("*",0);
  T->SetBranchStatus("bb.tdctrig.tdc",1);
  T->SetBranchStatus("bb.sh.nclus",1);
  T->SetBranchStatus("bb.sh.e",1);
  T->SetBranchStatus("bb.sh.rowblk",1);
  T->SetBranchStatus("sbs.hcal.nclus",1);
  T->SetBranchStatus("sbs.hcal.rowblk",1);
  T->SetBranchStatus("Ndata.bb.tdctrig.tdcelemID",1);
  T->SetBranchStatus("bb.tdctrig.tdcelemID",1);
  T->SetBranchAddress("bb.tdctrig.tdc", &bb_tdctrig_tdc);
  T->SetBranchAddress("bb.tdctrig.tdcelemID", &bb_tdctrig_tdcelemID);
  T->SetBranchAddress("Ndata.bb.tdctrig.tdcelemID", &Ndata_bb_tdctrig_tdcelemID);
  T->SetBranchAddress("bb.sh.nclus", &bb_sh_nclus);
  T->SetBranchAddress("bb.sh.e", &bb_sh_e);
  T->SetBranchAddress("sbs.hcal.nclus", &sbs_hcal_nclus);
  T->SetBranchAddress("bb.sh.rowblk", &bb_sh_rowblk);
  T->SetBranchAddress("sbs.hcal.rowblk", &sbs_hcal_rowblk);
  
  // Loop through events
  Long64_t nevents = T->GetEntries();
  for(Long64_t nevent=0; nevent<nevents; nevent++){

    T->GetEntry(nevent);

    if(bb_sh_nclus==0 || sbs_hcal_nclus==0 || bb_sh_e<0.5) continue;

    Double_t bbcal_time=0., hcal_time=0.;
    for(Int_t ihit=0; ihit<Ndata_bb_tdctrig_tdcelemID; ihit++){
      if(bb_tdctrig_tdcelemID[ihit]==5) bbcal_time=bb_tdctrig_tdc[ihit];
      if(bb_tdctrig_tdcelemID[ihit]==0) hcal_time=bb_tdctrig_tdc[ihit];
    }
    Double_t diff = hcal_time - bbcal_time; 
    if(fabs(diff-510.)<20.){
      h2_bbcal_hcal_corr->Fill( bb_sh_rowblk+1, sbs_hcal_rowblk+1);
    }
  }
  
  gStyle->SetPalette(53);
  h2_bbcal_hcal_corr->SetMinimum(-0.1);
  h2_bbcal_hcal_corr->SetStats(0);
  h2_bbcal_hcal_corr->Draw("colz");

  cout << "Processed macro with " << nevents << " entries." << endl;

  st->Stop();
  cout << "CPU time = " << st->CpuTime() << " s " << " Real time = " << st->RealTime() << " s " << endl;
}
