void hcal_heatmap(){

  cout << "Processing macro.." << endl;

  TStopwatch *st = new TStopwatch();
  st->Start(kTRUE);

  Int_t Ndata_bb_tdctrig_tdcelemID = 0;
  Double_t sbs_hcal_colblk = 0., sbs_hcal_rowblk = 0.; 
  Double_t bb_sh_nclus = 0., sbs_hcal_nclus = 0.;
  Double_t e_kine_W2 = 0.;
  Double_t bb_tdctrig_tdc[6] = {0.}, bb_tdctrig_tdcelemID[6] = {0.};

  TH2F *h2_hcal_heatmap = new TH2F("h2_hcal_heatmap","HCal Elastic Coincidence Events",12,0,12,24,0,24);

  // Declare branches
  TTree *T = (TTree*) gDirectory->Get("T");
  T->SetBranchStatus("*",0);
  T->SetBranchStatus("bb.tdctrig.tdc",1);
  T->SetBranchStatus("bb.sh.nclus",1);
  T->SetBranchStatus("sbs.hcal.colblk",1);
  T->SetBranchStatus("sbs.hcal.nclus",1);
  T->SetBranchStatus("sbs.hcal.rowblk",1);
  T->SetBranchStatus("Ndata.bb.tdctrig.tdcelemID",1);
  T->SetBranchStatus("bb.tdctrig.tdcelemID",1);
  T->SetBranchStatus("e.kine.W2",1);
  
  T->SetBranchAddress("bb.tdctrig.tdc", &bb_tdctrig_tdc);
  T->SetBranchAddress("bb.tdctrig.tdcelemID", &bb_tdctrig_tdcelemID);
  T->SetBranchAddress("Ndata.bb.tdctrig.tdcelemID", &Ndata_bb_tdctrig_tdcelemID);
  T->SetBranchAddress("e.kine.W2", &e_kine_W2);
  T->SetBranchAddress("sbs.hcal.nclus", &sbs_hcal_nclus);
  T->SetBranchAddress("bb.sh.nclus", &bb_sh_nclus);
  T->SetBranchAddress("sbs.hcal.colblk", &sbs_hcal_colblk);
  T->SetBranchAddress("sbs.hcal.rowblk", &sbs_hcal_rowblk);
  
  // Acquire the number of entries
  Long64_t nevents = T->GetEntries();
  
  for(Long64_t nevent=0; nevent<nevents; nevent++){

    T->GetEntry(nevent);

    if(bb_sh_nclus==0 || sbs_hcal_nclus==0) continue;

    Double_t bbcal_time=0., hcal_time=0.;
    for(Int_t ihit=0; ihit<Ndata_bb_tdctrig_tdcelemID; ihit++){
      if(bb_tdctrig_tdcelemID[ihit]==5) bbcal_time=bb_tdctrig_tdc[ihit];
      if(bb_tdctrig_tdcelemID[ihit]==0) hcal_time=bb_tdctrig_tdc[ihit];

      //cout << hcal_time - bbcal_time << endl;

    }
    Double_t diff = hcal_time - bbcal_time; 
    if(fabs(diff-510.)<20.&&e_kine_W2<1.2&&e_kine_W2>0.5){
      h2_hcal_heatmap->Fill( sbs_hcal_colblk+1, sbs_hcal_rowblk+1);

      //cout << "Yes" << endl;

    }
  }
  
  gStyle->SetPalette(53);
  h2_hcal_heatmap->SetMinimum(-0.1);
  h2_hcal_heatmap->SetStats(0);
  h2_hcal_heatmap->Draw("colz");

  cout << "Processed macro with " << nevents << " entries." << endl;

  st->Stop();
  cout << "CPU time = " << st->CpuTime() << " s " << " Real time = " << st->RealTime() << " s " << endl;
}
