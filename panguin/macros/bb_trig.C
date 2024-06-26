void bb_trig(){

  cout << "Processing macro.." << endl;

  TStopwatch *st = new TStopwatch();
  st->Start(kTRUE);

  Int_t Ndata_bb_tdctrig_tdcelemID = 0;
  Double_t bb_sh_e = 0., sbs_hcal_rowblk = 0.; 
  Double_t bb_sh_nclus = 0., sbs_hcal_nclus = 0.;
  //Double_t bb_sh_e = 0.;
  Double_t bb_tdctrig_tdc[6] = {0.}, bb_tdctrig_tdcelemID[6] = {0.};

  TH1D *h1_bb_trig = new TH1D("h1_bb_trig","BBCal trigger; BBCal Trigger", 140, 260, 330);

  // Declare branches
  TTree *T = (TTree*) gDirectory->Get("T");
  T->SetBranchStatus("*",0);
  T->SetBranchStatus("bb.tdctrig.tdc",1);
  T->SetBranchStatus("bb.sh.nclus",1);
  T->SetBranchStatus("sbs.hcal.nclus",1);
  T->SetBranchStatus("bb.sh.e",1);
  T->SetBranchStatus("Ndata.bb.tdctrig.tdcelemID",1);
  T->SetBranchStatus("bb.tdctrig.tdcelemID",1);

  T->SetBranchAddress("bb.tdctrig.tdc", &bb_tdctrig_tdc);
  T->SetBranchAddress("bb.tdctrig.tdcelemID", &bb_tdctrig_tdcelemID);
  T->SetBranchAddress("Ndata.bb.tdctrig.tdcelemID", &Ndata_bb_tdctrig_tdcelemID);
  T->SetBranchAddress("bb.sh.nclus", &bb_sh_nclus);
  T->SetBranchAddress("sbs.hcal.nclus", &sbs_hcal_nclus);
  T->SetBranchAddress("bb.sh.e", &bb_sh_e);

  // Loop through events
  Long64_t nevents = T->GetEntries();
  for(Long64_t nevent=0; nevent<nevents; nevent++){

    T->GetEntry(nevent);

    if( bb_sh_nclus==0 || bb_sh_e<0.5 || sbs_hcal_nclus==0 ) continue;

    Double_t bb_time=0.;
    for(Int_t ihit=0; ihit<Ndata_bb_tdctrig_tdcelemID; ihit++){
      if(bb_tdctrig_tdcelemID[ihit]==5) bb_time=bb_tdctrig_tdc[ihit];
    }
    h1_bb_trig->Fill(bb_time);
  }
  
  h1_bb_trig->SetStats(0);
  h1_bb_trig->Draw();

  cout << "Processed macro with " << nevents << " entries." << endl;

  st->Stop();
  cout << "CPU time = " << st->CpuTime() << " s " << " Real time = " << st->RealTime() << " s " << endl;
}
