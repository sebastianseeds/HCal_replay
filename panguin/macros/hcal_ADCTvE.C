void hcal_ADCTvE(){

  cout << "Processing macro.." << endl;

  TStopwatch *st = new TStopwatch();
  st->Start(kTRUE);

  Double_t sbs_hcal_clus_e[30] = {0.}, sbs_hcal_clus_atime[30] = {0.}, sbs_hcal_nclus = 0., sbs_hcal_clus_a_amp_p[30]={0.}; 
  Double_t fEvtHdr_fTrigBits = -1;

  TH2D *h2_hcal_adctime_v_e = new TH2D("h2_hcal_adctime_v_e","HCal ADC Time vs Clus E; ns; GeV",80,0,160,00,1,6);

  // Declare branches
  TTree *T = (TTree*) gDirectory->Get("T");
  T->SetBranchStatus("*",0);
  T->SetBranchStatus("sbs.hcal.clus.atime",1);
  //T->SetBranchStatus("sbs.hcal.clus.a_amp_p",1);
  T->SetBranchStatus("sbs.hcal.clus.e",1);
  T->SetBranchStatus("sbs.hcal.nclus",1);
  T->SetBranchStatus("fEvtHdr.fTrigBits",1);

  T->SetBranchAddress("sbs.hcal.clus.atime", &sbs_hcal_clus_atime);
  //T->SetBranchAddress("sbs.hcal.clus.a_amp_p", &sbs_hcal_clus_a_amp_p);
  T->SetBranchAddress("sbs.hcal.clus.e", &sbs_hcal_clus_e);
  T->SetBranchAddress("sbs.hcal.nclus", &sbs_hcal_nclus);
  T->SetBranchAddress("fEvtHdr.fTrigBits", &fEvtHdr_fTrigBits);

  // Acquire the number of entries
  Long64_t nevents = T->GetEntries();
  
  for(Long64_t nevent=0; nevent<nevents; nevent++){

    T->GetEntry(nevent);
    //sbs_hcal_clus_e[0]<0.35
    if( sbs_hcal_clus_atime[0]<1 || fEvtHdr_fTrigBits==32 || sbs_hcal_nclus==0 ) continue;
    
    h2_hcal_adctime_v_e->Fill( sbs_hcal_clus_atime[0], sbs_hcal_clus_e[0] );

  }

  gStyle->SetPalette(53);
  h2_hcal_adctime_v_e->SetMinimum(-0.1);
  h2_hcal_adctime_v_e->SetStats(0);
  h2_hcal_adctime_v_e->Draw("colz");

  cout << "Processed macro with " << nevents << " entries." << endl;
  
  st->Stop();
  cout << "CPU time = " << st->CpuTime() << " s " << " Real time = " << st->RealTime() << " s " << endl;
}
