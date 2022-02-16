void hcal_ADCTime_NoLED(){

  cout << "Processing macro.." << endl;

  TStopwatch *st = new TStopwatch();
  st->Start(kTRUE);

  Double_t sbs_hcal_clus_e[30] = {0.}, sbs_hcal_clus_atime[30] = {0.}, sbs_hcal_nclus = 0., sbs_hcal_ledbit = -2.; 
  UInt_t fEvtHdr_fTrigBits = -1;

  TH1D *h1_hcal_adctime_noled = new TH1D("h1_hcal_adctime_noled","HCal ADC Time; HCal ADC time",80,0,160);

  // Declare branches
  TTree *T = (TTree*) gDirectory->Get("T");
  T->SetBranchStatus("*",0);
  T->SetBranchStatus("sbs.hcal.clus.atime",1);
  T->SetBranchStatus("sbs.hcal.clus.e",1);
  T->SetBranchStatus("sbs.hcal.nclus",1);
  T->SetBranchStatus("sbs.hcal.ledbit",1);
  T->SetBranchStatus("fEvtHdr.fTrigBits",1);

  T->SetBranchAddress("sbs.hcal.clus.atime", &sbs_hcal_clus_atime);
  T->SetBranchAddress("sbs.hcal.clus.e", &sbs_hcal_clus_e);
  T->SetBranchAddress("sbs.hcal.nclus", &sbs_hcal_nclus);
  T->SetBranchAddress("sbs.hcal.ledbit", &sbs_hcal_ledbit);
  T->SetBranchAddress("fEvtHdr.fTrigBits", &fEvtHdr_fTrigBits);

  // Acquire the number of entries
  Long64_t nevents = T->GetEntries();
  
  for(Long64_t nevent=0; nevent<nevents; nevent++){

    T->GetEntry(nevent);
    if( sbs_hcal_clus_atime[0]<1 || fEvtHdr_fTrigBits==32 || sbs_hcal_nclus==0 ) continue;
    
    h1_hcal_adctime_noled->Fill( sbs_hcal_clus_atime[0] );

  }
  
  h1_hcal_adctime_noled->SetMinimum(-0.1);
  h1_hcal_adctime_noled->SetStats(0);
  h1_hcal_adctime_noled->Draw();

  cout << "Processed macro with " << nevents << " entries." << endl;
  
  st->Stop();
  cout << "CPU time = " << st->CpuTime() << " s " << " Real time = " << st->RealTime() << " s " << endl;
}
