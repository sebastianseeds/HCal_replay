void hcal_ADCTime_LED(){

  cout << "Processing macro.." << endl;

  TStopwatch *st = new TStopwatch();
  st->Start(kTRUE);

  Double_t sbs_hcal_clus_e[30] = {0.}, sbs_hcal_clus_atime[30] = {0.}, sbs_hcal_nclus = 0.; 
  UInt_t fEvtHdr_fTrigBits = -1;

  Double_t sbs_hcal_ledbit = 0.;

  TH1D *h1_hcal_adctime_led = new TH1D("h1_hcal_adctime_led","HCal ADC Time; HCal ADC time",80,0,160);

  // Declare branches
  TTree *T = (TTree*) gDirectory->Get("T");
  T->SetBranchStatus("*",0);
  T->SetBranchStatus("sbs.hcal.clus.atime",1);
  T->SetBranchStatus("sbs.hcal.clus.e",1);
  T->SetBranchStatus("sbs.hcal.nclus",1);
  T->SetBranchStatus("fEvtHdr.fTrigBits",1);
  T->SetBranchStatus("sbs.hcal.ledbit",1);

  T->SetBranchAddress("sbs.hcal.clus.atime", &sbs_hcal_clus_atime);
  T->SetBranchAddress("sbs.hcal.clus.e", &sbs_hcal_clus_e);
  T->SetBranchAddress("sbs.hcal.nclus", &sbs_hcal_nclus);
  T->SetBranchAddress("fEvtHdr.fTrigBits", &fEvtHdr_fTrigBits);
  T->SetBranchAddress("sbs.hcal.ledbit", &sbs_hcal_ledbit);

  // Acquire the number of entries
  Long64_t nevents = T->GetEntries();
  
  for(Long64_t nevent=0; nevent<nevents; nevent++){

    T->GetEntry(nevent);
    //sbs_hcal_clus_e[0]<0.35
    if( sbs_hcal_clus_atime[0]<1 || sbs_hcal_nclus==0 ) continue;

    if(fEvtHdr_fTrigBits==32) {
    //if(sbs_hcal_ledbit!=-1){
      h1_hcal_adctime_led->Fill( sbs_hcal_clus_atime[0] );
    }
  }
  
  h1_hcal_adctime_led->SetMinimum(-0.1);
  h1_hcal_adctime_led->SetStats(0);
  h1_hcal_adctime_led->Draw();

  cout << "Processed macro with " << nevents << " entries." << endl;
  
  st->Stop();
  cout << "CPU time = " << st->CpuTime() << " s " << " Real time = " << st->RealTime() << " s " << endl;
}
