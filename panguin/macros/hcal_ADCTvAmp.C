void hcal_ADCTvAmp(){

  cout << "Processing macro.." << endl;

  TStopwatch *st = new TStopwatch();
  st->Start(kTRUE);

  Double_t sbs_hcal_a_amp_p[300] = {0.}, sbs_hcal_a_time[300] = {0.}, sbs_hcal_nclus = 0.; 
  //Double_t sbs_hcal_adcelemID[300] = {0.};
  Int_t Ndata_sbs_hcal_adcelemID = 0;
  Double_t sbs_hcal_clus_e[300] = {0.};
  Double_t fEvtHdr_fTrigBits = -1;

  TH2D *h2_hcal_adctime_v_amp = new TH2D("h2_hcal_adctime_v_amp","HCal ADC Time (ns) vs Amp (mV); ns; mV",80,0,160,40,0,20);

  // Declare branches
  TTree *T = (TTree*) gDirectory->Get("T");
  T->SetBranchStatus("*",0);
  T->SetBranchStatus("sbs.hcal.a_time",1);
  T->SetBranchStatus("sbs.hcal.a_amp_p",1);
  T->SetBranchStatus("sbs.hcal.clus.e",1);
  T->SetBranchStatus("sbs.hcal.nclus",1);
  T->SetBranchStatus("fEvtHdr.fTrigBits",1);
  T->SetBranchStatus("Ndata.sbs.hcal.adcelemID",1);

  T->SetBranchAddress("sbs.hcal.a_time", &sbs_hcal_a_time);
  T->SetBranchAddress("sbs.hcal.a_amp_p", &sbs_hcal_a_amp_p);
  T->SetBranchAddress("sbs.hcal.clus.e", &sbs_hcal_clus_e);
  T->SetBranchAddress("sbs.hcal.nclus", &sbs_hcal_nclus);
  T->SetBranchAddress("fEvtHdr.fTrigBits", &fEvtHdr_fTrigBits);
  T->SetBranchAddress("Ndata.sbs.hcal.adcelemID", &Ndata_sbs_hcal_adcelemID);

  // Acquire the number of entries
  Long64_t nevents = T->GetEntries();
  
  for(Long64_t nevent=0; nevent<nevents; nevent++){

    T->GetEntry(nevent);
    //sbs_hcal_clus_e[0]<0.35
    if( fEvtHdr_fTrigBits==32 || sbs_hcal_nclus==0 || sbs_hcal_clus_e[0]<0.35 ) continue;

    for( int i=0; i<Ndata_sbs_hcal_adcelemID; i++ ){

      if( sbs_hcal_a_time[i]<1 ) continue;
    
      h2_hcal_adctime_v_amp->Fill( sbs_hcal_a_time[i], sbs_hcal_a_amp_p[i] );
    }
  }

  gStyle->SetPalette(53);
  h2_hcal_adctime_v_amp->SetMinimum(-0.1);
  h2_hcal_adctime_v_amp->SetStats(0);
  h2_hcal_adctime_v_amp->Draw("colz");

  cout << "Processed macro with " << nevents << " entries." << endl;
  
  st->Stop();
  cout << "CPU time = " << st->CpuTime() << " s " << " Real time = " << st->RealTime() << " s " << endl;
}
