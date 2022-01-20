void hcal_timing_expert(){

  cout << "Processing macro.." << endl;

  //TStopwatch *st = new TStopwatch();
  //st->Start(kTRUE);

  Double_t sbs_hcal_tdc = 0., sbs_hcal_a_time = 0.; 

  TH2D *h2_hcal_tdc_vs_adc = new TH2D("h2_hcal_tdc_vs_adc","HCal Timing Check; HCal TDC; HCal ADC time",300,-300,300,80,0,160);

  // Declare branches
  TTree *T = (TTree*) gDirectory->Get("T");
  T->SetBranchStatus("*",0);
  T->SetBranchStatus("sbs.hcal.a_time",1);
  T->SetBranchStatus("sbs.hcal.tdc",1);

  T->SetBranchAddress("sbs.hcal.a_time", &sbs_hcal_a_time);
  T->SetBranchAddress("sbs.hcal.tdc", &sbs_hcal_tdc);

  // Acquire the number of entries
  Long64_t nevents = T->GetEntries();
  
  for(Long64_t nevent=0; nevent<nevents; nevent++){

    T->GetEntry(nevent);
    
    if( sbs_hcal_a_time<1 || sbs_hcal_tdc>1000 ) continue;
    
    h2_hcal_tdc_vs_adc->Fill( sbs_hcal_tdc, sbs_hcal_a_time );
    
    /*
    Double_t bbcal_time=0., hcal_time=0.;
    for(Int_t ihit=0; ihit<Ndata_bb_tdctrig_tdcelemID; ihit++){
      if( bb_tdctrig_tdcelemID[ihit]==5 ) bbcal_time=bb_tdctrig_tdc[ ihit ];
      if( bb_tdctrig_tdcelemID[ihit]==0 ) hcal_time=bb_tdctrig_tdc[ ihit ];

      //cout << hcal_time - bbcal_time << endl;

    }
    Double_t diff = hcal_time - bbcal_time; 
    if( fabs( diff-510. )<20.&&e_kine_W2<1.2&&e_kine_W2>0.5 ){
      h2_hcal_heatmap->Fill( sbs_hcal_colblk+1, sbs_hcal_rowblk+1 );

      //cout << "Yes" << endl; 
      
      }
  */
  }
  
  gStyle->SetPalette(53);
  h2_hcal_tdc_vs_adc->SetMinimum(-0.1);
  h2_hcal_tdc_vs_adc->SetStats(0);
  h2_hcal_tdc_vs_adc->Draw("colz");

  cout << "Processed macro with " << nevents << " entries." << endl;
  
  //st->Stop();
  //cout << "CPU time = " << st->CpuTime() << " s " << " Real time = " << st->RealTime() << " s " << endl;
}
