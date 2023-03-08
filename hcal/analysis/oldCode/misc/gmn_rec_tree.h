//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Oct 17 01:45:40 2021 by ROOT version 6.24/06
// from TTree T/Hall A Analyzer Output DST
// found on file: gmn_replayed_11251_stream0_seg0_0.root
//////////////////////////////////////////////////////////

#ifndef gmn_rec_tree_h
#define gmn_rec_tree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class gmn_rec_tree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           Ndata_bb_eps_over_etot;
   Double_t        bb_eps_over_etot[1];   //[Ndata.bb.eps_over_etot]
   Int_t           Ndata_bb_etot_over_p;
   Double_t        bb_etot_over_p[6];   //[Ndata.bb.etot_over_p]
   Int_t           Ndata_bb_gem_hit_ADCU;
   Double_t        bb_gem_hit_ADCU[21];   //[Ndata.bb.gem.hit.ADCU]
   Int_t           Ndata_bb_gem_hit_ADCV;
   Double_t        bb_gem_hit_ADCV[21];   //[Ndata.bb.gem.hit.ADCV]
   Int_t           Ndata_bb_gem_hit_ADCasym;
   Double_t        bb_gem_hit_ADCasym[21];   //[Ndata.bb.gem.hit.ADCasym]
   Int_t           Ndata_bb_gem_hit_ADCavg;
   Double_t        bb_gem_hit_ADCavg[21];   //[Ndata.bb.gem.hit.ADCavg]
   Int_t           Ndata_bb_gem_hit_ADCmaxsampU;
   Double_t        bb_gem_hit_ADCmaxsampU[21];   //[Ndata.bb.gem.hit.ADCmaxsampU]
   Int_t           Ndata_bb_gem_hit_ADCmaxsampUclust;
   Double_t        bb_gem_hit_ADCmaxsampUclust[21];   //[Ndata.bb.gem.hit.ADCmaxsampUclust]
   Int_t           Ndata_bb_gem_hit_ADCmaxsampV;
   Double_t        bb_gem_hit_ADCmaxsampV[21];   //[Ndata.bb.gem.hit.ADCmaxsampV]
   Int_t           Ndata_bb_gem_hit_ADCmaxsampVclust;
   Double_t        bb_gem_hit_ADCmaxsampVclust[21];   //[Ndata.bb.gem.hit.ADCmaxsampVclust]
   Int_t           Ndata_bb_gem_hit_ADCmaxstripU;
   Double_t        bb_gem_hit_ADCmaxstripU[21];   //[Ndata.bb.gem.hit.ADCmaxstripU]
   Int_t           Ndata_bb_gem_hit_ADCmaxstripV;
   Double_t        bb_gem_hit_ADCmaxstripV[21];   //[Ndata.bb.gem.hit.ADCmaxstripV]
   Int_t           Ndata_bb_gem_hit_BUILD_ALL_SAMPLES_U;
   Double_t        bb_gem_hit_BUILD_ALL_SAMPLES_U[21];   //[Ndata.bb.gem.hit.BUILD_ALL_SAMPLES_U]
   Int_t           Ndata_bb_gem_hit_BUILD_ALL_SAMPLES_V;
   Double_t        bb_gem_hit_BUILD_ALL_SAMPLES_V[21];   //[Ndata.bb.gem.hit.BUILD_ALL_SAMPLES_V]
   Int_t           Ndata_bb_gem_hit_CM_GOOD_U;
   Double_t        bb_gem_hit_CM_GOOD_U[21];   //[Ndata.bb.gem.hit.CM_GOOD_U]
   Int_t           Ndata_bb_gem_hit_CM_GOOD_V;
   Double_t        bb_gem_hit_CM_GOOD_V[21];   //[Ndata.bb.gem.hit.CM_GOOD_V]
   Int_t           Ndata_bb_gem_hit_ENABLE_CM_U;
   Double_t        bb_gem_hit_ENABLE_CM_U[21];   //[Ndata.bb.gem.hit.ENABLE_CM_U]
   Int_t           Ndata_bb_gem_hit_ENABLE_CM_V;
   Double_t        bb_gem_hit_ENABLE_CM_V[21];   //[Ndata.bb.gem.hit.ENABLE_CM_V]
   Int_t           Ndata_bb_gem_hit_Tavg;
   Double_t        bb_gem_hit_Tavg[21];   //[Ndata.bb.gem.hit.Tavg]
   Int_t           Ndata_bb_gem_hit_Utime;
   Double_t        bb_gem_hit_Utime[21];   //[Ndata.bb.gem.hit.Utime]
   Int_t           Ndata_bb_gem_hit_UtimeMaxStrip;
   Double_t        bb_gem_hit_UtimeMaxStrip[21];   //[Ndata.bb.gem.hit.UtimeMaxStrip]
   Int_t           Ndata_bb_gem_hit_Vtime;
   Double_t        bb_gem_hit_Vtime[21];   //[Ndata.bb.gem.hit.Vtime]
   Int_t           Ndata_bb_gem_hit_VtimeMaxStrip;
   Double_t        bb_gem_hit_VtimeMaxStrip[21];   //[Ndata.bb.gem.hit.VtimeMaxStrip]
   Int_t           Ndata_bb_gem_hit_ccor_clust;
   Double_t        bb_gem_hit_ccor_clust[21];   //[Ndata.bb.gem.hit.ccor_clust]
   Int_t           Ndata_bb_gem_hit_ccor_strip;
   Double_t        bb_gem_hit_ccor_strip[21];   //[Ndata.bb.gem.hit.ccor_strip]
   Int_t           Ndata_bb_gem_hit_deltat;
   Double_t        bb_gem_hit_deltat[21];   //[Ndata.bb.gem.hit.deltat]
   Int_t           Ndata_bb_gem_hit_eresidu;
   Double_t        bb_gem_hit_eresidu[21];   //[Ndata.bb.gem.hit.eresidu]
   Int_t           Ndata_bb_gem_hit_eresidv;
   Double_t        bb_gem_hit_eresidv[21];   //[Ndata.bb.gem.hit.eresidv]
   Int_t           Ndata_bb_gem_hit_isampmaxUclust;
   Double_t        bb_gem_hit_isampmaxUclust[21];   //[Ndata.bb.gem.hit.isampmaxUclust]
   Int_t           Ndata_bb_gem_hit_isampmaxUstrip;
   Double_t        bb_gem_hit_isampmaxUstrip[21];   //[Ndata.bb.gem.hit.isampmaxUstrip]
   Int_t           Ndata_bb_gem_hit_isampmaxVclust;
   Double_t        bb_gem_hit_isampmaxVclust[21];   //[Ndata.bb.gem.hit.isampmaxVclust]
   Int_t           Ndata_bb_gem_hit_isampmaxVstrip;
   Double_t        bb_gem_hit_isampmaxVstrip[21];   //[Ndata.bb.gem.hit.isampmaxVstrip]
   Int_t           Ndata_bb_gem_hit_layer;
   Double_t        bb_gem_hit_layer[21];   //[Ndata.bb.gem.hit.layer]
   Int_t           Ndata_bb_gem_hit_module;
   Double_t        bb_gem_hit_module[21];   //[Ndata.bb.gem.hit.module]
   Int_t           Ndata_bb_gem_hit_nstripu;
   Double_t        bb_gem_hit_nstripu[21];   //[Ndata.bb.gem.hit.nstripu]
   Int_t           Ndata_bb_gem_hit_nstripv;
   Double_t        bb_gem_hit_nstripv[21];   //[Ndata.bb.gem.hit.nstripv]
   Int_t           Ndata_bb_gem_hit_residu;
   Double_t        bb_gem_hit_residu[21];   //[Ndata.bb.gem.hit.residu]
   Int_t           Ndata_bb_gem_hit_residv;
   Double_t        bb_gem_hit_residv[21];   //[Ndata.bb.gem.hit.residv]
   Int_t           Ndata_bb_gem_hit_trackindex;
   Double_t        bb_gem_hit_trackindex[21];   //[Ndata.bb.gem.hit.trackindex]
   Int_t           Ndata_bb_gem_hit_u;
   Double_t        bb_gem_hit_u[21];   //[Ndata.bb.gem.hit.u]
   Int_t           Ndata_bb_gem_hit_umoment;
   Double_t        bb_gem_hit_umoment[21];   //[Ndata.bb.gem.hit.umoment]
   Int_t           Ndata_bb_gem_hit_usigma;
   Double_t        bb_gem_hit_usigma[21];   //[Ndata.bb.gem.hit.usigma]
   Int_t           Ndata_bb_gem_hit_ustriphi;
   Double_t        bb_gem_hit_ustriphi[21];   //[Ndata.bb.gem.hit.ustriphi]
   Int_t           Ndata_bb_gem_hit_ustriplo;
   Double_t        bb_gem_hit_ustriplo[21];   //[Ndata.bb.gem.hit.ustriplo]
   Int_t           Ndata_bb_gem_hit_ustripmax;
   Double_t        bb_gem_hit_ustripmax[21];   //[Ndata.bb.gem.hit.ustripmax]
   Int_t           Ndata_bb_gem_hit_v;
   Double_t        bb_gem_hit_v[21];   //[Ndata.bb.gem.hit.v]
   Int_t           Ndata_bb_gem_hit_vmoment;
   Double_t        bb_gem_hit_vmoment[21];   //[Ndata.bb.gem.hit.vmoment]
   Int_t           Ndata_bb_gem_hit_vsigma;
   Double_t        bb_gem_hit_vsigma[21];   //[Ndata.bb.gem.hit.vsigma]
   Int_t           Ndata_bb_gem_hit_vstriphi;
   Double_t        bb_gem_hit_vstriphi[21];   //[Ndata.bb.gem.hit.vstriphi]
   Int_t           Ndata_bb_gem_hit_vstriplo;
   Double_t        bb_gem_hit_vstriplo[21];   //[Ndata.bb.gem.hit.vstriplo]
   Int_t           Ndata_bb_gem_hit_vstripmax;
   Double_t        bb_gem_hit_vstripmax[21];   //[Ndata.bb.gem.hit.vstripmax]
   Int_t           Ndata_bb_gem_hit_xglobal;
   Double_t        bb_gem_hit_xglobal[21];   //[Ndata.bb.gem.hit.xglobal]
   Int_t           Ndata_bb_gem_hit_xlocal;
   Double_t        bb_gem_hit_xlocal[21];   //[Ndata.bb.gem.hit.xlocal]
   Int_t           Ndata_bb_gem_hit_yglobal;
   Double_t        bb_gem_hit_yglobal[21];   //[Ndata.bb.gem.hit.yglobal]
   Int_t           Ndata_bb_gem_hit_ylocal;
   Double_t        bb_gem_hit_ylocal[21];   //[Ndata.bb.gem.hit.ylocal]
   Int_t           Ndata_bb_gem_hit_zglobal;
   Double_t        bb_gem_hit_zglobal[21];   //[Ndata.bb.gem.hit.zglobal]
   Int_t           Ndata_bb_gem_n2Dhit_layer;
   Double_t        bb_gem_n2Dhit_layer[5];   //[Ndata.bb.gem.n2Dhit_layer]
   Int_t           Ndata_bb_gem_nclustu_layer;
   Double_t        bb_gem_nclustu_layer[5];   //[Ndata.bb.gem.nclustu_layer]
   Int_t           Ndata_bb_gem_nclustv_layer;
   Double_t        bb_gem_nclustv_layer[5];   //[Ndata.bb.gem.nclustv_layer]
   Int_t           Ndata_bb_gem_nstripsu_layer;
   Double_t        bb_gem_nstripsu_layer[5];   //[Ndata.bb.gem.nstripsu_layer]
   Int_t           Ndata_bb_gem_nstripsv_layer;
   Double_t        bb_gem_nstripsv_layer[5];   //[Ndata.bb.gem.nstripsv_layer]
   Int_t           Ndata_bb_gem_track_chi2ndf;
   Double_t        bb_gem_track_chi2ndf[7];   //[Ndata.bb.gem.track.chi2ndf]
   Int_t           Ndata_bb_gem_track_nhits;
   Double_t        bb_gem_track_nhits[7];   //[Ndata.bb.gem.track.nhits]
   Int_t           Ndata_bb_gem_track_x;
   Double_t        bb_gem_track_x[7];   //[Ndata.bb.gem.track.x]
   Int_t           Ndata_bb_gem_track_xp;
   Double_t        bb_gem_track_xp[7];   //[Ndata.bb.gem.track.xp]
   Int_t           Ndata_bb_gem_track_y;
   Double_t        bb_gem_track_y[7];   //[Ndata.bb.gem.track.y]
   Int_t           Ndata_bb_gem_track_yp;
   Double_t        bb_gem_track_yp[7];   //[Ndata.bb.gem.track.yp]
   Int_t           Ndata_bb_grinch_tdc_tdc;
   Double_t        bb_grinch_tdc_tdc[14];   //[Ndata.bb.grinch_tdc.tdc]
   Int_t           Ndata_bb_grinch_tdc_tdc_mult;
   Double_t        bb_grinch_tdc_tdc_mult[14];   //[Ndata.bb.grinch_tdc.tdc_mult]
   Int_t           Ndata_bb_grinch_tdc_tdc_te;
   Double_t        bb_grinch_tdc_tdc_te[14];   //[Ndata.bb.grinch_tdc.tdc_te]
   Int_t           Ndata_bb_grinch_tdc_tdc_tot;
   Double_t        bb_grinch_tdc_tdc_tot[14];   //[Ndata.bb.grinch_tdc.tdc_tot]
   Int_t           Ndata_bb_grinch_tdc_tdccol;
   Double_t        bb_grinch_tdc_tdccol[14];   //[Ndata.bb.grinch_tdc.tdccol]
   Int_t           Ndata_bb_grinch_tdc_tdcelemID;
   Double_t        bb_grinch_tdc_tdcelemID[14];   //[Ndata.bb.grinch_tdc.tdcelemID]
   Int_t           Ndata_bb_grinch_tdc_tdclayer;
   Double_t        bb_grinch_tdc_tdclayer[14];   //[Ndata.bb.grinch_tdc.tdclayer]
   Int_t           Ndata_bb_grinch_tdc_tdcrow;
   Double_t        bb_grinch_tdc_tdcrow[14];   //[Ndata.bb.grinch_tdc.tdcrow]
   Int_t           Ndata_bb_hodoadc_bar_adc_L_a;
   Double_t        bb_hodoadc_bar_adc_L_a[20];   //[Ndata.bb.hodoadc.bar.adc.L.a]
   Int_t           Ndata_bb_hodoadc_bar_adc_L_ac;
   Double_t        bb_hodoadc_bar_adc_L_ac[20];   //[Ndata.bb.hodoadc.bar.adc.L.ac]
   Int_t           Ndata_bb_hodoadc_bar_adc_L_ap;
   Double_t        bb_hodoadc_bar_adc_L_ap[20];   //[Ndata.bb.hodoadc.bar.adc.L.ap]
   Int_t           Ndata_bb_hodoadc_bar_adc_R_a;
   Double_t        bb_hodoadc_bar_adc_R_a[20];   //[Ndata.bb.hodoadc.bar.adc.R.a]
   Int_t           Ndata_bb_hodoadc_bar_adc_R_ac;
   Double_t        bb_hodoadc_bar_adc_R_ac[20];   //[Ndata.bb.hodoadc.bar.adc.R.ac]
   Int_t           Ndata_bb_hodoadc_bar_adc_R_ap;
   Double_t        bb_hodoadc_bar_adc_R_ap[20];   //[Ndata.bb.hodoadc.bar.adc.R.ap]
   Int_t           Ndata_bb_hodoadc_bar_adc_id;
   Double_t        bb_hodoadc_bar_adc_id[20];   //[Ndata.bb.hodoadc.bar.adc.id]
   Int_t           Ndata_bb_hodoadc_bar_adc_mean;
   Double_t        bb_hodoadc_bar_adc_mean[20];   //[Ndata.bb.hodoadc.bar.adc.mean]
   Int_t           Ndata_bb_hodotdc_bar_tdc_L_le;
   Double_t        bb_hodotdc_bar_tdc_L_le[66];   //[Ndata.bb.hodotdc.bar.tdc.L.le]
   Int_t           Ndata_bb_hodotdc_bar_tdc_L_leW;
   Double_t        bb_hodotdc_bar_tdc_L_leW[66];   //[Ndata.bb.hodotdc.bar.tdc.L.leW]
   Int_t           Ndata_bb_hodotdc_bar_tdc_L_te;
   Double_t        bb_hodotdc_bar_tdc_L_te[66];   //[Ndata.bb.hodotdc.bar.tdc.L.te]
   Int_t           Ndata_bb_hodotdc_bar_tdc_L_teW;
   Double_t        bb_hodotdc_bar_tdc_L_teW[66];   //[Ndata.bb.hodotdc.bar.tdc.L.teW]
   Int_t           Ndata_bb_hodotdc_bar_tdc_L_tot;
   Double_t        bb_hodotdc_bar_tdc_L_tot[66];   //[Ndata.bb.hodotdc.bar.tdc.L.tot]
   Int_t           Ndata_bb_hodotdc_bar_tdc_L_totW;
   Double_t        bb_hodotdc_bar_tdc_L_totW[66];   //[Ndata.bb.hodotdc.bar.tdc.L.totW]
   Int_t           Ndata_bb_hodotdc_bar_tdc_R_le;
   Double_t        bb_hodotdc_bar_tdc_R_le[66];   //[Ndata.bb.hodotdc.bar.tdc.R.le]
   Int_t           Ndata_bb_hodotdc_bar_tdc_R_leW;
   Double_t        bb_hodotdc_bar_tdc_R_leW[66];   //[Ndata.bb.hodotdc.bar.tdc.R.leW]
   Int_t           Ndata_bb_hodotdc_bar_tdc_R_te;
   Double_t        bb_hodotdc_bar_tdc_R_te[66];   //[Ndata.bb.hodotdc.bar.tdc.R.te]
   Int_t           Ndata_bb_hodotdc_bar_tdc_R_teW;
   Double_t        bb_hodotdc_bar_tdc_R_teW[66];   //[Ndata.bb.hodotdc.bar.tdc.R.teW]
   Int_t           Ndata_bb_hodotdc_bar_tdc_R_tot;
   Double_t        bb_hodotdc_bar_tdc_R_tot[66];   //[Ndata.bb.hodotdc.bar.tdc.R.tot]
   Int_t           Ndata_bb_hodotdc_bar_tdc_R_totW;
   Double_t        bb_hodotdc_bar_tdc_R_totW[66];   //[Ndata.bb.hodotdc.bar.tdc.R.totW]
   Int_t           Ndata_bb_hodotdc_bar_tdc_id;
   Double_t        bb_hodotdc_bar_tdc_id[66];   //[Ndata.bb.hodotdc.bar.tdc.id]
   Int_t           Ndata_bb_hodotdc_bar_tdc_meantime;
   Double_t        bb_hodotdc_bar_tdc_meantime[66];   //[Ndata.bb.hodotdc.bar.tdc.meantime]
   Int_t           Ndata_bb_hodotdc_bar_tdc_timediff;
   Double_t        bb_hodotdc_bar_tdc_timediff[66];   //[Ndata.bb.hodotdc.bar.tdc.timediff]
   Int_t           Ndata_bb_hodotdc_bar_tdc_timehitpos;
   Double_t        bb_hodotdc_bar_tdc_timehitpos[66];   //[Ndata.bb.hodotdc.bar.tdc.timehitpos]
   Int_t           Ndata_bb_hodotdc_bar_tdc_vpos;
   Double_t        bb_hodotdc_bar_tdc_vpos[66];   //[Ndata.bb.hodotdc.bar.tdc.vpos]
   Int_t           Ndata_bb_ps_clus_col;
   Double_t        bb_ps_clus_col[5];   //[Ndata.bb.ps.clus.col]
   Int_t           Ndata_bb_ps_clus_e;
   Double_t        bb_ps_clus_e[5];   //[Ndata.bb.ps.clus.e]
   Int_t           Ndata_bb_ps_clus_e_c;
   Double_t        bb_ps_clus_e_c[5];   //[Ndata.bb.ps.clus.e_c]
   Int_t           Ndata_bb_ps_clus_eblk;
   Double_t        bb_ps_clus_eblk[5];   //[Ndata.bb.ps.clus.eblk]
   Int_t           Ndata_bb_ps_clus_eblk_c;
   Double_t        bb_ps_clus_eblk_c[5];   //[Ndata.bb.ps.clus.eblk_c]
   Int_t           Ndata_bb_ps_clus_id;
   Double_t        bb_ps_clus_id[5];   //[Ndata.bb.ps.clus.id]
   Int_t           Ndata_bb_ps_clus_nblk;
   Double_t        bb_ps_clus_nblk[5];   //[Ndata.bb.ps.clus.nblk]
   Int_t           Ndata_bb_ps_clus_row;
   Double_t        bb_ps_clus_row[5];   //[Ndata.bb.ps.clus.row]
   Int_t           Ndata_bb_ps_clus_x;
   Double_t        bb_ps_clus_x[5];   //[Ndata.bb.ps.clus.x]
   Int_t           Ndata_bb_ps_clus_y;
   Double_t        bb_ps_clus_y[5];   //[Ndata.bb.ps.clus.y]
   Int_t           Ndata_bb_ps_clus_blk_col;
   Double_t        bb_ps_clus_blk_col[5];   //[Ndata.bb.ps.clus_blk.col]
   Int_t           Ndata_bb_ps_clus_blk_e;
   Double_t        bb_ps_clus_blk_e[5];   //[Ndata.bb.ps.clus_blk.e]
   Int_t           Ndata_bb_ps_clus_blk_e_c;
   Double_t        bb_ps_clus_blk_e_c[5];   //[Ndata.bb.ps.clus_blk.e_c]
   Int_t           Ndata_bb_ps_clus_blk_id;
   Double_t        bb_ps_clus_blk_id[5];   //[Ndata.bb.ps.clus_blk.id]
   Int_t           Ndata_bb_ps_clus_blk_row;
   Double_t        bb_ps_clus_blk_row[5];   //[Ndata.bb.ps.clus_blk.row]
   Int_t           Ndata_bb_ps_clus_blk_x;
   Double_t        bb_ps_clus_blk_x[5];   //[Ndata.bb.ps.clus_blk.x]
   Int_t           Ndata_bb_ps_clus_blk_y;
   Double_t        bb_ps_clus_blk_y[5];   //[Ndata.bb.ps.clus_blk.y]
   Int_t           Ndata_bb_sh_clus_col;
   Double_t        bb_sh_clus_col[5];   //[Ndata.bb.sh.clus.col]
   Int_t           Ndata_bb_sh_clus_e;
   Double_t        bb_sh_clus_e[5];   //[Ndata.bb.sh.clus.e]
   Int_t           Ndata_bb_sh_clus_e_c;
   Double_t        bb_sh_clus_e_c[5];   //[Ndata.bb.sh.clus.e_c]
   Int_t           Ndata_bb_sh_clus_eblk;
   Double_t        bb_sh_clus_eblk[5];   //[Ndata.bb.sh.clus.eblk]
   Int_t           Ndata_bb_sh_clus_eblk_c;
   Double_t        bb_sh_clus_eblk_c[5];   //[Ndata.bb.sh.clus.eblk_c]
   Int_t           Ndata_bb_sh_clus_id;
   Double_t        bb_sh_clus_id[5];   //[Ndata.bb.sh.clus.id]
   Int_t           Ndata_bb_sh_clus_nblk;
   Double_t        bb_sh_clus_nblk[5];   //[Ndata.bb.sh.clus.nblk]
   Int_t           Ndata_bb_sh_clus_row;
   Double_t        bb_sh_clus_row[5];   //[Ndata.bb.sh.clus.row]
   Int_t           Ndata_bb_sh_clus_x;
   Double_t        bb_sh_clus_x[5];   //[Ndata.bb.sh.clus.x]
   Int_t           Ndata_bb_sh_clus_y;
   Double_t        bb_sh_clus_y[5];   //[Ndata.bb.sh.clus.y]
   Int_t           Ndata_bb_sh_clus_blk_col;
   Double_t        bb_sh_clus_blk_col[12];   //[Ndata.bb.sh.clus_blk.col]
   Int_t           Ndata_bb_sh_clus_blk_e;
   Double_t        bb_sh_clus_blk_e[12];   //[Ndata.bb.sh.clus_blk.e]
   Int_t           Ndata_bb_sh_clus_blk_e_c;
   Double_t        bb_sh_clus_blk_e_c[12];   //[Ndata.bb.sh.clus_blk.e_c]
   Int_t           Ndata_bb_sh_clus_blk_id;
   Double_t        bb_sh_clus_blk_id[12];   //[Ndata.bb.sh.clus_blk.id]
   Int_t           Ndata_bb_sh_clus_blk_row;
   Double_t        bb_sh_clus_blk_row[12];   //[Ndata.bb.sh.clus_blk.row]
   Int_t           Ndata_bb_sh_clus_blk_x;
   Double_t        bb_sh_clus_blk_x[12];   //[Ndata.bb.sh.clus_blk.x]
   Int_t           Ndata_bb_sh_clus_blk_y;
   Double_t        bb_sh_clus_blk_y[12];   //[Ndata.bb.sh.clus_blk.y]
   Int_t           Ndata_bb_tr_beta;
   Double_t        bb_tr_beta[7];   //[Ndata.bb.tr.beta]
   Int_t           Ndata_bb_tr_chi2;
   Double_t        bb_tr_chi2[7];   //[Ndata.bb.tr.chi2]
   Int_t           Ndata_bb_tr_d_ph;
   Double_t        bb_tr_d_ph[7];   //[Ndata.bb.tr.d_ph]
   Int_t           Ndata_bb_tr_d_th;
   Double_t        bb_tr_d_th[7];   //[Ndata.bb.tr.d_th]
   Int_t           Ndata_bb_tr_d_x;
   Double_t        bb_tr_d_x[7];   //[Ndata.bb.tr.d_x]
   Int_t           Ndata_bb_tr_d_y;
   Double_t        bb_tr_d_y[7];   //[Ndata.bb.tr.d_y]
   Int_t           Ndata_bb_tr_dbeta;
   Double_t        bb_tr_dbeta[7];   //[Ndata.bb.tr.dbeta]
   Int_t           Ndata_bb_tr_dtime;
   Double_t        bb_tr_dtime[7];   //[Ndata.bb.tr.dtime]
   Int_t           Ndata_bb_tr_flag;
   Double_t        bb_tr_flag[7];   //[Ndata.bb.tr.flag]
   Int_t           Ndata_bb_tr_ndof;
   Double_t        bb_tr_ndof[7];   //[Ndata.bb.tr.ndof]
   Int_t           Ndata_bb_tr_p;
   Double_t        bb_tr_p[7];   //[Ndata.bb.tr.p]
   Int_t           Ndata_bb_tr_pathl;
   Double_t        bb_tr_pathl[7];   //[Ndata.bb.tr.pathl]
   Int_t           Ndata_bb_tr_ph;
   Double_t        bb_tr_ph[7];   //[Ndata.bb.tr.ph]
   Int_t           Ndata_bb_tr_px;
   Double_t        bb_tr_px[7];   //[Ndata.bb.tr.px]
   Int_t           Ndata_bb_tr_py;
   Double_t        bb_tr_py[7];   //[Ndata.bb.tr.py]
   Int_t           Ndata_bb_tr_pz;
   Double_t        bb_tr_pz[7];   //[Ndata.bb.tr.pz]
   Int_t           Ndata_bb_tr_r_ph;
   Double_t        bb_tr_r_ph[7];   //[Ndata.bb.tr.r_ph]
   Int_t           Ndata_bb_tr_r_th;
   Double_t        bb_tr_r_th[7];   //[Ndata.bb.tr.r_th]
   Int_t           Ndata_bb_tr_r_x;
   Double_t        bb_tr_r_x[7];   //[Ndata.bb.tr.r_x]
   Int_t           Ndata_bb_tr_r_y;
   Double_t        bb_tr_r_y[7];   //[Ndata.bb.tr.r_y]
   Int_t           Ndata_bb_tr_tg_dp;
   Double_t        bb_tr_tg_dp[7];   //[Ndata.bb.tr.tg_dp]
   Int_t           Ndata_bb_tr_tg_ph;
   Double_t        bb_tr_tg_ph[7];   //[Ndata.bb.tr.tg_ph]
   Int_t           Ndata_bb_tr_tg_th;
   Double_t        bb_tr_tg_th[7];   //[Ndata.bb.tr.tg_th]
   Int_t           Ndata_bb_tr_tg_x;
   Double_t        bb_tr_tg_x[7];   //[Ndata.bb.tr.tg_x]
   Int_t           Ndata_bb_tr_tg_y;
   Double_t        bb_tr_tg_y[7];   //[Ndata.bb.tr.tg_y]
   Int_t           Ndata_bb_tr_th;
   Double_t        bb_tr_th[7];   //[Ndata.bb.tr.th]
   Int_t           Ndata_bb_tr_time;
   Double_t        bb_tr_time[7];   //[Ndata.bb.tr.time]
   Int_t           Ndata_bb_tr_vx;
   Double_t        bb_tr_vx[7];   //[Ndata.bb.tr.vx]
   Int_t           Ndata_bb_tr_vy;
   Double_t        bb_tr_vy[7];   //[Ndata.bb.tr.vy]
   Int_t           Ndata_bb_tr_vz;
   Double_t        bb_tr_vz[7];   //[Ndata.bb.tr.vz]
   Int_t           Ndata_bb_tr_x;
   Double_t        bb_tr_x[7];   //[Ndata.bb.tr.x]
   Int_t           Ndata_bb_tr_y;
   Double_t        bb_tr_y[7];   //[Ndata.bb.tr.y]
   Int_t           Ndata_bb_x_bcp;
   Double_t        bb_x_bcp[1];   //[Ndata.bb.x_bcp]
   Int_t           Ndata_bb_x_fcp;
   Double_t        bb_x_fcp[1];   //[Ndata.bb.x_fcp]
   Int_t           Ndata_bb_y_bcp;
   Double_t        bb_y_bcp[1];   //[Ndata.bb.y_bcp]
   Int_t           Ndata_bb_y_fcp;
   Double_t        bb_y_fcp[1];   //[Ndata.bb.y_fcp]
   Int_t           Ndata_sbs_hcal_a_amp;
   Double_t        sbs_hcal_a_amp[16];   //[Ndata.sbs.hcal.a_amp]
   Int_t           Ndata_sbs_hcal_a_c;
   Double_t        sbs_hcal_a_c[16];   //[Ndata.sbs.hcal.a_c]
   Int_t           Ndata_sbs_hcal_a_p;
   Double_t        sbs_hcal_a_p[16];   //[Ndata.sbs.hcal.a_p]
   Int_t           Ndata_sbs_hcal_a_time;
   Double_t        sbs_hcal_a_time[16];   //[Ndata.sbs.hcal.a_time]
   Int_t           Ndata_sbs_hcal_adccol;
   Double_t        sbs_hcal_adccol[16];   //[Ndata.sbs.hcal.adccol]
   Int_t           Ndata_sbs_hcal_adcrow;
   Double_t        sbs_hcal_adcrow[16];   //[Ndata.sbs.hcal.adcrow]
   Int_t           Ndata_sbs_hcal_clus_col;
   Double_t        sbs_hcal_clus_col[10];   //[Ndata.sbs.hcal.clus.col]
   Int_t           Ndata_sbs_hcal_clus_e;
   Double_t        sbs_hcal_clus_e[10];   //[Ndata.sbs.hcal.clus.e]
   Int_t           Ndata_sbs_hcal_clus_eblk;
   Double_t        sbs_hcal_clus_eblk[10];   //[Ndata.sbs.hcal.clus.eblk]
   Int_t           Ndata_sbs_hcal_clus_id;
   Double_t        sbs_hcal_clus_id[10];   //[Ndata.sbs.hcal.clus.id]
   Int_t           Ndata_sbs_hcal_clus_nblk;
   Double_t        sbs_hcal_clus_nblk[10];   //[Ndata.sbs.hcal.clus.nblk]
   Int_t           Ndata_sbs_hcal_clus_row;
   Double_t        sbs_hcal_clus_row[10];   //[Ndata.sbs.hcal.clus.row]
   Int_t           Ndata_sbs_hcal_clus_blk_e;
   Double_t        sbs_hcal_clus_blk_e[6];   //[Ndata.sbs.hcal.clus_blk.e]
   Int_t           Ndata_sbs_hcal_clus_blk_id;
   Double_t        sbs_hcal_clus_blk_id[6];   //[Ndata.sbs.hcal.clus_blk.id]
   Double_t        BB_gold_beta;
   Double_t        BB_gold_dp;
   Double_t        BB_gold_index;
   Double_t        BB_gold_ok;
   Double_t        BB_gold_p;
   Double_t        BB_gold_ph;
   Double_t        BB_gold_px;
   Double_t        BB_gold_py;
   Double_t        BB_gold_pz;
   Double_t        BB_gold_th;
   Double_t        BB_gold_x;
   Double_t        BB_gold_y;
   Double_t        bb_gem_hit_ngoodhits;
   Double_t        bb_gem_nlayershit;
   Double_t        bb_gem_nlayershitu;
   Double_t        bb_gem_nlayershituv;
   Double_t        bb_gem_nlayershitv;
   Double_t        bb_gem_track_besttrack;
   Double_t        bb_gem_track_ntrack;
   Double_t        bb_grinch_tdc_ngoodADChits;
   Double_t        bb_grinch_tdc_ngoodTDChits;
   Double_t        bb_grinch_tdc_nhits;
   Double_t        bb_grinch_tdc_nrefhits;
   Double_t        bb_hodotdc_bar_ngoodbars;
   Double_t        bb_hodotdc_nclus;
   Double_t        bb_ps_e;
   Double_t        bb_ps_e_c;
   Double_t        bb_ps_nblk;
   Double_t        bb_ps_nclus;
   Double_t        bb_ps_ngoodADChits;
   Double_t        bb_ps_x;
   Double_t        bb_ps_y;
   Double_t        bb_sh_e;
   Double_t        bb_sh_e_c;
   Double_t        bb_sh_nblk;
   Double_t        bb_sh_nclus;
   Double_t        bb_sh_ngoodADChits;
   Double_t        bb_sh_x;
   Double_t        bb_sh_y;
   Double_t        bb_tr_n;
   Double_t        sbs_hcal_e;
   Double_t        sbs_hcal_nblk;
   Double_t        sbs_hcal_nclus;
   Double_t        sbs_hcal_ngoodADChits;
   Double_t        singletrack;
   Double_t        anytrack;
 //THaEvent        *Event_Branch;
   ULong64_t       fEvtHdr_fEvtTime;
   UInt_t          fEvtHdr_fEvtNum;
   UInt_t          fEvtHdr_fEvtType;
   UInt_t          fEvtHdr_fEvtLen;
   Int_t           fEvtHdr_fHelicity;
   Int_t           fEvtHdr_fTargetPol;
   UInt_t          fEvtHdr_fRun;

   // List of branches
   TBranch        *b_Ndata_bb_eps_over_etot;   //!
   TBranch        *b_bb_eps_over_etot;   //!
   TBranch        *b_Ndata_bb_etot_over_p;   //!
   TBranch        *b_bb_etot_over_p;   //!
   TBranch        *b_Ndata_bb_gem_hit_ADCU;   //!
   TBranch        *b_bb_gem_hit_ADCU;   //!
   TBranch        *b_Ndata_bb_gem_hit_ADCV;   //!
   TBranch        *b_bb_gem_hit_ADCV;   //!
   TBranch        *b_Ndata_bb_gem_hit_ADCasym;   //!
   TBranch        *b_bb_gem_hit_ADCasym;   //!
   TBranch        *b_Ndata_bb_gem_hit_ADCavg;   //!
   TBranch        *b_bb_gem_hit_ADCavg;   //!
   TBranch        *b_Ndata_bb_gem_hit_ADCmaxsampU;   //!
   TBranch        *b_bb_gem_hit_ADCmaxsampU;   //!
   TBranch        *b_Ndata_bb_gem_hit_ADCmaxsampUclust;   //!
   TBranch        *b_bb_gem_hit_ADCmaxsampUclust;   //!
   TBranch        *b_Ndata_bb_gem_hit_ADCmaxsampV;   //!
   TBranch        *b_bb_gem_hit_ADCmaxsampV;   //!
   TBranch        *b_Ndata_bb_gem_hit_ADCmaxsampVclust;   //!
   TBranch        *b_bb_gem_hit_ADCmaxsampVclust;   //!
   TBranch        *b_Ndata_bb_gem_hit_ADCmaxstripU;   //!
   TBranch        *b_bb_gem_hit_ADCmaxstripU;   //!
   TBranch        *b_Ndata_bb_gem_hit_ADCmaxstripV;   //!
   TBranch        *b_bb_gem_hit_ADCmaxstripV;   //!
   TBranch        *b_Ndata_bb_gem_hit_BUILD_ALL_SAMPLES_U;   //!
   TBranch        *b_bb_gem_hit_BUILD_ALL_SAMPLES_U;   //!
   TBranch        *b_Ndata_bb_gem_hit_BUILD_ALL_SAMPLES_V;   //!
   TBranch        *b_bb_gem_hit_BUILD_ALL_SAMPLES_V;   //!
   TBranch        *b_Ndata_bb_gem_hit_CM_GOOD_U;   //!
   TBranch        *b_bb_gem_hit_CM_GOOD_U;   //!
   TBranch        *b_Ndata_bb_gem_hit_CM_GOOD_V;   //!
   TBranch        *b_bb_gem_hit_CM_GOOD_V;   //!
   TBranch        *b_Ndata_bb_gem_hit_ENABLE_CM_U;   //!
   TBranch        *b_bb_gem_hit_ENABLE_CM_U;   //!
   TBranch        *b_Ndata_bb_gem_hit_ENABLE_CM_V;   //!
   TBranch        *b_bb_gem_hit_ENABLE_CM_V;   //!
   TBranch        *b_Ndata_bb_gem_hit_Tavg;   //!
   TBranch        *b_bb_gem_hit_Tavg;   //!
   TBranch        *b_Ndata_bb_gem_hit_Utime;   //!
   TBranch        *b_bb_gem_hit_Utime;   //!
   TBranch        *b_Ndata_bb_gem_hit_UtimeMaxStrip;   //!
   TBranch        *b_bb_gem_hit_UtimeMaxStrip;   //!
   TBranch        *b_Ndata_bb_gem_hit_Vtime;   //!
   TBranch        *b_bb_gem_hit_Vtime;   //!
   TBranch        *b_Ndata_bb_gem_hit_VtimeMaxStrip;   //!
   TBranch        *b_bb_gem_hit_VtimeMaxStrip;   //!
   TBranch        *b_Ndata_bb_gem_hit_ccor_clust;   //!
   TBranch        *b_bb_gem_hit_ccor_clust;   //!
   TBranch        *b_Ndata_bb_gem_hit_ccor_strip;   //!
   TBranch        *b_bb_gem_hit_ccor_strip;   //!
   TBranch        *b_Ndata_bb_gem_hit_deltat;   //!
   TBranch        *b_bb_gem_hit_deltat;   //!
   TBranch        *b_Ndata_bb_gem_hit_eresidu;   //!
   TBranch        *b_bb_gem_hit_eresidu;   //!
   TBranch        *b_Ndata_bb_gem_hit_eresidv;   //!
   TBranch        *b_bb_gem_hit_eresidv;   //!
   TBranch        *b_Ndata_bb_gem_hit_isampmaxUclust;   //!
   TBranch        *b_bb_gem_hit_isampmaxUclust;   //!
   TBranch        *b_Ndata_bb_gem_hit_isampmaxUstrip;   //!
   TBranch        *b_bb_gem_hit_isampmaxUstrip;   //!
   TBranch        *b_Ndata_bb_gem_hit_isampmaxVclust;   //!
   TBranch        *b_bb_gem_hit_isampmaxVclust;   //!
   TBranch        *b_Ndata_bb_gem_hit_isampmaxVstrip;   //!
   TBranch        *b_bb_gem_hit_isampmaxVstrip;   //!
   TBranch        *b_Ndata_bb_gem_hit_layer;   //!
   TBranch        *b_bb_gem_hit_layer;   //!
   TBranch        *b_Ndata_bb_gem_hit_module;   //!
   TBranch        *b_bb_gem_hit_module;   //!
   TBranch        *b_Ndata_bb_gem_hit_nstripu;   //!
   TBranch        *b_bb_gem_hit_nstripu;   //!
   TBranch        *b_Ndata_bb_gem_hit_nstripv;   //!
   TBranch        *b_bb_gem_hit_nstripv;   //!
   TBranch        *b_Ndata_bb_gem_hit_residu;   //!
   TBranch        *b_bb_gem_hit_residu;   //!
   TBranch        *b_Ndata_bb_gem_hit_residv;   //!
   TBranch        *b_bb_gem_hit_residv;   //!
   TBranch        *b_Ndata_bb_gem_hit_trackindex;   //!
   TBranch        *b_bb_gem_hit_trackindex;   //!
   TBranch        *b_Ndata_bb_gem_hit_u;   //!
   TBranch        *b_bb_gem_hit_u;   //!
   TBranch        *b_Ndata_bb_gem_hit_umoment;   //!
   TBranch        *b_bb_gem_hit_umoment;   //!
   TBranch        *b_Ndata_bb_gem_hit_usigma;   //!
   TBranch        *b_bb_gem_hit_usigma;   //!
   TBranch        *b_Ndata_bb_gem_hit_ustriphi;   //!
   TBranch        *b_bb_gem_hit_ustriphi;   //!
   TBranch        *b_Ndata_bb_gem_hit_ustriplo;   //!
   TBranch        *b_bb_gem_hit_ustriplo;   //!
   TBranch        *b_Ndata_bb_gem_hit_ustripmax;   //!
   TBranch        *b_bb_gem_hit_ustripmax;   //!
   TBranch        *b_Ndata_bb_gem_hit_v;   //!
   TBranch        *b_bb_gem_hit_v;   //!
   TBranch        *b_Ndata_bb_gem_hit_vmoment;   //!
   TBranch        *b_bb_gem_hit_vmoment;   //!
   TBranch        *b_Ndata_bb_gem_hit_vsigma;   //!
   TBranch        *b_bb_gem_hit_vsigma;   //!
   TBranch        *b_Ndata_bb_gem_hit_vstriphi;   //!
   TBranch        *b_bb_gem_hit_vstriphi;   //!
   TBranch        *b_Ndata_bb_gem_hit_vstriplo;   //!
   TBranch        *b_bb_gem_hit_vstriplo;   //!
   TBranch        *b_Ndata_bb_gem_hit_vstripmax;   //!
   TBranch        *b_bb_gem_hit_vstripmax;   //!
   TBranch        *b_Ndata_bb_gem_hit_xglobal;   //!
   TBranch        *b_bb_gem_hit_xglobal;   //!
   TBranch        *b_Ndata_bb_gem_hit_xlocal;   //!
   TBranch        *b_bb_gem_hit_xlocal;   //!
   TBranch        *b_Ndata_bb_gem_hit_yglobal;   //!
   TBranch        *b_bb_gem_hit_yglobal;   //!
   TBranch        *b_Ndata_bb_gem_hit_ylocal;   //!
   TBranch        *b_bb_gem_hit_ylocal;   //!
   TBranch        *b_Ndata_bb_gem_hit_zglobal;   //!
   TBranch        *b_bb_gem_hit_zglobal;   //!
   TBranch        *b_Ndata_bb_gem_n2Dhit_layer;   //!
   TBranch        *b_bb_gem_n2Dhit_layer;   //!
   TBranch        *b_Ndata_bb_gem_nclustu_layer;   //!
   TBranch        *b_bb_gem_nclustu_layer;   //!
   TBranch        *b_Ndata_bb_gem_nclustv_layer;   //!
   TBranch        *b_bb_gem_nclustv_layer;   //!
   TBranch        *b_Ndata_bb_gem_nstripsu_layer;   //!
   TBranch        *b_bb_gem_nstripsu_layer;   //!
   TBranch        *b_Ndata_bb_gem_nstripsv_layer;   //!
   TBranch        *b_bb_gem_nstripsv_layer;   //!
   TBranch        *b_Ndata_bb_gem_track_chi2ndf;   //!
   TBranch        *b_bb_gem_track_chi2ndf;   //!
   TBranch        *b_Ndata_bb_gem_track_nhits;   //!
   TBranch        *b_bb_gem_track_nhits;   //!
   TBranch        *b_Ndata_bb_gem_track_x;   //!
   TBranch        *b_bb_gem_track_x;   //!
   TBranch        *b_Ndata_bb_gem_track_xp;   //!
   TBranch        *b_bb_gem_track_xp;   //!
   TBranch        *b_Ndata_bb_gem_track_y;   //!
   TBranch        *b_bb_gem_track_y;   //!
   TBranch        *b_Ndata_bb_gem_track_yp;   //!
   TBranch        *b_bb_gem_track_yp;   //!
   TBranch        *b_Ndata_bb_grinch_tdc_tdc;   //!
   TBranch        *b_bb_grinch_tdc_tdc;   //!
   TBranch        *b_Ndata_bb_grinch_tdc_tdc_mult;   //!
   TBranch        *b_bb_grinch_tdc_tdc_mult;   //!
   TBranch        *b_Ndata_bb_grinch_tdc_tdc_te;   //!
   TBranch        *b_bb_grinch_tdc_tdc_te;   //!
   TBranch        *b_Ndata_bb_grinch_tdc_tdc_tot;   //!
   TBranch        *b_bb_grinch_tdc_tdc_tot;   //!
   TBranch        *b_Ndata_bb_grinch_tdc_tdccol;   //!
   TBranch        *b_bb_grinch_tdc_tdccol;   //!
   TBranch        *b_Ndata_bb_grinch_tdc_tdcelemID;   //!
   TBranch        *b_bb_grinch_tdc_tdcelemID;   //!
   TBranch        *b_Ndata_bb_grinch_tdc_tdclayer;   //!
   TBranch        *b_bb_grinch_tdc_tdclayer;   //!
   TBranch        *b_Ndata_bb_grinch_tdc_tdcrow;   //!
   TBranch        *b_bb_grinch_tdc_tdcrow;   //!
   TBranch        *b_Ndata_bb_hodoadc_bar_adc_L_a;   //!
   TBranch        *b_bb_hodoadc_bar_adc_L_a;   //!
   TBranch        *b_Ndata_bb_hodoadc_bar_adc_L_ac;   //!
   TBranch        *b_bb_hodoadc_bar_adc_L_ac;   //!
   TBranch        *b_Ndata_bb_hodoadc_bar_adc_L_ap;   //!
   TBranch        *b_bb_hodoadc_bar_adc_L_ap;   //!
   TBranch        *b_Ndata_bb_hodoadc_bar_adc_R_a;   //!
   TBranch        *b_bb_hodoadc_bar_adc_R_a;   //!
   TBranch        *b_Ndata_bb_hodoadc_bar_adc_R_ac;   //!
   TBranch        *b_bb_hodoadc_bar_adc_R_ac;   //!
   TBranch        *b_Ndata_bb_hodoadc_bar_adc_R_ap;   //!
   TBranch        *b_bb_hodoadc_bar_adc_R_ap;   //!
   TBranch        *b_Ndata_bb_hodoadc_bar_adc_id;   //!
   TBranch        *b_bb_hodoadc_bar_adc_id;   //!
   TBranch        *b_Ndata_bb_hodoadc_bar_adc_mean;   //!
   TBranch        *b_bb_hodoadc_bar_adc_mean;   //!
   TBranch        *b_Ndata_bb_hodotdc_bar_tdc_L_le;   //!
   TBranch        *b_bb_hodotdc_bar_tdc_L_le;   //!
   TBranch        *b_Ndata_bb_hodotdc_bar_tdc_L_leW;   //!
   TBranch        *b_bb_hodotdc_bar_tdc_L_leW;   //!
   TBranch        *b_Ndata_bb_hodotdc_bar_tdc_L_te;   //!
   TBranch        *b_bb_hodotdc_bar_tdc_L_te;   //!
   TBranch        *b_Ndata_bb_hodotdc_bar_tdc_L_teW;   //!
   TBranch        *b_bb_hodotdc_bar_tdc_L_teW;   //!
   TBranch        *b_Ndata_bb_hodotdc_bar_tdc_L_tot;   //!
   TBranch        *b_bb_hodotdc_bar_tdc_L_tot;   //!
   TBranch        *b_Ndata_bb_hodotdc_bar_tdc_L_totW;   //!
   TBranch        *b_bb_hodotdc_bar_tdc_L_totW;   //!
   TBranch        *b_Ndata_bb_hodotdc_bar_tdc_R_le;   //!
   TBranch        *b_bb_hodotdc_bar_tdc_R_le;   //!
   TBranch        *b_Ndata_bb_hodotdc_bar_tdc_R_leW;   //!
   TBranch        *b_bb_hodotdc_bar_tdc_R_leW;   //!
   TBranch        *b_Ndata_bb_hodotdc_bar_tdc_R_te;   //!
   TBranch        *b_bb_hodotdc_bar_tdc_R_te;   //!
   TBranch        *b_Ndata_bb_hodotdc_bar_tdc_R_teW;   //!
   TBranch        *b_bb_hodotdc_bar_tdc_R_teW;   //!
   TBranch        *b_Ndata_bb_hodotdc_bar_tdc_R_tot;   //!
   TBranch        *b_bb_hodotdc_bar_tdc_R_tot;   //!
   TBranch        *b_Ndata_bb_hodotdc_bar_tdc_R_totW;   //!
   TBranch        *b_bb_hodotdc_bar_tdc_R_totW;   //!
   TBranch        *b_Ndata_bb_hodotdc_bar_tdc_id;   //!
   TBranch        *b_bb_hodotdc_bar_tdc_id;   //!
   TBranch        *b_Ndata_bb_hodotdc_bar_tdc_meantime;   //!
   TBranch        *b_bb_hodotdc_bar_tdc_meantime;   //!
   TBranch        *b_Ndata_bb_hodotdc_bar_tdc_timediff;   //!
   TBranch        *b_bb_hodotdc_bar_tdc_timediff;   //!
   TBranch        *b_Ndata_bb_hodotdc_bar_tdc_timehitpos;   //!
   TBranch        *b_bb_hodotdc_bar_tdc_timehitpos;   //!
   TBranch        *b_Ndata_bb_hodotdc_bar_tdc_vpos;   //!
   TBranch        *b_bb_hodotdc_bar_tdc_vpos;   //!
   TBranch        *b_Ndata_bb_ps_clus_col;   //!
   TBranch        *b_bb_ps_clus_col;   //!
   TBranch        *b_Ndata_bb_ps_clus_e;   //!
   TBranch        *b_bb_ps_clus_e;   //!
   TBranch        *b_Ndata_bb_ps_clus_e_c;   //!
   TBranch        *b_bb_ps_clus_e_c;   //!
   TBranch        *b_Ndata_bb_ps_clus_eblk;   //!
   TBranch        *b_bb_ps_clus_eblk;   //!
   TBranch        *b_Ndata_bb_ps_clus_eblk_c;   //!
   TBranch        *b_bb_ps_clus_eblk_c;   //!
   TBranch        *b_Ndata_bb_ps_clus_id;   //!
   TBranch        *b_bb_ps_clus_id;   //!
   TBranch        *b_Ndata_bb_ps_clus_nblk;   //!
   TBranch        *b_bb_ps_clus_nblk;   //!
   TBranch        *b_Ndata_bb_ps_clus_row;   //!
   TBranch        *b_bb_ps_clus_row;   //!
   TBranch        *b_Ndata_bb_ps_clus_x;   //!
   TBranch        *b_bb_ps_clus_x;   //!
   TBranch        *b_Ndata_bb_ps_clus_y;   //!
   TBranch        *b_bb_ps_clus_y;   //!
   TBranch        *b_Ndata_bb_ps_clus_blk_col;   //!
   TBranch        *b_bb_ps_clus_blk_col;   //!
   TBranch        *b_Ndata_bb_ps_clus_blk_e;   //!
   TBranch        *b_bb_ps_clus_blk_e;   //!
   TBranch        *b_Ndata_bb_ps_clus_blk_e_c;   //!
   TBranch        *b_bb_ps_clus_blk_e_c;   //!
   TBranch        *b_Ndata_bb_ps_clus_blk_id;   //!
   TBranch        *b_bb_ps_clus_blk_id;   //!
   TBranch        *b_Ndata_bb_ps_clus_blk_row;   //!
   TBranch        *b_bb_ps_clus_blk_row;   //!
   TBranch        *b_Ndata_bb_ps_clus_blk_x;   //!
   TBranch        *b_bb_ps_clus_blk_x;   //!
   TBranch        *b_Ndata_bb_ps_clus_blk_y;   //!
   TBranch        *b_bb_ps_clus_blk_y;   //!
   TBranch        *b_Ndata_bb_sh_clus_col;   //!
   TBranch        *b_bb_sh_clus_col;   //!
   TBranch        *b_Ndata_bb_sh_clus_e;   //!
   TBranch        *b_bb_sh_clus_e;   //!
   TBranch        *b_Ndata_bb_sh_clus_e_c;   //!
   TBranch        *b_bb_sh_clus_e_c;   //!
   TBranch        *b_Ndata_bb_sh_clus_eblk;   //!
   TBranch        *b_bb_sh_clus_eblk;   //!
   TBranch        *b_Ndata_bb_sh_clus_eblk_c;   //!
   TBranch        *b_bb_sh_clus_eblk_c;   //!
   TBranch        *b_Ndata_bb_sh_clus_id;   //!
   TBranch        *b_bb_sh_clus_id;   //!
   TBranch        *b_Ndata_bb_sh_clus_nblk;   //!
   TBranch        *b_bb_sh_clus_nblk;   //!
   TBranch        *b_Ndata_bb_sh_clus_row;   //!
   TBranch        *b_bb_sh_clus_row;   //!
   TBranch        *b_Ndata_bb_sh_clus_x;   //!
   TBranch        *b_bb_sh_clus_x;   //!
   TBranch        *b_Ndata_bb_sh_clus_y;   //!
   TBranch        *b_bb_sh_clus_y;   //!
   TBranch        *b_Ndata_bb_sh_clus_blk_col;   //!
   TBranch        *b_bb_sh_clus_blk_col;   //!
   TBranch        *b_Ndata_bb_sh_clus_blk_e;   //!
   TBranch        *b_bb_sh_clus_blk_e;   //!
   TBranch        *b_Ndata_bb_sh_clus_blk_e_c;   //!
   TBranch        *b_bb_sh_clus_blk_e_c;   //!
   TBranch        *b_Ndata_bb_sh_clus_blk_id;   //!
   TBranch        *b_bb_sh_clus_blk_id;   //!
   TBranch        *b_Ndata_bb_sh_clus_blk_row;   //!
   TBranch        *b_bb_sh_clus_blk_row;   //!
   TBranch        *b_Ndata_bb_sh_clus_blk_x;   //!
   TBranch        *b_bb_sh_clus_blk_x;   //!
   TBranch        *b_Ndata_bb_sh_clus_blk_y;   //!
   TBranch        *b_bb_sh_clus_blk_y;   //!
   TBranch        *b_Ndata_bb_tr_beta;   //!
   TBranch        *b_bb_tr_beta;   //!
   TBranch        *b_Ndata_bb_tr_chi2;   //!
   TBranch        *b_bb_tr_chi2;   //!
   TBranch        *b_Ndata_bb_tr_d_ph;   //!
   TBranch        *b_bb_tr_d_ph;   //!
   TBranch        *b_Ndata_bb_tr_d_th;   //!
   TBranch        *b_bb_tr_d_th;   //!
   TBranch        *b_Ndata_bb_tr_d_x;   //!
   TBranch        *b_bb_tr_d_x;   //!
   TBranch        *b_Ndata_bb_tr_d_y;   //!
   TBranch        *b_bb_tr_d_y;   //!
   TBranch        *b_Ndata_bb_tr_dbeta;   //!
   TBranch        *b_bb_tr_dbeta;   //!
   TBranch        *b_Ndata_bb_tr_dtime;   //!
   TBranch        *b_bb_tr_dtime;   //!
   TBranch        *b_Ndata_bb_tr_flag;   //!
   TBranch        *b_bb_tr_flag;   //!
   TBranch        *b_Ndata_bb_tr_ndof;   //!
   TBranch        *b_bb_tr_ndof;   //!
   TBranch        *b_Ndata_bb_tr_p;   //!
   TBranch        *b_bb_tr_p;   //!
   TBranch        *b_Ndata_bb_tr_pathl;   //!
   TBranch        *b_bb_tr_pathl;   //!
   TBranch        *b_Ndata_bb_tr_ph;   //!
   TBranch        *b_bb_tr_ph;   //!
   TBranch        *b_Ndata_bb_tr_px;   //!
   TBranch        *b_bb_tr_px;   //!
   TBranch        *b_Ndata_bb_tr_py;   //!
   TBranch        *b_bb_tr_py;   //!
   TBranch        *b_Ndata_bb_tr_pz;   //!
   TBranch        *b_bb_tr_pz;   //!
   TBranch        *b_Ndata_bb_tr_r_ph;   //!
   TBranch        *b_bb_tr_r_ph;   //!
   TBranch        *b_Ndata_bb_tr_r_th;   //!
   TBranch        *b_bb_tr_r_th;   //!
   TBranch        *b_Ndata_bb_tr_r_x;   //!
   TBranch        *b_bb_tr_r_x;   //!
   TBranch        *b_Ndata_bb_tr_r_y;   //!
   TBranch        *b_bb_tr_r_y;   //!
   TBranch        *b_Ndata_bb_tr_tg_dp;   //!
   TBranch        *b_bb_tr_tg_dp;   //!
   TBranch        *b_Ndata_bb_tr_tg_ph;   //!
   TBranch        *b_bb_tr_tg_ph;   //!
   TBranch        *b_Ndata_bb_tr_tg_th;   //!
   TBranch        *b_bb_tr_tg_th;   //!
   TBranch        *b_Ndata_bb_tr_tg_x;   //!
   TBranch        *b_bb_tr_tg_x;   //!
   TBranch        *b_Ndata_bb_tr_tg_y;   //!
   TBranch        *b_bb_tr_tg_y;   //!
   TBranch        *b_Ndata_bb_tr_th;   //!
   TBranch        *b_bb_tr_th;   //!
   TBranch        *b_Ndata_bb_tr_time;   //!
   TBranch        *b_bb_tr_time;   //!
   TBranch        *b_Ndata_bb_tr_vx;   //!
   TBranch        *b_bb_tr_vx;   //!
   TBranch        *b_Ndata_bb_tr_vy;   //!
   TBranch        *b_bb_tr_vy;   //!
   TBranch        *b_Ndata_bb_tr_vz;   //!
   TBranch        *b_bb_tr_vz;   //!
   TBranch        *b_Ndata_bb_tr_x;   //!
   TBranch        *b_bb_tr_x;   //!
   TBranch        *b_Ndata_bb_tr_y;   //!
   TBranch        *b_bb_tr_y;   //!
   TBranch        *b_Ndata_bb_x_bcp;   //!
   TBranch        *b_bb_x_bcp;   //!
   TBranch        *b_Ndata_bb_x_fcp;   //!
   TBranch        *b_bb_x_fcp;   //!
   TBranch        *b_Ndata_bb_y_bcp;   //!
   TBranch        *b_bb_y_bcp;   //!
   TBranch        *b_Ndata_bb_y_fcp;   //!
   TBranch        *b_bb_y_fcp;   //!
   TBranch        *b_Ndata_sbs_hcal_a_amp;   //!
   TBranch        *b_sbs_hcal_a_amp;   //!
   TBranch        *b_Ndata_sbs_hcal_a_c;   //!
   TBranch        *b_sbs_hcal_a_c;   //!
   TBranch        *b_Ndata_sbs_hcal_a_p;   //!
   TBranch        *b_sbs_hcal_a_p;   //!
   TBranch        *b_Ndata_sbs_hcal_a_time;   //!
   TBranch        *b_sbs_hcal_a_time;   //!
   TBranch        *b_Ndata_sbs_hcal_adccol;   //!
   TBranch        *b_sbs_hcal_adccol;   //!
   TBranch        *b_Ndata_sbs_hcal_adcrow;   //!
   TBranch        *b_sbs_hcal_adcrow;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_col;   //!
   TBranch        *b_sbs_hcal_clus_col;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_e;   //!
   TBranch        *b_sbs_hcal_clus_e;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_eblk;   //!
   TBranch        *b_sbs_hcal_clus_eblk;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_id;   //!
   TBranch        *b_sbs_hcal_clus_id;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_nblk;   //!
   TBranch        *b_sbs_hcal_clus_nblk;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_row;   //!
   TBranch        *b_sbs_hcal_clus_row;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_blk_e;   //!
   TBranch        *b_sbs_hcal_clus_blk_e;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_blk_id;   //!
   TBranch        *b_sbs_hcal_clus_blk_id;   //!
   TBranch        *b_BB_gold_beta;   //!
   TBranch        *b_BB_gold_dp;   //!
   TBranch        *b_BB_gold_index;   //!
   TBranch        *b_BB_gold_ok;   //!
   TBranch        *b_BB_gold_p;   //!
   TBranch        *b_BB_gold_ph;   //!
   TBranch        *b_BB_gold_px;   //!
   TBranch        *b_BB_gold_py;   //!
   TBranch        *b_BB_gold_pz;   //!
   TBranch        *b_BB_gold_th;   //!
   TBranch        *b_BB_gold_x;   //!
   TBranch        *b_BB_gold_y;   //!
   TBranch        *b_bb_gem_hit_ngoodhits;   //!
   TBranch        *b_bb_gem_nlayershit;   //!
   TBranch        *b_bb_gem_nlayershitu;   //!
   TBranch        *b_bb_gem_nlayershituv;   //!
   TBranch        *b_bb_gem_nlayershitv;   //!
   TBranch        *b_bb_gem_track_besttrack;   //!
   TBranch        *b_bb_gem_track_ntrack;   //!
   TBranch        *b_bb_grinch_tdc_ngoodADChits;   //!
   TBranch        *b_bb_grinch_tdc_ngoodTDChits;   //!
   TBranch        *b_bb_grinch_tdc_nhits;   //!
   TBranch        *b_bb_grinch_tdc_nrefhits;   //!
   TBranch        *b_bb_hodotdc_bar_ngoodbars;   //!
   TBranch        *b_bb_hodotdc_nclus;   //!
   TBranch        *b_bb_ps_e;   //!
   TBranch        *b_bb_ps_e_c;   //!
   TBranch        *b_bb_ps_nblk;   //!
   TBranch        *b_bb_ps_nclus;   //!
   TBranch        *b_bb_ps_ngoodADChits;   //!
   TBranch        *b_bb_ps_x;   //!
   TBranch        *b_bb_ps_y;   //!
   TBranch        *b_bb_sh_e;   //!
   TBranch        *b_bb_sh_e_c;   //!
   TBranch        *b_bb_sh_nblk;   //!
   TBranch        *b_bb_sh_nclus;   //!
   TBranch        *b_bb_sh_ngoodADChits;   //!
   TBranch        *b_bb_sh_x;   //!
   TBranch        *b_bb_sh_y;   //!
   TBranch        *b_bb_tr_n;   //!
   TBranch        *b_sbs_hcal_e;   //!
   TBranch        *b_sbs_hcal_nblk;   //!
   TBranch        *b_sbs_hcal_nclus;   //!
   TBranch        *b_sbs_hcal_ngoodADChits;   //!
   TBranch        *b_singletrack;   //!
   TBranch        *b_anytrack;   //!
   TBranch        *b_Event_Branch_fEvtHdr_fEvtTime;   //!
   TBranch        *b_Event_Branch_fEvtHdr_fEvtNum;   //!
   TBranch        *b_Event_Branch_fEvtHdr_fEvtType;   //!
   TBranch        *b_Event_Branch_fEvtHdr_fEvtLen;   //!
   TBranch        *b_Event_Branch_fEvtHdr_fHelicity;   //!
   TBranch        *b_Event_Branch_fEvtHdr_fTargetPol;   //!
   TBranch        *b_Event_Branch_fEvtHdr_fRun;   //!

   gmn_rec_tree(TTree *tree=0);
   virtual ~gmn_rec_tree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef gmn_rec_tree_cxx
gmn_rec_tree::gmn_rec_tree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("gmn_replayed_11251_stream0_seg0_0.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("gmn_replayed_11251_stream0_seg0_0.root");
      }
      f->GetObject("T",tree);

   }
   Init(tree);
}

gmn_rec_tree::~gmn_rec_tree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t gmn_rec_tree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t gmn_rec_tree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void gmn_rec_tree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Ndata.bb.eps_over_etot", &Ndata_bb_eps_over_etot, &b_Ndata_bb_eps_over_etot);
   fChain->SetBranchAddress("bb.eps_over_etot", bb_eps_over_etot, &b_bb_eps_over_etot);
   fChain->SetBranchAddress("Ndata.bb.etot_over_p", &Ndata_bb_etot_over_p, &b_Ndata_bb_etot_over_p);
   fChain->SetBranchAddress("bb.etot_over_p", bb_etot_over_p, &b_bb_etot_over_p);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ADCU", &Ndata_bb_gem_hit_ADCU, &b_Ndata_bb_gem_hit_ADCU);
   fChain->SetBranchAddress("bb.gem.hit.ADCU", bb_gem_hit_ADCU, &b_bb_gem_hit_ADCU);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ADCV", &Ndata_bb_gem_hit_ADCV, &b_Ndata_bb_gem_hit_ADCV);
   fChain->SetBranchAddress("bb.gem.hit.ADCV", bb_gem_hit_ADCV, &b_bb_gem_hit_ADCV);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ADCasym", &Ndata_bb_gem_hit_ADCasym, &b_Ndata_bb_gem_hit_ADCasym);
   fChain->SetBranchAddress("bb.gem.hit.ADCasym", bb_gem_hit_ADCasym, &b_bb_gem_hit_ADCasym);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ADCavg", &Ndata_bb_gem_hit_ADCavg, &b_Ndata_bb_gem_hit_ADCavg);
   fChain->SetBranchAddress("bb.gem.hit.ADCavg", bb_gem_hit_ADCavg, &b_bb_gem_hit_ADCavg);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ADCmaxsampU", &Ndata_bb_gem_hit_ADCmaxsampU, &b_Ndata_bb_gem_hit_ADCmaxsampU);
   fChain->SetBranchAddress("bb.gem.hit.ADCmaxsampU", bb_gem_hit_ADCmaxsampU, &b_bb_gem_hit_ADCmaxsampU);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ADCmaxsampUclust", &Ndata_bb_gem_hit_ADCmaxsampUclust, &b_Ndata_bb_gem_hit_ADCmaxsampUclust);
   fChain->SetBranchAddress("bb.gem.hit.ADCmaxsampUclust", bb_gem_hit_ADCmaxsampUclust, &b_bb_gem_hit_ADCmaxsampUclust);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ADCmaxsampV", &Ndata_bb_gem_hit_ADCmaxsampV, &b_Ndata_bb_gem_hit_ADCmaxsampV);
   fChain->SetBranchAddress("bb.gem.hit.ADCmaxsampV", bb_gem_hit_ADCmaxsampV, &b_bb_gem_hit_ADCmaxsampV);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ADCmaxsampVclust", &Ndata_bb_gem_hit_ADCmaxsampVclust, &b_Ndata_bb_gem_hit_ADCmaxsampVclust);
   fChain->SetBranchAddress("bb.gem.hit.ADCmaxsampVclust", bb_gem_hit_ADCmaxsampVclust, &b_bb_gem_hit_ADCmaxsampVclust);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ADCmaxstripU", &Ndata_bb_gem_hit_ADCmaxstripU, &b_Ndata_bb_gem_hit_ADCmaxstripU);
   fChain->SetBranchAddress("bb.gem.hit.ADCmaxstripU", bb_gem_hit_ADCmaxstripU, &b_bb_gem_hit_ADCmaxstripU);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ADCmaxstripV", &Ndata_bb_gem_hit_ADCmaxstripV, &b_Ndata_bb_gem_hit_ADCmaxstripV);
   fChain->SetBranchAddress("bb.gem.hit.ADCmaxstripV", bb_gem_hit_ADCmaxstripV, &b_bb_gem_hit_ADCmaxstripV);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.BUILD_ALL_SAMPLES_U", &Ndata_bb_gem_hit_BUILD_ALL_SAMPLES_U, &b_Ndata_bb_gem_hit_BUILD_ALL_SAMPLES_U);
   fChain->SetBranchAddress("bb.gem.hit.BUILD_ALL_SAMPLES_U", bb_gem_hit_BUILD_ALL_SAMPLES_U, &b_bb_gem_hit_BUILD_ALL_SAMPLES_U);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.BUILD_ALL_SAMPLES_V", &Ndata_bb_gem_hit_BUILD_ALL_SAMPLES_V, &b_Ndata_bb_gem_hit_BUILD_ALL_SAMPLES_V);
   fChain->SetBranchAddress("bb.gem.hit.BUILD_ALL_SAMPLES_V", bb_gem_hit_BUILD_ALL_SAMPLES_V, &b_bb_gem_hit_BUILD_ALL_SAMPLES_V);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.CM_GOOD_U", &Ndata_bb_gem_hit_CM_GOOD_U, &b_Ndata_bb_gem_hit_CM_GOOD_U);
   fChain->SetBranchAddress("bb.gem.hit.CM_GOOD_U", bb_gem_hit_CM_GOOD_U, &b_bb_gem_hit_CM_GOOD_U);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.CM_GOOD_V", &Ndata_bb_gem_hit_CM_GOOD_V, &b_Ndata_bb_gem_hit_CM_GOOD_V);
   fChain->SetBranchAddress("bb.gem.hit.CM_GOOD_V", bb_gem_hit_CM_GOOD_V, &b_bb_gem_hit_CM_GOOD_V);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ENABLE_CM_U", &Ndata_bb_gem_hit_ENABLE_CM_U, &b_Ndata_bb_gem_hit_ENABLE_CM_U);
   fChain->SetBranchAddress("bb.gem.hit.ENABLE_CM_U", bb_gem_hit_ENABLE_CM_U, &b_bb_gem_hit_ENABLE_CM_U);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ENABLE_CM_V", &Ndata_bb_gem_hit_ENABLE_CM_V, &b_Ndata_bb_gem_hit_ENABLE_CM_V);
   fChain->SetBranchAddress("bb.gem.hit.ENABLE_CM_V", bb_gem_hit_ENABLE_CM_V, &b_bb_gem_hit_ENABLE_CM_V);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.Tavg", &Ndata_bb_gem_hit_Tavg, &b_Ndata_bb_gem_hit_Tavg);
   fChain->SetBranchAddress("bb.gem.hit.Tavg", bb_gem_hit_Tavg, &b_bb_gem_hit_Tavg);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.Utime", &Ndata_bb_gem_hit_Utime, &b_Ndata_bb_gem_hit_Utime);
   fChain->SetBranchAddress("bb.gem.hit.Utime", bb_gem_hit_Utime, &b_bb_gem_hit_Utime);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.UtimeMaxStrip", &Ndata_bb_gem_hit_UtimeMaxStrip, &b_Ndata_bb_gem_hit_UtimeMaxStrip);
   fChain->SetBranchAddress("bb.gem.hit.UtimeMaxStrip", bb_gem_hit_UtimeMaxStrip, &b_bb_gem_hit_UtimeMaxStrip);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.Vtime", &Ndata_bb_gem_hit_Vtime, &b_Ndata_bb_gem_hit_Vtime);
   fChain->SetBranchAddress("bb.gem.hit.Vtime", bb_gem_hit_Vtime, &b_bb_gem_hit_Vtime);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.VtimeMaxStrip", &Ndata_bb_gem_hit_VtimeMaxStrip, &b_Ndata_bb_gem_hit_VtimeMaxStrip);
   fChain->SetBranchAddress("bb.gem.hit.VtimeMaxStrip", bb_gem_hit_VtimeMaxStrip, &b_bb_gem_hit_VtimeMaxStrip);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ccor_clust", &Ndata_bb_gem_hit_ccor_clust, &b_Ndata_bb_gem_hit_ccor_clust);
   fChain->SetBranchAddress("bb.gem.hit.ccor_clust", bb_gem_hit_ccor_clust, &b_bb_gem_hit_ccor_clust);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ccor_strip", &Ndata_bb_gem_hit_ccor_strip, &b_Ndata_bb_gem_hit_ccor_strip);
   fChain->SetBranchAddress("bb.gem.hit.ccor_strip", bb_gem_hit_ccor_strip, &b_bb_gem_hit_ccor_strip);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.deltat", &Ndata_bb_gem_hit_deltat, &b_Ndata_bb_gem_hit_deltat);
   fChain->SetBranchAddress("bb.gem.hit.deltat", bb_gem_hit_deltat, &b_bb_gem_hit_deltat);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.eresidu", &Ndata_bb_gem_hit_eresidu, &b_Ndata_bb_gem_hit_eresidu);
   fChain->SetBranchAddress("bb.gem.hit.eresidu", bb_gem_hit_eresidu, &b_bb_gem_hit_eresidu);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.eresidv", &Ndata_bb_gem_hit_eresidv, &b_Ndata_bb_gem_hit_eresidv);
   fChain->SetBranchAddress("bb.gem.hit.eresidv", bb_gem_hit_eresidv, &b_bb_gem_hit_eresidv);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.isampmaxUclust", &Ndata_bb_gem_hit_isampmaxUclust, &b_Ndata_bb_gem_hit_isampmaxUclust);
   fChain->SetBranchAddress("bb.gem.hit.isampmaxUclust", bb_gem_hit_isampmaxUclust, &b_bb_gem_hit_isampmaxUclust);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.isampmaxUstrip", &Ndata_bb_gem_hit_isampmaxUstrip, &b_Ndata_bb_gem_hit_isampmaxUstrip);
   fChain->SetBranchAddress("bb.gem.hit.isampmaxUstrip", bb_gem_hit_isampmaxUstrip, &b_bb_gem_hit_isampmaxUstrip);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.isampmaxVclust", &Ndata_bb_gem_hit_isampmaxVclust, &b_Ndata_bb_gem_hit_isampmaxVclust);
   fChain->SetBranchAddress("bb.gem.hit.isampmaxVclust", bb_gem_hit_isampmaxVclust, &b_bb_gem_hit_isampmaxVclust);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.isampmaxVstrip", &Ndata_bb_gem_hit_isampmaxVstrip, &b_Ndata_bb_gem_hit_isampmaxVstrip);
   fChain->SetBranchAddress("bb.gem.hit.isampmaxVstrip", bb_gem_hit_isampmaxVstrip, &b_bb_gem_hit_isampmaxVstrip);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.layer", &Ndata_bb_gem_hit_layer, &b_Ndata_bb_gem_hit_layer);
   fChain->SetBranchAddress("bb.gem.hit.layer", bb_gem_hit_layer, &b_bb_gem_hit_layer);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.module", &Ndata_bb_gem_hit_module, &b_Ndata_bb_gem_hit_module);
   fChain->SetBranchAddress("bb.gem.hit.module", bb_gem_hit_module, &b_bb_gem_hit_module);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.nstripu", &Ndata_bb_gem_hit_nstripu, &b_Ndata_bb_gem_hit_nstripu);
   fChain->SetBranchAddress("bb.gem.hit.nstripu", bb_gem_hit_nstripu, &b_bb_gem_hit_nstripu);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.nstripv", &Ndata_bb_gem_hit_nstripv, &b_Ndata_bb_gem_hit_nstripv);
   fChain->SetBranchAddress("bb.gem.hit.nstripv", bb_gem_hit_nstripv, &b_bb_gem_hit_nstripv);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.residu", &Ndata_bb_gem_hit_residu, &b_Ndata_bb_gem_hit_residu);
   fChain->SetBranchAddress("bb.gem.hit.residu", bb_gem_hit_residu, &b_bb_gem_hit_residu);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.residv", &Ndata_bb_gem_hit_residv, &b_Ndata_bb_gem_hit_residv);
   fChain->SetBranchAddress("bb.gem.hit.residv", bb_gem_hit_residv, &b_bb_gem_hit_residv);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.trackindex", &Ndata_bb_gem_hit_trackindex, &b_Ndata_bb_gem_hit_trackindex);
   fChain->SetBranchAddress("bb.gem.hit.trackindex", bb_gem_hit_trackindex, &b_bb_gem_hit_trackindex);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.u", &Ndata_bb_gem_hit_u, &b_Ndata_bb_gem_hit_u);
   fChain->SetBranchAddress("bb.gem.hit.u", bb_gem_hit_u, &b_bb_gem_hit_u);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.umoment", &Ndata_bb_gem_hit_umoment, &b_Ndata_bb_gem_hit_umoment);
   fChain->SetBranchAddress("bb.gem.hit.umoment", bb_gem_hit_umoment, &b_bb_gem_hit_umoment);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.usigma", &Ndata_bb_gem_hit_usigma, &b_Ndata_bb_gem_hit_usigma);
   fChain->SetBranchAddress("bb.gem.hit.usigma", bb_gem_hit_usigma, &b_bb_gem_hit_usigma);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ustriphi", &Ndata_bb_gem_hit_ustriphi, &b_Ndata_bb_gem_hit_ustriphi);
   fChain->SetBranchAddress("bb.gem.hit.ustriphi", bb_gem_hit_ustriphi, &b_bb_gem_hit_ustriphi);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ustriplo", &Ndata_bb_gem_hit_ustriplo, &b_Ndata_bb_gem_hit_ustriplo);
   fChain->SetBranchAddress("bb.gem.hit.ustriplo", bb_gem_hit_ustriplo, &b_bb_gem_hit_ustriplo);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ustripmax", &Ndata_bb_gem_hit_ustripmax, &b_Ndata_bb_gem_hit_ustripmax);
   fChain->SetBranchAddress("bb.gem.hit.ustripmax", bb_gem_hit_ustripmax, &b_bb_gem_hit_ustripmax);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.v", &Ndata_bb_gem_hit_v, &b_Ndata_bb_gem_hit_v);
   fChain->SetBranchAddress("bb.gem.hit.v", bb_gem_hit_v, &b_bb_gem_hit_v);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.vmoment", &Ndata_bb_gem_hit_vmoment, &b_Ndata_bb_gem_hit_vmoment);
   fChain->SetBranchAddress("bb.gem.hit.vmoment", bb_gem_hit_vmoment, &b_bb_gem_hit_vmoment);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.vsigma", &Ndata_bb_gem_hit_vsigma, &b_Ndata_bb_gem_hit_vsigma);
   fChain->SetBranchAddress("bb.gem.hit.vsigma", bb_gem_hit_vsigma, &b_bb_gem_hit_vsigma);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.vstriphi", &Ndata_bb_gem_hit_vstriphi, &b_Ndata_bb_gem_hit_vstriphi);
   fChain->SetBranchAddress("bb.gem.hit.vstriphi", bb_gem_hit_vstriphi, &b_bb_gem_hit_vstriphi);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.vstriplo", &Ndata_bb_gem_hit_vstriplo, &b_Ndata_bb_gem_hit_vstriplo);
   fChain->SetBranchAddress("bb.gem.hit.vstriplo", bb_gem_hit_vstriplo, &b_bb_gem_hit_vstriplo);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.vstripmax", &Ndata_bb_gem_hit_vstripmax, &b_Ndata_bb_gem_hit_vstripmax);
   fChain->SetBranchAddress("bb.gem.hit.vstripmax", bb_gem_hit_vstripmax, &b_bb_gem_hit_vstripmax);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.xglobal", &Ndata_bb_gem_hit_xglobal, &b_Ndata_bb_gem_hit_xglobal);
   fChain->SetBranchAddress("bb.gem.hit.xglobal", bb_gem_hit_xglobal, &b_bb_gem_hit_xglobal);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.xlocal", &Ndata_bb_gem_hit_xlocal, &b_Ndata_bb_gem_hit_xlocal);
   fChain->SetBranchAddress("bb.gem.hit.xlocal", bb_gem_hit_xlocal, &b_bb_gem_hit_xlocal);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.yglobal", &Ndata_bb_gem_hit_yglobal, &b_Ndata_bb_gem_hit_yglobal);
   fChain->SetBranchAddress("bb.gem.hit.yglobal", bb_gem_hit_yglobal, &b_bb_gem_hit_yglobal);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ylocal", &Ndata_bb_gem_hit_ylocal, &b_Ndata_bb_gem_hit_ylocal);
   fChain->SetBranchAddress("bb.gem.hit.ylocal", bb_gem_hit_ylocal, &b_bb_gem_hit_ylocal);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.zglobal", &Ndata_bb_gem_hit_zglobal, &b_Ndata_bb_gem_hit_zglobal);
   fChain->SetBranchAddress("bb.gem.hit.zglobal", bb_gem_hit_zglobal, &b_bb_gem_hit_zglobal);
   fChain->SetBranchAddress("Ndata.bb.gem.n2Dhit_layer", &Ndata_bb_gem_n2Dhit_layer, &b_Ndata_bb_gem_n2Dhit_layer);
   fChain->SetBranchAddress("bb.gem.n2Dhit_layer", bb_gem_n2Dhit_layer, &b_bb_gem_n2Dhit_layer);
   fChain->SetBranchAddress("Ndata.bb.gem.nclustu_layer", &Ndata_bb_gem_nclustu_layer, &b_Ndata_bb_gem_nclustu_layer);
   fChain->SetBranchAddress("bb.gem.nclustu_layer", bb_gem_nclustu_layer, &b_bb_gem_nclustu_layer);
   fChain->SetBranchAddress("Ndata.bb.gem.nclustv_layer", &Ndata_bb_gem_nclustv_layer, &b_Ndata_bb_gem_nclustv_layer);
   fChain->SetBranchAddress("bb.gem.nclustv_layer", bb_gem_nclustv_layer, &b_bb_gem_nclustv_layer);
   fChain->SetBranchAddress("Ndata.bb.gem.nstripsu_layer", &Ndata_bb_gem_nstripsu_layer, &b_Ndata_bb_gem_nstripsu_layer);
   fChain->SetBranchAddress("bb.gem.nstripsu_layer", bb_gem_nstripsu_layer, &b_bb_gem_nstripsu_layer);
   fChain->SetBranchAddress("Ndata.bb.gem.nstripsv_layer", &Ndata_bb_gem_nstripsv_layer, &b_Ndata_bb_gem_nstripsv_layer);
   fChain->SetBranchAddress("bb.gem.nstripsv_layer", bb_gem_nstripsv_layer, &b_bb_gem_nstripsv_layer);
   fChain->SetBranchAddress("Ndata.bb.gem.track.chi2ndf", &Ndata_bb_gem_track_chi2ndf, &b_Ndata_bb_gem_track_chi2ndf);
   fChain->SetBranchAddress("bb.gem.track.chi2ndf", bb_gem_track_chi2ndf, &b_bb_gem_track_chi2ndf);
   fChain->SetBranchAddress("Ndata.bb.gem.track.nhits", &Ndata_bb_gem_track_nhits, &b_Ndata_bb_gem_track_nhits);
   fChain->SetBranchAddress("bb.gem.track.nhits", bb_gem_track_nhits, &b_bb_gem_track_nhits);
   fChain->SetBranchAddress("Ndata.bb.gem.track.x", &Ndata_bb_gem_track_x, &b_Ndata_bb_gem_track_x);
   fChain->SetBranchAddress("bb.gem.track.x", bb_gem_track_x, &b_bb_gem_track_x);
   fChain->SetBranchAddress("Ndata.bb.gem.track.xp", &Ndata_bb_gem_track_xp, &b_Ndata_bb_gem_track_xp);
   fChain->SetBranchAddress("bb.gem.track.xp", bb_gem_track_xp, &b_bb_gem_track_xp);
   fChain->SetBranchAddress("Ndata.bb.gem.track.y", &Ndata_bb_gem_track_y, &b_Ndata_bb_gem_track_y);
   fChain->SetBranchAddress("bb.gem.track.y", bb_gem_track_y, &b_bb_gem_track_y);
   fChain->SetBranchAddress("Ndata.bb.gem.track.yp", &Ndata_bb_gem_track_yp, &b_Ndata_bb_gem_track_yp);
   fChain->SetBranchAddress("bb.gem.track.yp", bb_gem_track_yp, &b_bb_gem_track_yp);
   fChain->SetBranchAddress("Ndata.bb.grinch_tdc.tdc", &Ndata_bb_grinch_tdc_tdc, &b_Ndata_bb_grinch_tdc_tdc);
   fChain->SetBranchAddress("bb.grinch_tdc.tdc", bb_grinch_tdc_tdc, &b_bb_grinch_tdc_tdc);
   fChain->SetBranchAddress("Ndata.bb.grinch_tdc.tdc_mult", &Ndata_bb_grinch_tdc_tdc_mult, &b_Ndata_bb_grinch_tdc_tdc_mult);
   fChain->SetBranchAddress("bb.grinch_tdc.tdc_mult", bb_grinch_tdc_tdc_mult, &b_bb_grinch_tdc_tdc_mult);
   fChain->SetBranchAddress("Ndata.bb.grinch_tdc.tdc_te", &Ndata_bb_grinch_tdc_tdc_te, &b_Ndata_bb_grinch_tdc_tdc_te);
   fChain->SetBranchAddress("bb.grinch_tdc.tdc_te", bb_grinch_tdc_tdc_te, &b_bb_grinch_tdc_tdc_te);
   fChain->SetBranchAddress("Ndata.bb.grinch_tdc.tdc_tot", &Ndata_bb_grinch_tdc_tdc_tot, &b_Ndata_bb_grinch_tdc_tdc_tot);
   fChain->SetBranchAddress("bb.grinch_tdc.tdc_tot", bb_grinch_tdc_tdc_tot, &b_bb_grinch_tdc_tdc_tot);
   fChain->SetBranchAddress("Ndata.bb.grinch_tdc.tdccol", &Ndata_bb_grinch_tdc_tdccol, &b_Ndata_bb_grinch_tdc_tdccol);
   fChain->SetBranchAddress("bb.grinch_tdc.tdccol", bb_grinch_tdc_tdccol, &b_bb_grinch_tdc_tdccol);
   fChain->SetBranchAddress("Ndata.bb.grinch_tdc.tdcelemID", &Ndata_bb_grinch_tdc_tdcelemID, &b_Ndata_bb_grinch_tdc_tdcelemID);
   fChain->SetBranchAddress("bb.grinch_tdc.tdcelemID", bb_grinch_tdc_tdcelemID, &b_bb_grinch_tdc_tdcelemID);
   fChain->SetBranchAddress("Ndata.bb.grinch_tdc.tdclayer", &Ndata_bb_grinch_tdc_tdclayer, &b_Ndata_bb_grinch_tdc_tdclayer);
   fChain->SetBranchAddress("bb.grinch_tdc.tdclayer", bb_grinch_tdc_tdclayer, &b_bb_grinch_tdc_tdclayer);
   fChain->SetBranchAddress("Ndata.bb.grinch_tdc.tdcrow", &Ndata_bb_grinch_tdc_tdcrow, &b_Ndata_bb_grinch_tdc_tdcrow);
   fChain->SetBranchAddress("bb.grinch_tdc.tdcrow", bb_grinch_tdc_tdcrow, &b_bb_grinch_tdc_tdcrow);
   fChain->SetBranchAddress("Ndata.bb.hodoadc.bar.adc.L.a", &Ndata_bb_hodoadc_bar_adc_L_a, &b_Ndata_bb_hodoadc_bar_adc_L_a);
   fChain->SetBranchAddress("bb.hodoadc.bar.adc.L.a", bb_hodoadc_bar_adc_L_a, &b_bb_hodoadc_bar_adc_L_a);
   fChain->SetBranchAddress("Ndata.bb.hodoadc.bar.adc.L.ac", &Ndata_bb_hodoadc_bar_adc_L_ac, &b_Ndata_bb_hodoadc_bar_adc_L_ac);
   fChain->SetBranchAddress("bb.hodoadc.bar.adc.L.ac", bb_hodoadc_bar_adc_L_ac, &b_bb_hodoadc_bar_adc_L_ac);
   fChain->SetBranchAddress("Ndata.bb.hodoadc.bar.adc.L.ap", &Ndata_bb_hodoadc_bar_adc_L_ap, &b_Ndata_bb_hodoadc_bar_adc_L_ap);
   fChain->SetBranchAddress("bb.hodoadc.bar.adc.L.ap", bb_hodoadc_bar_adc_L_ap, &b_bb_hodoadc_bar_adc_L_ap);
   fChain->SetBranchAddress("Ndata.bb.hodoadc.bar.adc.R.a", &Ndata_bb_hodoadc_bar_adc_R_a, &b_Ndata_bb_hodoadc_bar_adc_R_a);
   fChain->SetBranchAddress("bb.hodoadc.bar.adc.R.a", bb_hodoadc_bar_adc_R_a, &b_bb_hodoadc_bar_adc_R_a);
   fChain->SetBranchAddress("Ndata.bb.hodoadc.bar.adc.R.ac", &Ndata_bb_hodoadc_bar_adc_R_ac, &b_Ndata_bb_hodoadc_bar_adc_R_ac);
   fChain->SetBranchAddress("bb.hodoadc.bar.adc.R.ac", bb_hodoadc_bar_adc_R_ac, &b_bb_hodoadc_bar_adc_R_ac);
   fChain->SetBranchAddress("Ndata.bb.hodoadc.bar.adc.R.ap", &Ndata_bb_hodoadc_bar_adc_R_ap, &b_Ndata_bb_hodoadc_bar_adc_R_ap);
   fChain->SetBranchAddress("bb.hodoadc.bar.adc.R.ap", bb_hodoadc_bar_adc_R_ap, &b_bb_hodoadc_bar_adc_R_ap);
   fChain->SetBranchAddress("Ndata.bb.hodoadc.bar.adc.id", &Ndata_bb_hodoadc_bar_adc_id, &b_Ndata_bb_hodoadc_bar_adc_id);
   fChain->SetBranchAddress("bb.hodoadc.bar.adc.id", bb_hodoadc_bar_adc_id, &b_bb_hodoadc_bar_adc_id);
   fChain->SetBranchAddress("Ndata.bb.hodoadc.bar.adc.mean", &Ndata_bb_hodoadc_bar_adc_mean, &b_Ndata_bb_hodoadc_bar_adc_mean);
   fChain->SetBranchAddress("bb.hodoadc.bar.adc.mean", bb_hodoadc_bar_adc_mean, &b_bb_hodoadc_bar_adc_mean);
   fChain->SetBranchAddress("Ndata.bb.hodotdc.bar.tdc.L.le", &Ndata_bb_hodotdc_bar_tdc_L_le, &b_Ndata_bb_hodotdc_bar_tdc_L_le);
   fChain->SetBranchAddress("bb.hodotdc.bar.tdc.L.le", bb_hodotdc_bar_tdc_L_le, &b_bb_hodotdc_bar_tdc_L_le);
   fChain->SetBranchAddress("Ndata.bb.hodotdc.bar.tdc.L.leW", &Ndata_bb_hodotdc_bar_tdc_L_leW, &b_Ndata_bb_hodotdc_bar_tdc_L_leW);
   fChain->SetBranchAddress("bb.hodotdc.bar.tdc.L.leW", bb_hodotdc_bar_tdc_L_leW, &b_bb_hodotdc_bar_tdc_L_leW);
   fChain->SetBranchAddress("Ndata.bb.hodotdc.bar.tdc.L.te", &Ndata_bb_hodotdc_bar_tdc_L_te, &b_Ndata_bb_hodotdc_bar_tdc_L_te);
   fChain->SetBranchAddress("bb.hodotdc.bar.tdc.L.te", bb_hodotdc_bar_tdc_L_te, &b_bb_hodotdc_bar_tdc_L_te);
   fChain->SetBranchAddress("Ndata.bb.hodotdc.bar.tdc.L.teW", &Ndata_bb_hodotdc_bar_tdc_L_teW, &b_Ndata_bb_hodotdc_bar_tdc_L_teW);
   fChain->SetBranchAddress("bb.hodotdc.bar.tdc.L.teW", bb_hodotdc_bar_tdc_L_teW, &b_bb_hodotdc_bar_tdc_L_teW);
   fChain->SetBranchAddress("Ndata.bb.hodotdc.bar.tdc.L.tot", &Ndata_bb_hodotdc_bar_tdc_L_tot, &b_Ndata_bb_hodotdc_bar_tdc_L_tot);
   fChain->SetBranchAddress("bb.hodotdc.bar.tdc.L.tot", bb_hodotdc_bar_tdc_L_tot, &b_bb_hodotdc_bar_tdc_L_tot);
   fChain->SetBranchAddress("Ndata.bb.hodotdc.bar.tdc.L.totW", &Ndata_bb_hodotdc_bar_tdc_L_totW, &b_Ndata_bb_hodotdc_bar_tdc_L_totW);
   fChain->SetBranchAddress("bb.hodotdc.bar.tdc.L.totW", bb_hodotdc_bar_tdc_L_totW, &b_bb_hodotdc_bar_tdc_L_totW);
   fChain->SetBranchAddress("Ndata.bb.hodotdc.bar.tdc.R.le", &Ndata_bb_hodotdc_bar_tdc_R_le, &b_Ndata_bb_hodotdc_bar_tdc_R_le);
   fChain->SetBranchAddress("bb.hodotdc.bar.tdc.R.le", bb_hodotdc_bar_tdc_R_le, &b_bb_hodotdc_bar_tdc_R_le);
   fChain->SetBranchAddress("Ndata.bb.hodotdc.bar.tdc.R.leW", &Ndata_bb_hodotdc_bar_tdc_R_leW, &b_Ndata_bb_hodotdc_bar_tdc_R_leW);
   fChain->SetBranchAddress("bb.hodotdc.bar.tdc.R.leW", bb_hodotdc_bar_tdc_R_leW, &b_bb_hodotdc_bar_tdc_R_leW);
   fChain->SetBranchAddress("Ndata.bb.hodotdc.bar.tdc.R.te", &Ndata_bb_hodotdc_bar_tdc_R_te, &b_Ndata_bb_hodotdc_bar_tdc_R_te);
   fChain->SetBranchAddress("bb.hodotdc.bar.tdc.R.te", bb_hodotdc_bar_tdc_R_te, &b_bb_hodotdc_bar_tdc_R_te);
   fChain->SetBranchAddress("Ndata.bb.hodotdc.bar.tdc.R.teW", &Ndata_bb_hodotdc_bar_tdc_R_teW, &b_Ndata_bb_hodotdc_bar_tdc_R_teW);
   fChain->SetBranchAddress("bb.hodotdc.bar.tdc.R.teW", bb_hodotdc_bar_tdc_R_teW, &b_bb_hodotdc_bar_tdc_R_teW);
   fChain->SetBranchAddress("Ndata.bb.hodotdc.bar.tdc.R.tot", &Ndata_bb_hodotdc_bar_tdc_R_tot, &b_Ndata_bb_hodotdc_bar_tdc_R_tot);
   fChain->SetBranchAddress("bb.hodotdc.bar.tdc.R.tot", bb_hodotdc_bar_tdc_R_tot, &b_bb_hodotdc_bar_tdc_R_tot);
   fChain->SetBranchAddress("Ndata.bb.hodotdc.bar.tdc.R.totW", &Ndata_bb_hodotdc_bar_tdc_R_totW, &b_Ndata_bb_hodotdc_bar_tdc_R_totW);
   fChain->SetBranchAddress("bb.hodotdc.bar.tdc.R.totW", bb_hodotdc_bar_tdc_R_totW, &b_bb_hodotdc_bar_tdc_R_totW);
   fChain->SetBranchAddress("Ndata.bb.hodotdc.bar.tdc.id", &Ndata_bb_hodotdc_bar_tdc_id, &b_Ndata_bb_hodotdc_bar_tdc_id);
   fChain->SetBranchAddress("bb.hodotdc.bar.tdc.id", bb_hodotdc_bar_tdc_id, &b_bb_hodotdc_bar_tdc_id);
   fChain->SetBranchAddress("Ndata.bb.hodotdc.bar.tdc.meantime", &Ndata_bb_hodotdc_bar_tdc_meantime, &b_Ndata_bb_hodotdc_bar_tdc_meantime);
   fChain->SetBranchAddress("bb.hodotdc.bar.tdc.meantime", bb_hodotdc_bar_tdc_meantime, &b_bb_hodotdc_bar_tdc_meantime);
   fChain->SetBranchAddress("Ndata.bb.hodotdc.bar.tdc.timediff", &Ndata_bb_hodotdc_bar_tdc_timediff, &b_Ndata_bb_hodotdc_bar_tdc_timediff);
   fChain->SetBranchAddress("bb.hodotdc.bar.tdc.timediff", bb_hodotdc_bar_tdc_timediff, &b_bb_hodotdc_bar_tdc_timediff);
   fChain->SetBranchAddress("Ndata.bb.hodotdc.bar.tdc.timehitpos", &Ndata_bb_hodotdc_bar_tdc_timehitpos, &b_Ndata_bb_hodotdc_bar_tdc_timehitpos);
   fChain->SetBranchAddress("bb.hodotdc.bar.tdc.timehitpos", bb_hodotdc_bar_tdc_timehitpos, &b_bb_hodotdc_bar_tdc_timehitpos);
   fChain->SetBranchAddress("Ndata.bb.hodotdc.bar.tdc.vpos", &Ndata_bb_hodotdc_bar_tdc_vpos, &b_Ndata_bb_hodotdc_bar_tdc_vpos);
   fChain->SetBranchAddress("bb.hodotdc.bar.tdc.vpos", bb_hodotdc_bar_tdc_vpos, &b_bb_hodotdc_bar_tdc_vpos);
   fChain->SetBranchAddress("Ndata.bb.ps.clus.col", &Ndata_bb_ps_clus_col, &b_Ndata_bb_ps_clus_col);
   fChain->SetBranchAddress("bb.ps.clus.col", bb_ps_clus_col, &b_bb_ps_clus_col);
   fChain->SetBranchAddress("Ndata.bb.ps.clus.e", &Ndata_bb_ps_clus_e, &b_Ndata_bb_ps_clus_e);
   fChain->SetBranchAddress("bb.ps.clus.e", bb_ps_clus_e, &b_bb_ps_clus_e);
   fChain->SetBranchAddress("Ndata.bb.ps.clus.e_c", &Ndata_bb_ps_clus_e_c, &b_Ndata_bb_ps_clus_e_c);
   fChain->SetBranchAddress("bb.ps.clus.e_c", bb_ps_clus_e_c, &b_bb_ps_clus_e_c);
   fChain->SetBranchAddress("Ndata.bb.ps.clus.eblk", &Ndata_bb_ps_clus_eblk, &b_Ndata_bb_ps_clus_eblk);
   fChain->SetBranchAddress("bb.ps.clus.eblk", bb_ps_clus_eblk, &b_bb_ps_clus_eblk);
   fChain->SetBranchAddress("Ndata.bb.ps.clus.eblk_c", &Ndata_bb_ps_clus_eblk_c, &b_Ndata_bb_ps_clus_eblk_c);
   fChain->SetBranchAddress("bb.ps.clus.eblk_c", bb_ps_clus_eblk_c, &b_bb_ps_clus_eblk_c);
   fChain->SetBranchAddress("Ndata.bb.ps.clus.id", &Ndata_bb_ps_clus_id, &b_Ndata_bb_ps_clus_id);
   fChain->SetBranchAddress("bb.ps.clus.id", bb_ps_clus_id, &b_bb_ps_clus_id);
   fChain->SetBranchAddress("Ndata.bb.ps.clus.nblk", &Ndata_bb_ps_clus_nblk, &b_Ndata_bb_ps_clus_nblk);
   fChain->SetBranchAddress("bb.ps.clus.nblk", bb_ps_clus_nblk, &b_bb_ps_clus_nblk);
   fChain->SetBranchAddress("Ndata.bb.ps.clus.row", &Ndata_bb_ps_clus_row, &b_Ndata_bb_ps_clus_row);
   fChain->SetBranchAddress("bb.ps.clus.row", bb_ps_clus_row, &b_bb_ps_clus_row);
   fChain->SetBranchAddress("Ndata.bb.ps.clus.x", &Ndata_bb_ps_clus_x, &b_Ndata_bb_ps_clus_x);
   fChain->SetBranchAddress("bb.ps.clus.x", bb_ps_clus_x, &b_bb_ps_clus_x);
   fChain->SetBranchAddress("Ndata.bb.ps.clus.y", &Ndata_bb_ps_clus_y, &b_Ndata_bb_ps_clus_y);
   fChain->SetBranchAddress("bb.ps.clus.y", bb_ps_clus_y, &b_bb_ps_clus_y);
   fChain->SetBranchAddress("Ndata.bb.ps.clus_blk.col", &Ndata_bb_ps_clus_blk_col, &b_Ndata_bb_ps_clus_blk_col);
   fChain->SetBranchAddress("bb.ps.clus_blk.col", bb_ps_clus_blk_col, &b_bb_ps_clus_blk_col);
   fChain->SetBranchAddress("Ndata.bb.ps.clus_blk.e", &Ndata_bb_ps_clus_blk_e, &b_Ndata_bb_ps_clus_blk_e);
   fChain->SetBranchAddress("bb.ps.clus_blk.e", bb_ps_clus_blk_e, &b_bb_ps_clus_blk_e);
   fChain->SetBranchAddress("Ndata.bb.ps.clus_blk.e_c", &Ndata_bb_ps_clus_blk_e_c, &b_Ndata_bb_ps_clus_blk_e_c);
   fChain->SetBranchAddress("bb.ps.clus_blk.e_c", bb_ps_clus_blk_e_c, &b_bb_ps_clus_blk_e_c);
   fChain->SetBranchAddress("Ndata.bb.ps.clus_blk.id", &Ndata_bb_ps_clus_blk_id, &b_Ndata_bb_ps_clus_blk_id);
   fChain->SetBranchAddress("bb.ps.clus_blk.id", bb_ps_clus_blk_id, &b_bb_ps_clus_blk_id);
   fChain->SetBranchAddress("Ndata.bb.ps.clus_blk.row", &Ndata_bb_ps_clus_blk_row, &b_Ndata_bb_ps_clus_blk_row);
   fChain->SetBranchAddress("bb.ps.clus_blk.row", bb_ps_clus_blk_row, &b_bb_ps_clus_blk_row);
   fChain->SetBranchAddress("Ndata.bb.ps.clus_blk.x", &Ndata_bb_ps_clus_blk_x, &b_Ndata_bb_ps_clus_blk_x);
   fChain->SetBranchAddress("bb.ps.clus_blk.x", bb_ps_clus_blk_x, &b_bb_ps_clus_blk_x);
   fChain->SetBranchAddress("Ndata.bb.ps.clus_blk.y", &Ndata_bb_ps_clus_blk_y, &b_Ndata_bb_ps_clus_blk_y);
   fChain->SetBranchAddress("bb.ps.clus_blk.y", bb_ps_clus_blk_y, &b_bb_ps_clus_blk_y);
   fChain->SetBranchAddress("Ndata.bb.sh.clus.col", &Ndata_bb_sh_clus_col, &b_Ndata_bb_sh_clus_col);
   fChain->SetBranchAddress("bb.sh.clus.col", bb_sh_clus_col, &b_bb_sh_clus_col);
   fChain->SetBranchAddress("Ndata.bb.sh.clus.e", &Ndata_bb_sh_clus_e, &b_Ndata_bb_sh_clus_e);
   fChain->SetBranchAddress("bb.sh.clus.e", bb_sh_clus_e, &b_bb_sh_clus_e);
   fChain->SetBranchAddress("Ndata.bb.sh.clus.e_c", &Ndata_bb_sh_clus_e_c, &b_Ndata_bb_sh_clus_e_c);
   fChain->SetBranchAddress("bb.sh.clus.e_c", bb_sh_clus_e_c, &b_bb_sh_clus_e_c);
   fChain->SetBranchAddress("Ndata.bb.sh.clus.eblk", &Ndata_bb_sh_clus_eblk, &b_Ndata_bb_sh_clus_eblk);
   fChain->SetBranchAddress("bb.sh.clus.eblk", bb_sh_clus_eblk, &b_bb_sh_clus_eblk);
   fChain->SetBranchAddress("Ndata.bb.sh.clus.eblk_c", &Ndata_bb_sh_clus_eblk_c, &b_Ndata_bb_sh_clus_eblk_c);
   fChain->SetBranchAddress("bb.sh.clus.eblk_c", bb_sh_clus_eblk_c, &b_bb_sh_clus_eblk_c);
   fChain->SetBranchAddress("Ndata.bb.sh.clus.id", &Ndata_bb_sh_clus_id, &b_Ndata_bb_sh_clus_id);
   fChain->SetBranchAddress("bb.sh.clus.id", bb_sh_clus_id, &b_bb_sh_clus_id);
   fChain->SetBranchAddress("Ndata.bb.sh.clus.nblk", &Ndata_bb_sh_clus_nblk, &b_Ndata_bb_sh_clus_nblk);
   fChain->SetBranchAddress("bb.sh.clus.nblk", bb_sh_clus_nblk, &b_bb_sh_clus_nblk);
   fChain->SetBranchAddress("Ndata.bb.sh.clus.row", &Ndata_bb_sh_clus_row, &b_Ndata_bb_sh_clus_row);
   fChain->SetBranchAddress("bb.sh.clus.row", bb_sh_clus_row, &b_bb_sh_clus_row);
   fChain->SetBranchAddress("Ndata.bb.sh.clus.x", &Ndata_bb_sh_clus_x, &b_Ndata_bb_sh_clus_x);
   fChain->SetBranchAddress("bb.sh.clus.x", bb_sh_clus_x, &b_bb_sh_clus_x);
   fChain->SetBranchAddress("Ndata.bb.sh.clus.y", &Ndata_bb_sh_clus_y, &b_Ndata_bb_sh_clus_y);
   fChain->SetBranchAddress("bb.sh.clus.y", bb_sh_clus_y, &b_bb_sh_clus_y);
   fChain->SetBranchAddress("Ndata.bb.sh.clus_blk.col", &Ndata_bb_sh_clus_blk_col, &b_Ndata_bb_sh_clus_blk_col);
   fChain->SetBranchAddress("bb.sh.clus_blk.col", bb_sh_clus_blk_col, &b_bb_sh_clus_blk_col);
   fChain->SetBranchAddress("Ndata.bb.sh.clus_blk.e", &Ndata_bb_sh_clus_blk_e, &b_Ndata_bb_sh_clus_blk_e);
   fChain->SetBranchAddress("bb.sh.clus_blk.e", bb_sh_clus_blk_e, &b_bb_sh_clus_blk_e);
   fChain->SetBranchAddress("Ndata.bb.sh.clus_blk.e_c", &Ndata_bb_sh_clus_blk_e_c, &b_Ndata_bb_sh_clus_blk_e_c);
   fChain->SetBranchAddress("bb.sh.clus_blk.e_c", bb_sh_clus_blk_e_c, &b_bb_sh_clus_blk_e_c);
   fChain->SetBranchAddress("Ndata.bb.sh.clus_blk.id", &Ndata_bb_sh_clus_blk_id, &b_Ndata_bb_sh_clus_blk_id);
   fChain->SetBranchAddress("bb.sh.clus_blk.id", bb_sh_clus_blk_id, &b_bb_sh_clus_blk_id);
   fChain->SetBranchAddress("Ndata.bb.sh.clus_blk.row", &Ndata_bb_sh_clus_blk_row, &b_Ndata_bb_sh_clus_blk_row);
   fChain->SetBranchAddress("bb.sh.clus_blk.row", bb_sh_clus_blk_row, &b_bb_sh_clus_blk_row);
   fChain->SetBranchAddress("Ndata.bb.sh.clus_blk.x", &Ndata_bb_sh_clus_blk_x, &b_Ndata_bb_sh_clus_blk_x);
   fChain->SetBranchAddress("bb.sh.clus_blk.x", bb_sh_clus_blk_x, &b_bb_sh_clus_blk_x);
   fChain->SetBranchAddress("Ndata.bb.sh.clus_blk.y", &Ndata_bb_sh_clus_blk_y, &b_Ndata_bb_sh_clus_blk_y);
   fChain->SetBranchAddress("bb.sh.clus_blk.y", bb_sh_clus_blk_y, &b_bb_sh_clus_blk_y);
   fChain->SetBranchAddress("Ndata.bb.tr.beta", &Ndata_bb_tr_beta, &b_Ndata_bb_tr_beta);
   fChain->SetBranchAddress("bb.tr.beta", bb_tr_beta, &b_bb_tr_beta);
   fChain->SetBranchAddress("Ndata.bb.tr.chi2", &Ndata_bb_tr_chi2, &b_Ndata_bb_tr_chi2);
   fChain->SetBranchAddress("bb.tr.chi2", bb_tr_chi2, &b_bb_tr_chi2);
   fChain->SetBranchAddress("Ndata.bb.tr.d_ph", &Ndata_bb_tr_d_ph, &b_Ndata_bb_tr_d_ph);
   fChain->SetBranchAddress("bb.tr.d_ph", bb_tr_d_ph, &b_bb_tr_d_ph);
   fChain->SetBranchAddress("Ndata.bb.tr.d_th", &Ndata_bb_tr_d_th, &b_Ndata_bb_tr_d_th);
   fChain->SetBranchAddress("bb.tr.d_th", bb_tr_d_th, &b_bb_tr_d_th);
   fChain->SetBranchAddress("Ndata.bb.tr.d_x", &Ndata_bb_tr_d_x, &b_Ndata_bb_tr_d_x);
   fChain->SetBranchAddress("bb.tr.d_x", bb_tr_d_x, &b_bb_tr_d_x);
   fChain->SetBranchAddress("Ndata.bb.tr.d_y", &Ndata_bb_tr_d_y, &b_Ndata_bb_tr_d_y);
   fChain->SetBranchAddress("bb.tr.d_y", bb_tr_d_y, &b_bb_tr_d_y);
   fChain->SetBranchAddress("Ndata.bb.tr.dbeta", &Ndata_bb_tr_dbeta, &b_Ndata_bb_tr_dbeta);
   fChain->SetBranchAddress("bb.tr.dbeta", bb_tr_dbeta, &b_bb_tr_dbeta);
   fChain->SetBranchAddress("Ndata.bb.tr.dtime", &Ndata_bb_tr_dtime, &b_Ndata_bb_tr_dtime);
   fChain->SetBranchAddress("bb.tr.dtime", bb_tr_dtime, &b_bb_tr_dtime);
   fChain->SetBranchAddress("Ndata.bb.tr.flag", &Ndata_bb_tr_flag, &b_Ndata_bb_tr_flag);
   fChain->SetBranchAddress("bb.tr.flag", bb_tr_flag, &b_bb_tr_flag);
   fChain->SetBranchAddress("Ndata.bb.tr.ndof", &Ndata_bb_tr_ndof, &b_Ndata_bb_tr_ndof);
   fChain->SetBranchAddress("bb.tr.ndof", bb_tr_ndof, &b_bb_tr_ndof);
   fChain->SetBranchAddress("Ndata.bb.tr.p", &Ndata_bb_tr_p, &b_Ndata_bb_tr_p);
   fChain->SetBranchAddress("bb.tr.p", bb_tr_p, &b_bb_tr_p);
   fChain->SetBranchAddress("Ndata.bb.tr.pathl", &Ndata_bb_tr_pathl, &b_Ndata_bb_tr_pathl);
   fChain->SetBranchAddress("bb.tr.pathl", bb_tr_pathl, &b_bb_tr_pathl);
   fChain->SetBranchAddress("Ndata.bb.tr.ph", &Ndata_bb_tr_ph, &b_Ndata_bb_tr_ph);
   fChain->SetBranchAddress("bb.tr.ph", bb_tr_ph, &b_bb_tr_ph);
   fChain->SetBranchAddress("Ndata.bb.tr.px", &Ndata_bb_tr_px, &b_Ndata_bb_tr_px);
   fChain->SetBranchAddress("bb.tr.px", bb_tr_px, &b_bb_tr_px);
   fChain->SetBranchAddress("Ndata.bb.tr.py", &Ndata_bb_tr_py, &b_Ndata_bb_tr_py);
   fChain->SetBranchAddress("bb.tr.py", bb_tr_py, &b_bb_tr_py);
   fChain->SetBranchAddress("Ndata.bb.tr.pz", &Ndata_bb_tr_pz, &b_Ndata_bb_tr_pz);
   fChain->SetBranchAddress("bb.tr.pz", bb_tr_pz, &b_bb_tr_pz);
   fChain->SetBranchAddress("Ndata.bb.tr.r_ph", &Ndata_bb_tr_r_ph, &b_Ndata_bb_tr_r_ph);
   fChain->SetBranchAddress("bb.tr.r_ph", bb_tr_r_ph, &b_bb_tr_r_ph);
   fChain->SetBranchAddress("Ndata.bb.tr.r_th", &Ndata_bb_tr_r_th, &b_Ndata_bb_tr_r_th);
   fChain->SetBranchAddress("bb.tr.r_th", bb_tr_r_th, &b_bb_tr_r_th);
   fChain->SetBranchAddress("Ndata.bb.tr.r_x", &Ndata_bb_tr_r_x, &b_Ndata_bb_tr_r_x);
   fChain->SetBranchAddress("bb.tr.r_x", bb_tr_r_x, &b_bb_tr_r_x);
   fChain->SetBranchAddress("Ndata.bb.tr.r_y", &Ndata_bb_tr_r_y, &b_Ndata_bb_tr_r_y);
   fChain->SetBranchAddress("bb.tr.r_y", bb_tr_r_y, &b_bb_tr_r_y);
   fChain->SetBranchAddress("Ndata.bb.tr.tg_dp", &Ndata_bb_tr_tg_dp, &b_Ndata_bb_tr_tg_dp);
   fChain->SetBranchAddress("bb.tr.tg_dp", bb_tr_tg_dp, &b_bb_tr_tg_dp);
   fChain->SetBranchAddress("Ndata.bb.tr.tg_ph", &Ndata_bb_tr_tg_ph, &b_Ndata_bb_tr_tg_ph);
   fChain->SetBranchAddress("bb.tr.tg_ph", bb_tr_tg_ph, &b_bb_tr_tg_ph);
   fChain->SetBranchAddress("Ndata.bb.tr.tg_th", &Ndata_bb_tr_tg_th, &b_Ndata_bb_tr_tg_th);
   fChain->SetBranchAddress("bb.tr.tg_th", bb_tr_tg_th, &b_bb_tr_tg_th);
   fChain->SetBranchAddress("Ndata.bb.tr.tg_x", &Ndata_bb_tr_tg_x, &b_Ndata_bb_tr_tg_x);
   fChain->SetBranchAddress("bb.tr.tg_x", bb_tr_tg_x, &b_bb_tr_tg_x);
   fChain->SetBranchAddress("Ndata.bb.tr.tg_y", &Ndata_bb_tr_tg_y, &b_Ndata_bb_tr_tg_y);
   fChain->SetBranchAddress("bb.tr.tg_y", bb_tr_tg_y, &b_bb_tr_tg_y);
   fChain->SetBranchAddress("Ndata.bb.tr.th", &Ndata_bb_tr_th, &b_Ndata_bb_tr_th);
   fChain->SetBranchAddress("bb.tr.th", bb_tr_th, &b_bb_tr_th);
   fChain->SetBranchAddress("Ndata.bb.tr.time", &Ndata_bb_tr_time, &b_Ndata_bb_tr_time);
   fChain->SetBranchAddress("bb.tr.time", bb_tr_time, &b_bb_tr_time);
   fChain->SetBranchAddress("Ndata.bb.tr.vx", &Ndata_bb_tr_vx, &b_Ndata_bb_tr_vx);
   fChain->SetBranchAddress("bb.tr.vx", bb_tr_vx, &b_bb_tr_vx);
   fChain->SetBranchAddress("Ndata.bb.tr.vy", &Ndata_bb_tr_vy, &b_Ndata_bb_tr_vy);
   fChain->SetBranchAddress("bb.tr.vy", bb_tr_vy, &b_bb_tr_vy);
   fChain->SetBranchAddress("Ndata.bb.tr.vz", &Ndata_bb_tr_vz, &b_Ndata_bb_tr_vz);
   fChain->SetBranchAddress("bb.tr.vz", bb_tr_vz, &b_bb_tr_vz);
   fChain->SetBranchAddress("Ndata.bb.tr.x", &Ndata_bb_tr_x, &b_Ndata_bb_tr_x);
   fChain->SetBranchAddress("bb.tr.x", bb_tr_x, &b_bb_tr_x);
   fChain->SetBranchAddress("Ndata.bb.tr.y", &Ndata_bb_tr_y, &b_Ndata_bb_tr_y);
   fChain->SetBranchAddress("bb.tr.y", bb_tr_y, &b_bb_tr_y);
   fChain->SetBranchAddress("Ndata.bb.x_bcp", &Ndata_bb_x_bcp, &b_Ndata_bb_x_bcp);
   fChain->SetBranchAddress("bb.x_bcp", bb_x_bcp, &b_bb_x_bcp);
   fChain->SetBranchAddress("Ndata.bb.x_fcp", &Ndata_bb_x_fcp, &b_Ndata_bb_x_fcp);
   fChain->SetBranchAddress("bb.x_fcp", bb_x_fcp, &b_bb_x_fcp);
   fChain->SetBranchAddress("Ndata.bb.y_bcp", &Ndata_bb_y_bcp, &b_Ndata_bb_y_bcp);
   fChain->SetBranchAddress("bb.y_bcp", bb_y_bcp, &b_bb_y_bcp);
   fChain->SetBranchAddress("Ndata.bb.y_fcp", &Ndata_bb_y_fcp, &b_Ndata_bb_y_fcp);
   fChain->SetBranchAddress("bb.y_fcp", bb_y_fcp, &b_bb_y_fcp);
   fChain->SetBranchAddress("Ndata.sbs.hcal.a_amp", &Ndata_sbs_hcal_a_amp, &b_Ndata_sbs_hcal_a_amp);
   fChain->SetBranchAddress("sbs.hcal.a_amp", sbs_hcal_a_amp, &b_sbs_hcal_a_amp);
   fChain->SetBranchAddress("Ndata.sbs.hcal.a_c", &Ndata_sbs_hcal_a_c, &b_Ndata_sbs_hcal_a_c);
   fChain->SetBranchAddress("sbs.hcal.a_c", sbs_hcal_a_c, &b_sbs_hcal_a_c);
   fChain->SetBranchAddress("Ndata.sbs.hcal.a_p", &Ndata_sbs_hcal_a_p, &b_Ndata_sbs_hcal_a_p);
   fChain->SetBranchAddress("sbs.hcal.a_p", sbs_hcal_a_p, &b_sbs_hcal_a_p);
   fChain->SetBranchAddress("Ndata.sbs.hcal.a_time", &Ndata_sbs_hcal_a_time, &b_Ndata_sbs_hcal_a_time);
   fChain->SetBranchAddress("sbs.hcal.a_time", sbs_hcal_a_time, &b_sbs_hcal_a_time);
   fChain->SetBranchAddress("Ndata.sbs.hcal.adccol", &Ndata_sbs_hcal_adccol, &b_Ndata_sbs_hcal_adccol);
   fChain->SetBranchAddress("sbs.hcal.adccol", sbs_hcal_adccol, &b_sbs_hcal_adccol);
   fChain->SetBranchAddress("Ndata.sbs.hcal.adcrow", &Ndata_sbs_hcal_adcrow, &b_Ndata_sbs_hcal_adcrow);
   fChain->SetBranchAddress("sbs.hcal.adcrow", sbs_hcal_adcrow, &b_sbs_hcal_adcrow);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.col", &Ndata_sbs_hcal_clus_col, &b_Ndata_sbs_hcal_clus_col);
   fChain->SetBranchAddress("sbs.hcal.clus.col", sbs_hcal_clus_col, &b_sbs_hcal_clus_col);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.e", &Ndata_sbs_hcal_clus_e, &b_Ndata_sbs_hcal_clus_e);
   fChain->SetBranchAddress("sbs.hcal.clus.e", sbs_hcal_clus_e, &b_sbs_hcal_clus_e);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.eblk", &Ndata_sbs_hcal_clus_eblk, &b_Ndata_sbs_hcal_clus_eblk);
   fChain->SetBranchAddress("sbs.hcal.clus.eblk", sbs_hcal_clus_eblk, &b_sbs_hcal_clus_eblk);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.id", &Ndata_sbs_hcal_clus_id, &b_Ndata_sbs_hcal_clus_id);
   fChain->SetBranchAddress("sbs.hcal.clus.id", sbs_hcal_clus_id, &b_sbs_hcal_clus_id);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.nblk", &Ndata_sbs_hcal_clus_nblk, &b_Ndata_sbs_hcal_clus_nblk);
   fChain->SetBranchAddress("sbs.hcal.clus.nblk", sbs_hcal_clus_nblk, &b_sbs_hcal_clus_nblk);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.row", &Ndata_sbs_hcal_clus_row, &b_Ndata_sbs_hcal_clus_row);
   fChain->SetBranchAddress("sbs.hcal.clus.row", sbs_hcal_clus_row, &b_sbs_hcal_clus_row);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus_blk.e", &Ndata_sbs_hcal_clus_blk_e, &b_Ndata_sbs_hcal_clus_blk_e);
   fChain->SetBranchAddress("sbs.hcal.clus_blk.e", sbs_hcal_clus_blk_e, &b_sbs_hcal_clus_blk_e);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus_blk.id", &Ndata_sbs_hcal_clus_blk_id, &b_Ndata_sbs_hcal_clus_blk_id);
   fChain->SetBranchAddress("sbs.hcal.clus_blk.id", sbs_hcal_clus_blk_id, &b_sbs_hcal_clus_blk_id);
   fChain->SetBranchAddress("BB.gold.beta", &BB_gold_beta, &b_BB_gold_beta);
   fChain->SetBranchAddress("BB.gold.dp", &BB_gold_dp, &b_BB_gold_dp);
   fChain->SetBranchAddress("BB.gold.index", &BB_gold_index, &b_BB_gold_index);
   fChain->SetBranchAddress("BB.gold.ok", &BB_gold_ok, &b_BB_gold_ok);
   fChain->SetBranchAddress("BB.gold.p", &BB_gold_p, &b_BB_gold_p);
   fChain->SetBranchAddress("BB.gold.ph", &BB_gold_ph, &b_BB_gold_ph);
   fChain->SetBranchAddress("BB.gold.px", &BB_gold_px, &b_BB_gold_px);
   fChain->SetBranchAddress("BB.gold.py", &BB_gold_py, &b_BB_gold_py);
   fChain->SetBranchAddress("BB.gold.pz", &BB_gold_pz, &b_BB_gold_pz);
   fChain->SetBranchAddress("BB.gold.th", &BB_gold_th, &b_BB_gold_th);
   fChain->SetBranchAddress("BB.gold.x", &BB_gold_x, &b_BB_gold_x);
   fChain->SetBranchAddress("BB.gold.y", &BB_gold_y, &b_BB_gold_y);
   fChain->SetBranchAddress("bb.gem.hit.ngoodhits", &bb_gem_hit_ngoodhits, &b_bb_gem_hit_ngoodhits);
   fChain->SetBranchAddress("bb.gem.nlayershit", &bb_gem_nlayershit, &b_bb_gem_nlayershit);
   fChain->SetBranchAddress("bb.gem.nlayershitu", &bb_gem_nlayershitu, &b_bb_gem_nlayershitu);
   fChain->SetBranchAddress("bb.gem.nlayershituv", &bb_gem_nlayershituv, &b_bb_gem_nlayershituv);
   fChain->SetBranchAddress("bb.gem.nlayershitv", &bb_gem_nlayershitv, &b_bb_gem_nlayershitv);
   fChain->SetBranchAddress("bb.gem.track.besttrack", &bb_gem_track_besttrack, &b_bb_gem_track_besttrack);
   fChain->SetBranchAddress("bb.gem.track.ntrack", &bb_gem_track_ntrack, &b_bb_gem_track_ntrack);
   fChain->SetBranchAddress("bb.grinch_tdc.ngoodADChits", &bb_grinch_tdc_ngoodADChits, &b_bb_grinch_tdc_ngoodADChits);
   fChain->SetBranchAddress("bb.grinch_tdc.ngoodTDChits", &bb_grinch_tdc_ngoodTDChits, &b_bb_grinch_tdc_ngoodTDChits);
   fChain->SetBranchAddress("bb.grinch_tdc.nhits", &bb_grinch_tdc_nhits, &b_bb_grinch_tdc_nhits);
   fChain->SetBranchAddress("bb.grinch_tdc.nrefhits", &bb_grinch_tdc_nrefhits, &b_bb_grinch_tdc_nrefhits);
   fChain->SetBranchAddress("bb.hodotdc.bar.ngoodbars", &bb_hodotdc_bar_ngoodbars, &b_bb_hodotdc_bar_ngoodbars);
   fChain->SetBranchAddress("bb.hodotdc.nclus", &bb_hodotdc_nclus, &b_bb_hodotdc_nclus);
   fChain->SetBranchAddress("bb.ps.e", &bb_ps_e, &b_bb_ps_e);
   fChain->SetBranchAddress("bb.ps.e_c", &bb_ps_e_c, &b_bb_ps_e_c);
   fChain->SetBranchAddress("bb.ps.nblk", &bb_ps_nblk, &b_bb_ps_nblk);
   fChain->SetBranchAddress("bb.ps.nclus", &bb_ps_nclus, &b_bb_ps_nclus);
   fChain->SetBranchAddress("bb.ps.ngoodADChits", &bb_ps_ngoodADChits, &b_bb_ps_ngoodADChits);
   fChain->SetBranchAddress("bb.ps.x", &bb_ps_x, &b_bb_ps_x);
   fChain->SetBranchAddress("bb.ps.y", &bb_ps_y, &b_bb_ps_y);
   fChain->SetBranchAddress("bb.sh.e", &bb_sh_e, &b_bb_sh_e);
   fChain->SetBranchAddress("bb.sh.e_c", &bb_sh_e_c, &b_bb_sh_e_c);
   fChain->SetBranchAddress("bb.sh.nblk", &bb_sh_nblk, &b_bb_sh_nblk);
   fChain->SetBranchAddress("bb.sh.nclus", &bb_sh_nclus, &b_bb_sh_nclus);
   fChain->SetBranchAddress("bb.sh.ngoodADChits", &bb_sh_ngoodADChits, &b_bb_sh_ngoodADChits);
   fChain->SetBranchAddress("bb.sh.x", &bb_sh_x, &b_bb_sh_x);
   fChain->SetBranchAddress("bb.sh.y", &bb_sh_y, &b_bb_sh_y);
   fChain->SetBranchAddress("bb.tr.n", &bb_tr_n, &b_bb_tr_n);
   fChain->SetBranchAddress("sbs.hcal.e", &sbs_hcal_e, &b_sbs_hcal_e);
   fChain->SetBranchAddress("sbs.hcal.nblk", &sbs_hcal_nblk, &b_sbs_hcal_nblk);
   fChain->SetBranchAddress("sbs.hcal.nclus", &sbs_hcal_nclus, &b_sbs_hcal_nclus);
   fChain->SetBranchAddress("sbs.hcal.ngoodADChits", &sbs_hcal_ngoodADChits, &b_sbs_hcal_ngoodADChits);
   fChain->SetBranchAddress("singletrack", &singletrack, &b_singletrack);
   fChain->SetBranchAddress("anytrack", &anytrack, &b_anytrack);
   fChain->SetBranchAddress("fEvtHdr.fEvtTime", &fEvtHdr_fEvtTime, &b_Event_Branch_fEvtHdr_fEvtTime);
   fChain->SetBranchAddress("fEvtHdr.fEvtNum", &fEvtHdr_fEvtNum, &b_Event_Branch_fEvtHdr_fEvtNum);
   fChain->SetBranchAddress("fEvtHdr.fEvtType", &fEvtHdr_fEvtType, &b_Event_Branch_fEvtHdr_fEvtType);
   fChain->SetBranchAddress("fEvtHdr.fEvtLen", &fEvtHdr_fEvtLen, &b_Event_Branch_fEvtHdr_fEvtLen);
   fChain->SetBranchAddress("fEvtHdr.fHelicity", &fEvtHdr_fHelicity, &b_Event_Branch_fEvtHdr_fHelicity);
   fChain->SetBranchAddress("fEvtHdr.fTargetPol", &fEvtHdr_fTargetPol, &b_Event_Branch_fEvtHdr_fTargetPol);
   fChain->SetBranchAddress("fEvtHdr.fRun", &fEvtHdr_fRun, &b_Event_Branch_fEvtHdr_fRun);
   Notify();
}

Bool_t gmn_rec_tree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void gmn_rec_tree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t gmn_rec_tree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef gmn_rec_tree_cxx
