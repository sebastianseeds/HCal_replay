# HCal_replay
SBS HCal Combined Analysis Platform
The purpose of this repository is to provide a software platform for evaluating the performance of the SBS Hadron Calorimeter, to process data collected from the calorimeter, and to generate calibration parameters for the calorimeter. A history of significant configurations are also included. Please reference each individual script or file for usage and details.

Please see [wiki](https://sbs.jlab.org/wiki/index.php/HOW_TOs#BigBite_Spectrometer_.28BB.29) for more instructions.
---
Directory and subdirectory descriptions follow:
## hcal
- BB_HCal_Coorelations: Contains scripts which produce elastic correlations between the HCal hadron and BigBite e' elastic correlations. Useful with limited statistics and insufficient optics to verify signal in HCal.
- clusterTiming: Produces plots comparing fADC time between cluster elements for evaluation of HCal DB parameter sbs.hcal.tmax for max timing cut on cluster element inclusion.
- config: Contains PMT serial numbers and layout in calorimeter.
- DB_old: Contains old and out-of-use database parameters with timestamps for reference or rollback.
- detection_efficiency: Produces plots which evaluate the relative detection efficiency between protons and neutrons in the calorimeter.
- displays: Produces interactive GUIs which allow for viewing of full waveform data from the calorimeter. Also contains similar display GUIs for clustering and the dedicated trigger signal fADC.
- energyCalibration: The central platform containing scripts necessary to produce fADC gain parameters (included in the HCal DB file db_sbs.hcal.dat) from in-beam signal. With these parameters can also produce target HV parameters for reset. NOTE: Changing HV should not be done unless absolutely necessary. Consult SME before making any changes to HV.
- GUI: Repo of all necessary GUI to control HCal functions for recording purposes. NOTE: These copies are not meant to be active for use.
- HV: Contains HV GUI configuration files for record.
- oldCode (Under Maintenance): Contains many out-of-use scripts used to produce simple comparisons and plots. Also contains alpha extraction, PMT quantum efficiency analysis, and cosmic fADC gain calibration scripts. These last three are under maintenance and will be moved up to hcal when complete.
- outFiles: Contains larger .root output analysis files which are kept for comparisons.
- parseElastics: Script to process several data files into a single large root file with tuned elastic cuts to reduce the file size and centralize elastic cut parameters.
- protonDisplay: Contains a script to make simple elastic cuts available without tracking information on LH2 data over sparse statistics for use during commissioning to find signal and evaluate functionality of HCal.
- misc: Has several miscillaneous scripts and files to evaluate parameters and modify the analysis environment. These include HCal BBCal trigger time difference evaluation and a gmn tree generator script.
- setFiles: Contains files which contain current HCal detector configurations and characteristics used for analysis throughout the platform.
- signal_centering: Simple script used to center fADC time in the ADC window with limited information from BigBite. Constructed in response to the sudden-latency-shift issues from GMn running 2021.
- TDC_alignment: Scripts to evaluate the timing resolution and produced offset parameters to align TDC/fADC time signals corresponding to elastic events.
- timingCalibration: Scripts to correct timing signals by event via observation of timewalk, TOF, and trigger jitter. 
- TOF (*temporary/under-construction*): Contains standalone script to evaluate the time-of-flight corrections to timing calibrations for HCal. Will be added to timingCalibration when complete
- trigger_analysis: Computes the timing difference between the various triggers including LED/Cosmic-paddles/Overlapping-regions/BBCal/etc triggers.

## HV
Contains many production HCal and GRINCH HV settings for record.

## Panguin
Contains functional panguin build for evaluation of production-level online monitoring script updates for HCal

## Replay
Contains replay scripts for use with standalone SBS-Offline and analyzer builds for evaluation of HCal tree variables not included in standard replay. Also can be used to evaluate changes to SBS-Offline when needed.
> Current Replay Script: replay_hcal_SAS_general.C

NOTE: Replay will reference any *.odef/*.cdef/*.dat contained in replay directory before $DB_DIR

## scripts
Contains several scripts to issue swif2 batch farm general and standalone replay jobs.