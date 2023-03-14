#!bin/sh

echo "\n"
echo "--------------------------------------------------------------------"
echo "Setup script run and environment parameters loaded for HCal analysis."
echo "Current JLab environment " ${JLAB_VERSION}
echo "--------------------------------------------------------------------"
echo "\n"

#Set environment variables in bash
setenv CONFIG_DIR /w/halla-scshelf2102/sbs/seeds/HCal_replay/hcal/hcalCalibration/config
# setenv GMN_DATA_DIR /w/halla-scshelf2102/sbs/sbs-gmn
# setenv PARSE_OUTPUT_DIR /lustre19/expphy/volatile/halla/sbs/seeds/parse
# setenv GMN_OUTPUT_DIR /lustre19/expphy/volatile/halla/sbs/seeds/gmn_analysis

# #Include path to header files
# setenv PATH /w/halla-scshelf2102/sbs/seeds/gmn_analysis/gmna/include:/w/halla-scshelf2102/sbs/seeds/gmn_analysis/gmna/physics:${PATH}
# #setenv LD_LIBRARY_PATH /w/halla-scshelf2102/sbs/seeds/gmn_analysis/gmna/include:${LD_LIBRARY_PATH}
