#!/bin/sh

#script to set up the environment for SBS-offline
export SBS=/adaqfs/home/a-onl/sbs/sbs_devel/install
export SBSOFFLINE=/adaqfs/home/a-onl/sbs/sbs_devel/install
#export SBS_REPLAY=/adaqfs/home/a-onl/sbs/sbs_devel/SBS-replay

if test "x$PATH" = "x" ; then
    export PATH=/adaqfs/home/a-onl/sbs/sbs_devel/install/bin
else
    export PATH=/adaqfs/home/a-onl/sbs/sbs_devel/install/bin:$PATH
fi

OS=`uname -s`


if [ "$OS" = "Darwin" ]
then # Mac OS: set DYLD_LIBRARY_PATH to library directory:
    if test "x$DYLD_LIBRARY_PATH" = "x"; then
	export DYLD_LIBRARY_PATH=/adaqfs/home/a-onl/sbs/sbs_devel/install/lib64
    else
	export DYLD_LIBRARY_PATH=/adaqfs/home/a-onl/sbs/sbs_devel/install/lib64:$DYLD_LIBRARY_PATH
    fi
fi

# set LD_LIBRARY_PATH regardless of OS:
if test "x$LD_LIBRARY_PATH" = "x"; then
    export LD_LIBRARY_PATH=/adaqfs/home/a-onl/sbs/sbs_devel/install/lib64
else
    export LD_LIBRARY_PATH=/adaqfs/home/a-onl/sbs/sbs_devel/install/lib64:$LD_LIBRARY_PATH
fi


