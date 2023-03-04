#!/bin/bash

# set up the analyzer for analysis needs

echo "Setting up the Analyzer.."
echo "Executing: module use /group/halla/modulefiles.."
module use /group/halla/modulefiles
echo "Executing: module load analyzer.."
module load analyzer
echo "--> Done!" 
