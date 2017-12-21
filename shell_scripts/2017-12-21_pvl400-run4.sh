#!/bin/bash

# untar your R installation
tar -xzf R.tar.gz
tar -xzf SLIBS.tar.gz
# make sure the script will use your R installation
export PATH=$(pwd)/R/bin:$PATH
export LD_LIBRARY_PATH=$(pwd)/SS:$LD_LIBRARY_PATH
# run R, with the name of your  R script
R CMD BATCH '--args argname='$1' nsnp='$2' s1='$3' index='$4' nboot_per_job='$5'' Rscript/2017-12-21_pvl400-run4.R ''$1'_pvl400-run4.Rout'

