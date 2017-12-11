#!/bin/bash

# untar your R installation
tar -xzf R.tar.gz
tar -xzf SLIBS.tar.gz
# make sure the script will use your R installation
export PATH=$(pwd)/R/bin:$PATH
export LD_LIBRARY_PATH=$(pwd)/SS:$LD_LIBRARY_PATH
# run R, with the name of your  R script
#R CMD BATCH '--args argname='$1'' Rscript/2017-09-25-boot.R ''$1'.Rout'
R CMD BATCH '--args argname='$1' nsnp='$2' s1='$3' index='$4' nboot_per_job='$5'' Rscript/2017-12-07_pvl400-typeI-error-run3.R '2017-12-07_'$1'_pvl400-run3.Rout'

