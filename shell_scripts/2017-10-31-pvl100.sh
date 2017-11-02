#!/bin/bash

# untar your R installation
tar -xzf R.tar.gz

# make sure the script will use your R installation
export PATH=$(pwd)/R/bin:$PATH

# run R, with the name of your  R script
#R CMD BATCH '--args argname='$1'' Rscript/2017-09-25-boot.R ''$1'.Rout'
R CMD BATCH '--args argname='$1' nsnp='$2' s1='$3' index='$4' nboot_per_job='$5'' Rscript/2017-10-31-pvl100.R ''$1'_pvl100.Rout'

