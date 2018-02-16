#!/bin/bash

# untar your R installation
tar -xzf R6.tar.gz
tar -xzf SLIBS.tar.gz
# make sure the script will use your R installation
export PATH=$(pwd)/R/bin:$PATH
export LD_LIBRARY_PATH=$(pwd)/SS:$LD_LIBRARY_PATH
# run R, with the name of your  R script
R CMD BATCH '--args argname='$1' run_num='$2' chr='$3' phe1_name='$4' phe2_name='$5' peak1='$6' peak2='$7' probs_file='$8'' Rscript/pvl-chr17.R 'pvl_run'$2'-'$1'.Rout'


#rm *.Rout
