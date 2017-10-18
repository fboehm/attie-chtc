## ------------------------------------------------------------------------
##First read in the arguments listed at the command line
args <- commandArgs(TRUE)
print(args)
##args is now a list of character vectors
## First check to see if arguments are passed.
## Then cycle through each element of the list and evaluate the expressions.
if(length(args)==0){
  print("No arguments supplied.")
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}
print(argname)
proc_num <- as.numeric(argname)
print(proc_num)
nboot <- as.numeric(nboot)
print(nboot)
snp1 <- as.numeric(snp1)
nn <- as.numeric(nn)
n_snp <- as.numeric(n_snp)


yid <- proc_num %% nboot

## ------------------------------------------------------------------------
library(dplyr)
library(qtl2pleio)
library(gemma2)
library(pleiotropy)
# load data
phe_nona <- as.matrix(read.csv(paste0("sim_data/2017-10-11-sim_pleio_", yid, ".csv")))
source("Rscript/setup.R")
source("Rscript/setup-chr2a.R")
## ------------------------------------------------------------------------
## set start POINT HERE
#snp1 <- 175



start_snp_i <- (i - 1)* 10 + snp1
start_snp_j <- (j - 1)* 10 + snp1


scan_pvl(probs = pp3, pheno = phe_nona, kinship = kinship, start_snp1 = start_snp_i, start_snp2 = start_snp_j, n_snp = 10) -> scan_out


fn <- paste0("boot-out_", proc_num, ".csv")

write.csv(scan_out, fn)

q("no")



