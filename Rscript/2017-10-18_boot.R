# translate command line args
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
nsnp <- as.numeric(nsnp)
print(nsnp)
s1 <- as.numeric(s1)
print(s1)
index <- as.numeric(index)

###############
# first, load setup.R
source("setup.R")
source("setup-chr17.R")
library(qtl2pleio)

# simulate a phenotype
X1 <- pp3[ , , index] #index is from command line args
pleiotropy::stagger_mats(X1, X1) -> X
set.seed(proc_num)
calc_covs(pheno = phe_nona, kinship = kinship_nona) -> cc_out
cc_out$Vg -> Vg
cc_out$Ve -> Ve
# calculate Sigma
calc_Sigma(Vg = Vg, Ve = Ve, K = kinship_nona) -> Sigma
solve(Sigma) -> Sigma_inv
# calc Bhat 
B <- calc_Bhat(X = X, Sigma_inv = Sigma_inv, Y = phe_nona)

sim1(X = X, B = B, Vg = Vg, Ve = Ve, kinship = kinship_nona) -> Ysim
scan_pvl(probs = pp3, pheno = Ysim, kinship = kinship_nona, start_snp1 = s1, n_snp = nsnp) -> loglik_mat
# in above call, s1 & nsnp come from command line args
calc_lrt(loglik_mat) -> lrt
fn <- paste0("2017-10-18_boot_", proc_num, "_", paste(phenames, collapse = "_"), ".txt")
write.table(lrt, fn)
q("no")


