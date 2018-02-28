library(qtl2pleio)
library(stringr)
library(dplyr)
library(readr)
args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
# translate command line args
print(args)
print(args$argname)
(proc_num <- as.numeric(args$argname))
(run_num <- as.numeric(args$run_num))
(nboot_per_job <- args$nboot_per_job) # think 1 boot per job
###############
library(qtl2)
file <- paste0("https://raw.githubusercontent.com/rqtl/",
               "qtl2data/master/DO_Recla/recla.zip")
recla <- read_cross2(file)
# make sex a covariate for use in pvl_scan
recla[[6]][ , 1, drop = FALSE] -> sex
# insert pseudomarkers
insert_pseudomarkers(recla, step = 0.1) -> pseudomap

## ------------------------------------------------------------------------
lapply(FUN = length, X = pseudomap) -> lens
sum(lens)

## ------------------------------------------------------------------------
probs <- calc_genoprob(recla, map = pseudomap)

## ------------------------------------------------------------------------
aprobs <- genoprob_to_alleleprob(probs)

## ------------------------------------------------------------------------
kinship <- calc_kinship(aprobs, "loco")

## ------------------------------------------------------------------------
recla$pheno -> ph
log(ph) -> lph
apply(FUN = broman::winsorize, X = lph, MARGIN = 2) -> wlph

## ------------------------------------------------------------------------
phe <- wlph[, c(7, 10)]
pp <- aprobs$`8`
gm <- pseudomap$`8`
## ------------------------------------------------------------------------
library(qtl2pleio)


###############
# simulate a phenotype
X1 <- pp[ , , pleio_peak_index]
cbind(X1, sex) -> Xpre
gemma2::stagger_mats(Xpre, Xpre) -> X
set.seed(proc_num)
calc_covs(pheno = phe, kinship = k) -> cc_out
(cc_out$Vg -> Vg)
(cc_out$Ve -> Ve)
# calculate Sigma
calc_Sigma(Vg = Vg, Ve = Ve, K = k) -> Sigma
solve(Sigma) -> Sigma_inv
# calc Bhat 
B <- calc_Bhat(X = X, Sigma_inv = Sigma_inv, Y = phe)
# Start loop
lrt <- numeric()
for (i in 1:nboot_per_job){
  sim1(X = X, B = B, Vg = Vg, Ve = Ve, kinship = k) -> foo
  matrix(foo, ncol = 2, byrow = FALSE) -> Ysim
  scan_pvl(probs = pp, pheno = Ysim, kinship = k, start_snp1 = s1, n_snp = nsnp) -> loglik
# in above call, s1 & nsnp come from command line args
  calc_lrt_tib(loglik) -> lrt[i]
}

fn_out <- paste0("recla-boot-run", run_num, "_", proc_num, ".txt")
write.table(lrt, fn_out)
q("no")


