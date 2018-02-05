## GOAL: run pvl scans - one for each pair of local expression traits near chr 11 hotspot

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
run_num <- as.numeric(run_num)
print(run_num)
############### There are 24 local genes of interest. So we put 24 in the modular division below
1 + (proc_num %% 24) -> local_id
1 + (proc_num %/% 24) -> nonlocal_id

# first, load setup2.R
source("Rscript/setup2.R")
#source("Rscript/setup-chr17-G83.R")
library(dplyr)
library(qtl2pleio)
load("probs_1.RData")
. -> pp
pm <- pmap$`1`
kinship <- K$`1`

## Determine which markers are shared between pmap & gmap
snp_g <- dimnames(pp)[[3]]
snp_p <- names(pm)
shared_snps <- intersect(snp_g, snp_p)
pp2 <- pp[ , , snp_g %in% shared_snps]
pm2 <- pm[snp_p %in% shared_snps]

# assign which columns of expression data object to analyze
# use proc_num & ncol(expr)
expr <- readRDS("data/chr1-34Mb-local.rds")
#ncol(expr) -> k
phe_local <- expr[, local_id]
readRDS("data/insulin_secretion.rds") -> nonlocal_traits
nonlocal_traits[ , nonlocal_id] -> phe_nonlocal
phenames <- c(colnames(expr)[local_id], colnames(nonlocal_traits)[nonlocal_id])
samples_to_drop <- c(360, 370, 268, 269, 309, 310) %>% as.character()
#### find samples that are listed in both phenotype objects. 
# note that expr has 378 mice, while insulin secr phenotypes object has 500 mice.
rownames(nonlocal_traits) -> rn_nonlocal
rownames(expr) -> rn_local
intersect(rn_local, rn_nonlocal) -> shared_phe_rn
phe_nonlocal2 <- phe_nonlocal[rn_nonlocal %in% shared_phe_rn]
phe_local2 <- phe_local[rn_local %in% shared_phe_rn]
# are the two phenotypes ordered the same way?
names(phe_local2) == names(phe_nonlocal2)
####
indic <- (is.na(phe_local2) | is.na(phe_nonlocal2))
phe_nona <- cbind(phe_local2, phe_nonlocal2)[!indic, ]
(phe <- phe_nona[!rownames(phe_nona) %in% samples_to_drop, ])
rownames(phe) <- stringr::str_replace(rownames(phe), "DO-", "") %>% as.numeric() 
rownames(pp2) %in% rownames(phe) -> pp_ind
length(pp_ind) == nrow(pp2)
pp2[pp_ind, , ] -> pp3
dim(pp3)
intersect(rownames(phe), rownames(pp3)) -> rn
pp3[rownames(pp3) %in% rn, , ] -> pp4
phe[rownames(phe) %in% rn, ] -> phe3
rownames(phe3) == rownames(pp4)
cbind(rownames(phe3), rownames(pp4))
pp4[order(rownames(pp4)), , ] -> pp5
phe3[order(rownames(phe3)), ] -> phe4
rownames(phe4) == rownames(pp5)
dim(pp5)
dim(phe4)

rownames(kinship) %in% rownames(phe4) -> kin_ind
k2 <- kinship[kin_ind, kin_ind]
k3 <- k2[order(rownames(k2)), order(rownames(k2))]
rownames(k3) == rownames(phe4)
# scan to find the pleiotropy peak
scan_pvl(probs = pp5, pheno = phe4, kinship = k3, start_snp1 = s1, n_snp = nsnp) -> out
#
fn_out <- paste0("pvl-run", run_num, "_", proc_num, "_", paste(phenames, collapse = "_"), ".rds")
#save(list = "out", file = fn)
saveRDS(object = out, file = fn_out)
q("no")
