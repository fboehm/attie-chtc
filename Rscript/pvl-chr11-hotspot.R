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
k <- as.integer(k)
print(k)
###############
# first, load setup2.R
source("Rscript/setup2.R")
#source("Rscript/setup-chr17-G83.R")
library(dplyr)
library(qtl2pleio)
load("probs_11.RData")
. -> pp
pm <- pmap$`11`
kinship <- K$`11`

## Determine which markers are shared between pmap & gmap
snp_g <- dimnames(pp)[[3]]
snp_p <- names(pm)
shared_snps <- intersect(snp_g, snp_p)
pp2 <- pp[ , , snp_g %in% shared_snps]
pm2 <- pm[snp_p %in% shared_snps]

# assign which columns of expression data object to analyze
# use proc_num & ncol(expr)
load("data/expr-chr11-hotspot-local.RData")
#ncol(expr) -> k
make_indices <- function(proc_num, k){
  stopifnot(proc_num < choose(k, 2))
  cumsum((k - 1):1) -> cuts1
  cuts2 <- c(0, cuts1)
  which(!(1 + proc_num > cuts2))[1] - 1 -> index1
  index2 <- index1 + 1 + proc_num - cuts2[index1]
  return(c(index1, index2))
}
make_indices(proc_num, k) -> expr_cols 
phe_pre <- expr[, expr_cols]

phe_pre[ , 1] -> phe1
phe_pre[ , 2] -> phe2
phenames <- colnames(phe_pre)
samples_to_drop <- c(360, 370, 268, 269, 309, 310) %>% as.character()

indic <- (is.na(phe1) | is.na(phe2))
phe_nona <- phe_pre[!indic, ]

(phe <- phe_nona[!rownames(phe_nona) %in% samples_to_drop, ])
rownames(phe) <- stringr::str_replace(rownames(phe), "DO-", "") %>% as.numeric() 

rownames(pp2) %in% rownames(phe) -> pp_ind
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
print(phe4)
rownames(kinship) %in% rownames(phe4) -> kin_ind
k2 <- kinship[kin_ind, kin_ind]
dim(k2)

k3 <- k2[order(rownames(k2)), order(rownames(k2))]

rownames(k3) == rownames(phe4)
# scan to find the pleiotropy peak
scan_pvl(probs = pp5, pheno = phe4, kinship = k3, start_snp1 = s1, n_snp = nsnp) -> out

fn_out <- paste0("pvl-run", run_num, "_", proc_num, "_", paste(phenames, collapse = "_"), ".txt")
#save(list = "out", file = fn)
write.table(out, fn_out, quote = FALSE)
q("no")
