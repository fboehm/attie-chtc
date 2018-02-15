library(qtl2pleio)

## ------------------------------------------------------------------------
library(stringr)
library(dplyr)
library(readr)
##First read in the arguments listed at the command line

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
print(args)
##args is now a list of character vectors
## First check to see if arguments are passed.
## Then cycle through each element of the list and evaluate the expressions.
#if(length(args)==0){
#  print("No arguments supplied.")
#}else{
#  for(i in 1:length(args)){
#    eval(parse(text=args[[i]]))
#  }
#}
print(args$argname)
proc_num <- as.numeric(args$argname)
print(proc_num)
#nsnp <- as.numeric(nsnp)
#print(nsnp)
#s1 <- as.numeric(s1)
#print(s1)
run_num <- as.numeric(args$run_num)
print(run_num)
chr <- as.integer(args$chr)
print(chr)
(peak1 <- as.numeric(args$peak1))
(peak2 <- as.numeric(args$peak2))
(probs_file <- args$probs_file)
(phe1_name <- args$phe1_name)
(phe2_name <- args$phe2_name)
###############


## ------------------------------------------------------------------------
load("data/pheno_clin_v6.RData")
## ------------------------------------------------------------------------
phenames <- c(phe1_name, phe2_name)
clin_phe <- pheno_clin[, colnames(pheno_clin) %in% phenames]

## ------------------------------------------------------------------------
# load allele probs file for appropriate chromosome
pp_fn <- paste0("attie_DO500_genoprobs_v5_chr", chr, ".rds")
readRDS(pp_fn) -> pp
# split mouse ids to remove "DO"
#str_split(rownames(clin_phe), pattern = "DO") -> splitted
#sapply(FUN = function(x)x[2], splitted) %>% as.numeric() -> phe_id
#rownames(clin_phe) <- phe_id
# subset genotypes
pp2 <- pp[rownames(pp) %in% rownames(clin_phe), , ]
dim(pp2)
clin_phe <- clin_phe[rownames(clin_phe) %in% rownames(pp2), ]
dim(clin_phe)
if (!check_dimnames(pp2, clin_phe)){
  arrange_by_rownames(pp2, clin_phe) -> pp2
}
check_dimnames(pp2, clin_phe)

## ------------------------------------------------------------------------
readRDS("data/kinship_loco_v5.rds") -> K
K[[chr]] -> kinship
rownames(kinship)
k2 <- kinship[rownames(kinship) %in% rownames(pp2), rownames(kinship) %in% rownames(pp2)]
arrange_by_rownames(k2, pp2) -> k2
## ------------------------------------------------------------------------
readRDS("data/grid_pmap.rds") -> pmap
pm <- pmap[[chr]]
which(pm > peak1)[1] -> p1
which(pm > peak2)[1] -> p2

(start_index <- max(min(p1, p2) - 50, 1))
(stop_index <- min(dim(pp)[3], max(p1, p2) + 50))
(n_snp <- stop_index - start_index + 1)

## ------------------------------------------------------------------------

## ---- cache = TRUE-------------------------------------------------------
# get covariates
load("data/covar.RData") # rownames DO NOT HAVE prefix "DO"
rownames(covar) <- rownames(clin_phe)
covariates <- covar[ rownames(covar) %in% rownames(pp2), c(1, 2, 8:11)]
arrange_by_rownames(covariates, pp2) -> covariates
dim(covariates)
check_dimnames(covariates, pp2) -> indicator
if (!indicator) stop()
# run scan
scan_out <- scan_pvl(probs = pp2, pheno = clin_phe, covariates = covariates, kinship = k2, start_snp1 = start_index, n_snp = n_snp)
## ------------------------------------------------------------------------
fn_out <- paste0("pvl-run", run_num, "_", proc_num, "_", paste(phenames, collapse = "_"), ".txt")
#save(list = "out", file = fn)
write.table(scan_out, fn_out, quote = FALSE)
q("no")
