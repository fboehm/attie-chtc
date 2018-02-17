library(qtl2pleio)
library(stringr)
library(dplyr)
library(readr)
args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
print(args)

print(args$argname)
proc_num <- as.numeric(args$argname)
print(proc_num)
run_num <- as.numeric(args$run_num)
print(run_num)
chr <- as.integer(args$chr)
print(chr)
(peak1 <- as.numeric(args$peak1))
(peak2 <- as.numeric(args$peak2))
(probs_file <- args$probs_file)
(phe1_name <- args$phe1_name)
(phe2_name <- args$phe2_name)

phenames <- c(phe1_name, phe2_name)

## ------------------------------------------------------------------------
load("data/pheno_clin_v6.RData")
## ------------------------------------------------------------------------
clin_phe1 <- pheno_clin[, colnames(pheno_clin) == phenames[1]]
clin_phe2 <- pheno_clin[, colnames(pheno_clin) == phenames[2]]
cbind(clin_phe1, clin_phe2) -> clin_phe

## ------------------------------------------------------------------------
# load allele probs file for appropriate chromosome
pp_fn <- paste0("attie_DO500_genoprobs_v5_chr", chr, ".rds")
readRDS(pp_fn) -> pp
# subset genotypes
pp2 <- pp[rownames(pp) %in% rownames(clin_phe), , ]
dim(pp2)

clin_phe <- clin_phe[rownames(clin_phe) %in% rownames(pp2), ]
dim(clin_phe)
if (!check_dimnames(pp2, clin_phe)){
  arrange_by_rownames(pp2, clin_phe) -> pp2
}
check_dimnames(pp2, clin_phe)

## ---- read_analyses_file-------------------------------------------------
readr::read_csv("data/analyses_clin.csv") -> anal
library(dplyr)
anal %>%
  filter(pheno == phenames[1]) -> anal_phe1
anal %>%
  filter(pheno == phenames[2]) -> anal_phe2

## ------------------------------------------------------------------------
# apply transformations to clin_phe as needed in analyses_clin.csv
foo <- do.call(anal_phe1$transf[1], list(clin_phe[ , 1] + anal_phe1$offset[1]))
if (anal_phe1$winsorize[1]){
  foo <- broman::winsorize(foo)
}
foo -> clin_phe[ , 1]

foo <- do.call(anal_phe2$transf[1], list(clin_phe[ , 2] + anal_phe2$offset[1]))
if (anal_phe2$winsorize[1]){
  foo <- broman::winsorize(foo)
}
foo -> clin_phe[ , 2]

## ------------------------------------------------------------------------
# define covariates needed
(colnames(anal_phe1)[9:20])[unlist(anal_phe1[1, 9:20])] -> phe1_cov
(colnames(anal_phe2)[9:20])[unlist(anal_phe2[1, 9:20])] -> phe2_cov
cov_names <- union(phe1_cov, phe2_cov)

## ----read_kinship--------------------------------------------------------
# we'll use both sex and diet days, but not wave indicators, below
## ------------------------------------------------------------------------
readRDS("data/kinship_loco_v5.rds") -> K
K[[chr]] -> kinship
rownames(kinship)
k2 <- kinship[rownames(kinship) %in% rownames(pp2), rownames(kinship) %in% rownames(pp2)]
arrange_by_rownames(k2, pp2) -> k2

## ----read_pmap-----------------------------------------------------------
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
if (length(cov_names) > 0){
  (covariates <- covar[ , colnames(covar) %in% cov_names])
  (apply(FUN = function(x)identical(x, rep(x[1], length(x))), X = covariates, MARGIN = 2) -> cov_cols)
# check covariate columns for all entries having the same value
# remove those columns that have all entries being a single value
  covariates <- covariates[ , !cov_cols]
#
  arrange_by_rownames(covariates, pp2) -> covariates
  dim(covariates)
  check_dimnames(covariates, pp2) -> indicator
  if (!indicator) stop()
} else {
  covariates <- NULL
}
## ---- cache = TRUE-------------------------------------------------------
# run scan
scan_out <- scan_pvl(probs = pp2, pheno = clin_phe, covariates = covariates, kinship = k2, start_snp1 = start_index, n_snp = n_snp, max_iter = 10000, max_prec = 0.00001)
## ------------------------------------------------------------------------
fn_out <- paste0("pvl-run", run_num, "_", proc_num, "_", paste(phenames, collapse = "_"), ".txt")
save(list = "out", file = fn)

