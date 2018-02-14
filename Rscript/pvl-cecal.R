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
peak1 <- as.numeric(args$peak1)
peak2 <- as.numeric(args$peak2)
###############

## ------------------------------------------------------------------------
#mytab <- read_tsv("~/Box Sync/attie/traeger/annotation/toFred_summedKO_cecalLipids_mediation_reformat_restrictedToCommonMice_20180211.txt") %>%
#  select(- phenotype_2_chr) %>%
#  rename(chr = phenotype_1_chr)

## ------------------------------------------------------------------------
cecal <- readRDS("data/attie_cecum_lipids_zscore_normalized.lltmod.v2.rds")

## ------------------------------------------------------------------------
summed_ko <- readRDS("summedKO_DO_finalRSEM_log2transformed_20171212.rds")

## ------------------------------------------------------------------------
phenames <- c(phe1, phe2)
cecal_phe <- cecal[, colnames(cecal) %in% phenames]
names(cecal_phe) <- rownames(cecal)
ko_phe <- summed_ko[, colnames(summed_ko) %in% phenames]
names(ko_phe) <- rownames(summed_ko)

## ------------------------------------------------------------------------
names(ko_phe) == names(cecal_phe)
intersect(names(ko_phe), names(cecal_phe)) -> shared_ids
ko_phe[names(ko_phe) %in% shared_ids] -> ko2
cecal_phe[names(cecal_phe) %in% shared_ids] -> ce2
ce2[match(x = names(ko2), table = names(ce2))] -> ce3
phe <- tibble(mouse_id = names(ce3), cecal = ce3, ko = ko2)

## ------------------------------------------------------------------------
pp_fn <- paste0("probs_", chr, ".RData")
load(pp_fn)
. -> pp
#rownames(probs$`2`)
str_split(phe$mouse_id, pattern = "DO") -> splitted
sapply(FUN = function(x)x[2], splitted) %>% as.numeric() -> phe_id
# subset genotypes data
#pp <- probs$`2`
pp2 <- pp[rownames(pp) %in% phe_id, , ]
tib2 <- phe %>%
  filter(phe_id %in% rownames(pp2))
t3 <- tib2 %>%
  mutate(mouse_id2 = phe_id[phe_id %in% rownames(pp2)])
# check to see if the orderings match for phe & geno 
t3$mouse_id2 == rownames(pp2)
t3[match(x = rownames(pp2), table = t3$mouse_id2), ] -> t4
t4$mouse_id2 == rownames(pp2)

## ------------------------------------------------------------------------
load("data/kinship_qtl2.RData")
K[[chr]] -> kinship
k2 <- kinship[rownames(kinship) %in% rownames(pp2), rownames(kinship) %in% rownames(pp2)]
rownames(k2) == rownames(pp2)
colnames(k2) == rownames(pp2)
rownames(k2) == t4$mouse_id2

## ------------------------------------------------------------------------
load("data/probs_pmap.RData")
pm <- pmap[[chr]]
which(pm > peak1)[1] -> p1
which(pm > peak2)[1] -> p2

start_index <- min(min(p1, p2) - 50, 1)
stop_index <- min(dim(pp)[3], max(p1, p2) + 50)
n_snp <- stop_index - start_index + 1


## ------------------------------------------------------------------------
library(qtl2pleio)

## ---- cache = TRUE-------------------------------------------------------
pheno <- as.matrix(t4[, 2:3])
rownames(pheno) <- t4$mouse_id2
colnames(pheno) <- colnames(t4)[2:3]
scan_out <- scan_pvl(probs = pp2, pheno = pheno, kinship = k2, start_snp1 = start_index, n_snp = n_snp)
## ------------------------------------------------------------------------
fn_out <- paste0("pvl-run", run_num, "_", proc_num, "_", paste(phenames, collapse = "_"), ".txt")
#save(list = "out", file = fn)
write.table(scan_out, fn_out, quote = FALSE)
q("no")
