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
### define nn
nn <- 3

i <- proc_num %/% nn + 1
j <- proc_num %% nn + 1
print(i)
print(j)

## ------------------------------------------------------------------------

library(regress)
library(pleiotropy)
library(dplyr)
# load packages
library(broman) # contains the winsorize() function

# load data


### genotype probabilities ("probs") in form used by R/qtl2
### see ../R/0_DOQTLprobs2qtl2.R for how they were converted
PATH_TO_DERIVED_DATA <- "data"
#PATH_TO_DERIVED_DATA <- "~/Box Sync/attie/attiedo"
#PATH_TO_DERIVED_DATA <- "~/attie"
#load(file.path(PATH_TO_DERIVED_DATA, "DerivedData/GM_Attie_allele_call_haploprobs_4qtl2_wave5.Rdata"))
load("probs_17.RData") 
. -> pp
### clinical phenotypes + phenotype dictionary
### ("pheno_clin" and "pheno_clin_dict")
load(file.path(PATH_TO_DERIVED_DATA, "pheno_clin.RData"))
### kinship matrices ("loco" method) ("K")
load(file.path(PATH_TO_DERIVED_DATA, "kinship.RData"))
### covariate matrix ("covar")
#load(file.path(PATH_TO_DERIVED_DATA, "DerivedData/covar.RData"))
### physical map of the markers in the probs array
load(file.path(PATH_TO_DERIVED_DATA, "probs_pmap.RData"))

## ------------------------------------------------------------------------

pm <- pmap$`17`
kinship <- K$`17`
## Determine which markers are shared between pmap & gmap
snp_g <- dimnames(pp)[[3]]
snp_p <- names(pm)
shared_snps <- intersect(snp_g, snp_p)
pp2 <- pp[ , , snp_g %in% shared_snps]
pm2 <- pm[snp_p %in% shared_snps]

## ------------------------------------------------------------------------
phenames <- c("Insulin at sac", "oGTT weight")
(g <- grep(phenames[1], pheno_clin_dict$name))
pheno_clin_dict[g, ]
foo <- pheno_clin[ , pheno_clin_dict[g, "short_name"]]
phe1 <- log(broman::winsorize(foo)) # take logs and winsorize (rather heavily!)
(g <- grep(phenames[2], pheno_clin_dict$name))
g <- 78 # use the first value here, after checking the two names!
pheno_clin_dict[g, ]
bar <- pheno_clin[ , pheno_clin_dict[g, "short_name"]]
phe2 <- log(broman::winsorize(bar)) 
# match phenotypes & genotypes ordering
foo <- cbind(phe1, phe2)
rownames(foo) <- rownames(pheno_clin)
rownames(foo) -> rn
match_out <- match(rownames(pp2), rn)
phe <- foo[match_out, ]
rownames(phe) == rownames(pp2)
# Remove NAs from phe
## first, make an indicator for missingness
missing <- is.na(phe[ , 1]) | is.na(phe[ , 2])
## now, filter phe to remove those with a missing phenotype value
phe_nona <- phe[!missing, ]
print(dim(phe_nona))
## filter kinship to remove those with a missing phenotype value

kinship_nona <- kinship[!missing, !missing]
print(dim(kinship_nona))
## filter pp2 to remove those with a missing phenotype
pp2_nona <- pp2[!missing, , ]
print(dim(pp2_nona))
## ------------------------------------------------------------------------
start_snp <- which(pm2 > 31.1)[1]
#stop_snp <- which(pm2 > 150)[1]
start_snp_i <- (i - 1)* 10 + 1 + start_snp 
stop_snp_i <- i * 10 + start_snp

start_snp_j <- (j - 1)* 10 + 1 + start_snp
stop_snp_j <- j * 10 + start_snp

## ------------------------------------------------------------------------
n_snp <- stop_snp_i - start_snp_i + 1
loglik <- matrix(NA, nrow = n_snp, ncol = n_snp)
#reg_obj <- list()
#reg_obj_i <- list()
i_count <- 0
j_count <- 0
for (i in start_snp_i:stop_snp_i){
  i_count <- i_count + 1
  print(i)
  for (j in start_snp_j:stop_snp_j){
    j_count <- j_count + 1
    out <- fit1_bvlmm(Y = phe_nona, 
                                  X1 = pp2_nona[ , , i], 
                                  X2 = pp2_nona[ , , j], 
                                  Kmat = kinship_nona
                                    )
    
    assemble_Ve(out$sigma) -> Ve
    assemble_Vg(out$sigma) -> Vg
    ll <- calc_ll_bvlmm(yvec = as.vector(phe_nona), 
                        Xmat = stagger_mats(pp2_nona[ , , i], pp2_nona[ , , j]),
                        betahat = out$beta, 
                        Vg = Vg, 
                        Ve = Ve,
                        Kmat = kinship_nona)
    print(c(i_count, j_count))
    ll -> loglik[i_count, ((j_count - 1) %% 10) + 1]
   # out -> reg_obj_i[[j - start_snp + 1]]
  }
  #reg_obj[[i - start_snp + 1]] <- reg_obj_i
}

save("loglik", file = paste0("out_", i, "_", j, ".RData"))
q("no")
