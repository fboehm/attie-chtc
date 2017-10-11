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
### define nn -- CHANGE THIS AS NEEDED
nn <- 15
## nn should be set to sqrt of total number of jobs.
## for 400 jobs, put nn = 20
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
load("probs_2.RData") 
. -> pp
### clinical phenotypes + phenotype dictionary
### ("pheno_clin" and "pheno_clin_dict")
#load(file.path(PATH_TO_DERIVED_DATA, "pheno_clin.RData"))
### kinship matrices ("loco" method) ("K")
load(file.path(PATH_TO_DERIVED_DATA, "kinship.RData"))
### covariate matrix ("covar")
#load(file.path(PATH_TO_DERIVED_DATA, "DerivedData/covar.RData"))
### physical map of the markers in the probs array
load(file.path(PATH_TO_DERIVED_DATA, "probs_pmap.RData"))
# read csv file containing simulated phenotype
library(readr)
read_csv("sim_data/2017-10-05-sim_link20.csv") -> phe
phe <- as.matrix(phe)


## ------------------------------------------------------------------------

pm <- pmap$`2`
kinship <- K$`2`
## Determine which markers are shared between pmap & gmap
snp_g <- dimnames(pp)[[3]]
snp_p <- names(pm)
shared_snps <- intersect(snp_g, snp_p)
pp2 <- pp[ , , snp_g %in% shared_snps]
pm2 <- pm[snp_p %in% shared_snps]

## ------------------------------------------------------------------------
phenames <- c("t1", "t2")

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

snp1 <- 180
#stop_snp <- which(pm2 > 150)[1]
start_snp_i <- (i - 1)* 10 + snp1 
start_snp_j <- (j - 1)* 10 + snp1

## ------------------------------------------------------------------------

library(qtl2pleio)
scan_pvl(probs = pp2_nona, 
         pheno = phe_nona, 
         kinship = kinship_nona, 
         start_snp1 = start_snp_i, 
         start_snp2 = start_snp_j, 
         n_snp = 10
         ) -> foo
fn <- paste0("scan_pvl-out-", i, "-", j, ".csv")
write.csv(foo, fn)

q("no")


