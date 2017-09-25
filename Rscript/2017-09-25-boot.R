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
set.seed(proc_num)
## ------------------------------------------------------------------------
library(dplyr)
library(qtl2pleio)
library(gemma2)
library(pleiotropy)
# load data

source("Rscript/setup.R")
source("Rscript/setup-chr5.R")
eigen2(kinship) -> e_out
e_out$vectors -> U
e_out$values -> eval
## ------------------------------------------------------------------------
n_mouse <- nrow(pp3)

X1pre <- rep(1, n_mouse) %>% as.matrix() %>% t()
X1 <- X1pre %*% U
Y <- t(phe_nona) %*% U

MphEM(Y = Y, X = X1, eval = eval, V_g = diag(2), V_e = diag(2)) -> g_out
g_out[[length(g_out)]][[2]] -> Vg
g_out[[length(g_out)]][[3]] -> Ve
snp1 <- which(pm2 > 147)[1]
scan_pvl(probs = pp3, 
         pheno = phe_nona, 
         kinship = kinship, 
         start_snp1 = snp1, 
         start_snp2 = snp1, 
         n_snp = 250
         ) -> foo1
tidy_scan_pvl(foo1, pmap$`5`) -> bar


bar2 <- bar %>%
  filter(trace == "pleio")
index <- which.max(bar2$lod) + snp1

# get Bhat
Sigma <- kinship %x% Vg + diag(n_mouse) %x% Ve
X1 <- pp3[ , , index]
X <- pleiotropy::stagger_mats(X1, X1)
Bhat <- calc_Bhat(X = X, Y = phe_nona, Sigma = Sigma)
B <- matrix(data = Bhat, byrow = FALSE, ncol = 2)
system.time(pvl_boot(X = X1, B = B, Vg_initial = Vg, Ve_initial = Ve, 
                     kinship = kinship, probs = pp3, start_snp = snp1, 
                     n_snp = 250, nboot = 1) -> boot_out
)

fn <- paste0("boot-out_", proc_num, ".csv")

write.csv(boot_out, fn)

q("no")



