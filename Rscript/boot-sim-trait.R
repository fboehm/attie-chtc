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
## process index from command line
as.numeric(index) -> snp_index


## ------------------------------------------------------------------------
library(dplyr)
library(qtl2pleio)
library(gemma2)
library(pleiotropy)
# load data
phe_nona <- as.matrix(read.csv("sim_data/2017-10-11-sim_pleio.csv"))
source("Rscript/setup.R")
source("Rscript/setup-chr2a.R")
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
## set start POINT HERE
snp1 <- 175

# get Bhat
Sigma <- kinship %x% Vg + diag(n_mouse) %x% Ve
X1 <- pp3[ , , snp_index]
X <- pleiotropy::stagger_mats(X1, X1)
phe_nona_vec <- as.vector(phe_nona)
Bhat <- calc_Bhat(X = X, Y = as.matrix(phe_nona_vec), Sigma = Sigma)
B <- matrix(data = Bhat, byrow = FALSE, ncol = 2)
system.time(pvl_boot(X = X1, B = B, Vg_initial = Vg, Ve_initial = Ve, 
                     kinship = kinship, probs = pp3, start_snp = snp1, 
                     n_snp = 50, nboot = 1) -> boot_out
)

fn <- paste0("boot-out_", proc_num, ".csv")

write.csv(boot_out, fn)

q("no")



