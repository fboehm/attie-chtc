pp <- probs
pm <- pmap$`17`
kinship <- K$`17`

## Determine which markers are shared between pmap & gmap
snp_g <- dimnames(pp)[[3]]
snp_p <- names(pm)
shared_snps <- intersect(snp_g, snp_p)
pp2 <- pp[ , , snp_g %in% shared_snps]
pm2 <- pm[snp_p %in% shared_snps]


#Now, we get the phenotype data.

phenames <- c("G8.3+GLP1 secretion", "weight change for week 11 vs week 1")
(g <- 59) # not sure why phenames[1] seems to not work with grep
pheno_clin_dict[g, ]
foo <- pheno_clin[ , pheno_clin_dict[g, "short_name"]]
phe1 <- log(broman::winsorize(foo)) # take logs and winsorize (rather heavily!)
(g <- grep(phenames[2], pheno_clin_dict$name))
# use the first value here, after checking the two names!
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
p1 <- as.vector(phe[, 1])
rownames(phe) -> names(p1)
p2 <- as.vector(phe[, 2])
rownames(phe) -> names(p2)

#We now make a phenotypes matrix without missing values.

na_ind <- is.na(phe[, 1]) | is.na(phe[, 2])
phe_nona <- phe[!na_ind, ]
pp3 <- pp2[!na_ind, , ]
kinship_nona <- kinship[!na_ind, !na_ind]
