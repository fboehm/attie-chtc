#pp <- probs$`2`
pp <- probs

pm <- pmap$`2`
kinship <- K$`2`

## Determine which markers are shared between pmap & gmap
snp_g <- dimnames(pp)[[3]]
snp_p <- names(pm)
shared_snps <- intersect(snp_g, snp_p)
pp3 <- pp[ , , snp_g %in% shared_snps]
pm2 <- pm[snp_p %in% shared_snps]


#Now, we get the phenotype data.

#We now make a phenotypes matrix without missing values.

n <- nrow(kinship)
