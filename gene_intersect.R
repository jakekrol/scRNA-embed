# get intersecting genes
f <- list.files(pattern="rds$")
n <- length(f)
f_1 <- f[1]
get_genes <- function(f) {
    d <- readRDS(f)
    genes <- rownames(d)
    return(genes)
}
genes <- get_genes(f_1)
print(paste("Num genes", length(genes)))
f <- f[-1]
for (file in f) {
    d <- readRDS(file)
    g_d <- rownames(d)
    genes <- intersect(genes, g_d)
    print(paste("Num genes", length(genes)))
}
writeLines(genes, "intersecting_genes.txt")
