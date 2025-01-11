#!/usr/bin/env Rscript --no-init-file

# Load required libraries
suppressMessages(library(Matrix))
suppressMessages(library(argparse))

# Argument parser setup
parser <- ArgumentParser(description = "Convert .rds file to Matrix Market (.mtx) format")
parser$add_argument("input_rds", type = "character", help = "Input .rds file")
parser$add_argument("output_prefix", type = "character", help = "Output prefix for .mtx, genes, and barcodes files")
args <- parser$parse_args()

# Read the input RDS file
input_rds <- args$input_rds
output_prefix <- args$output_prefix

cat("Reading RDS file:", input_rds, "\n")
data <- readRDS(input_rds)

# Validate data type
if (!inherits(data, "dgCMatrix")) {
  stop("The input RDS file does not contain a dgCMatrix object.")
}

# Extract matrix, genes, and barcodes
cat("Extracting matrix, genes, and barcodes...\n")
sparse_matrix <- data
sparse_matrix <- t(sparse_matrix)
genes <- colnames(sparse_matrix)
barcodes <- rownames(sparse_matrix)

if (is.null(genes) || is.null(barcodes)) {
  stop("Row or column names are missing in the matrix.")
}

# Write Matrix Market file (.mtx)
cat("Writing Matrix Market file...\n")
mtx_file <- paste0(output_prefix, ".mtx")
writeMM(sparse_matrix, file = mtx_file)

# Write genes file
cat("Writing genes file...\n")
genes_file <- paste0(output_prefix, "_genes.tsv")
write.table(genes, file = genes_file, quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)

# Write barcodes file
cat("Writing barcodes file...\n")
barcodes_file <- paste0(output_prefix, "_barcodes.tsv")
write.table(barcodes, file = barcodes_file, quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)

cat("Conversion complete.\n")

