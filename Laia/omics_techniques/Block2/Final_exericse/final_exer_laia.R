# Load packages
library(edgeR)
library(Glimma)
library(SummarizedExperiment)
library(factoextra)
library(pheatmap)
library(EnsDb.Hsapiens.v86)
library(GO.db)
library(org.Hs.eg.db)

################################################################################
# DATA
################################################################################

# Load the data
files_to_load <- c("GSE161731_counts.csv.gz", "GSE161731_counts_key.csv.gz")

# Download them
for (f in files_to_load) {
  url <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE161nnn/GSE161731/suppl/", f)
  download.file(url, destfile = f)
}

# Load the files in R
counts <- read.csv("GSE161731_counts.csv.gz", 
                   header = TRUE, check.names = FALSE, row.names = 1)
metadata <- read.csv("GSE161731_counts_key.csv.gz",
                     header = TRUE, check.names = FALSE, row.names = 1)

# Obtain gene annotations
genes <- genes(EnsDb.Hsapiens.v86)

# Build the SUMMARIZED EXPERIMENT object
comm.s <- intersect(colnames(counts), rownames(metadata))       # Get common samples (in 'counts' and 'metadata')
comm.g <- intersect(rownames(counts), genes$gene_id)            # Get common genes (in 'counts' and the annotation)
se <- SummarizedExperiment(assay = list("counts" = counts[comm.g, comm.s]),
                           colData = metadata[comm.s, ],
                           rowRanges = genes[comm.g])
se

