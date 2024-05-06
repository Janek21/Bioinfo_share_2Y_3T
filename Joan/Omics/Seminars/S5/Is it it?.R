library(edgeR)
library(EnsDb.Hsapiens.v86)
library(Glimma)
library(SummarizedExperiment)
library(factoextra)
library(pheatmap)

files_to_download <- c("GSE161731_counts.csv.gz",      # Gene expression (count matrix)
                       "GSE161731_counts_key.csv.gz")  # Metadata
# Download them
 for (f in files_to_download) {
   url <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE161nnn/GSE161731/suppl/", f)
   download.file(url, destfile = f)
 }
   
 # Load both files in R
 counts <- read.csv("GSE161731_counts.csv.gz", 
                    header = TRUE, check.names = FALSE, row.names = 1)
 metadata <- read.csv("GSE161731_counts_key.csv.gz",
                      header = TRUE, check.names = FALSE, row.names = 1)
 
 # Obtain gene annotations
 genes <- genes(EnsDb.Hsapiens.v86)

 
  comm.s <- intersect(colnames(counts), rownames(metadata))       # Get common samples (in 'counts' and 'metadata')
  comm.g <- intersect(rownames(counts), genes$gene_id)            # Get common genes (in 'counts' and the annotation)
  se <- SummarizedExperiment(assay = list("counts" = counts[comm.g, comm.s]),
                              colData = metadata[comm.s, ],
                              rowRanges = genes[comm.g])
  
   se 
   colData(se)   
   
   
   rowRanges(se)
   
   
   
   #3.2
   
   table(se$cohort)
   
   table(table(se$subject_id))
   
   se <- se[, !duplicated(se$subject_id)]
table(se$cohort)   


class(se$age)
