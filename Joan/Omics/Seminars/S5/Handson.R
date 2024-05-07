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

 
#3.2
 table(se$cohort) 

 table(table(se$subject_id))

 
 se <- se[, !duplicated(se$subject_id)]  # keep only the first (arbitrary)
 table(se$cohort)

 class(se$age) 

  se$age <- gsub(">89", "90", se$age)
  se$age <- as.numeric(as.character(se$age))
 class(se$age) 
 
 
 
 # Remove/modify some characters in the "cohort" and "race" columns
 cols <- c("cohort", "race")
 colData(se)[, cols] <- apply(colData(se)[, cols], 2, function(x){
   y <- gsub("-", "", x)   # Remove "-" 
   z <- gsub(" ", "_", y)  # Change " " (space) by "_"
   w <- gsub("/", "_", z)  # Change "/" by "_"
   return(w)
 })
   
 # Rename columns to individual id
 colnames(se) <- se$subject_id
 
 # Remove uninteresting columns
 colData(se)[, c("subject_id", "time_since_onset", "hospitalized")]  <- NULL 
   
 # Check out the SummarizedExperiment object:
 se 

 colData(se)

 
 se <- se[, se$cohort %in% c("Bacterial", "Influenza")]
 table(se$cohort) 

 
 se0 <- se  

 #3.4
 
 dgl <- calcNormFactors(se, method = "TMM") 
 dgl

 
 
 
  # Calculate CPM (Counts Per Million Reads)
  assays(se)$CPM <- cpm(se) 
 # Calculate CPM based on effective library size (using TMM normalization factors)
   assays(se)$TMM <- cpm(dgl, normalized.lib.sizes = TRUE)

boxplot(log2(assays(se)$counts + 1), las = 2) # Raw read counts


boxplot(log2(assays(se)$TMM + 1), las = 2) # CPM with efective library sizes computed by TMM

#4
#4.1

PCA <- prcomp(log2(t(assays(se)$TMM)+1))
 fviz_pca_ind(PCA, addEllipses = T,
              col.ind = se$cohort,
              # try others e.g. se$race, se$gender, se$age
              pointsize = 3) 

 #4.2
 
 pheatmap(log2(assays(se)$TMM + 1),
             show_rownames = FALSE, 
             annotation_col = as.data.frame(colData(se)))

 # Suppose you want to select only certain rows or columns from your data before creating the heatmap
 
 # Let's say you want to select the first 100 rows and all columns:
 subset_data <- assays(se)$TMM[1:10000, ]
 
 # Then you can apply the log transformation and add 1:
 subset_data_log <- log2(subset_data + 1)
 
 # Now you can create the heatmap with the subset of data:
 pheatmap(subset_data_log, show_rownames = FALSE, annotation_col = as.data.frame(colData(se)))
 

 #4.3
 
 se <- se0  # Initial (preprocessed) dataset
 se <- se[, colnames(se) != "896282"]

 keep <- filterByExpr(se, group = se$cohort)
 se <- se[keep, ]
 dgl <- calcNormFactors(se, method = "TMM") 
 assays(se)$CPM <- cpm(se) 
 assays(se)$TMM <- cpm(dgl, normalized.lib.sizes = TRUE)
 se

 #5
 
 cohort <- as.factor(se$cohort)
 design <- model.matrix(~0 + cohort)
 colnames(design) <- levels(cohort)
 design

 
 race <- as.factor(se$race)
 age <- se$age
 design <- model.matrix(~0 + cohort + age + race)
 colnames(design)[1:2] <- levels(cohort)
 design 

 
 contr.matrix <- makeContrasts(
   Bacterial - Influenza, 
   levels = colnames(design)
 )
 contr.matrix 
 