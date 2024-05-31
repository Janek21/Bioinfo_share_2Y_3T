# RNA-SEQ DATA ANALYSIS

################################################################################
#   R PACKAGES
################################################################################
BiocManager::install("edgeR")
BiocManager::install("Glimma")
BiocManager::install("SummarizedExperiment")
BiocManager::install("factoextra")
BiocManager::install("pheatmap")
BiocManager::install("GO.db")
BiocManager::install("org.Hs.eg.db")

library(edgeR)
library(Glimma)
library(SummarizedExperiment)
library(factoextra)
library(pheatmap)
library(EnsDb.Hsapiens.v86)
library(GO.db)
library(org.Hs.eg.db)

################################################################################
#   DATA
################################################################################

# DOWNLOAD RNA-SEQ DATA FROM GEO AND LOAD IT INTO R
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

# DATA CLEANING
# - Remove duplicates
# Count samples per condition
table(se$cohort)  

# Some individuals have replicated measurements
table(table(se$subject_id))

se <- se[, !duplicated(se$subject_id)]  # keep only the first (arbitrary)
table(se$cohort)

# - Ensure that the variables are the appropriate class
# Convert "age" to numeric
class(se$age)

se$age <- gsub(">89", "90", se$age)
se$age <- as.numeric(as.character(se$age))
class(se$age)

# - Metadata cleanup
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

# - Subset the conditions of interest --> select bacterial and influenza
se <- se[, se$cohort %in% c("Bacterial", "Influenza")]
table(se$cohort)

se0 <- se  # Make a copy of the SummarizedExperiment for later usage

# REMOVAL OF LOWLY EXPRESSED GENES
keep <- filterByExpr(se, group = se$cohort) # Check out ?filterByExpr
table(keep)

se <- se[keep, ]

# DATA NORMALIZATION
dgl <- calcNormFactors(se, method = "TMM") 
dgl

# Calculate CPM (Counts Per Million Reads)
assays(se)$CPM <- cpm(se) 
 
# Calculate CPM based on effective library size (using TMM normalization factors)
assays(se)$TMM <- cpm(dgl, normalized.lib.sizes = TRUE)

boxplot(log2(assays(se)$counts + 1), las = 2) # Raw read counts
boxplot(log2(assays(se)$CPM + 1), las = 2) # CPM

boxplot(log2(assays(se)$TMM + 1), las = 2) # CPM with efFective library sizes computed by TMM

################################################################################
#   EXPLORATORY DATA ANALYSIS
################################################################################
# PRINCIPAL COMPONENT ANALYSIS
# Perform PCA analysis on log-transformed TMM normalized counts
PCA <- prcomp(log2(t(assays(se)$TMM) + 1))

# Visualize PCA results with ellipses for groups
fviz_pca_ind(PCA, addEllipses = TRUE,
             col.ind = se$cohort,  # Color by 'cohort', can try 'se$race', 'se$gender', 'se$age'
             pointsize = 3)

# CLUSTERING / HEATMAP
# Generate a heatmap of log-transformed TMM normalized counts
pheatmap(log2(assays(se)$TMM + 1),
         show_rownames = FALSE, 
         annotation_col = as.data.frame(colData(se)))

# OUTLIER REMOVAL
se <- se0  # Initial (preprocessed) dataset
se <- se[, colnames(se) != "896282"]  

# Filter the SummarizedExperiment object to keep only highly expressed genes
keep <- filterByExpr(se, group = se$cohort)
se <- se[keep, ]

# Calculate normalization factors using the TMM method
dgl <- calcNormFactors(se, method = "TMM")

# Calculate and store Counts Per Million (CPM) for both unnormalized and normalized data
assays(se)$CPM <- cpm(se)
assays(se)$TMM <- cpm(dgl, normalized.lib.sizes = TRUE)

# Display the updated SummarizedExperiment object
se

################################################################################
#   DIFFERENTIAL GENE EXPRESSION ANALYSIS
################################################################################
# CREATING A DESIGN MATRIX AND CONTRASTS
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

# REMOVING MEAN-VARIANCE RELATIONSHIP FROM COUNT DATA
y <- voom(dgl, design, plot = TRUE)  # voom uses the TMM normalization factors available from 'dgl'

# FITTING LINEAR MODELS FOR COMPARISONS OF INTEREST
fit <- lmFit(y, design)

# Apply the specified contrasts to the fitted model
fit <- contrasts.fit(fit, contr.matrix)

# Apply Empirical Bayes method to shrink the standard errors
fit <- eBayes(fit)

# EXAMINING THE NUMBER OF DE GENES
results <- topTable(fit, sort.by = "P", n = Inf) # Get all results, ordered by P-value
head(results)

dim(results)

results.5fdr <-  subset(results, adj.P.Val < 0.05) 
nrow(results.5fdr)

dt <- decideTests(fit, adjust.method = "fdr", p.value = 0.05)
summary(dt)

# USEFUL GRAPHICAL REPRESENTATIONS OF DIFFERENTIAL EXPRESSION RESULTS
glMDPlot(fit, status=dt, side.main="ENSEMBL", counts = log2(assays(se)$TMM+1), groups = se$cohort)

dgl$sampls$group <- dgl$samples$cohort  # Allows coloring by cohort
glimmaVolcano(fit, dge = dgl)

# Create a data frame with cohort data and log-transformed expression values
df <- data.frame(colData(se), expression = log2(assays(se)$TMM["ENSG00000111335", ] + 1))

# Plot the data using ggplot2
ggplot(df, aes(x = cohort, y = expression, fill = cohort)) +
  geom_boxplot() +
  ylab("log2(CPM+1)") +
  ggtitle("Expression of ENSG00000111335 vs 'cohort'") +
  theme(plot.title = element_text(hjust = 0.5))


# Create a heatmap using pheatmap with selected genes and annotations
pheatmap(log2(assays(se)$TMM[rownames(results.5fdr), ] + 1),
         show_rownames = FALSE,
         annotation_col = as.data.frame(colData(se)))

################################################################################
#   BIOLOGICAL ANALYISIS
################################################################################
# GO ENRICHMENT ANALYSIS
dgl2 <- dgl
dgl2 <- dgl2[!is.na(dgl2$genes$entrezid), ]
dgl2 <- dgl2[!duplicated(dgl2$genes$entrezid), ]
dgl2 <- dgl2[!grepl(" ", dgl2$genes$entrezid), ]
rownames(dgl2) <- as.character(dgl2$genes$entrezid)

y2 <- voom(dgl2, design, plot = FALSE)
fit2 <- lmFit(y2, design)
fit2 <- contrasts.fit(fit2, contr.matrix)
fit2 <- eBayes(fit2)

go <- goana(fit2, species="Hs")
topGO(go, n=10)

# KEGG PATHWAY ANALYSIS
keg <- kegga(fit2, species="Hs")
topKEGG(keg, n=10)


