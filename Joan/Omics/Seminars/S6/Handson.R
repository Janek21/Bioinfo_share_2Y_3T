library(dplyr)
library(gplots)
library(pheatmap)

source("BEDtoolsR.R")

setwd("/home/joan/Desktop/Link to Joan/Omics/Seminars/S6")
metadata<- read.table("metadata.roadmap_clean.txt", sep = "\t")

metadata <- read.delim("metadata.roadmap_clean.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(metadata) <- c("ID", "Tissue", "Description")

tissue <- c("Brain", "Muscle", "Heart", "Digestive")
samples <- metadata[grep(paste(tissue, collapse = "|"), metadata$Tissue), ]

metadata$PERIOD <- ifelse(grepl("Fetal", metadata$Description), "Fetal", "Non-fetal")

roadmap <- list()

for(f in samlpes){
  if(file.exists(paste0("all.mnemonics.bedFiles/",f,"_18_core_K27ac_mnemonics.bed.gz"))){
    print(f)
    roadmap[[f]]<- read.table(gzfile(paste0("all.mnemonics.bedFiles/",f,"_18_core_K27ac_mnemonics.bed.gz")))
  }
}


names(roadmap) <- mapvalues(names(roadmap), from = selected_metadata$ID, to= selected_metadata$UniqueName )

state="TssA"
bed1 <- filter(roadmap[["Muscle_Fetal_B090"]], V4==state)
bed2 <- filter(roadmap[["Digestive_Adult_E102"]], V4==state)

bedTools.2

bedTools.2in(bed1 = filter(roadmap[["Muscle_Fetal_B090"]], )