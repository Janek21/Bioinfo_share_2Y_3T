library(dplyr)
library(gplots)
library(pheatmap)
library(plyr)

setwd('/home/laia/Desktop/omics_techniques/Roadmap')
source('BEDtoolsR.R')

metadata <- read.table("metadata.roadmap_clean.txt", sep = "\t")
View(metadata)

# See if there is a diff between fetal and adult tissue
# Make an extra column in which we store it the tissue is fetal or adult.
metadata$PERIOD <-  ifelse(grepl('fetal|Fetal', metadata$V3), "Fetal", "Adult")

# Store in another column the name of the shortname and the period
metadata$Name <- paste0(metadata$V2,'_', metadata$PERIOD)

# we can ave diff samples of the same tissue and period, we need a unique identifier, we'll have a descriptive unique identifier
metadata$NameUnique <- paste0(metadata$Name,'_', metadata$V1)

###############################################################################

# We are only working on tissues brain, muscle, heart, digestive
tissues <- c('Brain', 'Muscle', 'Heart', 'Digestive')
samples <- filter(metadata, V2 %in% tissues)$V1

tissues <- c('Brain','Muscle','Heart','Digestive')
samples <- dplyr::filter(metadata, V2 %in% tissues)$V1

# Store all the information
# Create list
roadmap <- list()

# Read the file 
for (f in samples){
  if (file.exists(paste0("all.mnemonics.bedFiles/", f,"_18_core_K27ac_mnemonics.bed.gz"))){
    print(f)
    roadmap[[f]] <- read.table(gzfile(paste0("all.mnemonics.bedFiles//", f,"_18_core_K27ac_mnemonics.bed.gz")))
  }
}

# change the names. Assign each name
names(roadmap) <- mapvalues(names(roadmap), from=metadata$V1, to=metadata$NameUnique)

roadmap[["Heart_Adult_E065"]] %>% head()

# compute the intersection
bedTools.2in(bed1 = roadmap[["Heart_Adult_E065"]], bed2 = roadmap[["Digestive_Adult_E101"]]) -> INT
head(INT)

# compute the intersection
state <- "1_TssA"
bed1 <- filter(roadmap[["Heart_Adult_E065"]], V4==state)
bed2 <- filter(roadmap[["Digestive_Adult_E101"]], V4==state)
int <- bedTools.2in(bed1 = bed1, bed2 = bed2)

head(bed1)

calculateJaccard <- function(a,b){
  intersection <- bedTools.2in(bed1 = a, bed2 = b)
  union <- 
  
}


m <- matrix(nrow = length(names(roadmap)), ncol = length(names(roadmap)), dimnames = list(names(roadmap), names(roadmap)))
for (i in names(roadmap)){
  for(j in names(roadmap)){
    print(c(i,j))
    a <- filter(roadmap[[i]], V4==state)
    b <- filter(roadmap[[j]], V4==state)
    z <- calculateJaccard(a,b)
    m[i,j] <- z
  }
}









