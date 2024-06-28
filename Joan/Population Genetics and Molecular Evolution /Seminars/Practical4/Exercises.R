efetch("NT_033777", db="nucleotide", rettype="fasta_cds_na", retmode="text", strand=2, seqstart=5976177, seqstop=5976817, outfile="./Dmel_obp83b_cds.fasta")
efetch("NC_052523", db="nucleotide", rettype="fasta_cds_na", retmode="text", strand=2, seqstart=2470507, seqstop=2471257, outfile="./Dsim_obp83b_cds.fasta")
efetch("NW_020825194", db="nucleotide", rettype="fasta_cds_na", retmode="text", strand=1, seqstart=26500029, seqstop=26500718, outfile="./Dere_obp83b_cds.fasta")
efetch("NW_024545863", db="nucleotide", rettype="fasta_cds_na", retmode="text", strand=2, seqstart=12746083, seqstop=12746802, outfile="./Dele_obp83b_cds.fasta")
efetch("NC_046679", db="nucleotide", rettype="fasta_cds_na", retmode="text", strand=1, seqstart=16336491, seqstop=16337172, outfile="./Dpseu_obp83b_cds.fasta")

library(ape)

## read the sequences
myseqs<-read.dna("./obp83b_cds.fasta", format="fasta")

## visualize the sequences
myseqs

## multiple sequence alignment
myalign <- clustal(myseqs)

## partial visualization of the generated alignment
myalign

## write a file (maintaining the gaps introduced in the msa step)
write.dna(myalign, file="./obp83b_cds_aln.fasta", format="fasta")


library(seqinr)

translate<-seqinr::translate

## function to translate a codon:
dna_to_aa <- function(codon) {
  dna<-s2c(codon)
  aa<-translate(dna, frame = 0, sens = "F", numcode = 1, NAstring = "X", ambiguous = FALSE)
  return (aa)
}  

## function to determine if a change between two codons is synonym
is_synonymous<-function(codon1, codon2) {
  return (dna_to_aa(codon1) == dna_to_aa(codon2))
}

## function to estimate the number of synonymous sites in a codon:
synpos<-function(codon) {
  pos<-s2c(codon)
  for (i in 1:length(pos)) {
    base = pos[i]
    BASES = c('A', 'G', 'T', 'C')
    other_bases = BASES[BASES!=base]
    syn=0
    for (new_base in other_bases) {
      new_pos=c(pos[0:(i-1)], new_base, pos[-(1:i)])
      new_codon=paste(new_pos, collapse="")
      s<-(is_synonymous(codon, new_codon))
      syn <- syn + length(s[s==TRUE])
    }
    nonsynp <- syn/3
  }
  return (nonsynp)
}

syn1<-synpos("ATGATGCAGAGTCTGTAA")
syn2<-synpos("ATGAGGCACAGTCTGTAA")

nonsyn1<-nchar("ATGATGCAGAGTCTGTAA")-syn1
nonsyn2<-nchar("ATGAGGCACAGTCTGTAA")-syn2


library(phylotools)

old_names<-get.fasta.name("./obp83b_cds_aln.fasta")
## new names must be in the same order than old names...
new_names<-c("Dmel_obp83b", 
             "Dsim_obp83b", 
             "Dere_obp83b", 
             "Dele_obp83b", 
             "Dpse_obp83b")
ref2 <- data.frame(old_names, new_names)
rename.fasta(infile = "./obp83b_cds_aln.fasta", ref_table = ref2, outfile = "./obp83b_cds_renamed_aln.fasta")

cds <- read.alignment("./obp83b_cds_renamed_aln.fasta",format="fasta")

## matrix with synonymous divergences (Ks)
kaks(cds)$ks

## matrix with nonsynonymous divergences (Ka)
kaks(cds)$ka

## example for D. pseudoobscura:
x<-efetch("NC_046679", db="nucleotide", rettype="gb", retmode="text", strand=1, seqstart=16336491, seqstop=16337172)
con <- content(x, as = "textConnection")
readLines(con)

## example for D.pseudobscura
## download the complete gene region
efetch("NC_046679", db="nucleotide", rettype="fasta", retmode="text", strand=1, seqstart=16336491,      seqstop=16337172, outfile="./Dpse_opb83.fasta")

## create an object with the sequence of the complete gene region
seq<-read.dna("./Dpse_opb83.fasta", format="fasta")

## create a new sequence with only intron sequences
## CDS = join(54..128,186..261,317..594)
## 5' region = 1-53
## intron 1 = 129-185
## intron2 = 262-316
seqnc<-seq
seqnc=seqnc[,c(1:53,129:185,262:316)]

## change species label
seqnc<-updateLabel(seqnc, labels(seqnc), "Dpse_obp83b")

## write a fasta file with noncoding sequences of the opb83b gene of ALL SPECIES (with append=TRUE, if your run this script for all species, all sequences will be added to the same file; remember to change the ids in each case)
write.dna(seqnc, format="fasta", file="./obp83b_nc.fasta", append=TRUE)




## read the noncoding sequences
myseqs<-read.dna("./obp83b_nc.fasta", format="fasta")


## visualize the sequences
myseqs

## multiple sequence alignment
myalign <- clustal(myseqs)

## partial visualization of the generated alignment
myalign

## write the alignment to file
write.dna(myalign, file="./obp83b_nc_aln.fasta", format="fasta")

## read the noncoding alignment
nc <- read.dna("./obp83b_nc_aln.fasta", format="fasta")

## noncoding divergences in the opb83b gene
dist.dna(nc)



library(phangorn)

cds <- read.alignment("./obp83b_cds_renamed_aln.fasta",format="fasta")

## NJ tree
syntree <- NJ(kaks(cds)$ka)
plot(syntree)
write.tree(syntree, file="./obp83b_Ks.tree")


## read the tree to a phylo class object
tree<-read.tree(file="./obp83b_Ks.tree")

## plot the unroted tree
plot(tree, main="Obp gene tree (number of nonsynonymous substitutions per site)")

## root the tree using Drosophila pseudoobscura as the outgroup
rooted<-root(tree, "Dpse_obp83b")
plot(rooted)
add.scale.bar(ask=TRUE)

## set calibration
## select the node corresponding to the ancestor of D. melanogater, D. simulans and D. erecta, and     specify a range of 6.1-18.9 MYA
cal <- makeChronosCalib(rooted, interactive = TRUE)

## fit the model (penalized likelihood)
chr<-chronos(rooted, model="clock", calibration = cal)

# plot de calibrated (ultrametric) tree
plot(chr, main="calibrated tree (million years)")
axisPhylo()
tiplabels()
nodelabels()
chr$edge
chr$edge.length