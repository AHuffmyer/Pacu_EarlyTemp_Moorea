# Title: P. acuta GFF adjustments
# Author: J. Ashey, E. Chille
# Date: 09/01/20

# Need to do some pacu gff adjustments so it can run properly in STAR. Here, I'll be adding transcript_id= to 'gene' column because we needs that label to do functional annotation

#Load libraries
library(tidyverse)

#Load  gene gff
pacu.gff <- read.csv(file="Pacu2021/Data/TagSeq/Pocillopora_acuta_HIv1.genes.gff3", header=FALSE, sep="\t", skip=1) 

#rename columns
colnames(pacu.gff) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "gene")
head(pacu.gff)

# Creating transcript id
pacu.gff$transcript_id <- sub(";.*", "", pacu.gff$gene)
pacu.gff$transcript_id <- gsub("ID=", "", pacu.gff$transcript_id) #remove ID= 

# Checking what kinds of ids are in gff
# WHY NO GENE???
unique(pacu.gff$id)
# [1] "CDS"         "exon"        "intron"      "stop_codon"  "start_codon"      

#If id == mRNA, exon, start_codon, stop_codon, CDS, tRNA, add ;transcript_id= <gene line ID without ID= stopping at first ; , else replace with original gene

#make a parent id to add in the gene id column 
pacu.gff$parent_id <- sub(".*Parent=", "", pacu.gff$gene)
pacu.gff$parent_id <- sub(";.*", "", pacu.gff$parent_id)

pacu.gff <- pacu.gff %>% 
  mutate(gene = ifelse(id != "gene", paste0(gene, "transcript_id=", pacu.gff$transcript_id, ";gene_id=", pacu.gff$parent_id),  paste0(gene)))
head(pacu.gff)

# Remove last col
pacu.gff <- pacu.gff[,-10]
pacu.gff <- pacu.gff[,-10]
head(pacu.gff)  

#save file
write.table(pacu.gff, file="Pacu2021/Data/TagSeq/Pacu.GFFannotations.fixed_transcript.gff3", sep="\t", col.names = FALSE, row.names=FALSE, quote=FALSE)
