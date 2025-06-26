
################  Load libraries needed  ################ 
library(phyloseq)
library(readxl)
library(dplyr)
library(tidyverse)
library(vegan)
library(stringr)
library(RColorBrewer)
library(tinytex)
library(knitr)
library(kableExtra)
library(patchwork)
library(ggrepel)
library(cowplot)
library(Biostrings)
library(data.table)
library(taxize)
library(traits)
library(formatR)

source("Functions/BuildPS.R")
source("Functions/CombinePS.R")


################# Import files and make phyloseq object #################    

PS.seq <- read.csv("18S_12S_COI_allclusters_seqtab.csv")
row.names(PS.seq) <- PS.seq$UniqueID
PS.seq$UniqueID <- NULL
PS.seq$X <- NULL

nrow(PS.seq) #check the number of rows in the seq table  
colSums(PS.seq) %>% sort() #Check sample sums  


################ Import taxonomy assignments  ################    
PS18S.tax <- read.delim("18S.qiimeblast.taxonomy.tsv") %>% BUILDTAXTAB(.)
PSCOI.tax <- read.delim("COI.qiimeblast.taxonomy.tsv") %>% BUILDTAXTAB(.)
PS12S.tax <- read.delim("12S.qiimeblast.taxonomy.tsv") %>% BUILDTAXTAB(.)

#Generate taxonomy table for unassigned sequences  
PSfastaBelow <- readDNAStringSet("clustersbelowmin10threshold.fasta")
PSfastaBelow.tax <- BELOWTAXTAB(PSfastaBelow)

write.csv(PSfastaBelow.tax,"clustersbelowthreshold_seqtab.csv")

PSfastaBelow.tax <- read.csv("clustersbelowthreshold_seqtab.csv",stringsAsFactors = FALSE,colClasses = c("character","character","character","character","character","character","character","character","integer"))
row.names(PSfastaBelow.tax) <- PSfastaBelow.tax$X
PSfastaBelow.tax$X <- NULL
 
PS.tax <- bind_rows(PSCOI.tax, PS18S.tax, PS12S.tax, PSfastaBelow.tax)
#PS.tax[is.na(PS.tax)] <- 0

# check number of rows in each component  
nrow(PSCOI.tax)
nrow(PS18S.tax)
nrow(PS12S.tax)
nrow(PSfastaBelow.tax)
nrow(PS.tax)

################ Import metadata   ################    
PS.metadata <- read_excel("eDNA_metadata.xlsx", na = ".")
PS.metadata$Barcode <- PS.metadata$Barcode %>% str_to_lower( locale = "en")
rownames(PS.metadata) <- PS.metadata$Barcode
PS.metadata <- as.data.frame(PS.metadata)
PS.metadata$FileName <- NULL

#remove barcodes that aren't in both the seq table and metadata  
keep <- list(colnames(PS.seq))
subsetPS.metadata <- PS.metadata[PS.metadata$Barcode %in% keep[[1]],] 

################ Make the phyloseq object  ################    
PS.seq <- as.matrix(PS.seq) %>% otu_table(., taxa_are_rows = TRUE)
PS.TAX <- as.matrix(PS.tax) %>% tax_table(.)
PS.META <- sample_data(subsetPS.metadata) 
PS.ps <- phyloseq(PS.seq, PS.TAX, PS.META)
print(PS.ps)

saveRDS(PS.ps,"PS.ps.rds")

