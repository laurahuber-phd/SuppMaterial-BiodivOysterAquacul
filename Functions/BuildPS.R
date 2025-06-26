BUILDTAXTAB <- function (TEMPLATE.tax) {
  TMP.tax <- TEMPLATE.tax %>% 
    mutate(SeqID = (word(Feature.ID, 2, sep=";") ), .before = 1)%>% 
    mutate(header = (word(Feature.ID, 3, sep=";") ), .before = 3) %>% 
    mutate_at(c('header'), ~na_if(., ''))
  #View(TMP.tax)
  TMP.tax <- TMP.tax %>% mutate(Duplicated = duplicated(TMP.tax$SeqID), .before = 4) %>% 
    filter(., Duplicated != "TRUE")
  TMP.tax$Feature.ID <- NULL
  TMP.tax$header <- NULL
  TMP.tax$Duplicated <- NULL
  TMP.tax <- separate(TMP.tax,Taxon,into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),sep = "[;]")
  TMP.tax$Kingdom <- gsub("k__","",as.character(TMP.tax$Kingdom))
  TMP.tax$Phylum <- gsub("p__","",as.character(TMP.tax$Phylum))
  TMP.tax$Class <- gsub("c__","",as.character(TMP.tax$Class))
  TMP.tax$Order <- gsub("o__","",as.character(TMP.tax$Order))
  TMP.tax$Family <- gsub("f__","",as.character(TMP.tax$Family))
  TMP.tax$Genus <- gsub("g__","",as.character(TMP.tax$Genus))
  TMP.tax$Species <- gsub("s__","",as.character(TMP.tax$Species))
  
  row.names(TMP.tax) <- TMP.tax$SeqID
  TMP.tax$SeqID <- NULL
  
  TMP.tax
}

BELOWTAXTAB <- function (TEMPLATEfastaBelow) { 
  TEMPLATEfastaBelownames <- names(TEMPLATEfastaBelow)
  TEMPLATEfastaBelowseqs <- paste(TEMPLATEfastaBelow)
  TEMPLATEfastaBelow.df <- data.frame(TEMPLATEfastaBelownames,TEMPLATEfastaBelowseqs)
  
  new.TEMPLATEfastaBelow <- TEMPLATEfastaBelow.df %>% mutate(barcode = (word(TEMPLATEfastaBelownames, 1,sep = ";")), .before = 1) %>%
    mutate(SeqID = (word(TEMPLATEfastaBelownames, 2, sep = ";")), .before = 2) %>%
    mutate(size = (word(TEMPLATEfastaBelownames, 3, sep = ";")), .before = 3)
  
  TEMPLATEfastaBelow.tax <- new.TEMPLATEfastaBelow %>% mutate(Barcode = gsub("barcodelabel=","",barcode)) %>%
    mutate(UniqueID = gsub("otu=","",SeqID)) %>%
    mutate(Abundance = as.numeric(gsub("size=","",size)))
  TEMPLATEfastaBelow.tax$barcode <- NULL
  TEMPLATEfastaBelow.tax$Barcode <- NULL
  TEMPLATEfastaBelow.tax$Abundance <- NULL
  TEMPLATEfastaBelow.tax$SeqID <- NULL
  TEMPLATEfastaBelow.tax$size <- NULL
  TEMPLATEfastaBelow.tax$TEMPLATEfastaBelownames <- NULL
  TEMPLATEfastaBelow.tax$TEMPLATEfastaBelowseqs <- NULL
  
  TEMPLATEfastaBelow.tax$Kingdom <- "Unassigned"
  TEMPLATEfastaBelow.tax$Phylum <- NA
  TEMPLATEfastaBelow.tax$Class <- NA
  TEMPLATEfastaBelow.tax$Order <- NA
  TEMPLATEfastaBelow.tax$Family <- NA
  TEMPLATEfastaBelow.tax$Genus <- NA
  TEMPLATEfastaBelow.tax$Species <- NA
  TEMPLATEfastaBelow.tax$Consensus <- 0
  
  class(TEMPLATEfastaBelow.tax$Phylum) <- "character"
  class(TEMPLATEfastaBelow.tax$Class) <- "character"
  class(TEMPLATEfastaBelow.tax$Order) <- "character"
  class(TEMPLATEfastaBelow.tax$Family) <- "character"
  class(TEMPLATEfastaBelow.tax$Genus) <- "character"
  class(TEMPLATEfastaBelow.tax$Species) <- "character"
  
  row.names(TEMPLATEfastaBelow.tax) <- TEMPLATEfastaBelow.tax$UniqueID
  TEMPLATEfastaBelow.tax$UniqueID <- NULL
  TEMPLATEfastaBelow.tax
}
