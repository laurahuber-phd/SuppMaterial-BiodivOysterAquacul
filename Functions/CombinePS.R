
CONSENSUSFILTER <- function(TEMPLATE.melt) {
  # this function takes a psmelt sequencing database and trims the sequence assignments based on consensus score
  # can keep anything that's a consensus of 1
  TEMPLATE.melt.con1 <- TEMPLATE.melt %>% filter(., Consensus == "1.0000000") %>% filter(., Abundance > 0) 
  # trim consensus 0.75 and 0.8 to Order
  TEMPLATE.melt.conord <- TEMPLATE.melt %>% filter(., Consensus == "0.7500000" | Consensus == "0.8000000") %>% filter(., Abundance > 0)
  TEMPLATE.melt.conord$Species <- NA
  TEMPLATE.melt.conord$Genus <- NA
  TEMPLATE.melt.conord$Family <- NA
  # trim consensus 0.6 / anything lower than 0.75 to Phylum
  TEMPLATE.melt.conphy <- TEMPLATE.melt %>% filter(., Consensus == "0.6000000" | Consensus == "0.6666667") %>% filter(., Abundance > 0)
  TEMPLATE.melt.conphy$Species <- NA
  TEMPLATE.melt.conphy$Genus <- NA
  TEMPLATE.melt.conphy$Family <- NA
  TEMPLATE.melt.conphy$Order <- NA
  TEMPLATE.melt.conphy$Class <- NA
  
  # put them back together 
  TEMPLATE.melt.trim <- rbind(TEMPLATE.melt.con1,TEMPLATE.melt.conord,TEMPLATE.melt.conphy)
  TEMPLATE.melt.trim
}

# this function has not been double checked 
CONSENSUSFILTERPS <- function(TEMPLATE.ps) {
  # subsets database 
  TEMPLATE.MET <- TEMPLATE.ps %>% 
    prune_samples(sample_sums(.)>=1000, .) %>% 
    subset_taxa(Kingdom == "Metazoa" | Kingdom == "Eukaryota")
  
  if (sum(str_detect(tax_table(TEMPLATE.MET)[,8],"1.0000000")) > 0 ){
  # make database for 100 percent consensus assignments 
  TEMPLATE.con1 <-  TEMPLATE.MET %>% 
    subset_taxa(., Consensus == "1.0000000") %>% 
    prune_taxa(taxa_sums(.) > 0, .) %>% 
    prune_samples(sample_sums(.)>0,.)
  }
  # make database for lower consensus assignment percentages 
  if (sum(str_detect(tax_table(TEMPLATE.MET)[,8],"0.8000000")) > 0 ){
  TEMPLATE.confam <- TEMPLATE.MET %>% 
    subset_taxa(., Consensus == "0.8000000") %>% 
    prune_taxa(taxa_sums(.) > 0, .) %>% 
    prune_samples(sample_sums(.)>0,.)
  
  tax_table(TEMPLATE.confam)[,"Species"] = NA
  tax_table(TEMPLATE.confam)[,"Genus"] = NA
  }
  # trim consensus at 0.75
  if (sum(str_detect(tax_table(TEMPLATE.MET)[,8],"0.7500000")) > 0 ){
  TEMPLATE.conord <- TEMPLATE.MET %>%
    subset_taxa(., Consensus == "0.7500000" ) %>% 
    prune_taxa(taxa_sums(.) > 0, .) %>% 
    prune_samples(sample_sums(.)>0,.)
  
  tax_table(TEMPLATE.conord)[,"Species"] = NA
  tax_table(TEMPLATE.conord)[,"Genus"] = NA
  tax_table(TEMPLATE.conord)[,"Family"] = NA
  }
  # trim consensus 0.6 / anything lower than 0.75 to Phylum
  if (sum(str_detect(tax_table(TEMPLATE.MET)[,8],"0.6000000") | str_detect(tax_table(TEMPLATE.MET)[,8],"0.6666667")) > 0 ){
  TEMPLATE.conphy <- TEMPLATE.MET %>% 
    subset_taxa(., Consensus == "0.6000000" | Consensus == "0.6666667") %>% 
    prune_taxa(taxa_sums(.) > 0, .) %>% 
    prune_samples(sample_sums(.)>0,.)
  
  tax_table(TEMPLATE.conphy)[,"Species"] = NA
  tax_table(TEMPLATE.conphy)[,"Genus"] = NA
  tax_table(TEMPLATE.conphy)[,"Family"] = NA
  tax_table(TEMPLATE.conphy)[,"Class"] = NA
  tax_table(TEMPLATE.conphy)[,"Order"] = NA
  }
  # put them back together 
  if ((exists("TEMPLATE.con1") & 
       exists("TEMPLATE.conord") & 
       exists("TEMPLATE.confam") & 
       exists("TEMPLATE.conphy")) == TRUE) {
    TEMPLATE.trim <- merge_phyloseq(TEMPLATE.con1,TEMPLATE.conord,TEMPLATE.confam,TEMPLATE.conphy)}
  else if ((exists("TEMPLATE.con1") & 
            exists("TEMPLATE.conord") & 
            #exists("TEMPLATE.confam") & 
            exists("TEMPLATE.conphy")) == TRUE & 
           exists("TEMPLATE.confam") == FALSE) {
    TEMPLATE.trim <- merge_phyloseq(TEMPLATE.con1,TEMPLATE.conord,TEMPLATE.conphy)
  }
  else if ((exists("TEMPLATE.con1") & 
            #exists("TEMPLATE.conord") & 
            #exists("TEMPLATE.confam") & 
            exists("TEMPLATE.conphy")) == TRUE & 
           (exists("TEMPLATE.conord") & 
            exists("TEMPLATE.confam")) == FALSE) {
    TEMPLATE.trim <- merge_phyloseq(TEMPLATE.con1,TEMPLATE.conphy)
  }
  else if ((exists("TEMPLATE.con1") & 
            #exists("TEMPLATE.conord") & 
            exists("TEMPLATE.confam") & 
            exists("TEMPLATE.conphy")) == TRUE & 
           exists("TEMPLATE.conord") == FALSE) {
             TEMPLATE.trim <- merge_phyloseq(TEMPLATE.con1,TEMPLATE.confam,TEMPLATE.conphy)
           }

  TEMPLATE.trim
}

PSSPGENGLOM <- function(TEMPLATE.subsetMET){
# This function will collapse a dataset by genera
# it does this by splitting the dataset into three groups 
# group one is the dataset with all taxonomic assignments from Kingdom to Family
  # these are not collapsed
# group two is the dataset that has all assignments that go to Genera - these are then collapsed
# group three is the dataset that has all Species level assignments - these are then collapsed
TEMPLATE.subsetMET.na.gen <- TEMPLATE.subsetMET %>% 
  subset_taxa(., is.na(Genus))  
TEMPLATE.subsetMET.genglom <- TEMPLATE.subsetMET %>% 
  subset_taxa(., is.na(Species) & !is.na(Genus)) %>% 
  tax_glom(., taxrank = "Genus") 
TEMPLATE.subsetMET.spglom <- TEMPLATE.subsetMET %>%
  subset_taxa(., !is.na(Species)) %>% 
  tax_glom(., taxrank = "Species") 

merge_phyloseq(TEMPLATE.subsetMET.na.gen,TEMPLATE.subsetMET.genglom,TEMPLATE.subsetMET.spglom)
}

TRANSFORM <- function(TEMPLATE.GLOM,TRANSFORMATION) {
  if (TRANSFORMATION == "Log") {
    TEMPLATE.GLOM %>% 
      #subset_samples(., Location != "CONTROL") %>% 
      #subset_taxa(., !is.na(Phylum) ) %>% 
      prune_samples(sample_sums(.) > 0,.) %>% 
      transform_sample_counts(., function(x) log(1 + x)) %>% 
      prune_samples(sample_sums(.) > 0,.)
    }
  else if (TRANSFORMATION == "Rank") {
    TEMPLATE.GLOM %>% 
      #subset_samples(., Location != "CONTROL") %>% 
      #subset_taxa(., !is.na(Phylum) ) %>% 
      prune_samples(sample_sums(.) > 0,.) %>% 
      transform_sample_counts(.,rank,ties.method="min") %>% 
      prune_samples(sample_sums(.) > 0,.)
  }
  else if (TRANSFORMATION == "Presence") {
    TEMPLATE.GLOM %>% 
      #subset_samples(., Location != "CONTROL") %>% 
      #subset_taxa(., !is.na(Phylum) ) %>% 
      prune_samples(sample_sums(.) > 0,.) %>% 
      transform_sample_counts(., function(abund) 1*(abund > 0)) %>% 
      prune_samples(sample_sums(.) > 0,.)
  }
  else if (TRANSFORMATION == "Relative") {
    TEMPLATE.GLOM %>% 
      #subset_samples(., Location != "CONTROL") %>% 
      #subset_taxa(., !is.na(Phylum) ) %>% 
      prune_samples(sample_sums(.) > 0,.) %>% 
      transform_sample_counts(., function(x) x/sum(x)) %>% 
                                prune_samples(sample_sums(.) > 0,.) 
  }
  else { print("Available transformations include 'Log', 'Rank', 'Relative' (aka Relative abundance), and 'Presence' (aka Presence/Absence)") }
  
}

COMBINEMELT <- function(TEMPLATE.transformed,IDENTIFIER) {
  IDENTIFIER.melt.df <- psmelt(TEMPLATE.transformed)
  IDENTIFIER.melt.df$OTU <- sub("^", IDENTIFIER, IDENTIFIER.melt.df$OTU ) # add run identifier to otu names so otus aren't accidentally combined  
  IDENTIFIER.melt.df$SampleName <- gsub(" ", "_",IDENTIFIER.melt.df$SampleName) # replace spaces 
  #IDENTIFIER.melt.df$Primers %>% unique() # check unique values of primers
  # change Primer names as needed (especially important for COI/12S)
  IDENTIFIER.melt.df$Sample=ifelse(IDENTIFIER.melt.df$Primers == "18S",paste0(IDENTIFIER.melt.df$SampleName,'_18S'),
                                  ifelse(IDENTIFIER.melt.df$Primers == "nCOI",paste0(IDENTIFIER.melt.df$SampleName,'_COI'),
                                         ifelse(IDENTIFIER.melt.df$Primers == "COI",paste0(IDENTIFIER.melt.df$SampleName,'_COI'),
                                         ifelse(IDENTIFIER.melt.df$Primers == "MiFish",paste0(IDENTIFIER.melt.df$SampleName,'_12S'),
                                                ifelse(IDENTIFIER.melt.df$Primers == "12S",paste0(IDENTIFIER.melt.df$SampleName,'_12S'),0)
                                                )
                                         )
                                         )
                                  ) 
  IDENTIFIER.melt.df
}


COMBINEpostMELT <- function(IDENTIFIER.melt.df,IDENTIFIER) {
  IDENTIFIER.melt.df$OTU <- sub("^", IDENTIFIER, IDENTIFIER.melt.df$OTU ) # add run identifier to otu names so otus aren't accidentally combined  
  IDENTIFIER.melt.df$SampleName <- gsub(" ", "_",IDENTIFIER.melt.df$SampleName) # replace spaces 
  #IDENTIFIER.melt.df$Primers %>% unique() # check unique values of primers
  # change Primer names as needed (especially important for COI/12S)
  IDENTIFIER.melt.df$Sample=ifelse(IDENTIFIER.melt.df$Primers == "18S",paste0(IDENTIFIER.melt.df$SampleName,'_18S'),
                                   ifelse(IDENTIFIER.melt.df$Primers == "nCOI",paste0(IDENTIFIER.melt.df$SampleName,'_COI'),
                                          ifelse(IDENTIFIER.melt.df$Primers == "COI",paste0(IDENTIFIER.melt.df$SampleName,'_COI'),
                                                 ifelse(IDENTIFIER.melt.df$Primers == "MiFish",paste0(IDENTIFIER.melt.df$SampleName,'_12S'),
                                                        ifelse(IDENTIFIER.melt.df$Primers == "12S",paste0(IDENTIFIER.melt.df$SampleName,'_12S'),0)
                                                 )
                                          )
                                   )
  ) 
  IDENTIFIER.melt.df
}


COMBINEMELTsample <- function(TEMPLATE.transformed,IDENTIFIER) {
  IDENTIFIER.melt.df <- psmelt(TEMPLATE.transformed)
  IDENTIFIER.melt.df$OTU <- sub("^", IDENTIFIER, IDENTIFIER.melt.df$OTU ) # add run identifier to otu names so otus aren't accidentally combined  
  IDENTIFIER.melt.df$SampleName <- gsub(" ", "_",IDENTIFIER.melt.df$SampleName) # replace spaces 
  #IDENTIFIER.melt.df$Primers %>% unique() # check unique values of primers
  # change Primer names as needed (especially important for COI/12S)
  IDENTIFIER.melt.df$OTU_new=ifelse(IDENTIFIER.melt.df$Primers == "18S",paste0(IDENTIFIER.melt.df$OTU,'_18S'),
                                    ifelse(IDENTIFIER.melt.df$Primers == "nCOI",paste0(IDENTIFIER.melt.df$OTU,'_COI'),
                                           ifelse(IDENTIFIER.melt.df$Primers == "COI",paste0(IDENTIFIER.melt.df$OTU,'_COI'),
                                                  ifelse(IDENTIFIER.melt.df$Primers == "MiFish",paste0(IDENTIFIER.melt.df$OTU,'_12S'),
                                                         ifelse(IDENTIFIER.melt.df$Primers == "12S",paste0(IDENTIFIER.melt.df$OTU,'_12S'),0)
                                                  )
                                          )
                                   )
  ) 
  IDENTIFIER.melt.df
}


COMBINEpostMELTsample <- function(IDENTIFIER.melt.df,IDENTIFIER) {
  IDENTIFIER.melt.df$OTU <- sub("^", IDENTIFIER, IDENTIFIER.melt.df$OTU ) # add run identifier to otu names so otus aren't accidentally combined  
  IDENTIFIER.melt.df$SampleName <- gsub(" ", "_",IDENTIFIER.melt.df$SampleName) # replace spaces 
  #IDENTIFIER.melt.df$Primers %>% unique() # check unique values of primers
  # change Primer names as needed (especially important for COI/12S)
  IDENTIFIER.melt.df$OTU_new=ifelse(IDENTIFIER.melt.df$Primers == "18S",paste0(IDENTIFIER.melt.df$OTU,'_18S'),
                                   ifelse(IDENTIFIER.melt.df$Primers == "nCOI",paste0(IDENTIFIER.melt.df$OTU,'_COI'),
                                          ifelse(IDENTIFIER.melt.df$Primers == "COI",paste0(IDENTIFIER.melt.df$OTU,'_COI'),
                                                 ifelse(IDENTIFIER.melt.df$Primers == "MiFish",paste0(IDENTIFIER.melt.df$OTU,'_12S'),
                                                        ifelse(IDENTIFIER.melt.df$Primers == "12S",paste0(IDENTIFIER.melt.df$OTU,'_12S'),0)
                                                 )
                                          )
                                   )
  ) 
  IDENTIFIER.melt.df
}





COMBINEMELTshort <- function(TEMPLATE.transformed,IDENTIFIER) {
  IDENTIFIER.melt.df <- psmelt(TEMPLATE.transformed)
  IDENTIFIER.melt.df$OTU <- sub("^", IDENTIFIER, IDENTIFIER.melt.df$OTU ) # add run identifier to otu names so otus aren't accidentally combined  
  IDENTIFIER.melt.df$SampleName <- gsub(" ", "_",IDENTIFIER.melt.df$SampleName) # replace spaces 
  #IDENTIFIER.melt.df$Primers %>% unique() # check unique values of primers
  # change Primer names as needed (especially important for COI/12S)
  #IDENTIFIER.melt.df$Sample=ifelse(IDENTIFIER.melt.df$Primers == "18S",paste0(IDENTIFIER.melt.df$SampleName,'_18S'),
  #                                ifelse(IDENTIFIER.melt.df$Primers == "nCOI",paste0(IDENTIFIER.melt.df$SampleName,'_COI'),
  #                                       ifelse(IDENTIFIER.melt.df$Primers == "COI",paste0(IDENTIFIER.melt.df$SampleName,'_COI'),
  #                                       ifelse(IDENTIFIER.melt.df$Primers == "MiFish",paste0(IDENTIFIER.melt.df$SampleName,'_12S'),
  #                                              ifelse(IDENTIFIER.melt.df$Primers == "12S",paste0(IDENTIFIER.melt.df$SampleName,'_12S'),0)
  #                                              )
  #                                       )
  #                                       )
  #                                ) 
  IDENTIFIER.melt.df
}

COMBINErenameSAMPLE <- function(IDENTIFIER.melt.df) {
  #IDENTIFIER.melt.df <- psmelt(TEMPLATE.transformed)
  #IDENTIFIER.melt.df$OTU <- sub("^", IDENTIFIER, IDENTIFIER.melt.df$OTU ) # add run identifier to otu names so otus aren't accidentally combined  
  #IDENTIFIER.melt.df$SampleName <- gsub(" ", "_",IDENTIFIER.melt.df$SampleName) # replace spaces 
  IDENTIFIER.melt.df$OTU_new=ifelse(IDENTIFIER.melt.df$Primers == "18S",paste0(IDENTIFIER.melt.df$OTU,'_18S'),
                                    ifelse(IDENTIFIER.melt.df$Primers == "nCOI",paste0(IDENTIFIER.melt.df$OTU,'_COI'),
                                           ifelse(IDENTIFIER.melt.df$Primers == "COI",paste0(IDENTIFIER.melt.df$OTU,'_COI'),
                                                  ifelse(IDENTIFIER.melt.df$Primers == "MiFish",paste0(IDENTIFIER.melt.df$OTU,'_12S'),
                                                         ifelse(IDENTIFIER.melt.df$Primers == "12S",paste0(IDENTIFIER.melt.df$OTU,'_12S'),0)
                                                  )
                                           )
                                    )
  ) 
  IDENTIFIER.melt.df
}


