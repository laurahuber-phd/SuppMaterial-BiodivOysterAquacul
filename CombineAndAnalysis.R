
################  Load libraries  ################    
library(phyloseq)
library(webshot2)
library(ggpubr)
library(readxl)
library(dplyr)
library(tidyverse)
library(vegan)
library(stringr)
library(RColorBrewer)
library(kableExtra)
library(tinytex)
library(ggrepel)
library(patchwork)
library(cowplot)
library(traits)
library(taxize)
library(lme4)
library(nlme)
library(ape)
library(ggpattern)
library(pheatmap)
library(gghighlight)
library(pairwiseAdonis)
library(VennDiagram)
library(flextable)
library(formatR)
library(knitr)
library(indicspecies)

# set seed 
set.seed(100)


################  Load data  ################    
AF.ps <- readRDS("R_objects/AF.spgen.glom.bla.ps.RDS")
AG.ps <- readRDS("R_objects/AG.spgen.glom.bla.ps.RDS")
AH.ps <- readRDS("R_objects/AH.spgen.glom.bla.ps.RDS")
AP.ps <- readRDS("R_objects/AP.spgen.glom.bla.ps.RDS")
AWR.ps <- readRDS("R_objects/AWR.spgen.glom.bla.ps.RDS")
AWF.ps <- readRDS("R_objects/AWF.spgen.glom.bla.ps.RDS")

################   Combine all sites  ################  
sum(length(sample_names(AG.ps)),
    length(sample_names(AH.ps)),
    length(sample_names(AP.ps)),
    length(sample_names(AF.ps)),
    length(sample_names(AWR.ps)),
    length(sample_names(AWF.ps))
)
barcodelist<- sprintf("barcode%s",seq(01:239))

sample_names(AG.ps) <- barcodelist[01:44] #length(sample_names(AG.ps))
sample_names(AH.ps) <- barcodelist[45:87] #43 samples
sample_names(AP.ps) <- barcodelist[88:130] #43 samples
sample_names(AF.ps) <- barcodelist[131:173] #43 samples
sample_names(AWR.ps) <- barcodelist[174:204] #31 samples
sample_names(AWF.ps) <- barcodelist[205:239] #35 samples

taxa_names(AG.ps) <- paste("AG_",taxa_names(AG.ps),sep="")
taxa_names(AH.ps) <- paste("AH_",taxa_names(AH.ps),sep="")
taxa_names(AP.ps) <- paste("AP_",taxa_names(AP.ps),sep="")
taxa_names(AF.ps) <- paste("AF_",taxa_names(AF.ps),sep="")
taxa_names(AWR.ps) <- paste("AWR_",taxa_names(AWR.ps),sep="")
taxa_names(AWF.ps) <- paste("AWF_",taxa_names(AWF.ps),sep="")

sample_data(AG.ps)$Run <- "AG"
sample_data(AH.ps)$Run <- "AH"
sample_data(AP.ps)$Run <- "AP"
sample_data(AF.ps)$Run <- "AF"
sample_data(AWR.ps)$Run <- "AWR"
sample_data(AWF.ps)$Run <- "AWF"

NCeDNA.ps <- merge_phyloseq(AG.ps,
                            AH.ps,
                            AP.ps,
                            AF.ps,
                            AWR.ps,
                            AWF.ps)


################  Transformation  ################    
NCeDNA.ps <- NCeDNA.ps %>% 
  subset_taxa(., is.na(Genus) | Genus != "Sus" & Genus != "Homo" & Genus != "Sigmodon" & Genus != "Bos" & Genus != "Crassostrea" ) %>% subset_taxa(Order != "Carnivora")

sample_data(NCeDNA.ps)$Season <- factor(sample_data(NCeDNA.ps)$Season,levels=c("Spring","Summer","Fall"))
sample_data(NCeDNA.ps)$SeasonNew <- factor(sample_data(NCeDNA.ps)$SeasonNew,levels=c("Spring","Summer","Fall"))
sample_data(NCeDNA.ps)$Location <- factor(sample_data(NCeDNA.ps)$Location,levels=c("A1","A2","A3","R1","R2","R3","CONTROL"))
sample_data(NCeDNA.ps)$Month <- factor(sample_data(NCeDNA.ps)$Month,levels=c("Apr","May","Jun","Jul","Aug","Sep","Oct","Nov"))

################  Convert to presence/absence  ################   
NCeDNA.famglom.pa.ps <- NCeDNA.ps %>% 
  subset_samples(!str_detect(SampleName,"Control") & !str_detect(Location,"CONTROL") ) %>% 
  prune_samples(sample_sums(.)>1000,.) %>% 
  prune_taxa(taxa_sums(.)>0,.) %>% 
  tax_glom("Family",NArm=FALSE)
otu_table(NCeDNA.famglom.pa.ps) <- otu_table(decostand(otu_table(NCeDNA.famglom.pa.ps),method="pa"),taxa_are_rows = TRUE)

fam.dist.df <- NCeDNA.famglom.pa.ps %>% psmelt() %>% 
  filter(., !is.na(Family)) %>% 
  mutate(across(where(is.factor), str_remove_all, pattern = fixed(" ")))

ReefFam <- fam.dist.df %>% filter(., !is.na(Abundance) & Abundance >0 & SiteType=="Reef") %>% .$Family %>% unique()
ReefSpringFam <- fam.dist.df %>% filter(., !is.na(Abundance) & Abundance >0 & SiteType=="Reef" & Season=="Spring") %>% .$Family %>% unique()
ReefSummerFam <- fam.dist.df %>% filter(., !is.na(Abundance) & Abundance >0 & SiteType=="Reef" & Season=="Summer") %>% .$Family %>% unique()
ReefFallFam <- fam.dist.df %>% filter(., !is.na(Abundance) & Abundance >0 & SiteType=="Reef" & Season=="Fall") %>% .$Family %>% unique()

AquacultureFam <- fam.dist.df %>% filter(., !is.na(Abundance) & Abundance >0 & SiteType=="Aquaculture") %>% .$Family %>% unique()
AquacultureSpringFam <- fam.dist.df %>% filter(., !is.na(Abundance) & Abundance >0 & SiteType=="Aquaculture" & Season=="Spring") %>% .$Family %>% unique()
AquacultureSummerFam <- fam.dist.df %>% filter(., !is.na(Abundance) & Abundance >0 & SiteType=="Aquaculture" & Season=="Summer") %>% .$Family %>% unique()
AquacultureFallFam <- fam.dist.df %>% filter(., !is.na(Abundance) & Abundance >0 & SiteType=="Aquaculture" & Season=="Fall") %>% .$Family %>% unique()

#estimate richness  
divtable <- estimate_richness(NCeDNA.famglom.pa.ps,measures=c("Shannon","Simpson","InvSimpson","Observed"))
plot_richness(NCeDNA.pa.ps,measures=c("Shannon","Simpson","InvSimpson","Observed"))

shapiro.test(divtable$Observed)
shapiro.test(divtable$InvSimpson)
shapiro.test(divtable$Shannon)

################  Assumption checks  ################   
#check dispersions  
jaccard_seppa_fam <- distance((NCeDNA.famglom.pa.ps), "jaccard",binary=TRUE) #calculate distance matrix
df_seppa_fam <- as(sample_data(NCeDNA.famglom.pa.ps), "data.frame") #make metadata -> needs to be in the same order as the distance matrix
groups_seppa_fam <- df_seppa_fam[["Site"]] #only going to check dispersion of sites first
#groups
mod_seppa_fam <- betadisper(jaccard_seppa_fam, groups_seppa_fam)  #are the communities are homogenously dispersed
anova(mod_seppa_fam)  #perform an anova on the dispersions; if the ANOVA's p-value is not significant it means that group dispersions are homogenous
TukeyHSD(mod_seppa_fam)  #A Tukey's test can be done to see if and which groups differ in relation to their variances

beta.deco.st <- betadisper(jaccard_combpa_fam, df_combpa_fam$SiteType)
beta.deco.sea <- betadisper(jaccard_combpa_fam, df_combpa_fam$Season)
beta.deco.sal <- betadisper(jaccard_combpa_fam, df_combpa_fam$Sal_ppt)
beta.deco.ph <- betadisper(jaccard_combpa_fam, df_combpa_fam$pH)
beta.deco.tid <- betadisper(jaccard_combpa_fam, df_combpa_fam$Tide)
beta.deco.si <- betadisper(jaccard_combpa_fam, df_combpa_fam$Site)
beta.deco.temp <- betadisper(jaccard_combpa_fam, df_combpa_fam$Temp_C)
permutest(beta.deco.st) 
permutest(beta.deco.sea)
permutest(beta.deco.sal) 
permutest(beta.deco.ph) 
permutest(beta.deco.tid) 
permutest(beta.deco.si)
permutest(beta.deco.temp) 

#Mantel Test  
latlongtable <- data.frame(Site=c("Hoop Pole Creek",
                                  "Gloucester",
                                  "Wilmington Reef",
                                  "Duke Aquafarm",
                                  "Carteret Community College",
                                  "Wilmington Farm"),
                           Lat=c(34.7068954,34.7243564,34.141367,34.7017321,34.7226472,34.14047),
                           Long=c(-76.7510587,-76.5493131,-77.857867,-76.7046435,-76.7604147,-77.86231))

Sites.List <- sample_data(NCeDNA.famglom.pa.ps)$SampleName #if using the un-collapsed dataset
latlongtable <- data.frame(Site=Sites.List,
                           Lat=ifelse(str_detect(Sites.List,"AF"),34.7068954,
                                      ifelse(str_detect(Sites.List,"AG"),34.7243564,
                                             ifelse(str_detect(Sites.List,"AF"),34.7017321,
                                                    ifelse(str_detect(Sites.List,"AP"),34.7226472,
                                                           ifelse(str_detect(Sites.List,"AWR"),34.141367,
                                                                  ifelse(str_detect(Sites.List,"AWF"),34.14047,0)))))),
                           Long=ifelse(str_detect(Sites.List,"AF"),-76.7510587,
                                       ifelse(str_detect(Sites.List,"AG"),-76.5493131,
                                              ifelse(str_detect(Sites.List,"AF"),-76.7046435,
                                                     ifelse(str_detect(Sites.List,"AP"),-76.7604147,
                                                            ifelse(str_detect(Sites.List,"AWR"),-77.857867,
                                                                   ifelse(str_detect(Sites.List,"AWF"),-77.86231,0)))))))

latlong.dist <- as.matrix(dist(cbind(latlongtable$Lat,latlongtable$Long)))
spatial.auto.test <- mantel(latlong.dist,(distance((NCeDNA.famglom.pa.ps), "jaccard",binary=TRUE)))

# add locations to metadata  
rownames(divtable) <- gsub("X", "", rownames(divtable))
meta <- bind_cols(Sample=rownames(divtable),as(sample_data(NCeDNA.pa.ps),"data.frame"),divtable) 

meta.m <- meta %>%
  #mutate(Month=word(DateShort,end=1,sep=" ")) %>% 
  mutate(
    Lat=ifelse(str_detect(Location,"R2"),34.7068954,
               ifelse(str_detect(Location,"R1"),34.7243564,
                      ifelse(str_detect(Location,"A1"),34.7017321,
                             ifelse(str_detect(Location,"A2"),34.7226472,
                                    ifelse(str_detect(Location,"R3"),34.141367,
                                           ifelse(str_detect(Location,"A3"),34.14047,0)))))),
    Long=ifelse(str_detect(Location,"R2"),-76.7510587,
                ifelse(str_detect(Location,"R1"),-76.5493131,
                       ifelse(str_detect(Location,"A1"),-76.7046435,
                              ifelse(str_detect(Location,"A2"),-76.7604147,
                                     ifelse(str_detect(Location,"R3"),-77.857867,
                                            ifelse(str_detect(Location,"A3"),-77.86231,0))))))) %>% 
  filter(Sample != "111121.AP")

#check tidal effect  
tidecheck.meta <- meta.m
tidecheck.meta$TideNum <- case_when(tidecheck.meta$Tide =="Rising" ~ 1,
                                    tidecheck.meta$Tide =="High" ~ 2,
                                    tidecheck.meta$Tide =="Falling" ~ 3,
                                    tidecheck.meta$Tide =="Low" ~ 4,
                                    .default=NA
)

cor.test(tidecheck.meta[!is.na(tidecheck.meta$TideNum),]$TideNum,tidecheck.meta[!is.na(tidecheck.meta$TideNum),]$Observed,
         alternative="two.sided",method="pearson") 
vegan::adonis2(tidecheck.meta[!is.na(tidecheck.meta$TideNum),]$Observed~TideNum,tidecheck.meta[!is.na(tidecheck.meta$TideNum),])
summary(aov(Observed ~ TideNum, data = tidecheck.meta))

################  Analysis  ################   
h <- fam.dist.df %>% 
  select(c(Season,Site,SiteType,Family,Abundance)) %>%
  group_by(Season,Site) %>%
  group_split(.keep=T) %>% 
  lapply(.,function(x) subset(x,Abundance>0 & !duplicated(Family)) )

SharedSpringFams <- lapply(h,function(x){ y <- data.frame(Site = x$Site, Family = x$Family, Season=x$Season);y[y$Season=="Spring",]}) %>% 
  Filter(function(x) dim(x)[1]>0,.) %>% 
  do.call("rbind",.) %>% group_by(Family) %>%  
  mutate(shared = n_distinct(Site) == n_distinct(.$Site)) %>%
  filter(.,shared==TRUE) %>% select(Family) %>% unique()
SharedSpringFams <- SharedSpringFams$Family %>% as.character()

SharedSummerFams <- lapply(h,function(x){ y <- data.frame(Site = x$Site, Family = x$Family, Season=x$Season);y[y$Season=="Summer",]}) %>% 
  Filter(function(x) dim(x)[1]>0,.) %>% 
  do.call("rbind",.) %>% group_by(Family) %>%  
  mutate(shared = n_distinct(Site) == n_distinct(.$Site)) %>%
  filter(.,shared==TRUE) %>% select(Family) %>% unique()
SharedSummerFams <- SharedSummerFams$Family %>% as.character()

SharedFallFams <- lapply(h,function(x){ y <- data.frame(Site = x$Site, Family = x$Family, Season=x$Season);y[y$Season=="Fall",]}) %>% 
  Filter(function(x) dim(x)[1]>0,.) %>% 
  do.call("rbind",.) %>% group_by(Family) %>%  
  mutate(shared = n_distinct(Site) == n_distinct(.$Site)) %>%
  filter(.,shared==TRUE) %>% select(Family) %>% unique()
SharedFallFams <- SharedFallFams$Family %>% as.character()

k <- fam.dist.df %>% 
  select(c(Season,Sample,Site,SiteType,Family,Abundance)) %>%
  group_by(SiteType) %>%
  group_split(.keep=T) %>% #makes a list for each intersection (Spring/Farm,Spring/Reef,Summer/Farm,etc)
  lapply(.,function(x) subset(x,Abundance>0)) #makes a list of all >0 Families found in that intersection 

#is there a statistical difference in the species richness of samples across sites, site types, or seasons?  
r.kru.obsit <- kruskal.test(Observed ~ Site, data = meta.m)
r.dun.obsit <- dunn.test::dunn.test(meta.m$Observed,meta.m$Site,method="bonferroni",kw=TRUE,altp=TRUE)

r.kru.obst <- kruskal.test(Observed ~ SiteType, data = meta.m)
r.dun.obst <- dunn.test::dunn.test(meta.m$Observed,meta.m$SiteType,method="bonferroni",kw=TRUE,altp=TRUE)

r.kru.obsea <- kruskal.test(Observed ~ Season, data = meta.m)
r.dun.obsea <- dunn.test::dunn.test(meta.m$Observed,meta.m$Season,method="bonferroni",kw=TRUE,altp=TRUE)

## richness and evenness  
stat.alphadiv <- fam.dist.df %>%
  select(c(Family,Abundance,Sample,Location,Sound,Month,SeasonNew,SiteType)) %>% 
  reframe(Abundance=sum(Abundance),Season=SeasonNew,SiteType=SiteType,Sound=Sound,.by=c(Family,Location,Month)) %>% 
  distinct() %>% mutate(NewSample = paste0(Location,sep="_",Month))
stat.alphadiv.df <- stat.alphadiv %>% select(-c(Abundance,Family)) %>% distinct()

stat.alphadiv.otu <- stat.alphadiv %>% select(c(NewSample,Abundance,Family)) %>% distinct() %>% pivot_wider(names_from=Family,values_from=Abundance) %>% as.data.frame()
rownames(stat.alphadiv.otu) <- stat.alphadiv.otu$NewSample
stat.alphadiv.otu$NewSample <- NULL

stat.shannondiv <- diversity(stat.alphadiv.otu,index="shannon")
stat.invsimpsondiv <- diversity(stat.alphadiv.otu,index="invsimpson")

perm.AD <-how(nperm=999)
setBlocks(perm.AD) <- with(stat.alphadiv.df,Sound)

## Inverse Simpson  
kruskal.test(InvSimpson ~  SiteType, data=meta.subset)
kruskal.test(InvSimpson ~  Season, data=meta.subset)
kruskal.test(InvSimpson ~  Tide, data=meta.subset)
kruskal.test(InvSimpson ~  Temp_C, data=meta.subset)
## Shannon  
kruskal.test(Shannon ~  SiteType, data=meta.subset)
kruskal.test(Shannon ~  Season, data=meta.subset)
kruskal.test(Shannon ~  Tide, data=meta.subset)
kruskal.test(Shannon ~  Temp_C, data=meta.subset)


################ Combine all primers  ################    
NCeDNA.df.zeros <- psmelt(NCeDNA.ps) %>% filter(Abundance > 0)

NCeDNA.bysample.otu <- NCeDNA.df.zeros %>% 
  select(., c(SampleName,OTU,Abundance))  %>% 
  reframe(., Abundance=sum(Abundance), .by=c(SampleName,OTU)) %>% 
  filter(., !str_detect(SampleName,"Control") & !str_detect(SampleName,"EB")) %>% 
  pivot_wider(., id_cols = "OTU", names_from = "SampleName", values_from = "Abundance", values_fill = 0) %>%
  as.data.frame() # %>% as.matrix()
row.names(NCeDNA.bysample.otu) <- NCeDNA.bysample.otu$OTU
NCeDNA.bysample.otu$OTU <- NULL

NCeDNA.bysample.tax <- NCeDNA.df.zeros %>% 
  filter(., !str_detect(SampleName,"Control") & !str_detect(SampleName,"EB")) %>% 
  select(., c(OTU,Kingdom,Phylum,Class,Order,Family,Genus,Species,Consensus)) %>% 
  unique() %>% 
  as.data.frame()
row.names(NCeDNA.bysample.tax) <- NCeDNA.bysample.tax$OTU
NCeDNA.bysample.tax$OTU <- NULL

NCeDNA.bysample.sd <- NCeDNA.df.zeros %>% 
  filter(., !str_detect(SampleName,"Control") & !str_detect(SampleName,"EB") ) %>% 
  select(., -c(OTU,Primers,Sample,Abundance,Kingdom,Phylum,Class,Order,Family,Genus,Species,Consensus)) %>% 
  as.data.frame() %>% distinct()
row.names(NCeDNA.bysample.sd) <- NCeDNA.bysample.sd$SampleName
NCeDNA.bysample.sd$SampleName <-NULL

NCeDNA.bysample.OTU <- as.matrix(NCeDNA.bysample.otu)%>% otu_table(., taxa_are_rows = TRUE)
NCeDNA.bysample.TAX <- as.matrix(NCeDNA.bysample.tax) %>% tax_table(.)
NCeDNA.bysample.SD <- sample_data(NCeDNA.bysample.sd)
NCeDNA.bysample.ps <- phyloseq(NCeDNA.bysample.OTU,NCeDNA.bysample.TAX,NCeDNA.bysample.SD)

df_combpa_fam <- as(sample_data(NCeDNA.bysample.famglom.pa.ps), "data.frame") #make metadata -> needs to be in the same order as the distance matrix
jaccard_combpa_fam <- distance((NCeDNA.bysample.famglom.pa.ps), "jaccard",binary=TRUE) #calculate distance matrix
groups_combpa_fam <- df_combpa_fam[["Site"]] #only going to check dispersion of sites first
mod_combpa_fam <- betadisper(jaccard_combpa_fam, groups_combpa_fam)  #are the communities are homogenously dispersed
anova(mod_combpa_fam)  #perform an anova on the dispersions; if the ANOVA's p-value is not significant it means that group dispersions are homogenous
TukeyHSD(mod_combpa_fam)  #A Tukey's test can be done to see if and which groups differ in relation to their variances
beta.deco.st <- betadisper(jaccard_combpa_fam, df_combpa_fam$SiteType)
beta.deco.sea <- betadisper(jaccard_combpa_fam, df_combpa_fam$Season)
beta.deco.sal <- betadisper(jaccard_combpa_fam, df_combpa_fam$Sal_ppt)
beta.deco.ph <- betadisper(jaccard_combpa_fam, df_combpa_fam$pH)
beta.deco.tid <- betadisper(jaccard_combpa_fam, df_combpa_fam$Tide)
beta.deco.si <- betadisper(jaccard_combpa_fam, df_combpa_fam$Site)
beta.deco.temp <- betadisper(jaccard_combpa_fam, df_combpa_fam$Temp_C)
permutest(beta.deco.st) 
permutest(beta.deco.sea)
permutest(beta.deco.sal) 
permutest(beta.deco.ph) 
permutest(beta.deco.tid) 
permutest(beta.deco.si)
permutest(beta.deco.temp)

########## PERMANOVA   ################     
df_combpa_fam$Temp_C <- scale(df_combpa_fam$Temp_C,center=TRUE,scale=TRUE)
df_combpa_fam$Sal_ppt <- scale(df_combpa_fam$Sal_ppt,center=TRUE,scale=TRUE)

permcomb <-how(nperm=999)
setBlocks(permcomb) <- with(df_combpa_fam,Sound)

NCeDNA.comb.pa.jaccard.adonis <- vegan::adonis2(jaccard_combpa_fam~SiteType*Month,permutations=permcomb,df_combpa_fam)
adonis_jaccard_locmo <- vegan::adonis2(jaccard_combpa_fam~Lo*Month,permutations=permcomb,df_combpa_fam)
pwintx_jaccard_monst <- pairwise.adonis2(jaccard_combpa_fam~Month*SiteType,strata="Sound",df_combpa_fam)
pwintx_jaccard_sitemon <- pairwise.adonis2(jaccard_combpa_fam~Location*Month,strata="Sound",df_combpa_fam)

NCeDNA.comb.pa.raup.adonis <- vegan::adonis2(raup_combpa_fam~SiteType*Month,permutations=permcomb,df_combpa_fam)
pwintx_raup_monst <- pairwise.adonis2(raup_combpa_fam~Month*SiteType,strata="Sound",df_combpa_fam)
adonis_raup_locmo <- vegan::adonis2(raup_combpa_fam~Location*Month,permutations=permcomb,df_combpa_fam)
pwintx_raup_sitemon <- pairwise.adonis2(raup_combpa_fam~Location*Month,strata="Sound",df_combpa_fam)

## RDA  
VP.fam.step <- varpart(t(otu_table(NCeDNA.bysample.famglom.pa.ps)
),~SiteType,~Season,~pH,data=df_combpa_fam)

VP.fam.step.table <- VP.fam.step$part[2][[1]]
VP.fam.step.table

plot(VP.fam.step)

##
NCeDNA.bysample.famglom.pa.ps <- NCeDNA.bysample.ps %>% 
  prune_samples(sample_sums(.)>1000,.) %>% 
  tax_glom(.,"Family",NArm=FALSE) %>% prune_taxa(taxa_sums(.)>0,.)
otu_table(NCeDNA.bysample.famglom.pa.ps) <- otu_table(decostand(otu_table(NCeDNA.bysample.famglom.pa.ps),method="pa"),taxa_are_rows = TRUE)
min_num = 4
min_samples = 3

pruned_min <- NCeDNA.bysample.famglom.pa.ps %>% 
  subset_taxa(!is.na(Family)) %>% 
  prune_taxa(taxa_sums(.)>3,.) %>% 
  filter_taxa(function(x) sum(x > 0) >= min_samples,prune=TRUE)

toptaxa.25 <- cbind(tax_table(pruned_min)[names(sort(rowSums(otu_table(pruned_min)),decreasing = TRUE)[1:25]),2:6], sort(rowSums(otu_table(pruned_min)),decreasing = TRUE)[1:25]) %>% .[,"Family"] %>% as.vector()

toptaxa.aquaculture.all25 <- cbind(tax_table(pruned_min %>% subset_samples(., SiteType == "Aquaculture"))[names(sort(rowSums(otu_table(pruned_min %>% subset_samples(., SiteType == "Aquaculture"))),decreasing = TRUE)[1:25]),2:6], sort(rowSums(otu_table(pruned_min %>% subset_samples(., SiteType == "Aquaculture"))),decreasing = TRUE)[1:25])%>% .[,"Family"] %>% as.vector()

toptaxa.reef.all25 <- cbind(tax_table(pruned_min %>% subset_samples(., SiteType == "Reef"))[names(sort(rowSums(otu_table(pruned_min %>% subset_samples(., SiteType == "Reef"))),decreasing = TRUE)[1:25]),2:6], sort(rowSums(otu_table(pruned_min %>% subset_samples(., SiteType == "Reef"))),decreasing = TRUE)[1:25]) %>% .[,"Family"] %>% as.vector()

######### Indicator species  
indic.sp.season <- fam.dist.df %>% filter(Kingdom=="Metazoa") %>% 
  select(c(Family,Abundance,Sample,Location,Sound,Month,Season,SiteType)) %>% reframe(Abundance=sum(Abundance),Season=Season,SiteType=SiteType,Sound=Sound,.by=c(Family,Location,Season)) %>% distinct() %>% mutate(NewSample = paste0(Location,sep="_",Season))

indic.sp.df <- indic.sp.season %>% select(-c(Abundance,Family)) %>% distinct()

indic.sp.otu <- indic.sp.season %>% select(c(NewSample,Abundance,Family)) %>% distinct() %>% pivot_wider(names_from=Family,values_from=Abundance) %>% as.data.frame()
rownames(indic.sp.otu) <- indic.sp.otu$NewSample
indic.sp.otu$NewSample <- NULL

perm_indic <-how(nperm=999)
setBlocks(perm_indic) <- with(indic.sp.df,Sound)

indval <- multipatt(indic.sp.otu, indic.sp.df$SiteType, 
                    control = perm_indic) 
summary(indval)
sig.indic.spp <- indval$sign[(indval$sign$p.value<0.05)&(!is.na(indval$sign$p.value)),]
rownames(sig.indic.spp) 

save.image("Analysis.RData")