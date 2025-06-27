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

# set seed 
set.seed(100)

#Load data  
load("Analysis.RData")



################# Figure 1 ################   
#### Richness for each site over time  
highlight.df <- fam.dist.df %>% 
  filter(., !is.na(Abundance) & Abundance > 0) %>% 
  group_by(Sample) %>% 
  mutate(., Date=as.Date(Date),
         NFam=(length(unique(Family)))
  ) %>% ungroup() %>% filter(.,!duplicated(Sample)) 

glou.highlight <- highlight.df %>%
  ggplot(., aes(x = Date, y = NFam,
                color=Location,group=Location)) +
  geom_line(stat = "identity",linewidth=3) +
  scale_x_date(date_breaks="1 month",date_labels = "%b") + 
  gghighlight(Location=="R1",#"Gloucester",
              use_group_by = FALSE,
              use_direct_label = FALSE,
              unhighlighted_params = list(linewidth = 1, colour = alpha("grey",0.4))) + 
  theme_pubclean() + theme(text=element_text(size=16)) + rotate.xlabels + 
  labs(x="2021",y="Number of Families") + noleg

hoop.highlight <- highlight.df %>%
  ggplot(., aes(x = Date, y = NFam,
                color=Location,group=Location)) +
  geom_line(stat = "identity",linewidth=3) +
  scale_x_date(date_breaks="1 month",date_labels = "%b") + 
  gghighlight(Location=="R2",#"Hoop Pole",
              use_group_by = FALSE,
              use_direct_label = FALSE,
              unhighlighted_params = list(linewidth = 1, colour = alpha("grey",0.4))) + 
  theme_pubclean() + theme(text=element_text(size=16)) + rotate.xlabels + 
  labs(x="2021",y="Number of Families") + noleg

dukefarm.highlight <- highlight.df %>%
  ggplot(., aes(x = Date, y = NFam,
                color=Location,group=Location)) +
  geom_line(stat = "identity",linewidth=3) +
  scale_x_date(date_breaks="1 month",date_labels = "%b") + 
  gghighlight(Location=="A1",#"Duke Aquafarm",
              use_group_by = FALSE,
              use_direct_label = FALSE,
              unhighlighted_params = list(linewidth = 1, colour = alpha("grey",0.4))) + 
  theme_pubclean() + theme(text=element_text(size=16)) + rotate.xlabels + 
  labs(x="2021",y="Number of Families") + noleg

collegefarm.highlight <- highlight.df %>%
  ggplot(., aes(x = Date, y = NFam,
                color=Location,group=Location)) +
  geom_line(stat = "identity",linewidth=3) +
  scale_x_date(date_breaks="1 month",date_labels = "%b") + 
  gghighlight(Location=="A2",#"Carteret Community College",
              use_group_by = FALSE,
              use_direct_label = FALSE,
              unhighlighted_params = list(linewidth = 1, colour = alpha("grey",0.4))) + 
  theme_pubclean() + theme(text=element_text(size=16)) + rotate.xlabels + 
  labs(x="2021",y="Number of Families") + noleg

hatcheryfarm.highlight <- highlight.df %>%
  ggplot(., aes(x = Date, y = NFam,
                color=Location,group=Location)) +
  geom_line(stat = "identity",linewidth=3) +
  scale_x_date(date_breaks="1 month",date_labels = "%b") + 
  gghighlight(Location=="A3",#Wilmington Hatchery",
              use_group_by = FALSE,
              use_direct_label = FALSE,
              unhighlighted_params = list(linewidth = 1, colour = alpha("grey",0.4))) + 
  theme_pubclean() + theme(text=element_text(size=16)) + rotate.xlabels + 
  labs(x="2021",y="Number of Families") + noleg

hatcheryreef.highlight <- highlight.df %>%
  ggplot(., aes(x = Date, y = NFam,
                color=Location,group=Location)) +
  geom_line(stat = "identity",linewidth=3) +
  scale_x_date(date_breaks="1 month",date_labels = "%b") + 
  gghighlight(Location=="R3",#"Wilmington Reef",
              use_group_by = FALSE,
              use_direct_label = FALSE,
              unhighlighted_params = list(linewidth = 1, colour = alpha("grey",0.4))) + 
  theme_pubclean() + theme(text=element_text(size=16)) + rotate.xlabels + 
  labs(x="2021",y="Number of Families") + noleg

family.number.bymonth.line.plot  <- wrap_plots(
  dukefarm.highlight, #A1
  collegefarm.highlight, #A2
  hatcheryfarm.highlight, #A3
  glou.highlight, #R1
  hoop.highlight, #R2
  hatcheryreef.highlight #R3
) +
  plot_layout(nrow=2) + 
  plot_annotation(tag_levels=list(c("A1","A2","A3","R1","R2","R3")))

col_grid <- rgb(235, 235, 235, 100, maxColorValue = 255)

#### Shared each season  
rich.sites <- fam.dist.df %>% 
  filter(., !is.na(Abundance) & Abundance > 0) %>%
  select(c(OTU,Sample,Abundance,Run,Location,Site,
           SiteType,Season,DateShort,
           Kingdom,Phylum,Class,Order,Family)) %>% 
  group_by(SiteType,Sample) %>% #for each sample, how many of the families are there? 
  mutate(., BySample_NFam=(length(unique(Family)))) %>% 
  ungroup() %>%
  group_by(Site,Season) %>% #for each season at each site, how many families are there? 
  mutate(., BySite_NFamFOUNDATTypeSeason=(length(unique(Family)))) %>% ungroup() %>% 
  group_by(.,SiteType,Season) %>% #for each season at each site type, how many families are there?
  mutate(.,ByType_NFamFOUNDINSeason = (length(unique(Family))),
         ByType_PercSeasonFamOFTotal = (ifelse(SiteType=="Aquaculture",(length(unique(Family))/length(AquacultureFam)*100),ifelse(SiteType=="Reef",(length(unique(Family))/length(ReefFam)*100),0))) #what percentage of the total number of families at a site type are found in that season? 
  ) %>% ungroup() %>%
  filter(., !duplicated(Sample)) 

Lolli.df <- rich.sites %>% 
  select(c(SiteType,Season,ByType_NFamFOUNDINSeason)) %>% 
  distinct() %>% 
  rbind(.,data.frame(SiteType=c("Shared","Shared","Shared"),
                     Season=c("Spring","Fall","Summer"),
                     ByType_NFamFOUNDINSeason=c(
                       length(ReefSpringFam[ReefSpringFam%in%AquacultureSpringFam]),
                       length(ReefSummerFam[ReefSummerFam%in%AquacultureSummerFam]),
                       length(ReefFallFam[ReefFallFam%in%AquacultureFallFam])) ) )

season.shared.fam.lolli.plot <- Lolli.df %>% filter(.,SiteType!="Shared") %>% 
  ggplot(aes(y=Season,x=ByType_NFamFOUNDINSeason,fill=SiteType)) + 
  geom_linerange( aes(y=Season, xmin=0,
                      xmax=ByType_NFamFOUNDINSeason),
                  position=position_dodge(width=0.3),
                  size=2)  + 
  geom_point(aes(fill=SiteType,color=SiteType),
             position = position_dodge(width=0.3),
             size=8) +
  geom_linerange(data=(Lolli.df %>%
                         filter(.,SiteType=="Shared")),
                 aes(y=Season,
                     xmin=0,xmax=ByType_NFamFOUNDINSeason,
                     color=SiteType),size=12) + 
  theme_pubr() + 
  scale_x_continuous(expand = c(0, NA), 
                     limits = c(0, 150)) + 
  theme(legend.position = "bottom",
        text=element_text(size=18,face="bold"),
        axis.text.y = element_text(angle=270,hjust=0.5),
        panel.grid.major.x = element_line(color=col_grid),
        panel.grid.minor.x = element_line(color=col_grid),
        plot.margin = margin(5,15,5,5)) +
  labs(x="Number of Families",y=NULL,title="Seasonal patterns") +
  scale_fill_manual(name=" ",
                    labels=c("Aquaculture","Reef","Shared"),
                    values=c(Aquaculture="#2C8905",SiteType.Palette.duke[2],Shared="#5b1e5b")) + 
  scale_color_manual(name=" ",
                     labels=c("Aquaculture","Reef","Shared"),
                     values=c(Aquaculture="#2C8905",SiteType.Palette.duke[2],Shared="#5b1e5b")) 

##### Venn Diagram  
v <- venn.diagram(list("Aquaculture" = AquacultureFam, "Reef" = ReefFam),
                  filename=NULL,
                  alpha=c(0.8,0.8),
                  euler.d = TRUE,
                  scaled = FALSE,
                  cat.fontfamily = "Times",
                  cat.pos=c(290, 70),cat.dist=c(-0.1,-0.09),
                  cat.fontface="italic",
                  cat.cex = 1,
                  lwd = c(0.3, 0.3),#vector giving the width of each circle's circumference
                  fontfamily = "Times",
                  fontface="bold",
                  cex = 2,
                  width = 10,
                  height = 10,
                  units = "cm",
                  main="", 
                  main.fontface = "bold",
                  main.fontfamily="Times",
                  main.cex=NULL,
                  imagetype = "png",
                  margin=0.01
)
grid.newpage()
grid.draw(v)



################# Figure 2 Cleveland plot  ################    
pruned_min <- NCeDNA.bysample.famglom.pa.ps %>% 
  subset_taxa(!is.na(Family)) %>% 
  prune_taxa(taxa_sums(.)>3,.) %>% 
  filter_taxa(function(x) sum(x > 0) >= min_samples,prune=TRUE)

pruned_min_psmelt <- psmelt(pruned_min)

pruned_min_sepsample <- NCeDNA.famglom.pa.ps %>% subset_samples(!str_detect(SampleName,"Control") & !str_detect(Location,"CONTROL") 
) %>% 
  subset_taxa(!is.na(Family)) %>% 
  prune_taxa(taxa_sums(.)>3,.) %>%
  filter_taxa(function(x) sum(x > 0) >= min_samples,prune=TRUE) 

pruned_min_sepsample_psmelt <- psmelt(pruned_min_sepsample) 

pruned_min_sepsample_psmelt$Month <- factor(pruned_min_sepsample_psmelt$Month, levels=c("Apr","May","Jun","Jul","Aug","Sep","Oct","Nov"))
pruned_min_sepsample_psmelt$Site <- factor(pruned_min_sepsample_psmelt$Site,  levels=c("Duke Aquafarm",
                                                                                       "Carteret Community College",
                                                                                       "Wilmington Hatchery",
                                                                                       "Gloucester",
                                                                                       "Hoop Pole",
                                                                                       "Wilmington Reef"))
cleveland.df <- pruned_min_sepsample_psmelt %>% 
  filter(.,Family%in%toptaxa.25) %>% 
  filter(Abundance > 0) %>%
  group_by(SiteType,Primers,Month,Phylum,Family) %>% 
  summarize(Primers=unique(Primers),
            Family=unique(Family),
            Phylum=unique(Phylum),
            Kingdom=unique(Kingdom),
            FamAbund=sum(Abundance)) %>% 
  mutate(FarmAbund=case_when(SiteType=="Aquaculture" ~ FamAbund,
                             .default=NA),
         ReefAbund=case_when(SiteType=="Reef" ~ FamAbund,
                             .default=NA)) %>%
  ungroup() %>% 
  group_by(SiteType) %>% 
  arrange(Phylum) %>% 
  mutate(Family=factor(Family,levels=unique(Family)),
         Phylum=factor(Phylum,levels=unique(Phylum)))


cleveland.met.plot <-   ggplot(cleveland.df %>% filter(Kingdom=="Metazoa") ) + 
  geom_point(aes(x=Family,y=FarmAbund,color=SiteType,size=FarmAbund)) + 
  geom_point(aes(x=Family,y=ReefAbund,color=SiteType,size=ReefAbund)) + 
  coord_flip() + 
  facet_grid(Primers~Month,scales="free_y") + 
  scale_color_manual(name="Site Type") + 
  guides(color=guide_legend(override.aes = list(size=6))) +
  scale_size(name="Occurrences \n per month",range=c(2,8)) +
  guides(color=theme(legend.key.size = unit(50,"pt"))) +
  theme_pubclean() + theme(legend.position="right") +
  labs(y="Number of samples a family was present in a month",title="Metazoans") + 
  theme(panel.border = element_rect(fill=NA,color="black",size=0.1),
        plot.margin = margin(c(0,0,0,0),unit="lines"),
        legend.background=element_rect(fill = 'transparent',color='transparent'),
        plot.background=element_rect(fill = 'transparent'),
        rect = element_rect(fill='transparent',color='transparent'),
        panel.grid.major.y =  element_line(color='black',linetype=3), 
        panel.grid.major.x =  element_line(color='grey',linetype=3), 
        axis.line = element_line(colour = "black"),
        text = element_text(size = 18,family = "Times"))

cleveland.cel.plot <-   ggplot(cleveland.df %>% filter(Kingdom=="Chromista") ) + 
  geom_point(aes(x=Family,y=FarmAbund,color=SiteType,size=FarmAbund)) + 
  geom_point(aes(x=Family,y=ReefAbund,color=SiteType,size=ReefAbund)) + 
  coord_flip() + 
  facet_grid(Primers~Month,scales="free_y") + 
  scale_color_manual(name="Site Type") + 
  guides(color=guide_legend(override.aes = list(size=6))) +
  scale_size(name="Occurrences \n per month",range=c(2,8)) +
  guides(color=theme(legend.key.size = unit(50,"pt"))) +
  theme_pubclean() + theme(legend.position="right") +
  labs(y="Number of samples a family was present in a month",title="Eukaryotes") + 
  theme(panel.border = element_rect(fill=NA,color="black",size=0.1),
        plot.margin = margin(c(0,0,0,0),unit="lines"),
        legend.background=element_rect(fill = 'transparent',color='transparent'),
        plot.background=element_rect(fill = 'transparent'),
        rect = element_rect(fill='transparent',color='transparent'),
        panel.grid.major.y =  element_line(color='black',linetype=3), 
        panel.grid.major.x =  element_line(color='grey',linetype=3), 
        axis.line = element_line(colour = "black"),
        text = element_text(size = 18,family = "Times"))

cleveland.cel.12.plot <-   ggplot(cleveland.df %>% filter(Kingdom=="Chromista" & Primers=="12S") ) + 
  geom_point(aes(x=Family,y=FarmAbund,color=SiteType,size=FarmAbund)) + 
  geom_point(aes(x=Family,y=ReefAbund,color=SiteType,size=ReefAbund)) + 
  coord_flip() + 
  facet_grid(Primers~Month,scales="free_y") + 
  scale_color_manual(name="Site Type") + 
  guides(color=guide_legend(override.aes = list(size=6))) +
  scale_size(name="Occurrences \n per month",range=c(2,8)) +
  guides(color=theme(legend.key.size = unit(50,"pt"))) +
  theme_pubclean() + theme(legend.position="right") +
  labs(y="Number of samples a family was present in a month",title="Eukaryotes") + 
  theme(panel.border = element_rect(fill=NA,color="black",size=0.1),
        plot.margin = margin(c(0,0,-3,0),unit="lines"),
        legend.background=element_rect(fill = 'transparent',color='transparent'),
        plot.background=element_rect(fill = 'transparent'),
        rect = element_rect(fill='transparent',color='transparent'),
        panel.grid.major.y =  element_line(color='black',linetype=3), 
        panel.grid.major.x =  element_line(color='grey',linetype=3), 
        axis.line.y = element_line(colour = "black"),
        text = element_text(size = 18,family = "Times"))

cleveland.cel.18coi.plot <-   ggplot(cleveland.df %>% filter(Kingdom=="Chromista" & Primers!="12S") ) + 
  geom_point(aes(x=Family,y=FarmAbund,color=SiteType,size=FarmAbund)) + 
  geom_point(aes(x=Family,y=ReefAbund,color=SiteType,size=ReefAbund)) + 
  coord_flip() + 
  facet_grid(Primers~Month,scales="free_y") + 
  scale_color_manual(name="Site Type") + 
  guides(color=guide_legend(override.aes = list(size=6))) +
  scale_size(name="Occurrences \n per month",range=c(2,8)) +
  guides(color=theme(legend.key.size = unit(50,"pt"))) +
  theme_pubclean() + theme(legend.position="none") +
  labs(y="Number of samples a family was present in a month") + 
  theme(panel.border = element_rect(fill=NA,color="black",size=0.1),
        legend.background=element_rect(fill = 'transparent',color='transparent'),
        plot.background=element_rect(fill = 'transparent'),
        rect = element_rect(fill='transparent',color='transparent'),
        panel.grid.major.y =  element_line(color='black',linetype=3), 
        panel.grid.major.x =  element_line(color='grey',linetype=3), 
        axis.line = element_line(colour = "black"),
        text = element_text(size = 18,family = "Times"))

design5="
AAAAA
AAAAA
AAAAA
BBBBB
CCCCC
CCCCC
CCCCC
"
fig2.barplot <- wrap_plots((cleveland.met.plot + theme(legend.position="none")
),
( cleveland.cel.12.plot + theme(axis.title.y=element_blank(),
                                #legend.position = "none",
                                axis.title.x = element_blank(),
                                axis.text.x = element_blank(),
                                axis.ticks.x = element_blank()
)
),
( cleveland.cel.18coi.plot + theme(
  legend.position = "none",
  strip.text.x = element_blank()
)
),
guides="collect") + 
  plot_layout(design=design5) +
  plot_annotation(tag_levels=list(c('A.','B.',' ')))



################# Figure 3 NMDS ################      
NCeDNA.famglom.pa.jaccard.nmds <- ordinate(
  physeq = NCeDNA.bysample.famglom.pa.ps, 
  k=3,try=500,
  method = "NMDS", 
  distance = "jaccard"
)

NCeDNA.famglom.pa.raup.nmds <- ordinate(
  physeq = NCeDNA.bysample.famglom.pa.ps, 
  k=3,try=500,
  method = "NMDS", 
  #method="cca",
  distance = "raup"
)

NMDS.fam.pa.sitediff <- plot_ordination(NCeDNA.bysample.famglom.pa.ps,NCeDNA.famglom.pa.jaccard.nmds,
                                        type="samples",color="Location",shape="Location") + 
  geom_point(size=3)+
  theme_bw()

NMDS.fam.pa.monthdiff <-plot_ordination(NCeDNA.bysample.famglom.pa.ps,NCeDNA.famglom.pa.jaccard.nmds,
                                        type="samples",color="SiteType",shape="Location") + 
  geom_point(size=3) + facet_wrap(~Month,nrow=1) + stat_ellipse(aes(group=SiteType,fill=SiteType)) +
  theme_bw() + theme(legend.position="bottom")
  
NMDS.raup.fam.pa.sitediff <- plot_ordination(NCeDNA.bysample.famglom.pa.ps,NCeDNA.famglom.pa.raup.nmds,
                                               type="samples",color="Location",shape="Location") + 
  geom_point(size=3)+
  theme_bw()

rev.fig3.nmds <-  wrap_plots(
  (NMDS.fam.pa.sitediff+
     stat_ellipse(aes(group=SiteType,fill=SiteType))),
  (NMDS.raup.fam.pa.sitediff+
     guides(shape=FALSE,color=FALSE)+
     stat_ellipse(aes(group=SiteType,fill=SiteType))),
  (NMDS.fam.pa.monthdiff+
     guides(shape=FALSE)+
     theme(legend.position="right"))
) + 
  plot_annotation(tag_level="A") +
  plot_layout(design="
AABB
CCCC",guides="collect")



################# Figure 4 Heatmap ################    
make_bold_names <- function(mat, rc_fun, rc_names) {
  bold_names <- rc_fun(mat)
  ids <- rc_names %>% match(rc_fun(mat))
  ids %>%
    walk(
      function(i)
        bold_names[i] <<-
        bquote(bold(.(rc_fun(mat)[i]))) %>%
        as.expression()
    )
  bold_names
}

f.3 <- pruned_min %>% 
  prune_taxa(taxa_sums(.)>0,.) %>% 
  subset_taxa(., !is.na(Family) 
  ) %>% 
  psmelt() %>% 
  group_by(., SiteType,Season,Sample) %>% 
  mutate(., RelAbund = Abundance/sum(Abundance)) %>% 
  ungroup()

hm.fam.pa.otu <- f.3 %>% 
  filter(Kingdom=="Metazoa") %>%
  group_by(SiteType,Season,Family) %>% 
  summarise(abundance = sum(Abundance)) %>%
  pivot_wider(.,names_from=c(SiteType,Season),
              names_sep=" ",
              values_from=abundance,
              values_fill=0) %>% 
  ungroup() %>% as.data.frame()

row.names(hm.fam.pa.otu) <- hm.fam.pa.otu$Family
hm.fam.pa.otu$Family <- NULL 

test <- t(decostand(hm.fam.pa.otu,method="pa")%>%filter(rowSums(.)<6)) %>% as.data.frame()
test$Site <- ifelse(str_detect(rownames(test),"Aqua"), "Aquaculture",ifelse(str_detect(rownames(test),"Reef"), "Reef",NA))

A <- test %>% filter(Site=="Aquaculture") %>% select(-Site) %>% .[,colSums(.)>1] %>% colnames()
R <- test %>% filter(Site=="Reef") %>% select(-Site) %>% .[,colSums(.)>1] %>% colnames()
B <- c(A,R) %>% unique()

phy <- f.3 %>% filter(Family%in%B&!duplicated(Family)) %>% 
  arrange(Family) %>% 
  select(c(Family,Phylum)) %>% 
  as.data.frame() %>% 
  mutate(across(where(is.factor), str_remove_all, pattern = fixed(" ")))
rownames(phy) <- phy$Family
phy$Family <- NULL

king <- f.3 %>% filter(!duplicated(Family)) %>% arrange(Family) %>% select(c(Family,Kingdom)) %>% filter(Family%in%rownames(test)) %>% as.data.frame() %>% mutate(across(where(is.factor), str_remove_all, pattern = fixed(" ")))
rownames(king) <- king$Family
king$Family <- NULL

my.breaks <- seq(0, 1, by=0.5)
my.colors <- c("#337FAD","#E0BC6A")

ann_colors <- list(Phylum=NCeDNA.Phylum.Palette.duke[names(NCeDNA.Phylum.Palette.duke)%in%phy$Phylum])
heatmap.bolded.plot <- pheatmap(test[colnames(test)%in%B],
                                width=12,height=5,fontsize=12,
                                annotation_col=phy,
                                annotation_colors=ann_colors,
                                breaks=my.breaks,
                                color=my.colors,
                                legend_breaks = c(0,1),
                                legend = F,
                                treeheight_row = 0,
                                treeheight_col = 0,
                                cutree_rows = 2,
                                labels_col = make_bold_names(test[colnames(test)%in%B],colnames, rownames(sig.indic.spp)[rownames(sig.indic.spp)%in%colnames(test)])
)
grid.draw(heatmap.bolded.plot)



################# Table 2 Reads per primer set for each site  ################    

ah.reads.2assgn.df = data.frame(
  reads18S_assigned_king=sum(sample_sums(AH.ps %>% subset_samples(Primers=="18S") %>%
                                           prune_taxa(taxa_sums(.)>9,.) %>%
                                           subset_taxa(Kingdom!="Unassigned")
  )),
  reads18S_assigned_fam=sum(sample_sums(AH.ps %>% subset_samples(Primers=="18S") %>%
                                          prune_taxa(taxa_sums(.)>9,.) %>%
                                          subset_taxa(!is.na(Family))
  )),
  reads12S_assigned_king=sum(sample_sums(AH.ps %>% subset_samples(Primers=="12S") %>%
                                           prune_taxa(taxa_sums(.)>9,.) %>%
                                           subset_taxa(Kingdom!="Unassigned")
  )),
  reads12S_assigned_fam=sum(sample_sums(AH.ps %>% subset_samples(Primers=="12S") %>%
                                          prune_taxa(taxa_sums(.)>9,.) %>%
                                          subset_taxa(!is.na(Family))
  )),
  readsCOI_assigned_king=sum(sample_sums(AH.ps %>% subset_samples(Primers=="nCOI") %>%
                                           prune_taxa(taxa_sums(.)>9,.) %>%
                                           subset_taxa(Kingdom!="Unassigned")
  )),
  readsCOI_assigned_fam=sum(sample_sums(AH.ps %>% subset_samples(Primers=="nCOI") %>%
                                          prune_taxa(taxa_sums(.)>9,.) %>%
                                          subset_taxa(!is.na(Family))
  ))
)

ag.reads.2assgn.df = data.frame(
  reads18S_assigned_king=sum(sample_sums(AG.ps %>% subset_samples(Primers=="18S") %>%
                                           prune_taxa(taxa_sums(.)>9,.) %>%
                                           subset_taxa(Kingdom!="Unassigned")
  )),
  reads18S_assigned_fam=sum(sample_sums(AG.ps %>% subset_samples(Primers=="18S") %>%
                                          prune_taxa(taxa_sums(.)>9,.) %>%
                                          subset_taxa(!is.na(Family))
  )),
  reads12S_assigned_king=sum(sample_sums(AG.ps %>% subset_samples(Primers=="MiFish") %>%
                                           prune_taxa(taxa_sums(.)>9,.) %>%
                                           subset_taxa(Kingdom!="Unassigned")
  )),
  reads12S_assigned_fam=sum(sample_sums(AG.ps %>% subset_samples(Primers=="MiFish") %>%
                                          prune_taxa(taxa_sums(.)>9,.) %>%
                                          subset_taxa(!is.na(Family))
  )),
  readsCOI_assigned_king=sum(sample_sums(AG.ps %>% subset_samples(Primers=="nCOI") %>%
                                           prune_taxa(taxa_sums(.)>9,.) %>%
                                           subset_taxa(Kingdom!="Unassigned")
  )),
  readsCOI_assigned_fam=sum(sample_sums(AG.ps %>% subset_samples(Primers=="nCOI") %>%
                                          prune_taxa(taxa_sums(.)>9,.) %>%
                                          subset_taxa(!is.na(Family))
  ))
)

awr.reads.2assgn.df = data.frame(
  reads18S_assigned_king=sum(sample_sums(AWR.ps %>% subset_samples(Primers=="18S") %>%
                                           prune_taxa(taxa_sums(.)>9,.) %>%
                                           subset_taxa(Kingdom!="Unassigned")
  )),
  reads18S_assigned_fam=sum(sample_sums(AWR.ps %>% subset_samples(Primers=="18S") %>%
                                          prune_taxa(taxa_sums(.)>9,.) %>%
                                          subset_taxa(!is.na(Family))
  )),
  reads12S_assigned_king=sum(sample_sums(AWR.ps %>% subset_samples(Primers=="12S") %>%
                                           prune_taxa(taxa_sums(.)>9,.) %>%
                                           subset_taxa(Kingdom!="Unassigned")
  )),
  reads12S_assigned_fam=sum(sample_sums(AWR.ps %>% subset_samples(Primers=="12S") %>%
                                          prune_taxa(taxa_sums(.)>9,.) %>%
                                          subset_taxa(!is.na(Family))
  )),
  readsCOI_assigned_king=sum(sample_sums(AWR.ps %>% subset_samples(Primers=="nCOI") %>%
                                           prune_taxa(taxa_sums(.)>9,.) %>%
                                           subset_taxa(Kingdom!="Unassigned")
  )),
  readsCOI_assigned_fam=sum(sample_sums(AWR.ps %>% subset_samples(Primers=="nCOI") %>%
                                          prune_taxa(taxa_sums(.)>9,.) %>%
                                          subset_taxa(!is.na(Family))
  ))
)

awf.reads.2assgn.df = data.frame(
  reads18S_assigned_king=sum(sample_sums(AWF.ps %>% subset_samples(Primers=="18S") %>%
                                           prune_taxa(taxa_sums(.)>9,.) %>%
                                           subset_taxa(Kingdom!="Unassigned")
  )),
  reads18S_assigned_fam=sum(sample_sums(AWF.ps %>% subset_samples(Primers=="18S") %>%
                                          prune_taxa(taxa_sums(.)>9,.) %>%
                                          subset_taxa(!is.na(Family))
  )),
  reads12S_assigned_king=sum(sample_sums(AWF.ps %>% subset_samples(Primers=="12S") %>%
                                           prune_taxa(taxa_sums(.)>9,.) %>%
                                           subset_taxa(Kingdom!="Unassigned")
  )),
  reads12S_assigned_fam=sum(sample_sums(AWF.ps %>% subset_samples(Primers=="12S") %>%
                                          prune_taxa(taxa_sums(.)>9,.) %>%
                                          subset_taxa(!is.na(Family))
  )),
  readsCOI_assigned_king=sum(sample_sums(AWF.ps %>% subset_samples(Primers=="nCOI") %>%
                                           prune_taxa(taxa_sums(.)>9,.) %>%
                                           subset_taxa(Kingdom!="Unassigned")
  )),
  readsCOI_assigned_fam=sum(sample_sums(AWF.ps %>% subset_samples(Primers=="nCOI") %>%
                                          prune_taxa(taxa_sums(.)>9,.) %>%
                                          subset_taxa(!is.na(Family))
  ))
)

ap.reads.2assgn.df = data.frame(
  reads18S_assigned_king=sum(sample_sums(AP.ps %>% subset_samples(Primers=="18S" | is.na(Primers)) %>%
                                           prune_taxa(taxa_sums(.)>9,.) %>%
                                           subset_taxa(Kingdom!="Unassigned")
  )),
  reads18S_assigned_fam=sum(sample_sums(AP.ps %>% subset_samples(Primers=="18S" | is.na(Primers)) %>%
                                          prune_taxa(taxa_sums(.)>9,.) %>%
                                          subset_taxa(!is.na(Family))
  )),
  reads12S_assigned_king=sum(sample_sums(AP.ps %>% subset_samples(Primers=="12S") %>%
                                           prune_taxa(taxa_sums(.)>9,.) %>%
                                           subset_taxa(Kingdom!="Unassigned")
  )),
  reads12S_assigned_fam=sum(sample_sums(AP.ps %>% subset_samples(Primers=="12S") %>%
                                          prune_taxa(taxa_sums(.)>9,.) %>%
                                          subset_taxa(!is.na(Family))
  )),
  readsCOI_assigned_king=sum(sample_sums(AP.ps %>% subset_samples(Primers=="nCOI") %>%
                                           prune_taxa(taxa_sums(.)>9,.) %>%
                                           subset_taxa(Kingdom!="Unassigned")
  )),
  readsCOI_assigned_fam=sum(sample_sums(AP.ps %>% subset_samples(Primers=="nCOI") %>%
                                          prune_taxa(taxa_sums(.)>9,.) %>%
                                          subset_taxa(!is.na(Family))
  ))
)

af.reads.2assgn.df = data.frame(
  reads18S_assigned_king=sum(sample_sums(AF.ps %>% subset_samples(Primers=="18S") %>%
                                           prune_taxa(taxa_sums(.)>9,.) %>%
                                           subset_taxa(Kingdom!="Unassigned")
  )),
  reads18S_assigned_fam=sum(sample_sums(AF.ps %>% subset_samples(Primers=="18S") %>%
                                          prune_taxa(taxa_sums(.)>9,.) %>%
                                          subset_taxa(!is.na(Family))
  )),
  reads12S_assigned_king=sum(sample_sums(AF.ps %>% subset_samples(Primers=="12S") %>%
                                           prune_taxa(taxa_sums(.)>9,.) %>%
                                           subset_taxa(Kingdom!="Unassigned")
  )),
  reads12S_assigned_fam=sum(sample_sums(AF.ps %>% subset_samples(Primers=="12S") %>%
                                          prune_taxa(taxa_sums(.)>9,.) %>%
                                          subset_taxa(!is.na(Family))
  )),
  readsCOI_assigned_king=sum(sample_sums(AF.ps %>% subset_samples(Primers=="nCOI") %>%
                                           prune_taxa(taxa_sums(.)>9,.) %>%
                                           subset_taxa(Kingdom!="Unassigned")
  )),
  readsCOI_assigned_fam=sum(sample_sums(AF.ps %>% subset_samples(Primers=="nCOI") %>%
                                          prune_taxa(taxa_sums(.)>9,.) %>%
                                          subset_taxa(!is.na(Family))
  ))
)

reads.2assgn.names <- c("18S: Number of reads assigned to Kingdom",
                        "18S: Number of reads assigned to Family",
                        "12S: Number of reads assigned to Kingdom",
                        "12S: Number of reads assigned to Family",
                        "COI: Number of reads assigned to Kingdom",
                        "COI: Number of reads assigned to Family")
reads.2assgn.table <- bind_rows(
  "A1"=af.reads.2assgn.df,
  "A2"=ap.reads.2assgn.df,
  "A3"=awf.reads.2assgn.df,
  "R1"=ag.reads.2assgn.df,
  "R2"=ah.reads.2assgn.df,
  "R3"=awr.reads.2assgn.df) %>% t() %>% 
  as.data.frame() %>% 
  rename("A1"="V1","A2"="V2","A3"="V3",
         "R1"="V4","R2"="V5","R3"="V6") %>% bind_cols(Reads=reads.2assgn.names,.)

ST.reads.2assgn.table <- reads.2assgn.table %>% 
  flextable(.) %>% bold(part="header") %>%
  set_caption("Supplementary Table: Number of reads in consensus sequences assigned to taxa at Kingdom and Family level, separated by primer set")



################# Table 3 Jaccard PERMANOVA  ################    
border_sep <- officer::fp_border(color="black",style="thick",width=2)
border_alt <- officer::fp_border(color="black",style="dotted")

permanova.jacc.stmo.df <- NCeDNA.comb.pa.jaccard.adonis %>% as.matrix() %>% as.data.frame() %>% rownames_to_column(var="Variable")

permanova.jacc.stmo.table <- permanova.jacc.stmo.df %>% 
  mutate(across(where(is.double), ~round(.,digits=3))) %>% 
  flextable() %>% bold(part="header")  %>%
  hline(i=5,part="body",border=border_alt) %>% #dotted sep for residuals
  hline(i=4,part="body",border=border_alt) %>% #dotted sep for residuals 
  hline(i=3,part = "body",border=border_sep)  %>% #thick border to sep total 
  bold(i = ~  permanova.jacc.stmo.df[6] < 0.05 & !is.na(permanova.jacc.stmo.df[6]),j="Pr(>F)", bold=TRUE) 



#################  Table 4 Raup-Crick PERMANOVA ################    
permanova.raup.stmo.df <- NCeDNA.comb.pa.raup.adonis %>% as.matrix() %>% as.data.frame() %>% rownames_to_column(var="Variable")

permanova.raup.stmo.table <- permanova.raup.stmo.df %>%  mutate(across(where(is.double), ~round(.,digits=3))) %>% flextable() %>% bold(part="header")  %>%
  hline(i=5,part="body",border=border_alt) %>% #dotted sep for residuals
  hline(i=4,part="body",border=border_alt) %>% #dotted sep for residuals 
  hline(i=3,part = "body",border=border_sep)  %>% #thick border to sep total 
  bold(i = ~  permanova.raup.stmo.df[6] < 0.05 & !is.na(permanova.raup.stmo.df[6]),j="Pr(>F)", bold=TRUE) 



################# S3 primer demultiplex to consensus ################    

ah.reads.2con.df = data.frame(
  reads18S_ntaxa=sum(ntaxa(AH.ps %>% subset_samples(Primers=="18S") %>%
                             prune_taxa(taxa_sums(.)>9,.)
  )),
  reads18S_sumseq=sum(sample_sums(AH.ps %>% subset_samples(Primers=="18S") %>%
                                    prune_taxa(taxa_sums(.)>9,.) 
  )),
  reads12S_ntaxa=sum(ntaxa(AH.ps %>% subset_samples(Primers=="12S") %>%
                             prune_taxa(taxa_sums(.)>9,.) 
  )),
  reads12S_sumseq=sum(sample_sums(AH.ps %>% subset_samples(Primers=="12S") %>%
                                    prune_taxa(taxa_sums(.)>9,.) 
  )),
  readsCOI_ntaxa=sum(ntaxa(AH.ps %>% subset_samples(Primers=="nCOI") %>%
                             prune_taxa(taxa_sums(.)>9,.)
  )),
  readsCOI_sumseq=sum(sample_sums(AH.ps %>% subset_samples(Primers=="nCOI") %>%
                                    prune_taxa(taxa_sums(.)>9,.)
  ))
)

ag.reads.2con.df = data.frame(
  reads18S_ntaxa=sum(ntaxa(AG.ps %>% subset_samples(Primers=="18S") %>%
                             prune_taxa(taxa_sums(.)>9,.) 
  )),
  reads18S_sumseq=sum(sample_sums(AG.ps %>% subset_samples(Primers=="18S") %>%
                                    prune_taxa(taxa_sums(.)>9,.) 
  )),
  reads12S_ntaxa=sum(ntaxa(AG.ps %>% subset_samples(Primers=="MiFish") %>%
                             prune_taxa(taxa_sums(.)>9,.)
  )),
  reads12S_sumseq=sum(sample_sums(AG.ps %>% subset_samples(Primers=="MiFish") %>%
                                    prune_taxa(taxa_sums(.)>9,.) 
  )),
  readsCOI_ntaxa=sum(ntaxa(AG.ps %>% subset_samples(Primers=="nCOI") %>%
                             prune_taxa(taxa_sums(.)>9,.)
  )),
  readsCOI_sumseq=sum(sample_sums(AG.ps %>% subset_samples(Primers=="nCOI") %>%
                                    prune_taxa(taxa_sums(.)>9,.) 
  ))
)


awr.reads.2con.df = data.frame(
  reads18S_ntaxa=sum(ntaxa(AWR.ps %>% subset_samples(Primers=="18S") %>%
                             prune_taxa(taxa_sums(.)>9,.)
  )),
  reads18S_sumseq=sum(sample_sums(AWR.ps %>% subset_samples(Primers=="18S") %>%
                                    prune_taxa(taxa_sums(.)>9,.) 
  )),
  reads12S_ntaxa=sum(ntaxa(AWR.ps %>% subset_samples(Primers=="12S") %>%
                             prune_taxa(taxa_sums(.)>9,.) 
  )),
  reads12S_sumseq=sum(sample_sums(AWR.ps %>% subset_samples(Primers=="12S") %>%
                                    prune_taxa(taxa_sums(.)>9,.) 
  )),
  readsCOI_ntaxa=sum(ntaxa(AWR.ps %>% subset_samples(Primers=="nCOI") %>%
                             prune_taxa(taxa_sums(.)>9,.) 
  )),
  readsCOI_sumseq=sum(sample_sums(AWR.ps %>% subset_samples(Primers=="nCOI") %>%
                                    prune_taxa(taxa_sums(.)>9,.) 
  ))
)


awf.reads.2con.df = data.frame(
  reads18S_ntaxa=sum(ntaxa(AWF.ps %>% subset_samples(Primers=="18S") %>%
                             prune_taxa(taxa_sums(.)>9,.) 
  )),
  reads18S_sumseq=sum(sample_sums(AWF.ps %>% subset_samples(Primers=="18S") %>%
                                    prune_taxa(taxa_sums(.)>9,.) 
  )),
  reads12S_ntaxa=sum(ntaxa(AWF.ps %>% subset_samples(Primers=="12S") %>%
                             prune_taxa(taxa_sums(.)>9,.) 
  )),
  reads12S_sumseq=sum(sample_sums(AWF.ps %>% subset_samples(Primers=="12S") %>%
                                    prune_taxa(taxa_sums(.)>9,.) 
  )),
  readsCOI_ntaxa=sum(ntaxa(AWF.ps %>% subset_samples(Primers=="nCOI") %>%
                             prune_taxa(taxa_sums(.)>9,.) 
  )),
  readsCOI_sumseq=sum(sample_sums(AWF.ps %>% subset_samples(Primers=="nCOI") %>%
                                    prune_taxa(taxa_sums(.)>9,.) 
  ))
)


ap.reads.2con.df = data.frame(
  reads18S_ntaxa=sum(ntaxa(AP.ps %>% subset_samples(Primers=="18S" | is.na(Primers)) %>%
                             prune_taxa(taxa_sums(.)>9,.) 
  )),
  reads18S_sumseq=sum(sample_sums(AP.ps %>% subset_samples(Primers=="18S" | is.na(Primers)) %>%
                                    prune_taxa(taxa_sums(.)>9,.) 
  )),
  reads12S_ntaxa=sum(ntaxa(AP.ps %>% subset_samples(Primers=="12S") %>%
                             prune_taxa(taxa_sums(.)>9,.) 
  )),
  reads12S_sumseq=sum(sample_sums(AP.ps %>% subset_samples(Primers=="12S") %>%
                                    prune_taxa(taxa_sums(.)>9,.) 
  )),
  readsCOI_ntaxa=sum(ntaxa(AP.ps %>% subset_samples(Primers=="nCOI") %>%
                             prune_taxa(taxa_sums(.)>9,.)
  )),
  readsCOI_sumseq=sum(sample_sums(AP.ps %>% subset_samples(Primers=="nCOI") %>%
                                    prune_taxa(taxa_sums(.)>9,.) 
  ))
)


af.reads.2con.df = data.frame(
  reads18S_ntaxa=sum(ntaxa(AF.ps %>% subset_samples(Primers=="18S") %>%
                             prune_taxa(taxa_sums(.)>9,.) 
  )),
  reads18S_sumseq=sum(sample_sums(AF.ps %>% subset_samples(Primers=="18S") %>%
                                    prune_taxa(taxa_sums(.)>9,.) 
  )),
  reads12S_ntaxa=sum(ntaxa(AF.ps %>% subset_samples(Primers=="12S") %>%
                             prune_taxa(taxa_sums(.)>9,.) 
  )),
  reads12S_sumseq=sum(sample_sums(AF.ps %>% subset_samples(Primers=="12S") %>%
                                    prune_taxa(taxa_sums(.)>9,.) 
  )),
  readsCOI_ntaxa=sum(ntaxa(AF.ps %>% subset_samples(Primers=="nCOI") %>%
                             prune_taxa(taxa_sums(.)>9,.) 
  )),
  readsCOI_sumseq=sum(sample_sums(AF.ps %>% subset_samples(Primers=="nCOI") %>%
                                    prune_taxa(taxa_sums(.)>9,.)
  ))
)

reads.2con.names <- c("18S: Number of consensus sequences",
                      "18S: Number of reads in consensus sequences",
                      "12S: Number of consensus sequences",
                      "12S: Number of reads in consensus sequences",
                      "COI: Number of consensus sequences",
                      "COI: Number of reads in consensus sequences")

reads.2con.table <- bind_rows(
  "A1"=af.reads.2con.df,
  "A2"=ap.reads.2con.df,
  "A3"=awf.reads.2con.df,
  "R1"=ag.reads.2con.df,
  "R2"=ah.reads.2con.df,
  "R3"=awr.reads.2con.df) %>% t() %>% 
  as.data.frame() %>% 
  rename("A1"="V1","A2"="V2","A3"="V3",
         "R1"="V4","R2"="V5","R3"="V6") %>% bind_cols(Reads=reads.2con.names,.)

demuxfilt <- read.csv("Primers_Demux_PassFilter.csv")

ST.demuxfilt <- demuxfilt %>% rename("A1"="AF","A2"="AP","A3"="AWF",
                                     "R1"="AG","R2"="AH","R3"="AWR") %>% 
  mutate(Reads=paste0(X.1,": ",X) %>% word(.,start=1,sep="\\(")) %>% 
  select(c("Reads","A1","A2","A3","R1","R2","R3")) %>% 
  bind_rows(reads.2con.table) %>%
  mutate(order=case_when(str_detect(Reads,"Assigned") ~ "1",
                         str_detect(Reads,"passed") ~ "2",
                         str_detect(Reads,"Number of reads") ~ "4",
                         str_detect(Reads, "Number of consensus") ~ "3"),
         primer=case_when(str_detect(Reads, "COI") ~ "1",
                          str_detect(Reads, "18S") ~ "2",
                          str_detect(Reads, "12S") ~ "3")) %>% 
  group_by(primer) %>% 
  arrange(order,.by_group = TRUE) %>% ungroup() %>%
  select(-c("order","primer"))

ST.demuxfilt.table <- ST.demuxfilt %>% 
  flextable(.) %>% bold(part="header") %>%
  bg(i=c(4,8,12),bg="grey") %>%
  hline(i=c(4,8,12),part = "body",border=border_sep) %>%
  set_caption("Supplementary Table:Progression of sequencing read numbers for each primer set (COI, 18S, and 12S). Each read was basecalled using the Oxford Nanopore Technology basecalling software and demultiplexed. Sequences that passed quality and length filtering were then clustered into consensus sequences and polished before assigning taxonomy.")

################# S4-6 Negative control distribution table ################    
nc.r1 <- c("Number of reads",
           paste0(round(mean(blank.df$reads),digits=3),"±",round(sd(blank.df$reads),digits=3)),
           paste0(range(blank.df$reads)[1],"-",round(range(blank.df$reads)[2],digits=0))
)
nc.r2 <- c("Number of consensus sequences",
           paste0(round(mean(blank.df$ncons),digits=3),"±",round(sd(blank.df$ncons),digits=3)),
           paste0(range(blank.df$ncons)[1],"-",round(range(blank.df$ncons)[2],digits=0))
)
nc.r3 <- c("Percent of reads in consensus sequences",
           paste0(round(mean(blank.df$perc_incon),digits=3),"±",round(sd(blank.df$perc_incon),digits=3),"%"),
           paste0(range(blank.df$perc_incon)[1],"-",round(range(blank.df$perc_incon)[2],digits=0),"%")
)
nc.r4 <- c("Percent of reads assigned",
           paste0(round(mean(blank.df$perc_assign),digits=3),"±",round(sd(blank.df$perc_assign),digits=3),"%"),
           paste0(range(blank.df$perc_assign)[1],"-",round(range(blank.df$perc_assign)[2],digits=0),"%")
)
negcontrol.table <- as.data.frame(rbind(nc.r1,nc.r2,nc.r3,nc.r4)) %>% rename(" " = "V1",
                                                                             "Mean±SD" = "V2",
                                                                             "Range" = "V3")
ST.negcontrol.table <- flextable(negcontrol.table) %>% 
  bold(part = "header") %>% vline() %>%
  set_caption("Supplementary Table: Measures of central tendency for negative control and PCR blank samples")

ST.negcontrol.A3.table <- negcontrol.A3.df %>% select(c(NewFam,NewPerc)) %>% 
  arrange(desc(NewPerc)) %>% 
  replace_na(list(NewFam="Not assigned to Family level") ) %>% 
  rename("Family"="NewFam",
         "Percent of reads"="NewPerc") %>% 
  flextable(.) %>% bold(part="header") %>%
  set_caption("Supplementary Table: Percent composition of taxonomic assignments for negative control sample in A3")

ST.negcontrol.R2.table <- negcontrol.R2.df %>% select(c(NewFam,NewPerc)) %>% 
  arrange(desc(NewPerc)) %>% mutate(Fc=as.character(NewFam),
                                    Family=ifelse(is.na(Fc),"Not assigned",Fc),.before=NewFam) %>% select(-c(NewFam,Fc)) %>%
  rename("Percent of reads"="NewPerc") %>% 
  flextable(.) %>% bold(part="header") %>%
  set_caption("Supplementary Table: Percent composition of taxonomic assignments for negative control sample in R2")



################# S7 Reads sequenced by site group ################     
ah.all.df = data.frame(reads=sample_sums(AH.ps), #number of reads per sample (same as sum(taxa_sums()))
                       ntaxa=colSums(otu_table(AH.ps)>0), #number of taxa per sample (same as ntaxa if you prune taxa>0)
                       reads_ncons=sample_sums(otu_table(AH.ps %>% #number of reads in consensus sequences 
                                                           prune_taxa(taxa_sums(.)>9,.) 
                       )),
                       ncons=colSums(otu_table(AH.ps %>% #number of consensus sequences  
                                                 prune_taxa(taxa_sums(.)>9,.)
                       )>0),  
                       
                       reads_assigned=sample_sums(AH.ps %>% 
                                                    prune_taxa(taxa_sums(.)>9,.) %>% 
                                                    subset_taxa(Kingdom!="Unassigned")), #number of assigned reads
                       ncons_assigned=colSums(otu_table(AH.ps %>% 
                                                          prune_taxa(taxa_sums(.)>9,.) %>% 
                                                          subset_taxa(Kingdom!="Unassigned")
                       )>0), #number of assigned consensus sequences
                       sample_data(AH.ps))
ag.all.df  = data.frame(reads=sample_sums(AG.ps), #number of reads per sample (same as sum(taxa_sums()))
                        ntaxa=colSums(otu_table(AG.ps)>0), #number of taxa per sample (same as ntaxa if you prune taxa>0)
                        reads_ncons=sample_sums(otu_table(AG.ps %>% #number of reads in consensus sequences 
                                                            prune_taxa(taxa_sums(.)>9,.) 
                        )),
                        ncons=colSums(otu_table(AG.ps %>% #number of consensus sequences  
                                                  prune_taxa(taxa_sums(.)>9,.)
                        )>0),  
                        
                        reads_assigned=sample_sums(AG.ps %>% 
                                                     prune_taxa(taxa_sums(.)>9,.) %>% 
                                                     subset_taxa(Kingdom!="Unassigned")), #number of assigned reads
                        ncons_assigned=colSums(otu_table(AG.ps %>% 
                                                           prune_taxa(taxa_sums(.)>9,.) %>% 
                                                           subset_taxa(Kingdom!="Unassigned")
                        )>0), #number of assigned consensus sequences
                        sample_data(AG.ps))
awr.all.df  = data.frame(reads=sample_sums(AWR.ps), #number of reads per sample (same as sum(taxa_sums()))
                         ntaxa=colSums(otu_table(AWR.ps)>0), #number of taxa per sample (same as ntaxa if you prune taxa>0)
                         reads_ncons=sample_sums(otu_table(AWR.ps %>% #number of reads in consensus sequences 
                                                             prune_taxa(taxa_sums(.)>9,.) 
                         )),
                         ncons=colSums(otu_table(AWR.ps %>% #number of consensus sequences  
                                                   prune_taxa(taxa_sums(.)>9,.)
                         )>0),  
                         
                         reads_assigned=sample_sums(AWR.ps %>% 
                                                      prune_taxa(taxa_sums(.)>9,.) %>% 
                                                      subset_taxa(Kingdom!="Unassigned")), #number of assigned reads
                         ncons_assigned=colSums(otu_table(AWR.ps %>% 
                                                            prune_taxa(taxa_sums(.)>9,.) %>% 
                                                            subset_taxa(Kingdom!="Unassigned")
                         )>0), #number of assigned consensus sequences
                         sample_data(AWR.ps))
ap.all.df = data.frame(reads=sample_sums(AP.ps), #number of reads per sample (same as sum(taxa_sums()))
                       ntaxa=colSums(otu_table(AP.ps)>0), #number of taxa per sample (same as ntaxa if you prune taxa>0)
                       reads_ncons=sample_sums(otu_table(AP.ps %>% #number of reads in consensus sequences 
                                                           prune_taxa(taxa_sums(.)>9,.) 
                       )),
                       ncons=colSums(otu_table(AP.ps %>% #number of consensus sequences  
                                                 prune_taxa(taxa_sums(.)>9,.)
                       )>0),  
                       
                       reads_assigned=sample_sums(AP.ps %>% 
                                                    prune_taxa(taxa_sums(.)>9,.) %>% 
                                                    subset_taxa(Kingdom!="Unassigned")), #number of assigned reads
                       ncons_assigned=colSums(otu_table(AP.ps %>% 
                                                          prune_taxa(taxa_sums(.)>9,.) %>% 
                                                          subset_taxa(Kingdom!="Unassigned")
                       )>0), #number of assigned consensus sequences
                       sample_data(AP.ps))
ap.all.df$Primers[is.na(ap.all.df$Primers)] <- "18S"
af.all.df = data.frame(reads=sample_sums(AF.ps), #number of reads per sample (same as sum(taxa_sums()))
                       ntaxa=colSums(otu_table(AF.ps)>0), #number of taxa per sample (same as ntaxa if you prune taxa>0)
                       reads_ncons=sample_sums(otu_table(AF.ps %>% #number of reads in consensus sequences 
                                                           prune_taxa(taxa_sums(.)>9,.) 
                       )),
                       ncons=colSums(otu_table(AF.ps %>% #number of consensus sequences  
                                                 prune_taxa(taxa_sums(.)>9,.)
                       )>0),  
                       
                       reads_assigned=sample_sums(AF.ps %>% 
                                                    prune_taxa(taxa_sums(.)>9,.) %>% 
                                                    subset_taxa(Kingdom!="Unassigned")), #number of assigned reads
                       ncons_assigned=colSums(otu_table(AF.ps %>% 
                                                          prune_taxa(taxa_sums(.)>9,.) %>% 
                                                          subset_taxa(Kingdom!="Unassigned")
                       )>0), #number of assigned consensus sequences
                       sample_data(AF.ps))
awf.all.df = data.frame(reads=sample_sums(AWF.ps), #number of reads per sample (same as sum(taxa_sums()))
                        ntaxa=colSums(otu_table(AWF.ps)>0), #number of taxa per sample (same as ntaxa if you prune taxa>0)
                        reads_ncons=sample_sums(otu_table(AWF.ps %>% #number of reads in consensus sequences 
                                                            prune_taxa(taxa_sums(.)>9,.) 
                        )),
                        ncons=colSums(otu_table(AWF.ps %>% #number of consensus sequences  
                                                  prune_taxa(taxa_sums(.)>9,.)
                        )>0),  
                        
                        reads_assigned=sample_sums(AWF.ps %>% 
                                                     prune_taxa(taxa_sums(.)>9,.) %>% 
                                                     subset_taxa(Kingdom!="Unassigned")), #number of assigned reads
                        ncons_assigned=colSums(otu_table(AWF.ps %>% 
                                                           prune_taxa(taxa_sums(.)>9,.) %>% 
                                                           subset_taxa(Kingdom!="Unassigned")
                        )>0), #number of assigned consensus sequences
                        sample_data(AWF.ps))

all.df <- bind_rows(ah.all.df, #r2
                    ag.all.df, #r1
                    awr.all.df, #r3
                    af.all.df, #a1
                    ap.all.df, #a2
                    awf.all.df, #a3
                    .id="id") 

all.df <- all.df %>% mutate(Primers=case_when(Primers=="18S" ~ "18S",
                                              Primers=="nCOI" ~ "COI",
                                              Primers=="12S" ~ "12S",
                                              Primers=="MiFish" ~ "12S"),
                            id=case_when(id=="1" ~ "R2",
                                         id=="2" ~ "R1",
                                         id=="3" ~ "R3",
                                         id=="4" ~ "A1",
                                         id=="5" ~ "A2",
                                         id=="6" ~ "A3",
                            ),
                            Location = case_when(Location == "G" ~ "R1",
                                                 Location == "H" ~ "R2",
                                                 Location == "WR" ~ "R3",
                                                 Location == "F" ~ "A1",
                                                 Location == "P" ~ "A2",
                                                 Location == "WF" ~ "A3",
                                                 TRUE ~ as.character(Location)))

NumberOfReadsPerSamplePlot <- all.df %>% filter(Location!="CONTROL") %>% 
  ggplot(., aes(x = Location, y = reads, color = reads_ncons)) +
  geom_boxplot(color="black") + theme_bw() +
  geom_jitter(width=.2, height=0)  + 
  facet_wrap(~Primers) +
  theme(axis.text.x = element_text(angle = 90,hjust = 1), text=element_text(family="Times")) +
  scale_color_continuous(name="Assigned consensus sequences",breaks=c(50000,100000,150000)) +
  labs(x="Site",y="Sequenced Reads",title="Reads sequenced by site group",subtitle="and the number of taxonomically assigned consensus sequences") 



################# S8 Negative control read distribution and tax assignment   ################    
NumberOfReadsPerControlPlot <-  all.df %>% 
  filter(Location=="CONTROL" & !str_detect(SampleName,"Positive")) %>% 
  mutate(Label=case_when(reads>20000 ~ id) ) %>% 
  ggplot(., aes(y = reads, x = Primers, color = ncons_assigned)) + 
  geom_boxplot(color="black") + theme_bw() +
  geom_point(aes(color = ncons_assigned), alpha = 1,position=position_jitter(width=0.15,height=0),size=3) +
  scale_color_continuous(name="Assigned \n consensus sequences")  + 
  theme(axis.text.x = element_text(angle = 90), text=element_text(family="Times"),plot.margin = unit(c(0,0,2,0),"line")) +
  labs(title="Distribution of reads in negative control samples",subtitle="and the number of assigned consensus sequences",y="Sequenced Reads") +
  geom_label_repel(aes(label=Label),direction=c("y"),point.padding=0.25)

NumberOfReadsPerControlPlot 
negcontrol.R2.df <- AH.ps %>% subset_samples(str_detect(SampleName,"EB") & Primers=="12S") %>% 
  prune_taxa(taxa_sums(.)>0,.) %>% 
  tax_glom("Family",NArm=FALSE) %>% psmelt() %>%
  select(c(OTU,SampleName,Abundance,Kingdom,Phylum,Class,Order,Family)) %>% 
  mutate(Perc=Abundance/sum(Abundance),
         NewFam=case_when(Perc>0.01~Family,Perc<0.01 ~ "Taxa below 1% of sequences")) %>% group_by(NewFam) %>% reframe(NewPerc=sum(Perc),PerAbun=as.factor(NewPerc),NewAbun=sum(Abundance))  

negcontrol.R2.df$NewFam <- factor(negcontrol.R2.df$NewFam, 
                                  levels=c(NA," Naviculaceae"," Hominidae"," Mugilidae"," Fundulidae","Taxa below 1% of sequences"))

negcontrol.R2.plot <- negcontrol.R2.df %>%
  ggplot(.,aes(fill=NewFam,y=NewFam,x=NewAbun)) + geom_bar(stat="identity",position="stack") + theme_bw() +
  theme(text=element_text(family="Times"))  + labs(y="",x="Abundance",title="R2 Negative control",subtitle="Family-level taxonomic assignments") + theme(legend.position="none") + scale_x_continuous(expand=expansion(mult=c(0,0.1)))

negcontrol.A3.df <- AWF.ps %>% subset_samples(Location=="CONTROL" & Primers=="12S") %>% 
  prune_taxa(taxa_sums(.)>0,.) %>% 
  tax_glom("Family",NArm=FALSE) %>% psmelt() %>%
  select(c(OTU,SampleName,Abundance,Kingdom,Phylum,Class,Order,Family)) %>% 
  mutate(Perc=Abundance/sum(Abundance),
         NewFam=case_when(Perc>0.01~Family,Perc<0.01 ~ "Taxa below 1% of sequences")) %>% group_by(NewFam) %>% reframe(NewPerc=sum(Perc),PerAbun=as.factor(NewPerc),NewAbun=sum(Abundance))  

negcontrol.A3.plot <- negcontrol.A3.df %>%
  ggplot(.,aes(fill=NewFam,y=fct_reorder(NewFam,NewAbun),x=NewAbun)) + geom_bar(stat="identity",position="stack") + theme_bw() +
  theme(text=element_text(family="Times"))  + labs(y="",x="Abundance",title="A3 Negative control",subtitle="Family-level taxonomic assignments") + theme(legend.position="none") + scale_x_continuous(expand=expansion(mult=c(0,0.1)))



#################  S14 Top 25 families by season   ################     
df1 = data.frame(col = sort(SharedSpringFams), Spring=sort(SharedSpringFams))
df2 = data.frame(col = sort(sort(SharedSummerFams)), Summer=sort(sort(SharedSummerFams)))
df3 = data.frame(col = sort(SharedFallFams), Fall=sort(SharedFallFams))
shared.df <- Reduce(function(x, y) merge(x, y, all=TRUE), list(df1, df2, df3))

ST.sharedfams <- flextable(shared.df %>% rename("Family" = "col")) %>% 
  bold(part = "header") %>% 
  vline() %>% hline() %>% 
  bg(i =~ is.na(Spring),j =~ Spring,bg="grey") %>% 
  bg(i =~ is.na(Summer),j =~ Summer,bg="grey") %>% 
  bg(i =~ is.na(Fall),j =~ Fall,bg="grey") %>% 
  set_caption("Supplementary Table: The number of clusters in each size group identified by CD-HIT for each site. ")



################# S15 Top 25 families by site type   ################    
df6 = data.frame(col = sort(toptaxa.aquaculture.all25), Aquaculture=sort(toptaxa.aquaculture.all25))
df5 = data.frame(col = sort(toptaxa.reef.all25), Reef=sort(toptaxa.reef.all25))

abundanttaxa.df <- Reduce(function(x, y) merge(x, y, all=TRUE), list(df5, df6))

abundanttaxa.df[2:3] = replace(data.frame(lapply(abundanttaxa.df[2:3], as.character), stringsAsFactors = FALSE),
                               !is.na(abundanttaxa.df[2:3]), "x")

ST.abundanttaxa <- flextable(abundanttaxa.df %>% 
                               rename("Family" = "col")) %>% 
  bold(part = "header") %>% vline() %>% hline() %>% 
  bg(i =~ is.na(Reef),j =~ Reef,bg="grey") %>%
  bg(i =~ is.na(Aquaculture),j =~ Aquaculture,bg="grey") %>%
  set_caption("Supplementary Table: The top 25 most abundant families identified at oyster reefs and at floating oyster aquaculture facilities.")



################# S16 NMDS for individual primers  ################    
NCeDNA.jaccard.nmds <- ordinate(
  physeq=NCeDNA.famglom.pa.ps,
  k=3,try=500,
  method = "NMDS", 
  distance = "jaccard"
)

primers.nmds.jacc.plot <- plot_ordination(NCeDNA.pa.ps,
                                          NCeDNA.jaccard.nmds,
                                          type="samples",shape="SiteType",color="SiteType") + facet_wrap(~Primers) + labs(title="NMDS on Jaccard distance for each individual primer set") + theme_bw() +
  geom_jitter(width=2, height=0.5)  + 
  theme(axis.text.x = element_text(angle = 90,hjust = 1), text=element_text(family="Times")) + scale_color_manual(values=SiteType.Palette.clrblnd)




#################  S17 Pairwise PERMANOVA Jaccard Sites  ################    
df <- pwintx_jaccard_sitemon[[2]] %>% as.matrix() %>% as.data.frame() %>% rownames_to_column(var="Variable") %>% mutate(Comparison= names(pwintx_jaccard_sitemon[2]),.before=Variable )
for (i in 3:length(pwintx_jaccard_sitemon)) {
  new <- pwintx_jaccard_sitemon[[i]] %>% as.matrix() %>% as.data.frame() %>% rownames_to_column(var="Variable") %>% mutate(Comparison= names(pwintx_jaccard_sitemon[i]),.before=Variable )
  df <- bind_rows(df,new)
}

pw.permanova.jacc.site <- df  %>%  mutate(across(where(is.double), ~round(.,digits=3))) %>% flextable() %>% bold(part="header")  %>%
  bg(.,i=seq(5,75,by=5), bg="#ccd7dd") %>%  #highlight total
  hline(i=seq(3,73,by=5),part="body",border=border_alt) %>% #dotted sep for residuals
  hline(i=seq(4,74,by=5),part="body",border=border_alt) %>% #dotted sep for residuals 
  hline(i=seq(5,75,by=5),part = "body",border=border_sep)  %>% #thick border to sep total 
  bold(i = ~  df[7] < 0.05 & !is.na(df[7]),j="Pr(>F)", bold=TRUE) 




#################  S18 Pairwise PERMANOVA Raup-Crick Sites    ################    
dfr <- pwintx_raup_sitemon[[2]] %>% as.matrix() %>% as.data.frame() %>% rownames_to_column(var="Variable") %>% mutate(Comparison= names(pwintx_raup_sitemon[2]),.before=Variable )
for (i in 3:length(pwintx_raup_sitemon)) {
  new <- pwintx_raup_sitemon[[i]] %>% as.matrix() %>% as.data.frame() %>% rownames_to_column(var="Variable") %>% mutate(Comparison= names(pwintx_raup_sitemon[i]),.before=Variable )
  dfr <- bind_rows(dfr,new)
}

pw.permanova.raup.site <- dfr  %>%  mutate(across(where(is.double), ~round(.,digits=3))) %>% flextable() %>% bold(part="header")  %>%
  bg(.,i=seq(5,75,by=5), bg="#ccd7dd") %>%  #highlight total
  hline(i=seq(3,73,by=5),part="body",border=border_alt) %>% #dotted sep for residuals
  hline(i=seq(4,74,by=5),part="body",border=border_alt) %>% #dotted sep for residuals 
  hline(i=seq(5,75,by=5),part = "body",border=border_sep) %>% #thick border to sep total 
  bold(i = ~  dfr[7] < 0.05 & !is.na(dfr[7]),j="Pr(>F)", bold=TRUE) 



#################   S19 Pairwise PERMANOVA Jaccard Months   ################    
pw.apr.may <- pwintx_jaccard_monst$May_vs_Apr %>% as.matrix() %>% as.data.frame() %>% rownames_to_column(var="Variable")
pw.may.jun <- pwintx_jaccard_monst$Jun_vs_May %>% as.matrix() %>% as.data.frame()%>% rownames_to_column(var="Variable")
pw.jun.jul <- pwintx_jaccard_monst$Jun_vs_Jul %>% as.matrix() %>% as.data.frame()%>% rownames_to_column(var="Variable")
pw.jul.aug <- pwintx_jaccard_monst$Aug_vs_Jul %>% as.matrix() %>% as.data.frame()%>% rownames_to_column(var="Variable")
pw.aug.sep <- pwintx_jaccard_monst$Sep_vs_Aug %>% as.matrix() %>% as.data.frame()%>% rownames_to_column(var="Variable")
pw.sep.oct <- pwintx_jaccard_monst$Sep_vs_Oct %>% as.matrix() %>% as.data.frame()%>% rownames_to_column(var="Variable")

pw.permanova.jacc.month.df <- 
  bind_rows(pw.apr.may,pw.may.jun,
            pw.jun.jul,pw.jul.aug,
            pw.aug.sep,pw.sep.oct,
            .id="Comparison") %>%
  mutate(Comparison=case_when(
    Comparison=="1" ~ "Apr-May",
    Comparison=="2" ~ "May-Jun",
    Comparison=="3" ~ "Jun-Jul",
    Comparison=="4" ~ "Jul-Aug",
    Comparison=="5" ~ "Aug-Sep",
    Comparison=="6" ~ "Sep-Oct"
  )) 
pw.permanova.jacc.month <- pw.permanova.jacc.month.df  %>%  mutate(across(where(is.double), ~round(.,digits=3))) %>% flextable() %>% bold(part="header")  %>%
  bg(.,i=seq(5,30,by=5), bg="#ccd7dd") %>%  #highlight total
  hline(i=seq(3,27,by=5),part="body",border=border_alt) %>% #dotted sep for residuals
  hline(i=seq(4,29,by=5),part="body",border=border_alt) %>% #dotted sep for residuals 
  hline(i=seq(5,30,by=5),part = "body",border=border_sep)  %>% #thick border to sep total 
  bold(i = ~  pw.permanova.jacc.month.df[7] < 0.05 & !is.na(pw.permanova.jacc.month.df[7]),j="Pr(>F)", bold=TRUE)



#################   S20 Pairwise PERMANOVA Raup-Crick Months    ################    
pw.raup.apr.may <- pwintx_raup_monst$May_vs_Apr %>% as.matrix() %>% as.data.frame() %>% rownames_to_column(var="Variable")
pw.raup.may.jun <- pwintx_raup_monst$Jun_vs_May %>% as.matrix() %>% as.data.frame()%>% rownames_to_column(var="Variable")
pw.raup.jun.jul <- pwintx_raup_monst$Jun_vs_Jul %>% as.matrix() %>% as.data.frame()%>% rownames_to_column(var="Variable")
pw.raup.jul.aug <- pwintx_raup_monst$Aug_vs_Jul %>% as.matrix() %>% as.data.frame()%>% rownames_to_column(var="Variable")
pw.raup.aug.sep <- pwintx_raup_monst$Sep_vs_Aug %>% as.matrix() %>% as.data.frame()%>% rownames_to_column(var="Variable")
pw.raup.sep.oct <- pwintx_raup_monst$Sep_vs_Oct %>% as.matrix() %>% as.data.frame()%>% rownames_to_column(var="Variable")

pw.permanova.raup.month.df <- 
  bind_rows(pw.raup.apr.may,pw.raup.may.jun,
            pw.raup.jun.jul,pw.raup.jul.aug,
            pw.raup.aug.sep,pw.raup.sep.oct,
            .id="Comparison") %>%
  mutate(Comparison=case_when(
    Comparison=="1" ~ "Apr-May",
    Comparison=="2" ~ "May-Jun",
    Comparison=="3" ~ "Jun-Jul",
    Comparison=="4" ~ "Jul-Aug",
    Comparison=="5" ~ "Aug-Sep",
    Comparison=="6" ~ "Sep-Oct"
  )) 

pw.permanova.raup.month <-pw.permanova.raup.month.df %>%  mutate(across(where(is.double), ~round(.,digits=3))) %>% flextable() %>% bold(part="header")  %>%
  bg(.,i=seq(5,30,by=5), bg="#ccd7dd") %>%  #highlight total
  hline(i=seq(3,27,by=5),part="body",border=border_alt) %>% #dotted sep for residuals
  hline(i=seq(4,29,by=5),part="body",border=border_alt) %>% #dotted sep for residuals 
  hline(i=seq(5,30,by=5),part = "body",border=border_sep) %>% #thick border to sep total 
  bold(i = ~  pw.permanova.raup.month.df[7] < 0.05 & !is.na(pw.permanova.raup.month.df[7]),j="Pr(>F)", bold=TRUE) 
