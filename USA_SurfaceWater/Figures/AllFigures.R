##########################FIG 1
##Upset Plot
library(circlize)
library(ComplexHeatmap)

allnames<-read.csv("GROWdb_with_vars_20220725_upsetPlot.csv",header=T,check.names = F)
row.names(allnames)=allnames$Number
m<-make_comb_mat(allnames)
UpSet(m,pt_size=unit(5,'mm'),lwd=3,comb_col=c("red","blue","black"))

library(ComplexUpset)
#sediment
sediment<-read.csv("GROWdb_with_vars_20220725_upsetPlot_v2.csv",header=T)
data = colnames(sediment)[4:7]

sed.plot=upset(sediment,data,sort_sets=FALSE, base_annotations=list('Intersection size'=intersection_size(counts=TRUE, mapping=aes(fill=Sample))),width_ratio=0.1)
sed.plot

##Map
library(sf)
library(mapview)
library(tidyverse)
library(tigris)
library(ggrepel)
library(scales)
library(ggsflabel)
library(extrafont)


sf_use_s2(FALSE)

# Shift and rescale but position AK/HI/PR outside the continental US rather than below
us_states_outside <- states(cb = TRUE, resolution = "20m") %>%
  shift_geometry(position = "below") %>%
  filter(!STUSPS %in% c("HI","PR"))

contus <- st_read('contus_xy.shp')

ak <- st_read('alaska_xy.shp')


grow <- read_csv('GROWdb_with_vars_20220725.csv') %>%
  mutate(MetaGcase=ifelse(!is.na(MetaG_JGI_accession),"MetaG",""),
         MetaTcase=ifelse(has_MetaT=="Y","MetaT","")) %>%
  mutate(Available_Meta=paste0(MetaGcase," ",MetaTcase)) %>%
  select(SampleName, Sample,Number,Latitude, Longitude, Available_Meta) %>%
  filter(SampleName %in% c(contus$SampleName, ak$SampleName)) %>%
  st_as_sf(coords=c("Longitude","Latitude"),crs=4326) %>%
  shift_geometry(position = "below")

nhd <- st_read('rivers_huc2.shp') %>% st_zm() %>%
  filter(!HUC2%in% c("21")) %>%
  shift_geometry(position = "below")

both <- dplyr::filter(grow, Available_Meta == "MetaG MetaT")

metaG <- filter(grow, Available_Meta == "MetaG ")

metaT <- filter(grow, Available_Meta == " MetaT")
name <- (grow)

theme_map <- function(...) {
  theme_minimal() +
    theme(
      text = element_text( color = "#22211d"),
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      # panel.grid.minor = element_line(color = "#ebebe5", size = 0.2),
      panel.grid.major = element_line(color = "white", size = 0.2),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "white", color = NA), 
      panel.background = element_rect(fill = "white", color = NA), 
      legend.background = element_rect(fill = "white", color = NA),
      panel.border = element_blank(),
      ...
    )
}


colors <- c('#4B0082', '#8A2BE2', '#FF00FF', '#DB7093', '#800000', 
            '#FF0000', '#FFA07A', '#FF8C00', '#FFD700', '#AFAA30',
            '#40E0D0', '#006400', '#3CB371', '#4B91CC', '#A1C0DE',
            '#94DEEC', '#936CDC', '#2B1776', '#191919') %>%          
  #grDevices::colors()[grep('gr(a|e)y|white', grDevices::colors(), invert = T)] %>%
  sample(.,19)

  
colors <- c("#FFA07A", "#936CDC", "#2B1776", "#FF00FF", "#8A2BE2", 
"#FFD700", "#40E0D0", "#006400", "#FF0000", "#191919", "#4B91CC", 
"#3CB371", "#FF8C00", "#4B0082", "#AFAA30", "#A1C0DE", "#94DEEC",
"#DB7093", "#800000")

ggplot() +
  theme_map() +
  geom_sf(data = us_states_outside, fill=NA, color="black", alpha = 0.01, size=0.25) +
  geom_sf(data = nhd, aes(color = HUC2), fill=NA, size = 0.20) +
  scale_color_manual(values = colors) +
  geom_sf(data=metaT, fill="orange", pch=21, color="black", size=5, alpha=0.5) +
  geom_sf(data=both, fill="blue", pch=21, color="black", size=5, alpha=0.5) +
  geom_sf(data=metaG, fill="red", pch=21, color="black", size=5, alpha=0.5) +
  geom_sf_label_repel(data=name, aes(label=Number),force = 50, nudge_x = -2, seed = 10, max.overlaps=60,size=3)
  
ggsave('maps/meta_map_092222.pdf')

##########################FIG 2
##See SRA analysis files 

#riverine metagenome compared to GROWdb stripchart Gbp
boxplot=read.csv("SRA_GROW_comparison_depth.csv")
#make plot
p <- ggplot(boxplot, aes(x=Source,y = Gbp)) + 
  geom_point(position=position_jitterdodge(jitter.width=0.75, dodge.width = 0),aes(color=Source,size=1))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90)) 

p

##Cladogram made on https://app.rawgraphs.io/

##########################FIG 3

library(ggplot2)
library(gridExtra)
library(reshape)
library("RColorBrewer")
#core figures genus
CoreScatter=read_csv("CoreFig.csv")
p=ggplot(CoreScatter,aes(x=MetaG_Mean,y=MetaG_Occupancyperc))+geom_point(aes(color=Phylum,size=.1))
p

####ribbon plot for metaG phyla

brewer.pal(n = 27, name = "Phylum")
[1] "#427638" "#FC8D62" "#8DA0CB" "#E78AC3" "#A6D854" "#FFD92F" "#E5C494" "#B3B3B3"

library(ggsankey)
library(tidyverse)
Phyla_Ribbon=read_csv("strict_mapping_table_95id_3x_60%cov_relabundance_phylum.csv")

Phyla_Ribbon=as.data.frame(Phyla_Ribbon)
Phyla_Ribbon.melt=melt(Phyla_Ribbon)
write.csv(Phyla_Ribbon.melt,file="Phyla_Ribbon_melt.csv")


p=ggplot(Phyla_Ribbon.melt, aes(x = Sample,
               node = variable,
               fill = variable,
               value = value)) +
  geom_sankey_bump(space = 0, type = "alluvial", color = "transparent", smooth = 15) +
  theme_sankey_bump(base_size = 16) +
  scale_color_manual(values=c("#427638","#F57F20")) +
  labs(x = NULL,
       y = "Relative Abundance",
       fill = "Phylum",
       color = NULL) +
  theme(legend.position = "bottom") 


p

pdf()
ggsave("filename.tiff")

####Genus metaT top 25 

top25=read_csv("RelAbund_norm.counts.rpk_edger.bins_mean_atLeast20_genustop25.csv")
top25=as.data.frame(top25)
top.25.melt=melt(top25)
p <- ggplot(top.25.melt, aes(x=Genus,y = value)) + geom_point(position=position_jitterdodge(jitter.width=0.55, dodge.width = 1),aes(color=Phylum,size=.05))+ geom_boxplot(position=position_dodge(2))+theme(axis.text.x = element_text(angle = 90)) 

p

##########################FIG4
library(vegan)
#read in gene table
Gene_table=read.csv("geTMM_norm.counts.rpk_edger_genes.csv")
Gene_table=as.data.frame(Gene_table)
#renaming genes names as rownames
rownames(Gene_table)<-Gene_table[,1]

Gene_table<-Gene_table[,-1]
#transpose
Gene_table=t(Gene_table)
#Read in MAG table 
MAG_table=read.csv("geTMM_norm.counts.rpk_edger_genome_atLeast20.csv")
#renaming samples as row names
rownames(MAG_table)<-MAG_table[,1]
MAG_table<-MAG_table[,-1]

######Diversity Stats
#Find Shannon diversity of MAGs
MAG_diversity=diversity(MAG_table, index="shannon", MARGIN=1, base= exp(1))
#write table to paste into excel
write.csv(MAG_diversity, file="MAG_diversity.csv")
#Find richness of all samples
MAG_richness=specnumber(MAG_table)
#write table for richness
write.csv(MAG_richness, file="MAG_richness.csv")

#Find Shannon diversity of Genes
Gene_diversity=diversity(Gene_table, index="shannon", MARGIN=1, base= exp(1))
#write table to paste into excel
write.csv(Gene_diversity, file="Gene_diversity.csv")
#Find richness of all samples
Gene_richness=specnumber(Gene_table)
#write table for richness
write.csv(Gene_richness, file="Gene_richness.csv")


#read in two tables to be joined, geospatial and geochem/microbial
Micro_Chem=read.csv("Correlation_MetaT57sample.csv")
row.names(Micro_Chem)=Micro_Chem[,1]
Micro_Chem=Micro_Chem[,-1]

Geospat=read.csv("GROWdb_with_vars_20220715_GeospatialForCorr.csv")
row.names(Geospat)=Geospat[,1]
Geospat=Geospat[,-1]

#join tables for corplot
merged_table=merge(Micro_Chem, Geospat,by = 'row.names')
####corr plot
library(RColorBrewer)
library(corrplot)
library(ggplot2)
col = (brewer.pal(n = 8, name = "RdGy"))

row.names(merged_table)=merged_table[,1]
merged_table=merged_table[,-1]
#correlation
cor=cor(Metadata_43, use='pairwise.complete.obs')
#make corrlation plot
corrplot(cor, type="upper")

#combine with a significance test
p.mat <- cor.mtest(cor, conf.level = .95)$p
#rerun correlation plot
corrplot(cor,  type="upper", p.mat = p.mat, sig.level = 0.05, insig = "blank", col=col, method='square')
ggsave("corrplot.pdf")
write.csv(cor,file="corrplot_43.csv")

#mantel tests to community data 
library(vegan)
library(dplyr)
######Mantel Test at 3 levels

###[1] MAG metagenomes

MAGs_142=read.csv("strict_mapping_table_95id_3x_60%cov_relabundance_142.csv")
row.names(MAGs_142)=MAGs_142[,1]
MAGs_142=MAGs_142[,-1]
Metadata_142=read.csv("GROWdb_with_vars_20220715_142_metaG.csv")
row.names(Metadata_142)=Metadata_142[,1]
Metadata_142=Metadata_142[,-1]
Metadata_142=Metadata_142[,-1]
MAGs_142_dist = dist(MAGs_142, method = "euclidean")

mantel_mags_sw = data.frame(matrix(nrow = ncol(Metadata_142),ncol = 3)) %>%
  rename(predictor = X1,mantelr = X2, pval= X3)

for (i in 1:ncol(Metadata_142)){
  out.clr = mantel(MAGs_142_dist ~ dist(Metadata_142[,i]))
  mantel_mags_sw[i,1] = colnames(Metadata_142[i]) 
  mantel_mags_sw[i,2] = out.clr[1]
  mantel_mags_sw[i,3] = out.clr[2] 
}
mantel_mags_sw

MAGs_115=read.csv("strict_mapping_table_95id_3x_60%cov_relabundance_115.csv")
row.names(MAGs_115)=MAGs_115[,1]
MAGs_115=MAGs_115[,-1]
Metadata_115=read.csv("GROWdb_with_vars_20220715_115_metaG.csv")
row.names(Metadata_115)=Metadata_115[,1]
Metadata_115=Metadata_115[,-1]
Metadata_115=Metadata_115[,-1]
Metadata_115_dist = dist(Metadata_115, method = "euclidean")

mantel_mags_sw = data.frame(matrix(nrow = ncol(Metadata_115),ncol = 3)) %>%
  rename(predictor = X1,mantelr = X2, pval= X3)

for (i in 1:ncol(Metadata_115)){
  out.clr = mantel(Metadata_115_dist ~ dist(Metadata_115[,i]), nperm = 10000, mrank = T)
  mantel_mags_sw[i,1] = colnames(Metadata_115[i]) 
  mantel_mags_sw[i,2] = out.clr[1]
  mantel_mags_sw[i,3] = out.clr[2] 
}
mantel_mags_sw

MAGs_50=read.csv("strict_mapping_table_95id_3x_60%cov_relabundance_50.csv")
row.names(MAGs_50)=MAGs_50[,1]
MAGs_50=MAGs_50[,-1]
Metadata_50=read.csv("GROWdb_with_vars_20220715_50_metaG.csv")
row.names(Metadata_50)=Metadata_50[,1]
Metadata_50=Metadata_50[,-1]
Metadata_50=Metadata_50[,-1]
Metadata_50_dist = dist(Metadata_50, method = "euclidean")

mantel_mags_sw = data.frame(matrix(nrow = ncol(Metadata_50),ncol = 3)) %>%
  rename(predictor = X1,mantelr = X2, pval= X3)

for (i in 1:ncol(Metadata_50)){
  out.clr = mantel(Metadata_50_dist ~ dist(Metadata_50[,i]), nperm = 10000, mrank = T)
  mantel_mags_sw[i,1] = colnames(Metadata_50[i]) 
  mantel_mags_sw[i,2] = out.clr[1]
  mantel_mags_sw[i,3] = out.clr[2] 
}
mantel_mags_sw

###[2] MAG MetaT

MAGs_52=read.csv("geTMM_norm.counts.rpk_edger_genome_relabund_52.csv")
row.names(MAGs_52)=MAGs_52[,1]
MAGs_52=MAGs_52[,-1]
Metadata_52=read.csv("corrplot_mergedtable_52.csv")
row.names(Metadata_52)=Metadata_52[,1]
Metadata_52=Metadata_52[,-1]
MAGs_52_dist = dist(MAGs_52, method = "euclidean")

mantel_mags_sw = data.frame(matrix(nrow = ncol(Metadata_52),ncol = 3)) %>%
  rename(predictor = X1,mantelr = X2, pval= X3)

for (i in 1:ncol(Metadata_52)){
  out.clr = mantel(MAGs_52_dist ~ dist(Metadata_52[,i]), nperm = 10000, mrank = T)
  mantel_mags_sw[i,1] = colnames(Metadata_52[i]) 
  mantel_mags_sw[i,2] = out.clr[1]
  mantel_mags_sw[i,3] = out.clr[2] 
}
mantel_mags_sw

MAGs_43=read.csv("geTMM_norm.counts.rpk_edger_genome_relabund_43.csv")

row.names(MAGs_43)=MAGs_43[,1]
MAGs_43=MAGs_43[,-1]
Metadata_43=read.csv("corrplot_mergedtable_43.csv")
row.names(Metadata_43)=Metadata_43[,1]
Metadata_43=Metadata_43[,-1]
MAGs_43_dist = dist(MAGs_43, method = "euclidean")

mantel_mags_sw = data.frame(matrix(nrow = ncol(Metadata_43),ncol = 3)) %>%
  rename(predictor = X1,mantelr = X2, pval= X3)

for (i in 1:ncol(Metadata_43)){
  out.clr = mantel(MAGs_43_dist ~ dist(Metadata_43[,i]), nperm = 10000, mrank = T)
  mantel_mags_sw[i,1] = colnames(Metadata_43[i]) 
  mantel_mags_sw[i,2] = out.clr[1]
  mantel_mags_sw[i,3] = out.clr[2] 
}
mantel_mags_sw

#stream order NMDS
OTU_table=read.csv("strict_mapping_table_95id_3x_60%cov_relabundance_105.csv")
#your OTU table should have samples as rows and otus as columns

#renaming your sample names as rownames
rownames(OTU_table)<-OTU_table[,1]

OTU_table<-OTU_table[,-1]

#Read in factors, for samples
Factors=read.csv("GROWdb_with_vars_20220715_105_metaG.csv")
#renaming your samples as row names
#note the order of this list need to be the same as the order in the OTU table
rownames(Factors)<-Factors[,1]
Factors<-Factors[,-1]

#NMDS
# NMDS streamorder and omernik with metaG 105
ord <- metaMDS(OTU_table, distance = "bray", k = 2,
               trymax = 200)
stressplot(ord)
nmds1 <- ord$points[,1]


nmds2 <- ord$points[,2]
ord_plot <- as.data.frame(cbind(nmds1,nmds2))

plt.nmds_ord = ggplot(ord_plot, aes(nmds1, nmds2)) +
    geom_point(aes(color = Factors$Omernik2),
             size = 3) +
  theme_classic() 
plt.nmds_ord

mrpp<-mrpp(OTU_table, Factors$Omernik2, permutations=999, distance="bray")

summary(mrpp)


mrpp

##########################FIG5
#see repo https://github.com/jmikayla1991/Genome-Resolved-Open-Watersheds-database-GROWdb/tree/main/USA_SurfaceWater/GenomeBasedCarbonAnalyses


  

