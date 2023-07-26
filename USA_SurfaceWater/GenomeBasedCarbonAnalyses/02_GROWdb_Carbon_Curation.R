library(readxl)
library(writexl)
library(dplyr)
library(tidyr)
library(FSA)
library(ggplot2)
library(ggalluvial)

metaG_adj= read_excel("Adjectives_manuallycalled.xlsx",sheet="for_metaG")
metaT_adj= read_excel("Adjectives_manuallycalled.xlsx", sheet="metaT_genes")
geTMM=read.delim("metaT_output_feature_counts_revstranded_07192022_97perc_sample_filterd_geTMM.csv",sep=",",header=TRUE)
tidy_geTMM = geTMM %>% gather(-X,key="fasta",value="geTMM")%>%filter(geTMM>0) %>% separate(fasta,into=c("river","code"),sep="_WHONDRS.")%>%separate(code,into=c("code","junk"),sep="_MT_bowtie.reformat97_sorted.bam")%>%select(-junk)
river_info = read.delim("Merged_Metadata3_allData_continuousOnly.csv",sep=",",header = TRUE)

# filter to carbon fixation pathways >=70% complete
c_fix_metaG= metaG_adj %>% select(-Taxonomy,-Photosynthetic,-Methanotroph,-Nitrifier,-`S Oxidizer`)%>%gather(-genome,key="c_fix_path",value="comp")%>%filter(comp>=0.70)
metaG_adj=metaG_adj %>% select(genome,Taxonomy,Photosynthetic,Methanotroph,Nitrifier,`S Oxidizer`)

metaT_fx=metaT_adj%>%separate(gene,into=c("MAG","junk"),remove=FALSE,sep=c("_scaff|_k121|_Ga0"))%>%
  select(-junk)%>%
  left_join(.,metaG_adj,by=c("MAG"="genome"))%>%
  select(gene,MAG,final_ID,module,c,c_type,photo,Guild,Taxonomy,Photosynthetic,Methanotroph,Nitrifier,`S Oxidizer`)%>%
  left_join(.,c_fix_metaG,by=c("MAG"="genome","module"="c_fix_path"))%>%
  mutate(guild2=ifelse(Guild=="c_fixation",ifelse(is.na(comp)==0,"central_c","c_fixation"),Guild))%>%select(-comp)%>%
  mutate(guild3=ifelse(is.na(guild2)==1,"other",Guild))%>%
  select(-Photosynthetic,-Methanotroph,-Nitrifier,-`S Oxidizer`)
  
MAG_metaT_guild=metaT_fx%>%left_join(.,tidy_geTMM,by=c("gene"="X"))%>%group_by(MAG,Taxonomy,guild3,photo,river,code)%>%summarise(sum=sum(geTMM))%>%ungroup()%>%select(-Taxonomy, -photo)%>%group_by(MAG,river,code,guild3)%>%summarise(sum=sum(sum))%>%
  mutate(life=ifelse(guild3=="photo","photo",
                     ifelse(guild3=="chemolithoautotroph","litho",
                            ifelse(guild3=="heterotroph-aromatic","organic",
                                   ifelse(guild3=="heterotroph-polymer","organic",
                                          ifelse(guild3=="heterotroph-sugar","organic",
                                                 ifelse(guild3=="heterotroph-ch4","organic",
                                                        ifelse(guild3=="methyl_c1","organic",
                                                               ifelse(guild3=="c_fixation","inorganic","central")))))))))%>%
  select(-guild3)%>%
  group_by(MAG,river,code,life)%>%
  summarise(total=sum(sum))%>%
  spread(key=life,value=total)%>%
  mutate_at(c(4,5,6,7,8),~replace_na(.,0))%>%
  mutate(central=ifelse(central>0,"central",""))%>%
  mutate(inorganic=ifelse(inorganic>0,"auto",""))%>%
  mutate(litho=ifelse(litho>0,"litho",""))%>%
  mutate(organic=ifelse(organic>0,"organo",""))%>%
  mutate(photo=ifelse(photo>0,"photo",""))%>%
  mutate(lifestyle=paste(photo,litho,inorganic,organic,central,sep=""))%>%
  mutate(final_lifestyle= ifelse(lifestyle=="organocentral"|lifestyle=="organo","heterotroph",
                                 ifelse(lifestyle=="central","central_only",
                                        ifelse(lifestyle=="autoorganocentral"|
                                                 lifestyle=="lithoorganocentral"|
                                                 lifestyle=="lithoautoorganocentral", "mixotroph",
                                               ifelse(lifestyle=="autocentral"|
                                                        lifestyle=="lithocentral"|
                                                        lifestyle=="litho"|
                                                        lifestyle=="lithoautocentral"|
                                                        lifestyle=="auto","chemolithoautotroph",
                                                      ifelse(lifestyle=="photoautocentral"|
                                                               lifestyle=="photoauto"|
                                                               lifestyle=="photoautoorganocentral","photoautotroph",
                                                             ifelse(lifestyle=="photo"|lifestyle=="photocentral","photo_only","photoheterotroph")))))))

MAG_metaT_guild=MAG_metaT_guild%>%ungroup()%>%select(MAG,river,code,final_lifestyle)%>%left_join(.,c_fix_metaG,by=c("MAG"="genome"))%>%ungroup()%>%
  group_by(MAG,river,code,final_lifestyle)%>%summarise(sum=sum(comp))%>%mutate(c_fix=ifelse(is.na(sum)==1,FALSE,TRUE))%>%select(-sum)%>%
  left_join(.,metaG_adj,by=c("MAG"="genome"))%>%mutate(litho=ifelse(Nitrifier==TRUE | `S Oxidizer`==TRUE,TRUE,FALSE))%>%
  mutate(final_lifestyle2=
           ifelse(final_lifestyle=="photo_only", ifelse(c_fix==TRUE,"photoautotroph","photoheterotroph"),
                  ifelse(final_lifestyle=="central_only" & Photosynthetic==TRUE,ifelse(c_fix==TRUE,"photoautotroph", "photoheterotroph"),
                         ifelse(final_lifestyle=="central_only", ifelse(c_fix==TRUE, ifelse(litho==TRUE,"chemolithoautotroph","heterotroph"), "heterotroph"),final_lifestyle))))


totals=tidy_geTMM %>%separate(X,into=c("MAG","scaff"),remove=FALSE,sep="_Ga0|_k121|_scaff") %>%group_by(MAG,river,code)%>%summarise(n=n(),sum=sum(geTMM))%>%left_join(MAG_metaT_guild,.,by=c("MAG","river","code")) %>% ungroup()%>%group_by(river,code)%>%summarise(total=sum(sum)) 


## pulling carbon functions from heterotrophs
carbon_gene=metaT_adj%>%separate(gene,into=c("MAG","scaff"),remove=FALSE,sep="_Ga0|_k121|_scaff")%>%left_join(.,tidy_geTMM,by=c("gene"="X"))%>%ungroup()%>%select(gene,MAG,river,code,geTMM,final_ID,module,Guild,c)%>%left_join(.,MAG_metaT_guild,by=c("MAG","river","code"))%>%select(-final_lifestyle,-c_fix,-Photosynthetic,-Methanotroph,-Nitrifier,-`S Oxidizer`,-litho)%>%filter(final_lifestyle2!="chemolithoautotroph",final_lifestyle2!="photoautotroph")%>%filter(c=="yes")%>%filter(Guild!="central_c",Guild!="c_fixation")

carbon_gene_totals=carbon_gene%>%ungroup()%>%select(river,code,geTMM)%>%group_by(river,code)%>%summarise(sum=sum(geTMM))%>%left_join(.,totals,by=c("river","code"))%>%left_join(.,river_info,by=c("code"="Sample"))%>%mutate(river_code=paste(StreamUpdated,river,code,sep="_"))


###########
## plotting fig 5A
###########
ave_carbon_gene=carbon_gene%>%left_join(.,carbon_gene_totals,by=c("river","code"))%>%ungroup()%>%group_by(river_code,Guild,sum,StreamUpdated)%>%summarise(guild_sum=sum(geTMM))%>%filter(StreamUpdated!="NA")%>%ungroup()%>%group_by(StreamUpdated,Guild)%>%summarise(ave_abund=mean(guild_sum/sum))
ave_carbon_gene_total=ave_carbon_gene%>%ungroup()%>%group_by(StreamUpdated)%>%summarise(total=sum(ave_abund))  
ave_carbon_gene%>%left_join(.,ave_carbon_gene_total,by="StreamUpdated")%>%group_by(StreamUpdated,Guild)%>%summarise(abund=ave_abund/total)%>%
  ggplot(aes(x=as.numeric(StreamUpdated),
             y=abund,
             stratum=ordered(Guild,levels=c("heterotroph-polymer","heterotroph-sugar","heterotroph-aromatic","heterotroph-ch4","methyl_c1","SCFA  and alcohol conversions","chemolithoautotroph")),
             alluvium=ordered(Guild,levels=c("heterotroph-polymer","heterotroph-sugar","heterotroph-aromatic","heterotroph-ch4","methyl_c1","SCFA  and alcohol conversions","chemolithoautotroph")),
             fill=ordered(Guild,levels=c("heterotroph-polymer","heterotroph-sugar","heterotroph-aromatic","heterotroph-ch4","methyl_c1","SCFA  and alcohol conversions","chemolithoautotroph"))))+
  geom_flow(stat = "alluvium", lode.guidance = "frontback",color = "darkgray") +
  geom_stratum()+
    scale_fill_manual(values=c("#00441b","#238b45","#74c476","#810f7c","#8c6bb1","#8c96c6","#9ebcda"))+
  theme_classic()+theme(legend.title = element_blank())

###########
## plotting fig 5B
###########
carbon_gene%>%left_join(.,carbon_gene_totals,by=c("river","code"))%>%ungroup()%>%group_by(river_code,Guild,sum,StreamUpdated)%>%summarise(guild_sum=sum(geTMM))%>%filter(StreamUpdated!="NA")%>%mutate(group=ifelse(StreamUpdated>=7,"river",ifelse(StreamUpdated>=4,"midsize","headwater")))%>%
  ggplot()+
  geom_boxplot(aes(x=group,y=guild_sum/sum,color=ordered(Guild,levels=c("heterotroph-polymer","heterotroph-sugar","heterotroph-aromatic","heterotroph-ch4","methyl_c1","SCFA  and alcohol conversions","chemolithoautotroph"))))+
  geom_jitter(aes(x=group,y=guild_sum/sum,color=ordered(Guild,levels=c("heterotroph-polymer","heterotroph-sugar","heterotroph-aromatic","heterotroph-ch4","methyl_c1","SCFA  and alcohol conversions","chemolithoautotroph"))))+
  scale_color_manual(values=c("#00441b","#238b45","#74c476","#810f7c","#8c6bb1","#8c96c6","#9ebcda"))+
  facet_wrap(~Guild,scales="free")+
  theme_classic()+theme(legend.title = element_blank())


###########
## stats to analyze differences in carbon use by stream order
###########


# CO ox
lifestyle=carbon_gene%>%left_join(.,carbon_gene_totals,by=c("river","code"))%>%ungroup()%>%group_by(river_code,Guild,sum,StreamUpdated)%>%summarise(guild_sum=sum(geTMM))%>%filter(StreamUpdated!="NA")%>%mutate(group=ifelse(StreamUpdated>=7,"river",ifelse(StreamUpdated>=4,"midsize","headwater")))%>%mutate(abund=guild_sum/sum)%>%filter(Guild=="chemolithoautotroph")
kruskal.test(abund~group,data=lifestyle) # not sig, p=0.3796

# polymer
lifestyle=carbon_gene%>%left_join(.,carbon_gene_totals,by=c("river","code"))%>%ungroup()%>%group_by(river_code,Guild,sum,StreamUpdated)%>%summarise(guild_sum=sum(geTMM))%>%filter(StreamUpdated!="NA")%>%mutate(group=ifelse(StreamUpdated>=7,"river",ifelse(StreamUpdated>=4,"midsize","headwater")))%>%mutate(abund=guild_sum/sum)%>%filter(Guild=="heterotroph-polymer")
kruskal.test(abund~group,data=lifestyle) # p-value = 0.01938
dunnTest(abund ~ group,
         data=lifestyle,
         method="bh") # headwater and river are significantly different, p.adj=0.022

# sugar
lifestyle=carbon_gene%>%left_join(.,carbon_gene_totals,by=c("river","code"))%>%ungroup()%>%group_by(river_code,Guild,sum,StreamUpdated)%>%summarise(guild_sum=sum(geTMM))%>%filter(StreamUpdated!="NA")%>%mutate(group=ifelse(StreamUpdated>=7,"river",ifelse(StreamUpdated>=4,"midsize","headwater")))%>%mutate(abund=guild_sum/sum)%>%filter(Guild=="heterotroph-sugar")
kruskal.test(abund~group,data=lifestyle) # p-value = 0.0005307
dunnTest(abund ~ group,
         data=lifestyle,
         method="bh") # headwater and river (p=0.00056), midsize and river  (p.adj=0.0439) are significantly different

# aromatic
lifestyle=carbon_gene%>%left_join(.,carbon_gene_totals,by=c("river","code"))%>%ungroup()%>%group_by(river_code,Guild,sum,StreamUpdated)%>%summarise(guild_sum=sum(geTMM))%>%filter(StreamUpdated!="NA")%>%mutate(group=ifelse(StreamUpdated>=7,"river",ifelse(StreamUpdated>=4,"midsize","headwater")))%>%mutate(abund=guild_sum/sum)%>%filter(Guild=="heterotroph-aromatic")
kruskal.test(abund~group,data=lifestyle) # p-value = 0.03062
dunnTest(abund ~ group,
         data=lifestyle,
         method="bh") # headwater and river (p.adj=0.02752059), midsize and river (p.adj=0.04589838) are significantly different

# methanotroph
lifestyle=carbon_gene%>%left_join(.,carbon_gene_totals,by=c("river","code"))%>%ungroup()%>%group_by(river_code,Guild,sum,StreamUpdated)%>%summarise(guild_sum=sum(geTMM))%>%filter(StreamUpdated!="NA")%>%mutate(group=ifelse(StreamUpdated>=7,"river",ifelse(StreamUpdated>=4,"midsize","headwater")))%>%mutate(abund=guild_sum/sum)%>%filter(Guild=="heterotroph-ch4")
kruskal.test(abund~group,data=lifestyle) # p-value = 0.002032
dunnTest(abund ~ group,
         data=lifestyle,
         method="bh") # headwater and river are significantly different (p.adj=0.002685785)

# methyl-c1
lifestyle=carbon_gene%>%left_join(.,carbon_gene_totals,by=c("river","code"))%>%ungroup()%>%group_by(river_code,Guild,sum,StreamUpdated)%>%summarise(guild_sum=sum(geTMM))%>%filter(StreamUpdated!="NA")%>%mutate(group=ifelse(StreamUpdated>=7,"river",ifelse(StreamUpdated>=4,"midsize","headwater")))%>%mutate(abund=guild_sum/sum)%>%filter(Guild=="methyl_c1")
kruskal.test(abund~group,data=lifestyle) # p-value = 0.0005856
dunnTest(abund ~ group,
         data=lifestyle,
         method="bh") # headwater and river (p.adj=0.001114524), headwater and midsize (p.adj=0.014989628) are significantly different

# SCFA
lifestyle=carbon_gene%>%left_join(.,carbon_gene_totals,by=c("river","code"))%>%ungroup()%>%group_by(river_code,Guild,sum,StreamUpdated)%>%summarise(guild_sum=sum(geTMM))%>%filter(StreamUpdated!="NA")%>%mutate(group=ifelse(StreamUpdated>=7,"river",ifelse(StreamUpdated>=4,"midsize","headwater")))%>%mutate(abund=guild_sum/sum)%>%filter(Guild=="SCFA  and alcohol conversions")
kruskal.test(abund~group,data=lifestyle) # p-value = 0.002426
dunnTest(abund ~ group,
         data=lifestyle,
         method="bh") # headwater and river are significantly different (p.adj=0.001703097)

