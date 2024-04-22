## ---------------------------
##
## Script: GROW_ML_explore_data
##
## Purpose: familarize myself with the Grow data and devise ml approach
##
## Author: Ikaia Leleiwi
##
## Date Created: Friday May 19th 2023
##
## Copyright (c) Ikaia Leleiwi, 2023
## Email: ileleiwi@gmail.com
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

## set working directory

setwd(paste0("GROW_ML_randomforest"))

## ---------------------------

##Libraries

library(tidyverse)
library(readxl)
library(vegan)
library(mapview)

##Data

#metadata
tax <- read_xlsx("data/SOMfiles/GROW_SuppFile2.xlsx")
som <- tax <- read_xlsx("data/SOMfiles/GROW_SuppFile1.xlsx")

#metag
metag_115_ra <- read_csv(paste0("data/metag/",
                                "strict_mapping_table_95id_3x_60%cov_relabundance_115.csv"))

metag_115_vars <- read_csv(paste0("data/metag/",
                                "GROWdb_with_vars_20220715_115_metaG.csv"))

#metat
metat_52 <- read_csv("data/metat/corrplot_mergedtable_52.csv") %>%
  rename("Sample" = "...1")


#combined metag metat
metat_metag <- metat_52 %>%
  left_join(metag_115_vars, by = "Sample") %>%
  drop_na() 

write_csv(metat_metag, "data/metat_metag.csv")

metat_metag_target_function <- metat_52 %>%
  left_join(metag_115_vars, by = "Sample") %>%
  drop_na() %>%
  pivot_longer(cols = c("Aerobic", "Microaerophillic", "Light_driven", 
                        "N_Reducer", "DNRA", "Methanotroph"),
               names_to = "function_class",
               values_to = "funct_cls_gene_cts") %>%
  select(function_class, funct_cls_gene_cts, everything())

metat_metag_target_landuse <- metat_52 %>%
  left_join(metag_115_vars, by = "Sample") %>%
  drop_na() %>%
  pivot_longer(cols = c(starts_with("Pct"), "PopDen2010Ws"),
               names_to = "landuse_class",
               values_to = "landuse_cts") %>%
  select(landuse_class, landuse_cts, everything())

pct <- metat_metag %>%
  select(Sample, starts_with("Pct")) %>%
  mutate(total = rowSums(across(where(is.numeric))))

fct_landuse <- metat_metag_target_function %>%
  select(function_class, funct_cls_gene_cts, starts_with("Pct")) %>%
  mutate(funct_cls_gene_cts = funct_cls_gene_cts + 1e-9,
         across(.cols = starts_with("Pct"), ~ ./funct_cls_gene_cts))



#functions
richness <- function(x){
  
  sum(x>0) 
  
}

shannon <- function(x){
  
  rabund <- x[x>0]/sum(x)
  -sum(rabund * log(rabund))
  
}


#clean metag relative abundance table
ra <- metag_115_ra %>%
  column_to_rownames(var = "Sample") %>%
  t() %>%
  as.data.frame() %>%
  mutate(total = rowSums(across(where(is.numeric)))) %>%
  filter(total > 0) #remove bins with 0 relative abundance in all samples

ra_pa <- ra %>%
  rownames_to_column(var = "bins") %>%
  imap_dfc(~if(is.numeric(.x)){ifelse(.x > 0, 1, 0)} else(.x))  #convert rel abund to present absent

alpha_div <- ra_pa %>%
  select(-total) %>%
  pivot_longer(cols = -bins,
               values_to = "pa",
               names_to = "sample") %>%
  group_by(sample) %>%
  summarise(richness = richness(pa),
            shannon = shannon(pa))

metag_ra_alpha <- ra %>%
  select(-total) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample") %>%
  left_join(alpha_div, by = "sample")

#plot metat samples on map
mapview(metat_52, xcol = "Long", ycol = "Lat", crs = 4269, grid = FALSE)


