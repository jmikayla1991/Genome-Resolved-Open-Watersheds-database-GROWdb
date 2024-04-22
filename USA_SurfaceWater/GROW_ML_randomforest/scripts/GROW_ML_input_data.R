## ---------------------------
##
## Script: GROW_ML_input_data
##
## Purpose: make the dataframe for input to ML models
##
## Author: Ikaia Leleiwi
##
## Date Created: May 23rd 2023
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
#taxonomy
tax <- read_xlsx("data/SOMfiles/GROW_SuppFile2.xlsx") 

#metag
metag_115_vars <- read_csv(paste0("data/metag/",
                                  "GROWdb_with_vars_20220715_115_metaG.csv"))


#metat
metat_52 <- read_csv("data/metat/corrplot_mergedtable_52.csv") %>%
  rename("Sample" = "...1")

bin_counts <- read_csv("data/metat/geTMM_norm.counts.rpk_edger_genome_relabund_52.csv")



#combined metag metat
metat_metag <- metat_52 %>%
  left_join(metag_115_vars, by = "Sample") %>%
  drop_na() 

write_csv(metat_metag, "data/metat_metag.csv")

#make target variable of diff function classes
metat_metag_target_function <- metat_52 %>%
  left_join(metag_115_vars, by = "Sample") %>%
  drop_na() %>%
  mutate_at(c("Aerobic", "Microaerophillic", "Light_driven", 
              "N_Reducer", "DNRA", "Methanotroph"),
            scale) %>% #scale function class counts across samples
  pivot_longer(cols = c("Aerobic", "Microaerophillic", "Light_driven", 
                        "N_Reducer", "DNRA", "Methanotroph"),
               names_to = "function_class",
               values_to = "funct_cls_gene_cts") %>%
  mutate(funct_cls_gene_cts = as.double(funct_cls_gene_cts)) %>%
  select(function_class, 
         funct_cls_gene_cts, starts_with("Pct"), 
         DamDensWs, PopDen2010Ws) #keep only target column, gene counts column, and relevant predictor columns

#visualize precictor columns
plot_hist <- function(var, name){
  
  ggplot(metat_metag_target_function, aes(x = (var))) +
    geom_histogram(bins = 100) +
    labs(x = name,
         title = name)
  
}

hist_list <- metat_metag_target_function %>%
  select(-funct_cls_gene_cts) %>%
  select_if(is.numeric) %>%
  map2(.y = names(.),
       ~plot_hist(.x, .y)) 

map(hist_list, ~print(.x))


#add pseudo count to predictor columns
#log transform
#scaled each predictor variable by the scaled function class counts
#drop gene counts variable from df
fct_landuse_scaled <- metat_metag_target_function %>%
  mutate(across(.cols = !function_class & !funct_cls_gene_cts,~.+1e9),
         across(.cols = !function_class & !funct_cls_gene_cts,log10),
         across(.cols = c(starts_with("Pct"), DamDensWs, PopDen2010Ws), 
                ~ ./funct_cls_gene_cts)) %>%
  select(-funct_cls_gene_cts)

write_csv(fct_landuse_scaled, "data/ml_input_scaled.csv")

#make ml input with taxa as target

clip_d <- function(x){
  return(str_remove(x,"^.__"))
}

taxa_landuse <- bin_counts %>%
  pivot_longer(cols = -Sample,
               names_to = "bin",
               values_to = "relative_abundance") %>%
  left_join(metag_115_vars, by = "Sample") %>%
  drop_na() %>%
  filter(relative_abundance > 0) %>%
  group_by(bin) %>%
  summarise_if(is.numeric, mean) %>%
  select(bin, starts_with("Pct"), DamDensWs, PopDen2010Ws) %>%
  left_join(tax, by = c("bin" = "user_genome")) %>%
  mutate_at(.vars = c("Domain", "Phylum", "Class", "Order", 
                      "Family", "Genus", "Species"),
            .funs = clip_d) %>%
  mutate(lowest_tax = case_when(Genus == "" ~ paste(Family, "(family)"),
                                Species == "" ~ paste(Genus, "(genus)"),
                                T ~ paste(Species, "(species)")))

write_csv(taxa_landuse, "data/ml_input_taxa.csv")


