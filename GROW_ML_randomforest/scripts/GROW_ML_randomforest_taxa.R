#!/usr/bin/env Rscript
## ---------------------------
##
## Script: GROW_ML_randomforest_taxa
##
## Purpose: predict spceific lineages with landuse
##
## Author: Ikaia Leleiwi
##
## Date Created: 
##
## Copyright (c) Ikaia Leleiwi, 2023
## Email: ileleiwi@gmail.com
##
## ---------------------------
##
## Notes: Script takes overnight to complete, best to run on server
##   
##
## ---------------------------

## set working directory

setwd(paste0("GROW_ML_randomforest"))

## ---------------------------

##Libraries

library(tidyverse)
library(caret)
library(pROC)
library(ranger) #importance_pvalues
#set.seed(123)

##Data

target_col <- c("bin", "lowest_tax", "Domain", "Phylum", "Class", "Order", 
                "Family", "Genus", "Species")
df <- read_csv("data/ml_input_taxa.csv") %>%
  select(target_col[5], starts_with("Pct"), DamDensWs, PopDen2010Ws) %>%
  rename("target" = target_col[5])

#find classes with very few observations
sparce_target <- df %>%
  group_by(target) %>%
  summarise(n = n()) %>%
  filter(n <= 3) %>%  
  pull(target)

#check for all 0's or sparce (only 1)
df_clean <- df %>%
  filter(!(target %in% sparce_target)) %>% #remove low observation classes from prediction dataset
  mutate_if(is.numeric,function(x) x + 1e-9) %>% #add pseudocount
  mutate_if(is.numeric,log10) #log transform


#look at data
plot_hist <- function(var, name){
  
  ggplot(df_clean, aes(x = (var))) +
    geom_histogram(bins = 100) +
    labs(x = name,
         title = name)
  
}

hist_list <- df_clean %>%
  select_if(is.numeric) %>%
  map2(.y = names(.),
       ~plot_hist(.x, .y)) 

map(hist_list, ~print(.x))

#####################bootstrap models############################
model_list <- list()
cm_list <- list()
varimp_list <- list()

for(s in seq(1:1000)){
  #split data
  set.seed(s)
  training_partition <- sample(nrow(df_clean), size = nrow(df_clean)*0.75)
  
  train <- df_clean[training_partition,]
  test <- df_clean[-training_partition,]
  
  #resampling method
  fit_control <- trainControl(## 10-fold CV
    method = "cv",
    number = 10,
    search = 'grid', 
    savePredictions = TRUE)
  
  #tuning paramaters
  rang_grid <- expand.grid(mtry = seq(2,8,2),
                           splitrule = c("gini", "extratrees"),
                           min.node.size = c(1, 3, 5, 6))
  
  rf_grid <- expand.grid(.mtry = c(sqrt(ncol(df_clean))))
  
  #fit a random forest model (using ranger)
  rang_fit <- train(as.factor(target) ~ .,
                    data = train,
                    method = "ranger",
                    trControl = fit_control,
                    tuneGrid = rang_grid,
                    importance = 'impurity',
                    metric = "Kappa")
  
  model_list[[length(model_list)+1]] <- rang_fit
  
  
  # predict the outcome on a test set
  rang_pred <- predict(rang_fit, test)
  
  # compare predicted outcome and true outcome
  cm <- confusionMatrix(rang_pred, factor(test$target, levels = levels(rang_pred)))
  cm_list[[length(cm_list)+1]] <- cm
 
  #variable importances
  #Conditional=True, adjusts for correlations between predictors.
  i_scores <- varImp(rang_fit, conditional=TRUE, scale = FALSE)
  
  #Gathering rownames in 'var'  and converting it to the factor
  i_scores_p <- i_scores$importance %>% 
    rownames_to_column("var") 
  
  scores_levels <- i_scores_p %>%
    arrange(Overall) %>%
    pull(var) %>%
    unique()
  
  i_scores_p <- i_scores_p %>%
    mutate(var = factor(var, levels = scores_levels))
  varimp_list[[length(varimp_list)+1]] <- i_scores_p
}

#collect all stats and variables from all the models and write out in dfs
overall_stats_list <- map(.x = seq(1:1000), .f = function(.x){
  
  out_df <- as.data.frame(cm_list[[.x]]["overall"]) %>%
    rownames_to_column(var = "stat") %>%
    mutate(model_name = paste0("model_", !!.x))
  return(out_df)
}) 
overall_stats_df <- do.call(rbind, overall_stats_list) %>%
  pivot_wider(names_from = "stat",
              values_from = "overall",
              values_fill = NaN)
write_csv(overall_stats_df, "rf_taxa_overall_stats_df.csv")

byClass_list <- map(.x = seq(1:1000), .f = function(.x){
  
  out_df <- as.data.frame(cm_list[[.x]]["byClass"]) %>%
    rownames_to_column(var = "Class") %>%
    mutate(model_name = paste0("model_", !!.x))
  return(out_df)
}) 
byClass_df <- do.call(rbind, byClass_list) 
write_csv(byClass_df, "rf_taxa_byClass_df.csv")


varimp_df <- map_dfr(.x = seq(1:1000), .f = function(.x){
  
  out_df <- as.data.frame(varimp_list[[.x]]) %>%
    mutate(model_name = paste0("model_", !!.x))
  return(out_df)
}) 
byClass_df <- do.call(rbind, byClass_list) 

write_csv(varimp_df, "rf_taxa_varimp_df.csv")

# #Plotting the bar and polar charts for comparing variables
# i_scores_p %>%
#   ggplot(aes(x = var, y = Overall)) +
#   geom_bar(stat = "identity") +
#   coord_flip()

