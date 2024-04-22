## ---------------------------
##
## Script: GRWO_ML_randomforest
##
## Purpose: train random forest model to predict function categories from landuse
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
library(caret)
library(pROC)
library(ranger) #importance_pvalues
set.seed(123)

##Data

df <- read_csv("data/ml_input_scaled.csv")

#check for all 0's or sparce (only 1)
df_clean <- df %>%
  select(where(~ any(. !=0 ))) %>% #remove columns with all 0's and all FALSE
  select(where(~ any(. !=1 ))) %>%
  rename(target = function_class) 

#find fct classes with very few observations
sparce_target <- df_clean %>%
  mutate(total = rowSums(across(where(is.numeric))))


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

#split data
training_partition <- df_clean$target %>%
  createDataPartition(p=0.75, list = FALSE)

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

rang_fit


# predict the outcome on a test set
rang_pred <- predict(rang_fit, test)

# compare predicted outcome and true outcome
confusionMatrix(rang_pred, as.factor(test$target))


#roc
rang_roc_obj <- multiclass.roc(test$target, as.numeric(rang_pred))
auc(rang_roc_obj)

rs <- rang_roc_obj[['rocs']]
plot.roc(rs[[1]])
sapply(2:length(rs),function(i) lines.roc(rs[[i]],col=i))


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
