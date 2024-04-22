## ---------------------------
##
## Script: GRWO_ML_randomforest_landusereg
##
## Purpose: train random forest regression model to predict landuse by individual functions
##
## Author: Ikaia Leleiwi
##
## Date Created: May 26th 2023
##
## Copyright (c) Ikaia Leleiwi, 2023
## Email: ileleiwi@gmail.com
##
## ---------------------------
##
## Notes:
#       Random Forest Regression R^2 Results
#       
#       
#       
#       
#       
#       
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

pred_cols <- c("Aerobic", "Microaerophillic", "Light_driven", "N_Reducer",
                "DNRA", "Methanotroph")

target_col <- c("DamDensWs", "PctFire2010Ws", "PctOw2016Ws", "PctIce2016Ws", 
                "PctUrbOp2016Ws", "PctUrbLo2016Ws", "PctUrbMd2016Ws", 
                "PctUrbHi2016Ws", "PctBl2016Ws", "PctDecid2016Ws", 
                "PctConif2016Ws", "PctMxFst2016Ws", "PctShrb2016Ws", 
                "PctGrs2016Ws", "PctHay2016Ws", "PctCrop2016Ws", 
                "PctWdWet2016Ws", "PctHbWet2016Ws", "PopDen2010Ws")

# df <- read_csv("data/metat_metag.csv") 
#   select(target_col[7], starts_with("Pct"), DamDensWs, PopDen2010Ws) %>%
#   rename("target" = target_col[7]) 

df <- read_csv("data/metat_metag.csv") %>%
  mutate(row_max = pmap_chr(across(starts_with("Pct")), ~ names(c(...)[which.max(c(...))])))  #gets column name of highest landuse cat. in each row

#check sparcity
df %>% 
  select(all_of(target_col)) %>%
  imap_dfc(~if(is.numeric(.x)){ifelse(.x > 0, 1, 0)} else(.x)) %>%
  colSums()

#drop PctFire2010Ws and PctIce2016Ws only 3 and 4 observations respectively
df1 <- df %>%
  select(-PctFire2010Ws, -PctIce2016Ws)

target_col <- c("DamDensWs", "PctOw2016Ws",  
                "PctUrbOp2016Ws", "PctUrbLo2016Ws", "PctUrbMd2016Ws", 
                "PctUrbHi2016Ws", "PctBl2016Ws", "PctDecid2016Ws", 
                "PctConif2016Ws", "PctMxFst2016Ws", "PctShrb2016Ws", 
                "PctGrs2016Ws", "PctHay2016Ws", "PctCrop2016Ws", 
                "PctWdWet2016Ws", "PctHbWet2016Ws", "PopDen2010Ws")

r2 <- list()
for(i in 1:length(target_col)){
  set.seed(123)
  df2 <- df1 %>%
  select(target_col[i], all_of(pred_cols)) %>%
  rename("target" = target_col[i])

#check for all 0's or sparce (only 1)
df_clean <- df2 %>%
  select(where(~ any(. !=0 ))) %>% #remove columns with all 0's and all FALSE
  select(where(~ any(. !=1 ))) %>%
  mutate_at(pred_cols, function(x) x + 1e-9) %>% #add pseudocount
  mutate_at(pred_cols, log10) #log transform


  
# #look at data
# plot_hist <- function(var, name){
# 
#   ggplot(df_clean, aes(x = (var))) +
#     geom_histogram(bins = 100) +
#     labs(x = name,
#          title = name)
# 
# }
# 
# hist_list <- df_clean %>%
#   select_if(is.numeric) %>%
#   map2(.y = names(.),
#        ~plot_hist(.x, .y))
# 
# map(hist_list, ~print(.x))

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
rang_grid <- expand.grid(mtry = seq(1,6,1),
                         splitrule = c("variance", "extratrees", "maxstat"),
                         min.node.size = c(1, 3, 5, 6))

rf_grid <- expand.grid(.mtry = c(sqrt(ncol(df_clean))))

#fit a random forest model (using ranger)
rang_fit <- train(target ~ .,
                  data = train,
                  method = "ranger",
                  trControl = fit_control,
                  tuneGrid = rang_grid,
                  importance = 'impurity',
                  metric = "Rsquared")

rang_fit


# predict the outcome on a test set
rang_pred <- predict(rang_fit, test)

plot(test$target ~ rang_pred)

#metrics
d <- test$target - rang_pred
mse = mean((d)^2)
mae = mean(abs(d))
rmse = sqrt(mse)
R2 = 1-(sum((d)^2)/sum((test$target-mean(test$target))^2))
cat(" MAE:", mae, "\n", "MSE:", mse, "\n", 
    "RMSE:", rmse, "\n", "R-squared:", R2)

r2[[length(r2)+1]] <- R2

}
names(r2) <- target_col
r2

######## test top landuse type as categorical target ########################################################################
df2 <- df1 %>%
  select(row_max, all_of(pred_cols)) %>%
  rename("target" = "row_max")

#check for all 0's or sparce (only 1)
df_clean <- df2 %>%
  select(where(~ any(. !=0 ))) %>% #remove columns with all 0's and all FALSE
  select(where(~ any(. !=1 ))) %>%
  mutate_at(pred_cols, function(x) x + 1e-9) %>% #add pseudocount
  mutate_at(pred_cols, log10) #log transform



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
rang_grid <- expand.grid(mtry = seq(1,6,1),
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
confusionMatrix(rang_pred, factor(test$target, levels = levels(rang_pred)))


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