####################
## Description:   (this was run with bartMachine v1.3.2)
##  - In this file we:
##    - Fit a BART PS model to all variables.
##    - Perform variable selection on PS model and refit model.
##    - Refit BART PS model
#################### 


## increase memory usage to 100GB of RAM (needs to be run before library(bartMachine))
options(java.parameters = "-Xmx100g")

## Load libraries
library(tidyverse)
library(bartMachine)

## Set up directory path to save files (stagered to ensure folders are created)

dir.create("Samples")
dir.create("Samples/SGLT2-GLP1")

output_path <- "Samples/SGLT2-GLP1/Aurum"
dir.create(output_path)

dir.create(paste0(output_path, "/ps_model"))

dir.create("Plots")


###############################################################################
###############################################################################
############################### Read Data In ##################################
###############################################################################
###############################################################################

## Load functions required

source("01.slade_aurum_functions.R")
source("02.slade_aurum_set_data.R")

## Load dataset
ps.model.train <- set_up_data_sglt2_glp1(dataset.type = "ps.model.train")

###############################################################################
###############################################################################
################################ FUNCTIONS ####################################
###############################################################################
###############################################################################


########
### Fit a propensity score model to all variables in dataset
########


# Check the optimal number of trees (does not work for classification)


if (class(try(
  
  bart_ps_model <- readRDS(paste0(output_path, "/ps_model/bart_ps_model.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  bart_ps_model <- bartMachine::bartMachine(X = ps.model.train %>%
                                              select(-patid,
                                                     -pated,
                                                     -drugclass),
                                            y = ps.model.train[,"drugclass"] %>%
                                              unlist(),
                                            num_trees = 50,
                                            use_missing_data = TRUE,
                                            num_burn_in = 15000,
                                            num_iterations_after_burn_in = 10000,
                                            serialize = TRUE)
  
  saveRDS(bart_ps_model, paste0(output_path, "/ps_model/bart_ps_model.rds"))
  
}


########
### Variable selection of propensity score model
########

if (class(try(
  
  vs_bart_ps_model <- readRDS(paste0(output_path, "/ps_model/vs_bart_ps_model.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  pdf(file = "Plots/11.04.prop_model_vs.pdf", width = 18, height = 11)
  # error with cv
  vs_bart_ps_model <- var_selection_by_permute(bart_ps_model)
  dev.off()
  
  ## Variables selected
  # [1] "prebmi"      "yrdrugstart" "preegfr"     "prehba1c"    "drugline"
  # [6] "ncurrtx"     "ethnicity"
  
  
  saveRDS(vs_bart_ps_model, paste0(output_path, "/ps_model/vs_bart_ps_model.rds"))
  
}

variables_chosen <- unique(gsub("_.*", "", vs_bart_ps_model$important_vars_local_names))

########
### Refit PS model with selected vars
########

## Fit initial model using all the available variables to estimate HbA1c outcome
if (class(try(
  
  bart_ps_model_final <- readRDS(paste0(output_path, "/ps_model/bart_ps_model_final.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  bart_ps_model_final <- bartMachine::bartMachine(X = ps.model.train %>%
                                                    select(
                                                      all_of(variables_chosen)
                                                    ),
                                                  y = ps.model.train[,"drugclass"] %>%
                                                    unlist(),
                                                  num_trees = 50,
                                                  use_missing_data = TRUE,
                                                  num_burn_in = 15000,
                                                  num_iterations_after_burn_in = 10000,
                                                  serialize = TRUE)
  
  saveRDS(bart_ps_model_final, paste0(output_path, "/ps_model/bart_ps_model_final.rds"))
  
}


## Checking ROC values in development dataset


# predictions <- bart_ps_model_final$p_hat_train
# 
# df <- as.data.frame(cbind(predictions, ps.model.train["drugclass"] %>%
#                             mutate(drugclass = ifelse(drugclass == "SGLT2", 0, 1)))) %>%
#   set_names(c("predictions", "labels"))
# 
# library(ROCR)
# pred <- prediction(df$predictions, df$labels)
# perf <- performance(pred,"tpr","fpr")
# pdf()
# plot(perf,colorize=TRUE)
# dev.off()
# 
# library(pROC)
# pROC_obj <- roc(df$labels,df$predictions,
#                 smoothed = TRUE,
#                 # arguments for ci
#                 ci=TRUE, ci.alpha=0.9, stratified=FALSE,
#                 # arguments for plot
#                 plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
#                 print.auc=TRUE, show.thres=TRUE)
# 
# 
# 
# values <- coords(pROC_obj, "best", ret=c("threshold", "specificity", "sensitivity", "accuracy",
#                                          "precision", "recall"), transpose = FALSE)
# 
# 
# confusion.matrix.dev <- as.data.frame(cbind(predictions, ps.model.train[,"drugclass"])) %>%
#   set_names(c("predictions", "labels")) %>%
#   mutate(predictions = ifelse(predictions > values$threshold, 1, 2)) %>% table()



###############################################################################
###############################################################################
###################### Propensity scores for patients #########################
###############################################################################
###############################################################################


# Prop scores for train dataset
patient_prop_scores <- ps.model.train %>%
  select(patid, pated) %>%
  cbind(prop.score = bart_ps_model_final$p_hat_train)


# Prop scores for test dataset

ps.model.test <- set_up_data_sglt2_glp1(dataset.type = "ps.model.test")

# calculate prop score
if (class(try(
  
  prop_score_testing_data <- readRDS(paste0(output_path, "/ps_model/prop_score_testing_data.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  prop_score_testing_data <- predict(bart_ps_model_final, ps.model.test %>%
                                       select(
                                         colnames(bart_ps_model_final$X)
                                       ))
  
  saveRDS(prop_score_testing_data, paste0(output_path, "/ps_model/prop_score_testing_data.rds"))
  
}



## Checking ROC values in validation dataset


# predictions <- prop_score_testing_data
# 
# df <- as.data.frame(cbind(predictions, ps.model.test["drugclass"] %>%
#                             mutate(drugclass = ifelse(drugclass == "SGLT2", 0, 1)))) %>%
#   set_names(c("predictions", "labels"))
# 
# library(ROCR)
# pred <- prediction(df$predictions, df$labels)
# perf <- performance(pred,"tpr","fpr")
# pdf()
# plot(perf,colorize=TRUE)
# dev.off()
# 
# library(pROC)
# pROC_obj <- roc(df$labels,df$predictions,
#                 smoothed = TRUE,
#                 # arguments for ci
#                 ci=TRUE, ci.alpha=0.9, stratified=FALSE,
#                 # arguments for plot
#                 plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
#                 print.auc=TRUE, show.thres=TRUE)
# 
# 
# 
# # values <- coords(pROC_obj, "best", ret=c("threshold", "specificity", "sensitivity", "accuracy",
# #                                          "precision", "recall"), transpose = FALSE)
# 
# 
# confusion.matrix.val <- as.data.frame(cbind(predictions, ps.model.test[,"drugclass"])) %>%
#   set_names(c("predictions", "labels")) %>%
#   mutate(predictions = ifelse(predictions > values$threshold, 1, 2)) %>% table()




patient_prop_scores <- patient_prop_scores %>%
  rbind(
    ps.model.test %>%
      select(patid, pated) %>%
      cbind(prop.score = prop_score_testing_data)
  ) %>%
  as.data.frame()

saveRDS(patient_prop_scores, paste0(output_path, "/ps_model/patient_prop_scores.rds"))





#:--------------------------------------------------------------------------------------------------
#:--------------------------------------------------------------------------------------------------
#:--------------------------------------------------------------------------------------------------

# ############### Redoing this in the newest version of the bartMachine for American dataset validation - bartMachine v1.3.3.1
# # 
# 
# 
# ## increase memory usage to 100GB of RAM (needs to be run before library(bartMachine))
# options(java.parameters = "-Xmx100g")
# 
# ## Load libraries
# library(tidyverse)
# library(bartMachine)
# 
# ## Set up directory path to save files (stagered to ensure folders are created)
# 
# dir.create("Samples")
# dir.create("Samples/SGLT2-GLP1")
# 
# output_path <- "Samples/SGLT2-GLP1/Aurum"
# dir.create(output_path)
# 
# dir.create(paste0(output_path, "/ps_model"))
# 
# dir.create("Plots")
# 
# ## Load functions required
# 
# source("01.slade_aurum_functions.R")
# source("02.slade_aurum_set_data.R")
# 
# ## Load dataset
# ps.model.train <- set_up_data_sglt2_glp1(dataset.type = "ps.model.train")
# 
# 
# ########
# ### Fit a propensity score model to all variables selected from the dataset
# ########
# 
# if (class(try(
#   
#   vs_bart_ps_model <- readRDS(paste0(output_path, "/ps_model/vs_bart_ps_model.rds"))
#   
#   , silent = TRUE)) == "try-error") {
#   
#   pdf(file = "Plots/11.04.prop_model_vs.pdf", width = 18, height = 11)
#   # error with cv
#   vs_bart_ps_model <- var_selection_by_permute(bart_ps_model)
#   dev.off()
#   
#   ## Variables selected
#   # [1] "prebmi"      "yrdrugstart" "preegfr"     "prehba1c"    "drugline"
#   # [6] "ncurrtx"     "ethnicity"
#   
#   
#   saveRDS(vs_bart_ps_model, paste0(output_path, "/ps_model/vs_bart_ps_model.rds"))
#   
# }
# 
# variables_chosen <- unique(gsub("_.*", "", vs_bart_ps_model$important_vars_local_names))
# 
# 
# 
# ## Fit initial model using all the available variables to estimate HbA1c outcome
# bart_ps_model_final <- bartMachine::bartMachine(X = ps.model.train %>%
#                                                   select(
#                                                     all_of(variables_chosen)
#                                                   ),
#                                                 y = ps.model.train[,"drugclass"] %>%
#                                                   unlist(),
#                                                 num_trees = 50,
#                                                 use_missing_data = TRUE,
#                                                 num_burn_in = 15000,
#                                                 num_iterations_after_burn_in = 10000,
#                                                 serialize = TRUE)
# 
# 
# bart_ps_model_test_delect <- readRDS("bart_ps_model_test_delect.rds")
# 
# bart_ps_model_final$X <- bart_ps_model_test_delect$X
# bart_ps_model_final$y <- bart_ps_model_test_delect$y
# 
# 
# saveRDS(bart_ps_model_final, "bart_ps_american_model_version.rds")
# 



