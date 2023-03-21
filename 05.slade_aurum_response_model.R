####################
## Description:    (this file uses BCF v.2.0.1)
##  - In this file we:
##    - Fit a SparseBCF model to choose the variables to include as control and moderators.
##    - Fit a BCF model with the specific chosen variables used in control and in moderators.
##    - Validate model
####################


## Load libraries
library(tidyverse)
# library(SparseBCF)    # This package can only be loaded when BCF isn't loaded. They conflict with each other.
require(bcf)
library(cobalt)

## Set up directory path to save files (stagered to ensure folders are created)

dir.create("Samples")
dir.create("Samples/SGLT2-GLP1")

output_path <- "Samples/SGLT2-GLP1/Aurum"
dir.create(output_path)

dir.create(paste0(output_path, "/response_model_bcf"))

dir.create(paste0(output_path, "/response_model_bcf/trees_prop"))

dir.create(paste0(output_path, "/response_model_bcf/trees_no_prop"))

dir.create(paste0(output_path, "/response_model_bcf/assessment"))

dir.create("Plots")


###############################################################################
###############################################################################
############################### Read Data In ##################################
###############################################################################
###############################################################################

## Load functions required

source("11.01.slade_aurum_functions.R")
source("11.02.slade_aurum_set_data.R")

## Load dataset
hba1c.train <- set_up_data_sglt2_glp1(dataset.type = "hba1c.train")

## Collect propensity score values
patient_prop_scores <- readRDS(paste0(output_path, "/ps_model/patient_prop_scores.rds"))

## Append values
hba1c.train <- hba1c.train %>%
  left_join(patient_prop_scores, by = c("patid", "pated"))

## Prepare dataset to be used in SparseBCF
hba1c.train.complete <- hba1c.train %>%
  # drop the variables with the most missingness (>40%)
  select(-preacr, -preast, -prehaematocrit, -prehaemoglobin, -pretriglyceride) %>%
  # only complete cases
  drop_na() %>%
  as.data.frame()


#:----------------------------------------------------------------------------------------------
# Fit 2 chains of SparseBCF

if (class(try(
  
  sparsebcf_chain_1 <- readRDS(paste0(output_path, "/response_model_bcf/sparsebcf_chain_1.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  
  sparsebcf_chain_1 <- SparseBCF(y = hba1c.train.complete$posthba1cfinal,
                                 z = hba1c.train.complete %>%
                                   select(drugclass) %>%
                                   mutate(drugclass = ifelse(drugclass == "GLP1", 0, 1)) %>%
                                   unlist(),
                                 x_control = hba1c.train.complete %>%
                                   select(
                                     -patid,
                                     -pated,
                                     -drugclass,
                                     -posthba1cfinal,
                                     -prop.score) %>%
                                   mutate_all(funs(as.numeric(.))) %>%
                                   as.matrix(),
                                 pihat = 1-hba1c.train.complete$prop.score,
                                 OOB = F,
                                 sparse = T,
                                 update_interval = 250,
                                 ntree_control = 200,
                                 sd_control = 2 * sd(hba1c.train.complete$posthba1cfinal),
                                 base_control = 0.95,
                                 power_control = 2,
                                 ntree_moderate = 200,
                                 sd_moderate = 2 * sd(hba1c.train.complete$posthba1cfinal),
                                 base_moderate = 0.95,
                                 power_moderate = 2,
                                 nburn = 50000,
                                 nsim = 200000,
                                 use_muscale = FALSE,
                                 use_tauscale = FALSE
  )
  
  
  saveRDS(sparsebcf_chain_1, paste0(output_path, "/response_model_bcf/sparsebcf_chain_1.rds"))
  
}

if (class(try(
  
  sparsebcf_chain_2 <- readRDS(paste0(output_path, "/response_model_bcf/sparsebcf_chain_2.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  
  sparsebcf_chain_2 <- SparseBCF(y = hba1c.train.complete$posthba1cfinal,
                                 z = hba1c.train.complete %>%
                                   select(drugclass) %>%
                                   mutate(drugclass = ifelse(drugclass == "GLP1", 0, 1)) %>%
                                   unlist(),
                                 x_control = hba1c.train.complete %>%
                                   select(
                                     -patid,
                                     -pated,
                                     -drugclass,
                                     -posthba1cfinal,
                                     -prop.score) %>%
                                   mutate_all(funs(as.numeric(.))) %>%
                                   as.matrix(),
                                 pihat = 1-hba1c.train.complete$prop.score,
                                 OOB = F,
                                 sparse = T,
                                 update_interval = 250,
                                 ntree_control = 200,
                                 sd_control = 2 * sd(hba1c.train.complete$posthba1cfinal),
                                 base_control = 0.95,
                                 power_control = 2,
                                 ntree_moderate = 200,
                                 sd_moderate = 2 * sd(hba1c.train.complete$posthba1cfinal),
                                 base_moderate = 0.95,
                                 power_moderate = 2,
                                 nburn = 50000,
                                 nsim = 200000,
                                 use_muscale = FALSE,
                                 use_tauscale = FALSE
  )
  
  
  saveRDS(sparsebcf_chain_2, paste0(output_path, "/response_model_bcf/sparsebcf_chain_2.rds"))
  
}


## Check trace plots for the only parameter than can be directly checked.
plot_sigma_sparsebcf <- ggplot() +
  geom_path(aes(x = 1:length(c(sparsebcf_chain_1$sigma, sparsebcf_chain_2$sigma)), y = c(sparsebcf_chain_1$sigma, sparsebcf_chain_2$sigma))) +
  ggtitle("Sigma") +
  ylab("Sigma") +
  theme(axis.title.x = element_blank())


#:----------------------------------------------------------------------------------------------
# Selecting variables

## Collect variables used in the modelling of SparseBCF
variables <- hba1c.train.complete %>% select(-patid, -pated, -drugclass, -posthba1cfinal, -prop.score) %>% colnames()

# Calculate the mean values of tau inclusion proportions for all chains
if (class(try(
  
  variables_tau_original <- readRDS(paste0(output_path, "/response_model_bcf/assessment/variables_tau_original.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  variables_tau_original <- colMeans(rbind(sparsebcf_chain_1$varprb_tau, sparsebcf_chain_2$varprb_tau)) %>% t() %>% as.data.frame()
  colnames(variables_tau_original) <- variables
  
  saveRDS(variables_tau_original, paste0(output_path, "/response_model_bcf/assessment/variables_tau_original.rds"))
  
}

if (class(try(
  
  variables_tau <- readRDS(paste0(output_path, "/response_model_bcf/variables_tau.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  # select variables with tau inclusion proportions above 1/n (n is number of variables)
  variables_tau <- colnames(variables_tau_original)[variables_tau_original > 0.023]
  
  saveRDS(variables_tau, file = paste0(output_path, "/response_model_bcf/variables_tau.rds"))
  
}

# Calculate the mean values of tau inclusion proportions for each chain
if (class(try(
  
  variables_tau_original_chain_1 <- readRDS(paste0(output_path, "/response_model_bcf/assessment/variables_tau_original_chain_1.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  variables_tau_original_chain_1 <- colMeans(sparsebcf_chain_1$varprb_tau) %>% t() %>% as.data.frame()
  colnames(variables_tau_original_chain_1) <- variables
  
  saveRDS(variables_tau_original_chain_1, paste0(output_path, "/response_model_bcf/assessment/variables_tau_original_chain_1.rds"))
  
}


if (class(try(
  
  variables_tau_original_chain_2 <- readRDS(paste0(output_path, "/response_model_bcf/assessment/variables_tau_original_chain_2.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  variables_tau_original_chain_2 <- colMeans(sparsebcf_chain_2$varprb_tau) %>% t() %>% as.data.frame()
  colnames(variables_tau_original_chain_2) <- variables
  
  saveRDS(variables_tau_original_chain_2, paste0(output_path, "/response_model_bcf/assessment/variables_tau_original_chain_2.rds"))
  
}

# Plot tau inclusion proportions
plot_variables_tau <- variables_tau_original %>%
  as.data.frame() %>%
  gather(key, value) %>%
  cbind(
    variables_tau_original_chain_1 %>%
      as.data.frame() %>%
      gather(key, chain_1) %>%
      select(chain_1)
    ,
    variables_tau_original_chain_2 %>%
      as.data.frame() %>%
      gather(key, chain_2) %>%
      select(chain_2)
  ) %>%
  arrange(desc(value)) %>%
  mutate(key = factor(key),
         colour = ifelse(value > 0.023, "Above", "Below")) %>%
  ggplot() +
  geom_hline(aes(yintercept = 0.023), colour = "red") +
  geom_point(aes(x = forcats::fct_reorder(key, value), y = chain_1, alpha = 0.3)) +
  geom_point(aes(x = forcats::fct_reorder(key, value), y = chain_2, alpha = 0.3)) +
  geom_point(aes(x = forcats::fct_reorder(key, value), y = value, colour = colour)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ggtitle("tau model") +
  ylab("Posterior splitting probabilities") +
  scale_colour_manual(values = c("Above" = "green", "Below" = "red"))


# Calculate the mean values of mu inclusion proportions for all chains
if (class(try(
  
  variables_mu_original <- readRDS(paste0(output_path, "/response_model_bcf/assessment/variables_mu_original.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  variables_mu_original <- colMeans(rbind(sparsebcf_chain_1$varprb_mu, sparsebcf_chain_2$varprb_mu)) %>% t() %>% as.data.frame()
  colnames(variables_mu_original) <- c(variables, "propensity score")
  
  saveRDS(variables_mu_original, paste0(output_path, "/response_model_bcf/assessment/variables_mu_original.rds"))
  
}


if (class(try(
  
  variables_mu <- readRDS(paste0(output_path, "/response_model_bcf/variables_mu.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  # select variables with mu inclusion proportions above 1/n (n is number of variables)
  dataset.variables.mu <- variables_mu_original %>%
    select(-`propensity score`)
  
  variables_mu <- colnames(dataset.variables.mu)[dataset.variables.mu > 0.023]
  
  saveRDS(variables_mu, file = paste0(output_path, "/response_model_bcf/variables_mu.rds"))
  
}

# Calculate the mean values of mu inclusion proportions for each chain
if (class(try(
  
  variables_mu_original_chain_1 <- readRDS(paste0(output_path, "/response_model_bcf/assessment/variables_mu_original_chain_1.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  variables_mu_original_chain_1 <- colMeans(sparsebcf_chain_1$varprb_mu) %>% t() %>% as.data.frame()
  colnames(variables_mu_original_chain_1) <- c(variables, "propensity score")
  
  saveRDS(variables_mu_original_chain_1, paste0(output_path, "/response_model_bcf/assessment/variables_mu_original_chain_1.rds"))
  
}

if (class(try(
  
  variables_mu_original_chain_2 <- readRDS(paste0(output_path, "/response_model_bcf/assessment/variables_mu_original_chain_2.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  variables_mu_original_chain_2 <- colMeans(sparsebcf_chain_2$varprb_mu) %>% t() %>% as.data.frame()
  colnames(variables_mu_original_chain_2) <- c(variables, "propensity score")
  
  saveRDS(variables_mu_original_chain_2, paste0(output_path, "/response_model_bcf/assessment/variables_mu_original_chain_2.rds"))
  
}

# Plot mu inclusion proportions
plot_variables_mu <- variables_mu_original %>%
  as.data.frame() %>%
  gather(key, value) %>%
  cbind(
    variables_mu_original_chain_1 %>%
      as.data.frame() %>%
      gather(key, chain_1) %>%
      select(chain_1)
    ,
    variables_mu_original_chain_2 %>%
      as.data.frame() %>%
      gather(key, chain_2) %>%
      select(chain_2)
  ) %>%
  arrange(desc(value)) %>%
  mutate(key = factor(key),
         colour = ifelse(value > 0.023, "Above", "Below")) %>%
  ggplot() +
  geom_hline(aes(yintercept = 0.023), colour = "red") +
  geom_point(aes(x = forcats::fct_reorder(key, value), y = chain_1, alpha = 0.3)) +
  geom_point(aes(x = forcats::fct_reorder(key, value), y = chain_2, alpha = 0.3)) +
  geom_point(aes(x = forcats::fct_reorder(key, value), y = value, colour = colour)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ggtitle("mu model") +
  ylab("Posterior splitting probabilities") +
  scale_colour_manual(values = c("Above" = "green", "Below" = "red"))


# remove objects that are not needed anymore (objects too big for the session)
rm(sparsebcf_chain_1)
rm(sparsebcf_chain_2)



#:----------------------------------------------------------------------------------------------
## Prepare dataset by only selecting the variable chosen (in both prognostic and moderator sides of the model)

hba1c.train.complete.vs <- hba1c.train.complete %>%
  select(patid, pated, posthba1cfinal, drugclass, unique(c(variables_mu, variables_tau)), prop.score)

## Hold out cohort
hba1c.test <- set_up_data_sglt2_glp1(dataset.type = "hba1c.test")

# collect propensity score values
patient_prop_scores <- readRDS(paste0(output_path, "/ps_model/patient_prop_scores.rds"))

hba1c.test <- hba1c.test %>%
  left_join(patient_prop_scores, by = c("patid", "pated"))

# only keep entries with complete data
hba1c.test.complete.vs <- hba1c.test %>%
  # selected variables from SparseBCF
  select(patid, pated, posthba1cfinal, drugclass, unique(c(variables_mu, variables_tau)), prop.score) %>%
  # only complete cases
  drop_na() %>%
  as.data.frame()


#:----------------------------------------------------------------------------------------------------
# Fit BCF model with propensity score included

if (class(try(
  
  bcf_model_prop <- readRDS(paste0(output_path, "/response_model_bcf/bcf_model_prop.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  
  bcf_model_prop = bcf::bcf(y = hba1c.train.complete.vs$posthba1cfinal,
                            z = hba1c.train.complete.vs %>%
                              select(drugclass) %>%
                              mutate(drugclass = ifelse(drugclass == "GLP1", 0, 1)) %>%
                              unlist(),
                            x_control = hba1c.train.complete.vs %>%
                              select(
                                all_of(variables_mu)
                              ) %>%
                              mutate_all(funs(as.numeric(.))) %>%
                              as.matrix(),
                            x_moderate = hba1c.train.complete.vs %>%
                              select(
                                all_of(variables_tau)
                              ) %>%
                              mutate_all(funs(as.numeric(.))) %>%
                              as.matrix(),
                            pihat = 1-hba1c.train.complete.vs$prop.score,
                            nburn = 200000,
                            nsim = 25000,
                            nthin = 4,
                            n_chains = 2,
                            # n_threads was max((RcppParallel::defaultNumThreads() - 2)/n_cores, 1) (this uses all of the server)
                            n_threads = 4,
                            update_interval = 500,
                            ntree_control = 200,
                            sd_control = 2 * sd(hba1c.train.complete.vs$posthba1cfinal),
                            base_control = 0.95,
                            power_control = 2,
                            ntree_moderate = 200,
                            sd_moderate = 2 * sd(hba1c.train.complete.vs$posthba1cfinal),
                            base_moderate = 0.95,
                            power_moderate = 2,
                            use_muscale = FALSE,
                            use_tauscale = FALSE,
                            save_tree_directory = paste0(output_path, "/response_model_bcf/trees_prop"))
  
  saveRDS(bcf_model_prop, paste0(output_path, "/response_model_bcf/bcf_model_prop.rds"))
  
}


if (class(try(
  
  predictions.hba1c.test_prop <- readRDS(paste0(output_path, "/response_model_bcf/predictions.hba1c.test_prop.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  predictions.hba1c.test_prop <- predict(object = bcf_model_prop,
                                         x_predict_control = hba1c.test.complete.vs %>%
                                           select(
                                             all_of(variables_mu)
                                           ) %>%
                                           mutate_all(funs(as.numeric(.))) %>%
                                           as.matrix(),
                                         x_predict_moderate = hba1c.test.complete.vs %>%
                                           select(
                                             all_of(variables_tau)
                                           ) %>%
                                           mutate_all(funs(as.numeric(.))) %>%
                                           as.matrix(),
                                         pi_pred = 1-hba1c.test.complete.vs$prop.score,
                                         z_pred = hba1c.test.complete.vs %>%
                                           select(drugclass) %>%
                                           mutate(drugclass = ifelse(drugclass == "GLP1", 0, 1)) %>%
                                           unlist(),
                                         save_tree_directory = paste0(output_path, "/response_model_bcf/trees_prop"))
  
  
  saveRDS(predictions.hba1c.test_prop, paste0(output_path, "/response_model_bcf/predictions.hba1c.test_prop.rds"))
  
}

#:----------------------------------------------------------------------------------------------------
# Fit BCF model without propensity score included

if (class(try(
  
  bcf_model <- readRDS(paste0(output_path, "/response_model_bcf/bcf_model.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  
  bcf_model = bcf::bcf(y = hba1c.train.complete.vs$posthba1cfinal,
                       z = hba1c.train.complete.vs %>%
                         select(drugclass) %>%
                         mutate(drugclass = ifelse(drugclass == "GLP1", 0, 1)) %>%
                         unlist(),
                       x_control = hba1c.train.complete.vs %>%
                         select(
                           all_of(variables_mu)
                         ) %>%
                         mutate_all(funs(as.numeric(.))) %>%
                         as.matrix(),
                       x_moderate = hba1c.train.complete.vs %>%
                         select(
                           all_of(variables_tau)
                         ) %>%
                         mutate_all(funs(as.numeric(.))) %>%
                         as.matrix(),
                       pihat = 1-hba1c.train.complete.vs$prop.score,
                       nburn = 200000,
                       nsim = 25000,
                       nthin = 4,
                       n_chains = 2,
                       # n_threads was max((RcppParallel::defaultNumThreads() - 2)/n_cores, 1) (this uses all of the server)
                       n_threads = 4,
                       update_interval = 500,
                       ntree_control = 200,
                       sd_control = 2 * sd(hba1c.train.complete.vs$posthba1cfinal),
                       base_control = 0.95,
                       power_control = 2,
                       ntree_moderate = 200,
                       sd_moderate = 2 * sd(hba1c.train.complete.vs$posthba1cfinal),
                       base_moderate = 0.95,
                       power_moderate = 2,
                       use_muscale = FALSE,
                       use_tauscale = FALSE,
                       include_pi = "none",
                       save_tree_directory = paste0(output_path, "/response_model_bcf/trees_no_prop"))
  
  saveRDS(bcf_model, paste0(output_path, "/response_model_bcf/bcf_model.rds"))
  
}


if (class(try(
  
  predictions.hba1c.test <- readRDS(paste0(output_path, "/response_model_bcf/predictions.hba1c.test.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  predictions.hba1c.test <- predict(object = bcf_model,
                                    x_predict_control = hba1c.test.complete.vs %>%
                                      select(
                                        all_of(variables_mu)
                                      ) %>%
                                      mutate_all(funs(as.numeric(.))) %>%
                                      as.matrix(),
                                    x_predict_moderate = hba1c.test.complete.vs %>%
                                      select(
                                        all_of(variables_tau)
                                      ) %>%
                                      mutate_all(funs(as.numeric(.))) %>%
                                      as.matrix(),
                                    pi_pred = 1-hba1c.test.complete.vs$prop.score,
                                    z_pred = hba1c.test.complete.vs %>%
                                      select(drugclass) %>%
                                      mutate(drugclass = ifelse(drugclass == "GLP1", 0, 1)) %>%
                                      unlist(),
                                    save_tree_directory = paste0(output_path, "/response_model_bcf/trees_no_prop"))
  
  
  saveRDS(predictions.hba1c.test, paste0(output_path, "/response_model_bcf/predictions.hba1c.test.rds"))
  
}


#:----------------------------------------------------------------------------------------------------


# plot mu predictions from model BCF (prop score) vs BCF (no prop score)
if (class(try(
  
  prop.score.comparison.mu <- readRDS(paste0(output_path, "/response_model_bcf/prop.score.comparison.mu.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  prop.score.comparison.mu <- cbind(bcf_prop = c(colMeans(bcf_model_prop$mu), colMeans(predictions.hba1c.test_prop$mu)),
                                    bcf_no_prop = c(colMeans(bcf_model$mu), colMeans(predictions.hba1c.test$mu)))
  
  saveRDS(prop.score.comparison.mu, paste0(output_path, "/response_model_bcf/prop.score.comparison.mu.rds"))
  
}

# plot tau predictions from model BCF (prop score) vs BCF (no prop score)
if (class(try(
  
  prop.score.comparison <- readRDS(paste0(output_path, "/response_model_bcf/prop.score.comparison.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  prop.score.comparison <- cbind(bcf_prop = c(colMeans(bcf_model_prop$tau), colMeans(predictions.hba1c.test_prop$tau)),
                                 bcf_no_prop = c(colMeans(bcf_model$tau), colMeans(predictions.hba1c.test$tau)))
  
  saveRDS(prop.score.comparison, paste0(output_path, "/response_model_bcf/prop.score.comparison.rds"))
  
}


## Plot the differences between BCF with prop score vs BCF without prop score
plot_comparison_bcf_models <- patchwork::wrap_plots(
  
  # bar plot of the numbers in each group
  prop.score.comparison %>%
    as.data.frame() %>%
    mutate(effect = ifelse(bcf_prop > 0 & bcf_no_prop > 0, "BCF P > 0 / BCF NP > 0",
                           ifelse(bcf_prop < 0 & bcf_no_prop < 0, "BCF P < 0 / BCF NP < 0",
                                  ifelse(bcf_prop > 0 & bcf_no_prop < 0, "BCF P > 0 / BCF NP < 0",
                                         "BCF P < 0 / BCF NP > 0")))) %>%
    count(effect) %>%
    mutate(perc = n / nrow(prop.score.comparison) *100) %>%
    ggplot() +
    geom_bar(aes(x = effect, y = perc), stat = "identity") +
    geom_text(aes(label = n, y= ..prop.., x = effect), stat= "count", vjust = -.5, colour = "blue")
  
  ,
  
  patchwork::wrap_plots(
    
    # scatter plot of tau values
    prop.score.comparison %>%
      as.data.frame() %>%
      ggplot(aes(x = bcf_prop, y = bcf_no_prop)) +
      geom_vline(aes(xintercept = 0), colour = "red") +
      geom_hline(aes(yintercept = 0), colour = "red") +
      geom_point() +
      ggtitle("Tau") +
      xlab("BCF model with propensity score") +
      ylab("BCF model without propensity score")
    
    ,
    
    # scatter plot of mu values
    prop.score.comparison.mu %>%
      as.data.frame() %>%
      ggplot(aes(x = bcf_prop, y = bcf_no_prop)) +
      geom_point() +
      ggtitle("Mu") +
      xlab("BCF model with propensity score") +
      ylab("BCF model without propensity score")
    
    ,
    
    # density plot of residuals for tau (BCF without prop - BCF with prop)
    prop.score.comparison %>%
      as.data.frame() %>%
      mutate(effect.difference = bcf_no_prop - bcf_prop) %>%
      ggplot() +
      geom_density(aes(x = effect.difference)) +
      xlab("Predicted effect difference (BCF NP - BCF P)")
    
    ,
    
    # density plot of residuals for mu (BCF without prop - BCF with prop)
    prop.score.comparison.mu %>%
      as.data.frame() %>%
      mutate(effect.difference = bcf_no_prop - bcf_prop) %>%
      ggplot() +
      geom_density(aes(x = effect.difference)) +
      xlab("Predicted effect difference (BCF NP - BCF P)")
    
    , ncol = 2, nrow = 2
    
    
    
  ) , ncol = 1, nrow = 2
  
)


## Produce pdf with these plots
pdf("Plots/11.05.BCF_prop_vs_no_prop.pdf", width = 8, height = 12)
plot_comparison_bcf_models
dev.off()

# remove objects from memory as they are not needed anymore
rm(bcf_model_prop)
rm(predictions.hba1c.test_prop)


# Trace plot for the only parameter directly observable
plot_sigma_bcf <- ggplot() +
  geom_path(aes(x = 1:length(bcf_model$sigma), y = bcf_model$sigma)) +
  ggtitle("Sigma") +
  ylab("Sigma") +
  theme(axis.title.x = element_blank())


#:----------------------------------------------------------------------------------------------------

#### Residuals for the BCF no prop model

# Calculate prediction interval, residual interval and standardised residual interval for Development cohort
cred_pred_dev <- calc_resid(hba1c.train.complete.vs, bcf_model$mu, "posthba1cfinal")

# Calculate prediction interval, residual interval and standardised residual interval for Validation cohort
cred_pred_val <- calc_resid(hba1c.test.complete.vs, predictions.hba1c.test$mu, "posthba1cfinal")

# Plot the standardised residual plots
plot_residuals <- resid_plot(cred_pred_dev, cred_pred_val, "Standardised Residuals of BCF model")


#### R2/RMSE
if (class(try(
  
  assessment_dev <- readRDS(paste0(output_path, "/response_model_bcf/assessment/assessment_dev.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  assessment_dev <- calc_assessment(hba1c.train.complete.vs, bcf_model$mu, "posthba1cfinal")
  
  saveRDS(assessment_dev, paste0(output_path, "/response_model_bcf/assessment/assessment_dev.rds"))
  
}

if (class(try(
  
  assessment_val <- readRDS(paste0(output_path, "/response_model_bcf/assessment/assessment_val.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  assessment_val <- calc_assessment(hba1c.test.complete.vs, predictions.hba1c.test$mu, "posthba1cfinal")
  
  saveRDS(assessment_val, paste0(output_path, "/response_model_bcf/assessment/assessment_val.rds"))
  
}

#:----------------------------------------------------------------------------------------------------

#### Plot predicted treatment effects

data_dev <- cbind(mean = colMeans(bcf_model$tau)) %>%
  as.data.frame()

plot_effect_1 <- hist_plot(data_dev, "", -20, 30)

data_val <- cbind(mean = colMeans(predictions.hba1c.test$tau)) %>%
  as.data.frame()

plot_effect_2 <- hist_plot(data_val, "", -20, 30)

####### Strata effect
## male sex
sex_male_dev <- hba1c.train.complete.vs %>%
  select(sex) %>%
  cbind(mean = colMeans(bcf_model$tau)) %>%
  filter(sex == "Male") %>%
  select(-sex)

plot_effect_1_male <- hist_plot(sex_male_dev, "", -20, 30)

## female sex
sex_female_dev <- hba1c.train.complete.vs %>%
  select(sex) %>%
  cbind(mean = colMeans(bcf_model$tau)) %>%
  filter(sex == "Female") %>%
  select(-sex)

plot_effect_1_female <- hist_plot(sex_female_dev, "", -20, 30)


## male sex
sex_male_val <- hba1c.test.complete.vs %>%
  select(sex) %>%
  cbind(mean = colMeans(predictions.hba1c.test$tau)) %>%
  filter(sex == "Male") %>%
  select(-sex)

plot_effect_2_male <- hist_plot(sex_male_val, "", -20, 30)

## female sex
sex_female_val <- hba1c.test.complete.vs %>%
  select(sex) %>%
  cbind(mean = colMeans(predictions.hba1c.test$tau)) %>%
  filter(sex == "Female") %>%
  select(-sex)

plot_effect_2_female <- hist_plot(sex_female_val, "", -20, 30)


#:----------------------------------------------------------------------------------------------------

####### Validation

### Development cohort
predicted_observed_dev <- hba1c.train.complete.vs %>%
  cbind(hba1c_diff = colMeans(bcf_model$tau)) %>%
  mutate(bestdrug = ifelse(hba1c_diff < 0, "SGLT2", "GLP1"),
         hba1c_diff.q = ntile(hba1c_diff, 10))

#:--------- Overall average treatment effect
# PSM 1:1 posthba1c ~ drugclass
if (class(try(
  
  ATE_matching_1_1_all_dev <- readRDS(paste0(output_path, "/response_model_bcf/assessment/ATE_matching_1_1_all_dev.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_matching_1_1_all_dev <- calc_ATE(data = predicted_observed_dev%>%mutate(hba1c_diff.q=1), 
                                       validation_type = "PSM", variable = "posthba1cfinal", 
                                       quantile_var = "hba1c_diff.q", prop_scores = predicted_observed_dev$prop.score, 
                                       order = "largest", breakdown = unique(c(variables_tau, variables_mu)))
  
  saveRDS(ATE_matching_1_1_all_dev, paste0(output_path, "/response_model_bcf/assessment/ATE_matching_1_1_all_dev.rds"))
  
}

### plot love plots for the matching
plot_all_dev_love <- cobalt::love.plot(ATE_matching_1_1_all_dev$matching_outputs[[1]], binary = "std", thresholds = c(m = .1), sample.names = c("Unmatched", "Matched"), title = "Full dataset")

### plot calibration plot
# plot_all_dev <- ATE_plot(ATE_matching_1_1_all_dev[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci")

# PSM 1:1 + adjust posthba1c ~ drugclass + adjust
ATE_matching_1_1_adjust_all_dev <- calc_ATE(data = predicted_observed_dev%>%mutate(hba1c_diff.q=1), 
                                            validation_type = "PSM + adjust", variable = "posthba1cfinal", 
                                            quantile_var = "hba1c_diff.q", prop_scores = predicted_observed_dev$prop.score, 
                                            order = "largest", breakdown = unique(c(variables_tau, variables_mu)))

# PDF with the love plot for overall dataset
pdf("Plots/11.05.dev_full_love_plot.pdf")
plot_all_dev_love
dev.off()

# PSM 1:1 posthba1c ~ drugclass
if (class(try(
  
  ATE_matching_1_1_validation_dev <- readRDS(paste0(output_path, "/response_model_bcf/assessment/ATE_matching_1_1_validation_dev.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_matching_1_1_validation_dev <- calc_ATE(data = predicted_observed_dev,
                                              validation_type = "PSM", variable = "posthba1cfinal",
                                              quantile_var = "hba1c_diff.q", prop_scores = predicted_observed_dev$prop.score, 
                                              order = "largest", breakdown = unique(c(variables_tau, variables_mu)))
  
  
  saveRDS(ATE_matching_1_1_validation_dev, paste0(output_path, "/response_model_bcf/assessment/ATE_matching_1_1_validation_dev.rds"))
  
}

# Plot calibration of treatment effects
plot_ATE_dev_prop_score_matching_1_1 <- ATE_plot(ATE_matching_1_1_validation_dev[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci")

# PDF of love plots for all deciles
pdf("Plots/11.05.dev_PSM_1_1_dev.pdf")
cobalt::love.plot(ATE_matching_1_1_validation_dev$matching_outputs[[1]], binary = "std", thresholds = c(m = .1), sample.names = c("Unmatched", "Matched"), title = "Match 1")
cobalt::love.plot(ATE_matching_1_1_validation_dev$matching_outputs[[2]], binary = "std", thresholds = c(m = .1), sample.names = c("Unmatched", "Matched"), title = "Match 2")
cobalt::love.plot(ATE_matching_1_1_validation_dev$matching_outputs[[3]], binary = "std", thresholds = c(m = .1), sample.names = c("Unmatched", "Matched"), title = "Match 3")
cobalt::love.plot(ATE_matching_1_1_validation_dev$matching_outputs[[4]], binary = "std", thresholds = c(m = .1), sample.names = c("Unmatched", "Matched"), title = "Match 4")
cobalt::love.plot(ATE_matching_1_1_validation_dev$matching_outputs[[5]], binary = "std", thresholds = c(m = .1), sample.names = c("Unmatched", "Matched"), title = "Match 5")
cobalt::love.plot(ATE_matching_1_1_validation_dev$matching_outputs[[6]], binary = "std", thresholds = c(m = .1), sample.names = c("Unmatched", "Matched"), title = "Match 6")
cobalt::love.plot(ATE_matching_1_1_validation_dev$matching_outputs[[7]], binary = "std", thresholds = c(m = .1), sample.names = c("Unmatched", "Matched"), title = "Match 7")
cobalt::love.plot(ATE_matching_1_1_validation_dev$matching_outputs[[8]], binary = "std", thresholds = c(m = .1), sample.names = c("Unmatched", "Matched"), title = "Match 8")
cobalt::love.plot(ATE_matching_1_1_validation_dev$matching_outputs[[9]], binary = "std", thresholds = c(m = .1), sample.names = c("Unmatched", "Matched"), title = "Match 9")
cobalt::love.plot(ATE_matching_1_1_validation_dev$matching_outputs[[10]], binary = "std", thresholds = c(m = .1), sample.names = c("Unmatched", "Matched"), title = "Match 10")
plot_ATE_dev_prop_score_matching_1_1
dev.off()

# PSM 1:1 posthba1c ~ drugclass + adjust (all variables used in bcf)
if (class(try(
  
  ATE_matching_1_1_adjust_validation_dev <- readRDS(paste0(output_path, "/response_model_bcf/assessment/ATE_matching_1_1_adjust_validation_dev.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_matching_1_1_adjust_validation_dev <- calc_ATE(data = predicted_observed_dev,
                                                     validation_type = "PSM + adjust", variable = "posthba1cfinal",
                                                     quantile_var = "hba1c_diff.q", prop_scores = predicted_observed_dev$prop.score, 
                                                     order = "largest", breakdown = unique(c(variables_tau, variables_mu)))
  
  saveRDS(ATE_matching_1_1_adjust_validation_dev, paste0(output_path, "/response_model_bcf/assessment/ATE_matching_1_1_adjust_validation_dev.rds"))
  
}

# Plot calibration of treatment effects
plot_ATE_dev_prop_score_matching_1_1_adjust <- ATE_plot(ATE_matching_1_1_adjust_validation_dev[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci")

# posthba1c ~ drugclass + adjust (all variables used in bcf)
if (class(try(
  
  ATE_adjust_validation_dev <- readRDS(paste0(output_path, "/response_model_bcf/assessment/ATE_adjust_validation_dev.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_adjust_validation_dev <- calc_ATE(data = predicted_observed_dev,
                                        validation_type = "Adjust", variable = "posthba1cfinal",
                                        quantile_var = "hba1c_diff.q",
                                        order = "largest", breakdown = unique(c(variables_tau, variables_mu)))
  
  saveRDS(ATE_adjust_validation_dev, paste0(output_path, "/response_model_bcf/assessment/ATE_adjust_validation_dev.rds"))
}

# Plot calibration of treatment effects
plot_ATE_adjust_validation_dev <- ATE_plot(ATE_adjust_validation_dev[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci")


# # PSM k:1 posthba1c ~ drugclass
# ATE_matching_k_1_validation_dev <- calc_ATE(data = predicted_observed_dev,
#                                             validation_type = "PSM", variable = "posthba1cfinal",
#                                             quantile_var = "hba1c_diff.q", prop_scores = predicted_observed_dev$prop.score, 
#                                             order = "largest", breakdown = unique(c(variables_tau, variables_mu)), replace = TRUE)
# 
# plot_ATE_dev_prop_score_matching_k_1 <- ATE_plot(ATE_matching_k_1_validation_dev[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci")
# 
# # PSM k:1 posthba1c ~ drugclass + adjust (all variables used in bcf)
# ATE_matching_k_1_adjust_validation_dev <- calc_ATE(data = predicted_observed_dev,
#                                                    validation_type = "PSM + adjust", variable = "posthba1cfinal",
#                                                    quantile_var = "hba1c_diff.q", prop_scores = predicted_observed_dev$prop.score, 
#                                                    order = "largest", breakdown = unique(c(variables_tau, variables_mu)), replace = TRUE)
# 
# plot_ATE_dev_prop_score_matching_k_1_adjust <- ATE_plot(ATE_matching_k_1_adjust_validation_dev[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci")
# 
# # IPTW
# ATE_weighting_validation_dev <- calc_ATE_validation_inverse_prop_weighting(predicted_observed_dev, "posthba1cfinal", predicted_observed_dev$prop.score)
# 
# plot_ATE_dev_prop_score_weighting  <- ATE_plot(ATE_weighting_validation_dev[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci")
# 
# # IPTW stabilised
# ATE_stabilised_validation_dev <- calc_ATE_validation_inverse_prop_weighting_stabilised(predicted_observed_dev, "posthba1cfinal", predicted_observed_dev$prop.score)
# 
# plot_ATE_dev_prop_score_stabilised  <- ATE_plot(ATE_stabilised_validation_dev[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci")



### Validation cohort
predicted_observed_val <- hba1c.test.complete.vs %>%
  cbind(hba1c_diff = colMeans(predictions.hba1c.test$tau)) %>%
  mutate(bestdrug = ifelse(hba1c_diff < 0, "SGLT2", "GLP1"),
         hba1c_diff.q = ntile(hba1c_diff, 10))

#:--------- Overall average treatment effect
# PSM 1:1 posthba1c ~ drugclass
if (class(try(
  
  ATE_matching_1_1_all_val <- readRDS(paste0(output_path, "/response_model_bcf/assessment/ATE_matching_1_1_all_val.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_matching_1_1_all_val <- calc_ATE(data = predicted_observed_val%>%mutate(hba1c_diff.q=1), 
                                       validation_type = "PSM", variable = "posthba1cfinal", 
                                       quantile_var = "hba1c_diff.q", prop_scores = predicted_observed_val$prop.score, 
                                       order = "largest", breakdown = unique(c(variables_tau, variables_mu)))
  
  saveRDS(ATE_matching_1_1_all_val, paste0(output_path, "/response_model_bcf/assessment/ATE_matching_1_1_all_val.rds"))
  
}

### plot love plots for the matching
plot_all_val_love <- cobalt::love.plot(ATE_matching_1_1_all_val$matching_outputs[[1]], binary = "std", thresholds = c(m = .1), sample.names = c("Unmatched", "Matched"), title = "Full dataset")

### plot calibration plot
# plot_all_val <- ATE_plot(ATE_matching_1_1_all_val[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci")

# PSM 1:1 + adjust posthba1c ~ drugclass + adjust
ATE_matching_1_1_adjust_all_val <- calc_ATE(data = predicted_observed_val%>%mutate(hba1c_diff.q=1), 
                                            validation_type = "PSM + adjust", variable = "posthba1cfinal", 
                                            quantile_var = "hba1c_diff.q", prop_scores = predicted_observed_val$prop.score, 
                                            order = "largest", breakdown = unique(c(variables_tau, variables_mu)))

# PDF with the love plot for overall dataset
pdf("Plots/11.05.val_full_love_plot.pdf")
plot_all_val_love
dev.off()

# PSM 1:1 posthba1c ~ drugclass
if (class(try(
  
  ATE_matching_1_1_validation_val <- readRDS(paste0(output_path, "/response_model_bcf/assessment/ATE_matching_1_1_validation_val.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_matching_1_1_validation_val <- calc_ATE(data = predicted_observed_val,
                                              validation_type = "PSM", variable = "posthba1cfinal",
                                              quantile_var = "hba1c_diff.q", prop_scores = predicted_observed_val$prop.score, 
                                              order = "largest", breakdown = unique(c(variables_tau, variables_mu)))
  
  saveRDS(ATE_matching_1_1_validation_val, paste0(output_path, "/response_model_bcf/assessment/ATE_matching_1_1_validation_val.rds"))
  
}

# Plot calibration of treatment effects
plot_ATE_val_prop_score_matching_1_1 <- ATE_plot(ATE_matching_1_1_validation_val[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci")

# PDF of love plots for all deciles
pdf("Plots/11.05.val_PSM_1_1_val.pdf")
cobalt::love.plot(ATE_matching_1_1_validation_val$matching_outputs[[1]], binary = "std", thresholds = c(m = .1), sample.names = c("Unmatched", "Matched"), title = "Match 1")
cobalt::love.plot(ATE_matching_1_1_validation_val$matching_outputs[[2]], binary = "std", thresholds = c(m = .1), sample.names = c("Unmatched", "Matched"), title = "Match 2")
cobalt::love.plot(ATE_matching_1_1_validation_val$matching_outputs[[3]], binary = "std", thresholds = c(m = .1), sample.names = c("Unmatched", "Matched"), title = "Match 3")
cobalt::love.plot(ATE_matching_1_1_validation_val$matching_outputs[[4]], binary = "std", thresholds = c(m = .1), sample.names = c("Unmatched", "Matched"), title = "Match 4")
cobalt::love.plot(ATE_matching_1_1_validation_val$matching_outputs[[5]], binary = "std", thresholds = c(m = .1), sample.names = c("Unmatched", "Matched"), title = "Match 5")
cobalt::love.plot(ATE_matching_1_1_validation_val$matching_outputs[[6]], binary = "std", thresholds = c(m = .1), sample.names = c("Unmatched", "Matched"), title = "Match 6")
cobalt::love.plot(ATE_matching_1_1_validation_val$matching_outputs[[7]], binary = "std", thresholds = c(m = .1), sample.names = c("Unmatched", "Matched"), title = "Match 7")
cobalt::love.plot(ATE_matching_1_1_validation_val$matching_outputs[[8]], binary = "std", thresholds = c(m = .1), sample.names = c("Unmatched", "Matched"), title = "Match 8")
cobalt::love.plot(ATE_matching_1_1_validation_val$matching_outputs[[9]], binary = "std", thresholds = c(m = .1), sample.names = c("Unmatched", "Matched"), title = "Match 9")
cobalt::love.plot(ATE_matching_1_1_validation_val$matching_outputs[[10]], binary = "std", thresholds = c(m = .1), sample.names = c("Unmatched", "Matched"), title = "Match 10")
plot_ATE_val_prop_score_matching_1_1
dev.off()

# PSM 1:1 posthba1c ~ drugclass + adjust (all variables used in bcf)
if (class(try(
  
  ATE_matching_1_1_adjust_validation_val <- readRDS(paste0(output_path, "/response_model_bcf/assessment/ATE_matching_1_1_adjust_validation_val.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_matching_1_1_adjust_validation_val <- calc_ATE(data = predicted_observed_val,
                                                     validation_type = "PSM + adjust", variable = "posthba1cfinal",
                                                     quantile_var = "hba1c_diff.q", prop_scores = predicted_observed_val$prop.score, 
                                                     order = "largest", breakdown = unique(c(variables_tau, variables_mu)))
  
  saveRDS(ATE_matching_1_1_adjust_validation_val, paste0(output_path, "/response_model_bcf/assessment/ATE_matching_1_1_adjust_validation_val.rds"))
  
}

# Plot calibration of treatment effects
plot_ATE_val_prop_score_matching_1_1_adjust <- ATE_plot(ATE_matching_1_1_adjust_validation_val[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci")

# posthba1c ~ drugclass + adjust (all variables used in bcf)
if (class(try(
  
  ATE_adjust_validation_val <- readRDS(paste0(output_path, "/response_model_bcf/assessment/ATE_adjust_validation_val.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_adjust_validation_val <- calc_ATE(data = predicted_observed_val,
                                        validation_type = "Adjust", variable = "posthba1cfinal",
                                        quantile_var = "hba1c_diff.q",
                                        order = "largest", breakdown = unique(c(variables_tau, variables_mu)))
  
  saveRDS(ATE_adjust_validation_val, paste0(output_path, "/response_model_bcf/assessment/ATE_adjust_validation_val.rds"))
}

# Plot calibration of treatment effects
plot_ATE_adjust_validation_val <- ATE_plot(ATE_adjust_validation_val[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci")


# # PSM k:1 posthba1c ~ drugclass
# ATE_matching_k_1_validation_val <- calc_ATE(data = predicted_observed_val,
#                                             validation_type = "PSM", variable = "posthba1cfinal",
#                                             quantile_var = "hba1c_diff.q", prop_scores = predicted_observed_val$prop.score, 
#                                             order = "largest", breakdown = unique(c(variables_tau, variables_mu)), replace = TRUE)
# 
# plot_ATE_val_prop_score_matching_k_1 <- ATE_plot(ATE_matching_k_1_validation_val[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci")
# 
# # PSM k:1 posthba1c ~ drugclass + adjust (all variables used in bcf)
# ATE_matching_k_1_adjust_validation_val <- calc_ATE(data = predicted_observed_val,
#                                                    validation_type = "PSM + adjust", variable = "posthba1cfinal",
#                                                    quantile_var = "hba1c_diff.q", prop_scores = predicted_observed_val$prop.score, 
#                                                    order = "largest", breakdown = unique(c(variables_tau, variables_mu)), replace = TRUE)
# 
# plot_ATE_val_prop_score_matching_k_1_adjust <- ATE_plot(ATE_matching_k_1_adjust_validation_val[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci")
# 
# # IPTW
# ATE_weighting_validation_val <- calc_ATE_validation_inverse_prop_weighting(predicted_observed_val, "posthba1cfinal", predicted_observed_val$prop.score)
# 
# plot_ATE_val_prop_score_weighting  <- ATE_plot(ATE_weighting_validation_val[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci")
# 
# # IPTW stabilised
# ATE_stabilised_validation_val <- calc_ATE_validation_inverse_prop_weighting_stabilised(predicted_observed_val, "posthba1cfinal", predicted_observed_val$prop.score)
# 
# plot_ATE_val_prop_score_stabilised  <- ATE_plot(ATE_stabilised_validation_val[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci")



## PDF with all of the plots necessary for checking the model fitting and validation
pdf(width = 10, height = 8, file = "Plots/11.05.slade_aurum_bcf.pdf")

# variables for treatment effects
plot_variables_tau
# variables for response
plot_variables_mu

# comparison between BCF models with/without prop scores
plot_comparison_bcf_models

# Standardised residuals for the bcf model
plot_residuals

patchwork::wrap_plots(list(plot_effect_1, plot_effect_2)) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(tag_levels = "A", # labels A = development, B = validation
                             title = "Treatment effect heterogeneity (A-Developement)", # title of full plot
                             theme = theme(plot.title = element_text(hjust = 0.5),
                                           legend.position = "bottom")) # center title of full plot

patchwork::wrap_plots(list(plot_effect_1_male, plot_effect_1_female), ncol = 2) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(tag_levels = "A",
                             title = "Developement: Treatment effect heterogeneity (A-Male)", # title of full plot
                             theme = theme(plot.title = element_text(hjust = 0.5),
                                           legend.position = "bottom")) # center title of full plot


patchwork::wrap_plots(list(plot_effect_2_male, plot_effect_2_female), ncol = 2) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(tag_levels = "A",
                             title = "Validation: Treatment effect heterogeneity (A-Male)", # title of full plot
                             theme = theme(plot.title = element_text(hjust = 0.5),
                                           legend.position = "bottom")) # center title of full plot

# validation of treatment effects for development cohort
patchwork::wrap_plots(list(plot_ATE_dev_prop_score_matching_1_1 + ggtitle("Match 1:1"),
                           plot_ATE_dev_prop_score_matching_1_1_adjust + ggtitle("Match 1:1 adjusted"),
                           plot_ATE_adjust_validation_dev + ggtitle("Adjusted")), ncol = 3) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(tag_levels = "A",
                             title = "Training cohort: Treatment effect validation", # title of full plot
                             theme = theme(plot.title = element_text(hjust = 0.5),
                                           legend.position = "bottom")) # center title of full plot

# validation of treatment effects for validation cohort
patchwork::wrap_plots(list(plot_ATE_val_prop_score_matching_1_1 + ggtitle("Match 1:1"),
                           plot_ATE_val_prop_score_matching_1_1_adjust + ggtitle("Match 1:1 adjusted"),
                           plot_ATE_adjust_validation_val + ggtitle("Adjusted")), ncol = 3) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(tag_levels = "A",
                             title = "Hold-out cohort: Treatment effect validation", # title of full plot
                             theme = theme(plot.title = element_text(hjust = 0.5),
                                           legend.position = "bottom")) # center title of full plot


dev.off()

#:-------------------------------------------------------------------------------------
## Decision tree - pseudo variable importance

library(rpart)
library(rattle)
library(rpart.plot)

rpart.dataset <- predicted_observed_dev

fit <- rpart(hba1c_diff ~ agetx + sex + ncurrtx + prehba1c + prebmi + preegfr + preheartfailure + preihd + preneuropathy + prepad + preretinopathy, data = rpart.dataset)

rpart.dataset.strata <- predicted_observed_dev %>%
  mutate(hba1cmonth = round(hba1cmonth)) %>%
  filter(hba1cmonth > 6)

fit2 <- rpart(hba1c_diff ~ agetx + sex + ncurrtx + prehba1c + prebmi + preegfr + preheartfailure + preihd + preneuropathy + prepad + preretinopathy, data = rpart.dataset.strata)


pdf(width = 20, height = 8, file = "Plots/11.05.effect_decision_tree.pdf")

prp(fit, pal.thresh = 0, box.palette="BuGn", extra = "auto", main = "Decision tree for treatment effects using development cohort")

prp(fit2, pal.thresh = 0, box.palette="BuGn", extra = "auto", main = "Strata:hba1cmonth > 6 (rounded)")

dev.off()


#:-------------------------------------------------------------------------------------
# Prediction of treatment effects for all patients

# patients which already have a calculated effect
patient_effects <- hba1c.train.complete.vs %>%
  select(patid, pated) %>%
  cbind(effects = colMeans(bcf_model$tau)) %>%
  rbind(
    hba1c.test.complete.vs %>%
      select(patid, pated) %>%
      cbind(effects = colMeans(predictions.hba1c.test$tau))
  )


# patients not yet calculated
full.cohort <- set_up_data_sglt2_glp1(dataset.type="full.cohort")

interim.dataset <- full.cohort[!(full.cohort$pated %in% patient_effects$pated),] %>%
  # left_join propensity scores%>%
  left_join(patient_prop_scores, by = c("patid", "pated")) %>%
  # select variables to make prediction
  select(patid, pated, drugclass, unique(c(variables_mu, variables_tau)), prop.score) %>%
  # this results in a lot of missing hba1cmonth, should I set these to 12 months?
  mutate(hba1cmonth = ifelse(is.na(hba1cmonth), 12, hba1cmonth)) %>%
  drop_na()


if (class(try(
  
  patient_effects <- readRDS(paste0(output_path, "/response_model_bcf/patient_effects.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  predictions.interim <- predict(object = bcf_model,
                                 x_predict_control = interim.dataset %>%
                                   select(
                                     all_of(variables_mu)
                                   ) %>%
                                   mutate_all(funs(as.numeric(.))) %>%
                                   as.matrix(),
                                 x_predict_moderate = interim.dataset %>%
                                   select(
                                     all_of(variables_tau)
                                   ) %>%
                                   mutate_all(funs(as.numeric(.))) %>%
                                   as.matrix(),
                                 pi_pred = 1-interim.dataset$prop.score,
                                 z_pred = interim.dataset %>%
                                   select(drugclass) %>%
                                   mutate(drugclass = ifelse(drugclass == "GLP1", 0, 1)) %>%
                                   unlist(),
                                 save_tree_directory = paste0(output_path, "/response_model_bcf/trees_no_prop"))
  
  
  
  patient_effects <- patient_effects %>%
    rbind(
      interim.dataset %>%
        select(patid, pated) %>%
        cbind(effects = colMeans(predictions.interim$tau))
    )
  
  
  saveRDS(patient_effects, paste0(output_path, "/response_model_bcf/patient_effects.rds"))
  
}


#:-------------------------------------------------------------------------------------
# Create table with patid, pated, predicted SGLT2 HbA1c, predicted GLP1 HbA1c, treatment effects difference

# patients which already have a calculated effect
patient_predicted_outcomes <- hba1c.train.complete.vs %>%
  select(patid, pated, drugclass) %>%
  cbind(effects.diff = colMeans(bcf_model$tau)) %>%
  cbind(pred = colMeans(bcf_model$yhat)) %>%
  # add SGLT2 prediction
  mutate(pred.SGLT2 = ifelse(drugclass == "SGLT2", pred, NA)) %>%
  # add GLP1 prediction 
  mutate(pred.GLP1 = ifelse(drugclass == "GLP1", pred, NA)) %>%
  # add SGLT2 predictions for those NA
  mutate(pred.SGLT2 = ifelse(is.na(pred.SGLT2), pred.GLP1+effects.diff, pred.SGLT2)) %>%
  # add GLP1 predictions for those NA
  mutate(pred.GLP1 = ifelse(is.na(pred.GLP1), pred.SGLT2-effects.diff, pred.GLP1)) %>%
  # remove pred, drugclass %>%
  select(-pred, -drugclass) %>%
  rbind(
    hba1c.test.complete.vs %>%
      select(patid, pated, drugclass) %>%
      cbind(effects.diff = colMeans(predictions.hba1c.test$tau)) %>%
      cbind(pred = colMeans(predictions.hba1c.test$yhat)) %>%
      # add SGLT2 prediction
      mutate(pred.SGLT2 = ifelse(drugclass == "SGLT2", pred, NA)) %>%
      # add GLP1 prediction 
      mutate(pred.GLP1 = ifelse(drugclass == "GLP1", pred, NA)) %>%
      # add SGLT2 predictions for those NA
      mutate(pred.SGLT2 = ifelse(is.na(pred.SGLT2), pred.GLP1+effects.diff, pred.SGLT2)) %>%
      # add GLP1 predictions for those NA
      mutate(pred.GLP1 = ifelse(is.na(pred.GLP1), pred.SGLT2-effects.diff, pred.GLP1)) %>%
      # remove pred, drugclass %>%
      select(-pred, -drugclass)
  )


# patients not yet calculated
full.cohort <- set_up_data_sglt2_glp1(dataset.type="full.cohort")

interim.dataset <- full.cohort[!(full.cohort$pated %in% patient_predicted_outcomes$pated),] %>%
  # left_join propensity scores%>%
  left_join(patient_prop_scores, by = c("patid", "pated")) %>%
  # select variables to make prediction
  select(patid, pated, drugclass, unique(c(variables_mu, variables_tau)), prop.score) %>%
  # this results in a lot of missing hba1cmonth, should I set these to 12 months?
  mutate(hba1cmonth = ifelse(is.na(hba1cmonth), 12, hba1cmonth)) %>%
  drop_na()


if (class(try(
  
  patient_predicted_outcomes <- readRDS(paste0(output_path, "/response_model_bcf/patient_predicted_outcomes.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  predictions.interim <- predict(object = bcf_model,
                                 x_predict_control = interim.dataset %>%
                                   select(
                                     all_of(variables_mu)
                                   ) %>%
                                   mutate_all(funs(as.numeric(.))) %>%
                                   as.matrix(),
                                 x_predict_moderate = interim.dataset %>%
                                   select(
                                     all_of(variables_tau)
                                   ) %>%
                                   mutate_all(funs(as.numeric(.))) %>%
                                   as.matrix(),
                                 pi_pred = 1-interim.dataset$prop.score,
                                 z_pred = interim.dataset %>%
                                   select(drugclass) %>%
                                   mutate(drugclass = ifelse(drugclass == "GLP1", 0, 1)) %>%
                                   unlist(),
                                 save_tree_directory = paste0(output_path, "/response_model_bcf/trees_no_prop"))
  
  
  
  patient_predicted_outcomes <- patient_predicted_outcomes %>%
    rbind(
      interim.dataset %>%
        select(patid, pated, drugclass) %>%
        cbind(effects.diff = colMeans(predictions.interim$tau)) %>%
        cbind(pred = colMeans(predictions.interim$yhat)) %>%
        # add SGLT2 prediction
        mutate(pred.SGLT2 = ifelse(drugclass == "SGLT2", pred, NA)) %>%
        # add GLP1 prediction 
        mutate(pred.GLP1 = ifelse(drugclass == "GLP1", pred, NA)) %>%
        # add SGLT2 predictions for those NA
        mutate(pred.SGLT2 = ifelse(is.na(pred.SGLT2), pred.GLP1+effects.diff, pred.SGLT2)) %>%
        # add GLP1 predictions for those NA
        mutate(pred.GLP1 = ifelse(is.na(pred.GLP1), pred.SGLT2-effects.diff, pred.GLP1)) %>%
        # remove pred, drugclass %>%
        select(-pred, -drugclass)
    )
  
  
  saveRDS(patient_predicted_outcomes, paste0(output_path, "/response_model_bcf/patient_predicted_outcomes.rds"))
  
}











