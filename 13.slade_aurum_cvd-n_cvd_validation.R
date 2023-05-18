####################
## Description:
##  - In this file we validate the model in those with and without CVD.
####################


## Load libraries
library(tidyverse)
library(patchwork)

## Set up directory path to save files (stagered to ensure folders are created)

dir.create("Samples")
dir.create("Samples/SGLT2-GLP1")

output_path <- "Samples/SGLT2-GLP1/Aurum"
dir.create(output_path)

## make directory for outputs
dir.create("Plots")



###############################################################################
###############################################################################
############################### Read Data In ##################################
###############################################################################
###############################################################################

## Load functions required

source("01.slade_aurum_functions.R")
source("02.slade_aurum_set_data.R")

###############################################################################
###############################################################################
########################## General variables ##################################
###############################################################################
###############################################################################

# load in variables used in the model
variables_mu <- readRDS(paste0(output_path, "/response_model_bcf/variables_mu.rds"))

variables_tau <- readRDS(paste0(output_path, "/response_model_bcf/variables_tau.rds"))



variables_mu <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/variables_mu.rds")

variables_tau <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/variables_tau.rds")

hba1c.test <- set_up_data_sglt2_glp1(dataset.type = "hba1c.test") %>%
  left_join(readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/patient_effects.rds"), by = c("patid", "pated")) %>%
  left_join(readRDS("Samples/SGLT2-GLP1/Aurum/ps_model/patient_prop_scores.rds"), by = c("patid", "pated"))

hba1c.test.with.cvd <- hba1c.test %>%
  filter(preangina == "Yes" | preaf == "Yes" | prerevasc == "Yes" | preheartfailure == "Yes" |preihd == "Yes" | premyocardialinfarction == "Yes" | prestroke == "Yes" | pretia == "Yes") %>%
  mutate(hba1c_diff.q = ntile(effects, 10)) %>%
  rename("hba1c_diff" = "effects") %>%
  drop_na(hba1c_diff.q)

ATE_psm_1_1_hba1c.test.with.cvd <- calc_ATE(data = hba1c.test.with.cvd,
                                            validation_type = "PSM", variable = "posthba1cfinal",
                                            quantile_var = "hba1c_diff.q", prop_scores = hba1c.test.with.cvd$prop.score, 
                                            order = "largest", breakdown = unique(c(variables_tau, variables_mu)))

ATE_psm_1_1_adjusted_hba1c.test.with.cvd <- calc_ATE(data = hba1c.test.with.cvd,
                                                     validation_type = "PSM + adjust", variable = "posthba1cfinal",
                                                     quantile_var = "hba1c_diff.q", prop_scores = hba1c.test.with.cvd$prop.score, 
                                                     order = "largest", breakdown = unique(c(variables_tau, variables_mu)))

ATE_adjusted_hba1c.test.with.cvd <- calc_ATE(data = hba1c.test.with.cvd,
                                             validation_type = "Adjust", variable = "posthba1cfinal",
                                             quantile_var = "hba1c_diff.q",
                                             order = "largest", breakdown = unique(c(variables_tau, variables_mu)))


plot_ATE_psm_1_1_hba1c.test.with.cvd <- ATE_plot(ATE_psm_1_1_hba1c.test.with.cvd[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -12, 12) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(
    title = paste0("Cohort with CVD (n=", format(sum(ATE_psm_1_1_hba1c.test.with.cvd[["effects"]]$n_drug1)*2,big.mark=",",scientific=FALSE), ")"),
    subtitle = "Propensity score matching"
  )

plot_ATE_psm_1_1_adjusted_hba1c.test.with.cvd <- ATE_plot(ATE_psm_1_1_adjusted_hba1c.test.with.cvd[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -12, 12) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(
    title = paste0("Cohort with CVD (n=", format(sum(ATE_psm_1_1_adjusted_hba1c.test.with.cvd[["effects"]]$n_drug1)*2,big.mark=",",scientific=FALSE), ")"),
    subtitle = "Propensity score matching and estimate adjustment"
  )

plot_ATE_adjusted_hba1c.test.with.cvd <- ATE_plot(ATE_adjusted_hba1c.test.with.cvd[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -12, 12) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(
    title = paste0("Cohort with CVD (n=", format(sum(ATE_adjusted_hba1c.test.with.cvd[["effects"]]$N),big.mark=",",scientific=FALSE), ")"),
    subtitle = "Estimate adjustment"
  )


hba1c.test.without.cvd <- hba1c.test %>%
  filter(!pated %in% hba1c.test.with.cvd$pated) %>%
  mutate(hba1c_diff.q = ntile(effects, 10)) %>%
  rename("hba1c_diff" = "effects") %>%
  drop_na(hba1c_diff.q)

ATE_psm_1_1_hba1c.test.without.cvd <- calc_ATE(data = hba1c.test.without.cvd,
                                               validation_type = "PSM", variable = "posthba1cfinal",
                                               quantile_var = "hba1c_diff.q", prop_scores = hba1c.test.without.cvd$prop.score, 
                                               order = "largest", breakdown = unique(c(variables_tau, variables_mu)))

ATE_psm_1_1_adjusted_hba1c.test.without.cvd <- calc_ATE(data = hba1c.test.without.cvd,
                                                        validation_type = "PSM + adjust", variable = "posthba1cfinal",
                                                        quantile_var = "hba1c_diff.q", prop_scores = hba1c.test.without.cvd$prop.score, 
                                                        order = "largest", breakdown = unique(c(variables_tau, variables_mu)))

ATE_adjusted_hba1c.test.without.cvd <- calc_ATE(data = hba1c.test.without.cvd,
                                                validation_type = "Adjust", variable = "posthba1cfinal",
                                                quantile_var = "hba1c_diff.q",
                                                order = "largest", breakdown = unique(c(variables_tau, variables_mu)))


plot_ATE_psm_1_1_hba1c.test.without.cvd <- ATE_plot(ATE_psm_1_1_hba1c.test.without.cvd[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -12, 12) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(
    title = paste0("Cohort without CVD (n=", format(sum(ATE_psm_1_1_hba1c.test.without.cvd[["effects"]]$n_drug1)*2,big.mark=",",scientific=FALSE), ")"),
    subtitle = "Propensity score matching"
  )

plot_ATE_psm_1_1_adjusted_hba1c.test.without.cvd <- ATE_plot(ATE_psm_1_1_adjusted_hba1c.test.without.cvd[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -12, 12) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(
    title = paste0("Cohort without CVD (n=", format(sum(ATE_psm_1_1_adjusted_hba1c.test.without.cvd[["effects"]]$n_drug1)*2,big.mark=",",scientific=FALSE), ")"),
    subtitle = "Propensity score matching and estimate adjustment"
  )

plot_ATE_adjusted_hba1c.test.without.cvd <- ATE_plot(ATE_adjusted_hba1c.test.without.cvd[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -12, 12) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(
    title = paste0("Cohort without CVD (n=", format(sum(ATE_adjusted_hba1c.test.without.cvd[["effects"]]$N),big.mark=",",scientific=FALSE), ")"),
    subtitle = "Estimate adjustment"
  )

plot_all <- (plot_ATE_psm_1_1_hba1c.test.with.cvd | plot_ATE_psm_1_1_hba1c.test.without.cvd) / (plot_ATE_psm_1_1_adjusted_hba1c.test.with.cvd | plot_ATE_psm_1_1_adjusted_hba1c.test.without.cvd) / (plot_ATE_adjusted_hba1c.test.with.cvd | plot_ATE_adjusted_hba1c.test.without.cvd) +
  plot_annotation(tag_levels = list(c("A.1", "A.2", "B.1", "B.2", "C.1", "C.2")))


pdf(width = 10, height = 15, "Plots/11.13.cvd_no_cvd_validation.pdf")
plot_all
dev.off()

