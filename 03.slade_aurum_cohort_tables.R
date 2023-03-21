####################
## Description:
##  - In this file we:
##    - Produce tables to describe cohorts
##    - Produce tables to describe forest plots.
####################


### Load libraries required
library(tidyverse)
library(tableone)

### Load functions required
source("11.02.slade_aurum_set_data.R")
source("11.01.slade_aurum_functions.R")

#:--------------------------------------------------------------------------------------
## Start of the file

### Run through the inclusion criteria to make all necessary datasets.

set_up_data_sglt2_glp1(dataset.type = "diagnostics")

### Make tables for each of the datasets

# Variables to be used in the table

vars <- c("yrdrugstart","agetx", "t2dmduration", "sex", "ethnicity", "drug_canagliflozin",
          "drug_dapagliflozin", "drug_empagliflozin", "drug_ertugliflozin", "drug_dulaglutide",
          "drug_exenatide_short", "drug_exenatide_long", "drug_liraglutide", "drug_lixisenatide",
          "deprivation", "smoke", "drugline", "ncurrtx", "MFN", "SU", "DPP4", "SGLT2", "TZD",
          "GLP1", "prehba1c", "prehba1c_na", "prebmi", "prebmi_na", "preegfr", "preegfr_na", 
          "prehdl", "prehdl_na", "prealt", "prealt_na", "prealbuminblood", "prealbuminblood_na",
          "prebilirubin", "prebilirubin_na", "pretotalcholesterol", "pretotalcholesterol_na",
          "premap", "premap_na", "prediabeticnephropathy", "preneuropathy", "preretinopathy",
          "preangina","CV_problems", "microvascular_complications", "ASCVD", "preaf", 
          "prerevasc", "preheartfailure", "prehypertension", "preihd", "premyocardialinfarction",
          "prepad", "prestroke","pretia", "preckd", "precld", "posthba1cfinal",
          "posthba1cfinal_na", "hba1cmonth", "hba1cmonth_na")

#:-----------------------
############## Predicted treatment benefits
####### > 5 mmol for either drug

### HbA1c train - dataset with complete variables used for fitting the model (smaller than hba1c.train)
hba1c.train.complete.benefit.simple <- set_up_data_sglt2_glp1(dataset.type = "hba1c.train") %>%
  left_join(readRDS("Samples/SGLT2-GLP1/Aurum//ps_model/patient_prop_scores.rds"), by = c("patid", "pated")) %>%
  # drop the variables with the most missingness
  select(-preacr, -preast, -prehaematocrit, -prehaemoglobin, -pretriglyceride) %>%
  # only complete cases
  drop_na() %>%
  as.data.frame() %>%
  # selected variables from SparseBCF
  select(patid, pated, posthba1cfinal, drugclass, unique(c(readRDS("Samples/SGLT2-GLP1/Aurum//response_model_bcf/variables_mu.rds"), readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/variables_tau.rds"))), prop.score) %>%
  # only complete cases
  drop_na() %>%
  as.data.frame() %>%
  left_join(readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/patient_effects.rds"), by = c("patid", "pated")) %>%
  mutate(benefit = ifelse(effects < -5, "SGLT2i", ifelse(effects > 5, "GLP1-RA", NA_real_))) %>%
  select(patid, pated, benefit)

hba1c.train.complete.benefit <- new.vars(hba1c.train.complete.benefit.simple)

tab.benefits <- CreateTableOne(vars = c("drugclass", vars), factorVars = c("MFN", "SU", "DPP4", "SGLT2", "TZD", "GLP1"), includeNA = TRUE, strata = "benefit", data = hba1c.train.complete.benefit, test = FALSE)
## Show table with SMD
print(tab.benefits, smd = TRUE)


### HbA1c test - dataset with complete variables used for validation
hba1c.test.complete.benefit.simple <- set_up_data_sglt2_glp1(dataset.type = "hba1c.test") %>%
  left_join(readRDS("Samples/SGLT2-GLP1/Aurum//ps_model/patient_prop_scores.rds"), by = c("patid", "pated")) %>%
  # selected variables from SparseBCF
  select(patid, pated, posthba1cfinal, drugclass, unique(c(readRDS("Samples/SGLT2-GLP1/Aurum//response_model_bcf/variables_mu.rds"), readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/variables_tau.rds"))), prop.score) %>%
  # only complete cases
  drop_na() %>%
  as.data.frame() %>%
  left_join(readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/patient_effects.rds"), by = c("patid", "pated")) %>%
  mutate(benefit = ifelse(effects < -5, "SGLT2i", ifelse(effects > 5, "GLP1-RA", NA_real_))) %>%
  select(patid, pated, benefit)

hba1c.test.complete.benefit <- new.vars(hba1c.test.complete.benefit.simple)

tab.benefits <- CreateTableOne(vars = c("drugclass", vars), factorVars = c("MFN", "SU", "DPP4", "SGLT2", "TZD", "GLP1"), includeNA = TRUE, strata = "benefit", data = hba1c.test.complete.benefit, test = FALSE)
## Show table with SMD
print(tab.benefits, smd = TRUE)


### HbA1c train
hba1c.train.benefit.simple <- set_up_data_sglt2_glp1(dataset.type = "hba1c.train") %>%
  left_join(readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/patient_effects.rds"), by = c("patid", "pated")) %>%
  mutate(benefit = ifelse(effects < -5, "SGLT2i", ifelse(effects > 5, "GLP1-RA", NA_real_))) %>%
  select(patid, pated, benefit)

hba1c.train.benefit <- new.vars(hba1c.train.benefit.simple)

# Abstract values (results)
tab.benefits <- CreateTableOne(vars = c("drugclass", vars), factorVars = c("MFN", "SU", "DPP4", "SGLT2", "TZD", "GLP1"), includeNA = TRUE, strata = "benefit", data = hba1c.train.benefit, test = FALSE)
## Show table with SMD
print(tab.benefits, smd = TRUE)


### Full cohort
full.cohort.benefit.simple <- set_up_data_sglt2_glp1(dataset.type = "full.cohort") %>%
  left_join(readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/patient_effects.rds"), by = c("patid", "pated")) %>%
  mutate(benefit = ifelse(effects < -5, "SGLT2i", ifelse(effects > 5, "GLP1-RA", NA_real_))) %>%
  select(patid, pated, benefit)

full.cohort.benefit <- new.vars(full.cohort.benefit.simple)

tab.benefits <- CreateTableOne(vars = c("drugclass", vars), factorVars = c("MFN", "SU", "DPP4", "SGLT2", "TZD", "GLP1"), includeNA = TRUE, strata = "benefit", data = full.cohort.benefit, test = FALSE)
## Show table with SMD
print(tab.benefits, smd = TRUE)


####### any benefit for either drug

### Full cohort
full.cohort.any_benefit.simple <- set_up_data_sglt2_glp1(dataset.type = "full.cohort") %>%
  left_join(readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/patient_effects.rds"), by = c("patid", "pated")) %>%
  mutate(benefit = ifelse(effects < 0, "SGLT2i", ifelse(effects > 0, "GLP1-RA", NA_real_))) %>%
  select(patid, pated, benefit)

full.cohort.any_benefit <- new.vars(full.cohort.any_benefit.simple)

tab.benefits <- CreateTableOne(vars = c("drugclass", vars), factorVars = c("MFN", "SU", "DPP4", "SGLT2", "TZD", "GLP1"), includeNA = TRUE, strata = "benefit", data = full.cohort.any_benefit, test = FALSE)
## Show table with SMD
print(tab.benefits, smd = TRUE)


table(full.cohort.any_benefit.simple$sex, full.cohort.any_benefit.simple$benefit)




#:-----------------------------------------------------------------------------
############## Generic table description

### Full cohort
full.cohort.simple <- set_up_data_sglt2_glp1(dataset.type = "full.cohort") %>%
  select(patid, pated)

full.cohort <- new.vars(full.cohort.simple)

## Construct a table
tab.full.cohort <- CreateTableOne(vars = vars, factorVars = c("MFN", "SU", "DPP4", "SGLT2", "TZD", "GLP1"), includeNA = TRUE, strata = "drugclass", data = full.cohort, test = FALSE)
## Show table with SMD
print(tab.full.cohort, smd = TRUE)

summary(full.cohort%>%select(agetx)) # results values
sd(full.cohort%>%select(agetx)%>%unlist(), na.rm = TRUE)
summary(full.cohort%>%select(sex))
table(full.cohort%>%select(sex))/nrow(full.cohort)*100
summary(full.cohort%>%select(ethnicity))
table(full.cohort%>%select(ethnicity))/nrow(full.cohort)*100

#:----------------------------
### PS model train
ps.model.train.simple <- set_up_data_sglt2_glp1(dataset.type = "ps.model.train") %>%
  select(patid, pated)

ps.model.train <- new.vars(ps.model.train.simple)

## Construct a table
tab.ps.model.train <- CreateTableOne(vars = vars, factorVars = c("MFN", "SU", "DPP4", "SGLT2", "TZD", "GLP1"), includeNA = TRUE, strata = "drugclass", data = ps.model.train, test = FALSE)
## Show table with SMD
print(tab.ps.model.train, smd = TRUE)


#:----------------------------
### PS model test
ps.model.test.simple <- set_up_data_sglt2_glp1(dataset.type = "ps.model.test") %>%
  select(patid, pated)

ps.model.test <- new.vars(ps.model.test.simple)

## Construct a table
tab.ps.model.test <- CreateTableOne(vars = vars, factorVars = c("MFN", "SU", "DPP4", "SGLT2", "TZD", "GLP1"), includeNA = TRUE, strata = "drugclass", data = ps.model.test, test = FALSE)
## Show table with SMD
print(tab.ps.model.test, smd = TRUE)


#:----------------------------
### HbA1c model train
hba1c.train.simple <- set_up_data_sglt2_glp1(dataset.type = "hba1c.train") %>%
  select(patid, pated)

hba1c.train <- new.vars(hba1c.train.simple)

## Construct a table
tab.hba1c.train <- CreateTableOne(vars = vars, factorVars = c("MFN", "SU", "DPP4", "SGLT2", "TZD", "GLP1"), includeNA = TRUE, strata = "drugclass", data = hba1c.train, test = FALSE)
## Show table with SMD
print(tab.hba1c.train, smd = TRUE)

# abstract values (methods)
summary(hba1c.train %>% select(prehba1c, agetx))
sd(hba1c.train%>%select(prehba1c)%>%unlist());sd(hba1c.train%>%select(agetx)%>%unlist())
table(hba1c.train %>% select(sex))/nrow(hba1c.train)*100


#:----------------------------
### HbA1c model test
hba1c.test.simple <- set_up_data_sglt2_glp1(dataset.type = "hba1c.test") %>%
  select(patid, pated)

hba1c.test <- new.vars(hba1c.test.simple)

## Construct a table
tab.hba1c.test <- CreateTableOne(vars = vars, factorVars = c("MFN", "SU", "DPP4", "SGLT2", "TZD", "GLP1"), includeNA = TRUE, strata = "drugclass", data = hba1c.test, test = FALSE)
## Show table with SMD
print(tab.hba1c.test, smd = TRUE)







#:-------------------------------------------------------------------------------------
#:-------------------------------------------------------------------------------------
#:-------------------------------------------------------------------------------------
#:-------------------------------------------------------------------------------------

output_path <- "Samples/SGLT2-GLP1/Aurum"

interval_breaks <- c(-5, -3, 0, 3, 5)

## Read in propensity scores
patient_prop_scores <- readRDS(paste0(output_path, "/ps_model/patient_prop_scores.rds"))

## Read in treatment effects
treatment_effects <- readRDS(paste0(output_path, "/response_model_bcf/patient_effects.rds"))

# load in variables used in the model
variables_mu <- readRDS(paste0(output_path, "/response_model_bcf/variables_mu.rds"))

variables_tau <- readRDS(paste0(output_path, "/response_model_bcf/variables_tau.rds"))

#### Tables containing the numerical results for clinical utilities
### HbA1c predictions of change

model_variables <- unique(c(variables_mu, variables_tau))[which(unique(c(variables_mu, variables_tau)) != "sex")]

## CPRD
## HbA1c clinical subgroups - predicted hba1c benefit
full.dataset <- set_up_data_sglt2_glp1(dataset.type = "full.cohort") %>%
  left_join(patient_prop_scores, by = c("patid", "pated")) %>%
  left_join(treatment_effects, by = c("patid", "pated")) %>%
  mutate(hba1c.change = posthba1cfinal - prehba1c,
         drugclass = factor(drugclass, levels = c("GLP1", "SGLT2")))

# still population into subgroups
group.full.dataset <- group_values(data = full.dataset,
                                   variable = "effects",
                                   breaks = interval_breaks) %>%
  drop_na(intervals) %>%
  rename("hba1c_diff" = "effects")



# predictions for the hba1c adjusted model
predictions_hba1c_stan_adjusted_overall <- readRDS(paste0(output_path, "/additional_outcomes/predictions_hba1c_stan_adjusted_overall.rds")) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci))

# predictions for the hba1c adjusted model full
predictions_hba1c_stan_adjusted_full <- readRDS(paste0(output_path, "/additional_outcomes/predictions_hba1c_stan_adjusted_full.rds")) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci))


print(predictions_hba1c_stan_adjusted_full, digits = 4)
print(predictions_hba1c_stan_adjusted_overall, digits = 2)


### HbA1c treatment effect
# combine everyone that can be included in the HbA1c analysis
hba1c <- set_up_data_sglt2_glp1(dataset.type = "hba1c.train") %>%
  left_join(patient_prop_scores, by = c("patid", "pated")) %>%
  left_join(treatment_effects, by = c("patid", "pated")) %>%
  rbind(
    set_up_data_sglt2_glp1(dataset.type = "hba1c.test") %>%
      left_join(patient_prop_scores, by = c("patid", "pated")) %>%
      left_join(treatment_effects, by = c("patid", "pated"))
  )

# split into subgroups according to hba1c benefit predicted
group.hba1c.dataset <- group_values(data = hba1c,
                                    variable = "effects",
                                    breaks = interval_breaks) %>%
  drop_na(intervals) %>%
  rename("hba1c_diff" = "effects")

# Full dataset

#-- PSM

ATE_psm_1_1_hba1c_full <- calc_ATE_validation_prop_matching(group.hba1c.dataset%>%mutate(intervals = as.numeric(1)), "posthba1cfinal", group.hba1c.dataset$prop.score, quantile_var = "intervals", adjust = FALSE, order = "largest")

ATE_psm_1_1_hba1c <- calc_ATE_validation_prop_matching(group.hba1c.dataset%>%mutate(intervals = as.numeric(intervals)), "posthba1cfinal", group.hba1c.dataset$prop.score, quantile_var = "intervals", adjust = FALSE, order = "largest")

print(ATE_psm_1_1_hba1c_full[["effects"]], digits =3)
print(ATE_psm_1_1_hba1c[["effects"]], digits =3) # abstract values (results)

#-- PSM + adjust

ATE_psm_1_1_adjust_hba1c_full <- calc_ATE_validation_prop_matching(group.hba1c.dataset%>%mutate(intervals = as.numeric(1)), "posthba1cfinal", group.hba1c.dataset$prop.score, quantile_var = "intervals", breakdown = model_variables, adjust = TRUE, order = "largest")

ATE_psm_1_1_adjust_hba1c <- calc_ATE_validation_prop_matching(group.hba1c.dataset%>%mutate(intervals = as.numeric(intervals)), "posthba1cfinal", group.hba1c.dataset$prop.score, quantile_var = "intervals", breakdown = model_variables, adjust = TRUE, order = "largest")

print(ATE_psm_1_1_adjust_hba1c_full[["effects"]], digits =3)
print(ATE_psm_1_1_adjust_hba1c[["effects"]], digits =3)

#-- Adjust

ATE_adjust_hba1c_full <- calc_ATE_validation_adjust(group.hba1c.dataset %>% mutate(intervals = as.numeric(1)), "posthba1cfinal", quantile_var = "intervals", breakdown = model_variables, adjust = TRUE)

ATE_adjust_hba1c <- calc_ATE_validation_adjust(group.hba1c.dataset, "posthba1cfinal", quantile_var = "intervals", breakdown = model_variables, adjust = TRUE)

print(ATE_adjust_hba1c_full[["effects"]], digits =3)
print(ATE_adjust_hba1c[["effects"]], digits =3)

# Strata by sex
group.hba1c.dataset.male <- group.hba1c.dataset %>% filter(sex == "Male")
group.hba1c.dataset.female <- group.hba1c.dataset %>% filter(sex == "Female")


#-- PSM

#------ Male

ATE_psm_1_1_hba1c_male_full <- calc_ATE_validation_prop_matching(group.hba1c.dataset.male%>%mutate(intervals = as.numeric(1)), "posthba1cfinal", group.hba1c.dataset.male$prop.score, quantile_var = "intervals", adjust = FALSE, order = "largest")

ATE_psm_1_1_hba1c_male <- calc_ATE_validation_prop_matching(group.hba1c.dataset.male%>%mutate(intervals = as.numeric(intervals)), "posthba1cfinal", group.hba1c.dataset.male$prop.score, quantile_var = "intervals", adjust = FALSE, order = "largest")

print(ATE_psm_1_1_hba1c_male_full[["effects"]], digits =3)
print(ATE_psm_1_1_hba1c_male[["effects"]], digits =3)

#------ Female

ATE_psm_1_1_hba1c_female_full <- calc_ATE_validation_prop_matching(group.hba1c.dataset.female%>%mutate(intervals = as.numeric(1)), "posthba1cfinal", group.hba1c.dataset.female$prop.score, quantile_var = "intervals", adjust = FALSE, order = "largest")

ATE_psm_1_1_hba1c_female <- calc_ATE_validation_prop_matching(group.hba1c.dataset.female%>%mutate(intervals = as.numeric(intervals)), "posthba1cfinal", group.hba1c.dataset.female$prop.score, quantile_var = "intervals", adjust = FALSE, order = "largest")

print(ATE_psm_1_1_hba1c_female_full[["effects"]], digits =3)
print(ATE_psm_1_1_hba1c_female[["effects"]], digits =3)

#-- PSM + adjust

#------ Male

ATE_psm_1_1_adjust_hba1c_male_full <- calc_ATE_validation_prop_matching(group.hba1c.dataset.male%>%mutate(intervals = as.numeric(1)), "posthba1cfinal", group.hba1c.dataset.male$prop.score, quantile_var = "intervals", breakdown = model_variables, adjust = TRUE, order = "largest")

ATE_psm_1_1_adjust_hba1c_male <- calc_ATE_validation_prop_matching(group.hba1c.dataset.male%>%mutate(intervals = as.numeric(intervals)), "posthba1cfinal", group.hba1c.dataset.male$prop.score, quantile_var = "intervals", breakdown = model_variables, adjust = TRUE, order = "largest")

print(ATE_psm_1_1_adjust_hba1c_male_full[["effects"]], digits =3)
print(ATE_psm_1_1_adjust_hba1c_male[["effects"]], digits =3)

#------ Female

ATE_psm_1_1_adjust_hba1c_female_full <- calc_ATE_validation_prop_matching(group.hba1c.dataset.female%>%mutate(intervals = as.numeric(1)), "posthba1cfinal", group.hba1c.dataset.female$prop.score, quantile_var = "intervals", breakdown = model_variables, adjust = TRUE, order = "largest")

ATE_psm_1_1_adjust_hba1c_female <- calc_ATE_validation_prop_matching(group.hba1c.dataset.female%>%mutate(intervals = as.numeric(intervals)), "posthba1cfinal", group.hba1c.dataset.female$prop.score, quantile_var = "intervals", breakdown = model_variables, adjust = TRUE, order = "largest")

print(ATE_psm_1_1_adjust_hba1c_female_full[["effects"]], digits =3)
print(ATE_psm_1_1_adjust_hba1c_female[["effects"]], digits =3)

#-- Adjust

#------ Male

ATE_adjust_hba1c_male_full <- calc_ATE_validation_adjust(group.hba1c.dataset.male%>%mutate(intervals = as.numeric(1)), "posthba1cfinal", quantile_var = "intervals", breakdown = model_variables, adjust = TRUE)

ATE_adjust_hba1c_male <- calc_ATE_validation_adjust(group.hba1c.dataset.male, "posthba1cfinal", quantile_var = "intervals", breakdown = model_variables, adjust = TRUE)

print(ATE_adjust_hba1c_male_full[["effects"]], digits =3)
print(ATE_adjust_hba1c_male[["effects"]], digits =3)

#------ Female

ATE_adjust_hba1c_female_full <- calc_ATE_validation_adjust(group.hba1c.dataset.female%>%mutate(intervals = as.numeric(1)), "posthba1cfinal", quantile_var = "intervals", breakdown = model_variables, adjust = TRUE)

ATE_adjust_hba1c_female <- calc_ATE_validation_adjust(group.hba1c.dataset.female, "posthba1cfinal", quantile_var = "intervals", breakdown = model_variables, adjust = TRUE)

print(ATE_adjust_hba1c_female_full[["effects"]], digits =3)
print(ATE_adjust_hba1c_female[["effects"]], digits =3)



## Scotland

# THIS IS STILL MISSING, NEEDS TO BE SOLVED


## America Vets dataset

# THIS IS STILL MISSING, NEEDS TO BE SOLVED


### Weight

## CPRD

# Full dataset

## Read in data for weight
weight.dataset <- set_up_data_sglt2_glp1(dataset.type = "weight.dataset") %>%
  left_join(patient_prop_scores, by = c("patid", "pated")) %>%
  left_join(treatment_effects, by = c("patid", "pated")) %>%
  mutate(w.change = postweight - preweight)


group.weight.dataset <- group_values(data = weight.dataset,
                                     variable = "effects",
                                     breaks = interval_breaks) %>%
  drop_na(intervals)

#-- PSM
predictions_weight_stan_psm_1_1_full <- readRDS(paste0(output_path, "/additional_outcomes/predictions_weight_stan_psm_1_1_full.rds"))
predictions_weight_stan_psm_1_1_overall <- readRDS(paste0(output_path, "/additional_outcomes/predictions_weight_stan_psm_1_1_overall.rds"))


print(predictions_weight_stan_psm_1_1_full%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

print(predictions_weight_stan_psm_1_1_overall%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

#-- PSM + adjust

predictions_weight_stan_psm_1_1_adjusted_full <- readRDS(paste0(output_path, "/additional_outcomes/predictions_weight_stan_psm_1_1_adjusted_full.rds"))

predictions_weight_stan_psm_1_1_adjusted_overall <- readRDS(paste0(output_path, "/additional_outcomes/predictions_weight_stan_psm_1_1_adjusted_overall.rds"))


print(predictions_weight_stan_psm_1_1_adjusted_full%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

print(predictions_weight_stan_psm_1_1_adjusted_overall%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

#-- Adjust

predictions_weight_stan_adjusted_full <- readRDS(paste0(output_path, "/additional_outcomes/predictions_weight_stan_adjusted_full.rds"))

predictions_weight_stan_adjusted_overall <- readRDS(paste0(output_path, "/additional_outcomes/predictions_weight_stan_adjusted_overall.rds"))


print(predictions_weight_stan_adjusted_full%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

print(predictions_weight_stan_adjusted_overall%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

# Strata by sex

#-- PSM

#------ Male

predictions_weight_stan_psm_1_1_full

predictions_weight_stan_psm_1_1 <- readRDS(paste0(output_path, "/additional_outcomes/predictions_weight_stan_psm_1_1.rds"))


print(predictions_weight_stan_psm_1_1_full%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

print(predictions_weight_stan_psm_1_1%>%filter(sex=="Male")%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

#------ Female

print(predictions_weight_stan_psm_1_1_full%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

print(predictions_weight_stan_psm_1_1%>%filter(sex=="Female")%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

#-- PSM + adjust

#------ Male

predictions_weight_stan_psm_1_1_adjusted <- readRDS(paste0(output_path, "/additional_outcomes/predictions_weight_stan_psm_1_1_adjusted.rds"))


print(predictions_weight_stan_psm_1_1_adjusted_full%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

print(predictions_weight_stan_psm_1_1_adjusted%>%filter(sex=="Male")%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

#------ Female

print(predictions_weight_stan_psm_1_1_adjusted_full%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

print(predictions_weight_stan_psm_1_1_adjusted%>%filter(sex=="Female")%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

#-- Adjust

#------ Male

predictions_weight_stan_adjusted <- readRDS(paste0(output_path, "/additional_outcomes/predictions_weight_stan_adjusted.rds"))


print(predictions_weight_stan_adjusted_full%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

print(predictions_weight_stan_adjusted%>%filter(sex=="Male")%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

#------ Female

print(predictions_weight_stan_adjusted_full%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

print(predictions_weight_stan_adjusted%>%filter(sex=="Female")%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)


### eGFR

## CPRD

# Full dataset

#-- PSM
predictions_egfr_stan_psm_1_1_full <- readRDS(paste0(output_path, "/additional_outcomes/predictions_egfr_stan_psm_1_1_full.rds"))
predictions_egfr_stan_psm_1_1_overall <- readRDS(paste0(output_path, "/additional_outcomes/predictions_egfr_stan_psm_1_1_overall.rds"))


print(predictions_egfr_stan_psm_1_1_full%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

print(predictions_egfr_stan_psm_1_1_overall%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

#-- PSM + adjust

predictions_egfr_stan_psm_1_1_adjusted_full <- readRDS(paste0(output_path, "/additional_outcomes/predictions_egfr_stan_psm_1_1_adjusted_full.rds"))

predictions_egfr_stan_psm_1_1_adjusted_overall <- readRDS(paste0(output_path, "/additional_outcomes/predictions_egfr_stan_psm_1_1_adjusted_overall.rds"))


print(predictions_egfr_stan_psm_1_1_adjusted_full%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

print(predictions_egfr_stan_psm_1_1_adjusted_overall%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

#-- Adjust

predictions_egfr_stan_adjusted_full <- readRDS(paste0(output_path, "/additional_outcomes/predictions_egfr_stan_adjusted_full.rds"))

predictions_egfr_stan_adjusted_overall <- readRDS(paste0(output_path, "/additional_outcomes/predictions_egfr_stan_adjusted_overall.rds"))


print(predictions_egfr_stan_adjusted_full%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

print(predictions_egfr_stan_adjusted_overall%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

# Strata by sex

#-- PSM

#------ Male

predictions_egfr_stan_psm_1_1_full

predictions_egfr_stan_psm_1_1 <- readRDS(paste0(output_path, "/additional_outcomes/predictions_egfr_stan_psm_1_1.rds"))


print(predictions_egfr_stan_psm_1_1_full%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

print(predictions_egfr_stan_psm_1_1%>%filter(sex=="Male")%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

#------ Female

print(predictions_egfr_stan_psm_1_1_full%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

print(predictions_egfr_stan_psm_1_1%>%filter(sex=="Female")%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

#-- PSM + adjust

#------ Male

predictions_egfr_stan_psm_1_1_adjusted <- readRDS(paste0(output_path, "/additional_outcomes/predictions_egfr_stan_psm_1_1_adjusted.rds"))


print(predictions_egfr_stan_psm_1_1_adjusted_full%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

print(predictions_egfr_stan_psm_1_1_adjusted%>%filter(sex=="Male")%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

#------ Female

print(predictions_egfr_stan_psm_1_1_adjusted_full%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

print(predictions_egfr_stan_psm_1_1_adjusted%>%filter(sex=="Female")%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

#-- Adjust

#------ Male

predictions_egfr_stan_adjusted <- readRDS(paste0(output_path, "/additional_outcomes/predictions_egfr_stan_adjusted.rds"))


print(predictions_egfr_stan_adjusted_full%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

print(predictions_egfr_stan_adjusted%>%filter(sex=="Male")%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

#------ Female

print(predictions_egfr_stan_adjusted_full%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

print(predictions_egfr_stan_adjusted%>%filter(sex=="Female")%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)



### Discontinuation

## CPRD

discontinuation.dataset <- set_up_data_sglt2_glp1(dataset.type = "discontinuation.dataset") %>%
  left_join(patient_prop_scores, by = c("patid", "pated")) %>%
  left_join(treatment_effects, by = c("patid", "pated")) %>%
  mutate(stopdrug_6m_3mFU = factor(stopdrug_6m_3mFU))

group.discontinuation.dataset <- group_values(data = discontinuation.dataset,
                                              variable = "effects",
                                              breaks = interval_breaks) %>%
  drop_na(intervals)

# Full dataset

#-- PSM
predictions_discontinuation_stan_psm_1_1_full <- readRDS(paste0(output_path, "/additional_outcomes/predictions_discontinuation_stan_psm_1_1_full.rds"))
predictions_discontinuation_stan_psm_1_1_overall <- readRDS(paste0(output_path, "/additional_outcomes/predictions_discontinuation_stan_psm_1_1_overall.rds"))


print(predictions_discontinuation_stan_psm_1_1_full%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

print(predictions_discontinuation_stan_psm_1_1_overall%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

#-- PSM + adjust

predictions_discontinuation_stan_psm_1_1_adjusted_full <- readRDS(paste0(output_path, "/additional_outcomes/predictions_discontinuation_stan_psm_1_1_adjusted_full.rds"))

predictions_discontinuation_stan_psm_1_1_adjusted_overall <- readRDS(paste0(output_path, "/additional_outcomes/predictions_discontinuation_stan_psm_1_1_adjusted_overall.rds"))


print(predictions_discontinuation_stan_psm_1_1_adjusted_full%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

print(predictions_discontinuation_stan_psm_1_1_adjusted_overall%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

#-- Adjust

predictions_discontinuation_stan_adjusted_full <- readRDS(paste0(output_path, "/additional_outcomes/predictions_discontinuation_stan_adjusted_full.rds"))

predictions_discontinuation_stan_adjusted_overall <- readRDS(paste0(output_path, "/additional_outcomes/predictions_discontinuation_stan_adjusted_overall.rds"))


print(predictions_discontinuation_stan_adjusted_full%>%mutate(mean = as.numeric(mean)*100, lci = as.numeric(lci)*100, uci = as.numeric(uci)*100), digits = 3)

print(predictions_discontinuation_stan_adjusted_overall%>%mutate(mean = as.numeric(mean)*100, lci = as.numeric(lci)*100, uci = as.numeric(uci)*100), digits = 3)

# Strata by sex

#-- PSM

#------ Male

predictions_discontinuation_stan_psm_1_1_full

predictions_discontinuation_stan_psm_1_1 <- readRDS(paste0(output_path, "/additional_outcomes/predictions_discontinuation_stan_psm_1_1.rds"))


print(predictions_discontinuation_stan_psm_1_1_full%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

print(predictions_discontinuation_stan_psm_1_1%>%filter(sex=="Male")%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

#------ Female

print(predictions_discontinuation_stan_psm_1_1_full%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

print(predictions_discontinuation_stan_psm_1_1%>%filter(sex=="Female")%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

#-- PSM + adjust

#------ Male

predictions_discontinuation_stan_psm_1_1_adjusted <- readRDS(paste0(output_path, "/additional_outcomes/predictions_discontinuation_stan_psm_1_1_adjusted.rds"))


print(predictions_discontinuation_stan_psm_1_1_adjusted_full%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

print(predictions_discontinuation_stan_psm_1_1_adjusted%>%filter(sex=="Male")%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

#------ Female

print(predictions_discontinuation_stan_psm_1_1_adjusted_full%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

print(predictions_discontinuation_stan_psm_1_1_adjusted%>%filter(sex=="Female")%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

#-- Adjust

#------ Male

predictions_discontinuation_stan_adjusted <- readRDS(paste0(output_path, "/additional_outcomes/predictions_discontinuation_stan_adjusted.rds"))


print(predictions_discontinuation_stan_adjusted_full%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

print(predictions_discontinuation_stan_adjusted%>%filter(sex=="Male")%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

#------ Female

print(predictions_discontinuation_stan_adjusted_full%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

print(predictions_discontinuation_stan_adjusted%>%filter(sex=="Female")%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)



### CVD

## CPRD

# Full dataset

patient_prop_scores_qrisk <- readRDS(paste0(output_path, "/additional_outcomes/patient_prop_scores_qrisk.rds"))

no_co.dataset <- set_up_data_sglt2_glp1(dataset.type="no_co.dataset") %>%
  left_join(patient_prop_scores_qrisk, by = c("patid", "pated")) %>%
  left_join(treatment_effects, by = c("patid", "pated"))

group.no_co.dataset <- group_values(data = no_co.dataset,
                                    variable = "effects",
                                    breaks = interval_breaks) %>%
  drop_na(intervals)

#-- PSM

predictions_no_co_cvd_stan_psm_1_1_full <- readRDS(paste0(output_path, "/additional_outcomes/predictions_no_co_cvd_stan_psm_1_1_full.rds"))
predictions_no_co_cvd_stan_psm_1_1_overall <- readRDS(paste0(output_path, "/additional_outcomes/predictions_no_co_cvd_stan_psm_1_1_overall.rds"))


print(predictions_no_co_cvd_stan_psm_1_1_full%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

print(predictions_no_co_cvd_stan_psm_1_1_overall%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

#-- PSM + adjust

predictions_no_co_cvd_stan_psm_1_1_adjusted_full <- readRDS(paste0(output_path, "/additional_outcomes/predictions_no_co_cvd_stan_psm_1_1_adjusted_full.rds"))

predictions_no_co_cvd_stan_psm_1_1_adjusted_overall <- readRDS(paste0(output_path, "/additional_outcomes/predictions_no_co_cvd_stan_psm_1_1_adjusted_overall.rds"))


print(predictions_no_co_cvd_stan_psm_1_1_adjusted_full%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

print(predictions_no_co_cvd_stan_psm_1_1_adjusted_overall%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

#-- Adjust

predictions_no_co_cvd_stan_adjusted_full <- readRDS(paste0(output_path, "/additional_outcomes/predictions_no_co_cvd_stan_adjusted_full.rds"))

predictions_no_co_cvd_stan_adjusted_overall <- readRDS(paste0(output_path, "/additional_outcomes/predictions_no_co_cvd_stan_adjusted_overall.rds"))


print(predictions_no_co_cvd_stan_adjusted_full%>%mutate(mean = exp(as.numeric(mean)), lci = exp(as.numeric(lci)), uci = exp(as.numeric(uci))), digits = 3)

print(predictions_no_co_cvd_stan_adjusted_overall%>%mutate(mean = exp(as.numeric(mean)), lci = exp(as.numeric(lci)), uci = exp(as.numeric(uci))), digits = 3)


# Strata by sex

#-- PSM

#------ Male

predictions_no_co_cvd_stan_psm_1_1_full

predictions_no_co_cvd_stan_psm_1_1 <- readRDS(paste0(output_path, "/additional_outcomes/predictions_no_co_cvd_stan_psm_1_1.rds"))


print(predictions_no_co_cvd_stan_psm_1_1_full%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

print(predictions_no_co_cvd_stan_psm_1_1%>%filter(sex=="Male")%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

#------ Female

print(predictions_no_co_cvd_stan_psm_1_1_full%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

print(predictions_no_co_cvd_stan_psm_1_1%>%filter(sex=="Female")%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

#-- PSM + adjust

#------ Male

predictions_no_co_cvd_stan_psm_1_1_adjusted <- readRDS(paste0(output_path, "/additional_outcomes/predictions_no_co_cvd_stan_psm_1_1_adjusted.rds"))


print(predictions_no_co_cvd_stan_psm_1_1_adjusted_full%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

print(predictions_no_co_cvd_stan_psm_1_1_adjusted%>%filter(sex=="Male")%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

#------ Female

print(predictions_no_co_cvd_stan_psm_1_1_adjusted_full%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

print(predictions_no_co_cvd_stan_psm_1_1_adjusted%>%filter(sex=="Female")%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

#-- Adjust

#------ Male

predictions_no_co_cvd_stan_adjusted <- readRDS(paste0(output_path, "/additional_outcomes/predictions_no_co_cvd_stan_adjusted.rds"))


print(predictions_no_co_cvd_stan_adjusted_full%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

print(predictions_no_co_cvd_stan_adjusted%>%filter(sex=="Male")%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

#------ Female

print(predictions_no_co_cvd_stan_adjusted_full%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

print(predictions_no_co_cvd_stan_adjusted%>%filter(sex=="Female")%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

### Microvascular complications

## CPRD 

# Full dataset

patient_prop_scores_qrisk <- readRDS(paste0(output_path, "/additional_outcomes/patient_prop_scores_qrisk.rds"))

micro_comp.dataset <- set_up_data_sglt2_glp1(dataset.type="micro_comp.dataset") %>%
  left_join(patient_prop_scores_qrisk, by = c("patid", "pated")) %>%
  left_join(treatment_effects, by = c("patid", "pated"))

group.micro_comp.dataset <- group_values(data = micro_comp.dataset,
                                         variable = "effects",
                                         breaks = interval_breaks) %>%
  drop_na(intervals)


# predictions for the Microvascular complications outcomes in the population with no Microvascular complications
predictions_micro_comp_micro_comp_stan_adjusted_overall <- readRDS(paste0(output_path, "/additional_outcomes/predictions_micro_comp_micro_comp_stan_adjusted_overall.rds"))

# predictions for the Microvascular complications outcomes in the population with no Microvascular complications
predictions_micro_comp_micro_comp_stan_adjusted_full <- readRDS(paste0(output_path, "/additional_outcomes/predictions_micro_comp_micro_comp_stan_adjusted_full.rds"))

print(predictions_micro_comp_micro_comp_stan_adjusted_full%>%mutate(mean = exp(as.numeric(mean)), lci = exp(as.numeric(lci)), uci = exp(as.numeric(uci))), digits = 2)

print(predictions_micro_comp_micro_comp_stan_adjusted_overall%>%mutate(mean = exp(as.numeric(mean)), lci = exp(as.numeric(lci)), uci = exp(as.numeric(uci))), digits = 2)


### HF

## CPRD

# Full dataset

#-- PSM

predictions_no_co_hf_stan_psm_1_1_full <- readRDS(paste0(output_path, "/additional_outcomes/predictions_no_co_hf_stan_psm_1_1_full.rds"))
predictions_no_co_hf_stan_psm_1_1_overall <- readRDS(paste0(output_path, "/additional_outcomes/predictions_no_co_hf_stan_psm_1_1_overall.rds"))


print(predictions_no_co_hf_stan_psm_1_1_full%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

print(predictions_no_co_hf_stan_psm_1_1_overall%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

#-- PSM + adjust

predictions_no_co_hf_stan_psm_1_1_adjusted_full <- readRDS(paste0(output_path, "/additional_outcomes/predictions_no_co_hf_stan_psm_1_1_adjusted_full.rds"))

predictions_no_co_hf_stan_psm_1_1_adjusted_overall <- readRDS(paste0(output_path, "/additional_outcomes/predictions_no_co_hf_stan_psm_1_1_adjusted_overall.rds"))


print(predictions_no_co_hf_stan_psm_1_1_adjusted_full%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

print(predictions_no_co_hf_stan_psm_1_1_adjusted_overall%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

#-- Adjust

predictions_no_co_hf_stan_adjusted_full <- readRDS(paste0(output_path, "/additional_outcomes/predictions_no_co_hf_stan_adjusted_full.rds"))

predictions_no_co_hf_stan_adjusted_overall <- readRDS(paste0(output_path, "/additional_outcomes/predictions_no_co_hf_stan_adjusted_overall.rds"))


print(predictions_no_co_hf_stan_adjusted_full%>%mutate(mean = exp(as.numeric(mean)), lci = exp(as.numeric(lci)), uci = exp(as.numeric(uci))), digits = 2)

print(predictions_no_co_hf_stan_adjusted_overall%>%mutate(mean = exp(as.numeric(mean)), lci = exp(as.numeric(lci)), uci = exp(as.numeric(uci))), digits = 2)


# Strata by sex

#-- PSM

#------ Male

predictions_no_co_hf_stan_psm_1_1_full

predictions_no_co_hf_stan_psm_1_1 <- readRDS(paste0(output_path, "/additional_outcomes/predictions_no_co_hf_stan_psm_1_1.rds"))


print(predictions_no_co_hf_stan_psm_1_1_full%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

print(predictions_no_co_hf_stan_psm_1_1%>%filter(sex=="Male")%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

#------ Female

print(predictions_no_co_hf_stan_psm_1_1_full%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

print(predictions_no_co_hf_stan_psm_1_1%>%filter(sex=="Female")%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

#-- PSM + adjust

#------ Male

predictions_no_co_hf_stan_psm_1_1_adjusted <- readRDS(paste0(output_path, "/additional_outcomes/predictions_no_co_hf_stan_psm_1_1_adjusted.rds"))


print(predictions_no_co_hf_stan_psm_1_1_adjusted_full%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

print(predictions_no_co_hf_stan_psm_1_1_adjusted%>%filter(sex=="Male")%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

#------ Female

print(predictions_no_co_hf_stan_psm_1_1_adjusted_full%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

print(predictions_no_co_hf_stan_psm_1_1_adjusted%>%filter(sex=="Female")%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

#-- Adjust

#------ Male

predictions_no_co_hf_stan_adjusted <- readRDS(paste0(output_path, "/additional_outcomes/predictions_no_co_hf_stan_adjusted.rds"))


print(predictions_no_co_hf_stan_adjusted_full%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

print(predictions_no_co_hf_stan_adjusted%>%filter(sex=="Male")%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

#------ Female

print(predictions_no_co_hf_stan_adjusted_full%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)

print(predictions_no_co_hf_stan_adjusted%>%filter(sex=="Female")%>%mutate(mean = as.numeric(mean), lci = as.numeric(lci), uci = as.numeric(uci)), digits = 3)



### CKD

#-- Adjust

# predictions for the CKD - stage 3/4/5 outcomes in the population with no CVD/HF/CKD
predictions_no_co_egfr40_or_ckd5_stan_adjusted_overall <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_no_co_egfr40_or_ckd5_stan_adjusted_overall.rds")

# predictions for the CKD - stage 3/4/5 outcomes in the population with no CVD/HF/CKD
predictions_no_co_egfr40_or_ckd5_stan_adjusted_full <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_no_co_egfr40_or_ckd5_stan_adjusted_full.rds")


print(predictions_no_co_egfr40_or_ckd5_stan_adjusted_full%>%mutate(mean = exp(as.numeric(mean)), lci = exp(as.numeric(lci)), uci = exp(as.numeric(uci))), digits = 2)

print(predictions_no_co_egfr40_or_ckd5_stan_adjusted_overall%>%mutate(mean = exp(as.numeric(mean)), lci = exp(as.numeric(lci)), uci = exp(as.numeric(uci))), digits = 2)


