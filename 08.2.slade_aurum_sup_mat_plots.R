####################
## Description:
##  - In this file we create the plots used in the supplementary material of the paper
####################


## Load libraries
library(tidyverse)
library(patchwork)
library(scales)
library(forestplot)
library(rms)

## Load functions required
source("11.01.slade_aurum_functions.R")
source("11.02.slade_aurum_set_data.R")

## make directory for outputs
dir.create("Plots/Paper")

## make directory for outputs
dir.create("Plots/Paper/Sup_Mat")


#### Supplementary manuscript plots / stuff


## Flowchart 1

# Finding the number of individuals discarded for the training/validation dataset used in the model

# training dataset
hba1c.train <- set_up_data_sglt2_glp1(dataset.type = "hba1c.train")

hba1c.train.complete <- set_up_data_sglt2_glp1(dataset.type = "hba1c.train") %>%
  # drop the variables with the most missingness (>40%)
  select(-preacr, -preast, -prehaematocrit, -prehaemoglobin, -pretriglyceride) %>%
  # only complete cases
  drop_na() %>%
  as.data.frame()

# nrow(hba1c.train)
# nrow(hba1c.train.complete)

# validation dataset

hba1c.test <- set_up_data_sglt2_glp1(dataset.type = "hba1c.test")

hba1c.test.complete <- set_up_data_sglt2_glp1(dataset.type = "hba1c.test") %>%
  # selected variables from SparseBCF
  select(patid, pated, posthba1cfinal, drugclass, unique(c(readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/variables_mu.rds"), readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/variables_tau.rds")))) %>%
  # only complete cases
  drop_na() %>%
  as.data.frame()

# nrow(hba1c.test)
# nrow(hba1c.test.complete)



## Flowchart 2

# Finding the number of individuals discarded for the weight dataset (incomplete vars)

# weight dataset
weight.dataset <- set_up_data_sglt2_glp1(dataset.type = "weight.dataset")

weight.dataset.effects <- set_up_data_sglt2_glp1(dataset.type = "weight.dataset") %>%
  left_join(readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/patient_effects.rds"), by = c("patid", "pated")) %>%
  drop_na(effects)

# nrow(weight.dataset)
# nrow(weight.dataset.effects)

# no comorbidities
no_co.dataset <- set_up_data_sglt2_glp1(dataset.type = "no_co.dataset")

no_co.dataset.effects <- set_up_data_sglt2_glp1(dataset.type = "no_co.dataset") %>%
  left_join(readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/patient_effects.rds"), by = c("patid", "pated")) %>%
  drop_na(effects)

# nrow(no_co.dataset)
# nrow(no_co.dataset.effects)

# discontinuation
discontinuation.dataset <- set_up_data_sglt2_glp1(dataset.type = "discontinuation.dataset")

discontinuation.dataset.effects <- set_up_data_sglt2_glp1(dataset.type = "discontinuation.dataset") %>%
  left_join(readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/patient_effects.rds"), by = c("patid", "pated")) %>%
  drop_na(effects)

# nrow(discontinuation.dataset)
# nrow(discontinuation.dataset.effects)

# microvascular complications
micro_comp.dataset <- set_up_data_sglt2_glp1(dataset.type = "micro_comp.dataset")

micro_comp.dataset.effects <- set_up_data_sglt2_glp1(dataset.type = "micro_comp.dataset") %>%
  left_join(readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/patient_effects.rds"), by = c("patid", "pated")) %>%
  drop_na(effects)

# nrow(micro_comp.dataset)
# nrow(micro_comp.dataset.effects)

#:--------------------------------------------------------------------------------
### Variable selection from the Sparse BCF model


## Variable selection of sparcebcf for the moderator effect
variables_tau <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/assessment/variables_tau_original_chain_1.rds")

plot_variables_tau <- variables_tau %>%
  rename("Current age" = "agetx",
         "Sex" = "sex",
         "Duration of diabetes" = "t2dmduration",
         "Number of glucose-lowering drug classes ever prescribed" = "drugline",
         "Number of other current glucose-lowering drugs" = "ncurrtx",
         "Month of Hba1C measure" = "hba1cmonth",
         "HbA1c" = "prehba1c",
         "BMI" = "prebmi",
         "eGFR" = "preegfr",
         "Albumin" = "prealbuminblood",
         "Alanine transaminase" = "prealt",
         "Bilirubin" = "prebilirubin",
         "HDL" = "prehdl",
         "Mean arterial blood pressure" = "premap",
         "Total cholesterol" = "pretotalcholesterol",
         "Angina" = "preangina",
         "Chronic liver disease" = "precld",
         "Nephropathy" = "prediabeticnephropathy",
         "Heart failure" = "preheartfailure",
         "Hypertension" = "prehypertension",
         "Ischaemic heart disease" = "preihd",
         "Myocardial infarction" = "premyocardialinfarction",
         "Neuropathy" = "preneuropathy",
         "Peripheral arterial disease" = "prepad",
         "Retinopathy" = "preretinopathy",
         "Atherosclerotic cardiovascular disease" = "prerevasc",
         "Stroke" = "prestroke",
         "Transient ischaemic attack" = "pretia",
         "Atrial fibrillation" = "preaf") %>%
  as.data.frame() %>%
  gather(key, value) %>%
  arrange(desc(value)) %>%
  mutate(key = factor(key),
         colour = ifelse(value > 0.023, "Above", "Below"),
         value = value*100) %>%
  ggplot(aes(y = forcats::fct_reorder(key, value), x = value, colour = colour)) +
  geom_vline(aes(xintercept = 0.023 * 100), colour = "red") +
  geom_segment(aes(x = 0, xend = value, yend = forcats::fct_reorder(key, value)), linetype = "dashed") +
  geom_point(size = 2) +
  ggtitle("Moderator component") +
  xlab("Posterior inclusion proportions (%)") +
  scale_colour_manual(values = c("Above" = "black", "Below" = "grey")) +
  theme_light() +
  theme(axis.text.y = element_text(angle = 30, face = "bold", colour = c(rep("grey70", 18), rep("black", 12))),
        axis.title.y = element_blank(),
        legend.position = "none")


## Variable selection of sparcebcf for the prognostic effect
variables_mu <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/assessment/variables_mu_original_chain_1.rds")

plot_variables_mu <- variables_mu %>%
  rename("Current age" = "agetx",
         "Sex" = "sex",
         "Duration of diabetes" = "t2dmduration",
         "Number of glucose-lowering drug classes ever prescribed" = "drugline",
         "Number of other current glucose-lowering drugs" = "ncurrtx",
         "Month of Hba1C measure" = "hba1cmonth",
         "HbA1c" = "prehba1c",
         "BMI" = "prebmi",
         "eGFR" = "preegfr",
         "Albumin" = "prealbuminblood",
         "Alanine transaminase" = "prealt",
         "Bilirubin" = "prebilirubin",
         "HDL" = "prehdl",
         "Mean arterial blood pressure" = "premap",
         "Total cholesterol" = "pretotalcholesterol",
         "Angina" = "preangina",
         "Chronic liver disease" = "precld",
         "Nephropathy" = "prediabeticnephropathy",
         "Heart failure" = "preheartfailure",
         "Hypertension" = "prehypertension",
         "Ischaemic heart disease" = "preihd",
         "Myocardial infarction" = "premyocardialinfarction",
         "Neuropathy" = "preneuropathy",
         "Peripheral arterial disease" = "prepad",
         "Retinopathy" = "preretinopathy",
         "Atherosclerotic cardiovascular disease" = "prerevasc",
         "Stroke" = "prestroke",
         "Transient ischaemic attack" = "pretia",
         "Atrial fibrillation" = "preaf",
         "Propensity score" = "propensity score") %>%
  as.data.frame() %>%
  gather(key, value) %>%
  arrange(desc(value)) %>%
  mutate(key = factor(key),
         colour = ifelse(value > 0.023, "Above", "Below"),
         value = value*100) %>%
  ggplot(aes(y = forcats::fct_reorder(key, value), x = value, colour = colour)) +
  geom_vline(aes(xintercept = 0.023 * 100), colour = "red") +
  geom_segment(aes(x = 0, xend = value, yend = forcats::fct_reorder(key, value)), linetype = "dashed") +
  geom_point(size = 2) +
  ggtitle("Prognostic component") +
  xlab("Posterior inclusion proportions (%)") +
  scale_colour_manual(values = c("Above" = "black", "Below" = "grey")) +
  theme_light() +
  theme(axis.text.y = element_text(angle = 30, face = "bold", colour = c(rep("grey70", 21), rep("black", 9))),
        axis.title.y = element_blank(),
        legend.position = "none")

# combining plots
plot_sup_1 <- patchwork::wrap_plots(list(plot_variables_mu, plot_variables_tau), ncol = 2) +
  patchwork::plot_annotation(tag_levels = c('A'),
                             title = "BART model variable selection") # title of full plot

# pdf of plot
pdf(width = 14, height = 6, "Plots/Paper/Sup_Mat/11.08.plot_sup_1.pdf")
plot_sup_1
dev.off()


#:--------------------------------------------------------------------------------
### Sex difference in trials

## glp1 order: Liraglutide, albi, priba, CPRD, scotland

## sglt2 order: trials, cprd, scotland 

# create dataset
df <- cbind(
  dataset = c("CPRD (n=11,296)", "PRIBA (n=550)", "Harmony 7 RCT Liraglutide (n=389)", "Harmony RCT (pooled) Albiglutide (n=1,682)", "CPRD (n=40,915)", "Trials meta-analysis (n=7,119)", "Scottish dataset (n=1,837)", "Scottish dataset (n=415)"),
  drug = c("GLP-1 receptor agonists", "GLP-1 receptor agonists", "GLP-1 receptor agonists", "GLP-1 receptor agonists", "SGLT2-inhibitors", "SGLT2-inhibitors", "SGLT2-inhibitors", "GLP-1 receptor agonists"),
  n = c(11296, 550, 389, 1682, 40915, 7119, 1837, 415),
  beta = c(2.0257, 2.972881, 4.403, 1.029, -0.9015, -1.05, -2.22353, 0.65486),
  lci = c(1.4409514, 0.4652423, 2.2, 0.14, -1.1514, -1.63, -3.3612905, -2.4401361),
  uci = c(2.61047723, 5.480521, 6.3, 1.9, -0.6516, -0.47, -1.0857777, 3.7498463)
) %>%
  as.data.frame() %>%
  mutate(n = as.numeric(n),
         beta = as.numeric(beta),
         lci = as.numeric(lci),
         uci = as.numeric(uci))

# make plot
plot_sup_2.1 <- df %>%
  filter(drug == "SGLT2-inhibitors") %>%
  mutate(dataset = factor(dataset, levels = c("Scottish dataset (n=1,837)", "CPRD (n=40,915)", "Trials meta-analysis (n=7,119)"))) %>%
  ggplot(aes(x = beta, y = dataset)) +
  geom_vline(aes(xintercept = 0), colour = "black") +
  geom_point(colour = "#f1a340") +
  geom_errorbar(aes(xmin = lci, xmax = uci), colour = "#f1a340", width = 0.3) +
  ggtitle("SGLT2-inhibitors") +
  xlab("HbA1c response difference in males versus females (mmol/mol)") +
  xlim(min(df%>%select(beta, lci, uci)%>%unlist()), max(df%>%select(beta, lci, uci)%>%unlist())) +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_line(colour = "white"))

plot_sup_2.2 <- df %>%
  filter(drug == "GLP-1 receptor agonists") %>%
  mutate(dataset = factor(dataset, levels = c("Scottish dataset (n=415)", "CPRD (n=11,296)", "PRIBA (n=550)", "Harmony RCT (pooled) Albiglutide (n=1,682)", "Harmony 7 RCT Liraglutide (n=389)"))) %>%
  ggplot(aes(x = beta, y = dataset)) +
  geom_vline(aes(xintercept = 0), colour = "black") +
  geom_point(colour = "dodgerblue2") +
  geom_errorbar(aes(xmin = lci, xmax = uci), colour = "dodgerblue2", width = 0.3) +
  ggtitle("GLP-1 receptor agonist") +
  xlab("HbA1c response difference in males versus females (mmol/mol)") +
  xlim(min(df%>%select(beta, lci, uci)%>%unlist()), max(df%>%select(beta, lci, uci)%>%unlist())) +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_line(colour = "white"))

plot_sup_2 <- patchwork::wrap_plots(list(plot_sup_2.1, plot_sup_2.2), ncol = 1, nrow = 2) +
  plot_annotation(tag_levels = c("A"))


# pdf of plot
pdf(width = 8, height = 6, "Plots/Paper/Sup_Mat/11.08.plot_sup_2.pdf")
plot_sup_2
dev.off()





#:--------------------------------------------------------------------------------
### Figure of treatment effects validation for Ethnicity

interval_breaks <- c(-5, -3, 0, 3, 5)

variables_mu <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/variables_mu.rds")

variables_tau <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/variables_tau.rds")

# Create dataset with all the variables needed
ethnicity.dataset <- set_up_data_sglt2_glp1(dataset.type = "full.cohort") %>%
  left_join(readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/patient_effects.rds"), by = c("patid", "pated")) %>%
  select(all_of(c("patid", "pated", "ethnicity", "drugclass", "posthba1cfinal", unique(c(variables_mu, variables_tau)), "effects"))) %>%
  drop_na() %>%
  mutate(ethnicity = factor(ethnicity, levels = c("White", "South Asian", "Black", "Other", "Mixed"), labels = c("White", "South Asian", "Black", "Other/Mixed", "Other/Mixed")))

# Create variable corresponding to the clinical subgroups
ethnicity.dataset.grouped <- group_values(data = ethnicity.dataset,
                                          variable = "effects",
                                          breaks = interval_breaks) %>%
  drop_na(intervals) %>%
  rename("hba1c_diff" = "effects")

# Subsetting dataset ethnicity - White
predicted_observed_white <- ethnicity.dataset.grouped %>%
  filter(ethnicity == "White")

# Subsetting dataset ethnicity - South Asian
predicted_observed_asian <- ethnicity.dataset.grouped %>%
  filter(ethnicity == "South Asian")

# Subsetting dataset ethnicity - Black
predicted_observed_black <- ethnicity.dataset.grouped %>%
  filter(ethnicity == "Black")

# Subsetting dataset ethnicity - Other/Mixed
predicted_observed_mixed <- ethnicity.dataset.grouped %>%
  filter(ethnicity == "Other/Mixed")



# White
ATE_adjust_white <- calc_ATE(data = predicted_observed_white%>%mutate(intervals=as.numeric(intervals)), validation_type = "Adjust",
                             variable = "posthba1cfinal", quantile_var = "intervals",
                             order = "largest", breakdown = unique(c(variables_mu, variables_tau)))

# South Asian
ATE_adjust_asian <- calc_ATE(data = predicted_observed_asian%>%mutate(intervals=as.numeric(intervals)), validation_type = "Adjust",
                             variable = "posthba1cfinal", quantile_var = "intervals",
                             order = "largest", breakdown = unique(c(variables_mu, variables_tau)))

# Black
ATE_adjust_black <- calc_ATE(data = predicted_observed_black%>%mutate(intervals=as.numeric(intervals)), validation_type = "Adjust",
                             variable = "posthba1cfinal", quantile_var = "intervals",
                             order = "largest", breakdown = unique(c(variables_mu, variables_tau)))

# Mixed / Other
ATE_adjust_mixed <- calc_ATE(data = predicted_observed_mixed%>%mutate(intervals=as.numeric(intervals)), validation_type = "Adjust",
                             variable = "posthba1cfinal", quantile_var = "intervals",
                             order = "largest", breakdown = unique(c(variables_mu, variables_tau)))


### Overall analysis (1 group)
# White
ATE_adjust_white_overall <- calc_ATE(data = predicted_observed_white%>%mutate(intervals=as.numeric(1)), validation_type = "Adjust",
                                     variable = "posthba1cfinal", quantile_var = "intervals",
                                     order = "largest", breakdown = unique(c(variables_mu, variables_tau)))

# South Asian
ATE_adjust_asian_overall <- calc_ATE(data = predicted_observed_asian%>%mutate(intervals=as.numeric(1)), validation_type = "Adjust",
                                     variable = "posthba1cfinal", quantile_var = "intervals",
                                     order = "largest", breakdown = unique(c(variables_mu, variables_tau)))

# Black
ATE_adjust_black_overall <- calc_ATE(data = predicted_observed_black%>%mutate(intervals=as.numeric(1)), validation_type = "Adjust",
                                     variable = "posthba1cfinal", quantile_var = "intervals",
                                     order = "largest", breakdown = unique(c(variables_mu, variables_tau)))

# Mixed / Other
ATE_adjust_mixed_overall <- calc_ATE(data = predicted_observed_mixed%>%mutate(intervals=as.numeric(1)), validation_type = "Adjust",
                                     variable = "posthba1cfinal", quantile_var = "intervals",
                                     order = "largest", breakdown = unique(c(variables_mu, variables_tau)))

# Set up axis limits for the forest plots
hba1c_strata_axis_min <- plyr::round_any(floor(min(c(
  ATE_adjust_white[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
  ATE_adjust_asian[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
  ATE_adjust_black[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
  ATE_adjust_mixed[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
  ATE_adjust_white_overall[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
  ATE_adjust_asian_overall[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
  ATE_adjust_black_overall[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
  ATE_adjust_mixed_overall[["effects"]] %>% select(c("obs","lci","uci")) %>% min()
))), 2, f = floor)

hba1c_strata_axis_max <- plyr::round_any(ceiling(max(c(
  ATE_adjust_white[["effects"]] %>% select(c("obs","lci","uci")) %>% max(),
  ATE_adjust_asian[["effects"]] %>% select(c("obs","lci","uci")) %>% max(),
  ATE_adjust_black[["effects"]] %>% select(c("obs","lci","uci")) %>% max(),
  ATE_adjust_mixed[["effects"]] %>% select(c("obs","lci","uci")) %>% max(),
  ATE_adjust_white_overall[["effects"]] %>% select(c("obs","lci","uci")) %>% max(),
  ATE_adjust_asian_overall[["effects"]] %>% select(c("obs","lci","uci")) %>% max(),
  ATE_adjust_black_overall[["effects"]] %>% select(c("obs","lci","uci")) %>% max(),
  ATE_adjust_mixed_overall[["effects"]] %>% select(c("obs","lci","uci")) %>% max()
))), 2, f = ceiling)


# plot of forest plot
plot_sup_3 <- rbind(
  cbind(obs = 50, lci = 50, uci = 50, group = "White", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(obs = 50, lci = 50, uci = 50, group = "South Asian", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(obs = 50, lci = 50, uci = 50, group = "Black", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(obs = 50, lci = 50, uci = 50, group = "Other/Mixed", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(ATE_adjust_white[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(1:3), group = "White"),
  cbind(ATE_adjust_asian[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(1:3), group = "South Asian"),
  cbind(ATE_adjust_black[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(1:3), group = "Black"),
  cbind(ATE_adjust_mixed[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(1:3), group = "Other/Mixed"),
  cbind(obs = 50, lci = 50, uci = 50, group = "White", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(obs = 50, lci = 50, uci = 50, group = "South Asian", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(obs = 50, lci = 50, uci = 50, group = "Black", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(obs = 50, lci = 50, uci = 50, group = "Other/Mixed", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(ATE_adjust_white[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(4:6), group = "White"),
  cbind(ATE_adjust_asian[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(4:6), group = "South Asian"),
  cbind(ATE_adjust_black[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(4:6), group = "Black"),
  cbind(ATE_adjust_mixed[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(4:6), group = "Other/Mixed"),
  cbind(ATE_adjust_white_overall[["effects"]] %>% select(obs, lci, uci), group = "White", intervals = "Average treatment effect"),
  cbind(ATE_adjust_asian_overall[["effects"]] %>% select(obs, lci, uci), group = "South Asian", intervals = "Average treatment effect"),
  cbind(ATE_adjust_black_overall[["effects"]] %>% select(obs, lci, uci), group = "Black", intervals = "Average treatment effect"),
  cbind(ATE_adjust_mixed_overall[["effects"]] %>% select(obs, lci, uci), group = "Other/Mixed", intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(obs = as.numeric(obs),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == 1, ">5 mmol/mol",
                            ifelse(intervals == 2, "3-5 mmol/mol", 
                                   ifelse(intervals == 3, "0-3 mmol/mol",
                                          ifelse(intervals == 4, "0-3 mmol/mol", 
                                                 ifelse(intervals == 5, "3-5 mmol/mol",
                                                        ifelse(intervals == 6, ">5 mmol/mol", intervals)))))),
         group = factor(group)) %>%
  rename("lower" = "lci", "upper" = "uci", "mean" = "obs", "labeltext" = "intervals") %>%
  mutate(n_white = c(rep(NA, 4), rep(format(ATE_adjust_white[["effects"]]$N[1:3],big.mark=",",scientific=FALSE), 4),
                     rep(NA, 4), rep(format(ATE_adjust_white[["effects"]]$N[4:6],big.mark=",",scientific=FALSE), 4),
                     rep(format(sum(ATE_adjust_white[["effects"]]$N),big.mark=",",scientific=FALSE), 4)),
         n_asian = c(rep(NA, 4), rep(format(ATE_adjust_asian[["effects"]]$N[1:3],big.mark=",",scientific=FALSE), 4),
                     rep(NA, 4), rep(format(ATE_adjust_asian[["effects"]]$N[4:6],big.mark=",",scientific=FALSE), 4),
                     rep(format(sum(ATE_adjust_asian[["effects"]]$N),big.mark=",",scientific=FALSE), 4)),
         n_black = c(rep(NA, 4), rep(format(ATE_adjust_black[["effects"]]$N[1:3],big.mark=",",scientific=FALSE), 4),
                     rep(NA, 4), rep(format(ATE_adjust_black[["effects"]]$N[4:6],big.mark=",",scientific=FALSE), 4),
                     rep(format(sum(ATE_adjust_black[["effects"]]$N),big.mark=",",scientific=FALSE), 4)),
         n_mixed = c(rep(NA, 4), rep(format(ATE_adjust_mixed[["effects"]]$N[1:3],big.mark=",",scientific=FALSE), 4),
                     rep(NA, 4), rep(format(ATE_adjust_mixed[["effects"]]$N[4:6],big.mark=",",scientific=FALSE), 4),
                     rep(format(sum(ATE_adjust_mixed[["effects"]]$N),big.mark=",",scientific=FALSE), 4)),
         mean = as.numeric(mean)) %>%
  group_by(group) %>%
  forestplot(labeltext = c("labeltext", "n_white", "n_asian", "n_black", "n_mixed"),
             ci.vertices = TRUE,
             ci.vertices.height = 0.1,
             title = "Adjusted",
             clip = c(hba1c_strata_axis_min, hba1c_strata_axis_max),
             xticks = seq(hba1c_strata_axis_min, hba1c_strata_axis_max, 2),
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Average treatment effect (mmol/mol)") %>%
  fp_add_header(labeltext = paste0("Overall population (n=", format(ethnicity.dataset.grouped%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                n_white = "White",
                n_asian = "South Asian",
                n_black = "Black",
                n_mixed = "Other/Mixed") %>%
  fp_set_style(box = c("#E64B35B2","#4DBBD5B2", "#00A087B2", "#3C5488B2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE)) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = "black")


pdf(width = 15, height = 6, "Plots/Paper/Sup_Mat/11.08.plot_sup_3.pdf")
plot_sup_3
dev.off()


#:--------------------------------------------------------------------------------
### Figure of risk of developing CKD (40% drop of Stage 5)

interval_breaks <- c(-5, -3, 0, 3, 5)

patient_prop_scores_qrisk <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/patient_prop_scores_qrisk.rds")

no_co.dataset <- set_up_data_sglt2_glp1(dataset.type="no_co.dataset") %>%
  left_join(patient_prop_scores_qrisk, by = c("patid", "pated")) %>%
  left_join(readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/patient_effects.rds"), by = c("patid", "pated"))

group.no_co.dataset <- group_values(data = no_co.dataset,
                                    variable = "effects",
                                    breaks = interval_breaks) %>%
  drop_na(intervals)

breakdown_adjust <- unique(c(variables_mu, variables_tau))
# categorical variables in breakdown
factors <- sapply(group.no_co.dataset[,breakdown_adjust], is.factor)

matching_no_co <- MatchIt::matchit(
  formula = formula(paste0("drugclass ~  agetx + t2dmduration + prehba1c + preegfr + prealt + drugline + ncurrtx + sex + preneuropathy + preretinopathy")),
  data = group.no_co.dataset,
  method = "nearest",
  distance = group.no_co.dataset[,"prop.score"],
  replace = FALSE,
  m.order = "largest",
  caliper = 0.05,
  mahvars = NULL, estimand = "ATT", exact = NULL, antiexact = NULL, discard = "none", reestimate = FALSE, s.weights = NULL, std.caliper = TRUE, ratio = 1, verbose = FALSE, include.obj = FALSE,
)

group.no_co.dataset.matched <- group.no_co.dataset %>%
  slice(which(matching_no_co$weights == 1))

# predictions for the CKD - stage 3/4/5 outcomes in the population with no CVD/HF/CKD
predictions_no_co_egfr40_or_ckd5_stan_psm_1_1_overall <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_no_co_egfr40_or_ckd5_stan_psm_1_1_overall.rds")

# predictions for the CKD - stage 3/4/5 outcomes in the population with no CVD/HF/CKD
predictions_no_co_egfr40_or_ckd5_stan_psm_1_1_full <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_no_co_egfr40_or_ckd5_stan_psm_1_1_full.rds")

# predictions for the CKD - stage 3/4/5 outcomes in the population with no CVD/HF/CKD
predictions_no_co_egfr40_or_ckd5_stan_psm_1_1_adjusted_overall <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_no_co_egfr40_or_ckd5_stan_psm_1_1_adjusted_overall.rds")

# predictions for the CKD - stage 3/4/5 outcomes in the population with no CVD/HF/CKD
predictions_no_co_egfr40_or_ckd5_stan_psm_1_1_adjusted_full <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_no_co_egfr40_or_ckd5_stan_psm_1_1_adjusted_full.rds")

# predictions for the CKD - stage 3/4/5 outcomes in the population with no CVD/HF/CKD
predictions_no_co_egfr40_or_ckd5_stan_adjusted_overall <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_no_co_egfr40_or_ckd5_stan_adjusted_overall.rds")

# predictions for the CKD - stage 3/4/5 outcomes in the population with no CVD/HF/CKD
predictions_no_co_egfr40_or_ckd5_stan_adjusted_full <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_no_co_egfr40_or_ckd5_stan_adjusted_full.rds")

# set axis limits
axis_min <- exp(min(c(predictions_no_co_egfr40_or_ckd5_stan_psm_1_1_overall%>%select(-intervals)%>%unlist()%>%as.numeric(),
                  predictions_no_co_egfr40_or_ckd5_stan_psm_1_1_full%>%unlist()%>%as.numeric(),
                  predictions_no_co_egfr40_or_ckd5_stan_psm_1_1_adjusted_overall%>%select(-intervals)%>%unlist()%>%as.numeric(),
                  predictions_no_co_egfr40_or_ckd5_stan_psm_1_1_adjusted_full%>%unlist()%>%as.numeric(),
                  predictions_no_co_egfr40_or_ckd5_stan_adjusted_overall%>%select(-intervals)%>%unlist()%>%as.numeric(),
                  predictions_no_co_egfr40_or_ckd5_stan_adjusted_full%>%unlist()%>%as.numeric())))

axis_max <- exp(max(c(predictions_no_co_egfr40_or_ckd5_stan_psm_1_1_overall%>%select(-intervals)%>%unlist()%>%as.numeric(),
                      predictions_no_co_egfr40_or_ckd5_stan_psm_1_1_full%>%unlist()%>%as.numeric(),
                      predictions_no_co_egfr40_or_ckd5_stan_psm_1_1_adjusted_overall%>%select(-intervals)%>%unlist()%>%as.numeric(),
                      predictions_no_co_egfr40_or_ckd5_stan_psm_1_1_adjusted_full%>%unlist()%>%as.numeric(),
                      predictions_no_co_egfr40_or_ckd5_stan_adjusted_overall%>%select(-intervals)%>%unlist()%>%as.numeric(),
                      predictions_no_co_egfr40_or_ckd5_stan_adjusted_full%>%unlist()%>%as.numeric())))

# plot forest plot
## PSM
plot_psm <- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_no_co_egfr40_or_ckd5_stan_psm_1_1_overall %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = NA, lci = NA, uci = NA),
  predictions_no_co_egfr40_or_ckd5_stan_psm_1_1_overall %>% slice(4:6),
  predictions_no_co_egfr40_or_ckd5_stan_psm_1_1_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[1])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[1])%>%filter(postdrug_egfr40_or_ckd5_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                            ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[2])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[2])%>%filter(postdrug_egfr40_or_ckd5_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                   ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[3])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[3])%>%filter(postdrug_egfr40_or_ckd5_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                          ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[4])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[4])%>%filter(postdrug_egfr40_or_ckd5_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                 ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[5])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[5])%>%filter(postdrug_egfr40_or_ckd5_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                        ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[6])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[6])%>%filter(postdrug_egfr40_or_ckd5_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci") %>%
  mutate(mean = exp(mean),
         lower = exp(lower),
         upper = exp(upper)) %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             zero = 1,
             title = "Propensity score matching individuals",
             xlog = TRUE,
             clip = c(axis_min, axis_max),
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazards ratio (95% CI, log scale)") %>%
  fp_add_header(paste0("Overall population (n=", format(nrow(group.no_co.dataset.matched),big.mark=",",scientific=FALSE), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = "black")


## PSM + adjusted
plot_psm_adusted <- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_no_co_egfr40_or_ckd5_stan_psm_1_1_adjusted_overall %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = NA, lci = NA, uci = NA),
  predictions_no_co_egfr40_or_ckd5_stan_psm_1_1_adjusted_overall %>% slice(4:6),
  predictions_no_co_egfr40_or_ckd5_stan_psm_1_1_adjusted_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[1])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[1])%>%filter(postdrug_egfr40_or_ckd5_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                            ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[2])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[2])%>%filter(postdrug_egfr40_or_ckd5_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                   ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[3])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[3])%>%filter(postdrug_egfr40_or_ckd5_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                          ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[4])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[4])%>%filter(postdrug_egfr40_or_ckd5_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                 ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[5])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[5])%>%filter(postdrug_egfr40_or_ckd5_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                        ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[6])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[6])%>%filter(postdrug_egfr40_or_ckd5_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci") %>%
  mutate(mean = exp(mean),
         lower = exp(lower),
         upper = exp(upper)) %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             zero = 1,
             title = "Propensity score matching individuals and adjustment of estimates",
             xlog = TRUE,
             clip = c(axis_min, axis_max),
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazards ratio (95% CI, log scale)") %>%
  fp_add_header(paste0("Overall population (n=", format(nrow(group.no_co.dataset.matched),big.mark=",",scientific=FALSE), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = "black")

## Adjusted
plot_adusted <- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_no_co_egfr40_or_ckd5_stan_adjusted_overall %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = NA, lci = NA, uci = NA),
  predictions_no_co_egfr40_or_ckd5_stan_adjusted_overall %>% slice(4:6),
  predictions_no_co_egfr40_or_ckd5_stan_adjusted_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.no_co.dataset$intervals)[1], paste0(">5 mmol/mol (n=", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[1])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[1])%>%filter(postdrug_egfr40_or_ckd5_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                            ifelse(intervals == levels(group.no_co.dataset$intervals)[2], paste0("3-5 mmol/mol (n=", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[2])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[2])%>%filter(postdrug_egfr40_or_ckd5_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                   ifelse(intervals == levels(group.no_co.dataset$intervals)[3], paste0("0-3 mmol/mol (n=", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[3])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[3])%>%filter(postdrug_egfr40_or_ckd5_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                          ifelse(intervals == levels(group.no_co.dataset$intervals)[4], paste0("0-3 mmol/mol (n=", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[4])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[4])%>%filter(postdrug_egfr40_or_ckd5_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                 ifelse(intervals == levels(group.no_co.dataset$intervals)[5], paste0("3-5 mmol/mol (n=", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[5])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[5])%>%filter(postdrug_egfr40_or_ckd5_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                        ifelse(intervals == levels(group.no_co.dataset$intervals)[6], paste0(">5 mmol/mol (n=", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[6])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[6])%>%filter(postdrug_egfr40_or_ckd5_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci") %>%
  mutate(mean = exp(mean),
         lower = exp(lower),
         upper = exp(upper)) %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             zero = 1,
             title = "Adjustment of estimates",
             xlog = TRUE,
             clip = c(axis_min, axis_max),
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazards ratio (95% CI, log scale)") %>%
  fp_add_header(paste0("Overall population (n=", format(nrow(group.no_co.dataset),big.mark=",",scientific=FALSE), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = "black")



pdf(width = 7, height = 12, "Plots/Paper/Sup_Mat/11.08.plot_sup_4.pdf")

grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 3,
                                           ncol = 1, heights = unit(c(5, 5, 5), "null"))))

# first plot
pushViewport(viewport(layout.pos.row = 1,
                      layout.pos.col = 1))
plot_psm
upViewport()
# second plot
pushViewport(viewport(layout.pos.row = 2,
                      layout.pos.col = 1))
plot_psm_adusted
upViewport()
# third plot
pushViewport(viewport(layout.pos.row = 3,
                      layout.pos.col = 1))
plot_adusted
upViewport()

grid.text(expression(bold("A")), vp = viewport(layout.pos.row = 1, layout.pos.col = 1), hjust = 28.5, vjust = -15)
grid.text(expression(bold("B")), vp = viewport(layout.pos.row = 2, layout.pos.col = 1), hjust = 28.5, vjust = -15)
grid.text(expression(bold("C")), vp = viewport(layout.pos.row = 3, layout.pos.col = 1), hjust = 28.5, vjust = -15)

dev.off()


#:--------------------------------------------------------------------------------
### Compare predictions from BCF SGLT2 vs John et al. SGLT2

# read in predictions from BCF model
patient_predicted_outcomes <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/patient_predicted_outcomes.rds")

# load the object for the Linear model. This model object is the object used
#   in the SGLT2vsDPP4 model by John Dennis et al.
load("m1_hba1cmodel_SGLT2_DPP4.Rdata")

# Collect the full cohort from CPRD Aurum that can be used to fit the SGLT2vsDPP4 linear model
full.cohort.updated <- set_up_data_sglt2_glp1(dataset.type = "full.cohort") %>%
  select(patid, pated, prehba1c, drugclass, drugline, ncurrtx, preegfr, prealt, prebmi, agetx, hba1cmonth) %>%
  rename("prealtlog"="prealt",
         "prehba1cmmol"="prehba1c",
         "egfr_ckdepi"="preegfr") %>%
  mutate(prealtlog = log(prealtlog),
         ncurrtx = factor(ncurrtx, levels = c("1","2","3","4","5+"), labels = c("0","1","2","3","3")),
         drugline = factor(drugline, levels = c("2","3","4","5+"), labels = c("2","3","4","5")),
         hba1cmonth = ifelse(is.na(hba1cmonth), 12, hba1cmonth)) %>%
  drop_na() %>%
  # apply rules from the model
  filter(prehba1cmmol < 120) %>%
  filter(prehba1cmmol >= 53) %>%
  filter(egfr_ckdepi > 45)


# predict for SGLT2
dataset.sglt2 <- full.cohort.updated %>%
  select(-patid, -pated) %>%
  mutate(drugclass = factor("SGLT2", levels = c("SGLT2", "DPP4")))
predictions.sglt2 <- predict(m1, dataset.sglt2)



## Compare SGLT2 BCF vs Linear regression
interim.dataset <- full.cohort.updated %>%
  select(patid, pated) %>%
  cbind(pred.SGLT2.lm = predictions.sglt2) %>%
  left_join(patient_predicted_outcomes %>%
              select(patid, pated, pred.SGLT2) %>%
              rename("pred.SGLT2.bcf" = "pred.SGLT2"), by = c("patid", "pated")) %>%
  select(-patid, -pated)

# Plot the comparison of SGLT2 predictions from the linear SGLT2vsDPP4 model vs SGLT2vsGLP1 model
plot_sup_5.1 <- interim.dataset %>%
  ggplot(aes(y = pred.SGLT2.bcf, x = pred.SGLT2.lm)) +
  geom_point() +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red", linetype = "dashed") +
  stat_smooth() +
  ylab("SGLT2i predictions with BCF") +
  xlab("SGLT2i predictions with Linear Regression") +
  ylim(min(interim.dataset), max(interim.dataset)) +
  xlim(min(interim.dataset), max(interim.dataset)) +
  ggtitle(paste0("Comparison of SGLT2i predictions (mmol/mol) (n=", format(interim.dataset %>% nrow(),big.mark=",",scientific=FALSE), ")")) +
  theme_bw()

# Plot histogram of difference between both predictions
plot_sup_5.2 <- interim.dataset %>%
  mutate(diff = as.numeric(pred.SGLT2.bcf - pred.SGLT2.lm)) %>%
  ggplot(aes(x = diff)) +
  geom_density() +
  xlab("HbA1c difference (mmol/mol)") +
  ggtitle("Difference between HbA1c predictions (mmol/mol)") +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

plot_sup_5 <- patchwork::wrap_plots(list(plot_sup_5.1, plot_sup_5.2), ncol = 2, nrow = 1) +
  plot_annotation(tag_levels = "A")

## PDF with the plot for the comparison  
pdf(width = 11, height = 5, "Plots/Paper/Sup_Mat/11.08.plot_sup_5.pdf")
plot_sup_5
dev.off()


#:--------------------------------------------------------------------------------
### Optimal therapy from both models combined

# predict for DPP4
dataset.dpp4 <- full.cohort.updated %>%
  select(-patid, -pated) %>%
  mutate(drugclass = factor("DPP4", levels = c("SGLT2", "DPP4")))
predictions.dpp4 <- predict(m1, dataset.dpp4)

interim.dataset <- patient_predicted_outcomes %>%
  left_join(full.cohort.updated %>%
              select(patid, pated) %>%
              cbind(pred.DPP4 = predictions.dpp4), by = c("patid", "pated")) %>%
  drop_na() %>%
  mutate(best_drug = ifelse(pred.SGLT2 < pred.GLP1 & pred.SGLT2 < pred.DPP4, "SGLT2i",
                            ifelse(pred.GLP1 < pred.SGLT2 & pred.GLP1 < pred.DPP4, "GLP1-RA",
                                   ifelse(pred.DPP4 < pred.SGLT2 & pred.DPP4 < pred.GLP1, "DPP4i", NA))),
         best_drug = factor(best_drug))

# Plot bar plot of the numbers of best drugs
plot_sup_6 <- interim.dataset %>%
  select(best_drug) %>%
  table() %>%
  as.data.frame() %>%
  rename("best_drug" = ".") %>%
  ggplot(aes(x = best_drug, y = Freq, fill = best_drug)) +
  geom_bar(stat="identity") +
  geom_text(aes(label = format(Freq,big.mark=",",scientific=FALSE)), vjust = -0.5) +
  xlab("Optimal predicted therapy") +
  ylab("Number of patients") +
  # ggtitle(paste0("Predicted optimal therapy (n=", format(interim.dataset%>%nrow(),big.mark=",",scientific=FALSE), ")")) +
  scale_fill_manual(values = c("red", "dodgerblue2", "#f1a340")) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_y_continuous(label=comma)

## PDF with the plot for the best therapy
pdf(width = 7, height = 7, "Plots/Paper/Sup_Mat/11.08.plot_sup_6.pdf")
plot_sup_6
dev.off()


#:--------------------------------------------------------------------------------
### ROC AUC curve for propensity score model

require(pROC)

# load PS model
bart_ps_model_final <- readRDS("Samples/SGLT2-GLP1/Aurum/ps_model/bart_ps_model_final.rds")

# load training dataset
ps.model.train <- set_up_data_sglt2_glp1(dataset.type = "ps.model.train")

# combine probabilities with real values
predictions.train <- bart_ps_model_final$p_hat_train

df.train <- as.data.frame(cbind(predictions.train, ps.model.train["drugclass"] %>%
                            mutate(drugclass = ifelse(drugclass == "SGLT2", 0, 1)))) %>%
  set_names(c("predictions", "labels"))


performance.train <- pROC::roc(df.train$labels, df.train$predictions)

plot_sup_7.1 <- ggroc(performance.train, legacy.axes = "TRUE") +
  geom_abline(aes(intercept = 0, slope = 1), lty = "dashed", colour = "red") +
  theme_bw() +
  labs(
    x = "1 - Specificity",
    y = "Sensitivity",
    title = paste0("Derivation cohort (n=", format(nrow(ps.model.train),big.mark=",",scientific=FALSE), ")"),
    subtitle = "Sensitivity vs 1-Specificity"
  )

precision_recall.train <- performance.train %>%
  coords(ret = "all", transpose = FALSE) %>%
  select(precision, recall)

plot_sup_7.2 <- precision_recall.train %>%
  ggplot(aes(x = recall, y = precision)) +
  geom_hline(aes(yintercept = 0.5), lty = "dashed", colour = "grey") +
  geom_path() +
  theme_bw() +
  ylim(0, 1) +
  labs(
    x = "Recall",
    y = "Precision",
    title = paste0("Derivation cohort (n=", format(nrow(ps.model.train),big.mark=",",scientific=FALSE), ")"),
    subtitle = "Precision vs Recall"
  )



# load PS model predictions for testing data
prop_score_testing_data <- readRDS("Samples/SGLT2-GLP1/Aurum/ps_model/prop_score_testing_data.rds")

# load testing dataset
ps.model.test <- set_up_data_sglt2_glp1(dataset.type = "ps.model.test")

# combined probabilities with real values
predictions.test <- prop_score_testing_data

df.test <- as.data.frame(cbind(predictions.test, ps.model.test["drugclass"] %>%
                            mutate(drugclass = ifelse(drugclass == "SGLT2", 0, 1)))) %>%
  set_names(c("predictions", "labels"))


performance.test <- pROC::roc(df.test$labels, df.test$predictions)

plot_sup_7.3 <- ggroc(performance.test, legacy.axes = "TRUE") +
  geom_abline(aes(intercept = 0, slope = 1), lty = "dashed", colour = "red") +
  theme_bw() +
  labs(
    x = "1 - Specificity",
    y = "Sensitivity",
    title = paste0("Validation cohort (n=", format(nrow(ps.model.test),big.mark=",",scientific=FALSE), ")"),
    subtitle = "Sensitivity vs 1-Specificity"
  )

precision_recall.test <- performance.test %>%
  coords(ret = "all", transpose = FALSE) %>%
  select(precision, recall)

plot_sup_7.4 <- precision_recall.test %>%
  ggplot(aes(x = recall, y = precision)) +
  geom_hline(aes(yintercept = 0.5), lty = "dashed", colour = "grey") +
  geom_path() +
  theme_bw() +
  ylim(0, 1) +
  labs(
    x = "Recall",
    y = "Precision",
    title = paste0("Validation cohort (n=", format(nrow(ps.model.test),big.mark=",",scientific=FALSE), ")"),
    subtitle = "Precision vs Recall"
  )


plot_sup_7 <- patchwork::wrap_plots(list(plot_sup_7.1, plot_sup_7.2, plot_sup_7.3, plot_sup_7.4), ncol = 2, nrow = 2) +
    plot_annotation(tag_levels = list(c("A.1", "A.2", "B.1", "B.2")))


## PDF with the plot for the best therapy
pdf(width = 8, height = 8, "Plots/Paper/Sup_Mat/11.08.plot_sup_7.pdf")
plot_sup_7
dev.off()



#:--------------------------------------------------------------------------------
### Variable selection for PS model

# load proportions of variable selection
vs_bart_ps_model <- readRDS("Samples/SGLT2-GLP1/Aurum/ps_model/vs_bart_ps_model.rds")

# dataset with values necessary for the plot
df <- t(apply(vs_bart_ps_model$permute_mat, 2, range)) %>%
  as.data.frame() %>%
  rename("min" = "V1",
         "max" = "V2") %>%
  cbind(labels = factor(names(vs_bart_ps_model$var_true_props_avg), 
                        levels = names(vs_bart_ps_model$var_true_props_av), 
                        labels = c("BMI", "Year of drug start", 
                                   "eGFR", "HbA1c", 
                                   "Number drugs ever prescribed - 5+", "Current age",
                                   "Duration of diabetes", "Number drugs ever prescribed - 2", 
                                   "Number drugs ever prescribed - 3", "Number other current drugs - 3",
                                   "Number other current drugs - 5+", "Heart failure - No",
                                   "Number other current drugs - 4", "Number drugs ever prescribed - 4", 
                                   "Heart failure - Yes", "Ethnicity - White",
                                   "Hospitalisation - No", "Chronic liver disease - Yes",
                                   "Sex - Female", "Number other current drugs - 2",
                                   "Hospitalisation - Yes", "Diabetic Nephropathy - Yes",
                                   "Chronic liver disease - No", "Number other current drugs - 1",
                                   "Retinopathy - No", "Sex - Male",
                                   "Retinopathy - Yes", "Neuropathy - No",
                                   "Ethnicity - South Asian", "Angina - Yes",
                                   "Atrial fibrillation - No", "Diabetic Nephropathy - No",
                                   "Neuropathy - Yes", "Peripheral arterial disease - No",
                                   "Angina - No", "Ethnicity - Mixed",
                                   "Transient ischaemic attack - No", "Ethnicity - Black",
                                   "Atrial fibrillation - Yes", "Ethnicity - Other",
                                   "Myocardial infarction - Yes", "Deprivation - 9",
                                   "Stroke - No", "Ischaemic heart disease - Yes",
                                   "Ischaemic heart disease - No", "Smoking status - Active",
                                   "Deprivation - 6", "Peripheral arterial disease - Yes",
                                   "Transient ischaemic attack - Yes", "Smoking status - Non-smoker",
                                   "Cardiac revascularisation - Yes", "Cardiac revascularisation - No",
                                   "Stroke - Yes", "Myocardial infarction - No",
                                   "Smoking status - Ex-smoker", "Hypertension - Yes",
                                   "Hypertension - No", "Deprivation - 8",
                                   "Deprivation - 1", "Deprivation - 5",
                                   "Deprivation - 3", "Deprivation - 4",
                                   "Deprivation - 7", "Deprivation - 10",
                                   "Deprivation - 2")),
        true.prop = vs_bart_ps_model$var_true_props_avg) %>%
  mutate(min = as.numeric(min),
         max = as.numeric(max),
         mean = (min+max)/2,
         max.values = ifelse(labels == "agetx", max, mean),
         colour.toggle = factor(ifelse(true.prop > max.values, 1, 0)),
         min = 0) %>%
  gather(label, plot.values, -labels, -max, -true.prop, -colour.toggle) %>%
  filter(labels != "Hospitalisation - No")  %>%
  filter(labels != "Hospitalisation - Yes") 

# plot the variable selection
plot_sup_8 <- df %>%
  ggplot(aes(x = labels, y = plot.values)) +
  geom_boxplot(width = 0, colour = "chartreuse4") +
  geom_point(aes(y = true.prop, shape = colour.toggle), size = 3) +
  scale_shape_manual(values = c(1, 16)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, colour = c(rep("black", 6), "grey70", rep("black", 3), rep("grey70", 2), rep("black", 2), "grey70", "black", rep("grey70", 2), "black", rep("grey70", 44))),
        legend.position = "none",
        axis.title.x = element_blank()) +
  labs(
    y = "Proportion included"
  )


## PDF with the plot for the best therapy
pdf(width = 11, height = 5, "Plots/Paper/Sup_Mat/11.08.plot_sup_8.pdf")
plot_sup_8
dev.off()


#:--------------------------------------------------------------------------------
### Comparing predictions of mu and tau for BCF prop / no prop

# comparing mu predictions
prop.score.comparison.mu <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/prop.score.comparison.mu.rds")

# comparing tau predictions
prop.score.comparison.tau <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/prop.score.comparison.rds")

# plot comparisons
plot_sup_9 <- patchwork::wrap_plots(
  
  # scatter plot of mu values
  prop.score.comparison.mu %>%
    as.data.frame() %>%
    ggplot(aes(x = bcf_prop, y = bcf_no_prop)) +
    geom_abline(aes(intercept = 0, slope = 1), colour = "red", lty = "dashed") +
    geom_point() +
    theme_bw() +
    labs(
      x = "BCF model with propensity score",
      y = "BCF model without propensity score",
      title = "Predicted HbA1c outcome (mmol/mol)",
      subtitle = "BCF without propensity score vs BCF with propensity score"
    )
      
  ,
  
  # density plot of residuals for mu (BCF without prop - BCF with prop)
  prop.score.comparison.mu %>%
    as.data.frame() %>%
    mutate(effect.difference = bcf_no_prop - bcf_prop) %>%
    ggplot() +
    theme_bw() +
    geom_density(aes(x = effect.difference)) +
    labs(
      x = "Difference in predicted HbA1c outcome (mmol/mol)",
      y = "density",
      title = "Predicted HbA1c outcome (mmol/mol)",
      subtitle = "BCF without propensity score - BCF with propensity score"
    ) +
    xlim(-10, 10)
  
  ,
  
  # scatter plot of tau values
  prop.score.comparison.tau %>%
    as.data.frame() %>%
    ggplot(aes(x = bcf_prop, y = bcf_no_prop)) +
    theme_bw() +
    geom_abline(aes(intercept = 0, slope = 1), colour = "red", lty = "dashed") +
    geom_point() +
    labs(
      x = "BCF model with propensity score",
      y = "BCF model without propensity score",
      title = "Predicted conditional average treatment effects (CATE)",
      subtitle = "BCF without propensity score vs BCF with propensity score"
    )
  
  ,
  
  # density plot of residuals for tau (BCF without prop - BCF with prop)
  prop.score.comparison.tau %>%
    as.data.frame() %>%
    mutate(effect.difference = bcf_no_prop - bcf_prop) %>%
    ggplot() +
    theme_bw() +
    geom_density(aes(x = effect.difference)) +
    labs(
      x = "Difference in predicted CATE (mmol/mol)",
      y = "density",
      title = "Predicted conditional average treatment effects (CATE)",
      subtitle = "BCF without propensity score - BCF with propensity score"
    ) +
    xlim(-4, 4)
  
  , ncol = 2, nrow = 2) +
  plot_annotation(tag_levels = list(c("A.1", "A.2", "B.1", "B.2")))


## PDF with the plot for the best therapy
pdf(width = 11, height = 11, "Plots/Paper/Sup_Mat/11.08.plot_sup_9.pdf")
plot_sup_9
dev.off()





#:--------------------------------------------------------------------------------
### Calibration of the model using PSM and PSM + adjusted

# load development cohort model calibration
ATE_matching_1_1_validation_dev <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/assessment/ATE_matching_1_1_validation_dev.rds")

ATE_matching_1_1_validation_val <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/assessment/ATE_matching_1_1_validation_val.rds")

ATE_matching_1_1_adjust_validation_dev <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/assessment/ATE_matching_1_1_adjust_validation_dev.rds")

ATE_matching_1_1_adjust_validation_val <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/assessment/ATE_matching_1_1_adjust_validation_val.rds")

# plot development cohort calibration
plot_ATE_matching_1_1_validation_dev <- ATE_plot(ATE_matching_1_1_validation_dev[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -12, 12) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(
    title = paste0("Derivation cohort (n=", format(sum(ATE_matching_1_1_validation_dev[["effects"]]$n_drug1)*2,big.mark=",",scientific=FALSE), ")"),
    subtitle = "Propensity score matching"
  )

plot_ATE_matching_1_1_validation_val <- ATE_plot(ATE_matching_1_1_validation_val[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -12, 12) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(
    title = paste0("Validation cohort (n=", format(sum(ATE_matching_1_1_validation_val[["effects"]]$n_drug1)*2,big.mark=",",scientific=FALSE), ")"),
    subtitle = "Propensity score matching"
  )

plot_ATE_matching_1_1_adjust_validation_dev <- ATE_plot(ATE_matching_1_1_adjust_validation_dev[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -12, 12) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(
    title = paste0("Derivation cohort (n=", format(sum(ATE_matching_1_1_validation_dev[["effects"]]$n_drug1)*2,big.mark=",",scientific=FALSE), ")"),
    subtitle = "Propensity score matching and estimate adjustment"
  )

plot_ATE_matching_1_1_adjust_validation_val <- ATE_plot(ATE_matching_1_1_adjust_validation_val[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -12, 12) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(
    title = paste0("Validation cohort (n=", format(sum(ATE_matching_1_1_validation_val[["effects"]]$n_drug1)*2,big.mark=",",scientific=FALSE), ")"),
    subtitle = "Propensity score matching and estimate adjustment"
  )

plot_sup_10 <- (plot_ATE_matching_1_1_validation_dev | plot_ATE_matching_1_1_validation_val) / (plot_ATE_matching_1_1_adjust_validation_dev | plot_ATE_matching_1_1_adjust_validation_val) +
  plot_annotation(tag_levels = list(c("A.1", "A.2", "B.1", "B.2")))


## PDF with the plot for the best therapy
pdf(width = 10, height = 10, "Plots/Paper/Sup_Mat/11.08.plot_sup_10.pdf")
plot_sup_10
dev.off()



#:--------------------------------------------------------------------------------
### Clinical subgroup HbA1c outcome

## PSM

# levels needed
interval_breaks <- c(-5, -3, 0, 3, 5)

variables_mu <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/variables_mu.rds")

variables_tau <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/variables_tau.rds")

## Read in propensity scores
patient_prop_scores <- readRDS("Samples/SGLT2-GLP1/Aurum/ps_model/patient_prop_scores.rds")

## Read in treatment effects
treatment_effects <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/patient_effects.rds")

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


breakdown_adjust <- unique(c(variables_mu, variables_tau))
# categorical variables in breakdown
factors <- sapply(group.full.dataset[,breakdown_adjust], is.factor)

matching_hba1c <- MatchIt::matchit(
  formula = formula(paste0("drugclass ~ agetx + t2dmduration + prehba1c + preegfr + prealt +", paste(breakdown_adjust[factors], collapse = "+"))),
  data = group.full.dataset,
  method = "nearest",
  distance = group.full.dataset[,"prop.score"],
  replace = FALSE,
  m.order = "largest",
  caliper = 0.05,
  mahvars = NULL, estimand = "ATT", exact = NULL, antiexact = NULL, discard = "none", reestimate = FALSE, s.weights = NULL, std.caliper = TRUE, ratio = 1, verbose = FALSE, include.obj = FALSE,
)

# require(cobalt)
# cobalt::love.plot(matching_hba1c, binary = "std", thresholds = c(m = .1), sample.names = c("Unmatched", "Matched"), title = "Full dataset")

n_drug <- 2*sum(!is.na(matching_hba1c$match.matrix))

group.full.dataset.matched <- group.full.dataset %>%
  slice(which(matching_hba1c$weights == 1))


# predictions for the hba1c adjusted model
predictions_hba1c_stan_psm_1_1_overall <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_hba1c_stan_psm_1_1_overall.rds") %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci))

predictions_hba1c_stan_psm_1_1_full <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_hba1c_stan_psm_1_1_full.rds") %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci))

# calculate axis limits for the plots
hba1c_strata_axis_min <- plyr::round_any(floor(min(c(predictions_hba1c_stan_psm_1_1_overall %>% select(c("mean","lci","uci")) %>% min(),
                                                     predictions_hba1c_stan_psm_1_1_full %>% select(c("mean","lci","uci")) %>% min()))), 5, f = floor)

hba1c_strata_axis_max <- plyr::round_any(floor(max(c(predictions_hba1c_stan_psm_1_1_overall %>% select(c("mean","lci","uci")) %>% max(),
                                                     predictions_hba1c_stan_psm_1_1_full %>% select(c("mean","lci","uci")) %>% max()))), 5, f = ceiling)


# plot predicted hba1c change
plot_a <- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = 50, lci = 50, uci = 50, drugclass = "SGLT2"),
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = 50, lci = 50, uci = 50, drugclass = "GLP1"),
  predictions_hba1c_stan_psm_1_1_overall %>% select(intervals, drugclass, mean, lci, uci) %>% slice(1:6),
  cbind(intervals = "Predicted HbA1c benefit on GLP-1 RA", mean = 50, lci = 50, uci = 50, drugclass = "SGLT2"),
  cbind(intervals = "Predicted HbA1c benefit on GLP-1 RA", mean = 50, lci = 50, uci = 50, drugclass = "GLP1"),
  predictions_hba1c_stan_psm_1_1_overall %>% select(intervals, drugclass, mean, lci, uci) %>% slice(7:12),
  predictions_hba1c_stan_psm_1_1_full %>% select(drugclass, mean, lci, uci) %>% mutate(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.full.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", format(group.full.dataset.matched%>%filter(intervals == levels(group.full.dataset.matched$intervals)[1])%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                            ifelse(intervals == levels(group.full.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", format(group.full.dataset.matched%>%filter(intervals == levels(group.full.dataset.matched$intervals)[2])%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                   ifelse(intervals == levels(group.full.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", format(group.full.dataset.matched%>%filter(intervals == levels(group.full.dataset.matched$intervals)[3])%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                          ifelse(intervals == levels(group.full.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", format(group.full.dataset.matched%>%filter(intervals == levels(group.full.dataset.matched$intervals)[4])%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                 ifelse(intervals == levels(group.full.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", format(group.full.dataset.matched%>%filter(intervals == levels(group.full.dataset.matched$intervals)[5])%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                        ifelse(intervals == levels(group.full.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", format(group.full.dataset.matched%>%filter(intervals == levels(group.full.dataset.matched$intervals)[6])%>%nrow(),big.mark=",",scientific=FALSE), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  mutate(group = ifelse(group == "SGLT2", "SGLT2i", "GLP-1RA")) %>%
  forestplot(ci.vertices = TRUE,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices.height = 0.1,
             title = "Average HbA1c change",
             clip = c(hba1c_strata_axis_min, hba1c_strata_axis_max),
             xticks = seq(hba1c_strata_axis_min, hba1c_strata_axis_max, 5),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Predicted HbA1c change (mmol/mol)") %>%
  fp_add_header(paste0("Overall population (n=", format(group.full.dataset.matched%>%nrow(),big.mark=",",scientific=FALSE), ")")) %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE)) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = "black")


### B - Weight

## Read in data for weight
weight.dataset <- set_up_data_sglt2_glp1(dataset.type = "weight.dataset") %>%
  left_join(patient_prop_scores, by = c("patid", "pated")) %>%
  left_join(treatment_effects, by = c("patid", "pated")) %>%
  mutate(w.change = postweight - preweight)


group.weight.dataset <- group_values(data = weight.dataset,
                                     variable = "effects",
                                     breaks = interval_breaks) %>%
  drop_na(intervals)

breakdown_adjust <- unique(c(variables_mu, variables_tau))
# categorical variables in breakdown
factors <- sapply(group.weight.dataset[,breakdown_adjust], is.factor)

matching_weight <- MatchIt::matchit(
  formula = formula(paste0("drugclass ~ preweight + agetx + t2dmduration + prehba1c + preegfr + prealt +", paste(breakdown_adjust[factors], collapse = "+"))),
  data = group.weight.dataset,
  method = "nearest",
  distance = group.weight.dataset[,"prop.score"],
  replace = FALSE,
  m.order = "largest",
  caliper = 0.05,
  mahvars = NULL, estimand = "ATT", exact = NULL, antiexact = NULL, discard = "none", reestimate = FALSE, s.weights = NULL, std.caliper = TRUE, ratio = 1, verbose = FALSE, include.obj = FALSE,
)

# require(cobalt)
# cobalt::love.plot(matching_weight, binary = "std", thresholds = c(m = .1), sample.names = c("Unmatched", "Matched"), title = "Full dataset")

n_drug <- 2*sum(!is.na(matching_weight$match.matrix))

group.weight.dataset.matched <- group.weight.dataset %>%
  slice(which(matching_weight$weights == 1))


# predictions for the weight adjusted model
predictions_weight_stan_psm_1_1_overall <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_weight_stan_psm_1_1_overall.rds")

# predictions for the weight adjusted model full
predictions_weight_stan_psm_1_1_full <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_weight_stan_psm_1_1_full.rds")


weight_strata_axis_min <- plyr::round_any(floor(min(c(predictions_weight_stan_psm_1_1_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min(),
                                                      predictions_weight_stan_psm_1_1_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min()))), 1, f = floor)

weight_strata_axis_max <- plyr::round_any(ceiling(max(c(predictions_weight_stan_psm_1_1_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max(),
                                                        predictions_weight_stan_psm_1_1_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max()))), 2, f = ceiling)


plot_b <- rbind(
  cbind(mean = 50, lci = 50, uci = 50, drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(mean = 50, lci = 50, uci = 50, drugclass = "GLP1", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_weight_stan_psm_1_1_overall %>%
    slice(c(1:6)),
  cbind(mean = 50, lci = 50, uci = 50, drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(mean = 50, lci = 50, uci = 50, drugclass = "GLP1", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_weight_stan_psm_1_1_overall %>%
    slice(-c(1:6)),
  predictions_weight_stan_psm_1_1_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == levels(group.weight.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", format(group.weight.dataset.matched%>%filter(intervals==levels(group.weight.dataset.matched$intervals)[1])%>%nrow(),big.mark=",",scientific=FALSE),")"),
                            ifelse(intervals == levels(group.weight.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", format(group.weight.dataset.matched%>%filter(intervals==levels(group.weight.dataset.matched$intervals)[2])%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                   ifelse(intervals == levels(group.weight.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", format(group.weight.dataset.matched%>%filter(intervals==levels(group.weight.dataset.matched$intervals)[3])%>%nrow(),big.mark=",",scientific=FALSE),")"),
                                          ifelse(intervals == levels(group.weight.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", format(group.weight.dataset.matched%>%filter(intervals==levels(group.weight.dataset.matched$intervals)[4])%>%nrow(),big.mark=",",scientific=FALSE),")"),
                                                 ifelse(intervals == levels(group.weight.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", format(group.weight.dataset.matched%>%filter(intervals==levels(group.weight.dataset.matched$intervals)[5])%>%nrow(),big.mark=",",scientific=FALSE),")"),
                                                        ifelse(intervals == levels(group.weight.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", format(group.weight.dataset.matched%>%filter(intervals==levels(group.weight.dataset.matched$intervals)[6])%>%nrow(),big.mark=",",scientific=FALSE),")"), intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  forestplot(ci.vertices = TRUE,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices.height = 0.1,
             title = "Average weight change",
             clip = c(weight_strata_axis_min, weight_strata_axis_max),
             xticks = seq(weight_strata_axis_min, weight_strata_axis_max, 1),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Predicted weight change (kg)") %>%
  fp_add_header(paste0("Overall population (n=", format(group.weight.dataset.matched%>%nrow(),big.mark=",",scientific=FALSE), ")")) %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE)) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = "black")


### C - Discontinuation

discontinuation.dataset <- set_up_data_sglt2_glp1(dataset.type = "discontinuation.dataset") %>%
  left_join(patient_prop_scores, by = c("patid", "pated")) %>%
  left_join(treatment_effects, by = c("patid", "pated")) %>%
  mutate(stopdrug_6m_3mFU = factor(stopdrug_6m_3mFU))

group.discontinuation.dataset <- group_values(data = discontinuation.dataset,
                                              variable = "effects",
                                              breaks = interval_breaks) %>%
  drop_na(intervals)

breakdown_adjust <- unique(c(variables_mu, variables_tau))
# categorical variables in breakdown
factors <- sapply(group.discontinuation.dataset[,breakdown_adjust], is.factor)

matching_discontinuation <- MatchIt::matchit(
  formula = formula(paste0("drugclass ~ agetx + t2dmduration +  + prehba1c + preegfr + prealt +", paste(breakdown_adjust[factors], collapse = "+"))),
  data = group.discontinuation.dataset,
  method = "nearest",
  distance = group.discontinuation.dataset[,"prop.score"], 
  replace = FALSE,
  m.order = "largest",
  caliper = 0.05,
  mahvars = NULL, estimand = "ATT", exact = NULL, antiexact = NULL, discard = "none", reestimate = FALSE, s.weights = NULL, std.caliper = TRUE, ratio = 1, verbose = FALSE, include.obj = FALSE,
)

# require(cobalt)
# cobalt::love.plot(matching_discontinuation, binary = "std", thresholds = c(m = .1), sample.names = c("Unmatched", "Matched"), title = "Full dataset")

n_drug <- 2*sum(!is.na(matching_discontinuation$match.matrix))

group.discontinuation.dataset.matched <- group.discontinuation.dataset %>%
  slice(which(matching_discontinuation$weights == 1))


# predictions for the discontinuation adjusted model
predictions_discontinuation_stan_psm_1_1_overall <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_discontinuation_stan_psm_1_1_overall.rds")

# predictions for the discontinuation model sex strata
predictions_discontinuation_stan_psm_1_1_full <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_discontinuation_stan_psm_1_1_full.rds")

discontinuation_overall_axis_max <- plyr::round_any(ceiling(max(c(predictions_discontinuation_stan_psm_1_1_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() * 100,
                                                                  predictions_discontinuation_stan_psm_1_1_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() * 100))), 5, f = ceiling)

discontinuation_strata_axis_max <- plyr::round_any(ceiling(max(c(predictions_discontinuation_stan_psm_1_1_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() * 100,
                                                                 predictions_discontinuation_stan_psm_1_1_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() * 100))), 5, f = ceiling)


plot_c <- rbind(
  cbind(mean = 150, lci = 150, uci = 150, drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(mean = 150, lci = 150, uci = 150, drugclass = "GLP1", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_discontinuation_stan_psm_1_1_overall %>%
    slice(c(1:6)) %>%
    mutate(mean = as.numeric(mean)*100,
           lci = as.numeric(lci)*100,
           uci = as.numeric(uci)*100),
  cbind(mean = 150, lci = 150, uci = 150, drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(mean = 150, lci = 150, uci = 150, drugclass = "GLP1", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_discontinuation_stan_psm_1_1_overall %>%
    slice(-c(1:6)) %>%
    mutate(mean = as.numeric(mean)*100,
           lci = as.numeric(lci)*100,
           uci = as.numeric(uci)*100),
  predictions_discontinuation_stan_psm_1_1_full %>% cbind(intervals = "Average treatment effect") %>%
    mutate(mean = as.numeric(mean)*100,
           lci = as.numeric(lci)*100,
           uci = as.numeric(uci)*100)
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == levels(group.discontinuation.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", format(group.discontinuation.dataset.matched%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[1])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.discontinuation.dataset.matched%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[1])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                            ifelse(intervals == levels(group.discontinuation.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", format(group.discontinuation.dataset.matched%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[2])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.discontinuation.dataset.matched%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[2])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                   ifelse(intervals == levels(group.discontinuation.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", format(group.discontinuation.dataset.matched%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[3])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.discontinuation.dataset.matched%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[3])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                          ifelse(intervals == levels(group.discontinuation.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", format(group.discontinuation.dataset.matched%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[4])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.discontinuation.dataset.matched%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[4])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                 ifelse(intervals == levels(group.discontinuation.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", format(group.discontinuation.dataset.matched%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[5])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.discontinuation.dataset.matched%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[5])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                        ifelse(intervals == levels(group.discontinuation.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", format(group.discontinuation.dataset.matched%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[6])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.discontinuation.dataset.matched%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[6])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"), intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  forestplot(ci.vertices = TRUE,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices.height = 0.1,
             title = "Risk of discontinuation",
             clip = c(0, discontinuation_overall_axis_max),
             xticks = seq(0, discontinuation_overall_axis_max, 5),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Discontinuation (%)") %>%
  fp_add_header(paste0("Overall population (n=", format(group.discontinuation.dataset.matched%>%nrow(),big.mark=",",scientific=FALSE), ")")) %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE)) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = "black")


### D - Risk of developing microvascular complications


patient_prop_scores_qrisk <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/patient_prop_scores_qrisk.rds")

micro_comp.dataset <- set_up_data_sglt2_glp1(dataset.type="micro_comp.dataset") %>%
  left_join(patient_prop_scores_qrisk, by = c("patid", "pated")) %>%
  left_join(treatment_effects, by = c("patid", "pated"))

group.micro_comp.dataset <- group_values(data = micro_comp.dataset,
                                         variable = "effects",
                                         breaks = interval_breaks) %>%
  drop_na(intervals)

matching_micro_comp <- MatchIt::matchit(
  formula = formula(paste0("drugclass ~  agetx + t2dmduration + prehba1c + preegfr + prealt + drugline + ncurrtx + sex + preneuropathy + preretinopathy")),
  data = group.micro_comp.dataset,
  method = "nearest",
  distance = group.micro_comp.dataset[,"prop.score"],
  replace = FALSE,
  m.order = "largest",
  caliper = 0.05,
  mahvars = NULL, estimand = "ATT", exact = NULL, antiexact = NULL, discard = "none", reestimate = FALSE, s.weights = NULL, std.caliper = TRUE, ratio = 1, verbose = FALSE, include.obj = FALSE,
)


# require(cobalt)
# cobalt::love.plot(matching_micro_comp, binary = "std", thresholds = c(m = .1), sample.names = c("Unmatched", "Matched"), title = "Full dataset")

n_drug <- 2*sum(!is.na(matching_micro_comp$match.matrix))

group.micro_comp.dataset.matched <- group.micro_comp.dataset %>%
  slice(which(matching_micro_comp$weights == 1))


# predictions for the CVD outcomes in the population with no CVD/HF/CKD
predictions_micro_comp_micro_comp_stan_psm_1_1_overall <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_micro_comp_micro_comp_stan_psm_1_1_overall.rds")

# predictions for the CVD outcomes in the population with no CVD/HF/CKD
predictions_micro_comp_micro_comp_stan_psm_1_1_full <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_micro_comp_micro_comp_stan_psm_1_1_full.rds")

# plot of forest plot
plot_d <- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_micro_comp_micro_comp_stan_psm_1_1_overall %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = NA, lci = NA, uci = NA),
  predictions_micro_comp_micro_comp_stan_psm_1_1_overall %>% slice(4:6),
  predictions_micro_comp_micro_comp_stan_psm_1_1_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.micro_comp.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", format(group.micro_comp.dataset.matched%>%filter(intervals==levels(group.micro_comp.dataset.matched$intervals)[1])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.micro_comp.dataset.matched%>%filter(intervals==levels(group.micro_comp.dataset.matched$intervals)[1])%>%filter(postdrug_micro_comp_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                            ifelse(intervals == levels(group.micro_comp.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", format(group.micro_comp.dataset.matched%>%filter(intervals==levels(group.micro_comp.dataset.matched$intervals)[2])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.micro_comp.dataset.matched%>%filter(intervals==levels(group.micro_comp.dataset.matched$intervals)[2])%>%filter(postdrug_micro_comp_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                   ifelse(intervals == levels(group.micro_comp.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", format(group.micro_comp.dataset.matched%>%filter(intervals==levels(group.micro_comp.dataset.matched$intervals)[3])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.micro_comp.dataset.matched%>%filter(intervals==levels(group.micro_comp.dataset.matched$intervals)[3])%>%filter(postdrug_micro_comp_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                          ifelse(intervals == levels(group.micro_comp.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", format(group.micro_comp.dataset.matched%>%filter(intervals==levels(group.micro_comp.dataset.matched$intervals)[4])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.micro_comp.dataset.matched%>%filter(intervals==levels(group.micro_comp.dataset.matched$intervals)[4])%>%filter(postdrug_micro_comp_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                 ifelse(intervals == levels(group.micro_comp.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", format(group.micro_comp.dataset.matched%>%filter(intervals==levels(group.micro_comp.dataset.matched$intervals)[5])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.micro_comp.dataset.matched%>%filter(intervals==levels(group.micro_comp.dataset.matched$intervals)[5])%>%filter(postdrug_micro_comp_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                        ifelse(intervals == levels(group.micro_comp.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", format(group.micro_comp.dataset.matched%>%filter(intervals==levels(group.micro_comp.dataset.matched$intervals)[6])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.micro_comp.dataset.matched%>%filter(intervals==levels(group.micro_comp.dataset.matched$intervals)[6])%>%filter(postdrug_micro_comp_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci") %>%
  mutate(mean = exp(mean),
         lower = exp(lower),
         upper = exp(upper)) %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             zero = 1,
             title = "Microvascular complications",
             xlog = TRUE,
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazard ratio (95% CI, log scale)") %>%
  fp_add_header(paste0("Overall population (n=", format(nrow(group.micro_comp.dataset.matched),big.mark=",",scientific=FALSE), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = "black")


### E - Risk of developing MACE/CVD


patient_prop_scores_qrisk <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/patient_prop_scores_qrisk.rds")

no_co.dataset <- set_up_data_sglt2_glp1(dataset.type="no_co.dataset") %>%
  left_join(patient_prop_scores_qrisk, by = c("patid", "pated")) %>%
  left_join(treatment_effects, by = c("patid", "pated"))

group.no_co.dataset <- group_values(data = no_co.dataset,
                                    variable = "effects",
                                    breaks = interval_breaks) %>%
  drop_na(intervals)

breakdown_adjust <- unique(c(variables_mu, variables_tau))
# categorical variables in breakdown
factors <- sapply(group.no_co.dataset[,breakdown_adjust], is.factor)

matching_no_co <- MatchIt::matchit(
  formula = formula(paste0("drugclass ~  agetx + t2dmduration + prehba1c + preegfr + prealt + drugline + ncurrtx + sex + preneuropathy + preretinopathy")),
  data = group.no_co.dataset,
  method = "nearest",
  distance = group.no_co.dataset[,"prop.score"],
  replace = FALSE,
  m.order = "largest",
  caliper = 0.05,
  mahvars = NULL, estimand = "ATT", exact = NULL, antiexact = NULL, discard = "none", reestimate = FALSE, s.weights = NULL, std.caliper = TRUE, ratio = 1, verbose = FALSE, include.obj = FALSE,
)

# require(cobalt)
# cobalt::love.plot(matching_no_co, binary = "std", thresholds = c(m = .1), sample.names = c("Unmatched", "Matched"), title = "Full dataset")

n_drug <- 2*sum(!is.na(matching_no_co$match.matrix))

group.no_co.dataset.matched <- group.no_co.dataset %>%
  slice(which(matching_no_co$weights == 1))


# predictions for the CVD outcomes in the population with no CVD/HF/CKD
predictions_no_co_cvd_stan_psm_1_1_overall <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_no_co_cvd_stan_psm_1_1_overall.rds")

# predictions for the CVD outcomes in the population with no CVD/HF/CKD
predictions_no_co_cvd_stan_psm_1_1_full <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_no_co_cvd_stan_psm_1_1_full.rds")


# plot of forest plot
plot_e <- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_no_co_cvd_stan_psm_1_1_overall %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = NA, lci = NA, uci = NA),
  predictions_no_co_cvd_stan_psm_1_1_overall %>% slice(4:6),
  predictions_no_co_cvd_stan_psm_1_1_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[1])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[1])%>%filter(postdrug_mace_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                            ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[2])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[2])%>%filter(postdrug_mace_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                   ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[3])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[3])%>%filter(postdrug_mace_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                          ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[4])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[4])%>%filter(postdrug_mace_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                 ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[5])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[5])%>%filter(postdrug_mace_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                        ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[6])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[6])%>%filter(postdrug_mace_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci") %>%
  mutate(mean = exp(mean),
         lower = exp(lower),
         upper = exp(upper)) %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             zero = 1,
             title = "Major adverse cardiovascular events",
             xlog = TRUE,
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazard ratio (95% CI, log scale)") %>%
  fp_add_header(paste0("Overall population (n=", format(nrow(group.no_co.dataset.matched),big.mark=",",scientific=FALSE), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = "black")


### F - Risk of developing Heart Failure


# predictions for the HF outcomes in the population with no CVD/HF/CKD
predictions_no_co_hf_stan_psm_1_1_overall <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_no_co_hf_stan_psm_1_1_overall.rds")

# predictions for the CVD outcomes in the population with no CVD/HF/CKD
predictions_no_co_hf_stan_psm_1_1_full <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_no_co_hf_stan_psm_1_1_full.rds")

plot_f <- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_no_co_hf_stan_psm_1_1_overall %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = NA, lci = NA, uci = NA),
  predictions_no_co_hf_stan_psm_1_1_overall %>% slice(4:6),
  predictions_no_co_hf_stan_psm_1_1_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[1])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[1])%>%filter(postdrug_hf_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                            ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[2])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[2])%>%filter(postdrug_hf_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                   ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[3])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[3])%>%filter(postdrug_hf_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                          ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[4])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[4])%>%filter(postdrug_hf_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                 ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[5])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[5])%>%filter(postdrug_hf_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                        ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[6])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[6])%>%filter(postdrug_hf_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci") %>%
  mutate(mean = exp(mean),
         lower = exp(lower),
         upper = exp(upper)) %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             zero = 1,
             title = "Heart failure",
             xlog = TRUE,
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazard ratio (95% CI, log scale)") %>%
  fp_add_header(paste0("Overall population (n=", format(nrow(group.no_co.dataset.matched),big.mark=",",scientific=FALSE), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = "black")


pdf(width = 20, height = 10, "Plots/Paper/Sup_Mat/11.08.plot_sup_11.pdf")

grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 2,
                                           ncol = 3, heights = unit(c(5, 5), "null"))))
# 1 2 3
# 4 5 6
# first plot
pushViewport(viewport(layout.pos.row = 1,
                      layout.pos.col = 1))
plot_a
upViewport()

# second plot
pushViewport(viewport(layout.pos.row = 1,
                      layout.pos.col = 2))
plot_b
upViewport()

# third plot
pushViewport(viewport(layout.pos.row = 1,
                      layout.pos.col = 3))
plot_c
upViewport()

# forth plot
pushViewport(viewport(layout.pos.row = 2,
                      layout.pos.col = 1))
plot_d
upViewport()

# fifth plot
pushViewport(viewport(layout.pos.row = 2,
                      layout.pos.col = 2))
plot_e
upViewport()

# third plot
pushViewport(viewport(layout.pos.row = 2,
                      layout.pos.col = 3))
plot_f
upViewport()

# add legends
grid.text(expression(bold("A.1")), vp = viewport(layout.pos.row = 1, layout.pos.col = 1), hjust = 12, vjust = -18)
grid.text(expression(bold("A.2")), vp = viewport(layout.pos.row = 1, layout.pos.col = 2), hjust = 12, vjust = -18)
grid.text(expression(bold("A.3")), vp = viewport(layout.pos.row = 1, layout.pos.col = 3), hjust = 12, vjust = -17.5)
grid.text(expression(bold("B.1")), vp = viewport(layout.pos.row = 2, layout.pos.col = 1), hjust = 12, vjust = -18)
grid.text(expression(bold("B.2")), vp = viewport(layout.pos.row = 2, layout.pos.col = 2), hjust = 12, vjust = -18)
grid.text(expression(bold("B.3")), vp = viewport(layout.pos.row = 2, layout.pos.col = 3), hjust = 12, vjust = -17.5)

dev.off()



#:--------------------------------------------------------------------------------
### Clinical subgroup HbA1c outcome

## PSM + adjust

# levels needed
interval_breaks <- c(-5, -3, 0, 3, 5)

variables_mu <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/variables_mu.rds")

variables_tau <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/variables_tau.rds")

## Read in propensity scores
patient_prop_scores <- readRDS("Samples/SGLT2-GLP1/Aurum/ps_model/patient_prop_scores.rds")

## Read in treatment effects
treatment_effects <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/patient_effects.rds")

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


breakdown_adjust <- unique(c(variables_mu, variables_tau))
# categorical variables in breakdown
factors <- sapply(group.full.dataset[,breakdown_adjust], is.factor)

matching_hba1c <- MatchIt::matchit(
  formula = formula(paste0("drugclass ~ agetx + t2dmduration + prehba1c + preegfr + prealt +", paste(breakdown_adjust[factors], collapse = "+"))),
  data = group.full.dataset,
  method = "nearest",
  distance = group.full.dataset[,"prop.score"],
  replace = FALSE,
  m.order = "largest",
  caliper = 0.05,
  mahvars = NULL, estimand = "ATT", exact = NULL, antiexact = NULL, discard = "none", reestimate = FALSE, s.weights = NULL, std.caliper = TRUE, ratio = 1, verbose = FALSE, include.obj = FALSE,
)

# require(cobalt)
# cobalt::love.plot(matching_hba1c, binary = "std", thresholds = c(m = .1), sample.names = c("Unmatched", "Matched"), title = "Full dataset")

n_drug <- 2*sum(!is.na(matching_hba1c$match.matrix))

group.full.dataset.matched <- group.full.dataset %>%
  slice(which(matching_hba1c$weights == 1))


# predictions for the hba1c adjusted model
predictions_hba1c_stan_psm_1_1_adjusted_overall <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_hba1c_stan_psm_1_1_adjusted_overall.rds") %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci))

predictions_hba1c_stan_psm_1_1_adjusted_full <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_hba1c_stan_psm_1_1_adjusted_full.rds") %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci))

# calculate axis limits for the plots
hba1c_strata_axis_min <- plyr::round_any(floor(min(c(predictions_hba1c_stan_psm_1_1_adjusted_overall %>% select(c("mean","lci","uci")) %>% min(),
                                                     predictions_hba1c_stan_psm_1_1_adjusted_full %>% select(c("mean","lci","uci")) %>% min()))), 5, f = floor)

hba1c_strata_axis_max <- plyr::round_any(floor(max(c(predictions_hba1c_stan_psm_1_1_adjusted_overall %>% select(c("mean","lci","uci")) %>% max(),
                                                     predictions_hba1c_stan_psm_1_1_adjusted_full %>% select(c("mean","lci","uci")) %>% max()))), 5, f = ceiling)


# plot predicted hba1c change
plot_a <- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = 50, lci = 50, uci = 50, drugclass = "SGLT2"),
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = 50, lci = 50, uci = 50, drugclass = "GLP1"),
  predictions_hba1c_stan_psm_1_1_adjusted_overall %>% select(intervals, drugclass, mean, lci, uci) %>% slice(1:6),
  cbind(intervals = "Predicted HbA1c benefit on GLP-1 RA", mean = 50, lci = 50, uci = 50, drugclass = "SGLT2"),
  cbind(intervals = "Predicted HbA1c benefit on GLP-1 RA", mean = 50, lci = 50, uci = 50, drugclass = "GLP1"),
  predictions_hba1c_stan_psm_1_1_adjusted_overall %>% select(intervals, drugclass, mean, lci, uci) %>% slice(7:12),
  predictions_hba1c_stan_psm_1_1_adjusted_full %>% select(drugclass, mean, lci, uci) %>% mutate(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.full.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", format(group.full.dataset.matched%>%filter(intervals == levels(group.full.dataset.matched$intervals)[1])%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                            ifelse(intervals == levels(group.full.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", format(group.full.dataset.matched%>%filter(intervals == levels(group.full.dataset.matched$intervals)[2])%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                   ifelse(intervals == levels(group.full.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", format(group.full.dataset.matched%>%filter(intervals == levels(group.full.dataset.matched$intervals)[3])%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                          ifelse(intervals == levels(group.full.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", format(group.full.dataset.matched%>%filter(intervals == levels(group.full.dataset.matched$intervals)[4])%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                 ifelse(intervals == levels(group.full.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", format(group.full.dataset.matched%>%filter(intervals == levels(group.full.dataset.matched$intervals)[5])%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                        ifelse(intervals == levels(group.full.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", format(group.full.dataset.matched%>%filter(intervals == levels(group.full.dataset.matched$intervals)[6])%>%nrow(),big.mark=",",scientific=FALSE), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  mutate(group = ifelse(group == "SGLT2", "SGLT2i", "GLP-1RA")) %>%
  forestplot(ci.vertices = TRUE,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices.height = 0.1,
             title = "Average HbA1c change",
             clip = c(hba1c_strata_axis_min, hba1c_strata_axis_max),
             xticks = seq(hba1c_strata_axis_min, hba1c_strata_axis_max, 5),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Predicted HbA1c change (mmol/mol)") %>%
  fp_add_header(paste0("Overall population (n=", format(group.full.dataset.matched%>%nrow(),big.mark=",",scientific=FALSE), ")")) %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE)) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = "black")


### B - Weight

## Read in data for weight
weight.dataset <- set_up_data_sglt2_glp1(dataset.type = "weight.dataset") %>%
  left_join(patient_prop_scores, by = c("patid", "pated")) %>%
  left_join(treatment_effects, by = c("patid", "pated")) %>%
  mutate(w.change = postweight - preweight)


group.weight.dataset <- group_values(data = weight.dataset,
                                     variable = "effects",
                                     breaks = interval_breaks) %>%
  drop_na(intervals)

breakdown_adjust <- unique(c(variables_mu, variables_tau))
# categorical variables in breakdown
factors <- sapply(group.weight.dataset[,breakdown_adjust], is.factor)

matching_weight <- MatchIt::matchit(
  formula = formula(paste0("drugclass ~ preweight + agetx + t2dmduration + prehba1c + preegfr + prealt +", paste(breakdown_adjust[factors], collapse = "+"))),
  data = group.weight.dataset,
  method = "nearest",
  distance = group.weight.dataset[,"prop.score"],
  replace = FALSE,
  m.order = "largest",
  caliper = 0.05,
  mahvars = NULL, estimand = "ATT", exact = NULL, antiexact = NULL, discard = "none", reestimate = FALSE, s.weights = NULL, std.caliper = TRUE, ratio = 1, verbose = FALSE, include.obj = FALSE,
)

# require(cobalt)
# cobalt::love.plot(matching_weight, binary = "std", thresholds = c(m = .1), sample.names = c("Unmatched", "Matched"), title = "Full dataset")

n_drug <- 2*sum(!is.na(matching_weight$match.matrix))

group.weight.dataset.matched <- group.weight.dataset %>%
  slice(which(matching_weight$weights == 1))


# predictions for the weight adjusted model
predictions_weight_stan_psm_1_1_adjusted_overall <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_weight_stan_psm_1_1_adjusted_overall.rds")

# predictions for the weight adjusted model full
predictions_weight_stan_psm_1_1_adjusted_full <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_weight_stan_psm_1_1_adjusted_full.rds")


weight_strata_axis_min <- plyr::round_any(floor(min(c(predictions_weight_stan_psm_1_1_adjusted_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min(),
                                                      predictions_weight_stan_psm_1_1_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min()))), 1, f = floor)

weight_strata_axis_max <- plyr::round_any(ceiling(max(c(predictions_weight_stan_psm_1_1_adjusted_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max(),
                                                        predictions_weight_stan_psm_1_1_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max()))), 2, f = ceiling)


plot_b <- rbind(
  cbind(mean = 50, lci = 50, uci = 50, drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(mean = 50, lci = 50, uci = 50, drugclass = "GLP1", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_weight_stan_psm_1_1_adjusted_overall %>%
    slice(c(1:6)),
  cbind(mean = 50, lci = 50, uci = 50, drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(mean = 50, lci = 50, uci = 50, drugclass = "GLP1", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_weight_stan_psm_1_1_adjusted_overall %>%
    slice(-c(1:6)),
  predictions_weight_stan_psm_1_1_adjusted_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == levels(group.weight.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", format(group.weight.dataset.matched%>%filter(intervals==levels(group.weight.dataset.matched$intervals)[1])%>%nrow(),big.mark=",",scientific=FALSE),")"),
                            ifelse(intervals == levels(group.weight.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", format(group.weight.dataset.matched%>%filter(intervals==levels(group.weight.dataset.matched$intervals)[2])%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                   ifelse(intervals == levels(group.weight.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", format(group.weight.dataset.matched%>%filter(intervals==levels(group.weight.dataset.matched$intervals)[3])%>%nrow(),big.mark=",",scientific=FALSE),")"),
                                          ifelse(intervals == levels(group.weight.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", format(group.weight.dataset.matched%>%filter(intervals==levels(group.weight.dataset.matched$intervals)[4])%>%nrow(),big.mark=",",scientific=FALSE),")"),
                                                 ifelse(intervals == levels(group.weight.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", format(group.weight.dataset.matched%>%filter(intervals==levels(group.weight.dataset.matched$intervals)[5])%>%nrow(),big.mark=",",scientific=FALSE),")"),
                                                        ifelse(intervals == levels(group.weight.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", format(group.weight.dataset.matched%>%filter(intervals==levels(group.weight.dataset.matched$intervals)[6])%>%nrow(),big.mark=",",scientific=FALSE),")"), intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  forestplot(ci.vertices = TRUE,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices.height = 0.1,
             title = "Average weight change",
             clip = c(weight_strata_axis_min, weight_strata_axis_max),
             xticks = seq(weight_strata_axis_min, weight_strata_axis_max, 1),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Predicted weight change (kg)") %>%
  fp_add_header(paste0("Overall population (n=", format(group.weight.dataset.matched%>%nrow(),big.mark=",",scientific=FALSE), ")")) %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE)) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = "black")


### C - Discontinuation

discontinuation.dataset <- set_up_data_sglt2_glp1(dataset.type = "discontinuation.dataset") %>%
  left_join(patient_prop_scores, by = c("patid", "pated")) %>%
  left_join(treatment_effects, by = c("patid", "pated")) %>%
  mutate(stopdrug_6m_3mFU = factor(stopdrug_6m_3mFU))

group.discontinuation.dataset <- group_values(data = discontinuation.dataset,
                                              variable = "effects",
                                              breaks = interval_breaks) %>%
  drop_na(intervals)

breakdown_adjust <- unique(c(variables_mu, variables_tau))
# categorical variables in breakdown
factors <- sapply(group.discontinuation.dataset[,breakdown_adjust], is.factor)

matching_discontinuation <- MatchIt::matchit(
  formula = formula(paste0("drugclass ~ agetx + t2dmduration +  + prehba1c + preegfr + prealt +", paste(breakdown_adjust[factors], collapse = "+"))),
  data = group.discontinuation.dataset,
  method = "nearest",
  distance = group.discontinuation.dataset[,"prop.score"], 
  replace = FALSE,
  m.order = "largest",
  caliper = 0.05,
  mahvars = NULL, estimand = "ATT", exact = NULL, antiexact = NULL, discard = "none", reestimate = FALSE, s.weights = NULL, std.caliper = TRUE, ratio = 1, verbose = FALSE, include.obj = FALSE,
)

# require(cobalt)
# cobalt::love.plot(matching_discontinuation, binary = "std", thresholds = c(m = .1), sample.names = c("Unmatched", "Matched"), title = "Full dataset")

n_drug <- 2*sum(!is.na(matching_discontinuation$match.matrix))

group.discontinuation.dataset.matched <- group.discontinuation.dataset %>%
  slice(which(matching_discontinuation$weights == 1))


# predictions for the discontinuation adjusted model
predictions_discontinuation_stan_psm_1_1_adjusted_overall <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_discontinuation_stan_psm_1_1_adjusted_overall.rds")

# predictions for the discontinuation model sex strata
predictions_discontinuation_stan_psm_1_1_adjusted_full <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_discontinuation_stan_psm_1_1_adjusted_full.rds")



discontinuation_overall_axis_max <- plyr::round_any(ceiling(max(c(predictions_discontinuation_stan_psm_1_1_adjusted_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() * 100,
                                                                  predictions_discontinuation_stan_psm_1_1_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() * 100))), 5, f = ceiling)

discontinuation_strata_axis_max <- plyr::round_any(ceiling(max(c(predictions_discontinuation_stan_psm_1_1_adjusted_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() * 100,
                                                                 predictions_discontinuation_stan_psm_1_1_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() * 100))), 5, f = ceiling)


plot_c <- rbind(
  cbind(mean = 150, lci = 150, uci = 150, drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(mean = 150, lci = 150, uci = 150, drugclass = "GLP1", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_discontinuation_stan_psm_1_1_adjusted_overall %>%
    slice(c(1:6)) %>%
    mutate(mean = as.numeric(mean)*100,
           lci = as.numeric(lci)*100,
           uci = as.numeric(uci)*100),
  cbind(mean = 150, lci = 150, uci = 150, drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(mean = 150, lci = 150, uci = 150, drugclass = "GLP1", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_discontinuation_stan_psm_1_1_adjusted_overall %>%
    slice(-c(1:6)) %>%
    mutate(mean = as.numeric(mean)*100,
           lci = as.numeric(lci)*100,
           uci = as.numeric(uci)*100),
  predictions_discontinuation_stan_psm_1_1_adjusted_full %>% cbind(intervals = "Average treatment effect") %>%
    mutate(mean = as.numeric(mean)*100,
           lci = as.numeric(lci)*100,
           uci = as.numeric(uci)*100)
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == levels(group.discontinuation.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", format(group.discontinuation.dataset.matched%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[1])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.discontinuation.dataset.matched%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[1])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                            ifelse(intervals == levels(group.discontinuation.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", format(group.discontinuation.dataset.matched%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[2])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.discontinuation.dataset.matched%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[2])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                   ifelse(intervals == levels(group.discontinuation.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", format(group.discontinuation.dataset.matched%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[3])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.discontinuation.dataset.matched%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[3])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                          ifelse(intervals == levels(group.discontinuation.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", format(group.discontinuation.dataset.matched%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[4])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.discontinuation.dataset.matched%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[4])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                 ifelse(intervals == levels(group.discontinuation.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", format(group.discontinuation.dataset.matched%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[5])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.discontinuation.dataset.matched%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[5])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                        ifelse(intervals == levels(group.discontinuation.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", format(group.discontinuation.dataset.matched%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[6])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.discontinuation.dataset.matched%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[6])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"), intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  forestplot(ci.vertices = TRUE,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices.height = 0.1,
             title = "Risk of discontinuation",
             clip = c(0, discontinuation_overall_axis_max),
             xticks = seq(0, discontinuation_overall_axis_max, 5),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Discontinuation (%)") %>%
  fp_add_header(paste0("Overall population (n=", format(group.discontinuation.dataset.matched%>%nrow(),big.mark=",",scientific=FALSE), ")")) %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE)) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = "black")


### D - Risk of developing microvascular complications


patient_prop_scores_qrisk <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/patient_prop_scores_qrisk.rds")

micro_comp.dataset <- set_up_data_sglt2_glp1(dataset.type="micro_comp.dataset") %>%
  left_join(patient_prop_scores_qrisk, by = c("patid", "pated")) %>%
  left_join(treatment_effects, by = c("patid", "pated"))

group.micro_comp.dataset <- group_values(data = micro_comp.dataset,
                                         variable = "effects",
                                         breaks = interval_breaks) %>%
  drop_na(intervals)

matching_micro_comp <- MatchIt::matchit(
  formula = formula(paste0("drugclass ~  agetx + t2dmduration + prehba1c + preegfr + prealt + drugline + ncurrtx + sex + preneuropathy + preretinopathy")),
  data = group.micro_comp.dataset,
  method = "nearest",
  distance = group.micro_comp.dataset[,"prop.score"],
  replace = FALSE,
  m.order = "largest",
  caliper = 0.05,
  mahvars = NULL, estimand = "ATT", exact = NULL, antiexact = NULL, discard = "none", reestimate = FALSE, s.weights = NULL, std.caliper = TRUE, ratio = 1, verbose = FALSE, include.obj = FALSE,
)


# require(cobalt)
# cobalt::love.plot(matching_micro_comp, binary = "std", thresholds = c(m = .1), sample.names = c("Unmatched", "Matched"), title = "Full dataset")

n_drug <- 2*sum(!is.na(matching_micro_comp$match.matrix))

group.micro_comp.dataset.matched <- group.micro_comp.dataset %>%
  slice(which(matching_micro_comp$weights == 1))


# predictions for the CVD outcomes in the population with no CVD/HF/CKD
predictions_micro_comp_micro_comp_stan_psm_1_1_adjusted_overall <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_micro_comp_micro_comp_stan_psm_1_1_adjusted_overall.rds")

# predictions for the CVD outcomes in the population with no CVD/HF/CKD
predictions_micro_comp_micro_comp_stan_psm_1_1_adjusted_full <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_micro_comp_micro_comp_stan_psm_1_1_adjusted_full.rds")

# plot of forest plot
plot_d <- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_micro_comp_micro_comp_stan_psm_1_1_adjusted_overall %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = NA, lci = NA, uci = NA),
  predictions_micro_comp_micro_comp_stan_psm_1_1_adjusted_overall %>% slice(4:6),
  predictions_micro_comp_micro_comp_stan_psm_1_1_adjusted_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.micro_comp.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", format(group.micro_comp.dataset.matched%>%filter(intervals==levels(group.micro_comp.dataset.matched$intervals)[1])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.micro_comp.dataset.matched%>%filter(intervals==levels(group.micro_comp.dataset.matched$intervals)[1])%>%filter(postdrug_micro_comp_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                            ifelse(intervals == levels(group.micro_comp.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", format(group.micro_comp.dataset.matched%>%filter(intervals==levels(group.micro_comp.dataset.matched$intervals)[2])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.micro_comp.dataset.matched%>%filter(intervals==levels(group.micro_comp.dataset.matched$intervals)[2])%>%filter(postdrug_micro_comp_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                   ifelse(intervals == levels(group.micro_comp.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", format(group.micro_comp.dataset.matched%>%filter(intervals==levels(group.micro_comp.dataset.matched$intervals)[3])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.micro_comp.dataset.matched%>%filter(intervals==levels(group.micro_comp.dataset.matched$intervals)[3])%>%filter(postdrug_micro_comp_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                          ifelse(intervals == levels(group.micro_comp.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", format(group.micro_comp.dataset.matched%>%filter(intervals==levels(group.micro_comp.dataset.matched$intervals)[4])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.micro_comp.dataset.matched%>%filter(intervals==levels(group.micro_comp.dataset.matched$intervals)[4])%>%filter(postdrug_micro_comp_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                 ifelse(intervals == levels(group.micro_comp.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", format(group.micro_comp.dataset.matched%>%filter(intervals==levels(group.micro_comp.dataset.matched$intervals)[5])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.micro_comp.dataset.matched%>%filter(intervals==levels(group.micro_comp.dataset.matched$intervals)[5])%>%filter(postdrug_micro_comp_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                        ifelse(intervals == levels(group.micro_comp.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", format(group.micro_comp.dataset.matched%>%filter(intervals==levels(group.micro_comp.dataset.matched$intervals)[6])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.micro_comp.dataset.matched%>%filter(intervals==levels(group.micro_comp.dataset.matched$intervals)[6])%>%filter(postdrug_micro_comp_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci") %>%
  mutate(mean = exp(mean),
         lower = exp(lower),
         upper = exp(upper)) %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             zero = 1,
             title = "Microvascular complications",
             xlog = TRUE,
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazard ratio (95% CI, log scale)") %>%
  fp_add_header(paste0("Overall population (n=", format(nrow(group.micro_comp.dataset.matched),big.mark=",",scientific=FALSE), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = "black")


### E - Risk of developing MACE/CVD


patient_prop_scores_qrisk <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/patient_prop_scores_qrisk.rds")

no_co.dataset <- set_up_data_sglt2_glp1(dataset.type="no_co.dataset") %>%
  left_join(patient_prop_scores_qrisk, by = c("patid", "pated")) %>%
  left_join(treatment_effects, by = c("patid", "pated"))

group.no_co.dataset <- group_values(data = no_co.dataset,
                                    variable = "effects",
                                    breaks = interval_breaks) %>%
  drop_na(intervals)

breakdown_adjust <- unique(c(variables_mu, variables_tau))
# categorical variables in breakdown
factors <- sapply(group.no_co.dataset[,breakdown_adjust], is.factor)

matching_no_co <- MatchIt::matchit(
  formula = formula(paste0("drugclass ~  agetx + t2dmduration + prehba1c + preegfr + prealt + drugline + ncurrtx + sex + preneuropathy + preretinopathy")),
  data = group.no_co.dataset,
  method = "nearest",
  distance = group.no_co.dataset[,"prop.score"],
  replace = FALSE,
  m.order = "largest",
  caliper = 0.05,
  mahvars = NULL, estimand = "ATT", exact = NULL, antiexact = NULL, discard = "none", reestimate = FALSE, s.weights = NULL, std.caliper = TRUE, ratio = 1, verbose = FALSE, include.obj = FALSE,
)

# require(cobalt)
# cobalt::love.plot(matching_no_co, binary = "std", thresholds = c(m = .1), sample.names = c("Unmatched", "Matched"), title = "Full dataset")

n_drug <- 2*sum(!is.na(matching_no_co$match.matrix))

group.no_co.dataset.matched <- group.no_co.dataset %>%
  slice(which(matching_no_co$weights == 1))


# predictions for the CVD outcomes in the population with no CVD/HF/CKD
predictions_no_co_cvd_stan_psm_1_1_adjusted_overall <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_no_co_cvd_stan_psm_1_1_adjusted_overall.rds")

# predictions for the CVD outcomes in the population with no CVD/HF/CKD
predictions_no_co_cvd_stan_psm_1_1_adjusted_full <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_no_co_cvd_stan_psm_1_1_adjusted_full.rds")


# plot of forest plot
plot_e <- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_no_co_cvd_stan_psm_1_1_adjusted_overall %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = NA, lci = NA, uci = NA),
  predictions_no_co_cvd_stan_psm_1_1_adjusted_overall %>% slice(4:6),
  predictions_no_co_cvd_stan_psm_1_1_adjusted_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[1])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[1])%>%filter(postdrug_mace_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                            ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[2])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[2])%>%filter(postdrug_mace_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                   ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[3])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[3])%>%filter(postdrug_mace_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                          ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[4])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[4])%>%filter(postdrug_mace_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                 ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[5])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[5])%>%filter(postdrug_mace_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                        ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[6])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[6])%>%filter(postdrug_mace_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci") %>%
  mutate(mean = exp(mean),
         lower = exp(lower),
         upper = exp(upper)) %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             zero = 1,
             title = "Major adverse cardiovascular events",
             xlog = TRUE,
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazard ratio (95% CI, log scale)") %>%
  fp_add_header(paste0("Overall population (n=", format(nrow(group.no_co.dataset.matched),big.mark=",",scientific=FALSE), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = "black")


### F - Risk of developing Heart Failure


# predictions for the HF outcomes in the population with no CVD/HF/CKD
predictions_no_co_hf_stan_psm_1_1_adjusted_overall <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_no_co_hf_stan_psm_1_1_adjusted_overall.rds")

# predictions for the CVD outcomes in the population with no CVD/HF/CKD
predictions_no_co_hf_stan_psm_1_1_adjusted_full <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_no_co_hf_stan_psm_1_1_adjusted_full.rds")

plot_f <- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_no_co_hf_stan_psm_1_1_adjusted_overall %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = NA, lci = NA, uci = NA),
  predictions_no_co_hf_stan_psm_1_1_adjusted_overall %>% slice(4:6),
  predictions_no_co_hf_stan_psm_1_1_adjusted_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[1])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[1])%>%filter(postdrug_hf_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                            ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[2])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[2])%>%filter(postdrug_hf_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                   ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[3])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[3])%>%filter(postdrug_hf_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                          ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[4])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[4])%>%filter(postdrug_hf_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                 ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[5])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[5])%>%filter(postdrug_hf_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                        ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[6])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[6])%>%filter(postdrug_hf_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci") %>%
  mutate(mean = exp(mean),
         lower = exp(lower),
         upper = exp(upper)) %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             zero = 1,
             title = "Heart failure",
             xlog = TRUE,
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazard ratio (95% CI, log scale)") %>%
  fp_add_header(paste0("Overall population (n=", format(nrow(group.no_co.dataset.matched),big.mark=",",scientific=FALSE), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = "black")


pdf(width = 20, height = 10, "Plots/Paper/Sup_Mat/11.08.plot_sup_12.pdf")

grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 2,
                                           ncol = 3, heights = unit(c(5, 5), "null"))))
# 1 2 3
# 4 5 6
# first plot
pushViewport(viewport(layout.pos.row = 1,
                      layout.pos.col = 1))
plot_a
upViewport()

# second plot
pushViewport(viewport(layout.pos.row = 1,
                      layout.pos.col = 2))
plot_b
upViewport()

# third plot
pushViewport(viewport(layout.pos.row = 1,
                      layout.pos.col = 3))
plot_c
upViewport()

# forth plot
pushViewport(viewport(layout.pos.row = 2,
                      layout.pos.col = 1))
plot_d
upViewport()

# fifth plot
pushViewport(viewport(layout.pos.row = 2,
                      layout.pos.col = 2))
plot_e
upViewport()

# third plot
pushViewport(viewport(layout.pos.row = 2,
                      layout.pos.col = 3))
plot_f
upViewport()

# add legends
grid.text(expression(bold("A.1")), vp = viewport(layout.pos.row = 1, layout.pos.col = 1), hjust = 12, vjust = -18)
grid.text(expression(bold("A.2")), vp = viewport(layout.pos.row = 1, layout.pos.col = 2), hjust = 12, vjust = -18)
grid.text(expression(bold("A.3")), vp = viewport(layout.pos.row = 1, layout.pos.col = 3), hjust = 12, vjust = -17.5)
grid.text(expression(bold("B.1")), vp = viewport(layout.pos.row = 2, layout.pos.col = 1), hjust = 12, vjust = -18)
grid.text(expression(bold("B.2")), vp = viewport(layout.pos.row = 2, layout.pos.col = 2), hjust = 12, vjust = -18)
grid.text(expression(bold("B.3")), vp = viewport(layout.pos.row = 2, layout.pos.col = 3), hjust = 12, vjust = -17.5)

dev.off()



#:--------------------------------------------------------------------------------
### Validation hold-out cohort with/without CVD calibration

variables_mu <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/variables_mu.rds")

variables_tau <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/variables_tau.rds")

hba1c.test <- set_up_data_sglt2_glp1(dataset.type = "hba1c.test") %>%
  left_join(readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/patient_effects.rds"), by = c("patid", "pated")) %>%
  left_join(readRDS("Samples/SGLT2-GLP1/Aurum/ps_model/patient_prop_scores.rds"), by = c("patid", "pated"))

hba1c.test.with.cvd <- hba1c.test %>%
  filter(preangina == "Yes" | preaf == "Yes" | prerevasc == "Yes" | preheartfailure == "Yes" | prehypertension == "Yes" | preihd == "Yes" | premyocardialinfarction == "Yes" | prestroke == "Yes" | pretia == "Yes") %>%
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


plot_ATE_psm_1_1_hba1c.test.with.cvd <- ATE_plot(ATE_psm_1_1_hba1c.test.with.cvd[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -10, 10) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(
    title = paste0("Cohort with CVD (n=", format(sum(ATE_psm_1_1_hba1c.test.with.cvd[["effects"]]$n_drug1)*2,big.mark=",",scientific=FALSE), ")"),
    subtitle = "Propensity score matching"
  )

plot_ATE_psm_1_1_adjusted_hba1c.test.with.cvd <- ATE_plot(ATE_psm_1_1_adjusted_hba1c.test.with.cvd[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -10, 10) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(
    title = paste0("Cohort with CVD (n=", format(sum(ATE_psm_1_1_adjusted_hba1c.test.with.cvd[["effects"]]$n_drug1)*2,big.mark=",",scientific=FALSE), ")"),
    subtitle = "Propensity score matching and estimate adjustment"
  )

plot_ATE_adjusted_hba1c.test.with.cvd <- ATE_plot(ATE_adjusted_hba1c.test.with.cvd[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -10, 10) +
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


plot_ATE_psm_1_1_hba1c.test.without.cvd <- ATE_plot(ATE_psm_1_1_hba1c.test.without.cvd[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -10, 10) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(
    title = paste0("Cohort without CVD (n=", format(sum(ATE_psm_1_1_hba1c.test.without.cvd[["effects"]]$n_drug1)*2,big.mark=",",scientific=FALSE), ")"),
    subtitle = "Propensity score matching"
  )

plot_ATE_psm_1_1_adjusted_hba1c.test.without.cvd <- ATE_plot(ATE_psm_1_1_adjusted_hba1c.test.without.cvd[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -10, 10) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(
    title = paste0("Cohort without CVD (n=", format(sum(ATE_psm_1_1_adjusted_hba1c.test.without.cvd[["effects"]]$n_drug1)*2,big.mark=",",scientific=FALSE), ")"),
    subtitle = "Propensity score matching and estimate adjustment"
  )

plot_ATE_adjusted_hba1c.test.without.cvd <- ATE_plot(ATE_adjusted_hba1c.test.without.cvd[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -10, 10) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(
    title = paste0("Cohort without CVD (n=", format(sum(ATE_adjusted_hba1c.test.without.cvd[["effects"]]$N),big.mark=",",scientific=FALSE), ")"),
    subtitle = "Estimate adjustment"
  )

plot_sup_13 <- (plot_ATE_psm_1_1_hba1c.test.with.cvd | plot_ATE_psm_1_1_hba1c.test.without.cvd) / (plot_ATE_psm_1_1_adjusted_hba1c.test.with.cvd | plot_ATE_psm_1_1_adjusted_hba1c.test.without.cvd) / (plot_ATE_adjusted_hba1c.test.with.cvd | plot_ATE_adjusted_hba1c.test.without.cvd) +
  plot_annotation(tag_levels = list(c("A.1", "A.2", "B.1", "B.2", "C.1", "C.2")))


pdf(width = 10, height = 15, "Plots/Paper/Sup_Mat/11.08.plot_sup_13.pdf")
plot_sup_13
dev.off()




#:--------------------------------------------------------------------------------
### Variables ranges across HbA1c benefit subgroups 

interval_breaks <- c(-5, -3, 0, 3, 5)

# Full cohort for average values
full.cohort <- set_up_data_sglt2_glp1(dataset.type="full.cohort") %>%
  left_join(readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/patient_effects.rds"), by = c("patid", "pated"))

group.full.cohort <- group_values(data = full.cohort,
                                              variable = "effects",
                                              breaks = interval_breaks) %>%
  drop_na(intervals) 

group.full.cohort <- group.full.cohort %>%
  mutate(intervals = ifelse(intervals == levels(group.full.cohort$intervals)[1], "SGLT2i benefit >5 mmol/mol",
                            ifelse(intervals == levels(group.full.cohort$intervals)[2], "SGLT2i benefit 3-5 mmol/mol",
                                   ifelse(intervals == levels(group.full.cohort$intervals)[3], "SGLT2i benefit 0-3 mmol/mol",
                                          ifelse(intervals == levels(group.full.cohort$intervals)[4], "GLP1-RA benefit 0-3 mmol/mol",
                                                 ifelse(intervals == levels(group.full.cohort$intervals)[5], "GLP1-RA benefit 3-5 mmol/mol",
                                                        ifelse(intervals == levels(group.full.cohort$intervals)[6], "GLP1-RA benefit >5 mmol/mol", intervals))))))
  ) %>%
  mutate(intervals = factor(intervals, levels = c("SGLT2i benefit >5 mmol/mol", "SGLT2i benefit 3-5 mmol/mol", "SGLT2i benefit 0-3 mmol/mol", "GLP1-RA benefit 0-3 mmol/mol", "GLP1-RA benefit 3-5 mmol/mol", "GLP1-RA benefit >5 mmol/mol")))

levels(group.full.cohort$sex) <- c("Females", "Males")
levels(group.full.cohort$ncurrtx) <- c("0", "1", "2", "3", "4+")

# number of other current drugs
plot_ncurrtx_strata <- group.full.cohort %>%
  select(ncurrtx, intervals) %>%
  ggplot(aes(x = intervals)) +
  geom_bar(aes(fill = ncurrtx), position = "fill") +
  ggtitle("Number of other current\nglucose-lowering drugs") +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 70, vjust = 1, hjust = 1),
        strip.background = element_rect(fill = "white"),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))

plot_sex_strata <- group.full.cohort %>%
  select(sex, intervals) %>%
  ggplot(aes(x = intervals)) +
  geom_bar(aes(fill = sex), position = "fill") +
  ggtitle("Sex") +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 70, vjust = 1, hjust = 1),
        strip.background = element_rect(fill = "white"),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))

plot_prehba1c_strata <- group.full.cohort %>%
  select(intervals, prehba1c) %>%
  drop_na() %>%
  ggplot(aes(x = intervals, y = prehba1c)) +
  geom_boxplot(outlier.shape = NA, fill = c(rep("#f1a340", 3), rep("dodgerblue2", 3)), alpha = c(1, 0.6, 0.2, 0.2, 0.6, 1)) +
  ggtitle("HbA1c") +
  ylim(25, 150) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 70, vjust = 1, hjust = 1),
        strip.background = element_rect(fill = "white"),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))

plot_preegfr_strata <- group.full.cohort %>%
  select(intervals, preegfr) %>%
  drop_na() %>%
  ggplot(aes(x = intervals, y = preegfr)) +
  geom_boxplot(outlier.shape = NA, fill = c(rep("#f1a340", 3), rep("dodgerblue2", 3)), alpha = c(1, 0.6, 0.2, 0.2, 0.6, 1)) +
  ggtitle("eGFR") +
  ylim(20, 140) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 70, vjust = 1, hjust = 1),
        strip.background = element_rect(fill = "white"),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))

plot_agetx_strata <- group.full.cohort %>%
  select(intervals, agetx) %>%
  drop_na() %>%
  ggplot(aes(x = intervals, y = agetx)) +
  geom_boxplot(outlier.shape = NA, fill = c(rep("#f1a340", 3), rep("dodgerblue2", 3)), alpha = c(1, 0.6, 0.2, 0.2, 0.6, 1)) +
  ggtitle("Current age") +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 70, vjust = 1, hjust = 1),
        strip.background = element_rect(fill = "white"),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))

plot_prebmi_strata <- group.full.cohort %>%
  select(intervals, prebmi) %>%
  drop_na() %>%
  ggplot(aes(x = intervals, y = prebmi)) +
  geom_boxplot(outlier.shape = NA, fill = c(rep("#f1a340", 3), rep("dodgerblue2", 3)), alpha = c(1, 0.6, 0.2, 0.2, 0.6, 1)) +
  ggtitle("BMI") +
  ylim(12.5, 62.5) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 70, vjust = 1, hjust = 1),
        strip.background = element_rect(fill = "white"),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))
  

plot_sup_14 <- patchwork::wrap_plots(list(plot_ncurrtx_strata, plot_sex_strata, plot_prehba1c_strata, plot_preegfr_strata, plot_agetx_strata, plot_prebmi_strata), nrow = 1)


pdf(width = 17, height = 8, "Plots/Paper/Sup_Mat/11.08.plot_sup_14.pdf")
plot_sup_14
dev.off()




