####################
## Description:
##  - In this file we create the plots used in the Diabetes UK conference 2023
####################


## Load libraries
library(tidyverse)
library(patchwork)
library(scales)

## Load functions required
source("01.slade_aurum_functions.R")
source("02.slade_aurum_set_data.R")


## make directory for outputs
dir.create("Plots/DUK_2023")


#### Presentation plots - Diabetes UK


#:--------------------------------------------------------------------------------
# Histogram of predicted treatment effects for overall population

hba1c.train <- set_up_data_sglt2_glp1(dataset.type = "hba1c.train") %>%
  # drop the variables with the most missingness (>40%)
  select(-preacr, -preast, -prehaematocrit, -prehaemoglobin, -pretriglyceride) %>%
  # only complete cases
  drop_na() %>%
  left_join(readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/patient_effects.rds"), by = c("patid", "pated"))

# plot histogram of treatment effects
plot_duk_1 <- hist_plot(data = hba1c.train %>%
                          rename("mean" = "effects"),
                        title = paste0("Development cohort (n=", format(nrow(hba1c.train),big.mark=",",scientific=FALSE), ")"),
                        xmin = -15,
                        xmax = 15) +
  scale_y_continuous(label=comma) +
  scale_x_continuous(breaks = c(-15, -10, -5, 0, 5, 10, 15)) +
  geom_vline(aes(xintercept = -5), lty = "dashed", colour = "red") +
  geom_vline(aes(xintercept = 5), lty = "dashed", colour = "red") +
  theme(
    plot.title = element_text(size = 24),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 17),
    legend.text = element_text(size = 15)
    
  )

# PDF with plot
pdf(width = 7, height = 6, "Plots/DUK_2023/11.08.plot_duk_1.pdf")
plot_duk_1
dev.off()



# :-------------------------------------------------------------------------------
# Calibration plots for development and validation cohorts

# # load development cohort model calibration
# ATE_adjust_validation_dev <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/assessment/ATE_adjust_validation_dev.rds")
# 
# # plot development cohort calibration
# plot_ATE_adjust_validation_dev <- ATE_plot(ATE_adjust_validation_dev[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -12, 12) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   ggtitle(paste0("Development cohort (n=",format(sum(ATE_adjust_validation_dev[["effects"]]$N),big.mark=",",scientific=FALSE), ")"))

# load validation cohort model calibration
ATE_adjust_validation_val <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/assessment/ATE_adjust_validation_val.rds")

# plot validation cohort calibration
plot_ATE_adjust_validation_val <- ATE_plot(ATE_adjust_validation_val[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -12, 12) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle(paste0("Hold-out CPRD cohort (n=",format(sum(ATE_adjust_validation_val[["effects"]]$N),big.mark=",",scientific=FALSE), ")"))

# # Combine plots together
# plot_duk_2 <- patchwork::wrap_plots(list(plot_ATE_adjust_validation_dev, plot_ATE_adjust_validation_val), ncol = 2, nrow = 1) &
#   theme(
#     plot.title = element_text(size = 15),
#     axis.title = element_text(size = 11),
#     axis.text.x = element_text(size = 12),
#     axis.text.y = element_text(size = 13)
#   )

# PDF with plot
pdf(width = 4.8, height = 5, "Plots/DUK_2023/11.08.plot_duk_2.pdf")
plot_ATE_adjust_validation_val +
  xlab("Predicted treatment effect") +
  ylab("Average treatment effect") +
  theme(
    plot.title = element_text(size = 19),
    axis.title = element_text(size = 16),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 13)
  )
dev.off()



#:--------------------------------------------------------------------------------
# Forest plot of predicted treatment effect benefit vs average treatment effect

## Load libraries
require(forestplot)

interval_breaks <- c(-5, -3, 0, 3, 5)

## Read in propensity scores
patient_prop_scores <- readRDS("Samples/SGLT2-GLP1/Aurum/ps_model/patient_prop_scores.rds")

## Read in treatment effects
treatment_effects <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/patient_effects.rds")

# load in variables used in the model
variables_mu <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/variables_mu.rds")

variables_tau <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/variables_tau.rds")

model_variables <- unique(c(variables_mu, variables_tau))[which(unique(c(variables_mu, variables_tau)) != "sex")]


# Overall population
hba1c <- set_up_data_sglt2_glp1(dataset.type = "hba1c.train") %>%
  # drop the variables with the most missingness (>40%)
  select(-preacr, -preast, -prehaematocrit, -prehaemoglobin, -pretriglyceride) %>%
  # only complete cases
  drop_na() %>%
  left_join(patient_prop_scores, by = c("patid", "pated")) %>%
  left_join(treatment_effects, by = c("patid", "pated"))

# split population into groups
group.hba1c.dataset <- group_values(data = hba1c,
                                    variable = "effects",
                                    breaks = interval_breaks) %>%
  drop_na(intervals) %>%
  rename("hba1c_diff" = "effects")

# calculate average treatment effect for the subgroups
ATE_adjust_hba1c <- calc_ATE(data = group.hba1c.dataset, validation_type = "Adjust", variable = "posthba1cfinal",
                             quantile_var = "intervals", breakdown = model_variables)

ATE_adjust_hba1c_full <- calc_ATE(data = group.hba1c.dataset %>% mutate(intervals = as.numeric(1)), validation_type = "Adjust", variable = "posthba1cfinal",
                                  quantile_var = "intervals", breakdown = model_variables)

# calculate axis limits for the plots
hba1c_strata_axis_min <- plyr::round_any(floor(min(c(ATE_adjust_hba1c[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
                                                     ATE_adjust_hba1c_full[["effects"]] %>% select(c("obs","lci","uci")) %>% min()))), 2, f = floor)

hba1c_strata_axis_max <- plyr::round_any(ceiling(max(c(ATE_adjust_hba1c[["effects"]] %>% select(c("obs","lci","uci")) %>% max(),
                                                       ATE_adjust_hba1c_full[["effects"]] %>% select(c("obs","lci","uci")) %>% max()))), 2, f = ceiling)

# plot treatment effects benefit
plot_duk_3 <- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", obs = NA, lci = NA, uci = NA),
  ATE_adjust_hba1c[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP-1 RA", obs = NA, lci = NA, uci = NA),
  ATE_adjust_hba1c[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(4:6),
  ATE_adjust_hba1c_full[["effects"]] %>% select(intervals, obs, lci, uci) %>% mutate(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(obs = as.numeric(obs),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.hba1c.dataset$intervals)[1], paste0(">5 mmol/mol (n=", format(group.hba1c.dataset%>%filter(intervals == levels(group.hba1c.dataset$intervals)[1])%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                            ifelse(intervals == levels(group.hba1c.dataset$intervals)[2], paste0("3-5 mmol/mol (n=", format(group.hba1c.dataset%>%filter(intervals == levels(group.hba1c.dataset$intervals)[2])%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                   ifelse(intervals == levels(group.hba1c.dataset$intervals)[3], paste0("0-3 mmol/mol (n=", format(group.hba1c.dataset%>%filter(intervals == levels(group.hba1c.dataset$intervals)[3])%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                          ifelse(intervals == levels(group.hba1c.dataset$intervals)[4], paste0("0-3 mmol/mol (n=", format(group.hba1c.dataset%>%filter(intervals == levels(group.hba1c.dataset$intervals)[4])%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                 ifelse(intervals == levels(group.hba1c.dataset$intervals)[5], paste0("3-5 mmol/mol (n=", format(group.hba1c.dataset%>%filter(intervals == levels(group.hba1c.dataset$intervals)[5])%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                        ifelse(intervals == levels(group.hba1c.dataset$intervals)[6], paste0(">5 mmol/mol (n=", format(group.hba1c.dataset%>%filter(intervals == levels(group.hba1c.dataset$intervals)[6])%>%nrow(),big.mark=",",scientific=FALSE), ")"), intervals))))))) %>%
  rename("mean" = "obs", "lower" = "lci", "upper" = "uci") %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             xticks = seq(hba1c_strata_axis_min, hba1c_strata_axis_max, 2),
             ci.vertices.height = 0.3,
             boxsize = .4,
             lwd.ci = 3,
             txt_gp = fpTxtGp(ticks=gpar(cex=1.5), 
                              title=gpar(cex=2.5),
                              legend =gpar(cex=2),
                              label=gpar(cex=1.5),
                              xlab=gpar(cex=1.5)),
             xlab = "Average HbA1c treatment effect (mmol/mol)") %>%
  fp_add_header(paste0("Development cohort (n=", format(group.hba1c.dataset%>%nrow(),big.mark=",",scientific=FALSE), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = "black")

# PDF with plot
pdf(width = 10.5, height = 6, "Plots/DUK_2023/11.08.plot_duk_3.pdf")
plot_duk_3
dev.off()




#:--------------------------------------------------------------------------------
# Forest plot of predicted treatment effect benefit vs predicted hba1c change

interval_breaks <- c(-5, -3, 0, 3, 5)

## Read in propensity scores
patient_prop_scores <- readRDS("Samples/SGLT2-GLP1/Aurum/ps_model/patient_prop_scores.rds")

## Read in treatment effects
treatment_effects <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/patient_effects.rds")

# Overall population
hba1c <- set_up_data_sglt2_glp1(dataset.type = "hba1c.train") %>%
  # drop the variables with the most missingness (>40%)
  select(-preacr, -preast, -prehaematocrit, -prehaemoglobin, -pretriglyceride) %>%
  # only complete cases
  drop_na() %>%
  left_join(patient_prop_scores, by = c("patid", "pated")) %>%
  left_join(treatment_effects, by = c("patid", "pated"))

# split population into groups
group.hba1c.dataset <- group_values(data = hba1c,
                                    variable = "effects",
                                    breaks = interval_breaks) %>%
  drop_na(intervals) %>%
  rename("hba1c_diff" = "effects")

# predictions for the hba1c adjusted model
predictions_hba1c_stan_adjusted_overall <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_hba1c_stan_adjusted_overall.rds") %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci))

# predictions for the hba1c adjusted model full
predictions_hba1c_stan_adjusted_full <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_hba1c_stan_adjusted_full.rds") %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci))


# calculate axis limits for the plots
hba1c_strata_axis_min <- plyr::round_any(floor(min(c(predictions_hba1c_stan_adjusted_overall %>% select(c("mean","lci","uci")) %>% min(),
                                                     predictions_hba1c_stan_adjusted_full %>% select(c("mean","lci","uci")) %>% min()))), 5, f = floor)

hba1c_strata_axis_max <- plyr::round_any(floor(max(c(predictions_hba1c_stan_adjusted_overall %>% select(c("mean","lci","uci")) %>% max(),
                                                     predictions_hba1c_stan_adjusted_full %>% select(c("mean","lci","uci")) %>% max()))), 5, f = ceiling)


# plot predicted hba1c change
plot_duk_4 <- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = 50, lci = 50, uci = 50, drugclass = "SGLT2"),
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = 50, lci = 50, uci = 50, drugclass = "GLP1"),
  predictions_hba1c_stan_adjusted_overall %>% select(intervals, drugclass, mean, lci, uci) %>% slice(1:6),
  cbind(intervals = "Predicted HbA1c benefit on GLP-1 RA", mean = 50, lci = 50, uci = 50, drugclass = "SGLT2"),
  cbind(intervals = "Predicted HbA1c benefit on GLP-1 RA", mean = 50, lci = 50, uci = 50, drugclass = "GLP1"),
  predictions_hba1c_stan_adjusted_overall %>% select(intervals, drugclass, mean, lci, uci) %>% slice(7:12),
  predictions_hba1c_stan_adjusted_full %>% select(drugclass, mean, lci, uci) %>% mutate(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == "(-21,-5]", paste0(">5 mmol/mol (n=", format(group.hba1c.dataset%>%filter(intervals == levels(group.hba1c.dataset$intervals)[1])%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                            ifelse(intervals == levels(group.hba1c.dataset$intervals)[2], paste0("3-5 mmol/mol (n=", format(group.hba1c.dataset%>%filter(intervals == levels(group.hba1c.dataset$intervals)[2])%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                   ifelse(intervals == levels(group.hba1c.dataset$intervals)[3], paste0("0-3 mmol/mol (n=", format(group.hba1c.dataset%>%filter(intervals == levels(group.hba1c.dataset$intervals)[3])%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                          ifelse(intervals == levels(group.hba1c.dataset$intervals)[4], paste0("0-3 mmol/mol (n=", format(group.hba1c.dataset%>%filter(intervals == levels(group.hba1c.dataset$intervals)[4])%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                 ifelse(intervals == levels(group.hba1c.dataset$intervals)[5], paste0("3-5 mmol/mol (n=", format(group.hba1c.dataset%>%filter(intervals == levels(group.hba1c.dataset$intervals)[5])%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                        ifelse(intervals == "(5,40]", paste0(">5 mmol/mol (n=", format(group.hba1c.dataset%>%filter(intervals == levels(group.hba1c.dataset$intervals)[6])%>%nrow(),big.mark=",",scientific=FALSE), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  mutate(group = ifelse(group == "SGLT2", "SGLT2i", "GLP-1 RA")) %>%
  forestplot(ci.vertices = TRUE,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices.height = 0.1,
             clip = c(hba1c_strata_axis_min, hba1c_strata_axis_max),
             xticks = seq(hba1c_strata_axis_min, hba1c_strata_axis_max, 5),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Predicted HbA1c change (mmol/mol)") %>%
  fp_add_header(paste0("Development cohort (n=", format(group.hba1c.dataset%>%nrow(),big.mark=",",scientific=FALSE), ")")) %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE)) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = "black")

# PDF with plot
pdf(width = 7, height = 5, "Plots/DUK_2023/11.08.plot_duk_4.pdf")
plot_duk_4
dev.off()



#:--------------------------------------------------------------------------------
# Forest plot of predicted treatment effect benefit vs weight change
require(forestplot)

interval_breaks <- c(-5, -3, 0, 3, 5)

## Read in propensity scores
patient_prop_scores <- readRDS("Samples/SGLT2-GLP1/Aurum/ps_model/patient_prop_scores.rds")

## Read in treatment effects
treatment_effects <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/patient_effects.rds")

# Weight change population
weight.dataset <- set_up_data_sglt2_glp1(dataset.type = "weight.dataset") %>%
  left_join(patient_prop_scores, by = c("patid", "pated")) %>%
  left_join(treatment_effects, by = c("patid", "pated")) %>%
  mutate(w.change = postweight - preweight)

# Split population into subgroups
group.weight.dataset <- group_values(data = weight.dataset,
                                     variable = "effects",
                                     breaks = interval_breaks) %>%
  drop_na(intervals)

# predictions for the weight adjusted model
predictions_weight_stan_adjusted_overall <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_weight_stan_adjusted_overall.rds")

# predictions for the weight adjusted model full
predictions_weight_stan_adjusted_full <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_weight_stan_adjusted_full.rds")

# limits for the plot
weight_strata_axis_min <- plyr::round_any(floor(min(c(predictions_weight_stan_adjusted_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min(),
                                                      predictions_weight_stan_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min()))), 1, f = floor)

weight_strata_axis_max <- plyr::round_any(ceiling(max(c(predictions_weight_stan_adjusted_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max(),
                                                        predictions_weight_stan_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max()))), 2, f = ceiling)

# forest plot
plot_duk_5 <- rbind(
  cbind(mean = 50, lci = 50, uci = 50, drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(mean = 50, lci = 50, uci = 50, drugclass = "GLP1", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_weight_stan_adjusted_overall %>%
    slice(c(1:6)),
  cbind(mean = 50, lci = 50, uci = 50, drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on GLP-1 RA"),
  cbind(mean = 50, lci = 50, uci = 50, drugclass = "GLP1", intervals = "Predicted HbA1c benefit on GLP-1 RA"),
  predictions_weight_stan_adjusted_overall %>%
    slice(-c(1:6)),
  predictions_weight_stan_adjusted_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == levels(group.weight.dataset$intervals)[1], paste0(">5 mmol/mol (n=", format(group.weight.dataset%>%filter(intervals==levels(group.weight.dataset$intervals)[1])%>%nrow(),big.mark=",",scientific=FALSE),")"),
                            ifelse(intervals == levels(group.weight.dataset$intervals)[2], paste0("3-5 mmol/mol (n=", format(group.weight.dataset%>%filter(intervals==levels(group.weight.dataset$intervals)[2])%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                   ifelse(intervals == levels(group.weight.dataset$intervals)[3], paste0("0-3 mmol/mol (n=", format(group.weight.dataset%>%filter(intervals==levels(group.weight.dataset$intervals)[3])%>%nrow(),big.mark=",",scientific=FALSE),")"),
                                          ifelse(intervals == levels(group.weight.dataset$intervals)[4], paste0("0-3 mmol/mol (n=", format(group.weight.dataset%>%filter(intervals==levels(group.weight.dataset$intervals)[4])%>%nrow(),big.mark=",",scientific=FALSE),")"),
                                                 ifelse(intervals == levels(group.weight.dataset$intervals)[5], paste0("3-5 mmol/mol (n=", format(group.weight.dataset%>%filter(intervals==levels(group.weight.dataset$intervals)[5])%>%nrow(),big.mark=",",scientific=FALSE),")"),
                                                        ifelse(intervals == levels(group.weight.dataset$intervals)[6], paste0(">5 mmol/mol (n=", format(group.weight.dataset%>%filter(intervals==levels(group.weight.dataset$intervals)[6])%>%nrow(),big.mark=",",scientific=FALSE),")"), intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  mutate(group = ifelse(group == "SGLT2", "SGLT2i", "GLP-1 RA")) %>%
  forestplot(ci.vertices = TRUE,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices.height = 0.1,
             # title = "Weight change",
             clip = c(weight_strata_axis_min, weight_strata_axis_max),
             xticks = seq(weight_strata_axis_min, weight_strata_axis_max, 1),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=1.5), 
                              title=gpar(cex=2.5),
                              legend =gpar(cex=2),
                              label=gpar(cex=1.5),
                              xlab=gpar(cex=1.5)),
             xlab = "Predicted weight change (kg)") %>%
  fp_add_header(paste0("Overall population (n=", format(group.weight.dataset%>%nrow(),big.mark=",",scientific=FALSE), ")")) %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE)) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = "black")

# PDF with plots
pdf(width = 11, height = 6, "Plots/DUK_2023/11.08.plot_duk_5.pdf")
plot_duk_5
dev.off()


#:--------------------------------------------------------------------------------
# Forest plot of predicted treatment effect benefit vs discontinuation
require(forestplot)

interval_breaks <- c(-5, -3, 0, 3, 5)

## Read in propensity scores
patient_prop_scores <- readRDS("Samples/SGLT2-GLP1/Aurum/ps_model/patient_prop_scores.rds")

## Read in treatment effects
treatment_effects <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/patient_effects.rds")

# Discontinuation
discontinuation.dataset <- set_up_data_sglt2_glp1(dataset.type = "discontinuation.dataset") %>%
  left_join(patient_prop_scores, by = c("patid", "pated")) %>%
  left_join(treatment_effects, by = c("patid", "pated")) %>%
  mutate(stopdrug_6m_3mFU = factor(stopdrug_6m_3mFU))

# Split population into subgroups
group.discontinuation.dataset <- group_values(data = discontinuation.dataset,
                                              variable = "effects",
                                              breaks = interval_breaks) %>%
  drop_na(intervals)

# predictions for the discontinuation adjusted model
predictions_discontinuation_stan_adjusted_overall <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_discontinuation_stan_adjusted_overall.rds")

# predictions for the discontinuation model sex strata
predictions_discontinuation_stan_adjusted_full <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_discontinuation_stan_adjusted_full.rds")

# limits for the plot
discontinuation_overall_axis_max <- plyr::round_any(ceiling(max(c(predictions_discontinuation_stan_adjusted_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() * 100,
                                                                  predictions_discontinuation_stan_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() * 100))), 5, f = ceiling)

discontinuation_strata_axis_max <- plyr::round_any(ceiling(max(c(predictions_discontinuation_stan_adjusted_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() * 100,
                                                                 predictions_discontinuation_stan_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() * 100))), 5, f = ceiling)

# forest plot
plot_duk_6 <- rbind(
  cbind(mean = 150, lci = 150, uci = 150, drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(mean = 150, lci = 150, uci = 150, drugclass = "GLP1", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_discontinuation_stan_adjusted_overall %>%
    slice(c(1:6)) %>%
    mutate(mean = as.numeric(mean)*100,
           lci = as.numeric(lci)*100,
           uci = as.numeric(uci)*100),
  cbind(mean = 150, lci = 150, uci = 150, drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on GLP-1 RA"),
  cbind(mean = 150, lci = 150, uci = 150, drugclass = "GLP1", intervals = "Predicted HbA1c benefit on GLP-1 RA"),
  predictions_discontinuation_stan_adjusted_overall %>%
    slice(-c(1:6)) %>%
    mutate(mean = as.numeric(mean)*100,
           lci = as.numeric(lci)*100,
           uci = as.numeric(uci)*100),
  predictions_discontinuation_stan_adjusted_full %>% cbind(intervals = "Average treatment effect") %>%
    mutate(mean = as.numeric(mean)*100,
           lci = as.numeric(lci)*100,
           uci = as.numeric(uci)*100)
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == levels(group.discontinuation.dataset$intervals)[1], paste0(">5 mmol/mol (n=", format(group.discontinuation.dataset%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[1])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.discontinuation.dataset%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[1])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                            ifelse(intervals == levels(group.discontinuation.dataset$intervals)[2], paste0("3-5 mmol/mol (n=", format(group.discontinuation.dataset%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[2])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.discontinuation.dataset%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[2])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                   ifelse(intervals == levels(group.discontinuation.dataset$intervals)[3], paste0("0-3 mmol/mol (n=", format(group.discontinuation.dataset%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[3])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.discontinuation.dataset%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[3])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                          ifelse(intervals == levels(group.discontinuation.dataset$intervals)[4], paste0("0-3 mmol/mol (n=", format(group.discontinuation.dataset%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[4])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.discontinuation.dataset%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[4])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                 ifelse(intervals == levels(group.discontinuation.dataset$intervals)[5], paste0("3-5 mmol/mol (n=", format(group.discontinuation.dataset%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[5])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.discontinuation.dataset%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[5])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                        ifelse(intervals == levels(group.discontinuation.dataset$intervals)[6], paste0(">5 mmol/mol (n=", format(group.discontinuation.dataset%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[6])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.discontinuation.dataset%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[6])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"), intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  mutate(group = ifelse(group == "SGLT2", "SGLT2i", "GLP-1 RA")) %>%
  forestplot(ci.vertices = TRUE,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices.height = 0.1,
             # title = "Discontinuation",
             clip = c(0, discontinuation_overall_axis_max),
             xticks = seq(0, discontinuation_overall_axis_max, 5),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=1.5), 
                              title=gpar(cex=2.5),
                              legend =gpar(cex=2),
                              label=gpar(cex=1.5),
                              xlab=gpar(cex=1.5)),
             xlab = "Discontinuation (%)") %>%
  fp_add_header(paste0("Overall population (n=", format(group.discontinuation.dataset%>%nrow(),big.mark=",",scientific=FALSE), ")")) %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE)) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = "black")

# PDF with plots
pdf(width = 11, height = 6, "Plots/DUK_2023/11.08.plot_duk_6.pdf")
plot_duk_6
dev.off()


#:--------------------------------------------------------------------------------
# Forest plot of predicted treatment effect benefit vs microvascular complications
require(forestplot)

interval_breaks <- c(-5, -3, 0, 3, 5)

## Read in propensity scores
patient_prop_scores <- readRDS("Samples/SGLT2-GLP1/Aurum/ps_model/patient_prop_scores.rds")

## Read in treatment effects
treatment_effects <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/patient_effects.rds")

# load new propensity scores
patient_prop_scores_qrisk <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/patient_prop_scores_qrisk.rds")

# Microvascular complications dataset
micro_comp.dataset <- set_up_data_sglt2_glp1(dataset.type="micro_comp.dataset") %>%
  left_join(patient_prop_scores_qrisk, by = c("patid", "pated")) %>%
  left_join(treatment_effects, by = c("patid", "pated"))

# Split population into subgroups
group.micro_comp.dataset <- group_values(data = micro_comp.dataset,
                                         variable = "effects",
                                         breaks = interval_breaks) %>%
  drop_na(intervals)

# predictions for the CVD outcomes in the population with no CVD/HF/CKD
predictions_micro_comp_micro_comp_stan_adjusted_overall <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_micro_comp_micro_comp_stan_adjusted_overall.rds")

# predictions for the CVD outcomes in the population with no CVD/HF/CKD
predictions_micro_comp_micro_comp_stan_adjusted_full <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_micro_comp_micro_comp_stan_adjusted_full.rds")

# forest plot
plot_duk_7 <- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_micro_comp_micro_comp_stan_adjusted_overall %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP-1 RA", mean = NA, lci = NA, uci = NA),
  predictions_micro_comp_micro_comp_stan_adjusted_overall %>% slice(4:6),
  predictions_micro_comp_micro_comp_stan_adjusted_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.micro_comp.dataset$intervals)[1], paste0(">5 mmol/mol (n=", format(group.micro_comp.dataset%>%filter(intervals==levels(group.micro_comp.dataset$intervals)[1])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.micro_comp.dataset%>%filter(intervals==levels(group.micro_comp.dataset$intervals)[1])%>%filter(postdrug_micro_comp_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                            ifelse(intervals == levels(group.micro_comp.dataset$intervals)[2], paste0("3-5 mmol/mol (n=", format(group.micro_comp.dataset%>%filter(intervals==levels(group.micro_comp.dataset$intervals)[2])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.micro_comp.dataset%>%filter(intervals==levels(group.micro_comp.dataset$intervals)[2])%>%filter(postdrug_micro_comp_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                   ifelse(intervals == levels(group.micro_comp.dataset$intervals)[3], paste0("0-3 mmol/mol (n=", format(group.micro_comp.dataset%>%filter(intervals==levels(group.micro_comp.dataset$intervals)[3])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.micro_comp.dataset%>%filter(intervals==levels(group.micro_comp.dataset$intervals)[3])%>%filter(postdrug_micro_comp_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                          ifelse(intervals == levels(group.micro_comp.dataset$intervals)[4], paste0("0-3 mmol/mol (n=", format(group.micro_comp.dataset%>%filter(intervals==levels(group.micro_comp.dataset$intervals)[4])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.micro_comp.dataset%>%filter(intervals==levels(group.micro_comp.dataset$intervals)[4])%>%filter(postdrug_micro_comp_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                 ifelse(intervals == levels(group.micro_comp.dataset$intervals)[5], paste0("3-5 mmol/mol (n=", format(group.micro_comp.dataset%>%filter(intervals==levels(group.micro_comp.dataset$intervals)[5])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.micro_comp.dataset%>%filter(intervals==levels(group.micro_comp.dataset$intervals)[5])%>%filter(postdrug_micro_comp_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                        ifelse(intervals == levels(group.micro_comp.dataset$intervals)[6], paste0(">5 mmol/mol (n=", format(group.micro_comp.dataset%>%filter(intervals==levels(group.micro_comp.dataset$intervals)[6])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.micro_comp.dataset%>%filter(intervals==levels(group.micro_comp.dataset$intervals)[6])%>%filter(postdrug_micro_comp_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci") %>%
  mutate(mean = exp(mean),
         lower = exp(lower),
         upper = exp(upper)) %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             zero = 1,
             # title = "Microvascular complications",
             xlog = TRUE,
             ci.vertices.height = 0.3,
             boxsize = .4,
             lwd.ci = 3,
             txt_gp = fpTxtGp(ticks=gpar(cex=1.5), 
                              title=gpar(cex=2.5),
                              legend =gpar(cex=2),
                              label=gpar(cex=1.5),
                              xlab=gpar(cex=1.5)),
             xlab = "Hazard ratio (95% CI, log scale)") %>%
  fp_add_header(paste0("Overall population (n=", format(nrow(group.micro_comp.dataset),big.mark=",",scientific=FALSE), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = "black")

# PDF with plots
pdf(width = 11, height = 5, "Plots/DUK_2023/11.08.plot_duk_7.pdf")
plot_duk_7
dev.off()



#:--------------------------------------------------------------------------------
# Forest plot of predicted treatment effect benefit vs CVD
require(forestplot)

interval_breaks <- c(-5, -3, 0, 3, 5)

## Read in propensity scores
patient_prop_scores <- readRDS("Samples/SGLT2-GLP1/Aurum/ps_model/patient_prop_scores.rds")

## Read in treatment effects
treatment_effects <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/patient_effects.rds")

# load new propensity scores
patient_prop_scores_qrisk <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/patient_prop_scores_qrisk.rds")

# CVD dataset with no other commorbidities (CVD/CKD/HF)
no_co.dataset <- set_up_data_sglt2_glp1(dataset.type="no_co.dataset") %>%
  left_join(patient_prop_scores_qrisk, by = c("patid", "pated")) %>%
  left_join(treatment_effects, by = c("patid", "pated"))

# Split population into subgroups
group.no_co.dataset <- group_values(data = no_co.dataset,
                                    variable = "effects",
                                    breaks = interval_breaks) %>%
  drop_na(intervals)

# predictions for the CVD outcomes in the population with no CVD/HF/CKD
predictions_no_co_cvd_stan_adjusted_overall <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_no_co_cvd_stan_adjusted_overall.rds")

# predictions for the CVD outcomes in the population with no CVD/HF/CKD
predictions_no_co_cvd_stan_adjusted_full <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_no_co_cvd_stan_adjusted_full.rds")

# forest plot
plot_duk_8 <- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_no_co_cvd_stan_adjusted_overall %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP-1 RA", mean = NA, lci = NA, uci = NA),
  predictions_no_co_cvd_stan_adjusted_overall %>% slice(4:6),
  predictions_no_co_cvd_stan_adjusted_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.no_co.dataset$intervals)[1], paste0(">5 mmol/mol (n=", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[1])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[1])%>%filter(postdrug_mace_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                            ifelse(intervals == levels(group.no_co.dataset$intervals)[2], paste0("3-5 mmol/mol (n=", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[2])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[2])%>%filter(postdrug_mace_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                   ifelse(intervals == levels(group.no_co.dataset$intervals)[3], paste0("0-3 mmol/mol (n=", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[3])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[3])%>%filter(postdrug_mace_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                          ifelse(intervals == levels(group.no_co.dataset$intervals)[4], paste0("0-3 mmol/mol (n=", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[4])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[4])%>%filter(postdrug_mace_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                 ifelse(intervals == levels(group.no_co.dataset$intervals)[5], paste0("3-5 mmol/mol (n=", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[5])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[5])%>%filter(postdrug_mace_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                        ifelse(intervals == levels(group.no_co.dataset$intervals)[6], paste0(">5 mmol/mol (n=", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[6])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[6])%>%filter(postdrug_mace_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci") %>%
  mutate(mean = exp(mean),
         lower = exp(lower),
         upper = exp(upper)) %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             zero = 1,
             title = "CVD outcomes",
             xlog = TRUE,
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazard ratio (95% CI, log scale)") %>%
  fp_add_header(paste0("Overall population (n=", format(nrow(group.no_co.dataset),big.mark=",",scientific=FALSE), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = "black")

# PDF with plots
pdf(width = 7, height = 5, "Plots/DUK_2023/11.08.plot_duk_8.pdf")
plot_duk_8
dev.off()


#:--------------------------------------------------------------------------------
# Anthony outcomes

### Predicted HbA1c outcome

# predictions for the hba1c adjusted model sex strata
predictions_hba1c <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_hba1c_stan_adjusted.rds") %>%
  # only male results
  filter(sex == "Male") %>%
  # only the results for the best predicted HbA1c benefit
  slice(1:2) %>%
  select(mean, drugclass) %>%
  mutate(mean = as.numeric(mean))

plot_1 <- predictions_hba1c %>%
  ggplot() +
  ggtitle(paste0("Predicted HbA1c change (", round(diff(predictions_hba1c$mean), digits = 3), " mmol/mol)")) +
  geom_col(aes(x = mean, y = drugclass, fill = drugclass)) +
  geom_vline(aes(xintercept = 0), colour = "black") +
  geom_text(label = "SGLT2i", y = "SGLT2", x = -1, size = 9, hjust = "right", fontface = "bold") +
  geom_text(label = "GLP-1 RA", y = "GLP1", x = -1, size = 9, hjust = "right", fontface = "bold") +
  scale_x_continuous(breaks = seq(-25, 0, 5)) +
  scale_fill_manual(values = c("dodgerblue2", "#f1a340"))+
  theme(
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.line.x = element_line(colour = "black"),
    legend.position = "none",
    panel.background = element_rect(fill = "grey90"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "grey90")
  )

# predictions for the weight adjusted model
predictions_weight <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_weight_stan_adjusted_overall.rds") %>%
  # only the results for the best predicted HbA1c benefit
  slice(1:2) %>%
  select(mean, drugclass) %>%
  mutate(mean = as.numeric(mean))

plot_2 <- predictions_weight %>%
  ggplot() +
  ggtitle(paste0("Predicted weight change (", round(diff(predictions_weight$mean), digits = 3), " kg)")) +
  geom_col(aes(x = mean, y = drugclass, fill = drugclass)) +
  geom_vline(aes(xintercept = 0), colour = "black") +
  geom_text(label = "SGLT2i", y = "SGLT2", x = -0.1, size = 9, hjust = "right", fontface = "bold") +
  geom_text(label = "GLP-1 RA", y = "GLP1", x = -0.1, size = 9, hjust = "right", fontface = "bold") +
  scale_x_continuous(breaks = seq(-4, 0, 1)) +
  scale_fill_manual(values = c("dodgerblue2", "#f1a340"))+
  theme(
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.line.x = element_line(colour = "black"),
    legend.position = "none",
    panel.background = element_rect(fill = "grey90"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "grey90")
  )

# predictions for the discontinuation adjusted model
predictions_discontinuation <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_discontinuation_stan_adjusted_overall.rds") %>%
  # only the results for the best predicted HbA1c benefit
  slice(1:2) %>%
  select(mean, drugclass) %>%
  mutate(mean = as.numeric(mean)*100)

plot_3 <- predictions_discontinuation %>%
  ggplot() +
  ggtitle(paste0("Discontinuation (", round(diff(predictions_discontinuation$mean), digits = 3), " %)")) +
  geom_col(aes(x = mean, y = drugclass, fill = drugclass)) +
  geom_vline(aes(xintercept = 0), colour = "black") +
  geom_text(label = "SGLT2i", y = "SGLT2", x = 0.3, size = 9, hjust = "left", fontface = "bold") +
  geom_text(label = "GLP-1 RA", y = "GLP1", x = 0.3, size = 9, hjust = "left", fontface = "bold") +
  scale_x_continuous(breaks = seq(0, 20, 5)) +
  scale_fill_manual(values = c("dodgerblue2", "#f1a340"))+
  theme(
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.line.x = element_line(colour = "black"),
    legend.position = "none",
    panel.background = element_rect(fill = "grey90"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "grey90")
  )

# Combine plots
plot_duk_9 <- patchwork::wrap_plots(list(plot_1, plot_2, plot_3), ncol = 1) +
  patchwork::plot_annotation(
    title = "Anthony"
  )

# PDF with plots
pdf(width = 8, height = 6, "Plots/DUK_2023/11.08.plot_duk_9.pdf")
plot_duk_9
dev.off()



#:--------------------------------------------------------------------------------
# Beryl outcomes

### Predicted HbA1c outcome

# predictions for the hba1c adjusted model sex strata
predictions_hba1c <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_hba1c_stan_adjusted.rds") %>%
  # only male results
  filter(sex == "Female") %>%
  # only the results for the best predicted HbA1c benefit
  slice(11:12) %>%
  select(mean, drugclass) %>%
  mutate(mean = as.numeric(mean))

plot_1 <- predictions_hba1c %>%
  ggplot() +
  ggtitle(paste0("Predicted HbA1c change (", round(diff(predictions_hba1c$mean), digits = 3), " mmol/mol)")) +
  geom_col(aes(x = mean, y = drugclass, fill = drugclass)) +
  geom_vline(aes(xintercept = 0), colour = "black") +
  geom_text(label = "SGLT2i", y = "SGLT2", x = -1, size = 9, hjust = "right", fontface = "bold") +
  geom_text(label = "GLP-1 RA", y = "GLP1", x = -1, size = 9, hjust = "right", fontface = "bold") +
  scale_x_continuous(breaks = seq(-20, 0, 5)) +
  scale_fill_manual(values = c("dodgerblue2", "#f1a340"))+
  theme(
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.line.x = element_line(colour = "black"),
    legend.position = "none",
    panel.background = element_rect(fill = "grey90"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "grey90")
  )

# predictions for the weight adjusted model
predictions_weight <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_weight_stan_adjusted_overall.rds") %>%
  # only the results for the best predicted HbA1c benefit
  slice(11:12) %>%
  select(mean, drugclass) %>%
  mutate(mean = as.numeric(mean))

plot_2 <- predictions_weight %>%
  ggplot() +
  ggtitle(paste0("Predicted weight change (", round(diff(predictions_weight$mean), digits = 3), " kg)")) +
  geom_col(aes(x = mean, y = drugclass, fill = drugclass)) +
  geom_vline(aes(xintercept = 0), colour = "black") +
  geom_text(label = "SGLT2i", y = "SGLT2", x = -0.1, size = 9, hjust = "right", fontface = "bold") +
  geom_text(label = "GLP-1 RA", y = "GLP1", x = -0.1, size = 9, hjust = "right", fontface = "bold") +
  scale_x_continuous(breaks = seq(-5, 0, 1)) +
  scale_fill_manual(values = c("dodgerblue2", "#f1a340"))+
  theme(
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.line.x = element_line(colour = "black"),
    legend.position = "none",
    panel.background = element_rect(fill = "grey90"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "grey90")
  )

# predictions for the discontinuation adjusted model
predictions_discontinuation <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_discontinuation_stan_adjusted_overall.rds") %>%
  # only the results for the best predicted HbA1c benefit
  slice(11:12) %>%
  select(mean, drugclass) %>%
  mutate(mean = as.numeric(mean)*100)

plot_3 <- predictions_discontinuation %>%
  ggplot() +
  ggtitle(paste0("Discontinuation (", round(diff(predictions_discontinuation$mean), digits = 3), " %)")) +
  geom_col(aes(x = mean, y = drugclass, fill = drugclass)) +
  geom_vline(aes(xintercept = 0), colour = "black") +
  geom_text(label = "SGLT2i", y = "SGLT2", x = 0.3, size = 9, hjust = "left", fontface = "bold") +
  geom_text(label = "GLP-1 RA", y = "GLP1", x = 0.3, size = 9, hjust = "left", fontface = "bold") +
  scale_x_continuous(breaks = seq(0, 30, 5)) +
  scale_fill_manual(values = c("dodgerblue2", "#f1a340"))+
  theme(
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.line.x = element_line(colour = "black"),
    legend.position = "none",
    panel.background = element_rect(fill = "grey90"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "grey90")
  )

# Combine plots
plot_duk_10 <- patchwork::wrap_plots(list(plot_1, plot_2, plot_3), ncol = 1) +
  patchwork::plot_annotation(
    title = "Beryl"
  )

# PDF with plots
pdf(width = 8, height = 6, "Plots/DUK_2023/11.08.plot_duk_10.pdf")
plot_duk_10
dev.off()


#:--------------------------------------------------------------------------------
# Scottish data
require(forestplot)

# load results
ATE_adjust_validation_subgroups <- readRDS("Samples/SGLT2-GLP1/Scotland_Results/ATE_adjust_validation_subgroups.rds")
ATE_adjust_validation_subgroups_full <- readRDS("Samples/SGLT2-GLP1/Scotland_Results/ATE_adjust_validation_subgroups_full.rds")

# limits for the plot
hba1c_overall_axis_min <- plyr::round_any(ceiling(min(c(ATE_adjust_validation_subgroups %>% select(c("obs","lci","uci")) %>% min(),
                                                        ATE_adjust_validation_subgroups_full %>% select(c("obs","lci","uci")) %>% min()))), 2, f = floor)

hba1c_overall_axis_max <- plyr::round_any(ceiling(max(c(ATE_adjust_validation_subgroups %>% select(c("obs","lci","uci")) %>% max(),
                                                        ATE_adjust_validation_subgroups_full %>% select(c("obs","lci","uci")) %>% max()))), 2, f = ceiling)

# forest plot with results
plot_duk_11 <- rbind(
  cbind(hba1c_diff.q = "Predicted HbA1c benefit on SGLT2i", obs = NA, lci = NA, uci = NA),
  ATE_adjust_validation_subgroups %>% select(hba1c_diff.q, obs, lci, uci) %>% slice(1:3),
  cbind(hba1c_diff.q = "Predicted HbA1c benefit on GLP-1 RA", obs = NA, lci = NA, uci = NA),
  ATE_adjust_validation_subgroups %>% select(hba1c_diff.q, obs, lci, uci) %>% slice(4:6),
  ATE_adjust_validation_subgroups_full %>% select(obs, lci, uci) %>% mutate(hba1c_diff.q = "Average treatment effect")
) %>%
  as.data.frame() %>%
  rename("intervals" = "hba1c_diff.q") %>%
  mutate(obs = as.numeric(obs),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == 1, paste0(">5 mmol/mol (n=", format(ATE_adjust_validation_subgroups%>%select(N)%>%slice(1)%>%unlist(),big.mark=",",scientific=FALSE), ")"), 
                            ifelse(intervals == 2, paste0("3-5 mmol/mol (n=", format(ATE_adjust_validation_subgroups%>%select(N)%>%slice(2)%>%unlist(),big.mark=",",scientific=FALSE), ")"), 
                                   ifelse(intervals == 3, paste0("0-3 mmol/mol (n=", format(ATE_adjust_validation_subgroups%>%select(N)%>%slice(3)%>%unlist(),big.mark=",",scientific=FALSE), ")"), 
                                          ifelse(intervals == 4, paste0("0-3 mmol/mol (n=", format(ATE_adjust_validation_subgroups%>%select(N)%>%slice(4)%>%unlist(),big.mark=",",scientific=FALSE), ")"), 
                                                 ifelse(intervals == 5, paste0("3-5 mmol/mol (n=", format(ATE_adjust_validation_subgroups%>%select(N)%>%slice(5)%>%unlist(),big.mark=",",scientific=FALSE), ")"),
                                                        ifelse(intervals == 6, paste0(">5 mmol/mol (n=", format(ATE_adjust_validation_subgroups%>%select(N)%>%slice(6)%>%unlist(),big.mark=",",scientific=FALSE), ")"), intervals))))))) %>%
  rename("mean" = "obs", "lower" = "lci", "upper" = "uci") %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             # title = "Adjusted model",
             xticks = seq(hba1c_overall_axis_min, hba1c_overall_axis_max, 2),
             ci.vertices.height = 0.3,
             boxsize = .4,
             lwd.ci = 3,
             txt_gp = fpTxtGp(ticks=gpar(cex=1.5), 
                              title=gpar(cex=2.5),
                              legend =gpar(cex=2),
                              label=gpar(cex=1.5),
                              xlab=gpar(cex=1.5)),
             xlab = "Average HbA1c treatment effect (mmol/mol)") %>%
  fp_add_header(paste0("Overall population (n=", format(sum(ATE_adjust_validation_subgroups$N),big.mark=",",scientific=FALSE), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = "black")

# PDF with plot
pdf(width = 12, height = 6, "Plots/DUK_2023/11.08.plot_duk_11.pdf")
plot_duk_11
dev.off()




#:------------------------------------------------------------------------------
# patient characteristics and results treatment effect differences - poster

interval_breaks <- c(-5, -3, 0, 3, 5)



full.cohort <- set_up_data_sglt2_glp1(dataset.type = "hba1c.train") %>%
  # drop the variables with the most missingness (>40%)
  select(-preacr, -preast, -prehaematocrit, -prehaemoglobin, -pretriglyceride) %>%
  # only complete cases
  drop_na() %>%
  as.data.frame() %>%
  select(patid, pated, posthba1cfinal, drugclass, unique(c(readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/variables_mu.rds"), readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/variables_tau.rds")))) %>%
  rbind(
    set_up_data_sglt2_glp1(dataset.type = "hba1c.test") %>%
      # selected variables from SparseBCF
      select(patid, pated, posthba1cfinal, drugclass, unique(c(readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/variables_mu.rds"), readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/variables_tau.rds")))) %>%
      # only complete cases
      drop_na() %>%
      as.data.frame()
  ) %>%
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


plot_sex_strata <- group.full.cohort %>%
  select(sex, intervals) %>%
  ggplot(aes(x = intervals)) +
  geom_bar(aes(fill = sex), position = "fill") +
  ggtitle("Sex") +
  ylab("Percentage (%)") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Accent") +
  theme_bw() +
  annotate("point", x = 1, y = -0.1, shape = 15, colour = "#f1a340", size = 6, alpha = 1) +
  annotate("point", x = 2, y = -0.1, shape = 15, colour = "#f1a340", size = 6, alpha = 0.6) +
  annotate("point", x = 3, y = -0.1, shape = 15, colour = "#f1a340", size = 6, alpha = 0.2) +
  annotate("point", x = 4, y = -0.1, shape = 15, colour = "dodgerblue2", size = 6, alpha = 0.2) +
  annotate("point", x = 5, y = -0.1, shape = 15, colour = "dodgerblue2", size = 6, alpha = 0.6) +
  annotate("point", x = 6, y = -0.1, shape = 15, colour = "dodgerblue2", size = 6, alpha = 1) +
  coord_cartesian(ylim = c(0, 1), clip = "off") +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 30),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 30),
        legend.text = element_text(size= 30),
        legend.title = element_text(size = 30), 
        axis.title.x = element_text(vjust = -50, colour = "white"),
        strip.background = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5, size = 30)) +
  guides(fill=guide_legend(title="Sex", title.position="top", title.hjust = 0.5))

plot_prehba1c_strata <- group.full.cohort %>%
  select(intervals, prehba1c) %>%
  drop_na() %>%
  ggplot(aes(x = intervals, y = prehba1c)) +
  geom_point(aes(fill = intervals, colour = intervals), alpha = 0, shape = 15, size = 5) +
  geom_boxplot(outlier.shape = NA, aes(fill = intervals), alpha = c(1, 0.6, 0.2, 0.2, 0.6, 1), show.legend = F) +
  ggtitle("Baseline HbA1c") +
  ylab("Baseline HbA1c (mmol/mol)") +
  theme_bw() +
  annotate("point", x = 1, y = 40, shape = 15, colour = "#f1a340", size = 6, alpha = 1) +
  annotate("point", x = 2, y = 40, shape = 15, colour = "#f1a340", size = 6, alpha = 0.6) +
  annotate("point", x = 3, y = 40, shape = 15, colour = "#f1a340", size = 6, alpha = 0.2) +
  annotate("point", x = 4, y = 40, shape = 15, colour = "dodgerblue2", size = 6, alpha = 0.2) +
  annotate("point", x = 5, y = 40, shape = 15, colour = "dodgerblue2", size = 6, alpha = 0.6) +
  annotate("point", x = 6, y = 40, shape = 15, colour = "dodgerblue2", size = 6, alpha = 1) +
  coord_cartesian(ylim = c(50, 150), clip = "off") +
  scale_y_continuous(breaks = c(53, 64, 75, 86, 97, 108, 119, 150)) +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 30),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 30),
        legend.text = element_text(size= 30),
        legend.title = element_text(size = 30), 
        axis.title.x = element_text(vjust = -50, colour = "white"),
        strip.background = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5, size = 30),
        panel.grid.minor = element_blank()) +
  scale_fill_manual(
    name = paste0("HbA1c benefit (mmol/mol) [n=", format(nrow(group.full.cohort),big.mark=",",scientific=FALSE),"]"),
    labels = c(paste0("SGLT2i > 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[1]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("SGLT2i 3 - 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[2]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("SGLT2i 0 - 3 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[3]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("GLP1-RA > 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[4]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("GLP1-RA 3 - 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[5]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("GLP1-RA 0 - 3 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[6]) %>% nrow(),big.mark=",",scientific=FALSE), "]")),
    breaks = c(levels(group.full.cohort$intervals)[1], levels(group.full.cohort$intervals)[2], levels(group.full.cohort$intervals)[3], levels(group.full.cohort$intervals)[6], levels(group.full.cohort$intervals)[5], levels(group.full.cohort$intervals)[4]),
    values = c("SGLT2i benefit >5 mmol/mol" = "#f1a340",
               "SGLT2i benefit 3-5 mmol/mol" = "#f1a340",
               "SGLT2i benefit 0-3 mmol/mol" = "#f1a340",
               "GLP1-RA benefit 0-3 mmol/mol" = "dodgerblue2",
               "GLP1-RA benefit 3-5 mmol/mol" = "dodgerblue2",
               "GLP1-RA benefit >5 mmol/mol" = "dodgerblue2")
  ) +
  scale_colour_manual(
    name = paste0("HbA1c benefit (mmol/mol) [n=", format(nrow(group.full.cohort),big.mark=",",scientific=FALSE),"]"),
    labels = c(paste0("SGLT2i > 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[1]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("SGLT2i 3 - 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[2]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("SGLT2i 0 - 3 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[3]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("GLP1-RA > 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[4]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("GLP1-RA 3 - 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[5]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("GLP1-RA 0 - 3 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[6]) %>% nrow(),big.mark=",",scientific=FALSE), "]")),
    breaks = c(levels(group.full.cohort$intervals)[1], levels(group.full.cohort$intervals)[2], levels(group.full.cohort$intervals)[3], levels(group.full.cohort$intervals)[6], levels(group.full.cohort$intervals)[5], levels(group.full.cohort$intervals)[4]),
    values = c("SGLT2i benefit >5 mmol/mol" = "#f1a340",
               "SGLT2i benefit 3-5 mmol/mol" = "#f1a340",
               "SGLT2i benefit 0-3 mmol/mol" = "#f1a340",
               "GLP1-RA benefit 0-3 mmol/mol" = "dodgerblue2",
               "GLP1-RA benefit 3-5 mmol/mol" = "dodgerblue2",
               "GLP1-RA benefit >5 mmol/mol" = "dodgerblue2")
  ) +
  guides(
    fill = guide_legend(title.position="top", title.hjust = 0.5, override.aes = list(alpha = c(1, 0.6, 0.2, 1, 0.6, 0.2), size = 6), shape = c(15,15,15,15,15,15),  ncol = 2, nrow = 3),
    colour = guide_legend(title.position="top", title.hjust = 0.5, override.aes = list(alpha = c(1, 0.6, 0.2, 1, 0.6, 0.2), size = 6), shape = c(15,15,15,15,15,15),  ncol = 2, nrow = 3)
  )


plot_preegfr_strata <- group.full.cohort %>%
  select(intervals, preegfr) %>%
  drop_na() %>%
  ggplot(aes(x = intervals, y = preegfr)) +
  geom_point(aes(fill = intervals, colour = intervals), alpha = 0, shape = 15, size = 5) +
  geom_boxplot(outlier.shape = NA, aes(fill = intervals), alpha = c(1, 0.6, 0.2, 0.2, 0.6, 1), show.legend = F) +
  ggtitle("eGFR") +
  ylab(expression(paste("eGFR (mL/min/1.3", m^{2}, ")"))) +
  scale_y_continuous(breaks = c(30, 60, 90, 120, 150)) +
  theme_bw() +
  annotate("point", x = 1, y = 18, shape = 15, colour = "#f1a340", size = 6, alpha = 1) +
  annotate("point", x = 2, y = 18, shape = 15, colour = "#f1a340", size = 6, alpha = 0.6) +
  annotate("point", x = 3, y = 18, shape = 15, colour = "#f1a340", size = 6, alpha = 0.2) +
  annotate("point", x = 4, y = 18, shape = 15, colour = "dodgerblue2", size = 6, alpha = 0.2) +
  annotate("point", x = 5, y = 18, shape = 15, colour = "dodgerblue2", size = 6, alpha = 0.6) +
  annotate("point", x = 6, y = 18, shape = 15, colour = "dodgerblue2", size = 6, alpha = 1) +
  coord_cartesian(ylim = c(30, 150), clip = "off") +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 30),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 30),
        legend.text = element_text(size= 30),
        legend.title = element_text(size = 30), 
        axis.title.x = element_text(vjust = -50, colour = "white"),
        strip.background = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5, size = 30)) +
  scale_fill_manual(
    name = paste0("HbA1c benefit (mmol/mol) [n=", format(nrow(group.full.cohort),big.mark=",",scientific=FALSE),"]"),
    labels = c(paste0("SGLT2i > 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[1]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("SGLT2i 3 - 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[2]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("SGLT2i 0 - 3 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[3]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("GLP1-RA > 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[4]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("GLP1-RA 3 - 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[5]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("GLP1-RA 0 - 3 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[6]) %>% nrow(),big.mark=",",scientific=FALSE), "]")),
    breaks = c(levels(group.full.cohort$intervals)[1], levels(group.full.cohort$intervals)[2], levels(group.full.cohort$intervals)[3], levels(group.full.cohort$intervals)[6], levels(group.full.cohort$intervals)[5], levels(group.full.cohort$intervals)[4]),
    values = c("SGLT2i benefit >5 mmol/mol" = "#f1a340",
               "SGLT2i benefit 3-5 mmol/mol" = "#f1a340",
               "SGLT2i benefit 0-3 mmol/mol" = "#f1a340",
               "GLP1-RA benefit 0-3 mmol/mol" = "dodgerblue2",
               "GLP1-RA benefit 3-5 mmol/mol" = "dodgerblue2",
               "GLP1-RA benefit >5 mmol/mol" = "dodgerblue2")
  ) +
  scale_colour_manual(
    name = paste0("HbA1c benefit (mmol/mol) [n=", format(nrow(group.full.cohort),big.mark=",",scientific=FALSE),"]"),
    labels = c(paste0("SGLT2i > 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[1]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("SGLT2i 3 - 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[2]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("SGLT2i 0 - 3 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[3]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("GLP1-RA > 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[4]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("GLP1-RA 3 - 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[5]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("GLP1-RA 0 - 3 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[6]) %>% nrow(),big.mark=",",scientific=FALSE), "]")),
    breaks = c(levels(group.full.cohort$intervals)[1], levels(group.full.cohort$intervals)[2], levels(group.full.cohort$intervals)[3], levels(group.full.cohort$intervals)[6], levels(group.full.cohort$intervals)[5], levels(group.full.cohort$intervals)[4]),
    values = c("SGLT2i benefit >5 mmol/mol" = "#f1a340",
               "SGLT2i benefit 3-5 mmol/mol" = "#f1a340",
               "SGLT2i benefit 0-3 mmol/mol" = "#f1a340",
               "GLP1-RA benefit 0-3 mmol/mol" = "dodgerblue2",
               "GLP1-RA benefit 3-5 mmol/mol" = "dodgerblue2",
               "GLP1-RA benefit >5 mmol/mol" = "dodgerblue2")
  ) +
  guides(
    fill = guide_legend(title.position="top", title.hjust = 0.5, override.aes = list(alpha = c(1, 0.6, 0.2, 1, 0.6, 0.2), size = 6), shape = c(15,15,15,15,15,15),  ncol = 2, nrow = 3),
    colour = guide_legend(title.position="top", title.hjust = 0.5, override.aes = list(alpha = c(1, 0.6, 0.2, 1, 0.6, 0.2), size = 6), shape = c(15,15,15,15,15,15),  ncol = 2, nrow = 3)
  )



plot_agetx_strata <- group.full.cohort %>%
  select(intervals, agetx) %>%
  drop_na() %>%
  ggplot(aes(x = intervals, y = agetx)) +
  geom_point(aes(fill = intervals, colour = intervals), alpha = 0, shape = 15, size = 5) +
  geom_boxplot(outlier.shape = NA, aes(fill = intervals), alpha = c(1, 0.6, 0.2, 0.2, 0.6, 1), show.legend = F) +
  ggtitle("Current age") +
  ylab("Current age (years)") +
  theme_bw() +
  annotate("point", x = 1, y = 10, shape = 15, colour = "#f1a340", size = 6, alpha = 1) +
  annotate("point", x = 2, y = 10, shape = 15, colour = "#f1a340", size = 6, alpha = 0.6) +
  annotate("point", x = 3, y = 10, shape = 15, colour = "#f1a340", size = 6, alpha = 0.2) +
  annotate("point", x = 4, y = 10, shape = 15, colour = "dodgerblue2", size = 6, alpha = 0.2) +
  annotate("point", x = 5, y = 10, shape = 15, colour = "dodgerblue2", size = 6, alpha = 0.6) +
  annotate("point", x = 6, y = 10, shape = 15, colour = "dodgerblue2", size = 6, alpha = 1) +
  coord_cartesian(ylim = c(18, 100), clip = "off") +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 30),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 30),
        legend.text = element_text(size= 30),
        legend.title = element_text(size = 30), 
        axis.title.x = element_text(vjust = -50, colour = "white"),
        strip.background = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5, size = 30)) +
  scale_fill_manual(
    name = paste0("HbA1c benefit (mmol/mol) [n=", format(nrow(group.full.cohort),big.mark=",",scientific=FALSE),"]"),
    labels = c(paste0("SGLT2i > 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[1]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("SGLT2i 3 - 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[2]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("SGLT2i 0 - 3 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[3]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("GLP1-RA > 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[4]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("GLP1-RA 3 - 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[5]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("GLP1-RA 0 - 3 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[6]) %>% nrow(),big.mark=",",scientific=FALSE), "]")),
    breaks = c(levels(group.full.cohort$intervals)[1], levels(group.full.cohort$intervals)[2], levels(group.full.cohort$intervals)[3], levels(group.full.cohort$intervals)[6], levels(group.full.cohort$intervals)[5], levels(group.full.cohort$intervals)[4]),
    values = c("SGLT2i benefit >5 mmol/mol" = "#f1a340",
               "SGLT2i benefit 3-5 mmol/mol" = "#f1a340",
               "SGLT2i benefit 0-3 mmol/mol" = "#f1a340",
               "GLP1-RA benefit 0-3 mmol/mol" = "dodgerblue2",
               "GLP1-RA benefit 3-5 mmol/mol" = "dodgerblue2",
               "GLP1-RA benefit >5 mmol/mol" = "dodgerblue2")
  ) +
  scale_colour_manual(
    name = paste0("HbA1c benefit (mmol/mol) [n=", format(nrow(group.full.cohort),big.mark=",",scientific=FALSE),"]"),
    labels = c(paste0("SGLT2i > 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[1]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("SGLT2i 3 - 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[2]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("SGLT2i 0 - 3 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[3]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("GLP1-RA > 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[4]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("GLP1-RA 3 - 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[5]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("GLP1-RA 0 - 3 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[6]) %>% nrow(),big.mark=",",scientific=FALSE), "]")),
    breaks = c(levels(group.full.cohort$intervals)[1], levels(group.full.cohort$intervals)[2], levels(group.full.cohort$intervals)[3], levels(group.full.cohort$intervals)[6], levels(group.full.cohort$intervals)[5], levels(group.full.cohort$intervals)[4]),
    values = c("SGLT2i benefit >5 mmol/mol" = "#f1a340",
               "SGLT2i benefit 3-5 mmol/mol" = "#f1a340",
               "SGLT2i benefit 0-3 mmol/mol" = "#f1a340",
               "GLP1-RA benefit 0-3 mmol/mol" = "dodgerblue2",
               "GLP1-RA benefit 3-5 mmol/mol" = "dodgerblue2",
               "GLP1-RA benefit >5 mmol/mol" = "dodgerblue2")
  ) +
  guides(
    fill = guide_legend(title.position="top", title.hjust = 0.5, override.aes = list(alpha = c(1, 0.6, 0.2, 1, 0.6, 0.2), size = 6), shape = c(15,15,15,15,15,15),  ncol = 2, nrow = 3),
    colour = guide_legend(title.position="top", title.hjust = 0.5, override.aes = list(alpha = c(1, 0.6, 0.2, 1, 0.6, 0.2), size = 6), shape = c(15,15,15,15,15,15),  ncol = 2, nrow = 3)
  )


plot_prebmi_strata <- group.full.cohort %>%
  select(intervals, prebmi) %>%
  drop_na() %>%
  ggplot(aes(x = intervals, y = prebmi)) +
  geom_point(aes(fill = intervals, colour = intervals), alpha = 0, shape = 15, size = 5) +
  geom_boxplot(outlier.shape = NA, aes(fill = intervals), alpha = c(1, 0.6, 0.2, 0.2, 0.6, 1), show.legend = F) +
  ggtitle("BMI") +
  ylab(expression(paste("BMI (kg/", m^{2}, ")"))) +
  theme_bw() +
  annotate("point", x = 1, y = 10.9, shape = 15, colour = "#f1a340", size = 6, alpha = 1) +
  annotate("point", x = 2, y = 10.9, shape = 15, colour = "#f1a340", size = 6, alpha = 0.6) +
  annotate("point", x = 3, y = 10.9, shape = 15, colour = "#f1a340", size = 6, alpha = 0.2) +
  annotate("point", x = 4, y = 10.9, shape = 15, colour = "dodgerblue2", size = 6, alpha = 0.2) +
  annotate("point", x = 5, y = 10.9, shape = 15, colour = "dodgerblue2", size = 6, alpha = 0.6) +
  annotate("point", x = 6, y = 10.9, shape = 15, colour = "dodgerblue2", size = 6, alpha = 1) +
  coord_cartesian(ylim = c(15, 55), clip = "off") +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 30),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 30),
        legend.text = element_text(size= 30),
        legend.title = element_text(size = 30), 
        axis.title.x = element_text(vjust = -50, colour = "white"),
        strip.background = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5, size = 30)) +
  scale_fill_manual(
    name = paste0("HbA1c benefit (mmol/mol) [n=", format(nrow(group.full.cohort),big.mark=",",scientific=FALSE),"]"),
    labels = c(paste0("SGLT2i > 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[1]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("SGLT2i 3 - 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[2]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("SGLT2i 0 - 3 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[3]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("GLP1-RA > 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[4]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("GLP1-RA 3 - 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[5]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("GLP1-RA 0 - 3 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[6]) %>% nrow(),big.mark=",",scientific=FALSE), "]")),
    breaks = c(levels(group.full.cohort$intervals)[1], levels(group.full.cohort$intervals)[2], levels(group.full.cohort$intervals)[3], levels(group.full.cohort$intervals)[6], levels(group.full.cohort$intervals)[5], levels(group.full.cohort$intervals)[4]),
    values = c("SGLT2i benefit >5 mmol/mol" = "#f1a340",
               "SGLT2i benefit 3-5 mmol/mol" = "#f1a340",
               "SGLT2i benefit 0-3 mmol/mol" = "#f1a340",
               "GLP1-RA benefit 0-3 mmol/mol" = "dodgerblue2",
               "GLP1-RA benefit 3-5 mmol/mol" = "dodgerblue2",
               "GLP1-RA benefit >5 mmol/mol" = "dodgerblue2")
  ) +
  scale_colour_manual(
    name = paste0("HbA1c benefit (mmol/mol) [n=", format(nrow(group.full.cohort),big.mark=",",scientific=FALSE),"]"),
    labels = c(paste0("SGLT2i > 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[1]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("SGLT2i 3 - 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[2]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("SGLT2i 0 - 3 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[3]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("GLP1-RA > 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[4]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("GLP1-RA 3 - 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[5]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("GLP1-RA 0 - 3 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[6]) %>% nrow(),big.mark=",",scientific=FALSE), "]")),
    breaks = c(levels(group.full.cohort$intervals)[1], levels(group.full.cohort$intervals)[2], levels(group.full.cohort$intervals)[3], levels(group.full.cohort$intervals)[6], levels(group.full.cohort$intervals)[5], levels(group.full.cohort$intervals)[4]),
    values = c("SGLT2i benefit >5 mmol/mol" = "#f1a340",
               "SGLT2i benefit 3-5 mmol/mol" = "#f1a340",
               "SGLT2i benefit 0-3 mmol/mol" = "#f1a340",
               "GLP1-RA benefit 0-3 mmol/mol" = "dodgerblue2",
               "GLP1-RA benefit 3-5 mmol/mol" = "dodgerblue2",
               "GLP1-RA benefit >5 mmol/mol" = "dodgerblue2")
  ) +
  guides(
    fill = guide_legend(title.position="top", title.hjust = 0.5, override.aes = list(alpha = c(1, 0.6, 0.2, 1, 0.6, 0.2), size = 6), shape = c(15,15,15,15,15,15),  ncol = 2, nrow = 3),
    colour = guide_legend(title.position="top", title.hjust = 0.5, override.aes = list(alpha = c(1, 0.6, 0.2, 1, 0.6, 0.2), size = 6), shape = c(15,15,15,15,15,15),  ncol = 2, nrow = 3)
  )



plot_duk_12 <- patchwork::wrap_plots(list(
  plot_sex_strata,
  plot_prehba1c_strata,
  plot_preegfr_strata,
  plot_agetx_strata
), nrow = 2, ncol = 2) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18),
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 20)
  )

pdf(width = 34/2.54, height = 26.88/2.54, "Plots/DUK_2023/11.08.plot_duk_12.pdf")
plot_duk_12
dev.off()




#:--------------------------------------------------------------------------------
# Histogram of predicted treatment effects for overall population - poster

# load treatment effects of overall population
patient_effects <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/patient_effects.rds")

# plot histogram of treatment effects
plot_duk_13 <- hist_plot(data = patient_effects %>%
                           rename("mean" = "effects"),
                         title = paste0("Overall cohort (n=", format(nrow(patient_effects),big.mark=",",scientific=FALSE), ")"),
                         xmin = -15,
                         xmax = 20) +
  scale_y_continuous(label=comma) +
  scale_x_continuous(breaks = c(-15, -10, -5, 0, 5, 10, 15), limits = c(-15, 15)) +
  # geom_vline(aes(xintercept = -5), lty = "dashed", colour = "red") +
  # geom_vline(aes(xintercept = 5), lty = "dashed", colour = "red") +
  ggtitle("Development cohort (n=27,319)") +
  theme(
    plot.title = element_text(size =30),
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 30),
    axis.title.x = element_text(size = 32)
  )

# PDF with plot
pdf(width = 20.72/2.54, height = 18.36/2.54, "Plots/DUK_2023/11.08.plot_duk_13.pdf")
plot_duk_13
dev.off()



# :-------------------------------------------------------------------------------
# Calibration plots for development and validation cohorts - poster

# load validation cohort model calibration
ATE_adjust_validation_val <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/assessment/ATE_adjust_validation_val.rds")

# plot validation cohort calibration
plot_ATE_adjust_validation_val <- ATE_plot(ATE_adjust_validation_val[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -8, 8) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle(paste0("Validation cohort (n=",format(sum(ATE_adjust_validation_val[["effects"]]$N),big.mark=",",scientific=FALSE), ")"))

# Combine plots together
plot_duk_14 <- plot_ATE_adjust_validation_val +
  theme(
    plot.title = element_text(size = 30),
    axis.text = element_text(size = 30),
    axis.title = element_text(size = 16)
  )

# PDF with plot
pdf(width = 17.2/2.54, height = 18.36/2.54, "Plots/DUK_2023/11.08.plot_duk_14.pdf")
plot_duk_14
dev.off()



#:--------------------------------------------------------------------------------
# Forest plot of predicted treatment effect benefit vs weight change - poster
require(forestplot)

interval_breaks <- c(-5, -3, 0, 3, 5)

## Read in propensity scores
patient_prop_scores <- readRDS("Samples/SGLT2-GLP1/Aurum/ps_model/patient_prop_scores.rds")

## Read in treatment effects
treatment_effects <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/patient_effects.rds")

# Weight change population
weight.dataset <- set_up_data_sglt2_glp1(dataset.type = "weight.dataset") %>%
  left_join(patient_prop_scores, by = c("patid", "pated")) %>%
  left_join(treatment_effects, by = c("patid", "pated")) %>%
  mutate(w.change = postweight - preweight)

# Split population into subgroups
group.weight.dataset <- group_values(data = weight.dataset,
                                     variable = "effects",
                                     breaks = interval_breaks) %>%
  drop_na(intervals)

# predictions for the weight adjusted model
predictions_weight_stan_adjusted_overall <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_weight_stan_adjusted_overall.rds")

# predictions for the weight adjusted model full
predictions_weight_stan_adjusted_full <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_weight_stan_adjusted_full.rds")

# limits for the plot
weight_strata_axis_min <- plyr::round_any(floor(min(c(predictions_weight_stan_adjusted_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min(),
                                                      predictions_weight_stan_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min()))), 1, f = floor)

weight_strata_axis_max <- plyr::round_any(ceiling(max(c(predictions_weight_stan_adjusted_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max(),
                                                        predictions_weight_stan_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max()))), 2, f = ceiling)

# forest plot
plot_duk_15 <- rbind(
  cbind(mean = 50, lci = 50, uci = 50, drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(mean = 50, lci = 50, uci = 50, drugclass = "GLP1", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_weight_stan_adjusted_overall %>%
    slice(c(1:6)),
  cbind(mean = 50, lci = 50, uci = 50, drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on GLP-1 RA"),
  cbind(mean = 50, lci = 50, uci = 50, drugclass = "GLP1", intervals = "Predicted HbA1c benefit on GLP-1 RA"),
  predictions_weight_stan_adjusted_overall %>%
    slice(-c(1:6)),
  predictions_weight_stan_adjusted_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == levels(group.weight.dataset$intervals)[1], paste0(">5 (n=", format(group.weight.dataset%>%filter(intervals==levels(group.weight.dataset$intervals)[1])%>%nrow(),big.mark=",",scientific=FALSE),")"),
                            ifelse(intervals == levels(group.weight.dataset$intervals)[2], paste0("3-5 (n=", format(group.weight.dataset%>%filter(intervals==levels(group.weight.dataset$intervals)[2])%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                   ifelse(intervals == levels(group.weight.dataset$intervals)[3], paste0("0-3 (n=", format(group.weight.dataset%>%filter(intervals==levels(group.weight.dataset$intervals)[3])%>%nrow(),big.mark=",",scientific=FALSE),")"),
                                          ifelse(intervals == levels(group.weight.dataset$intervals)[4], paste0("0-3 (n=", format(group.weight.dataset%>%filter(intervals==levels(group.weight.dataset$intervals)[4])%>%nrow(),big.mark=",",scientific=FALSE),")"),
                                                 ifelse(intervals == levels(group.weight.dataset$intervals)[5], paste0("3-5 (n=", format(group.weight.dataset%>%filter(intervals==levels(group.weight.dataset$intervals)[5])%>%nrow(),big.mark=",",scientific=FALSE),")"),
                                                        ifelse(intervals == levels(group.weight.dataset$intervals)[6], paste0(">5 (n=", format(group.weight.dataset%>%filter(intervals==levels(group.weight.dataset$intervals)[6])%>%nrow(),big.mark=",",scientific=FALSE),")"), intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  mutate(group = ifelse(group == "SGLT2", "SGLT2i", "GLP-1 RA")) %>%
  forestplot(ci.vertices = TRUE,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices.height = 0.3,
             title = "Weight change",
             clip = c(weight_strata_axis_min, weight_strata_axis_max),
             xticks = seq(weight_strata_axis_min, weight_strata_axis_max, 1),
             boxsize = .4,
             lwd.ci = 3,
             txt_gp = fpTxtGp(ticks=gpar(cex=1.5), 
                              title=gpar(cex=2.5),
                              legend =gpar(cex=2),
                              label=gpar(cex=1.5),
                              xlab=gpar(cex=1.5)),
             xlab = "Predicted weight change (kg) (less is better)") %>%
  fp_add_header(paste0("Overall population (n=", format(group.weight.dataset%>%nrow(),big.mark=",",scientific=FALSE), ")")) %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE)) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = "black")

# PDF with plots
pdf(width = 27.07/2.54, height = 18.03/2.54, "Plots/DUK_2023/11.08.plot_duk_15.pdf")
plot_duk_15
dev.off()


#:--------------------------------------------------------------------------------
# Forest plot of predicted treatment effect benefit vs discontinuation - poster
require(forestplot)

interval_breaks <- c(-5, -3, 0, 3, 5)

## Read in propensity scores
patient_prop_scores <- readRDS("Samples/SGLT2-GLP1/Aurum/ps_model/patient_prop_scores.rds")

## Read in treatment effects
treatment_effects <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/patient_effects.rds")

# Discontinuation
discontinuation.dataset <- set_up_data_sglt2_glp1(dataset.type = "discontinuation.dataset") %>%
  left_join(patient_prop_scores, by = c("patid", "pated")) %>%
  left_join(treatment_effects, by = c("patid", "pated")) %>%
  mutate(stopdrug_6m_3mFU = factor(stopdrug_6m_3mFU))

# Split population into subgroups
group.discontinuation.dataset <- group_values(data = discontinuation.dataset,
                                              variable = "effects",
                                              breaks = interval_breaks) %>%
  drop_na(intervals)

# predictions for the discontinuation adjusted model
predictions_discontinuation_stan_adjusted_overall <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_discontinuation_stan_adjusted_overall.rds")

# predictions for the discontinuation model sex strata
predictions_discontinuation_stan_adjusted_full <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_discontinuation_stan_adjusted_full.rds")

# limits for the plot
discontinuation_overall_axis_max <- plyr::round_any(ceiling(max(c(predictions_discontinuation_stan_adjusted_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() * 100,
                                                                  predictions_discontinuation_stan_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() * 100))), 5, f = ceiling)

discontinuation_strata_axis_max <- plyr::round_any(ceiling(max(c(predictions_discontinuation_stan_adjusted_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() * 100,
                                                                 predictions_discontinuation_stan_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() * 100))), 5, f = ceiling)

# forest plot
plot_duk_16 <- rbind(
  cbind(mean = 150, lci = 150, uci = 150, drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(mean = 150, lci = 150, uci = 150, drugclass = "GLP1", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_discontinuation_stan_adjusted_overall %>%
    slice(c(1:6)) %>%
    mutate(mean = as.numeric(mean)*100,
           lci = as.numeric(lci)*100,
           uci = as.numeric(uci)*100),
  cbind(mean = 150, lci = 150, uci = 150, drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on GLP-1 RA"),
  cbind(mean = 150, lci = 150, uci = 150, drugclass = "GLP1", intervals = "Predicted HbA1c benefit on GLP-1 RA"),
  predictions_discontinuation_stan_adjusted_overall %>%
    slice(-c(1:6)) %>%
    mutate(mean = as.numeric(mean)*100,
           lci = as.numeric(lci)*100,
           uci = as.numeric(uci)*100),
  predictions_discontinuation_stan_adjusted_full %>% cbind(intervals = "Average treatment effect") %>%
    mutate(mean = as.numeric(mean)*100,
           lci = as.numeric(lci)*100,
           uci = as.numeric(uci)*100)
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == levels(group.discontinuation.dataset$intervals)[1], paste0(">5 (n=", format(group.discontinuation.dataset%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[1])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.discontinuation.dataset%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[1])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                            ifelse(intervals == levels(group.discontinuation.dataset$intervals)[2], paste0("3-5 (n=", format(group.discontinuation.dataset%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[2])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.discontinuation.dataset%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[2])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                   ifelse(intervals == levels(group.discontinuation.dataset$intervals)[3], paste0("0-3 (n=", format(group.discontinuation.dataset%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[3])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.discontinuation.dataset%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[3])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                          ifelse(intervals == levels(group.discontinuation.dataset$intervals)[4], paste0("0-3 (n=", format(group.discontinuation.dataset%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[4])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.discontinuation.dataset%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[4])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                 ifelse(intervals == levels(group.discontinuation.dataset$intervals)[5], paste0("3-5 (n=", format(group.discontinuation.dataset%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[5])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.discontinuation.dataset%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[5])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                        ifelse(intervals == levels(group.discontinuation.dataset$intervals)[6], paste0(">5 (n=", format(group.discontinuation.dataset%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[6])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.discontinuation.dataset%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[6])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"), intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  mutate(group = ifelse(group == "SGLT2", "SGLT2i", "GLP-1 RA")) %>%
  forestplot(ci.vertices = TRUE,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices.height = 0.1,
             title = "Discontinuation",
             clip = c(0, discontinuation_overall_axis_max),
             xticks = seq(0, discontinuation_overall_axis_max, 5),
             boxsize = .2,
             # lwd.ci = 3,
             txt_gp = fpTxtGp(ticks=gpar(cex=1.5), 
                              title=gpar(cex=2.5),
                              legend =gpar(cex=2),
                              label=gpar(cex=1.5),
                              xlab=gpar(cex=1.5)),
             xlab = "Discontinuation (%) (less is better)") %>%
  fp_add_header(paste0("Overall population (n=", format(group.discontinuation.dataset%>%nrow(),big.mark=",",scientific=FALSE), ")")) %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE)) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = "black")

# PDF with plots
pdf(width = 27.07/2.54, height = 18.03/2.54, "Plots/DUK_2023/11.08.plot_duk_16.pdf")
plot_duk_16
dev.off()




#:--------------------------------------------------------------------------------
# Forest plot of predicted treatment effect benefit vs microvascular complications - poster
require(forestplot)

interval_breaks <- c(-5, -3, 0, 3, 5)

## Read in propensity scores
patient_prop_scores <- readRDS("Samples/SGLT2-GLP1/Aurum/ps_model/patient_prop_scores.rds")

## Read in treatment effects
treatment_effects <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/patient_effects.rds")

# load new propensity scores
patient_prop_scores_qrisk <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/patient_prop_scores_qrisk.rds")

# Microvascular complications dataset
micro_comp.dataset <- set_up_data_sglt2_glp1(dataset.type="micro_comp.dataset") %>%
  left_join(patient_prop_scores_qrisk, by = c("patid", "pated")) %>%
  left_join(treatment_effects, by = c("patid", "pated"))

# Split population into subgroups
group.micro_comp.dataset <- group_values(data = micro_comp.dataset,
                                         variable = "effects",
                                         breaks = interval_breaks) %>%
  drop_na(intervals)

# predictions for the CVD outcomes in the population with no CVD/HF/CKD
predictions_micro_comp_micro_comp_stan_adjusted_overall <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_micro_comp_micro_comp_stan_adjusted_overall.rds")

# predictions for the CVD outcomes in the population with no CVD/HF/CKD
predictions_micro_comp_micro_comp_stan_adjusted_full <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_micro_comp_micro_comp_stan_adjusted_full.rds")

# forest plot
plot_duk_17 <- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_micro_comp_micro_comp_stan_adjusted_overall %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP-1 RA", mean = NA, lci = NA, uci = NA),
  predictions_micro_comp_micro_comp_stan_adjusted_overall %>% slice(4:6),
  predictions_micro_comp_micro_comp_stan_adjusted_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.micro_comp.dataset$intervals)[1], paste0(">5 (n=", format(group.micro_comp.dataset%>%filter(intervals==levels(group.micro_comp.dataset$intervals)[1])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.micro_comp.dataset%>%filter(intervals==levels(group.micro_comp.dataset$intervals)[1])%>%filter(postdrug_micro_comp_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                            ifelse(intervals == levels(group.micro_comp.dataset$intervals)[2], paste0("3-5 (n=", format(group.micro_comp.dataset%>%filter(intervals==levels(group.micro_comp.dataset$intervals)[2])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.micro_comp.dataset%>%filter(intervals==levels(group.micro_comp.dataset$intervals)[2])%>%filter(postdrug_micro_comp_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                   ifelse(intervals == levels(group.micro_comp.dataset$intervals)[3], paste0("0-3 (n=", format(group.micro_comp.dataset%>%filter(intervals==levels(group.micro_comp.dataset$intervals)[3])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.micro_comp.dataset%>%filter(intervals==levels(group.micro_comp.dataset$intervals)[3])%>%filter(postdrug_micro_comp_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                          ifelse(intervals == levels(group.micro_comp.dataset$intervals)[4], paste0("0-3 (n=", format(group.micro_comp.dataset%>%filter(intervals==levels(group.micro_comp.dataset$intervals)[4])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.micro_comp.dataset%>%filter(intervals==levels(group.micro_comp.dataset$intervals)[4])%>%filter(postdrug_micro_comp_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                 ifelse(intervals == levels(group.micro_comp.dataset$intervals)[5], paste0("3-5 (n=", format(group.micro_comp.dataset%>%filter(intervals==levels(group.micro_comp.dataset$intervals)[5])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.micro_comp.dataset%>%filter(intervals==levels(group.micro_comp.dataset$intervals)[5])%>%filter(postdrug_micro_comp_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                        ifelse(intervals == levels(group.micro_comp.dataset$intervals)[6], paste0(">5 (n=", format(group.micro_comp.dataset%>%filter(intervals==levels(group.micro_comp.dataset$intervals)[6])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.micro_comp.dataset%>%filter(intervals==levels(group.micro_comp.dataset$intervals)[6])%>%filter(postdrug_micro_comp_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"), intervals))))))) %>%
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
             ci.vertices.height = 0.3,
             boxsize = .4,
             lwd.ci = 3,
             txt_gp = fpTxtGp(ticks=gpar(cex=1.5), 
                              title=gpar(cex=2.5),
                              legend =gpar(cex=2),
                              label=gpar(cex=1.5),
                              xlab=gpar(cex=1.5)),
             xlab = "Hazard ratio (95% CI, log scale)") %>%
  fp_add_header(paste0("Overall population (n=", format(nrow(group.micro_comp.dataset),big.mark=",",scientific=FALSE), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = "black")

# PDF with plots
pdf(width = 27.07/2.54, height = 18.03/2.54, "Plots/DUK_2023/11.08.plot_duk_17.pdf")
plot_duk_17
dev.off()




#:------------------------------------------------------------------------------
# patient characteristics and results treatment effect differences - presentation

interval_breaks <- c(-5, -3, 0, 3, 5)



full.cohort <- set_up_data_sglt2_glp1(dataset.type = "hba1c.train") %>%
  # drop the variables with the most missingness (>40%)
  select(-preacr, -preast, -prehaematocrit, -prehaemoglobin, -pretriglyceride) %>%
  # only complete cases
  drop_na() %>%
  as.data.frame() %>%
  select(patid, pated, posthba1cfinal, drugclass, unique(c(readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/variables_mu.rds"), readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/variables_tau.rds")))) %>%
  rbind(
    set_up_data_sglt2_glp1(dataset.type = "hba1c.test") %>%
      # selected variables from SparseBCF
      select(patid, pated, posthba1cfinal, drugclass, unique(c(readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/variables_mu.rds"), readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/variables_tau.rds")))) %>%
      # only complete cases
      drop_na() %>%
      as.data.frame()
  ) %>%
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


plot_sex_strata <- group.full.cohort %>%
  select(sex, intervals) %>%
  ggplot(aes(x = intervals)) +
  geom_bar(aes(fill = sex), position = "fill") +
  ggtitle("Sex") +
  ylab("Percentage (%)") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Accent") +
  theme_bw() +
  annotate("point", x = 1, y = -0.1, shape = 15, colour = "#f1a340", size = 6, alpha = 1) +
  annotate("point", x = 2, y = -0.1, shape = 15, colour = "#f1a340", size = 6, alpha = 0.6) +
  annotate("point", x = 3, y = -0.1, shape = 15, colour = "#f1a340", size = 6, alpha = 0.2) +
  annotate("point", x = 4, y = -0.1, shape = 15, colour = "dodgerblue2", size = 6, alpha = 0.2) +
  annotate("point", x = 5, y = -0.1, shape = 15, colour = "dodgerblue2", size = 6, alpha = 0.6) +
  annotate("point", x = 6, y = -0.1, shape = 15, colour = "dodgerblue2", size = 6, alpha = 1) +
  coord_cartesian(ylim = c(0, 1), clip = "off") +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 30),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 30),
        legend.text = element_text(size= 30),
        legend.title = element_text(size = 30), 
        axis.title.x = element_text(vjust = -50, colour = "white"),
        strip.background = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5, size = 30)) +
  guides(fill=guide_legend(title="Sex", title.position="top", title.hjust = 0.5))

plot_prehba1c_strata <- group.full.cohort %>%
  select(intervals, prehba1c) %>%
  drop_na() %>%
  ggplot(aes(x = intervals, y = prehba1c)) +
  geom_point(aes(fill = intervals, colour = intervals), alpha = 0, shape = 15, size = 5) +
  geom_boxplot(outlier.shape = NA, aes(fill = intervals), alpha = c(1, 0.6, 0.2, 0.2, 0.6, 1), show.legend = F) +
  ggtitle("Baseline HbA1c") +
  ylab("Baseline HbA1c (mmol/mol)") +
  theme_bw() +
  annotate("point", x = 1, y = 43, shape = 15, colour = "#f1a340", size = 6, alpha = 1) +
  annotate("point", x = 2, y = 43, shape = 15, colour = "#f1a340", size = 6, alpha = 0.6) +
  annotate("point", x = 3, y = 43, shape = 15, colour = "#f1a340", size = 6, alpha = 0.2) +
  annotate("point", x = 4, y = 43, shape = 15, colour = "dodgerblue2", size = 6, alpha = 0.2) +
  annotate("point", x = 5, y = 43, shape = 15, colour = "dodgerblue2", size = 6, alpha = 0.6) +
  annotate("point", x = 6, y = 43, shape = 15, colour = "dodgerblue2", size = 6, alpha = 1) +
  coord_cartesian(ylim = c(50, 150), clip = "off") +
  scale_y_continuous(breaks = c(53, 64, 75, 86, 97, 108, 119, 150)) +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 30),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 30),
        legend.text = element_text(size= 30),
        legend.title = element_text(size = 30), 
        axis.title.x = element_text(vjust = -50, colour = "white"),
        strip.background = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5, size = 30),
        panel.grid.minor = element_blank()) +
  scale_fill_manual(
    name = paste0("HbA1c benefit (mmol/mol) [n=", format(nrow(group.full.cohort),big.mark=",",scientific=FALSE),"]"),
    labels = c(paste0("SGLT2i > 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[1]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("SGLT2i 3 - 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[2]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("SGLT2i 0 - 3 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[3]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("GLP1-RA > 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[4]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("GLP1-RA 3 - 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[5]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("GLP1-RA 0 - 3 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[6]) %>% nrow(),big.mark=",",scientific=FALSE), "]")),
    breaks = c(levels(group.full.cohort$intervals)[1], levels(group.full.cohort$intervals)[2], levels(group.full.cohort$intervals)[3], levels(group.full.cohort$intervals)[6], levels(group.full.cohort$intervals)[5], levels(group.full.cohort$intervals)[4]),
    values = c("SGLT2i benefit >5 mmol/mol" = "#f1a340",
               "SGLT2i benefit 3-5 mmol/mol" = "#f1a340",
               "SGLT2i benefit 0-3 mmol/mol" = "#f1a340",
               "GLP1-RA benefit 0-3 mmol/mol" = "dodgerblue2",
               "GLP1-RA benefit 3-5 mmol/mol" = "dodgerblue2",
               "GLP1-RA benefit >5 mmol/mol" = "dodgerblue2")
  ) +
  scale_colour_manual(
    name = paste0("HbA1c benefit (mmol/mol) [n=", format(nrow(group.full.cohort),big.mark=",",scientific=FALSE),"]"),
    labels = c(paste0("SGLT2i > 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[1]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("SGLT2i 3 - 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[2]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("SGLT2i 0 - 3 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[3]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("GLP1-RA > 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[4]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("GLP1-RA 3 - 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[5]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("GLP1-RA 0 - 3 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[6]) %>% nrow(),big.mark=",",scientific=FALSE), "]")),
    breaks = c(levels(group.full.cohort$intervals)[1], levels(group.full.cohort$intervals)[2], levels(group.full.cohort$intervals)[3], levels(group.full.cohort$intervals)[6], levels(group.full.cohort$intervals)[5], levels(group.full.cohort$intervals)[4]),
    values = c("SGLT2i benefit >5 mmol/mol" = "#f1a340",
               "SGLT2i benefit 3-5 mmol/mol" = "#f1a340",
               "SGLT2i benefit 0-3 mmol/mol" = "#f1a340",
               "GLP1-RA benefit 0-3 mmol/mol" = "dodgerblue2",
               "GLP1-RA benefit 3-5 mmol/mol" = "dodgerblue2",
               "GLP1-RA benefit >5 mmol/mol" = "dodgerblue2")
  ) +
  guides(
    fill = guide_legend(title.position="top", title.hjust = 0.5, override.aes = list(alpha = c(1, 0.6, 0.2, 1, 0.6, 0.2), size = 6), shape = c(15,15,15,15,15,15),  ncol = 2, nrow = 3),
    colour = guide_legend(title.position="top", title.hjust = 0.5, override.aes = list(alpha = c(1, 0.6, 0.2, 1, 0.6, 0.2), size = 6), shape = c(15,15,15,15,15,15),  ncol = 2, nrow = 3)
  )


plot_preegfr_strata <- group.full.cohort %>%
  select(intervals, preegfr) %>%
  drop_na() %>%
  ggplot(aes(x = intervals, y = preegfr)) +
  geom_point(aes(fill = intervals, colour = intervals), alpha = 0, shape = 15, size = 5) +
  geom_boxplot(outlier.shape = NA, aes(fill = intervals), alpha = c(1, 0.6, 0.2, 0.2, 0.6, 1), show.legend = F) +
  ggtitle("eGFR") +
  ylab(expression(paste("eGFR (mL/min/1.3", m^{2}, ")"))) +
  scale_y_continuous(breaks = c(30, 60, 90, 120, 150)) +
  theme_bw() +
  annotate("point", x = 1, y = 18, shape = 15, colour = "#f1a340", size = 6, alpha = 1) +
  annotate("point", x = 2, y = 18, shape = 15, colour = "#f1a340", size = 6, alpha = 0.6) +
  annotate("point", x = 3, y = 18, shape = 15, colour = "#f1a340", size = 6, alpha = 0.2) +
  annotate("point", x = 4, y = 18, shape = 15, colour = "dodgerblue2", size = 6, alpha = 0.2) +
  annotate("point", x = 5, y = 18, shape = 15, colour = "dodgerblue2", size = 6, alpha = 0.6) +
  annotate("point", x = 6, y = 18, shape = 15, colour = "dodgerblue2", size = 6, alpha = 1) +
  coord_cartesian(ylim = c(30, 150), clip = "off") +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 30),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 30),
        legend.text = element_text(size= 30),
        legend.title = element_text(size = 30), 
        axis.title.x = element_text(vjust = -50, colour = "white"),
        strip.background = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5, size = 30)) +
  scale_fill_manual(
    name = paste0("HbA1c benefit (mmol/mol) [n=", format(nrow(group.full.cohort),big.mark=",",scientific=FALSE),"]"),
    labels = c(paste0("SGLT2i > 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[1]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("SGLT2i 3 - 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[2]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("SGLT2i 0 - 3 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[3]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("GLP1-RA > 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[4]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("GLP1-RA 3 - 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[5]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("GLP1-RA 0 - 3 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[6]) %>% nrow(),big.mark=",",scientific=FALSE), "]")),
    breaks = c(levels(group.full.cohort$intervals)[1], levels(group.full.cohort$intervals)[2], levels(group.full.cohort$intervals)[3], levels(group.full.cohort$intervals)[6], levels(group.full.cohort$intervals)[5], levels(group.full.cohort$intervals)[4]),
    values = c("SGLT2i benefit >5 mmol/mol" = "#f1a340",
               "SGLT2i benefit 3-5 mmol/mol" = "#f1a340",
               "SGLT2i benefit 0-3 mmol/mol" = "#f1a340",
               "GLP1-RA benefit 0-3 mmol/mol" = "dodgerblue2",
               "GLP1-RA benefit 3-5 mmol/mol" = "dodgerblue2",
               "GLP1-RA benefit >5 mmol/mol" = "dodgerblue2")
  ) +
  scale_colour_manual(
    name = paste0("HbA1c benefit (mmol/mol) [n=", format(nrow(group.full.cohort),big.mark=",",scientific=FALSE),"]"),
    labels = c(paste0("SGLT2i > 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[1]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("SGLT2i 3 - 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[2]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("SGLT2i 0 - 3 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[3]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("GLP1-RA > 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[4]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("GLP1-RA 3 - 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[5]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("GLP1-RA 0 - 3 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[6]) %>% nrow(),big.mark=",",scientific=FALSE), "]")),
    breaks = c(levels(group.full.cohort$intervals)[1], levels(group.full.cohort$intervals)[2], levels(group.full.cohort$intervals)[3], levels(group.full.cohort$intervals)[6], levels(group.full.cohort$intervals)[5], levels(group.full.cohort$intervals)[4]),
    values = c("SGLT2i benefit >5 mmol/mol" = "#f1a340",
               "SGLT2i benefit 3-5 mmol/mol" = "#f1a340",
               "SGLT2i benefit 0-3 mmol/mol" = "#f1a340",
               "GLP1-RA benefit 0-3 mmol/mol" = "dodgerblue2",
               "GLP1-RA benefit 3-5 mmol/mol" = "dodgerblue2",
               "GLP1-RA benefit >5 mmol/mol" = "dodgerblue2")
  ) +
  guides(
    fill = guide_legend(title.position="top", title.hjust = 0.5, override.aes = list(alpha = c(1, 0.6, 0.2, 1, 0.6, 0.2), size = 6), shape = c(15,15,15,15,15,15),  ncol = 2, nrow = 3),
    colour = guide_legend(title.position="top", title.hjust = 0.5, override.aes = list(alpha = c(1, 0.6, 0.2, 1, 0.6, 0.2), size = 6), shape = c(15,15,15,15,15,15),  ncol = 2, nrow = 3)
  )



plot_agetx_strata <- group.full.cohort %>%
  select(intervals, agetx) %>%
  drop_na() %>%
  ggplot(aes(x = intervals, y = agetx)) +
  geom_point(aes(fill = intervals, colour = intervals), alpha = 0, shape = 15, size = 5) +
  geom_boxplot(outlier.shape = NA, aes(fill = intervals), alpha = c(1, 0.6, 0.2, 0.2, 0.6, 1), show.legend = F) +
  ggtitle("Current age") +
  ylab("Current age (years)") +
  theme_bw() +
  annotate("point", x = 1, y = 10, shape = 15, colour = "#f1a340", size = 6, alpha = 1) +
  annotate("point", x = 2, y = 10, shape = 15, colour = "#f1a340", size = 6, alpha = 0.6) +
  annotate("point", x = 3, y = 10, shape = 15, colour = "#f1a340", size = 6, alpha = 0.2) +
  annotate("point", x = 4, y = 10, shape = 15, colour = "dodgerblue2", size = 6, alpha = 0.2) +
  annotate("point", x = 5, y = 10, shape = 15, colour = "dodgerblue2", size = 6, alpha = 0.6) +
  annotate("point", x = 6, y = 10, shape = 15, colour = "dodgerblue2", size = 6, alpha = 1) +
  coord_cartesian(ylim = c(18, 100), clip = "off") +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 30),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 30),
        legend.text = element_text(size= 30),
        legend.title = element_text(size = 30), 
        axis.title.x = element_text(vjust = -50, colour = "white"),
        strip.background = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5, size = 30)) +
  scale_fill_manual(
    name = paste0("HbA1c benefit (mmol/mol) [n=", format(nrow(group.full.cohort),big.mark=",",scientific=FALSE),"]"),
    labels = c(paste0("SGLT2i > 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[1]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("SGLT2i 3 - 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[2]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("SGLT2i 0 - 3 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[3]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("GLP1-RA > 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[4]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("GLP1-RA 3 - 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[5]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("GLP1-RA 0 - 3 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[6]) %>% nrow(),big.mark=",",scientific=FALSE), "]")),
    breaks = c(levels(group.full.cohort$intervals)[1], levels(group.full.cohort$intervals)[2], levels(group.full.cohort$intervals)[3], levels(group.full.cohort$intervals)[6], levels(group.full.cohort$intervals)[5], levels(group.full.cohort$intervals)[4]),
    values = c("SGLT2i benefit >5 mmol/mol" = "#f1a340",
               "SGLT2i benefit 3-5 mmol/mol" = "#f1a340",
               "SGLT2i benefit 0-3 mmol/mol" = "#f1a340",
               "GLP1-RA benefit 0-3 mmol/mol" = "dodgerblue2",
               "GLP1-RA benefit 3-5 mmol/mol" = "dodgerblue2",
               "GLP1-RA benefit >5 mmol/mol" = "dodgerblue2")
  ) +
  scale_colour_manual(
    name = paste0("HbA1c benefit (mmol/mol) [n=", format(nrow(group.full.cohort),big.mark=",",scientific=FALSE),"]"),
    labels = c(paste0("SGLT2i > 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[1]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("SGLT2i 3 - 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[2]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("SGLT2i 0 - 3 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[3]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("GLP1-RA > 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[4]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("GLP1-RA 3 - 5 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[5]) %>% nrow(),big.mark=",",scientific=FALSE), "]"),
               paste0("GLP1-RA 0 - 3 [n=", format(group.full.cohort %>% filter(intervals == levels(group.full.cohort$intervals)[6]) %>% nrow(),big.mark=",",scientific=FALSE), "]")),
    breaks = c(levels(group.full.cohort$intervals)[1], levels(group.full.cohort$intervals)[2], levels(group.full.cohort$intervals)[3], levels(group.full.cohort$intervals)[6], levels(group.full.cohort$intervals)[5], levels(group.full.cohort$intervals)[4]),
    values = c("SGLT2i benefit >5 mmol/mol" = "#f1a340",
               "SGLT2i benefit 3-5 mmol/mol" = "#f1a340",
               "SGLT2i benefit 0-3 mmol/mol" = "#f1a340",
               "GLP1-RA benefit 0-3 mmol/mol" = "dodgerblue2",
               "GLP1-RA benefit 3-5 mmol/mol" = "dodgerblue2",
               "GLP1-RA benefit >5 mmol/mol" = "dodgerblue2")
  ) +
  guides(
    fill = guide_legend(title.position="top", title.hjust = 0.5, override.aes = list(alpha = c(1, 0.6, 0.2, 1, 0.6, 0.2), size = 6), shape = c(15,15,15,15,15,15),  ncol = 2, nrow = 3),
    colour = guide_legend(title.position="top", title.hjust = 0.5, override.aes = list(alpha = c(1, 0.6, 0.2, 1, 0.6, 0.2), size = 6), shape = c(15,15,15,15,15,15),  ncol = 2, nrow = 3)
  )



plot_duk_18 <- patchwork::wrap_plots(list(
  plot_prehba1c_strata,
  plot_preegfr_strata,
  plot_agetx_strata,
  plot_sex_strata
), nrow = 1) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "none",
    # legend.text = element_text(size = 18),
    # legend.title = element_text(size = 18),
    plot.title = element_text(hjust = 0.5, size = 45),
    axis.text.y = element_text(size = 40),
    axis.title.y = element_text(size = 40)
  )

pdf(width = 34/2.54*2, height = 26.88/2.54, "Plots/DUK_2023/11.08.plot_duk_18.pdf")
plot_duk_18
dev.off()


