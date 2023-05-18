####################
## Description:
##  - In this file we calculate the treatment effect / response for 
##      the average patient for a range of values at each covariate.
####################


## Load libraries
library(tidyverse)


## Set up directory path to save files (stagered to ensure folders are created)

dir.create("Samples")
dir.create("Samples/SGLT2-GLP1")

output_path <- "Samples/SGLT2-GLP1/Aurum"
dir.create(output_path)

## make directory for outputs
dir.create(paste0(output_path, "/differential_response"))

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

# variables chosen
variables_mu <- readRDS(paste0(output_path, "/response_model_bcf/variables_mu.rds"))
variables_tau <- readRDS(paste0(output_path, "/response_model_bcf/variables_tau.rds"))

# treatment effects
patient_effects <- readRDS(paste0(output_path, "/response_model_bcf/patient_effects.rds"))

# Full cohort for average values
hba1c.train <- set_up_data_sglt2_glp1(dataset.type="hba1c.train") %>%
  left_join(patient_effects, by = c("patid", "pated"))

levels(hba1c.train$sex) <- c("Females", "Males")
levels(hba1c.train$ncurrtx) <- c("0", "1", "2", "3", "4+")



#:------------------------------------------------------------------------
# Stratify predicted treatment effects by variables
# sex
plot_sex_strata <- hba1c.train %>%
  select(sex, effects) %>%
  ggplot(aes(x = sex, y = effects)) +
  geom_hline(aes(yintercept = 0), colour = "red") +
  geom_boxplot(outlier.shape = NA) +
  ylim(-10, 10) +
  ggtitle("Sex") +
  ylab("Predicted treatment effects (mmol/mol)") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5))

# number of other current drugs
plot_ncurrtx_strata <- hba1c.train %>%
  select(ncurrtx, effects, sex) %>%
  ggplot(aes(x = ncurrtx, y = effects)) +
  geom_hline(aes(yintercept = 0), colour = "red") +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~sex) +
  ylim(-10, 10) +
  ggtitle("Number of other current\nglucose-lowering drugs") +
  ylab("Predicted treatment effects (mmol/mol)") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(fill="white"))

# Peripheral arterial disease
plot_prepad_strata <- hba1c.train %>%
  select(prepad, effects, sex) %>%
  ggplot(aes(x = prepad, y = effects)) +
  geom_hline(aes(yintercept = 0), colour = "red") +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~sex) +
  ylim(-10, 10) +
  ggtitle("Peripheral arterial disease") +
  ylab("Predicted treatment effects (mmol/mol)") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(fill="white"))

# Ischaemic heart disease
plot_preihd_strata <- hba1c.train %>%
  select(preihd, effects, sex) %>%
  ggplot(aes(x = preihd, y = effects)) +
  geom_hline(aes(yintercept = 0), colour = "red") +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~sex) +
  ylim(-10, 10) +
  ggtitle("Ischaemic heart disease") +
  ylab("Predicted treatment effects (mmol/mol)") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(fill="white"))

# Neuropathy
plot_preneuropathy_strata <- hba1c.train %>%
  select(preneuropathy, effects, sex) %>%
  ggplot(aes(x = preneuropathy, y = effects)) +
  geom_hline(aes(yintercept = 0), colour = "red") +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~sex) +
  ylim(-10, 10) +
  ggtitle("Neuropathy") +
  ylab("Predicted treatment effects (mmol/mol)") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(fill="white"))

# Retinopathy
plot_preretinopathy_strata <- hba1c.train %>%
  select(preretinopathy, effects, sex) %>%
  ggplot(aes(x = preretinopathy, y = effects)) +
  geom_hline(aes(yintercept = 0), colour = "red") +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~sex) +
  ylim(-10, 10) +
  ggtitle("Retinopathy") +
  ylab("Predicted treatment effects (mmol/mol)") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(fill="white"))

# Heart Failure
plot_preheartfailure_strata <- hba1c.train %>%
  select(preheartfailure, effects, sex) %>%
  ggplot(aes(x = preheartfailure, y = effects)) +
  geom_hline(aes(yintercept = 0), colour = "red") +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~sex) +
  ylim(-10, 10) +
  ggtitle("Heart failure") +
  ylab("Predicted treatment effects (mmol/mol)") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(fill="white"))

# HbA1c at baseline
breaks_hba1c <- quantile(hba1c.train$prehba1c, probs = seq(0.2, 0.9, 0.2), na.rm = TRUE)

plot_prehba1c_strata <- group_values(data = hba1c.train,
                                     variable = "prehba1c",
                                     breaks = breaks_hba1c) %>%
  select(intervals, effects, sex, prehba1c) %>%
  drop_na() %>%
  group_by(intervals) %>%
  mutate(intervals_labels = round(mean(prehba1c, na.rm = TRUE)),
         intervals_labels = factor(intervals_labels)) %>%
  ungroup() %>%
  ggplot(aes(x = intervals_labels, y = effects)) +
  geom_hline(aes(yintercept = 0), colour = "red") +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~sex) +
  ylim(-10, 10) +
  ggtitle("HbA1c") +
  ylab("Predicted treatment effects (mmol/mol)") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(fill="white"))

# eGFR at baseline
breaks_egfr <- quantile(hba1c.train$preegfr, probs = seq(0.2, 0.9, 0.2), na.rm = TRUE)

plot_preegfr_strata <- group_values(data = hba1c.train,
                                    variable = "preegfr",
                                    breaks = breaks_egfr) %>%
  select(intervals, effects, sex, preegfr) %>%
  drop_na() %>%
  group_by(intervals) %>%
  mutate(intervals_labels = round(mean(preegfr, na.rm = TRUE)),
         intervals_labels = factor(intervals_labels)) %>%
  ungroup() %>%
  ggplot(aes(x = intervals_labels, y = effects)) +
  geom_hline(aes(yintercept = 0), colour = "red") +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~sex) +
  ylim(-10, 10) +
  ggtitle("eGFR") +
  ylab("Predicted treatment effects (mmol/mol)") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(fill="white"))

# Age at baseline
breaks_agetx <- quantile(hba1c.train$agetx, probs = seq(0.2, 0.9, 0.2), na.rm = TRUE)

plot_agetx_strata <- group_values(data = hba1c.train,
                                  variable = "agetx",
                                  breaks = breaks_agetx) %>%
  select(intervals, effects, sex, agetx) %>%
  drop_na() %>%
  group_by(intervals) %>%
  mutate(intervals_labels = round(mean(agetx, na.rm = TRUE)),
         intervals_labels = factor(intervals_labels)) %>%
  ungroup() %>%
  ggplot(aes(x = intervals_labels, y = effects)) +
  geom_hline(aes(yintercept = 0), colour = "red") +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~sex) +
  ylim(-10, 10) +
  ggtitle("Current age") +
  ylab("Predicted treatment effects (mmol/mol)") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(fill="white"))

# BMI at baseline
breaks_prebmi <- quantile(hba1c.train$prebmi, probs = seq(0.2, 0.9, 0.2), na.rm = TRUE)

plot_prebmi_strata <- group_values(data = hba1c.train,
                                   variable = "prebmi",
                                   breaks = breaks_prebmi) %>%
  select(intervals, effects, sex, prebmi) %>%
  drop_na() %>%
  group_by(intervals) %>%
  mutate(intervals_labels = round(mean(prebmi, na.rm = TRUE)),
         intervals_labels = factor(intervals_labels)) %>%
  ungroup() %>%
  ggplot(aes(x = intervals_labels, y = effects)) +
  geom_hline(aes(yintercept = 0), colour = "red") +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~sex) +
  ylim(-10, 10) +
  ggtitle("BMI") +
  ylab("Predicted treatment effects (mmol/mol)") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(fill="white"))


# Combine all
plot_strata <- patchwork::wrap_plots(
  list(
    plot_ncurrtx_strata,
    plot_sex_strata,
    plot_preegfr_strata,
    plot_agetx_strata,
    plot_prebmi_strata +
      theme(axis.title.y = element_text(size = 11)),
    plot_prehba1c_strata,
    plot_preretinopathy_strata,
    plot_prepad_strata,
    plot_preneuropathy_strata,
    plot_preihd_strata,
    plot_preheartfailure_strata
  )) +
  patchwork::plot_annotation(
    title = "Boxplot of treatment effects for covariate strata"
  )


# PDF containing the plot
pdf(width = 16, height = 10, "Plots/11.07.diff_treat_effect.pdf")
plot_strata
dev.off()


