####################
## Description:
##  - In this file we compare the predictions of SGLT2/GLP1 from the BCF model
##    to the predictions form the SGLT2/DPP4 linear model.
####################

## Load libraries
library(rms)
library(tidyverse)
library(scales)

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
source("02.slade_aurum_set_data.R")

###############################################################################
###############################################################################
########################## General variables ##################################
###############################################################################
###############################################################################

# read in predictions from BCF model
patient_predicted_outcomes <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/patient_predicted_outcomes.rds")

# load the object for the Linear model. This model object is the object used
#   in the SGLT2vsDPP4 model by John Dennis et al.
# Regression parameters can be extracted from the original paper.
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

# predict for DPP4
dataset.dpp4 <- full.cohort.updated %>%
  select(-patid, -pated) %>%
  mutate(drugclass = factor("DPP4", levels = c("SGLT2", "DPP4")))
predictions.dpp4 <- predict(m1, dataset.dpp4)

# treatment effects
effects <- predictions.sglt2-predictions.dpp4



## Compare SGLT2 BCF vs Linear regression
interim.dataset <- full.cohort.updated %>%
  select(patid, pated) %>%
  cbind(pred.SGLT2.lm = predictions.sglt2) %>%
  left_join(patient_predicted_outcomes %>%
              select(patid, pated, pred.SGLT2) %>%
              rename("pred.SGLT2.bcf" = "pred.SGLT2"), by = c("patid", "pated")) %>%
  select(-patid, -pated)

# Plot the comparison of SGLT2 predictions from the linear SGLT2vsDPP4 model vs SGLT2vsGLP1 model
plot_comparison <- interim.dataset %>%
  ggplot(aes(y = pred.SGLT2.bcf, x = pred.SGLT2.lm)) +
  geom_point() +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red", linetype = "dashed") +
  stat_smooth() +
  ylab("SGLT2 predictions with BCF") +
  xlab("SGLT2 predictions with Linear Regression") +
  ylim(min(interim.dataset), max(interim.dataset)) +
  xlim(min(interim.dataset), max(interim.dataset)) +
  ggtitle(paste0("Comparison of SGLT2 predictions using Linear Regression and BCF (n=", interim.dataset %>% nrow(), ")")) +
  theme_bw()

## PDF with the plot for the comparison  
pdf(width = 7, height = 7, "Plots/11.09.plot_1.pdf")
plot_comparison
dev.off()

## r_squared 

## Calculate assessments of prediction
rsq <- function (x, y) cor(x, y) ^ 2

r_squared_value <- rsq(interim.dataset$pred.SGLT2.lm, interim.dataset$pred.SGLT2.bcf)


## What drug is best:
#-----------------
# Using SGLT2 BCF
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
plot_bar <- interim.dataset %>%
  select(best_drug) %>%
  table() %>%
  as.data.frame() %>%
  rename("best_drug" = ".") %>%
  ggplot(aes(x = best_drug, y = Freq, fill = best_drug)) +
  geom_bar(stat="identity") +
  geom_text(aes(label = format(Freq,big.mark=",",scientific=FALSE)), vjust = -0.5) +
  xlab("Optimal predicted therapy") +
  ylab("Number of patients") +
  ggtitle(paste0("Predicted optimal therapy (n=", format(interim.dataset%>%nrow(),big.mark=",",scientific=FALSE), ")")) +
  scale_fill_manual(values = c("red", "dodgerblue2", "#f1a340")) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_y_continuous(label=comma)

## PDF with the plot for the best therapy
pdf(width = 7, height = 7, "Plots/11.09.plot_2.pdf")
plot_bar
dev.off()


#---------------
## Plot the differential treatment effect: best vs closest other
interim.dataset <- patient_predicted_outcomes %>%
  left_join(full.cohort.updated %>%
              select(patid, pated) %>%
              cbind(pred.DPP4 = predictions.dpp4), by = c("patid", "pated")) %>%
  drop_na() %>%
  mutate(best_drug = ifelse(pred.SGLT2 < pred.GLP1 & pred.SGLT2 < pred.DPP4, "Favours SGLT2i",
                            ifelse(pred.GLP1 < pred.SGLT2 & pred.GLP1 < pred.DPP4, "Favours GLP1-RA",
                                   ifelse(pred.DPP4 < pred.SGLT2 & pred.DPP4 < pred.GLP1, "Favours DPP4i", NA))),
         best_drug = factor(best_drug),
         effect = ifelse(best_drug == "Favours SGLT2i" & pred.GLP1 < pred.DPP4, pred.SGLT2 - pred.GLP1,
                         ifelse(best_drug == "Favours SGLT2i" & pred.DPP4 < pred.GLP1, pred.SGLT2 - pred.DPP4,
                                ifelse(best_drug == "Favours GLP1-RA" & pred.SGLT2 < pred.DPP4, pred.GLP1 - pred.SGLT2,
                                       ifelse(best_drug == "Favours GLP1-RA" & pred.DPP4 < pred.SGLT2, pred.GLP1 - pred.DPP4,
                                              ifelse(best_drug == "Favours DPP4i" & pred.SGLT2 < pred.GLP1, pred.DPP4 - pred.SGLT2,
                                                     ifelse(best_drug == "Favours DPP4i" & pred.GLP1 < pred.SGLT2, pred.DPP4 - pred.GLP1, NA))))))) 

# Plot histogram of predicted differential treatment effect: best vs closest other
plot_histogram <- interim.dataset %>%
  select(best_drug, effect) %>%
  ggplot(aes(x = effect, fill = best_drug)) +
  geom_histogram(position = "identity", alpha = 0.5, color = "black", breaks = seq(-15, 0, by = 1)) +
  geom_vline(aes(xintercept = 0), linetype = "dashed") +
  xlab("Predicted HbA1c benefit (mmol/mol)") +
  ggtitle(paste0("Predicted HbA1c benefit against the next best therapy (n=", interim.dataset %>% nrow(), ")")) +
  scale_fill_manual(values = c("red", "dodgerblue2", "#f1a340")) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.y = element_blank()) +
  facet_wrap(~best_drug, nrow = 1) +
  scale_y_continuous(label=comma)

## PDF with the plot for the differential treatment effects
pdf(width = 9, height = 4, "Plots/11.09.plot_3.pdf")
plot_histogram
dev.off()


#---------------
## Decision tree for the covariates used

# load in variables used in the model
variables_mu <- readRDS(paste0(output_path, "/response_model_bcf/variables_mu.rds"))

variables_tau <- readRDS(paste0(output_path, "/response_model_bcf/variables_tau.rds"))

interim.dataset <- patient_predicted_outcomes %>%
  left_join(full.cohort.updated %>%
              select(patid, pated) %>%
              cbind(pred.DPP4 = predictions.dpp4), by = c("patid", "pated")) %>%
  drop_na() %>%
  mutate(best_drug = ifelse(pred.SGLT2 < pred.GLP1 & pred.SGLT2 < pred.DPP4, "Favours SGLT2i",
                            ifelse(pred.GLP1 < pred.SGLT2 & pred.GLP1 < pred.DPP4, "Favours GLP1-RA",
                                   ifelse(pred.DPP4 < pred.SGLT2 & pred.DPP4 < pred.GLP1, "Favours DPP4i", NA))),
         best_drug = factor(best_drug),
         effect = ifelse(best_drug == "Favours SGLT2i" & pred.GLP1 < pred.DPP4, pred.SGLT2 - pred.GLP1,
                         ifelse(best_drug == "Favours SGLT2i" & pred.DPP4 < pred.GLP1, pred.SGLT2 - pred.DPP4,
                                ifelse(best_drug == "Favours GLP1-RA" & pred.SGLT2 < pred.DPP4, pred.GLP1 - pred.SGLT2,
                                       ifelse(best_drug == "Favours GLP1-RA" & pred.DPP4 < pred.SGLT2, pred.GLP1 - pred.DPP4,
                                              ifelse(best_drug == "Favours DPP4i" & pred.SGLT2 < pred.GLP1, pred.DPP4 - pred.SGLT2,
                                                     ifelse(best_drug == "Favours DPP4i" & pred.GLP1 < pred.SGLT2, pred.DPP4 - pred.GLP1, NA))))))) %>%
  select(patid, pated, best_drug)

full.cohort.decision <- set_up_data_sglt2_glp1(dataset.type = "full.cohort") %>%
  select(all_of(c("patid", "pated", unique(c(variables_mu, variables_tau))))) %>%
  left_join(interim.dataset, by = c("patid", "pated")) %>%
  drop_na(best_drug) %>%
  select(-hba1cmonth, -patid, -pated)

# Load libraries
library(rpart)
library(rattle)
library(rpart.plot)


fit <- rpart(best_drug ~ agetx + t2dmduration + prehba1c + preegfr + prealt + prepad + sex + prebmi + preheartfailure + preihd + preneuropathy + preretinopathy, data = full.cohort.decision)

pdf(width = 10, height = 8, file = "Plots/11.09.plot_4.pdf")

prp(fit, pal.thresh = 0, extra = "auto", main = "Decision tree for treatment effects using development cohort", box.palette = list("dodgerblue2", "#f1a340"))

dev.off()




#---------------
## Table of characteristics

table.best_drug <- tableone::CreateTableOne(vars = c("agetx", "t2dmduration", "prehba1c", "preegfr", "prealt", "prepad", "sex", "prebmi", "preheartfailure", "preihd", "preneuropathy", "preretinopathy"), includeNA = TRUE, strata = "best_drug", data = full.cohort.decision, test = FALSE)


# print(table.best_drug)



## merge pdfs

qpdf::pdf_combine(input = c("Plots/11.09.plot_1.pdf",
                            "Plots/11.09.plot_2.pdf",
                            "Plots/11.09.plot_3.pdf",
                            "Plots/11.09.plot_4.pdf"),
                  output = "Plots/11.09.comparison_SGLT2_GLP1_DPP4.pdf")

# delete old pdfs
file.remove(c("Plots/11.09.plot_1.pdf", "Plots/11.09.plot_2.pdf", "Plots/11.09.plot_3.pdf", "Plots/11.09.plot_4.pdf"))


