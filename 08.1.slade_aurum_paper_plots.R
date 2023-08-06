####################
## Description:
##  - In this file we create the plots used in the paper + extended figures
####################


## Load libraries
library(tidyverse)
library(patchwork)
library(scales)
library(forestplot)

## Load functions required
source("11.01.slade_aurum_functions.R")
source("11.02.slade_aurum_set_data.R")

## make directory for outputs
dir.create("Plots/Paper")

## make directory for outputs
dir.create("Plots/Paper/Manuscript")


#### Main manuscript plots

#:------------------------------------------------------------------------------
# Figure 1: 3 panel plot: 
#   -A. Histogram of treatment effects - development cohort
#   -B. Calibration plot for development cohort
#   -C. Calibration plot for validation cohort

### A

bcf_model <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/bcf_model.rds")

data_dev <- cbind(mean = colMeans(bcf_model$tau)) %>%
  as.data.frame()

plot_effect_dev <- hist_plot(data_dev, paste0("Development cohort (n=", format(nrow(data_dev),big.mark=",",scientific=FALSE), ")"), -15, 20) +
  scale_y_continuous(label=comma) +
  xlab("Predicted conditional average treatment effect (mmol/mol)")

### B

# load development cohort model calibration
ATE_adjust_validation_dev <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/assessment/ATE_adjust_validation_dev.rds")

# plot development cohort calibration
plot_ATE_adjust_validation_dev <- ATE_plot(ATE_adjust_validation_dev[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -12, 12) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle(paste0("Development cohort (n=", format(sum(ATE_adjust_validation_dev[["effects"]]$N),big.mark=",",scientific=FALSE), ")"))

### C

# load validation cohort model calibration
ATE_adjust_validation_val <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/assessment/ATE_adjust_validation_val.rds")

# plot validation cohort calibration
plot_ATE_adjust_validation_val <- ATE_plot(ATE_adjust_validation_val[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -12, 12) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle(paste0("Validation cohort (n=", format(sum(ATE_adjust_validation_val[["effects"]]$N),big.mark=",",scientific=FALSE), ")"))


plot_main_1.1 <- plot_effect_dev | plot_ATE_adjust_validation_dev | plot_ATE_adjust_validation_val

plot_main_1 <- plot_main_1.1 +
  plot_annotation(tag_levels = "A",
                  theme = theme(legend.position = "bottom"))

pdf(width = 15, height = 6, "Plots/Paper/Manuscript/11.08.plot_main_1.pdf")
plot_main_1
dev.off()



#:------------------------------------------------------------------------------
# Figure 2:
#:--------------------------------------------------------------------------------
### Variables ranges across HbA1c benefit subgroups 

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


plot_t2dmduration_strata <- group.full.cohort %>%
  select(intervals, t2dmduration) %>%
  drop_na() %>%
  ggplot(aes(x = intervals, y = t2dmduration)) +
  geom_boxplot(outlier.shape = NA, fill = c(rep("#f1a340", 3), rep("dodgerblue2", 3)), alpha = c(1, 0.6, 0.2, 0.2, 0.6, 1)) +
  ggtitle("Duration of diabetes") +
  ylim(0, 30) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 70, vjust = 1, hjust = 1),
        strip.background = element_rect(fill = "white"),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))

# plot_ncurrtx_strata <- group.full.cohort %>%
#   select(ncurrtx, intervals) %>%
#   ggplot(aes(x = intervals)) +
#   geom_bar(aes(fill = ncurrtx), position = "fill") +
#   ggtitle("Number of other current\nglucose-lowering drugs") +
#   scale_y_continuous(labels = scales::percent) +
#   theme_bw() +
#   theme(legend.position = "bottom",
#         axis.title = element_blank(),
#         axis.text.x = element_text(angle = 70, vjust = 1, hjust = 1),
#         strip.background = element_rect(fill = "white"),
#         legend.title = element_blank(),
#         plot.title = element_text(hjust = 0.5))

plot_sex_strata <- group.full.cohort %>%
  select(sex, intervals) %>%
  ggplot(aes(x = intervals)) +
  geom_bar(aes(fill = sex), position = "fill") +
  ggtitle("Sex") +
  ylab("Percentage (%)") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Accent") +
  theme_bw() +
  annotate("point", x = 1, y = -0.08, shape = 15, colour = "#f1a340", size = 5, alpha = 1) +
  annotate("point", x = 2, y = -0.08, shape = 15, colour = "#f1a340", size = 5, alpha = 0.6) +
  annotate("point", x = 3, y = -0.08, shape = 15, colour = "#f1a340", size = 5, alpha = 0.2) +
  annotate("point", x = 4, y = -0.08, shape = 15, colour = "dodgerblue2", size = 5, alpha = 0.2) +
  annotate("point", x = 5, y = -0.08, shape = 15, colour = "dodgerblue2", size = 5, alpha = 0.6) +
  annotate("point", x = 6, y = -0.08, shape = 15, colour = "dodgerblue2", size = 5, alpha = 1) +
  coord_cartesian(ylim = c(0, 1), clip = "off") +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 13),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 13),
        legend.text = element_text(size= 13 ),
        legend.title = element_text(size = 13), 
        axis.title.x = element_text(vjust = -5, colour = "white"),
        strip.background = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5)) +
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
  annotate("point", x = 1, y = 42, shape = 15, colour = "#f1a340", size = 5, alpha = 1) +
  annotate("point", x = 2, y = 42, shape = 15, colour = "#f1a340", size = 5, alpha = 0.6) +
  annotate("point", x = 3, y = 42, shape = 15, colour = "#f1a340", size = 5, alpha = 0.2) +
  annotate("point", x = 4, y = 42, shape = 15, colour = "dodgerblue2", size = 5, alpha = 0.2) +
  annotate("point", x = 5, y = 42, shape = 15, colour = "dodgerblue2", size = 5, alpha = 0.6) +
  annotate("point", x = 6, y = 42, shape = 15, colour = "dodgerblue2", size = 5, alpha = 1) +
  coord_cartesian(ylim = c(50, 150), clip = "off") +
  scale_y_continuous(breaks = c(53, 64, 75, 86, 97, 108, 119, 150)) +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 13),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 13),
        legend.text = element_text(size= 13 ),
        legend.title = element_text(size = 13), 
        axis.title.x = element_text(vjust = -5, colour = "white"),
        strip.background = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5),
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
    fill = guide_legend(title.position="top", title.hjust = 0.5, override.aes = list(alpha = c(1, 0.6, 0.2, 1, 0.6, 0.2)), shape = c(15,15,15,15,15,15),  ncol = 2, nrow = 3),
    colour = guide_legend(title.position="top", title.hjust = 0.5, override.aes = list(alpha = c(1, 0.6, 0.2, 1, 0.6, 0.2)), shape = c(15,15,15,15,15,15),  ncol = 2, nrow = 3)
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
  annotate("point", x = 1, y = 19.5, shape = 15, colour = "#f1a340", size = 5, alpha = 1) +
  annotate("point", x = 2, y = 19.5, shape = 15, colour = "#f1a340", size = 5, alpha = 0.6) +
  annotate("point", x = 3, y = 19.5, shape = 15, colour = "#f1a340", size = 5, alpha = 0.2) +
  annotate("point", x = 4, y = 19.5, shape = 15, colour = "dodgerblue2", size = 5, alpha = 0.2) +
  annotate("point", x = 5, y = 19.5, shape = 15, colour = "dodgerblue2", size = 5, alpha = 0.6) +
  annotate("point", x = 6, y = 19.5, shape = 15, colour = "dodgerblue2", size = 5, alpha = 1) +
  coord_cartesian(ylim = c(30, 150), clip = "off") +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 13),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 13),
        legend.text = element_text(size= 13 ),
        legend.title = element_text(size = 13), 
        axis.title.x = element_text(vjust = -5, colour = "white"),
        strip.background = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5)) +
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
    fill = guide_legend(title.position="top", title.hjust = 0.5, override.aes = list(alpha = c(1, 0.6, 0.2, 1, 0.6, 0.2)), shape = c(15,15,15,15,15,15),  ncol = 2, nrow = 3),
    colour = guide_legend(title.position="top", title.hjust = 0.5, override.aes = list(alpha = c(1, 0.6, 0.2, 1, 0.6, 0.2)), shape = c(15,15,15,15,15,15),  ncol = 2, nrow = 3)
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
  annotate("point", x = 1, y = 11, shape = 15, colour = "#f1a340", size = 5, alpha = 1) +
  annotate("point", x = 2, y = 11, shape = 15, colour = "#f1a340", size = 5, alpha = 0.6) +
  annotate("point", x = 3, y = 11, shape = 15, colour = "#f1a340", size = 5, alpha = 0.2) +
  annotate("point", x = 4, y = 11, shape = 15, colour = "dodgerblue2", size = 5, alpha = 0.2) +
  annotate("point", x = 5, y = 11, shape = 15, colour = "dodgerblue2", size = 5, alpha = 0.6) +
  annotate("point", x = 6, y = 11, shape = 15, colour = "dodgerblue2", size = 5, alpha = 1) +
  coord_cartesian(ylim = c(18, 100), clip = "off") +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 13),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 13),
        legend.text = element_text(size= 13 ),
        legend.title = element_text(size = 13), 
        axis.title.x = element_text(vjust = -5, colour = "white"),
        strip.background = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5)) +
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
    fill = guide_legend(title.position="top", title.hjust = 0.5, override.aes = list(alpha = c(1, 0.6, 0.2, 1, 0.6, 0.2)), shape = c(15,15,15,15,15,15),  ncol = 2, nrow = 3),
    colour = guide_legend(title.position="top", title.hjust = 0.5, override.aes = list(alpha = c(1, 0.6, 0.2, 1, 0.6, 0.2)), shape = c(15,15,15,15,15,15),  ncol = 2, nrow = 3)
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
  annotate("point", x = 1, y = 11.5, shape = 15, colour = "#f1a340", size = 5, alpha = 1) +
  annotate("point", x = 2, y = 11.5, shape = 15, colour = "#f1a340", size = 5, alpha = 0.6) +
  annotate("point", x = 3, y = 11.5, shape = 15, colour = "#f1a340", size = 5, alpha = 0.2) +
  annotate("point", x = 4, y = 11.5, shape = 15, colour = "dodgerblue2", size = 5, alpha = 0.2) +
  annotate("point", x = 5, y = 11.5, shape = 15, colour = "dodgerblue2", size = 5, alpha = 0.6) +
  annotate("point", x = 6, y = 11.5, shape = 15, colour = "dodgerblue2", size = 5, alpha = 1) +
  coord_cartesian(ylim = c(15, 55), clip = "off") +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 13),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 13),
        legend.text = element_text(size= 13 ),
        legend.title = element_text(size = 13), 
        axis.title.x = element_text(vjust = -5, colour = "white"),
        strip.background = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5)) +
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
    fill = guide_legend(title.position="top", title.hjust = 0.5, override.aes = list(alpha = c(1, 0.6, 0.2, 1, 0.6, 0.2)), shape = c(15,15,15,15,15,15),  ncol = 2, nrow = 3),
    colour = guide_legend(title.position="top", title.hjust = 0.5, override.aes = list(alpha = c(1, 0.6, 0.2, 1, 0.6, 0.2)), shape = c(15,15,15,15,15,15),  ncol = 2, nrow = 3)
  )



plot_main_2 <- patchwork::wrap_plots(list(
  plot_sex_strata +
    theme(
      legend.position = "bottom"
    ), 
  
  patchwork::wrap_plots(list(
    
    plot_prehba1c_strata, 
    plot_preegfr_strata, 
    plot_agetx_strata, 
    plot_prebmi_strata
    
  ), nrow = 1) +
    plot_layout(guides = "collect") &
    theme(
      legend.position = "bottom"
    )
), nrow = 1) +
  plot_layout(
    widths = c(1, 5)
  )

pdf(width = 17, height = 8, "Plots/Paper/Manuscript/11.08.plot_main_2.pdf")
plot_main_2
dev.off()


# pdf(width = 5, height = 7, "Plots/Paper/Sup_Mat/11.08.plot_sup_15.pdf")
# plot_t2dmduration_strata
# dev.off()


#####################


plot_sex_strata <- group.full.cohort %>%
  select(sex, intervals) %>%
  ggplot(aes(x = intervals)) +
  geom_bar(aes(fill = sex), position = "fill") +
  ggtitle("Sex") +
  ylab("Percentage (%)") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Accent") +
  theme_bw() +
  annotate("point", x = 1, y = -0.08, shape = 15, colour = "#f1a340", size = 5, alpha = 1) +
  annotate("point", x = 2, y = -0.08, shape = 15, colour = "#f1a340", size = 5, alpha = 0.6) +
  annotate("point", x = 3, y = -0.08, shape = 15, colour = "#f1a340", size = 5, alpha = 0.2) +
  annotate("point", x = 4, y = -0.08, shape = 15, colour = "dodgerblue2", size = 5, alpha = 0.2) +
  annotate("point", x = 5, y = -0.08, shape = 15, colour = "dodgerblue2", size = 5, alpha = 0.6) +
  annotate("point", x = 6, y = -0.08, shape = 15, colour = "dodgerblue2", size = 5, alpha = 1) +
  coord_cartesian(ylim = c(0, 1), clip = "off") +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 13),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 13),
        legend.text = element_text(size= 13 ),
        legend.title = element_text(size = 13), 
        axis.title.x = element_text(vjust = -5, colour = "white"),
        strip.background = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5)) +
  guides(fill=guide_legend(title="Sex", title.position="top", title.hjust = 0.5))

plot_prehba1c_strata <- group.full.cohort %>%
  select(intervals, prehba1c) %>%
  drop_na() %>%
  ggplot(aes(x = intervals, y = prehba1c)) +
  geom_point(aes(fill = intervals, colour = intervals), alpha = 0, shape = 15, size = 5) +
  geom_boxplot(aes(fill = intervals), alpha = c(1, 0.6, 0.2, 0.2, 0.6, 1), show.legend = F) +
  ggtitle("Baseline HbA1c") +
  ylab("Baseline HbA1c (mmol/mol)") +
  theme_bw() +
  annotate("point", x = 1, y = 40, shape = 15, colour = "#f1a340", size = 5, alpha = 1) +
  annotate("point", x = 2, y = 40, shape = 15, colour = "#f1a340", size = 5, alpha = 0.6) +
  annotate("point", x = 3, y = 40, shape = 15, colour = "#f1a340", size = 5, alpha = 0.2) +
  annotate("point", x = 4, y = 40, shape = 15, colour = "dodgerblue2", size = 5, alpha = 0.2) +
  annotate("point", x = 5, y = 40, shape = 15, colour = "dodgerblue2", size = 5, alpha = 0.6) +
  annotate("point", x = 6, y = 40, shape = 15, colour = "dodgerblue2", size = 5, alpha = 1) +
  coord_cartesian(ylim = c(50, 180), clip = "off") +
  scale_y_continuous(breaks = c(53, 64, 75, 86, 97, 108, 119, 150)) +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 13),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 13),
        legend.text = element_text(size= 13 ),
        legend.title = element_text(size = 13), 
        axis.title.x = element_text(vjust = -5, colour = "white"),
        strip.background = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5),
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
    fill = guide_legend(title.position="top", title.hjust = 0.5, override.aes = list(alpha = c(1, 0.6, 0.2, 1, 0.6, 0.2)), shape = c(15,15,15,15,15,15),  ncol = 2, nrow = 3),
    colour = guide_legend(title.position="top", title.hjust = 0.5, override.aes = list(alpha = c(1, 0.6, 0.2, 1, 0.6, 0.2)), shape = c(15,15,15,15,15,15),  ncol = 2, nrow = 3)
  )


plot_preegfr_strata <- group.full.cohort %>%
  select(intervals, preegfr) %>%
  drop_na() %>%
  ggplot(aes(x = intervals, y = preegfr)) +
  geom_point(aes(fill = intervals, colour = intervals), alpha = 0, shape = 15, size = 5) +
  geom_boxplot(aes(fill = intervals), alpha = c(1, 0.6, 0.2, 0.2, 0.6, 1), show.legend = F) +
  ggtitle("eGFR") +
  ylab(expression(paste("eGFR (mL/min/1.3", m^{2}, ")"))) +
  scale_y_continuous(breaks = c(30, 60, 90, 120, 150)) +
  theme_bw() +
  annotate("point", x = 1, y = 9, shape = 15, colour = "#f1a340", size = 5, alpha = 1) +
  annotate("point", x = 2, y = 9, shape = 15, colour = "#f1a340", size = 5, alpha = 0.6) +
  annotate("point", x = 3, y = 9, shape = 15, colour = "#f1a340", size = 5, alpha = 0.2) +
  annotate("point", x = 4, y = 9, shape = 15, colour = "dodgerblue2", size = 5, alpha = 0.2) +
  annotate("point", x = 5, y = 9, shape = 15, colour = "dodgerblue2", size = 5, alpha = 0.6) +
  annotate("point", x = 6, y = 9, shape = 15, colour = "dodgerblue2", size = 5, alpha = 1) +
  coord_cartesian(ylim = c(20, 160), clip = "off") +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 13),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 13),
        legend.text = element_text(size= 13 ),
        legend.title = element_text(size = 13), 
        axis.title.x = element_text(vjust = -5, colour = "white"),
        strip.background = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5)) +
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
    fill = guide_legend(title.position="top", title.hjust = 0.5, override.aes = list(alpha = c(1, 0.6, 0.2, 1, 0.6, 0.2)), shape = c(15,15,15,15,15,15),  ncol = 2, nrow = 3),
    colour = guide_legend(title.position="top", title.hjust = 0.5, override.aes = list(alpha = c(1, 0.6, 0.2, 1, 0.6, 0.2)), shape = c(15,15,15,15,15,15),  ncol = 2, nrow = 3)
  )




plot_agetx_strata <- group.full.cohort %>%
  select(intervals, agetx) %>%
  drop_na() %>%
  ggplot(aes(x = intervals, y = agetx)) +
  geom_point(aes(fill = intervals, colour = intervals), alpha = 0, shape = 15, size = 5) +
  geom_boxplot(aes(fill = intervals), alpha = c(1, 0.6, 0.2, 0.2, 0.6, 1), show.legend = F) +
  ggtitle("Current age") +
  ylab("Current age (years)") +
  theme_bw() +
  annotate("point", x = 1, y = 5, shape = 15, colour = "#f1a340", size = 5, alpha = 1) +
  annotate("point", x = 2, y = 5, shape = 15, colour = "#f1a340", size = 5, alpha = 0.6) +
  annotate("point", x = 3, y = 5, shape = 15, colour = "#f1a340", size = 5, alpha = 0.2) +
  annotate("point", x = 4, y = 5, shape = 15, colour = "dodgerblue2", size = 5, alpha = 0.2) +
  annotate("point", x = 5, y = 5, shape = 15, colour = "dodgerblue2", size = 5, alpha = 0.6) +
  annotate("point", x = 6, y = 5, shape = 15, colour = "dodgerblue2", size = 5, alpha = 1) +
  coord_cartesian(ylim = c(12, 100), clip = "off") +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 13),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 13),
        legend.text = element_text(size= 13 ),
        legend.title = element_text(size = 13), 
        axis.title.x = element_text(vjust = -5, colour = "white"),
        strip.background = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5)) +
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
    fill = guide_legend(title.position="top", title.hjust = 0.5, override.aes = list(alpha = c(1, 0.6, 0.2, 1, 0.6, 0.2)), shape = c(15,15,15,15,15,15),  ncol = 2, nrow = 3),
    colour = guide_legend(title.position="top", title.hjust = 0.5, override.aes = list(alpha = c(1, 0.6, 0.2, 1, 0.6, 0.2)), shape = c(15,15,15,15,15,15),  ncol = 2, nrow = 3)
  )

plot_prebmi_strata <- group.full.cohort %>%
  select(intervals, prebmi) %>%
  drop_na() %>%
  ggplot(aes(x = intervals, y = prebmi)) +
  geom_point(aes(fill = intervals, colour = intervals), alpha = 0, shape = 15, size = 5) +
  geom_boxplot(aes(fill = intervals), alpha = c(1, 0.6, 0.2, 0.2, 0.6, 1), show.legend = F) +
  ggtitle("BMI") +
  ylab(expression(paste("BMI (kg/", m^{2}, ")"))) +
  theme_bw() +
  annotate("point", x = 1, y = 8, shape = 15, colour = "#f1a340", size = 5, alpha = 1) +
  annotate("point", x = 2, y = 8, shape = 15, colour = "#f1a340", size = 5, alpha = 0.6) +
  annotate("point", x = 3, y = 8, shape = 15, colour = "#f1a340", size = 5, alpha = 0.2) +
  annotate("point", x = 4, y = 8, shape = 15, colour = "dodgerblue2", size = 5, alpha = 0.2) +
  annotate("point", x = 5, y = 8, shape = 15, colour = "dodgerblue2", size = 5, alpha = 0.6) +
  annotate("point", x = 6, y = 8, shape = 15, colour = "dodgerblue2", size = 5, alpha = 1) +
  coord_cartesian(ylim = c(15, 100), clip = "off") +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 13),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 13),
        legend.text = element_text(size= 13 ),
        legend.title = element_text(size = 13), 
        axis.title.x = element_text(vjust = -5, colour = "white"),
        strip.background = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5)) +
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
    fill = guide_legend(title.position="top", title.hjust = 0.5, override.aes = list(alpha = c(1, 0.6, 0.2, 1, 0.6, 0.2)), shape = c(15,15,15,15,15,15),  ncol = 2, nrow = 3),
    colour = guide_legend(title.position="top", title.hjust = 0.5, override.aes = list(alpha = c(1, 0.6, 0.2, 1, 0.6, 0.2)), shape = c(15,15,15,15,15,15),  ncol = 2, nrow = 3)
  )


plot_main_2.1 <- patchwork::wrap_plots(list(
  plot_sex_strata, 
  plot_prehba1c_strata, 
  plot_preegfr_strata, 
  plot_agetx_strata, 
  plot_prebmi_strata
), nrow = 1) +
  plot_layout(guides = 'collect')  & theme(legend.position = "bottom")




plot_main_2.1 <- patchwork::wrap_plots(list(
  plot_sex_strata +
    theme(
      legend.position = "bottom"
    ), 
  
  patchwork::wrap_plots(list(
    
    plot_prehba1c_strata, 
    plot_preegfr_strata, 
    plot_agetx_strata, 
    plot_prebmi_strata
    
  ), nrow = 1) +
    plot_layout(guides = "collect") &
    theme(
      legend.position = "bottom"
    )
), nrow = 1) +
  plot_layout(
    widths = c(1, 5)
  )


pdf(width = 17, height = 8, "Plots/Paper/Manuscript/11.08.plot_main_2.1.pdf")
plot_main_2.1
dev.off()





#:------------------------------------------------------------------------------
# Figure 3: 6 panel plot: 
#   -A. Predicted HbA1c change for each clinical subgroup
#   -B. Predicted weight change for each clinical subgroup
#   -C. Discontinuation levels for each clinical subgroup
#   -D. Risk of microvascular complications
#   -E. Risk of CVD/MACE
#   -F. Risk of heart failure

require(forestplot)
require(grid)

### A - HbA1c

# levels needed
interval_breaks <- c(-5, -3, 0, 3, 5)

## Read in propensity scores
patient_prop_scores <- readRDS("Samples/SGLT2-GLP1/Aurum/ps_model/patient_prop_scores.rds")

## Read in treatment effects
treatment_effects <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/patient_effects.rds")

# Overall population
hba1c <- set_up_data_sglt2_glp1(dataset.type = "full.cohort") %>%
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
plot_a <- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = 50, lci = 50, uci = 50, drugclass = "SGLT2"),
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = 50, lci = 50, uci = 50, drugclass = "GLP1"),
  predictions_hba1c_stan_adjusted_overall %>% select(intervals, drugclass, mean, lci, uci) %>% slice(1:6),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = 50, lci = 50, uci = 50, drugclass = "SGLT2"),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = 50, lci = 50, uci = 50, drugclass = "GLP1"),
  predictions_hba1c_stan_adjusted_overall %>% select(intervals, drugclass, mean, lci, uci) %>% slice(7:12),
  predictions_hba1c_stan_adjusted_full %>% select(drugclass, mean, lci, uci) %>% mutate(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.hba1c.dataset$intervals)[1], paste0(">5 mmol/mol (n=", format(group.hba1c.dataset%>%filter(intervals == levels(group.hba1c.dataset$intervals)[1])%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                            ifelse(intervals == levels(group.hba1c.dataset$intervals)[2], paste0("3-5 mmol/mol (n=", format(group.hba1c.dataset%>%filter(intervals == levels(group.hba1c.dataset$intervals)[2])%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                   ifelse(intervals == levels(group.hba1c.dataset$intervals)[3], paste0("0-3 mmol/mol (n=", format(group.hba1c.dataset%>%filter(intervals == levels(group.hba1c.dataset$intervals)[3])%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                          ifelse(intervals == levels(group.hba1c.dataset$intervals)[4], paste0("0-3 mmol/mol (n=", format(group.hba1c.dataset%>%filter(intervals == levels(group.hba1c.dataset$intervals)[4])%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                 ifelse(intervals == levels(group.hba1c.dataset$intervals)[5], paste0("3-5 mmol/mol (n=", format(group.hba1c.dataset%>%filter(intervals == levels(group.hba1c.dataset$intervals)[5])%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                        ifelse(intervals == levels(group.hba1c.dataset$intervals)[6], paste0(">5 mmol/mol (n=", format(group.hba1c.dataset%>%filter(intervals == levels(group.hba1c.dataset$intervals)[6])%>%nrow(),big.mark=",",scientific=FALSE), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  mutate(group = ifelse(group == "SGLT2", "SGLT2i", "GLP1-RA")) %>%
  forestplot(ci.vertices = TRUE,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices.height = 0.1,
             title = "Average HbA1c change",
             clip = c(hba1c_strata_axis_min, hba1c_strata_axis_max),
             xticks = seq(hba1c_strata_axis_min, hba1c_strata_axis_max, 5),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Average HbA1c change (mmol/mol)") %>%
  fp_add_header(paste0("Overall population (n=", format(group.hba1c.dataset%>%nrow(),big.mark=",",scientific=FALSE), ")")) %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE)) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = "black")


### B - Weight

weight.dataset <- set_up_data_sglt2_glp1(dataset.type = "weight.dataset") %>%
  left_join(patient_prop_scores, by = c("patid", "pated")) %>%
  left_join(treatment_effects, by = c("patid", "pated")) %>%
  mutate(w.change = postweight - preweight)


group.weight.dataset <- group_values(data = weight.dataset,
                                     variable = "effects",
                                     breaks = interval_breaks) %>%
  drop_na(intervals)

# predictions for the weight adjusted model
predictions_weight_stan_adjusted_overall <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_weight_stan_adjusted_overall.rds")

# predictions for the weight adjusted model full
predictions_weight_stan_adjusted_full <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_weight_stan_adjusted_full.rds")


weight_strata_axis_min <- plyr::round_any(floor(min(c(predictions_weight_stan_adjusted_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min(),
                                                      predictions_weight_stan_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min()))), 1, f = floor)

weight_strata_axis_max <- plyr::round_any(ceiling(max(c(predictions_weight_stan_adjusted_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max(),
                                                        predictions_weight_stan_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max()))), 2, f = ceiling)


plot_b <- rbind(
  cbind(mean = 50, lci = 50, uci = 50, drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(mean = 50, lci = 50, uci = 50, drugclass = "GLP1", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_weight_stan_adjusted_overall %>%
    slice(c(1:6)),
  cbind(mean = 50, lci = 50, uci = 50, drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(mean = 50, lci = 50, uci = 50, drugclass = "GLP1", intervals = "Predicted HbA1c benefit on GLP1-RA"),
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
  mutate(group = ifelse(group == "SGLT2", "SGLT2i", "GLP1-RA")) %>%
  forestplot(ci.vertices = TRUE,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices.height = 0.1,
             title = "Average weight change",
             clip = c(weight_strata_axis_min, weight_strata_axis_max),
             xticks = seq(weight_strata_axis_min, weight_strata_axis_max, 1),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Average weight change (kg)") %>%
  fp_add_header(paste0("Overall population (n=", format(group.weight.dataset%>%nrow(),big.mark=",",scientific=FALSE), ")")) %>%
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


# predictions for the discontinuation adjusted model
predictions_discontinuation_stan_adjusted_overall <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_discontinuation_stan_adjusted_overall.rds")

# predictions for the discontinuation model sex strata
predictions_discontinuation_stan_adjusted_full <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_discontinuation_stan_adjusted_full.rds")

discontinuation_overall_axis_max <- plyr::round_any(ceiling(max(c(predictions_discontinuation_stan_adjusted_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() * 100,
                                                                  predictions_discontinuation_stan_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() * 100))), 5, f = ceiling)

discontinuation_strata_axis_max <- plyr::round_any(ceiling(max(c(predictions_discontinuation_stan_adjusted_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() * 100,
                                                                 predictions_discontinuation_stan_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() * 100))), 5, f = ceiling)


plot_c <- rbind(
  cbind(mean = 150, lci = 150, uci = 150, drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(mean = 150, lci = 150, uci = 150, drugclass = "GLP1", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_discontinuation_stan_adjusted_overall %>%
    slice(c(1:6)) %>%
    mutate(mean = as.numeric(mean)*100,
           lci = as.numeric(lci)*100,
           uci = as.numeric(uci)*100),
  cbind(mean = 150, lci = 150, uci = 150, drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(mean = 150, lci = 150, uci = 150, drugclass = "GLP1", intervals = "Predicted HbA1c benefit on GLP1-RA"),
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
  mutate(group = ifelse(group == "SGLT2", "SGLT2i", "GLP1-RA")) %>%
  forestplot(ci.vertices = TRUE,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices.height = 0.1,
             title = "Risk of discontinuation",
             clip = c(0, discontinuation_overall_axis_max),
             xticks = seq(0, discontinuation_overall_axis_max, 5),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Discontinuation (%)") %>%
  fp_add_header(paste0("Overall population (n=", format(group.discontinuation.dataset%>%nrow(),big.mark=",",scientific=FALSE), ")")) %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE)) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = "black")


### D - Risk of developing microvascular complications

micro_comp.dataset <- set_up_data_sglt2_glp1(dataset.type="micro_comp.dataset") %>%
  left_join(treatment_effects, by = c("patid", "pated"))

group.micro_comp.dataset <- group_values(data = micro_comp.dataset,
                                         variable = "effects",
                                         breaks = interval_breaks) %>%
  drop_na(intervals)

# predictions for the CVD outcomes in the population with no CVD/HF/CKD
predictions_micro_comp_micro_comp_stan_adjusted_overall <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_micro_comp_micro_comp_stan_adjusted_overall.rds")

# predictions for the CVD outcomes in the population with no CVD/HF/CKD
predictions_micro_comp_micro_comp_stan_adjusted_full <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_micro_comp_micro_comp_stan_adjusted_full.rds")

# plot of forest plot
plot_d <- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_micro_comp_micro_comp_stan_adjusted_overall %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = NA, lci = NA, uci = NA),
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
             title = "Microvascular complications",
             xlog = TRUE,
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazard ratio (95% CI, log scale)") %>%
  fp_add_header(paste0("Overall population (n=", format(nrow(group.micro_comp.dataset),big.mark=",",scientific=FALSE), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = "black")


### E - Risk of developing MACE/CVD

no_co.dataset <- set_up_data_sglt2_glp1(dataset.type="no_co.dataset") %>%
  left_join(treatment_effects, by = c("patid", "pated"))

group.no_co.dataset <- group_values(data = no_co.dataset,
                                    variable = "effects",
                                    breaks = interval_breaks) %>%
  drop_na(intervals)

# predictions for the CVD outcomes in the population with no CVD/HF/CKD
predictions_no_co_cvd_stan_adjusted_overall <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_no_co_cvd_stan_adjusted_overall.rds")

# predictions for the CVD outcomes in the population with no CVD/HF/CKD
predictions_no_co_cvd_stan_adjusted_full <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_no_co_cvd_stan_adjusted_full.rds")

# plot of forest plot
plot_e <- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_no_co_cvd_stan_adjusted_overall %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = NA, lci = NA, uci = NA),
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
             title = "Major adverse cardiovascular events",
             xlog = TRUE,
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazard ratio (95% CI, log scale)") %>%
  fp_add_header(paste0("Overall population (n=", format(nrow(group.no_co.dataset),big.mark=",",scientific=FALSE), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = "black")


### F - Risk of developing Heart Failure

no_co.dataset <- set_up_data_sglt2_glp1(dataset.type="no_co.dataset") %>%
  left_join(treatment_effects, by = c("patid", "pated"))

group.no_co.dataset <- group_values(data = no_co.dataset,
                                    variable = "effects",
                                    breaks = interval_breaks) %>%
  drop_na(intervals)

# predictions for the HF outcomes in the population with no CVD/HF/CKD
predictions_no_co_hf_stan_adjusted_overall <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_no_co_hf_stan_adjusted_overall.rds")

# predictions for the CVD outcomes in the population with no CVD/HF/CKD
predictions_no_co_hf_stan_adjusted_full <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_no_co_hf_stan_adjusted_full.rds")

plot_f <- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_no_co_hf_stan_adjusted_overall %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = NA, lci = NA, uci = NA),
  predictions_no_co_hf_stan_adjusted_overall %>% slice(4:6),
  predictions_no_co_hf_stan_adjusted_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.no_co.dataset$intervals)[1], paste0(">5 mmol/mol (n=", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[1])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[1])%>%filter(postdrug_hf_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                            ifelse(intervals == levels(group.no_co.dataset$intervals)[2], paste0("3-5 mmol/mol (n=", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[2])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[2])%>%filter(postdrug_hf_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                   ifelse(intervals == levels(group.no_co.dataset$intervals)[3], paste0("0-3 mmol/mol (n=", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[3])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[3])%>%filter(postdrug_hf_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                          ifelse(intervals == levels(group.no_co.dataset$intervals)[4], paste0("0-3 mmol/mol (n=", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[4])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[4])%>%filter(postdrug_hf_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                 ifelse(intervals == levels(group.no_co.dataset$intervals)[5], paste0("3-5 mmol/mol (n=", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[5])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[5])%>%filter(postdrug_hf_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                        ifelse(intervals == levels(group.no_co.dataset$intervals)[6], paste0(">5 mmol/mol (n=", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[6])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[6])%>%filter(postdrug_hf_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"), intervals))))))) %>%
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
  fp_add_header(paste0("Overall population (n=", format(nrow(group.no_co.dataset),big.mark=",",scientific=FALSE), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = "black")



pdf(width = 20, height = 10, "Plots/Paper/Manuscript/11.08.plot_main_3.pdf")

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



#:-------------------------------------------------------------------
#
# Extended figures
#

#:------------------------------------------------------------------------------
# Extended Figure 2: 3 panel plot: 
#   -A. Variable importance linear projection
#   -B. Overall cohort histogram split by sex
#   -C. Differential treatment effect

require(grid)

### A - Variable importance

# load variables used in the BCF model
variables_tau <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/variables_tau.rds")
variables_mu <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/variables_mu.rds")

# load BCF model
bcf_model <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/bcf_model.rds")

# load development dataset
hba1c.train.cleaned_up <- set_up_data_sglt2_glp1(dataset.type = "hba1c.train") %>%
  # drop the variables with the most missingness
  select(-preacr, -preast, -prehaematocrit, -prehaemoglobin, -pretriglyceride) %>%
  # only complete cases
  drop_na() %>%
  as.data.frame() %>%
  select(patid, pated, posthba1cfinal, drugclass, unique(c(variables_mu, variables_tau))) %>%
  cbind(effects = colMeans(bcf_model$tau))


m1 <- rms::ols(effects ~ rms::rcs(agetx, 3) + sex + ncurrtx + rms::rcs(prehba1c, 3) + rms::rcs(prebmi, 3) + rms::rcs(preegfr, 3) + preheartfailure + preihd + preneuropathy + prepad + preretinopathy, data = hba1c.train.cleaned_up, x = TRUE, y = TRUE)

values <- plot(anova(m1), what = 'proportion R2')

plot_a <- as.data.frame(values) %>%
  cbind(variable = c("Number of other current glucose-lowering drugs", "Sex", "eGFR", "Current age", "BMI", "HbA1c", "Retinopathy", "Peripheral arterial disease", "Neuropathy", "Ischaemic heart disease", "Heart failure")) %>%
  mutate(variable = factor(variable),
         values = values * 100) %>%
  ggplot(aes(y = forcats::fct_reorder(variable, values), x = values)) +
  geom_segment(aes(x = 0, xend = values, yend = forcats::fct_reorder(variable, values)), linetype = "dashed") +
  geom_point(size = 2, colour = "black") +
  ggtitle("Relative importance for treatment effect heterogeneity") +
  xlab("Relative Importance (%)") +
  theme_bw() +
  theme(axis.text.y = element_text(angle = 45, face = "bold"),
        axis.title.y = element_blank())

# theme(axis.text = element_text(face="bold"))


### B - Overall cohort histogram split by sex

hba1c.heterogeneity <- set_up_data_sglt2_glp1(dataset.type = "full.cohort") %>%
  left_join(readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/patient_effects.rds"), by = c("patid", "pated")) %>%
  select(sex, effects) %>%
  drop_na() %>%
  rename("mean" = "effects")


percent_var = round(((hba1c.heterogeneity %>% filter(sex == "Female" & mean < 0) %>% nrow()) / hba1c.heterogeneity %>% filter(sex == "Female") %>% nrow())*100)

plot_b_1 <- hist_plot(hba1c.heterogeneity %>% filter(sex == "Female"), "Females", -15, 18) +
  theme(legend.direction = "horizontal",
        legend.position = "bottom",
        legend.box = "horizontal") +
  labs(subtitle = paste0("SGLT2i = ",
                         percent_var,
                         "%, GLP1-RA = ",
                         100-percent_var,
                         "%"))

percent_var = round(((hba1c.heterogeneity %>% filter(sex == "Male" & mean < 0) %>% nrow()) / hba1c.heterogeneity %>% filter(sex == "Male") %>% nrow())*100)

plot_b_2 <- hist_plot(hba1c.heterogeneity %>% filter(sex == "Male"), "Males", -15, 18) +
  theme(legend.direction = "horizontal",
        legend.position = "bottom",
        legend.box = "horizontal") +
  labs(subtitle = paste0("SGLT2i = ",
                         percent_var,
                         "%, GLP1-RA = ",
                         100-percent_var,
                         "%"))

### C - Differential treatment effects

# Full cohort for average values
full.cohort <- set_up_data_sglt2_glp1(dataset.type="hba1c.train") %>%
  left_join(readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/patient_effects.rds"), by = c("patid", "pated"))

levels(full.cohort$sex) <- c("Females", "Males")
levels(full.cohort$ncurrtx) <- c("0", "1", "2", "3", "4+")
# Stratify predicted treatment effects by variables
# sex
plot_sex_strata <- full.cohort %>%
  select(sex, effects) %>%
  ggplot(aes(x = sex, y = effects)) +
  geom_hline(aes(yintercept = 0), colour = "red") +
  geom_violin(width = 0.3) +
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.3) +
  ylim(-10, 10) +
  ggtitle("Sex") +
  ylab("Predicted treatment effects (mmol/mol)") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5))

# number of other current drugs
plot_ncurrtx_strata <- full.cohort %>%
  select(ncurrtx, effects, sex) %>%
  ggplot(aes(x = ncurrtx, y = effects)) +
  geom_hline(aes(yintercept = 0), colour = "red") +
  geom_violin() +
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar") +
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
plot_prepad_strata <- full.cohort %>%
  select(prepad, effects, sex) %>%
  ggplot(aes(x = prepad, y = effects)) +
  geom_hline(aes(yintercept = 0), colour = "red") +
  geom_violin(width = 0.5) +
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.5) +
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
plot_preihd_strata <- full.cohort %>%
  select(preihd, effects, sex) %>%
  ggplot(aes(x = preihd, y = effects)) +
  geom_hline(aes(yintercept = 0), colour = "red") +
  geom_violin(width = 0.5) +
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.5) +
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
plot_preneuropathy_strata <- full.cohort %>%
  select(preneuropathy, effects, sex) %>%
  ggplot(aes(x = preneuropathy, y = effects)) +
  geom_hline(aes(yintercept = 0), colour = "red") +
  geom_violin(width = 0.5) +
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.5) +
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
plot_preretinopathy_strata <- full.cohort %>%
  select(preretinopathy, effects, sex) %>%
  ggplot(aes(x = preretinopathy, y = effects)) +
  geom_hline(aes(yintercept = 0), colour = "red") +
  geom_violin(width = 0.5) +
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.5) +
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
plot_preheartfailure_strata <- full.cohort %>%
  select(preheartfailure, effects, sex) %>%
  ggplot(aes(x = preheartfailure, y = effects)) +
  geom_hline(aes(yintercept = 0), colour = "red") +
  geom_violin(width = 0.5) +
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.5) +
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
breaks_hba1c <- quantile(full.cohort$prehba1c, probs = seq(0.2, 0.9, 0.2), na.rm = TRUE)

plot_prehba1c_strata <- group_values(data = full.cohort,
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
  geom_violin() +
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar") +
  facet_wrap(~sex) +
  ylim(-10, 10) +
  ggtitle("HbA1c (mmol/mol)") +
  ylab("Predicted treatment effects (mmol/mol)") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(fill="white"))

# eGFR at baseline
breaks_egfr <- quantile(full.cohort$preegfr, probs = seq(0.2, 0.9, 0.2), na.rm = TRUE)

plot_preegfr_strata <- group_values(data = full.cohort,
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
  geom_violin() +
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar") +
  facet_wrap(~sex) +
  ylim(-10, 10) +
  ggtitle(expression(paste("eGFR (mL/min per 1.3",m^{2},")"))) +
  ylab("Predicted treatment effects (mmol/mol)") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(fill="white"))


# Age at baseline
breaks_agetx <- quantile(full.cohort$agetx, probs = seq(0.2, 0.9, 0.2), na.rm = TRUE)

plot_agetx_strata <- group_values(data = full.cohort,
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
  geom_violin() +
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar") +
  facet_wrap(~sex) +
  ylim(-10, 10) +
  ggtitle("Current age (years)") +
  ylab("Predicted treatment effects (mmol/mol)") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(fill="white"))

# BMI at baseline
breaks_prebmi <- quantile(full.cohort$prebmi, probs = seq(0.2, 0.9, 0.2), na.rm = TRUE)

plot_prebmi_strata <- group_values(data = full.cohort,
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
  geom_violin() +
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar") +
  facet_wrap(~sex) +
  ylim(-10, 10) +
  ggtitle(expression(paste("BMI (kg/", m^{2}, ")"))) +
  ylab("Predicted treatment effects (mmol/mol)") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(fill="white"))


pdf(width = 17, height = 12, "Plots/Paper/Manuscript/11.08.plot_extended_2.pdf")

grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 5,
                                           ncol = 12, heights = unit(c(1.5, 0.2, 1, 1, 1), "null"))))

# first plot
pushViewport(viewport(layout.pos.row = 1,
                      layout.pos.col = 1:6))
grid.draw(ggplotGrob(plot_a))
upViewport()

# second plot
pushViewport(viewport(layout.pos.row = 1,
                      layout.pos.col = 7:9))
grid.draw(ggplotGrob(plot_b_1))
upViewport()

# third plot
pushViewport(viewport(layout.pos.row = 1,
                      layout.pos.col = 10:12))
grid.draw(ggplotGrob(plot_b_2))
upViewport()

# forth plot
pushViewport(viewport(layout.pos.row = 3,
                      layout.pos.col = 1:3))
grid.draw(ggplotGrob(plot_ncurrtx_strata +
                       ylab("") +
                       theme(axis.title.y = element_text())))
upViewport()

# fifth plot
pushViewport(viewport(layout.pos.row = 3,
                      layout.pos.col = 4:6))
grid.draw(ggplotGrob(plot_sex_strata))
upViewport()

# sixth plot
pushViewport(viewport(layout.pos.row = 3,
                      layout.pos.col = 7:9))
grid.draw(ggplotGrob(plot_preegfr_strata))
upViewport()

# seventh plot
pushViewport(viewport(layout.pos.row = 3,
                      layout.pos.col = 10:12))
grid.draw(ggplotGrob(plot_agetx_strata))
upViewport()

# nineth plot
pushViewport(viewport(layout.pos.row = 4,
                      layout.pos.col = 4:6))
grid.draw(ggplotGrob(plot_prehba1c_strata))
upViewport()

# tenth plot
pushViewport(viewport(layout.pos.row = 4,
                      layout.pos.col = 7:9))
grid.draw(ggplotGrob(plot_preretinopathy_strata))
upViewport()

# eleventh plot
pushViewport(viewport(layout.pos.row = 4,
                      layout.pos.col = 10:12))
grid.draw(ggplotGrob(plot_prepad_strata))
upViewport()

# twelfth plot
pushViewport(viewport(layout.pos.row = 5,
                      layout.pos.col = 1:3))
grid.draw(ggplotGrob(plot_preneuropathy_strata +
                       ylab("") +
                       theme(axis.title.y = element_text())))
upViewport()

# thirteenth plot
pushViewport(viewport(layout.pos.row = 5,
                      layout.pos.col = 4:6))
grid.draw(ggplotGrob(plot_preihd_strata))
upViewport()

# fourteenth plot
pushViewport(viewport(layout.pos.row = 5,
                      layout.pos.col = 7:9))
grid.draw(ggplotGrob(plot_preheartfailure_strata))
upViewport()

# eighth plot
pushViewport(viewport(layout.pos.row = 4,
                      layout.pos.col = 1:3))
grid.draw(ggplotGrob(plot_prebmi_strata +
                       theme(axis.title.y = element_text(size = 12))))
upViewport()

# legend
grid.text(expression(bold("A")), vp = viewport(layout.pos.row = 1, layout.pos.col = 1), hjust = 5, vjust = -14.7)
grid.text(expression(bold("B.1")), vp = viewport(layout.pos.row = 1, layout.pos.col = 7), hjust = 3, vjust = -14.7)
grid.text(expression(bold("B.2")), vp = viewport(layout.pos.row = 1, layout.pos.col = 10), hjust = 3, vjust = -14.7)
grid.text(expression(bold("C")), vp = viewport(layout.pos.row = 2, layout.pos.col = 1), hjust = 5, vjust = 2)

dev.off()



#:--------------------------------------------------------------------------------
# Extended Figure 3:
### Sex difference in trials

## glp1 order: Liraglutide, albi, priba, CPRD, scotland

## sglt2 order: trials, cprd, scotland 

# create dataset
df <- cbind(
  dataset = c("CPRD (n=11,296)", "PRIBA (n=550)", "Harmony 7 RCT Liraglutide (n=389)", "Harmony RCT (pooled) Albiglutide (n=1,682)", "CPRD (n=40,915)", "Trials meta-analysis (n=7,119)", "Tayside & Fife (n=1,837)", "Tayside & Fife (n=415)"),
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
plot_extended_3.1 <- df %>%
  filter(drug == "SGLT2-inhibitors") %>%
  mutate(dataset = factor(dataset, levels = c("Tayside & Fife (n=1,837)", "CPRD (n=40,915)", "Trials meta-analysis (n=7,119)"))) %>%
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

plot_extended_3.2 <- df %>%
  filter(drug == "GLP-1 receptor agonists") %>%
  mutate(dataset = factor(dataset, levels = c("Tayside & Fife (n=415)", "CPRD (n=11,296)", "PRIBA (n=550)", "Harmony RCT (pooled) Albiglutide (n=1,682)", "Harmony 7 RCT Liraglutide (n=389)"))) %>%
  ggplot(aes(x = beta, y = dataset)) +
  geom_vline(aes(xintercept = 0), colour = "black") +
  geom_point(colour = "dodgerblue2") +
  geom_errorbar(aes(xmin = lci, xmax = uci), colour = "dodgerblue2", width = 0.3) +
  ggtitle("GLP1-receptor agonists") +
  xlab("HbA1c response difference in males versus females (mmol/mol)") +
  xlim(min(df%>%select(beta, lci, uci)%>%unlist()), max(df%>%select(beta, lci, uci)%>%unlist())) +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_line(colour = "white"))

plot_extended_3 <- patchwork::wrap_plots(list(plot_extended_3.1, plot_extended_3.2), ncol = 1, nrow = 2) +
  plot_annotation(tag_levels = c("A"))


# pdf of plot
pdf(width = 8, height = 2, "Plots/Paper/Sup_Mat/11.08plot_extended_3.1.pdf")
plot_extended_3.1
dev.off()

pdf(width = 8, height = 3, "Plots/Paper/Sup_Mat/11.08.plot_extended_3.2.pdf")
plot_extended_3.2
dev.off()



#:--------------------------------------------------------------------------------
# Extended Figure 4:
### Figure of treatment effects validation for Ethnicity

# intervals
interval_breaks <- c(-5, -3, 0, 3, 5)

## Read in propensity scores
patient_prop_scores <- readRDS(paste0(output_path, "/ps_model/patient_prop_scores.rds"))

## Read in treatment effects
treatment_effects <- readRDS(paste0(output_path, "/response_model_bcf/patient_effects.rds"))

# load in variables used in the model
variables_mu <- readRDS(paste0(output_path, "/response_model_bcf/variables_mu.rds"))

variables_tau <- readRDS(paste0(output_path, "/response_model_bcf/variables_tau.rds"))

# Create dataset with all the variables needed
ethnicity.dataset <- set_up_data_sglt2_glp1(dataset.type = "hba1c.train") %>%
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
  left_join(patient_prop_scores, by = c("patid", "pated")) %>%
  left_join(treatment_effects, by = c("patid", "pated")) %>%
  left_join(set_up_data_sglt2_glp1(dataset.type = "full.cohort") %>%
              select(patid, pated, ethnicity), by = c("patid", "pated")) %>%
  select(c("patid", "pated", "ethnicity", "drugclass", "posthba1cfinal", unique(c(variables_mu, variables_tau)), "prop.score", "effects")) %>%
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

# Subsetting data ethnicity - Non-white (South Asian, Black, Other/Mixed)
predicted_observed_non_white <- ethnicity.dataset.grouped %>%
  filter(ethnicity != "White")


#### ---------------------- Propensity score matching (too low numbers)

#### ---------------------- Propensity score matching + adjusted (too low numbers)

#### ---------------------- Adjusted
# White
ATE_adjust_white <- calc_ATE(data = predicted_observed_white%>%mutate(intervals=as.numeric(intervals)), validation_type = "Adjust",
                             variable = "posthba1cfinal", quantile_var = "intervals",
                             order = "largest", breakdown = unique(c(variables_mu, variables_tau)))

# White
ATE_adjust_non_white <- calc_ATE(data = predicted_observed_non_white%>%mutate(intervals=as.numeric(intervals)), validation_type = "Adjust",
                                 variable = "posthba1cfinal", quantile_var = "intervals",
                                 order = "largest", breakdown = unique(c(variables_mu, variables_tau)))


### Overall analysis (1 group)
# White
ATE_adjust_white_overall <- calc_ATE(data = predicted_observed_white%>%mutate(intervals=as.numeric(1)), validation_type = "Adjust",
                                     variable = "posthba1cfinal", quantile_var = "intervals",
                                     order = "largest", breakdown = unique(c(variables_mu, variables_tau)))

# Non-White
ATE_adjust_non_white_overall <- calc_ATE(data = predicted_observed_non_white%>%mutate(intervals=as.numeric(1)), validation_type = "Adjust",
                                         variable = "posthba1cfinal", quantile_var = "intervals",
                                         order = "largest", breakdown = unique(c(variables_mu, variables_tau)))

#:-----------------------------------------------------------------------------------
# Forest plots for calibration variable

# Set up axis limits for the forest plots
hba1c_strata_axis_min <- plyr::round_any(floor(min(c(
  ATE_adjust_white[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
  ATE_adjust_non_white[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
  ATE_adjust_white_overall[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
  ATE_adjust_non_white_overall[["effects"]] %>% select(c("obs","lci","uci")) %>% min()
))), 2, f = floor)

hba1c_strata_axis_max <- plyr::round_any(ceiling(max(c(
  ATE_adjust_white[["effects"]] %>% select(c("obs","lci","uci")) %>% max(),
  ATE_adjust_non_white[["effects"]] %>% select(c("obs","lci","uci")) %>% max(),
  ATE_adjust_white_overall[["effects"]] %>% select(c("obs","lci","uci")) %>% max(),
  ATE_adjust_non_white_overall[["effects"]] %>% select(c("obs","lci","uci")) %>% max()
))), 2, f = ceiling)



plot_extended_4 <- rbind(
  cbind(obs = 50, lci = 50, uci = 50, group = "White", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(obs = 50, lci = 50, uci = 50, group = "Non-White", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(ATE_adjust_white[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(1:3), group = "White"),
  cbind(ATE_adjust_non_white[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(1:3), group = "Non-White"),
  cbind(obs = 50, lci = 50, uci = 50, group = "White", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(obs = 50, lci = 50, uci = 50, group = "Non-White", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(ATE_adjust_white[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(4:6), group = "White"),
  cbind(ATE_adjust_non_white[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(4:6), group = "Non-White"),
  cbind(ATE_adjust_white_overall[["effects"]] %>% select(obs, lci, uci), group = "White", intervals = "Average treatment effect"),
  cbind(ATE_adjust_non_white_overall[["effects"]] %>% select(obs, lci, uci), group = "Non-White", intervals = "Average treatment effect")
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
  mutate(n_white = c(rep(NA, 2), rep(ATE_adjust_white[["effects"]]$N[1:3], 2),
                     rep(NA, 2), rep(ATE_adjust_white[["effects"]]$N[4:6], 2),
                     rep(sum(ATE_adjust_white[["effects"]]$N), 2)),
         n_non_white = c(rep(NA, 2), rep(ATE_adjust_non_white[["effects"]]$N[1:3], 2),
                         rep(NA, 2), rep(ATE_adjust_non_white[["effects"]]$N[4:6], 2),
                         rep(sum(ATE_adjust_non_white[["effects"]]$N), 2)),
         mean = as.numeric(mean)) %>%
  group_by(group) %>%
  forestplot(labeltext = c("labeltext", "n_white", "n_non_white"),
             ci.vertices = TRUE,
             ci.vertices.height = 0.1,
             # title = "Adjusted",
             clip = c(hba1c_strata_axis_min, hba1c_strata_axis_max),
             xticks = seq(hba1c_strata_axis_min, hba1c_strata_axis_max, 2),
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Average treatment effect (mmol/mol)") %>%
  fp_add_header(labeltext = paste0("Overall population (n=", format(ethnicity.dataset.grouped%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                n_white = "White",
                n_non_white = "Non-White") %>%
  fp_set_style(box = c("#E64B35B2", "#3C5488B2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE)) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = "black")



pdf(width = 10, height = 6, "Plots/Paper/Sup_Mat/11.08.plot_extended_4.pdf")
plot_extended_4
dev.off()



